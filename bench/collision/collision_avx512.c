/*================================================================
   The JCP 2016 collision kernel (src/Collision_JCP2016.inc) written with
   AVX-512 intrinsics: one row of 8 cells per 512-bit vector, every
   multiply-add an explicit FMA, mask registers for the interface normal.

   The same body is instantiated four times (always_inline + constant
   arguments): with P, Ux, Uy as input or computed from g in the same pass
   (moments), and with exact divisions / square root or with the 14-bit
   reciprocal estimates of AVX-512 refined by two Newton steps (approx).

   Called from Fortran through BIND(C) interfaces (Collision_Variants_Mod):
   every argument by reference, arrays as in the Fortran kernel,

       C(-1:nx+2,-1:ny+2), P(nx,ny), Ux(nx,ny), Uy(nx,ny),
       h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

   par = [ Rhol, Rhoh, Gy, taul, tauh, M0, W, wc, Beta, kappa, grad2_min ].
   Compile with -DAMR_BLOCK_SIZE=<nx> (a multiple of 8) and -march=<a CPU
   with AVX-512F>; without AVX-512 only collision_avx512_available() = 0
   is built with stubs.
  ================================================================*/
#include <stddef.h>

#ifndef AMR_BLOCK_SIZE
#error "AMR_BLOCK_SIZE must be defined"
#endif
#define NX AMR_BLOCK_SIZE
#define NY AMR_BLOCK_SIZE

#if defined(__AVX512F__) && (NX % 8 == 0)
#include <immintrin.h>

int collision_avx512_available(void) { return 1; }

typedef __m512d V;
#define SET(x)       _mm512_set1_pd(x)
#define ADD(a,b)     _mm512_add_pd(a,b)
#define SUB(a,b)     _mm512_sub_pd(a,b)
#define MUL(a,b)     _mm512_mul_pd(a,b)
#define DIV(a,b)     _mm512_div_pd(a,b)
#define MAX(a,b)     _mm512_max_pd(a,b)
#define MIN(a,b)     _mm512_min_pd(a,b)
#define FMA(a,b,c)   _mm512_fmadd_pd(a,b,c)      /*   a*b + c  */
#define FMS(a,b,c)   _mm512_fmsub_pd(a,b,c)      /*   a*b - c  */
#define FNMA(a,b,c)  _mm512_fnmadd_pd(a,b,c)     /* -(a*b) + c */
#define LOAD(p)      _mm512_loadu_pd(p)
#define STORE(p,v)   _mm512_storeu_pd(p,v)

/* strides of the fields */
#define LDC  (NX+4)                 /* C: two ghost layers */
#define LDH  (NX+2)                 /* h, g: one ghost layer */
#define PLANE ((NX+2)*(NY+2))       /* one direction of h, g */

#define INLINE static inline __attribute__((always_inline))

/* 1/x and 1/sqrt(x): exact, or the 14-bit estimate with two Newton steps */
INLINE V recip(V x, int approx)
{
    if (!approx) return DIV(SET(1.0), x);
    if (approx == 2) return FMA(SET(1.e-3), x, SET(0.125));    /* no division at all: cost probe, wrong results */
    V y = _mm512_rcp14_pd(x);
    y = MUL(y, FNMA(x, y, SET(2.0)));
    y = MUL(y, FNMA(x, y, SET(2.0)));
    return y;
}
INLINE V rsqrt(V x, int approx)
{
    if (!approx) return DIV(SET(1.0), _mm512_sqrt_pd(x));
    if (approx == 2) return ADD(SET(1.e-3), x);
    V y  = _mm512_rsqrt14_pd(x);
    V hx = MUL(SET(0.5), x);
    y = MUL(y, FNMA(MUL(hx, y), y, SET(1.5)));
    y = MUL(y, FNMA(MUL(hx, y), y, SET(1.5)));
    return y;
}

/* a pair of opposite directions kp, km (1-based) with weight w, e.u = +-v:
     Gamma+- = Gs +- Ga,  h relaxes to C (Gamma+- +- sq),
     g+- = Pw + Gamma+- r+- + (Gamma+- amu - w dRho3) (B+- + D0)          */
#define PAIR(kp, km, w, w3, wdr, tw, v, vsq, nrm, PWM, rP, rM, Bp, Bm)              \
    {                                                                               \
        V Gs = MUL(w, FMA(c45, vsq, base0));                                        \
        V Ga = MUL(w3, v);                                                          \
        V Gp = ADD(Gs, Ga), Gm = SUB(Gs, Ga);                                       \
        V sq = MUL(tCH, MUL(tw, nrm));                                              \
        double *hp = hh + ((kp)-1)*PLANE, *hm = hh + ((km)-1)*PLANE;                \
        STORE(hp, FMA(omwc, LOAD(hp), MUL(cH, ADD(Gp, sq))));                       \
        STORE(hm, FMA(omwc, LOAD(hm), MUL(cH, SUB(Gm, sq))));                       \
        V Pw = PWM;                                                                 \
        STORE(gg + ((kp)-1)*PLANE, FMA(FMS(Gp, amu, wdr), ADD(Bp, D0), FMA(Gp, rP, Pw))); \
        STORE(gg + ((km)-1)*PLANE, FMA(FMS(Gm, amu, wdr), ADD(Bm, D0), FMA(Gm, rM, Pw))); \
    }

INLINE void kernel(const int moments, const int approx, const double *restrict inv_dx_p,
                   const double *restrict par, const double *restrict Cf, double *restrict Pf,
                   double *restrict Uxf, double *restrict Uyf, double *restrict hf, double *restrict gf)
{
    const double inv_dx = *inv_dx_p;
    const double Rhol = par[0], Rhoh = par[1], Gy = par[2], taul = par[3], tauh = par[4],
                 M0 = par[5], W = par[6], wc_ = par[7], Beta = par[8], kappa = par[9], g2min = par[10];
    const double w0_ = 4.0/9.0, ws_ = 1.0/9.0, wd_ = 1.0/36.0;
    const double dRho_ = Rhoh - Rhol, dRho3_ = dRho_/3.0;

    /* per block: broadcast constants (the compiler keeps the hot ones in registers
       and uses the others as embedded-broadcast memory operands) */
    const V cC = SET(0.5*inv_dx), cM = SET(0.25*inv_dx);
    const V one = SET(1.0), two = SET(2.0), three = SET(3.0), four = SET(4.0), half = SET(0.5);
    const V c13 = SET(1.0/3.0), c23 = SET(2.0/3.0), c16 = SET(1.0/6.0), c112 = SET(1.0/12.0);
    const V c15 = SET(1.5), c45 = SET(4.5), quarter = SET(0.25), twenty = SET(20.0), zero = _mm512_setzero_pd();
    const V rhol = SET(Rhol), rhoh = SET(Rhoh), dRho = SET(dRho_), dRho3 = SET(dRho3_), dRho6 = SET(dRho_/6.0);
    const V gy = SET(Gy), hGy = SET(0.5*Gy);
    const V invTaul = SET(1.0/taul), dInvTau = SET(1.0/tauh - 1.0/taul);
    const V tws = SET(12.0*M0/W*ws_), twd = SET(12.0*M0/W*wd_);
    const V wc = SET(wc_), omwc = SET(1.0 - wc_);
    const V w0 = SET(w0_), ws = SET(ws_), wd = SET(wd_);
    const V w0d = SET(w0_*dRho3_), wsd = SET(ws_*dRho3_), wdd = SET(wd_*dRho3_);
    const V ws3 = SET(3.0*ws_), wd3 = SET(3.0*wd_), ws45 = SET(4.5*ws_), wd18 = SET(18.0*wd_), wd12 = SET(12.0*wd_);
    const V c4B = SET(4.0*Beta), cK = SET(kappa*inv_dx*inv_dx/6.0), gmin = SET(g2min);

    for (int j = 1; j <= NY; j++)
    for (int i = 1; i <= NX; i += 8) {
        const double *c = Cf + (i+1) + LDC*(j+1);           /* C(i,j)   */
        double *hh = hf + i + LDH*j;                        /* h(i,j,1) */
        double *gg = gf + i + LDH*j;                        /* g(i,j,1) */
        const ptrdiff_t ip = (i-1) + NX*(j-1);              /* P(i,j)   */
#define Cat(di,dj) LOAD(c + (di) + LDC*(dj))

        const V C0  = Cat(0,0);
        const V Cp0 = Cat(+1,0), Cm0 = Cat(-1,0), C0p = Cat(0,+1), C0m = Cat(0,-1);
        const V Cpp = Cat(+1,+1), Cmm = Cat(-1,-1), Cmp = Cat(-1,+1), Cpm = Cat(+1,-1);

        /* central differences along e2, e3, e6, e7 */
        const V eC2 = MUL(SUB(Cp0, Cm0), cC), eC3 = MUL(SUB(C0p, C0m), cC);
        const V eC6 = MUL(SUB(Cpp, Cmm), cC), eC7 = MUL(SUB(Cmp, Cpm), cC);

        /* biased differences B_k = (4 C(x+e) - 3 C(x) - C(x+2e))/(4 dx) */
        const V C3 = MUL(three, C0);
#define BIAS(C1, di, dj) MUL(SUB(FMS(four, C1, C3), Cat(di,dj)), cM)
        const V B2 = BIAS(Cp0, +2, 0), B4 = BIAS(Cm0, -2, 0), B3 = BIAS(C0p, 0, +2), B5 = BIAS(C0m, 0, -2);
        const V B6 = BIAS(Cpp, +2, +2), B8 = BIAS(Cmm, -2, -2), B7 = BIAS(Cmp, -2, +2), B9 = BIAS(Cpm, +2, -2);

        /* isotropic gradients of the central and the mixed differences */
        const V DcDx = FMA(c23, eC2, MUL(c16, SUB(eC6, eC7)));
        const V DcDy = FMA(c23, eC3, MUL(c16, ADD(eC6, eC7)));
        const V t68 = ADD(SUB(B6, B8), eC6), t79 = ADD(SUB(B7, B9), eC7);
        const V DmDx = FMA(c13, ADD(SUB(B2, B4), eC2), MUL(c112, SUB(t68, t79)));
        const V DmDy = FMA(c13, ADD(SUB(B3, B5), eC3), MUL(c112, ADD(t68, t79)));

        /* unit normal, zero where |grad C|^2 <= grad2_min (zero-masking) */
        const V g2 = FMA(DcDx, DcDx, MUL(DcDy, DcDy));
        const __mmask8 kint = _mm512_cmp_pd_mask(g2, gmin, _CMP_GT_OQ);
        const V rn = _mm512_maskz_mov_pd(kint, rsqrt(MAX(g2, gmin), approx));
        const V ni = MUL(DcDx, rn), nj = MUL(DcDy, rn);

        /* chemical potential */
        V lap = FMA(four, ADD(ADD(C0m, Cm0), ADD(Cp0, C0p)), ADD(ADD(Cmm, Cpm), ADD(Cmp, Cpp)));
        lap = FNMA(twenty, C0, lap);
        const V xmu = FNMA(cK, lap, MUL(MUL(MUL(c4B, C0), SUB(C0, one)), SUB(C0, half)));

        const V Rho = FMA(C0, dRho, rhol);
        const V hF  = MUL(SUB(rhoh, Rho), hGy);             /* half the buoyancy */

        const V g2_ = LOAD(gg + 1*PLANE), g3_ = LOAD(gg + 2*PLANE), g4_ = LOAD(gg + 3*PLANE), g5_ = LOAD(gg + 4*PLANE);
        const V g6_ = LOAD(gg + 5*PLANE), g7_ = LOAD(gg + 6*PLANE), g8_ = LOAD(gg + 7*PLANE), g9_ = LOAD(gg + 8*PLANE);
        const V s24 = ADD(g2_, g4_), s35 = ADD(g3_, g5_), s68 = ADD(g6_, g8_), s79 = ADD(g7_, g9_);

        V vx, vy, Pc, uDc;
        if (moments) {
            /* pressure and velocity from g and the forces */
            const V S0 = ADD(LOAD(gg), ADD(ADD(s24, s35), ADD(s68, s79)));
            const V Sx = ADD(SUB(g2_, g4_), ADD(SUB(g6_, g7_), SUB(g9_, g8_)));
            const V Sy = ADD(SUB(g3_, g5_), SUB(ADD(g6_, g7_), ADD(g8_, g9_)));
            const V rRho = recip(Rho, approx);
            vx  = MUL(FMA(three, Sx, MUL(half, MUL(xmu, DcDx))), rRho);
            vy  = MUL(FMA(three, Sy, FMA(half, MUL(xmu, DcDy), hF)), rRho);
            uDc = FMA(vx, DcDx, MUL(vy, DcDy));
            Pc  = FMA(dRho6, uDc, S0);
            STORE(Uxf + ip, vx);
            STORE(Uyf + ip, vy);
            STORE(Pf  + ip, Pc);
        } else {
            vx  = LOAD(Uxf + ip);
            vy  = LOAD(Uyf + ip);
            Pc  = LOAD(Pf  + ip);
            uDc = FMA(vx, DcDx, MUL(vy, DcDy));
        }

        /* per-cell quantities */
        const V vxx = MUL(vx, vx), vyy = MUL(vy, vy);
        const V base0 = FNMA(c15, ADD(vxx, vyy), one);
        const V D0   = FNMA(vy, DmDy, FNMA(vx, DmDx, MUL(half, uDc)));
        const V rho3 = MUL(Rho, c13);
        const V amu  = ADD(dRho3, xmu);
        const V Pr   = SUB(Pc, rho3);
        const V r0   = FNMA(hF, vy, rho3), rp = ADD(r0, hF), rm = SUB(r0, hF);
        const V cH   = MUL(wc, C0), tCH = SUB(one, C0);

        /* (1 - s9)/4 = (2 - t)/(2 + t)/4, t = 1/tau_fluid interpolated in the clamped C */
        const V t  = FMA(MIN(MAX(C0, zero), one), dInvTau, invTaul);
        const V km = MUL(FNMA(quarter, t, half), recip(ADD(two, t), approx));

        /* relaxed stress moments of g - geq (pair sums of geq in closed form) */
        const V Q  = FMA(amu, uDc, MUL(two, FMA(hF, vy, rho3)));
        const V e8 = FMA(MUL(ws45, SUB(vxx, vyy)), Q,
                         MUL(ws3, FNMA(amu, FMS(vx, eC2, MUL(vy, eC3)), MUL(two, MUL(vy, hF)))));
        const V vp = ADD(vx, vy), vm = SUB(vy, vx);
        const V e9 = FNMA(MUL(wd3, amu), FMS(vp, eC6, MUL(vm, eC7)),
                          FMS(MUL(wd18, MUL(vx, vy)), Q, MUL(wd12, MUL(hF, vx))));
        const V m8 = MUL(km, SUB(SUB(s24, s35), e8));
        const V m9 = MUL(km, SUB(SUB(s68, s79), e9));

        /* rest direction */
        {
            const V Gp = MUL(w0, base0);
            STORE(hh, FMA(omwc, LOAD(hh), MUL(cH, Gp)));
            STORE(gg, FMA(FMS(Gp, amu, w0d), D0, FMA(Gp, r0, MUL(Pr, w0))));
        }
        PAIR(2, 4, ws, ws3, wsd, tws, vx, vxx,          ni,          FMA(Pr, ws, m8),  r0, r0, B2, B4)
        PAIR(3, 5, ws, ws3, wsd, tws, vy, vyy,          nj,          FMS(Pr, ws, m8),  rp, rm, B3, B5)
        PAIR(6, 8, wd, wd3, wdd, twd, vp, MUL(vp, vp),  ADD(ni, nj), FMA(Pr, wd, m9),  rp, rm, B6, B8)
        PAIR(7, 9, wd, wd3, wdd, twd, vm, MUL(vm, vm),  SUB(nj, ni), FMS(Pr, wd, m9),  rp, rm, B7, B9)
    }
    (void)gy;
}

#define EXPORT(name, moments, approx)                                                              \
    void name(const double *inv_dx, const double *par, const double *C, double *P, double *Ux,     \
              double *Uy, double *h, double *g)                                                    \
    { kernel(moments, approx, inv_dx, par, C, P, Ux, Uy, h, g); }

#else   /* no AVX-512: stubs, never called */

int collision_avx512_available(void) { return 0; }
#define EXPORT(name, moments, approx)                                                              \
    void name(const double *inv_dx, const double *par, const double *C, double *P, double *Ux,     \
              double *Uy, double *h, double *g)                                                    \
    { (void)inv_dx; (void)par; (void)C; (void)P; (void)Ux; (void)Uy; (void)h; (void)g; }

#endif

EXPORT(collision_jcp2016_avx512,             0, 0)
EXPORT(moments_collision_jcp2016_avx512,     1, 0)
EXPORT(collision_jcp2016_avx512_rcp,         0, 1)
EXPORT(moments_collision_jcp2016_avx512_rcp, 1, 1)
EXPORT(collision_jcp2016_avx512_nodiv,         0, 2)
EXPORT(moments_collision_jcp2016_avx512_nodiv, 1, 2)

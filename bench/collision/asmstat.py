#!/usr/bin/env python3
"""Instruction statistics of the functions in an assembly file (AT&T syntax, gfortran/ifx/flang -S).
   usage: asmstat.py variants.s [name-substring ...]
   Per function: instructions, packed/scalar FP instructions by register width, divisions and
   square roots, vector loads/stores to the stack (spills and reloads), gathers, calls."""
import re, sys
fn = sys.argv[1]; pats = [p.lower() for p in sys.argv[2:]]
funcs = {}; cur = None
for line in open(fn):
    s = line.strip()
    m = re.match(r'^([A-Za-z_.$][\w.$]*):', s)
    if m and not m.group(1).startswith('.'):
        cur = m.group(1); funcs[cur] = []
        continue
    if cur is None or not s or s.startswith(('.', '#')) or s.endswith(':'): continue
    funcs[cur].append(s)
print('%-44s %6s %6s %6s %6s %6s %5s %6s %6s %6s %5s' %
      ('function', 'instr', 'zmm', 'ymm', 'xmm_pd', 'scal', 'div', 'sqrt', 'spill', 'reload', 'call'))
for name, ins in funcs.items():
    if pats and not any(p in name.lower() for p in pats): continue
    if len(ins) < 30: continue
    c = dict(zmm=0, ymm=0, xmmpd=0, scal=0, div=0, sqrt=0, spill=0, reload=0, call=0, gather=0)
    for s in ins:
        op = s.split()[0]
        if op.startswith('call'): c['call'] += 1
        if not op.startswith('v') and not re.match(r'(add|sub|mul|div|sqrt|mov[au]p|movs)[sp]d', op): continue
        if 'zmm' in s: c['zmm'] += 1
        elif 'ymm' in s: c['ymm'] += 1
        elif re.search(r'pd\b|pd$', op) or op.endswith('pd'): c['xmmpd'] += 1
        elif op.endswith('sd') or 'sd' in op: c['scal'] += 1
        if 'div' in op: c['div'] += 1
        if 'sqrt' in op: c['sqrt'] += 1
        if 'gather' in op: c['gather'] += 1
        stack = re.search(r'\((%rsp|%rbp)[,)]', s)
        if stack and re.match(r'v?mov[au]p[ds]|vmovdq|v?movsd', op):
            ops = s.split(None, 1)[1]
            if re.search(r'\((%rsp|%rbp)[^)]*\)\s*$', ops) or re.search(r'\((%rsp|%rbp)[^)]*\)\{', ops): c['spill'] += 1
            else: c['reload'] += 1
        elif stack and op.startswith('v'):
            c['reload'] += 1          # stack operand used directly
    print('%-44s %6d %6d %6d %6d %6d %5d %6d %6d %6d %5d%s' % (name[:44], len(ins), c['zmm'], c['ymm'], c['xmmpd'],
          c['scal'], c['div'], c['sqrt'], c['spill'], c['reload'], c['call'], '  gather %d' % c['gather'] if c['gather'] else ''))

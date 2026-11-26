#!/bin/bash
# Run built configurations and tabulate best MLUPS:  ./run.sh <pool_MiB> <tag> [<tag> ...]
[ -f ~/fortran-env.sh ] && source ~/fortran-env.sh; set +u   # compiler run-time libraries (this sandbox)
# (single thread pinned to a core unless OMP_NUM_THREADS is set). Results go to results/<tag>_<pool>.txt
here=$(cd "$(dirname "$0")" && pwd); mkdir -p $here/results
mb=$1; shift
sweeps=$(python3 -c "print(max(4,int(240/$mb)))")
for tag in "$@"; do
  if [ -z "$OMP_NUM_THREADS" ]; then pin="taskset -c 2"; export OMP_NUM_THREADS=1; unset_after=1; else pin=""; unset_after=0; fi
  out=$here/results/${tag}_${mb}MiB_${OMP_NUM_THREADS}t.txt
  $pin $here/build/$tag/bench_collision $mb $sweeps ${TRIALS:-15} > $out 2>&1
  [ $unset_after = 1 ] && unset OMP_NUM_THREADS
done
python3 - "$mb" "${OMP_NUM_THREADS:-1}" "$@" <<'PY'
import sys, os
mb, nt, tags = sys.argv[1], sys.argv[2], sys.argv[3:]
here = os.path.dirname(os.path.abspath(sys.argv[0])) if False else os.getcwd()
res = {}; names = []
for t in tags:
    for line in open('results/%s_%sMiB_%st.txt' % (t, mb, nt)):
        if line.startswith('R '):
            f = line.split()
            res[(t, f[1])] = (float(f[3]), float(f[5]))
            if f[1] not in names: names.append(f[1])
print('%-18s' % ('pool %s MiB, %s thr' % (mb, nt)) + ''.join('%11s' % t for t in tags) + '   max err')
for n in names:
    print('%-18s' % n + ''.join('%11.1f' % res[(t, n)][0] if (t, n) in res else '%11s' % '-' for t in tags)
          + '   %.1e' % max(res[(t, n)][1] for t in tags if (t, n) in res))
PY

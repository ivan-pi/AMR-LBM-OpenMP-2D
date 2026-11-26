#!/bin/bash
# Run built configurations and tabulate best MLUPS:  ./run.sh <tag> [<tag> ...]
[ -f ~/fortran-env.sh ] && source ~/fortran-env.sh; set +u   # compiler run-time libraries (this sandbox)
# (single thread pinned to a core unless OMP_NUM_THREADS is set). Results go to results/<tag>_<pool>.txt
here=$(cd "$(dirname "$0")" && pwd); mkdir -p $here/results
mb=mesh
for tag in "$@"; do
  if [ -z "$OMP_NUM_THREADS" ]; then pin="taskset -c 2"; export OMP_NUM_THREADS=1; unset_after=1; else pin=""; unset_after=0; fi
  out=$here/results/${tag}_${mb}_${OMP_NUM_THREADS}t.txt
  $pin $here/build/$tag/bench_ghost ${SWEEPS:-20} ${TRIALS:-11} > $out 2>&1
  [ $unset_after = 1 ] && unset OMP_NUM_THREADS
done
python3 - "$mb" "${OMP_NUM_THREADS:-1}" "$@" <<'PY'
import sys, os
mb, nt, tags = sys.argv[1], sys.argv[2], sys.argv[3:]
here = os.path.dirname(os.path.abspath(sys.argv[0])) if False else os.getcwd()
res = {}; names = []
for t in tags:
    for line in open('results/%s_%s_%st.txt' % (t, mb, nt)):
        if line.startswith('R '):
            f = line.split()
            res[(t, f[1])] = (float(f[2]), float(f[5]), float(f[4]))
            if f[1] not in names: names.append(f[1])
print('fill | collision+fill, MLUPS, %s thr' % nt)
print('%-18s' % '' + ''.join('%15s' % t for t in tags) + '   max err')
for n in names:
    print('%-18s' % n + ''.join('%8.1f |%5.1f' % (res[(t, n)][0], res[(t, n)][2]) if (t, n) in res else '%15s' % '-' for t in tags)
          + '   %.1e' % max(res[(t, n)][1] for t in tags if (t, n) in res))
PY

import os
import sys
import h5py
from numpy import *

for ep in ["MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
    os.environ[ep] = "1"
try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass

if not len(sys.argv) == 5:
    print('usage:: python pid.bb.xyz2h5.py N xxx.ex.xyz episoid_start episoid_end')
    sys.exit()
N  = int(sys.argv[1])
fx = str(sys.argv[2])
e0 = int(sys.argv[3])
e1 = int(sys.argv[4]) + 1
lb = fx[:fx.find('.e', fx.find('kd'))]

ts = []
for edx in range(e0, e1):
    fx = lb + ".e%d.xyz" % (edx)
    if not os.path.isfile(fx):
        print('Cannot find ' + fx)
        sys.exit()
    else:
        with open(fx) as f:
            for line in f:
                if not line[0] == '#':
                    lt = line.strip().split()
                    tv = float(lt[0])
                    if not tv in ts:
                        ts.append(tv)
nf = len(ts)
N2 = int(N*2)

fw = h5py.File(lb + '.eX.h5', 'w'); fw.create_dataset('T', data=array(ts), dtype='f')
xs = fw.create_dataset('X', (nf, N2, 3), dtype='f')
ss = fw.create_dataset('S', (N2,), dtype='i')
cs = fw.create_dataset('C', (N2,), dtype='i')

for edx in range(e0, e1):
    print("edx: %2d/%2d" % (edx, e1-e0))
    fx = lb + ".e%d.xyz" % (edx)
    if not os.path.isfile(fx):
        print('Cannot find ' + fx)
        sys.exit()
    else:
        adx = 0
        with open(fx) as f:
            for line in f:
                if not line[0] == '#':
                    lt = line.strip().split()
                    fdx = ts.index( float(lt[0]) )
                    for k in range(0, 3):
                        xs[fdx, adx, k] = float(lt[k+1])
                    if fdx == 0:
                        ss[adx] = int(lt[4])
                        cs[adx] = int(lt[5])
                    adx = (adx+1)%N2
fw.close()


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

def pij2ps(pij, N):
    ps = zeros(N)
    ps[0] = 1
    for s in range(1, N):
        for i in range(0, N-s):
            ps[s] += pij[i,i+s]
        ps[s] = ps[s]/float(N-s)
    return ps

if not len(sys.argv) == 6:
    print('usage:: python pid.xb.h5px.dist.py N xxx.h5 frameStart frameStep fq[=0/1:gaussian/step]')
    sys.exit()
N  =   int(sys.argv[1])
fx =   str(sys.argv[2])
f0 =   int(sys.argv[3])
fs =   int(sys.argv[4])
rc =   arange(1,5)*0.5
fq =   int(sys.argv[5])

lb = fx[:-3]
rcn= len(rc)
rc2= rc**2
N4 = int(N*4)

P = zeros((N4, N4, rcn))
nsample = 0

if True:
    if not os.path.isfile(fx):
        print('Cannot find '+ fx)
        sys.exit()
    else:
        fr = h5py.File(fx, 'r')
        nf = fr.get('X')[()].shape[0]

        for fdx in range(f0, nf, fs):
            X = fr.get('X')[fdx]
            for i in range(0, N4):
                vs = X - X[i]
                r2 = sum(vs**2, axis=1)

                for qdx in range(0, rcn):
                    if fq == 0:
                        P[i,:,qdx] += exp(-1.5*r2/rc2[qdx])
                    else:
                        P[i,:,qdx] += array(r2<rc2[qdx], dtype=float)
            nsample += 1
        fr.close()
P = P/float(nsample)

for qdx in range(0, rcn):
    fo = lb + ".rc%3.1f_fq%d" % (rc[qdx], fq)

    fw = open(fo+'.px', 'w')
    fw.write("#NSamples: %d\n" % (nsample))
    for i in range(0, N4):
        lt = ''
        for j in range(0, N4):
            lt += "%11.5e " % (P[i,j,qdx])
        fw.write(lt+'\n')
    fw.close()

    fw = open(fo+'.ps-cis', 'w')
    for c in range(0, 4):
        ia = c*N
        ib = ia + N
        ps = pij2ps( P[ia:ib, ia:ib, qdx], N)
        fw.write("#chain:%d\n#s p(s)\n"%(c))
        for s in range(0, N):
            lt = "%3d %11.5e " % (s, ps[s])
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()

    fw = open(fo+'.ps-mp', 'w')
    for c in range(0, 2):
        ia = c*N
        ib = int(ia + 2*N)
        fw.write("#chain:%d-%d\n#i p_{mp}(i)\n"%(c, c+2))
        for i in range(0, N):
            lt = "%3d %11.5e " % (i, P[ia+i, ib+i, qdx])
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()


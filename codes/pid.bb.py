import os
import sys
from numpy import *
try:
    from numba import jit
except ImportError:
    print('Please install Numba')
    sys.exit()

for ep in ["MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
    os.environ[ep] = "1"
try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass

def initialize(N, sdx):
    N2 = int(N*2)
    X = zeros((N2, 3))
    for i in range(1, N):
        v = random.normal(0,1,3)
        u = v/sqrt(sum(v**2))
        X[i] = X[i-1] + u
        X[i+N] = X[i-1+N] + u

    ss = []
    with open(fq) as f:
        for line in f:
            if not line[0] == '#':
                lt = list(line.strip())
                ss.append( list(map(int, lt)) )
    S = ss[sdx]
    S+= S
    S = array(S, dtype=int)

    C = zeros(N2, dtype=int)
    for c in range(0, 2):
        C[c*N:c*N+N] = c
    return X, S, C

def readCheckpoint(ckp):
    if not os.path.isfile(ckp):
        print('Cannot find '+ckp)
        sys.exit()
    else:
        X = []
        S = []
        C = []
        tz = nan
        with open(ckp) as f:
            for line in f:
                if not line[0] == '#':
                    lt = line.strip().split()
                    X.append( array(list(map(float, lt[1:4]))) )
                    S.append( int(lt[4]) )
                    C.append( int(lt[5]) )
                    if isnan(tz):
                        tz = float(lt[0])
        X = array(X)
        S = array(S, dtype=int)
        C = array(C, dtype=int)
        return X, S, C, tz

def readX(X, S, C, head):
    lt = ''
    for i in range(0, len(X)):
        lt += head
        for k in range(0, 3):
            lt += "%+12.5e " % (X[i][k])
        lt += "%d " % (S[i])
        lt += "%d " % (C[i])
        lt += '\n'
    return lt

def X2ree(X, N):
    ree = zeros(2)
    for c in range(0, 2):
        idx = c*N
        jdx = idx + N - 1
        v = X[idx] - X[jdx]
        ree[c] = sqrt(sum(v**2))
    return ree

@jit(nopython=True)
def ppc_brownianDynamics(N, N2, X, S, kbond, Dt, Dsq, NoS, rsd):
    random.seed(rsd)
    step = 0
    while step < NoS:
        F = zeros((N2, 3)) # Deterministic force
        for c in range(0, 2):
            for i in range(1, N):
                idx = c*N + i
                v = X[idx] - X[idx-1]
                r = sqrt(sum(v**2))
                F[idx]   -= kbond*(r - 1.0)*v/r
                F[idx-1] += kbond*(r - 1.0)*v/r
        for i in range(0, N):
            if S[i] == 1:
                F[i] += F[i+N] # The buttom sites on the 2nd chain are virtual

        R = random.normal(0, 1, (N2, 3)) # Stochastic force
        for i in range(0, N):
            X[i] += Dt[i]*F[i] + Dsq[i]*R[i]
        for i in range(N, N2):
            if S[i] == 1:
                X[i] = X[i-N]
            else:
                X[i] += Dt[i]*F[i] + Dsq[i]*R[i]
        step += 1
    return X

# Brownian dynamics simulation of bubble-chain model
#
if not len(sys.argv) == 7:
    print('usage:: python pid.bb.py xxx.seq seq-idx[0/1] k_{bond} D_{pair}[0.7937] replica-index episoide')
    sys.exit()
fq =    str(sys.argv[1])
sdx =   int(sys.argv[2])
kbond=float(sys.argv[3]) # k_{0} (Eq. 3)
dpair=float(sys.argv[4]) # D_{pair} = D*d_{pair}
rdx  =  int(sys.argv[5])
edx  =  int(sys.argv[6])

lb = "bb.%s.s%d_kd%04.1f_dp%6.4f.R%d.e%d" % (fq[:-4], sdx, kbond, dpair, rdx, edx)

if not os.path.isfile(fq):
    print('Cannot find ' + fq)
    sys.exit()
else:
    with open(fq) as f:
        for line in f:
            if not line[0] == '#':
                N  = len(line.strip()) # Number of beads per chain
                break
    N2 = int(2*N) # Total number of beads

fx = 'mdRandomSeeds'
if os.path.isfile(fx):
    ldx = 0
    with open(fx) as f:
        for line in f:
            lt = line.strip()
            if ldx == rdx:
                randomSeed = int(lt)
                break
            ldx += 1
    random.seed(randomSeed)
else:
    print('Cannot find '+fx)
    sys.exit()

# Main entrance
if True:
    if edx == 0:
        X, S, C = initialize(N, sdx)
        tz = 0
    else:
        ckp = "bb.%s.s%d_kd%04.1f_dp%6.4f.R%d.e%d.ckp" % (fq[:-4], sdx, kbond, dpair, rdx, edx-1)
        X, S, C, tz = readCheckpoint(ckp)

    fw = open(lb+'.xyz', 'w')
    if edx == 0:
        fw.write("#time x y z species\n"); fw.write( readX(X, S, C, "%11.5e "%(0)) )
    fw.close()
    fw = open(lb+'.log', 'w')
    if edx == 0:
        ree = X2ree(X, N)
        fw.write("#time ree(0) ree(1)\n")
        log  = "%11.5e "%(0)
        for k in range(0, 2):
            log += "%11.5e " % (ree[k])
        fw.write(log+'\n')
    fw.close()

    D = ones(N2)
    for i in range(0, N2):
        if S[i] == 1:
            D[i] *= dpair
    dt = 0.01 # 0.01
    numberOfSteps = int(2e+0*(N**2)/dt) # 1e+2
    numberOfBlocks = 2000
    numberOfStepsPerBlock = int(numberOfSteps/numberOfBlocks)
    timePerBlock = dt*numberOfStepsPerBlock
    Dt = D*dt
    Dsq = sqrt(2.0*D*dt)

    indexOfBlock = 0
    while indexOfBlock < numberOfBlocks:
        rsd = randomSeed + numberOfBlocks*edx + indexOfBlock
        X = ppc_brownianDynamics(N, N2, X, S, kbond, Dt, Dsq, numberOfStepsPerBlock, rsd)
        indexOfBlock += 1

        st = indexOfBlock*timePerBlock
        fw = open(lb+'.xyz', 'a')
        fw.write( readX(X, S, C, "%11.5e "%(st+tz)) )
        fw.close()

        fw = open(lb+'.log', 'a')
        ree = X2ree(X, N)
        log = "%11.5e "%(st+tz)
        for k in range(0, 2):
            log += "%11.5e " % (ree[k])
        fw.write(log+'\n')
        fw.close()

        if indexOfBlock == numberOfBlocks:
            fw = open(lb+'.ckp', 'w')
            fw.write("#time x y z species c\n")
            fw.write( readX(X, S, C, "%11.5e "%(st+tz)) )
            fw.close()
        print("%5d / %5d "%(indexOfBlock, numberOfBlocks))


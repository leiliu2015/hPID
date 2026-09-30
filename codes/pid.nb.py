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

def initialize(N):
    X = zeros((4*N, 3))
    for c in range(0, 4):
        for i in range(1, N):
            idx = c*N + i
            v = random.normal(0,1,3)
            u = v/sqrt(sum(v**2))
            X[idx] = X[idx-1] + u

    S = []
    with open(fq) as f:
        for line in f:
            if not line[0] == '#':
                lt = list(line.strip())
                S += list(map(int, lt))
    S+= S
    S = array(S, dtype=int)

    C = zeros(4*N, dtype=int)
    for c in range(0, 4):
        C[c*N:c*N+N] = c
    PL = zeros(4*N, dtype=int) - 1 # Pairing list
    return X, S, C, PL

def readCheckpoint(ckp):
    if not os.path.isfile(ckp):
        print('Cannot find '+ckp)
        sys.exit()
    else:
        X = []
        S = []
        C = []
        PL = []
        tz = nan
        with open(ckp) as f:
            for line in f:
                if not line[0] == '#':
                    lt = line.strip().split()
                    X.append( array(list(map(float, lt[1:4]))) )
                    S.append( int(lt[4]) )
                    C.append( int(lt[5]) )
                    PL.append(int(lt[6]) )
                    if isnan(tz):
                        tz = float(lt[0])
        X = array(X)
        S = array(S, dtype=int)
        C = array(C, dtype=int)
        PL= array(PL,dtype=int)
        return X, S, C, PL, tz

def readX(X, S, C, PL, head):
    lt = ''
    for i in range(0, len(X)):
        lt += head
        for k in range(0, 3):
            lt += "%+12.5e " % (X[i][k])
        lt += "%d %d " % (S[i], C[i])
        lt += "%4d " % (PL[i])
        lt += '\n'
    return lt

def PL2fhp(PL, N, nob):
    N2 = int(2*N)
    hp = zeros((N, 2))
    for c in range(0, 2):
        for i in range(0, N):
            idx = c*N + i
            jdx = idx + N2
            if PL[idx] == jdx:
                hp[i][c] = 1.0
    return sum(hp, axis=0)/nob

def X2ree(X, N):
    ree = zeros(4)
    for c in range(0, 4):
        idx = c*N
        jdx = idx + N - 1
        v = X[idx] - X[jdx]
        ree[c] = sqrt(sum(v**2))
    return ree

@jit(nopython=True)
def wm_brownianDynamics(N, N4, X, S, C, rt, rc, kbond, kconf, kbend, PL, pu, PQ, D, dt, NoS, rsd):
    random.seed(rsd)
    Dt = D*dt
    Dsq = sqrt(2.0*D*dt)
    T2 = rt**2
    rcap2 = rc**2

    t2 = 2.0**(-2.0/3)
    t1 = 2.0**(-1.0/3)

    step = 0
    while step < NoS:

        if PQ and (pu>0):
            for i in range(0, N4):
                j = PL[i]
                if (not j == -1) and (random.random() < pu):
                    PL[i] = -1
                    PL[j] = -1

        F = zeros((N4, 3)) # Deterministic force
        for c in range(0, 4):
            idx = c*N
            v = X[idx]
            r2 = sum(v**2)
            if r2 > T2:
                r = sqrt(r2)
                F[idx] -= kconf*(r - rt)*v/r

            for i in range(1, N):
                idx = c*N + i
                v = X[idx] - X[idx-1]
                r = sqrt(sum(v**2))
                F[idx]   -= kbond*(r - 1.0)*v/r
                F[idx-1] += kbond*(r - 1.0)*v/r

            if kbend > 0:
                for i in range(1, N-1):
                    idx = c*N + i
                    v12 = X[idx] - X[idx-1]
                    v23 = X[idx+1] - X[idx]
                    r12 = sqrt(sum(v12**2))
                    r23 = sqrt(sum(v23**2))
                    dpt = sum(v12*v23)
                    f1 = kbend * (-v23/r12 + v12/(r12**3)*dpt) / r23
                    f3 = kbend * ( v12/r23 - v23/(r23**3)*dpt) / r12
                    F[idx-1] += f1
                    F[idx+1] += f3
                    F[idx] -= (f1+f3)

        if PQ:
            R = zeros((N4, 3)) # Random force
            for i in range(0, N4):
                if PL[i] == -1:
                    R[i] = random.normal(0, 1.0, 3)
                else:
                    R[i] = random.normal(0, t2, 3)
            for i in range(0, N4):
                if PL[i] == -1:
                    X[i] += Dt*F[i] + Dsq*R[i]
                else:
                    j = PL[i]
                    X[i] += Dt*F[i]*t1 + Dsq*R[i]
                    X[i] += Dt*F[j]*t1 + Dsq*R[j]
        else:
            R = random.normal(0, 1, (N4, 3))
            for i in range(0, N4):
                X[i] += Dt*F[i] + Dsq*R[i]

        if PQ:
            for i in range(0, N4):
                for j in range(i+1, N4):
                    if (S[i] == 1) and (S[j] == 1) and (PL[i] == -1) and (PL[j] == -1):
                        v = X[j] -  X[i]
                        r2 = sum(v**2)
                        if r2 < rcap2:
                            PL[i] = j
                            PL[j] = i
        step += 1
    return X, PL

# Brownian dynamics simulation of nonspecific-button model
#
if not len(sys.argv) == 8:
    print('usage:: python pid.nb.py xxx.seq k_{bond} k_{bend} p_{unpair} R_{tether} replica-index episoide')
    sys.exit()
fq =    str(sys.argv[1])
kbond=float(sys.argv[2]) # k_{0} (Eq. 3)
kb =  float(sys.argv[3]) # k_{b} (Eq. 4)
pu =  float(sys.argv[4]) # p_{u}
kconf = 35.              # k_{c} (Eq. 5)
rt =  float(sys.argv[5]) # R (Eq. 5)
rdx  =  int(sys.argv[6])
edx  =  int(sys.argv[7])

rc = 0.5 # Radius of capture
lb = "nb.%s.kd%04.1f_kb%3.1f_pu%5.0e_rt%3.1f.R%d.e%d" % (fq[:-4], kbond, kb, pu, rt, rdx, edx)

if not os.path.isfile(fq):
    print('Cannot find ' + fq)
    sys.exit()
else:
    with open(fq) as f:
        for line in f:
            if not line[0] == '#':
                N  = len(line.strip()) # Number of beads per chain
                break
    L  = int(N-1) # Number of links per chain
    N4 = int(4*N) # Total number of beads

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
        X, S, C, PL = initialize(N)
        tz = 0
        nob = array([sum(S[0:N]), sum(S[N:2*N])]) # number of buttons per chain
    else:
        ckp = "nb.%s.kd%04.1f_kb%3.1f_pu%5.0e_rt%3.1f.R%d.e%d.ckp" % (fq[:-4], kbond, kb, pu, rt, rdx, edx-1)
        X, S, C, PL, tz = readCheckpoint(ckp)
        nob = array([sum(S[0:N]), sum(S[N:2*N])])

    fw = open(lb+'.xyz', 'w')
    if edx == 0:
        fw.write("#time x y z species c pairing\n"); fw.write( readX(X, S, C, PL, "%11.5e "%(0)) )
    fw.close()
    fw = open(lb+'.log', 'w')
    if edx == 0:
        fhp = PL2fhp(PL, N, nob)
        ree = X2ree(X, N)

        fw.write("#time ree(0) ree(1) ree(2) ree(3) fhp(0-2) fhp(1-3)\n")
        log = "%11.5e "%(0)
        for k in range(0, 4):
            log += "%11.5e " % (ree[k])
        for k in range(0, 2):
            log += "%11.5e " % (fhp[k])
        fw.write(log+'\n')
    fw.close()

    D = 1.0
    dt = 0.01 # 0.01
    numberOfSteps = int(2e+0*(N**2)/dt); # 1e+2
    numberOfBlocks = 2000
    numberOfStepsPerBlock = int(numberOfSteps/numberOfBlocks)
    timePerBlock = dt*numberOfStepsPerBlock
    if edx == 0:
        PQ = False
    else:
        PQ = True

    indexOfBlock = 0
    while indexOfBlock < numberOfBlocks:
        rsd = randomSeed + numberOfBlocks*edx + indexOfBlock
        X, PL = wm_brownianDynamics(N, N4, X, S, C, rt, rc, kbond, kconf, kb, PL, pu, PQ, D, dt, numberOfStepsPerBlock, rsd)
        indexOfBlock += 1

        st = indexOfBlock*timePerBlock
        fw = open(lb+'.xyz', 'a')
        fw.write( readX(X, S, C, PL, "%11.5e "%(st+tz)) )
        fw.close()

        fw = open(lb+'.log', 'a')
        fhp = PL2fhp(PL, N, nob)
        ree = X2ree(X, N)
        log = "%11.5e "%(st+tz)
        for k in range(0, 4):
            log += "%11.5e " % (ree[k])
        for k in range(0, 2):
            log += "%11.5e " % (fhp[k])
        fw.write(log+'\n')
        fw.close()

        if indexOfBlock == numberOfBlocks:
            fw = open(lb+'.ckp', 'w')
            fw.write("#time x y z species c pairing\n")
            fw.write( readX(X, S, C, PL, "%11.5e "%(st+tz)) )
            fw.close()
        print("%5d / %5d "%(indexOfBlock, numberOfBlocks))


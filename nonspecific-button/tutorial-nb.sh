#!/bin/bash

N=90
fq=example30
cp ../codes/${fq}.seq ./
cp ../codes/mdRandomSeeds ./

kb=0.0
numberOfReplicas=1

# Simulations
for ((rdx=0; "$rdx" < ${numberOfReplicas}; rdx++))
do
    for ((edx=0; "$edx" < 2; edx++))
    do
        python ../codes/pid.nb.py ${fq}.seq 25.0 ${kb} 0.001 2.0 ${rdx} ${edx}
    done
done

# Analyses
python ../codes/pid.xb.xyz2h5.py ${N} nb.${fq}.kd25.0_kb${kb}_pu1e-03_rt2.0.R0.e1.xyz 1 1
python ../codes/pid.xb.h5px.dist.py ${N} nb.${fq}.kd25.0_kb${kb}_pu1e-03_rt2.0.R0.eX.h5 1000 1 0




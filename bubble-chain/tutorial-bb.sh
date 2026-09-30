#!/bin/bash

N=90
fq=example30
cp ../codes/${fq}.seq ./
cp ../codes/mdRandomSeeds ./

Dpair=0.7937
numberOfReplicas=1

# sdx=0/1: M1-P1/M2-P2
for ((sdx=0; "$sdx" < 2; sdx++))
do
    # Simulations
    for ((rdx=0; "$rdx" < ${numberOfReplicas}; rdx++))
    do
        for ((edx=0; "$edx" < 1; edx++))
        do
            python ../codes/pid.bb.py ${fq}.seq ${sdx} 25.0 ${Dpair} ${rdx} ${edx}
        done
    done

    # Analyses
    python ../codes/pid.bb.xyz2h5.py ${N} bb.example30.s${sdx}_kd25.0_dp${Dpair}.R0.e0.xyz 0 0
    python ../codes/pid.bb.h5px.dist.py ${N} bb.example30.s${sdx}_kd25.0_dp${Dpair}.R0.eX.h5 1000 1 0
done


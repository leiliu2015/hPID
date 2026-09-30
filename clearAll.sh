#/bin/bash

fds=(nonspecific-button specific-button bubble-chain)

for ((fdx=0; "$fdx" < 3; fdx++))
do
    fd=${fds[$fdx]}
    echo ${fd}

    cd ./${fd}
    if [ -f mdRandomSeeds ]
    then
        rm mdRandomSeeds
    fi

    if [ -f example30.seq ]
    then
        rm example30.seq
    fi

    if stat *.example* >/dev/null 2>&1
    then
        rm *.example*
    fi
    cd ../
done


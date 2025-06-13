#!/bin/bash

# this is an example of a parallel run setup, suitable for running on a single
# multi-core node. It is meant to illustrate the sequence of the typical parallel run.
# Appropriate scripts should be setup for use with batch systems, depending upon
# the specific implementation.


> Timings.txt

# number of cores to be used for a given parallel run
ncores=8

# number of parallel runs; better be a multiple of ncores.
nprocesses=2

# number of iterations for the calculation of xgrid at parallelstage 1. It is the old itmx1
nxgriditeration=2

# The pwhg_main executable should be in the ../ directory
prg=../pwhg_main


# Only do a single stage
STAGE="0"
while [[ $# -gt 0 ]]; do
KEY="$1"
case $KEY in
    --st1)
        STAGE="1"
        shift
        ;;
     --st2)
        STAGE="2"
        shift
        ;;
     --st3)
        STAGE="3"
        shift
        ;;
     --st4)
        STAGE="4"
        shift
        ;;
esac
done

# the following function limits the number of subprocesses
# to be not larger than the number of cores specified by the
# user
function limit_procs {
    while [ `jobs -p | wc -w` -ge $ncores ]
    do
	sleep 1
    done
}


#-----------------------------------------------------------------
# PARALLEL STAGE 1
#-----------------------------------------------------------------
# Compute the x-grids for integration
if [ $STAGE == "0" ] || [ $STAGE == "1" ]
then
parstage=1
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

# nxgriditeration stages of importance sampling grid calculation
for igrid in $(seq 1 $nxgriditeration)	     
do
    echo "     xgriditeration " $igrid
    (echo -n st$parstage xg$igrid ' ' ; date ) >> Timings.txt
    cat powheg.input-save | sed "s/xgriditeration.*/xgriditeration $igrid/ ; s/parallelstage.*/parallelstage $parstage/" > powheg.input
    cp powheg.input powheg.input-$parstage-$igrid
    
    for i in `seq $nprocesses`
    do
	echo $i | $prg > run-st1-xg$igrid-$i.log 2>&1 &
	limit_procs
    done
    wait    
done
(echo -n end st$parstage ' ' ; date ) >> Timings.txt
fi

#-----------------------------------------------------------------
# PARALLEL STAGE 2
#-----------------------------------------------------------------
# compute NLO and upper bounding envelope for underlying Born configurations

if [ $STAGE == "0" ] || [ $STAGE == "2" ]
then
parstage=2
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/" > powheg.input
cp powheg.input powheg.input-$parstage
(echo -n beg st$parstage ' ' ; date ) >> Timings.txt
for i in `seq $nprocesses`
do
    echo $i | $prg > run-st2-$i.log 2>&1 &
    limit_procs
done
# this stage should be fully completed before the next stage is started.
# on batch systems, failed jobs can be abandoned, if their number is not too large.
wait
(echo -n end st$parstage ' ' ; date ) >> Timings.txt
fi
#-----------------------------------------------------------------
# PARALLEL STAGE 3
#-----------------------------------------------------------------
# compute upper bounding coefficients for radiation


if [ $STAGE == "0" ] || [ $STAGE == "3" ]
then
parstage=3
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/" > powheg.input
cp powheg.input powheg.input-$parstage
(echo -n beg st$parstage ' ' ; date ) >> Timings.txt
for i in `seq $nprocesses`
do
    echo $i | $prg > run-st3-$i.log 2>&1 &
    limit_procs
done
# this stage should be fully completed before the next stage is started.
# on batch systems, failed jobs can be abandoned, if their number is not too large.
wait
(echo -n end st$parstage ' ' ; date ) >> Timings.txt
fi


if [ $STAGE == "0" ] || [ $STAGE == "4" ]
then
parstage=4
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/" > powheg.input
cp powheg.input powheg.input-$parstage
(echo -n beg st$parstage ' ' ; date ) >> Timings.txt
for i in `seq $nprocesses`
do
    echo $i | $prg > run-st4-$i.log 2>&1 &
    limit_procs
done
# this stage should be fully completed before the next stage is started.
# on batch systems, failed jobs can be abandoned, if their number is not too large.
wait
(echo -n end st$parstage ' ' ; date ) >> Timings.txt
fi

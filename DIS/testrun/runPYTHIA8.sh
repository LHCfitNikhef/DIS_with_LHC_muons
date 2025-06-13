#!/bin/bash

#specify path to current directory so that all files are found when submitting job to cluster
source $PBS_O_HOME/.bashrc


> Timings.txt

ncores=32
nprocesses=32
# First compile the pwhg_main executable in the ../ directory
#

# the following function limits the number of subprocesses
# to be not larger than the number of cores specified by the
# user

function limit_procs {
    while [ `jobs -p | wc -w` -gt $ncores ]
    do
	sleep 1
    done
}

    
outfile=pwgPOWHEG+PYTHIA8-output
    
for file in pwgevents-0*.lhe
do
    which=PYTHIA8-`echo $file | sed 's/pwgevents-//; s/lhe/log/ '`
    echo $file | ../main-PYTHIA8-lhef > $which &
done

wait

labels=( 'dummy1' '11' '12' '21' '22' '1H' 'H1' 'HH' )
# adjust number of weights accordingly
for i in {2..8}
do
    mergedata 1 pwgPOWHEG+PYTHIA8-output-00*W$i.top -o $outfile-${labels[$i-1]}.top
done

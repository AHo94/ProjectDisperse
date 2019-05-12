#!/bin/bash
#echo "Spawning worker processes"
#Nproc=50
Nproc=$1
for ((i=0; i<= Nproc; i++))
do
    (nice -n13 python SpeedProfileWorker.py &)
done

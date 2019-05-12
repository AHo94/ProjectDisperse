#!/bin/bash
#echo "Spawning worker processes"
#Nproc=50
Nproc=$1
RunEuclid=$2
for ((i=0; i<= Nproc; i++))
do
    ( nice -n13 python ParticlesPerFilament.py -euclid21 $((RunEuclid)) &)
done

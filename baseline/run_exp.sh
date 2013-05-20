#!/bin/bash

meas=./main

lpgen_bench=../lpbench

rm *.dat;

for i in {1..10}; do
    for j in 0008 0016 0032 0064 0128 0256 0512 1024; do
        $meas $lpgen_bench/lpb_${j}_$i.dlp $lpgen_bench/lpb_${j}_$i.lp >> $j.dat &
    done
done 

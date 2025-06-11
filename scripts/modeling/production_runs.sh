#!/bin/bash

for i in `seq 15`; do
    folder_name="$i"
    mpirun -np 8 ~/imp-clean/build/setup_environment.sh python ../scripts/modeling/mic10_dimer_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" ../inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;


for i in `seq 16 30`; do
    folder_name="$i"
    mpirun -np 8 ~/imp-clean/build/setup_environment.sh python ../scripts/modeling/mic10_dimer_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" ../inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;

for i in `seq 31 45`; do
    folder_name="$i"
    mpirun -np 8 ~/imp-clean/build/setup_environment.sh python ../scripts/modeling/mic10_dimer_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" ../inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;

for i in `seq 46 50`; do
    folder_name="$i"
    mpirun -np 8 ~/imp-clean/build/setup_environment.sh python ../scripts/modeling/mic10_dimer_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" ../inputs/data 2> err_$folder_name.log &
done
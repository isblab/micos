#!/bin/bash

for i in `seq 10`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/projects/micos/github/scripts/modeling/mic10_mic13_rb_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" /home/muskaan/projects/micos/github/inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;


for i in `seq 11 21`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/projects/micos/github/scripts/modeling/mic10_mic13_rb_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" /home/muskaan/projects/micos/github/inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;

for i in `seq 22 32`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/projects/micos/github/scripts/modeling/mic10_mic13_rb_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" /home/muskaan/projects/micos/github/inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;

for i in `seq 33 43`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/projects/micos/github/scripts/modeling/mic10_mic13_rb_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" /home/muskaan/projects/micos/github/inputs/data 2> err_$folder_name.log &
done

wait
sleep 600 ;

for i in `seq 44 50`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/projects/micos/github/scripts/modeling/mic10_mic13_rb_mic19_1full_1N_1C_independent_fixed.py "prod" "$folder_name" /home/muskaan/projects/micos/github/inputs/data 2> err_$folder_name.log &
done

exit
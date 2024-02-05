#!/bin/bash

for i in `seq 24`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/Documents/modeling_micos_complex/micos/scripts/modeling/mic19_1full_1N_1C_fixed_apdimer.py "prod" "$folder_name" /home/muskaan/Documents/modeling_micos_complex/micos/inputs/data 2> err_$folder_name.log &
done
sleep 600 ;

mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/Documents/modeling_micos_complex/micos/scripts/modeling/mic19_1full_1N_1C_fixed_apdimer.py "prod" "25" /home/muskaan/Documents/modeling_micos_complex/micos/inputs/data 2> "err_25.log" &

for i in `seq 26 50`; do
    folder_name="$i"
    mpirun -np 8 /home/muskaan/imp-clean/build/setup_environment.sh python /home/muskaan/Documents/modeling_micos_complex/micos/scripts/modeling/mic19_1full_1N_1C_fixed_apdimer.py "prod" "$folder_name" /home/muskaan/Documents/modeling_micos_complex/micos/inputs/data 2> err_$folder_name.log &
done

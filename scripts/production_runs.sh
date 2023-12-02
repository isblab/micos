#!/bin/bash

for i in `seq 24`; do
    folder_name="$i"
    ~/imp-clean/build/setup_environment.sh mpirun -np 8 python ~muskaan/Documents/modeling_micos_complex/micos/scripts/mic19_1full_1N_1C.py "prod" "$folder_name" ~muskaan/Documents/modeling_micos_complex/micos/inputs/data 2> err_$folder_name.log &
done
sleep 600 ;
~/imp-clean/build/setup_environment.sh mpirun -np 8 python ~muskaan/Documents/modeling_micos_complex/micos/scripts/mic19_1full_1N_1C.py "prod" "25" ~muskaan/Documents/modeling_micos_complex/micos/inputs/data 2> "err_output_25.log" &

for i in `seq 26 50`; do
    folder_name="$i"
    ~/imp-clean/build/setup_environment.sh mpirun -np 8 python ~muskaan/Documents/modeling_micos_complex/micos/scripts/mic19_1full_1N_1C.py "prod" "$folder_name" ~muskaan/Documents/modeling_micos_complex/micos/inputs/data 2> err_$folder_name.log &
done

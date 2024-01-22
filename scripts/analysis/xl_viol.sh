#!/bin/bash
sampling
input_files=("BDP_PIR_mouse.csv" "BDP_PIR_human.csv" "DHSO_DSSO_human.csv" "DSS_BS3_human.csv" "DSSO_mouse.csv")
parameters=("52" "52" "30" "35" "30")

for ((i=0; i<${#input_files[@]}; i++))
do
    input_file="../../micos/inputs/data/crosslinks/human/sampling_${input_files[$i]}"
    parameter=${parameters[$i]}

    echo "$input_file with parameter $parameter"
    $imp ~/imp-clean/build/setup_environment.sh python ~/IMP_Toolbox/analysis/fit_data/get_xlviol_val_set_v2.py sampcon_0_extracted.rmf3 "$input_file" "$parameter" &
done
#
# #evaluation
input_files=("BDP_PIR" "DHSO_DSSO" "DSS_BS3")
parameters=("52" "30" "35")

for ((i=0; i<${#input_files[@]}; i++))
do
    input_file="../../micos/inputs/data/crosslinks/human/evicalc_${input_files[$i]}"
    parameter=${parameters[$i]}

    echo "$input_file with parameter $parameter"
    $imp ~/imp-clean/build/setup_environment.sh python ~/IMP_Toolbox/analysis/fit_data/get_xlviol_val_set_v2.py sampcon_0_extracted.rmf3 "$input_file" "$parameter" &
done

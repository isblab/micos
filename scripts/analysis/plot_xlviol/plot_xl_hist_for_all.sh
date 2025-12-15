# !/bin/bash
# sampling
input_files=("BDP_PIR_mouse.csv.csv" "BDP_PIR_human.csv.csv" "DHSO_DSSO_human.csv.csv" "DSS_BS3_human.csv.csv" "DSSO_mouse.csv.csv")
parameters=("52" "52" "30" "35" "30")

for ((i=0; i<${#input_files[@]}; i++))
do
    parameter="${parameters[$i]}"

    python ../github/scripts/analysis/plot_xlviol/xl_distance_hist_plotter.py \
        "xl_viol/sampling/xl_violation_sampling_${input_files[$i]}" \
        "${input_files[$i]}" \
        "${parameter}"
done


mkdir -p sampling
mv *.png sampling

# evaluation
input_files=("BDP_PIR_mouse.csv.csv" "DHSO_DSSO_human.csv.csv" "DHSO_DSSO_mouse.csv.csv" "DSS_BS3_human.csv.csv")
parameters=("52" "30" "30" "35")

for ((i=0; i<${#input_files[@]}; i++))
do
        parameter="${parameters[$i]}"

    python ../github/scripts/analysis/plot_xlviol/xl_distance_hist_plotter.py \
        "xl_viol/evicalc/xl_violation_evicalc_${input_files[$i]}" \
        "${input_files[$i]}" \
        "${parameter}"
done

mkdir -p evaluation
mv *.png evaluation

rm -rf hist
mkdir hist
mv sampling evaluation hist
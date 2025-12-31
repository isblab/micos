#! /bin/bash

# python ../github/scripts/analysis/plot_xlviol/get_color_by_min_distance.py xl_viol/sampling/xl_violation_sampling_BDP_PIR_mouse.csv.csv sampling_BDP_PIR_mouse.csv
# python ../github/scripts/analysis/plot_xlviol/plot_circos.py sampling_BDP_PIR_mouse.csv

files=("BDP_PIR_mouse.csv.csv" "BDP_PIR_human.csv.csv" "DHSO_DSSO_human.csv.csv" "DSS_BS3_human.csv.csv" "DSSO_mouse.csv.csv") 
parameters=("52" "52" "30" "35" "30")

for ((i=0; i<${#files[@]}; i++))
do
    parameter="${parameters[$i]}"
    python ~/projects/micos/github/scripts/analysis/plot_xlviol/get_color_by_min_distance.py "xl_viol/sampling/xl_violation_sampling_${files[$i]}" sampling_${files[$i]} "${parameter}"
    python ~/projects/micos/github/scripts/analysis/plot_xlviol/plot_circos.py sampling_${files[$i]}

done

mkdir -p sampling
mv *.png *.csv sampling



files=("BDP_PIR_mouse.csv.csv" "DHSO_DSSO_human.csv.csv" "DHSO_DSSO_mouse.csv.csv" "DSS_BS3_human.csv.csv") 
parameters=("52" "30" "30" "35")

for ((i=0; i<${#files[@]}; i++))
do
    parameter="${parameters[$i]}"
    python ~/projects/micos/github/scripts/analysis/plot_xlviol/get_color_by_min_distance.py "xl_viol/evicalc/xl_violation_evicalc_${files[$i]}" evicalc_${files[$i]} "${parameter}"
    python ~/projects/micos/github/scripts/analysis/plot_xlviol/plot_circos.py evicalc_${files[$i]}

done

mkdir -p evaluation
mv *.png *.csv evaluation

rm -rf circos
mkdir circos
mv sampling evaluation circos
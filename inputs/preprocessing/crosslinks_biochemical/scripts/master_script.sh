#!/bin/bash

########## Map all data to human and yeast species  ###################
input_dir="../mapping_to_homologs/inputs/"
output_dir="../mapping_to_homologs/outputs/"
input_dir_biochem="../biochemical/"

##########----Step 1----###################

########## MAP all crosslinks including those from XlinkDB and other papers  ###################

########## mouse to human ###################
for input_file in "${input_dir}"mouse_PIR_xl.csv "${input_dir}"mouse_BDP_xl.csv "${input_dir}"mouse_DSSO_xl.csv
do
  base_filename="$(basename "$input_file")"
  output_file="${output_dir}${base_filename%.*}_XLinkDB_mouse_to_human"
  python 4_map_crosslinks_biochemical_to_homologs.py "$input_file" "$output_file" "MOUSE" "HUMAN"
done

########### yeast to human ##############

for data_ in "${input_dir}"yeast_DSSO_xl.csv
do
  base_filename="$(basename "$data_")"
  output_file="${output_dir}${base_filename%.*}_XLinkDB_yeast_to_human"
  python 4_map_crosslinks_biochemical_to_homologs.py "$data_" "$output_file" "YEAST" "HUMAN"
done

############ yeast to human crosslink papers ###############
for data_ in "${input_dir}"linden_BS3.csv
do
  base_filename="$(basename "$data_")"
  output_file="${output_dir}${base_filename%.*}_yeast_to_human"
  python 4_map_crosslinks_biochemical_to_homologs.py "$data_" "$output_file" "YEAST" "HUMAN"
done

########## mouse to yeast ################

for data in "${input_dir}"mouse_PIR_xl.csv "${input_dir}"mouse_BDP_xl.csv "${input_dir}"mouse_DSSO_xl.csv
do
  base_filename="$(basename "$data")"
  output_file="${output_dir}${base_filename%.*}_XLinkDB_mouse_to_yeast"
  python 4_map_crosslinks_biochemical_to_homologs.py "$data" "$output_file" "MOUSE" "YEAST"
done

############ human to yeast #############

for files in "${input_dir}"human_PIR_xl.csv "${input_dir}"human_BDP_xl.csv "${input_dir}"human_DSSO_xl.csv
do
  base_filename="$(basename "$files")"
  output_file="${output_dir}${base_filename%.*}_XLinkDB_human_to_yeast"
  python 4_map_crosslinks_biochemical_to_homologs.py "$files" "$output_file" "HUMAN" "YEAST"
done

############ human to yeast crosslink papers ###############
for files in "${input_dir}"ryl_BS3.csv "${input_dir}"yu_DSSO.csv "${input_dir}"sun_DSS.csv "${input_dir}"DSSO_Bartolec.csv "${input_dir}"DHSO_Bartolec.csv
do
  base_filename="$(basename "$files")"
  output_file="${output_dir}${base_filename%.*}_human_to_yeast"
  python 4_map_crosslinks_biochemical_to_homologs.py "$files" "$output_file" "HUMAN" "YEAST"
done

############ human to human to map mic25 to mic19 (only DSSO_Bartolec) ###############
for files in "${input_dir}"DSSO_Bartolec.csv
do
  base_filename="$(basename "$files")"
  output_file="${output_dir}${base_filename%.*}_human_to_human"
  python 4_map_crosslinks_biochemical_to_homologs.py "$files" "$output_file" "HUMAN" "HUMAN"
done

##########----Step 2----###################

########## Merge the crosslinks for human only ##################

output_dir_for_xlinks="../../../data/crosslinks/human/"

python 5_merge_files.py "${output_dir_for_xlinks}BDP_PIR_human.csv" "${input_dir}"human_BDP_xl.csv  "${input_dir}"human_PIR_xl.csv
python 5_merge_files.py "${output_dir_for_xlinks}BDP_PIR_mouse.csv" "${output_dir}"mouse_PIR_xl_XLinkDB_mouse_to_human  "${output_dir}"mouse_BDP_xl_XLinkDB_mouse_to_human
python 5_merge_files.py "${output_dir_for_xlinks}DHSO_DSSO.csv" "${input_dir}"human_DSSO_xl.csv  "${output_dir}"DSSO_Bartolec_human_to_human "${input_dir}"DHSO_Bartolec.csv "${input_dir}"yu_DSSO.csv
python 5_merge_files.py "${output_dir_for_xlinks}DSS_BS3.csv" "${input_dir}"ryl_BS3.csv  "${input_dir}"sun_DSS.csv
python 5_merge_files.py "${output_dir_for_xlinks}DSSO_XLinkDB_mouse.csv" "${output_dir}"mouse_DSSO_xl_XLinkDB_mouse_to_human

##########----Step 3----###################

########## Map biochemical data to human and yeast specieswise ###################

########### yeast to human ##############

for data_ in "${input_dir_biochem}"pairwise_biochemical_data.csv
do
  base_filename="$(basename "$data_")"
  output_file="${output_dir}${base_filename%.*}_biochemical_data_yeast_to_human"
  python 4_map_crosslinks_biochemical_to_homologs.py "$data_" "$output_file" "YEAST" "HUMAN"
done

############ human to yeast #############
for files in "${input_dir_biochem}"pairwise_biochemical_data.csv
do
  base_filename="$(basename "$files")"
  output_file="${output_dir}${base_filename%.*}_biochemical_data_human_to_yeast"
  python 4_map_crosslinks_biochemical_to_homologs.py "$files" "$output_file" "HUMAN" "YEAST"
done

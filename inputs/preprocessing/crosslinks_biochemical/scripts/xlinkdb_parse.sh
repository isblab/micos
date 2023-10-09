#!/bin/bash

######### Get XLinkDB crosslinks by type and species ##############

# Separate crosslinks based on the xlinker type
python 2_get_xlinkdb_crosslinks_by_xlinker.py

# Separate XlinkDB crosslinks by species
input_dir="../mapping_to_homologs/inputs/"

for file in "${input_dir}"DSSO.csv "${input_dir}"BDP.csv "${input_dir}"PIR.csv
do
  base_filename="$(basename "$file")"
  python 3_sort_crosslinks_by_species.py "$file" "${input_dir}human_$base_filename" "${input_dir}mouse_$base_filename" "${input_dir}yeast_$base_filename"
done

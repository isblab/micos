# #! /bin/bash

input_directory="../../../data/crosslinks/"

for file in $input_directory/*; do
  input_base_name=$(basename "$file")
  input_base_name_no_extension="${input_base_name%.*}"
  python filtering_mic19_xlinks.py "$file" "$input_base_name_no_extension"
done


directories=("../self_links" "../monomeric_links")

for directory in "${directories[@]}";do
  if [[ $directory == *"self_links"* ]]; then
      flag="self"
  elif [[ $directory == *"monomeric_links"* ]]; then
      flag="monomeric"
  else
      flag=""
  fi

  for input_file in $directory/*; do
      if [[ $input_file == *BDP* ]]; then
          length=52
      elif [[ $input_file == *BS3* ]]; then
        length=35
      else
          length=30
      fi
      input_base_name=$(basename "$input_file")
      input_base_name_no_extension="${input_base_name%.*}"

      python crosslink_satisfaction.py "$input_file" "${input_base_name_no_extension}" $length "AF_ap.pdb" "$flag"
      python crosslink_satisfaction.py "$input_file" "${input_base_name_no_extension}_CCP" $length "cccp_parallel.pdb" "$flag"
      python crosslink_satisfaction.py "$input_file" "${input_base_name_no_extension}_CAP" $length "cccp_antiparallel.pdb" "$flag"
  done
done

#!/bin/bash

files=(
    "../contact_maps/binarized_distance_matrices/MIC10.0-MIC13.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC19.0-MIC10.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC19.2-MIC10.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.0-MIC10.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.1-MIC10.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.2-MIC10.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.3-MIC10.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.0-MIC19.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.1-MIC19.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.2-MIC19.0_binarized_distances.csv"
    "../contact_maps/binarized_distance_matrices/MIC60.3-MIC19.0_binarized_distances.csv"
)

for file in "${files[@]}"; do
    # Extract protein names
    protein1=$(echo "$file" | awk -F'[-/]' '{print $4}')
    protein2=$(echo "$file" | awk -F'[-/]' '{print $5}')

    if [[ "$protein1" == "MIC10.0" || "$protein2" == "MIC13.0" ]]; then
        param1=1
    elif [[ "$protein1" == "MIC60.0" || "$protein1" == "MIC60.1" || "$protein1" == "MIC60.2" || "$protein1" == "MIC60.3" ]]; then
        param1=410
    elif [[ "$protein1" == "MIC19.0" ]]; then
        param1=1
    elif [[ "$protein1" == "MIC19.2" ]]; then
        param1=186
    fi

    param2=1

    python ../../../scripts/analysis/correcting_res_number_in_cm.py "$file" "$param1" "$param2"
done

#! /bin/bash

 ~/imp-clean/build/setup_environment.sh python ../github/scripts/analysis/extract_sampcon.py sampcon_0_extracted.rmf3 model_analysis/A_models_clust6.rmf3 sampcon/cluster.0.sample_A.txt model_analysis/B_models_clust6.rmf3 sampcon/cluster.0.sample_B.txt

~/imp-clean/build/setup_environment.sh python ~/IMP_Toolbox/analysis/surface_distance_maps_v2/contact_maps_surface_v2.py -rf sampcon_0_extracted.rmf3 -p 120 -t 10 -pt 0.25

~/imp-clean/build/setup_environment.sh python ../github/scripts/analysis/fit_to_binding_data.py sampcon_0_extracted.rmf3 ../github/inputs/data/biochemical/af3_predictions_biochemical.csv ../github/inputs/data/biochemical/validation_biochemical.csv

~/imp-clean/build/setup_environment.sh sh ../github/scripts/analysis/xl_viol.sh

~/imp-clean/build/setup_environment.sh python ~/prism/src/main.py  --input sampcon/cluster.0.prism.npz --input_type npz --output output/ --voxel_size 4 --return_spread --classes 3 --cores 8 --models 1.0 --n_breaks 50

~/imp-clean/build/setup_environment.sh python ~/prism/src/color_precision.py --resolution 1 --annotations_file output/annotations_cl2.txt --input sampcon/cluster.0/cluster_center_model.rmf3 --output output/micos_patch_colored_cluster_center_model.rmf3

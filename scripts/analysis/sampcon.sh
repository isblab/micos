#! /bin/bash

$imp  python  ~/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py -n micos -a -m cpu_omp -cc 2 -pr -d model_analysis/density.txt -gp -g 2.0  -sa model_analysis/A_models_clust1.txt -sb model_analysis/B_models_clust1.txt  -ra model_analysis/A_models_clust1.rmf3 -rb model_analysis/B_models_clust1.rmf3

echo "Done!"

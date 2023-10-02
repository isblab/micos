step1
go to xlinkdb site
search protein name in Search Name/Uniprot Accession and download all crosslink data files corresponding to this, saved in xlinkdb-Mic/
run the script get_unique_datasets_with_micos_xlinks.py to get the list of the unique datasets

step2
download these datasets from xlinkdb 4 and save in xlinkdb-datasets/
datasets ChemBiol17AAG_Bruce and Liu2015NatureMethods_Heck are not considered further
run sorting_dataset_for_micos_proteins_xlinker.py

step3
separate the crosslinks based on the species
run species_wise_sorting.py

setp4
map the crosslinks from one species to another
run mapping_crosslinks_biochemical_data_in_homologs.py

all these scripts are added in master_script.sh

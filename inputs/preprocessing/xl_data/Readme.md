### Preprocessing the crosslinks from XLinkDB 

step1
go to xlinkdb site `https://xlinkdb.gs.washington.edu/xlinkdb/`
search protein name in Search Name/Uniprot Accession and download all crosslink data files corresponding to this, saved in xlinkdb-Mic/
```
run the script python scripts/get_unique_datasets_with_micos_xlinks.py to get the list of the unique datasets
```

step2
download these datasets from xlinkdb 4 and save in xlinkdb-datasets/
datasets ChemBiol17AAG_Bruce and Liu2015NatureMethods_Heck are not considered further
```
run python scripts/sorting_dataset_for_micos_proteins_xlinker.py
```

step3
separate the crosslinks based on the species
```
run python scripts/species_wise_sorting.py ../inputs/DSSO_xl.csv human_xl.csv mouse_xl.csv yeast_xl.csv
```

setp4
map the crosslinks from one species to another
```
run python scripts/mapping_crosslinks_biochemical_data_in_homologs.py mouse_PIR_xl.csv mouse_PIR_xl_XLinkDB_mouse_to_human 'MOUSE' 'HUMAN'
```

all these scripts are added in `scripts/master_script.sh`

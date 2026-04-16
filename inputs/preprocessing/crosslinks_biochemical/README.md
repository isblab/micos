### Preprocessing the crosslinks and biochemical data

#### Biochemical data

`biochemical` This contains pairwise biochemical data studied in human and yeast. We are using only human data. No need to map yeast data to human.

#### Crosslink data

1. Filter the MICOS crosslinks from:
    a. The studies in `crosslinks/crosslink_papers` 
    b. The crosslink datasets in `crosslinks/xlinkdb`. Using 1-3 scripts in the `scripts` directory for this parsing the crosslinks.

For details, go to the respective directories mentioned above and follow the README for parsing the crosslinks from various studies. 


#### Mapping data across species and getting the inputs for modeling

2. Map the MICOS crosslinks and biochemical data from human to yeast and vice versa. Use `scripts/4_map_crosslinks_biochemical_to_homologs.py` for this.

3. Merge and save the human crosslinks (mapped to human and directly observed in human) in  `../../data/crosslinks/` to be used for modeling. Use `scripts/5_merge_files.py` for this.

Run the following script to run 4-5 scripts: 

```
scripts/mapping_to_homologs_master_script.sh
```


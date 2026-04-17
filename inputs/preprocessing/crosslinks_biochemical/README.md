### Preprocessing the crosslinks and biochemical data

#### Biochemical data

`biochemical` This contains pairwise biochemical data studied in human and yeast. We are using only human data. No need to map yeast data to human.

#### Crosslink data

1. First filter the MICOS crosslinks from various studies and from XlinkDB.
   
    a. The individual studies are in `crosslinks/crosslink_papers`

Go to the respective directories mentioned above and follow the README for filtering MICOS crosslinks from each study.

    b. The XlinkDB crosslink datasets are in `crosslinks/xlinkdb`.

Scripts `scripts/[1-3]_*` along with the wrapper script `xlinkdb_parse.sh` in the `scripts` sub-directory are for this.

#### Mapping data across species and getting the inputs for modeling

2. Map the MICOS crosslinks and biochemical data from human to yeast and vice versa. `scripts/4_map_crosslinks_biochemical_to_homologs.py` does this.

3. Merge and save the human crosslinks (mapped to human and directly observed in human) in  `../../data/crosslinks/` to be used for modeling. `scripts/5_merge_files.py` does this.

Run the following master script to run both above scripts `scripts/[4-5]_*`: 

```
scripts/mapping_to_homologs_master_script.sh
```


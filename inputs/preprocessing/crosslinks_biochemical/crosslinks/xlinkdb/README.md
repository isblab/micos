### Preprocessing the crosslinks from XLinkDB

#### Step1
Go to xlinkdb site `https://xlinkdb.gs.washington.edu/xlinkdb/`.
Search protein name in Search Name/Uniprot Accession and download all crosslink data files corresponding to the protein and save in `xlinkdb_Mic/`.
Get the names of all the datasets that have crosslinks for MICOS proteins using:
```
../../scripts/1_get_unique_xlinkdb_datasets_with_micos_xlinks.py
```

#### Step 2
Download all the crosslink datasets from xlinkdb which have MICOS proteins and save in `xlinkdb_datasets/`.

#### Step 3
Separate the crosslinks based on the xlinker type and the species using:

```
../../scripts/xlinkdb_parse.sh
```

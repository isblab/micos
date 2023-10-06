### Preprocessing the crosslinks and biochemical data

#### Biochemical data

`biochemical`
Contains pairwise biochemical data studied in human and yeast in pairwise_biochemical_data.csv


#### Crosslink data

1. Filter the micos crosslinks from the studies in `crosslinks/crosslink_papers` and the crosslink datasets in `crosslinks/xlinkdb`.

2. Map the micos crosslinks from human to yeast and vice versa.

3. Save the crosslinks mapped to human and studied in human in  `../../data/crosslinks/human` to be used in modeling.

For step 2 and 3, run the script: 

```
master_script.sh
```


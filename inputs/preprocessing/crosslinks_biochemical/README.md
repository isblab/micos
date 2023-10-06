### Preprocessing the crosslinks and biochemical data

#### Biochemical data

`biochemical`

This contains pairwise biochemical data studied in human and yeast.


#### Crosslink data

1. Filter the micos crosslinks from the studies in `crosslinks/crosslink_papers` and the crosslink datasets in `crosslinks/xlinkdb`.

2. Map the micos crosslinks from human to yeast and vice versa.

3. Save the crosslinks mapped and studied in human in  `../../data/crosslinks/human` to be used for modeling.

Run the following script to map the biochemical and crosslink data: 

```
master_script.sh
```


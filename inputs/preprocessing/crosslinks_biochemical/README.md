### Preprocessing the crosslinks and biochemical data

#### Biochemical data

`biochemical` This contains pairwise biochemical data studied in human and yeast.

#### Crosslink data

Filter the micos crosslinks from the studies in `crosslinks/crosslink_papers` and the crosslink datasets in `crosslinks/xlinkdb`.
For this, go to the respective directories mentioned above and follow the README for parsing the crosslinks from various studies. 

#### Mapping data across species and getting the inputs for modeling

Map the micos crosslinks and biochemical data from human to yeast and vice versa. Save the human crosslinks (mapped to human and directly observed in human) in  `../../data/crosslinks/human` to be used for modeling.

Run the following script do the above: 

```
master_script.sh
```


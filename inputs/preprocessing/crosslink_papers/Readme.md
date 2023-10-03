### Preprocessing the crosslinks from papers
 
#### These studies are done on human cells

`bartolec_2023`
download supplementary table 1 excel file and 4th and 5th sheet as csv files corresponding to DSSO and DHSO crosslinks, respectively.
```
run bartolec_filtering_data.py to get two csv files with micos crosslinks
```
remove 3 crosslinks corresponding to MICA and MICU proteins manually from the DSSO crosslink dataset.

`ryl_2020`
download supplementary table and save S2B sheet.
```
run ryl_filtering_data.py
```
`yugandhar_2020`
download supplementary zip file, extract 156740_0_supp_420642_q05ccw.xlsx and save 'supplementary table 4' sheet as csv file.
```
run yugandhar_filtering_data.py
```

`sun_2022`
download supplementary table 3 excel file and save 3rd and 8th sheet as csv files.
manually added the three croslinks in a csv file.

```
Linden_2020 study done in yeast, while aligning to the human, they map to gaps. Not considered further.
```

### Preprocessing the crosslinks from papers

#### These studies are done on human cells

`bartolec_2023`

Download supplementary table 1 excel file and 4th and 5th sheet as csv files corresponding to DSSO and DHSO crosslinks, respectively.
Run the script:
```
bartolec_filtering_data.py
```
Remove 3 crosslinks corresponding to MICA and MICU proteins manually from the DSSO crosslinks file.

`ryl_2020`

Download supplementary table and save S2B sheet.
Run the script:
```
ryl_filtering_data.py
```

`yugandhar_2020`

Download supplementary zip file, extract 156740_0_supp_420642_q05ccw.xlsx and save 'supplementary table 4' sheet as csv file.
Run the script:
```
run yugandhar_filtering_data.py
```

`sun_2022`

Download supplementary table 3 excel file and search manually for micos crosslinks. There are two crosslinks in 3rd sheet and one crosslink in 8th sheet, which are manually saved in a csv file.

####### linden_2020 study was done in yeast. Crosslink residues were aligned to gap in human species, so not considered further.

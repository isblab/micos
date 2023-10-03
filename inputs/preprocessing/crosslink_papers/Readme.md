### Preprocessing the crosslinks from papers
 
#### these studies are done on human cells

`Bartolec_2023`
download supplementary table 1 excel file and 4th and 5th sheet as csv files corresponding to DSSO and DHSO crosslinks, respectively.
run bartolec_filtering_data.py to get two csv files with micos crosslinks
remove 3 crosslinks corresponding to MICA and MICU proteins manually from the DSSO crosslink dataset.

`Ryl_2020`
download supplementary table and save S2B sheet.
run ryl_filtering_data.py

`Yugandhar_2020`
download supplementary zip file, extract 156740_0_supp_420642_q05ccw.xlsx and save 'supplementary table 4' sheet as csv file.
run yugandhar_filtering_data.py

`Sun_2022`
download supplementary table 3 excel file and save 3rd and 8th sheet as csv files.
manually added the three croslinks in a csv file.

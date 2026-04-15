### Preprocessing the crosslinks from papers

** Note: ** The filtered crosslinks from these studies (1-4) are saved in `../mapping_to_homologs/inputs` for the next steps.

#### Studies on human proteins

1. `bartolec_2023`

Download supplementary table 1 excel file and 4th and 5th sheet as csv files corresponding to DSSO and DHSO crosslinks, respectively.
Run the script:
```
bartolec_2023/bartolec_filtering_data.py
```
Remove 3 crosslinks corresponding to MICA and MICU proteins manually from the DSSO crosslinks file.

2. `ryl_2020`

Download supplementary table and save S2B sheet.
Run the script:
```
ryl_2020/ryl_filtering_data.py
```

3. `yugandhar_2020`

Download supplementary zip file, extract 156740_0_supp_420642_q05ccw.xlsx and save 'supplementary table 4' sheet as csv file.
Run the script:
```
yugandhar_2020/yugandhar_filtering_data.py
```

4. `sun_2022`

Download supplementary table 3 excel file and search manually for micos crosslinks. There are two crosslinks in 3rd sheet and one crosslink in 8th sheet, which are manually saved in a csv file.

5. `zhu_2024`

Download the supplementary table and parse through 'Supp Data 1_Combined (DSSO)' sheet to sort MICOS and non-MICOS crosslinks.
Run the script:
```
zhu_2024/zhu_filtering_data.py
```
These crosslinks were used only for validation. 

##### Studies on yeast (not used)
`linden_2020` study was done in yeast. These crosslink residues were aligned to gap in MSA in human species.  Therefore they were not considered further.

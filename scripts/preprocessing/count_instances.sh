#!/bin/bash

# # Output direc

### Yeast
for file in *_to_yeast; do
  count=$(awk -F',' 'NR > 1 && $1 != $3' "$file" | wc -l)

  echo "$file: $count"
done

# bartolec_DHSO_human_to_yeast: 0
# bartolec_DSSO_human_to_yeast: 6
# human_BDP_XLinkDB_human_to_yeast: 0
# human_DSSO_XLinkDB_human_to_yeast: 1
# human_PIR_XLinkDB_human_to_yeast: 4
# mouse_BDP_XLinkDB_mouse_to_yeast: 4
# mouse_DSSO_XLinkDB_mouse_to_yeast: 0
# mouse_PIR_XLinkDB_mouse_to_yeast: 8
# pairwise_biochemical_data_human_to_yeast: 11
# ryl_BS3_human_to_yeast: 0
# sun_DSS_human_to_yeast: 0
# yu_DSSO_human_to_yeast: 1


### Human
for file in *_to_human; do
  count=$(awk -F',' 'NR > 1 && $1 != $3' "$file" | wc -l)

  echo "$file: $count"
done

# bartolec_DSSO_human_to_human: 7
# linden_BS3_yeast_to_human: 0
# mouse_BDP_XLinkDB_mouse_to_human: 10
# mouse_DSSO_XLinkDB_mouse_to_human: 0
# mouse_PIR_XLinkDB_mouse_to_human: 10
# pairwise_biochemical_data_yeast_to_human: 11
# yeast_DSSO_XLinkDB_yeast_to_human: 0

# Input direc
for file in *.csv; do
  count=$(awk -F',' 'NR > 1 && $1 != $3' "$file" | wc -l)

  echo "$file: $count"
done
### Analysis of parallel and antiparallel configuration of Mic19 coiled-coil domain

The Mic19-Mic19 crosslinks in 59-174 region from BDP_PIR_mouse, DHSO_DSSO_human, DSS_BS3_human and DSSO_XLinkDB_mouse datasets are analysed on anti-parallel coiled-coil dimer predicted from AF-multimer and parallel and anti-parallel coiled-coil dimer modeled using CCCP server. The BDP_PIR_human dataset doesnot have Mic19-Mic19 crosslink in 59-174 region.

1. `distances` contains the satisfied and violated distances.

2. `monomeric_links` contains the crosslinks which are on different residues i.e., Protein1 Residue1 Protein1 Residue2.

3. `self_links` contains the crosslinks which are present on same residue i.e., Protein1 Residue1 Protein1 Residue1.

4. `servers` contains the outputs from COILS, MARCOIL, PCOILS and Psipred servers and coiled-coil dimers predicted by CCCP and AF-multimer.

5. `scripts` contains the script `filtering_mic19_xlinks.py` to filter Mic19-Mic19 crosslinks from the datasets and `crosslink_satisfaction.py` to calculate the percentage of crosslinks satisfied. 




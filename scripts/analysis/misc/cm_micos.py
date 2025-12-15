import warnings
import numpy as np
import pandas as pd
from af_pipeline.utils import get_patches_from_matrix
import glob, os

def get_interacting_patches(
        contact_map: np.array,
        region_of_interest: dict,
    ):
        """This is a dirty implementation to get the interacting patches. \n
        This is a temporary solution until we find a better way to get interacting
        patches for the given contact map.

        Args:
            contact_map (np.array): binary contact map.
            region_of_interest (dict): region of interest for the protein pair.

        Returns:
            patches (dict): interacting patches for the given region of interest of the protein pair.
        """

        patches = {}

        chain1, chain2 = region_of_interest.keys()
        p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

        if np.unique(contact_map).tolist() == [0]: # No interactions found.
            warnings.warn(
                f"No interacting patches found for {chain1}:{p1_region} and {chain2}:{p2_region}."
            )
            return patches

        patches_df = get_patches_from_matrix(
            matrix=contact_map,
            chain1=chain1,
            chain2=chain2
        )

        for patch_idx, patch in patches_df.iterrows():

            ch1_patch = patch[chain1]
            ch2_patch = patch[chain2]

            ch1_patch = sorted([int(x) for x in ch1_patch])
            ch2_patch = sorted([int(x) for x in ch2_patch])

            ch1_patch = np.array(ch1_patch) + region_of_interest[chain1][0]
            ch2_patch = np.array(ch2_patch) + region_of_interest[chain2][0]

            patches[patch_idx] = {
                chain1: np.array(ch1_patch),
                chain2: np.array(ch2_patch),
            }

        return patches, chain1, chain2 # In the OMG's version, function returns only patches.


def get_filtered_patch(p1, p2, patches, domain_list):
     
     # For example, Mic10 = [11,12,13,14,15,16] and Mic13 = [19, 20, 21, 22, 23] are making a contact
        # But it is inconsistent with membrane topology.       
        # superset1 is [IMS,TM] and superset2 is [TM]
        # Then we find the common domain = [TM] and trim the protein regions according to the TM region.
        # So we will keep Mic10 = [13, 14,15,16] and Mic13 =[19, 20, 21, 22, 23] as contact.
     

     to_remove = [] # To store the patches inconsistent with membrane topology or has less than 5 residues of any protein in a patch
     
     for patch in patches:
        ph1 = patches[patch][p1]
        ph2 = patches[patch][p2]

        p_domains1, p_domains2 = domain_list

        # Find all the domains of each protein in patch        
        superset1 = [
            (k)
            for k, v in p_domains1.items()
            for b in v
            if any(x in range(b[0], b[1] + 1) for x in ph1)
        ]

        superset2 = [
            (k)
            for k, v in p_domains2.items()
            for b in v
            if any(x in range(b[0], b[1] + 1) for x in ph2)
        ]
        
        # Get the common domains. 
        common_domains = [(s) for s in superset1 if s in superset2]
                

        def trim_patch(ph, domains):
            trimmed = []

            for ranges in domains:
                start, end = ranges[0], ranges[1]
                overlap = [int(x) for x in ph if start <= x <= end]
                if overlap:
                    trimmed.extend(overlap)

            return sorted(set(trimmed))

        if len(common_domains) == 0:
            # Remove those patches
            to_remove.append(patch)

        else:
            # Trim patches
            ph1_trimmed = trim_patch(ph1, p_domains1[common_domains[0]])
            ph2_trimmed = trim_patch(ph2, p_domains2[common_domains[0]])

            patches[patch][p1] = ph1_trimmed
            patches[patch][p2] = ph2_trimmed
            
            # # Remove patches with less than 5 residues in both the proteins. Not using this. 
            # if len(ph1_trimmed) > 5 or len(ph2_trimmed) > 5:
            #     patches[patch][p1] = ph1_trimmed
            #     patches[patch][p2] = ph2_trimmed
            # else:
            #     to_remove.append(patch)

     patches = {k: v for k, v in patches.items() if k not in to_remove}

     return patches
     


#### To obtain contacts filtered by membrane regions ###

# Modeled regions
regions_of_interest = {'MIC10': [1, 78], 'MIC13': [1, 118], 'MIC60': [410, 758], 'MIC19': [1, 227]}

# Membrane regions of each modeled protein
domains = {'MIC10': {'IMS': [[1,12], [61, 78]], 'TM': [[13,36], [40, 60]], 'Matrix': [[37, 39]] }, 
            'MIC13': {'IMS': [[24, 118]], 'TM': [[8, 23]], 'Matrix': [[1,7]]} , 
            'MIC19': {'IMS': [[1,227]]},
            'MIC60': {'IMS': [[410,758]]}}

cm_folder = '/home/muskaan/projects/micos/results/contact_maps/interprotein/distance_matrices'

for csv in glob.glob(os.path.join(cm_folder, "*")):
    
    df = pd.read_csv(csv, index_col=False) 
    if 'Unnamed: 0' in df.columns:
        df = df.drop('Unnamed: 0', axis=1)
    mean_dist = df.to_numpy()

    bin_csv = np.where(mean_dist <= float(10), 1, 0)
    
    # Get protein names
    p = []
    if 'MIC10' in csv:
        p.append('MIC10')
    if 'MIC13' in csv:
        p.append('MIC13')
    if 'MIC60' in csv:
        p.append('MIC60')
    if 'MIC19' in csv:
        p.append('MIC19')
    
    assert len(p) == 2 ; 'More than two protein names in the contact map file'

    # create the dict of protein and their region of interest
    roi = { p[0]: regions_of_interest[p[0]], 
            p[1]: regions_of_interest[p[1]],}
        
    patches, p1, p2 = get_interacting_patches(bin_csv, roi)
    domain_list = [domains[p1], domains[p2]]
    filtered_patch = get_filtered_patch(p1, p2, patches, domain_list)

    df = pd.DataFrame.from_dict(filtered_patch, orient='index')

    for col in df.columns:
        df[col] = df[col].apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)

    df.to_excel(f"/home/muskaan/projects/micos/results/contact_maps/{p1}_{p2}_patches.xlsx", index=True)
    exit()
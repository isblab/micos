# Originally written by Shreyas Arvindekar. Modified by stochastic13 (Satwik)
'''
Algorithm
Examine a set of multipliers (mean, mean-0.25 std, .. ) for set of data restraints
Choose the largest number of models (lowest multiplier)
    [such that nA and nB are each less than 15k: earlier ] about the same as : nA +nB < 20k
    It passes the KS test on total score and (A,B) looks similar on score distribution
    If you dont find a multiplier even after the narrower search in point 1, [take a random subset of the nearest multiplier that passes the KS test] OR take a single lenient cutoff on EV (less than mean) along with score multipliers on the other restraints.
'''

import os
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tabulate import tabulate
from scipy.stats import ks_2samp


# This function will take the individual dataframes and compare them with the mean-multiplier*std of common dataframe and output the dataframe of the models that satisfy the filter
def variable_filter(std_multiplier, df, score_list):
    temp = None
    for i in score_list:
        if temp is None:
            temp = df[i] <= (common_df_mean[i] + std_multiplier * common_df_std[i])
        else:
            temp = temp & (df[i] <= (common_df_mean[i] + std_multiplier * common_df_std[i]))
    return df[temp]


parser = argparse.ArgumentParser()
# parser.add_argument('-e', '--evr', action='store_true', default=False,
# help='Shall I use the Excluded Volume Restraint in the filter?')
parser.add_argument('-c', '--cluster_num', default=0, help='On which cluster shall I run the filter?')
parser.add_argument('-lc', '--lowest_cutoff', default=-2,
                    help='What is the standard deviation multiplier for the most stringent cutoff?')
parser.add_argument('-hc', '--highest_cutoff', default=3,
                    help='What is the standard deviation multiplier for the most lenient cutoff?')
parser.add_argument('-ss', '--step_size', default=0.01, help='What step size shall I use?')
parser.add_argument('-n', '--num_models', default=30000, type=int,
                    help='What is the maximum number of models you want?')
parser.add_argument('-g', '--gsmsel', default=os.getcwd() + "/model_analysis/", type=str,
                    help='Where is the gsm_sel directory')

# print("If you are using the EVR flag, consider reducing the step size")

args = parser.parse_args()

cluster_num = str(args.cluster_num)
lowest_cutoff = float(args.lowest_cutoff)
highest_cutoff = float(args.highest_cutoff)
step_size = float(args.step_size)
num_models = int(args.num_models)
gsm_sel_dir = str(args.gsmsel)
# print(args)
sys_name = 'ignore'

data_restraint_names = []

cluster_csv_fileA = gsm_sel_dir + 'selected_models_A_cluster' + cluster_num + '_detailed.csv'
cluster_csv_fileB = gsm_sel_dir + 'selected_models_B_cluster' + cluster_num + '_detailed.csv'

std_mult_dtrst = []
for i in np.arange(highest_cutoff, lowest_cutoff - step_size, -step_size):
    std_mult_dtrst.append(round(i, 3))
print(f'Std-multiplier list: {std_mult_dtrst}')

# Reading the CSV files and combining them to apply a common cutoff
columns_to_ignore = ["traj", "rmf3_file", "half", "cluster", "frame_RMF3"]
dfA = pd.read_csv(cluster_csv_fileA, usecols=lambda column: column not in columns_to_ignore )
dfB = pd.read_csv(cluster_csv_fileB, usecols=lambda column: column not in columns_to_ignore )
print('Loaded the csv files')
df_list = [dfA, dfB]
common_df = pd.concat(df_list, ignore_index=True)
print('Calculating mean')
common_df_mean = common_df.mean()
print('Calculating standard deviation')
common_df_std = common_df.std()
print('Running variable filter')
out = []
out_str = f'Cluster number: {cluster_num} \nLowest cutoff: {lowest_cutoff} \nHighest cutoff: {highest_cutoff} \nStep size: {step_size} \nMaximum number of models to be selected: {num_models} \n\n'
mult_found = False

dfA = pd.read_csv(cluster_csv_fileA)
dfB = pd.read_csv(cluster_csv_fileB)

default_restraints = [
    'EV_sum',
    'XLs_sum',
    'ILR_sum',
    'ZAR_sum',
    'TMR_sum',
    'MLR_sum',
    'MPDBR_sum',
    'Total_Score'
]

best_multiplier = None
best_dfA = None
best_dfB = None
best_total = 0
best_ksd = None
best_ksp = None

for multiplier in std_mult_dtrst:
    sel_dfA = variable_filter(
        multiplier, dfA, default_restraints
    )

    sel_dfB = variable_filter(
        multiplier, dfB, default_restraints
    )

    nModelsA = len(sel_dfA)
    nModelsB = len(sel_dfB)
    nModelsT = nModelsA + nModelsB

    # Cannot run KS test if either set is empty
    if nModelsA == 0 or nModelsB == 0:
        print(
            f"{multiplier:.2f}: "
            f"A={nModelsA}, B={nModelsB}, "
            f"total={nModelsT} -- empty set"
        )
        continue

    scoresA = sel_dfA['Total_Score'].values
    scoresB = sel_dfB['Total_Score'].values

    ksd_pval = ks_2samp(scoresA, scoresB)
    ksd = ksd_pval.statistic
    ksp = ksd_pval.pvalue

    ks_pass = (
        ksp > 0.05
        or
        (ksp <= 0.05 and ksd < 0.3)
    )

    print(
        f"{multiplier:.2f}: "
        f"A={nModelsA}, B={nModelsB}, "
        f"total={nModelsT}, "
        f"KS-D={ksd:.3f}, KS-p={ksp:.4g}, "
        f"KS={'PASS' if ks_pass else 'FAIL'}"
    )

    # Save the current passing multiplier
    if nModelsT >= num_models and ks_pass:

        best_multiplier = multiplier
        best_dfA = sel_dfA.copy()
        best_dfB = sel_dfB.copy()
        best_total = nModelsT
        best_ksd = ksd
        best_ksp = ksp

    # Since we are going from high -> low,
    # the first failure after a passing multiplier
    # means the previous multiplier was the most stringent passing one.
    elif best_multiplier is not None:

        print(
            f"\nStopping search at multiplier {multiplier:.2f}. "
            f"Most stringent passing multiplier is {best_multiplier:.2f}."
        )
        break


if best_multiplier is None:

    print(
        "\nNo multiplier found that provides "
        f"at least {num_models} models AND passes KS."
    )

else:

    print("\nOptimal multiplier found.")
    print(f"Multiplier: {best_multiplier}")
    print(f"Available A models: {len(best_dfA)}")
    print(f"Available B models: {len(best_dfB)}")
    print(f"Available total: {best_total}")
    print(f"KS-D: {best_ksd}")
    print(f"KS-p: {best_ksp}")

    # ------------------------------------------
    # Downsample to exactly num_models
    # ------------------------------------------

    target_A = round(
        num_models * len(best_dfA) / best_total
    )

    target_B = num_models - target_A

    best_dfA = best_dfA.sample(
        n=target_A,
        random_state=42
    )

    best_dfB = best_dfB.sample(
        n=target_B,
        random_state=42
    )

    print("\nFinal selection:")
    print(f"A: {len(best_dfA)}")
    print(f"B: {len(best_dfB)}")
    print(f"Total: {len(best_dfA) + len(best_dfB)}")

    best_dfA.to_csv(
        cluster_csv_fileA.replace(
            'selected_models',
            'good_scoring_models'
        ),
        index=False
    )

    best_dfB.to_csv(
        cluster_csv_fileB.replace(
            'selected_models',
            'good_scoring_models'
        ),
        index=False
    )
    
    with open('var_filt_out.log', 'w') as outf:
        outf.write(out_str)
        outf.write(
            f'\nOptimal filter found.\n'
            f'Extracted at {best_multiplier}\n'
            f'A models: {len(best_dfA)}\n'
            f'B models: {len(best_dfB)}\n'
            f'Total models: {len(best_dfA) + len(best_dfB)}\n'
            f'KS-D: {best_ksd}\n'
            f'KS-p: {best_ksp}\n'
        )

    # Plot score distributions
    scoresA = best_dfA['Total_Score'].values
    scoresB = best_dfB['Total_Score'].values

    nBins = int(
        max(scoresA.max(), scoresB.max())
        -
        min(scoresA.min(), scoresB.min())
    )

    plt.figure()
    plt.hist(scoresA, bins=nBins, histtype='step', label='ScoresA')
    plt.hist(scoresB, bins=nBins, histtype='step', label='ScoresB')
    plt.title('Scores of sampleA and sampleB')
    plt.xlabel('Total Score')
    plt.ylabel('nModels')
    plt.legend()
    plt.savefig('var_filt_out.png')
    plt.show()

# The old mult_found block is no longer needed because
# best_multiplier is used above.
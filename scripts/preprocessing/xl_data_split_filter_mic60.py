import pandas as pd
import glob
import os

folder = "."

evicalc_files = glob.glob(os.path.join(folder, "evicalc_*.csv"))
sampling_files = glob.glob(os.path.join(folder, "sampling_*.csv"))

names = {
    os.path.basename(f).split("_", 1)[1]
    for f in evicalc_files + sampling_files
}


def filter_df(df):
    """Remove MIC60 residues < 451 from a crosslink DataFrame."""
    remove = (
        ((df["Protein1"] == "MIC60") & (df["Residue1"] < 451)) |
        ((df["Protein2"] == "MIC60") & (df["Residue2"] < 451))
    )
    removed = remove.sum()
    df = df[~remove].copy()
    return df, removed


def filter_file(filepath):
    """Read, filter, and overwrite a crosslink file."""
    df = pd.read_csv(filepath)
    df, removed = filter_df(df)
    df.to_csv(filepath, index=False)
    return df, removed

def count_interprotein(df):
    """Count crosslinks between different proteins."""
    return (df["Protein1"] != df["Protein2"]).sum()

for name in sorted(names):

    evicalc_file = os.path.join(folder, "evicalc_" + name)
    sampling_file = os.path.join(folder, "sampling_" + name)

    evicalc_exists = os.path.exists(evicalc_file)
    sampling_exists = os.path.exists(sampling_file)

    # --------------------------------------------------
    # Both sampling and evicalc files exist
    # Merge -> filter -> split 25/75
    # --------------------------------------------------

    if evicalc_exists and sampling_exists:

        df1 = pd.read_csv(evicalc_file)
        df2 = pd.read_csv(sampling_file)

        df = pd.concat([df1, df2], ignore_index=True)

        # Use the same filtering function
        df, removed = filter_df(df)

        # Shuffle
        df = df.sample(
            frac=1,
            random_state=42
        ).reset_index(drop=True)

        # Split 25/75
        n_evicalc = int(len(df) * 0.25)

        evicalc = df.iloc[:n_evicalc]
        sampling = df.iloc[n_evicalc:]

        evicalc.to_csv(evicalc_file, index=False)
        sampling.to_csv(sampling_file, index=False)

        print(name)
        print("  Both files present")
        print("  Removed:", removed)
        print("  Total:", len(df))
        print("  Evicalc:", len(evicalc))
        print("  Sampling:", len(sampling))
        print("  Evicalc inter-protein:", count_interprotein(evicalc))
        print("  Sampling inter-protein:", count_interprotein(sampling))

    else:
        # Filter it, but do NOT split. This is for evicalc_dhso_dsso_mouse as there are 2 xlinks, 
        # sampling bdp_pir_human has 13 xlinks and sampling_dsso_mouse has only intra-mic19 xlinks
        
        filepath = evicalc_file if evicalc_exists else sampling_file
        file_type = "evicalc" if evicalc_exists else "sampling"

        df, removed = filter_file(filepath)

        print(name)
        print(f"  Only {file_type} present")
        print("  Removed:", removed)
        print("  Remaining:", len(df))
        print("  Inter-protein:", count_interprotein(df))
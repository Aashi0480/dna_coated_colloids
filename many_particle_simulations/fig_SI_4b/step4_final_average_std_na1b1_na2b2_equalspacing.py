import pandas as pd
import numpy as np

# Correct z ranges
def group_atoms_by_z_range(z):
    if z < 7:
        return "z<7"
    elif 7 <= z < 15:
        return "7<=z<15"
    elif 15 <= z < 23:
        return "15<=z<23"
    elif 23 <= z < 31:
        return "23<=z<31"
    elif 31 <= z < 39:
        return "31<=z<39"
    elif 39 <= z < 47:
        return "39<=z<47"
    else:
        return "z>47"


# Load the data into a DataFrame and process
def process_file(file_path, output_file):

    N = 20   # number of snapshots

    # Read the input file
    data = pd.read_csv(
        file_path,
        sep=r"\s+",
        header=None,
        names=["timestep", "z", "a1b1_a2b2_sur_2sbrg"]
    )

    # Ensure numeric data types
    data["z"] = pd.to_numeric(data["z"], errors="coerce")
    data["a1b1_a2b2_sur_2sbrg"] = pd.to_numeric(
        data["a1b1_a2b2_sur_2sbrg"], errors="coerce"
    )

    # Remove invalid rows
    data.dropna(subset=["z", "a1b1_a2b2_sur_2sbrg"], inplace=True)

    # Assign z-range group
    data["z_range"] = data["z"].apply(group_atoms_by_z_range)

    # Compute mean and std
    result = data.groupby("z_range").agg({
        "z": ["mean", "std"],
        "a1b1_a2b2_sur_2sbrg": ["mean", "std"]
    }).reset_index()

    # Flatten column names
    result.columns = [
        "z_range",
        "avg_z",
        "std_z",
        "avg_a1b1_a2b2_sur_2sbrg",
        "std_a1b1_a2b2_sur_2sbrg"
    ]

    # Replace NaN std (single-point group case)
    result["std_z"] = result["std_z"].fillna(0)
    result["std_a1b1_a2b2_sur_2sbrg"] = result[
        "std_a1b1_a2b2_sur_2sbrg"
    ].fillna(0)

    # Compute error bars
    result["err_z"] = result["std_z"] / np.sqrt(N)
    result["err_a1b1_a2b2_sur_2sbrg"] = (
        result["std_a1b1_a2b2_sur_2sbrg"] / np.sqrt(N)
    )

    # Sort by avg_z
    result = result.sort_values(by="avg_z", ascending=True)

    # Save output
    result.to_csv(output_file, index=False, sep="\t")

    print(f"Processed data saved to {output_file}")


# run
process_file(
    "na1b1_na2b2_averages.txt",
    "na1b1_na2b2_sum_averages.txt"
)

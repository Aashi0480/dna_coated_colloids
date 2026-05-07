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

    # Read file
    data = pd.read_csv(
        file_path,
        sep=r"\s+",
        header=None,
        names=["timestep", "z", "b1b2a1_b1b2a2"]
    )

    # Ensure numeric data types
    data["z"] = pd.to_numeric(data["z"], errors="coerce")
    data["b1b2a1_b1b2a2"] = pd.to_numeric(
        data["b1b2a1_b1b2a2"], errors="coerce"
    )

    # Remove NaN rows
    data.dropna(subset=["z", "b1b2a1_b1b2a2"], inplace=True)

    # Assign z range groups
    data["z_range"] = data["z"].apply(group_atoms_by_z_range)

    # Compute mean and std
    result = data.groupby("z_range").agg({
        "z": ["mean", "std"],
        "b1b2a1_b1b2a2": ["mean", "std"]
    }).reset_index()

    # Flatten column names
    result.columns = [
        "z_range",
        "avg_z",
        "std_z",
        "avg_b1b2a1_b1b2a2",
        "std_b1b2a1_b1b2a2"
    ]

    # Replace NaN std (single value case)
    result["std_z"] = result["std_z"].fillna(0)
    result["std_b1b2a1_b1b2a2"] = result[
        "std_b1b2a1_b1b2a2"
    ].fillna(0)

    # Error bars = sigma / sqrt(20)
    result["err_z"] = result["std_z"] / np.sqrt(N)
    result["err_b1b2a1_b1b2a2"] = (
        result["std_b1b2a1_b1b2a2"] / np.sqrt(N)
    )

    # Sort by avg_z
    result = result.sort_values(by="avg_z", ascending=True)

    # Save output
    result.to_csv(output_file, index=False, sep="\t")

    print(f"Processed data saved to {output_file}")


# run
process_file(
    "b1b2a1_b1b2a2_averages.txt",
    "b1b2a1_b1b2a2_sum_averages.txt"
)

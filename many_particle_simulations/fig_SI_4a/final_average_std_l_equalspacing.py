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

    # Read input file
    data = pd.read_csv(
        file_path,
        sep=r"\s+",
        header=None,
        names=["timestep", "z", "l"]
    )

    # Ensure numeric values
    data["z"] = pd.to_numeric(data["z"], errors="coerce")
    data["l"] = pd.to_numeric(data["l"], errors="coerce")

    # Remove invalid rows
    data.dropna(subset=["z", "l"], inplace=True)

    # Assign z ranges
    data["z_range"] = data["z"].apply(group_atoms_by_z_range)

    # Calculate mean and std
    result = data.groupby("z_range").agg({
        "z": ["mean", "std"],
        "l": ["mean", "std"]
    }).reset_index()

    # Flatten columns
    result.columns = [
        "z_range",
        "avg_z",
        "std_z",
        "avg_l",
        "std_l"
    ]

    # Replace NaN std (single-value group case)
    result["std_z"] = result["std_z"].fillna(0)
    result["std_l"] = result["std_l"].fillna(0)

    # Calculate error bars
    result["err_z"] = result["std_z"] / np.sqrt(N)
    result["err_l"] = result["std_l"] / np.sqrt(N)

    # Sort by z
    result = result.sort_values(by="avg_z", ascending=True)

    # Save output
    result.to_csv(output_file, index=False, sep="\t")

    print(f"Processed data saved to {output_file}")


# run
process_file(
    "l_averages.txt",
    "l_sum_averages.txt"
)

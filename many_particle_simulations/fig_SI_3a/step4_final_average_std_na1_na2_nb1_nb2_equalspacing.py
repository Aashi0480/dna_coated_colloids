import pandas as pd
import numpy as np

# Define the ranges
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
        names=["timestep", "z", "na1_na2_nb1_nb2"]
    )

    # Ensure numeric data types
    data["z"] = pd.to_numeric(data["z"], errors="coerce")
    data["na1_na2_nb1_nb2"] = pd.to_numeric(
        data["na1_na2_nb1_nb2"], errors="coerce"
    )

    # Remove NaN rows
    data.dropna(subset=["z", "na1_na2_nb1_nb2"], inplace=True)

    # Assign z ranges
    data["z_range"] = data["z"].apply(group_atoms_by_z_range)

    # Mean and standard deviation
    result = data.groupby("z_range").agg({
        "z": ["mean", "std"],
        "na1_na2_nb1_nb2": ["mean", "std"]
    }).reset_index()

    # Flatten multi-level columns
    result.columns = [
        "z_range",
        "avg_z",
        "std_z",
        "avg_na1_na2_nb1_nb2",
        "std_na1_na2_nb1_nb2"
    ]

    # Replace NaN std (single value group case)
    result["std_z"] = result["std_z"].fillna(0)
    result["std_na1_na2_nb1_nb2"] = result[
        "std_na1_na2_nb1_nb2"
    ].fillna(0)

    # Error bars
    result["err_z"] = result["std_z"] / np.sqrt(N)
    result["err_na1_na2_nb1_nb2"] = (
        result["std_na1_na2_nb1_nb2"] / np.sqrt(N)
    )

    # Sort by avg_z
    result = result.sort_values(by="avg_z", ascending=True)

    # Save output
    result.to_csv(output_file, index=False, sep="\t")

    print(f"Processed data saved to {output_file}")


# run
process_file(
    "na1_na2_nb1_nb2_averages.txt",
    "na1_na2_nb1_nb2_sum_averages.txt"
)

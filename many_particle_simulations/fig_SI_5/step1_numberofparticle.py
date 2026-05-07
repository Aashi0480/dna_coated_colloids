import pandas as pd
import matplotlib.pyplot as plt
import os

# Function to read the dump file and extract relevant data
def read_dump_file(dump_file):
    data = []
    timestep_count = {}
    with open(dump_file, 'r') as file:
        lines = file.readlines()
        
        timestep = None
        in_atoms_section = False
        columns = None

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("ITEM: TIMESTEP"):
                i += 1
                if i < len(lines):
                    timestep = float(lines[i].strip())  # Use float instead of int
                    timestep = round(timestep)  # Convert to nearest integer
                    if timestep not in timestep_count:
                        timestep_count[timestep] = 0
            elif line.startswith("ITEM: ATOMS"):
                columns = line.split()[2:]  # Extract column names
                in_atoms_section = True
            elif line.startswith("ITEM:") and in_atoms_section:
                in_atoms_section = False
            elif in_atoms_section:
                atom_data = list(map(float, line.split()))
                atom_data.append(timestep)
                data.append(atom_data)
                timestep_count[timestep] += 1  # Count atoms for the current timestep
            i += 1
    
    # Create a DataFrame with extracted data
    df = pd.DataFrame(data, columns=columns + ["timestep"])
    return df, timestep_count

# Read the dump file

#dump_file = 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-06_1000_modified.dump'
#dump_file = 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na200_modified.dump'
dump_file = 'long_Gillespie_BD_negno_test18_AB_DGT0_DG_-20_rho_ideal_1.0e-05_1000.dump'
#dump_file = 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_modified.dump'
#dump_file = 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-08_1000.dump'
#dump_file = 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na50_modified.dump'


df, timestep_count = read_dump_file(dump_file)

# Extract the timestep and atom counts
timesteps_sorted = sorted(timestep_count.keys())
atom_counts = [timestep_count[t] for t in timesteps_sorted]

# Generate output file name based on the input file
base_name = os.path.splitext(dump_file)[0]  # Remove the file extension
output_file = f"{base_name}_timestep_vs_atom_counts.txt"

# Save to a .txt file
output_df = pd.DataFrame({'Timestep': timesteps_sorted, 'Number of Atoms': atom_counts})
output_df.to_csv(output_file, sep='\t', index=False)  # Save as a tab-separated file

# Plot the number of atoms as a function of timestep
plt.figure(figsize=(12, 8))
plt.plot(timesteps_sorted, atom_counts, marker='None', linestyle='-', color='blue', label='Number of Atoms')
plt.xlabel('Timestep')
plt.ylabel('Number of Atoms')
plt.title('Number of Atoms vs. Timestep')
plt.legend()
plt.grid(False)
plt.show()

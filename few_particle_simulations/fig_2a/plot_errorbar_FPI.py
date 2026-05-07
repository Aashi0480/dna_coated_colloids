import numpy as np
import matplotlib.pyplot as plt

# Initialize lists for storing potential values and errors from each file
potential_values_all = []

# Loop over the 20 input files
for i in range(1, 21):
    file_path = f'distance_sur_A_BDon_10k_2nd_{i}.txt'  # Updated file path
    positions = np.loadtxt(file_path)

    start_distance = 5.3
    end_distance = 7
    filtered_positions = positions[(positions >= start_distance) & (positions <= end_distance)]

    num_bins = 20
    hist, bin_edges = np.histogram(filtered_positions, bins=num_bins)
    bin_midpoints = (bin_edges[1:] + bin_edges[:-1]) / 2
    last_bin_midpoint = bin_midpoints[-1]
    rmax_squared = last_bin_midpoint ** 2
    potential_values = -np.log(hist / bin_midpoints ** 2) + np.log(hist[-1] / rmax_squared) - 118
    #potential_values = -np.log(hist / bin_midpoints ** 2) + np.log(hist[-1] / rmax_squared)

    potential_values_all.append(potential_values)

# Convert potential_values_all to a numpy array
potential_values_all = np.array(potential_values_all)

# Calculate the error bars
error_bars = np.sqrt(np.mean(potential_values_all ** 2, axis=0) - np.mean(potential_values_all, axis=0) ** 2) / np.sqrt(20)

# Create a figure with a specific size
plt.figure(figsize=(10, 8))

# Calculate the mean of potential values
mean_potential = np.mean(potential_values_all, axis=0)

# Plot the shaded error region
plt.fill_between(bin_midpoints, mean_potential - error_bars, mean_potential + error_bars, color='black', alpha=0.5, label='Error in kinetic')

# Plot the mean values
plt.plot(bin_midpoints, mean_potential, color='black', label='Kinetic', linestyle='None', marker='o', markersize=20)

# Load data from the second file
fpi_data = np.loadtxt("FPIenergy_surA_-10-2.txt")

# Plot FPI data
plt.plot(fpi_data[:, 0], fpi_data[:, 1], color='red', label='Thermodynamic', linestyle='-.', marker='s', markersize=20)

# Add labels and title with larger font sizes
plt.xlabel(r'$r_\mathrm{AS}$', fontsize=52)
plt.ylabel(r'$\beta F$', fontsize=52)
#plt.title('Potential and FPI Data Comparison', fontsize=20)
plt.legend(fontsize=22)

# Set larger tick label sizes
plt.tick_params(axis='both', labelsize=36)

# Adjust layout to prevent clipping of tick-labels
plt.tight_layout()

# Save the plot as an image (e.g., PNG)
plt.savefig('combined_plot_sur_A_BDon_FPI.png')

# Print a message indicating the save location
print('Combined plot saved as combined_plot_sur_A_BDon_FPI.png')

# Save the potential values and error bars to a text file
output_file_path = 'energy_mean_and_error_bars_sur_A_BDon_10k_2nd_test.txt'
with open(output_file_path, 'w') as output_file:
    for i in range(len(bin_midpoints)):
        output_file.write(f'{bin_midpoints[i]} {np.mean(potential_values_all, axis=0)[i]} {error_bars[i]}\n')

# Print a message indicating the save location for the text file
print(f'Distance, Mean_V_BI values, and errors saved to {output_file_path}')

# Show the plot (optional)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SMALL_SIZE = 25
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rcParams["font.family"] = "Arial"

# List of input files (replace with your actual filenames)
files = [
    'long_Gillespie_BD_negno_test18_AB_DGT0_DG_-20_rho_ideal_1.0e-05_1000_timestep_vs_atom_counts.txt',
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-04_1000_modified_timestep_vs_atom_counts.txt',
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_modified_timestep_vs_atom_counts.txt',
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na200_modified.dump',
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na200_modified_timestep_vs_atom_counts.txt',
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_modified_timestep_vs_atom_counts.txt',
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-06_1000_modified_timestep_vs_atom_counts.txt', 
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-07_1000_modified_timestep_vs_atom_counts.txt', 
    #'/media/dell/Expansion/Pritam_sir_work/28_june/long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-08_1000_modified_timestep_vs_atom_counts.txt', 
    'long_Gillespie_BD_negno_test18_AB_DGT-4_DG_-20_rho_ideal_1.0e-05_1000_timestep_vs_atom_counts.txt'
]

# Corresponding labels
labels = [
    r'$\Delta G_T=0$',
    #r'$\rho_0 = DGT{-2}$',
    #r'$\rho_0 = DGT{-2}$',
    #r'$\rho_0 = DGT{-2}$',
    #r'$\rho_0 = DGT{-2}$',
    r'$\Delta G_T=-2$',
    #r'$\rho_0 = DGT{-2}$',
    #r'$\rho_0 = DGT{-2}$',
    #r'$\rho_0 = DGT{-2}$',
    r'$\Delta G_T=-4$'
]

# Plot styles
markers = ['o', 's', '^']
line_styles = ['-', '--', '-.', ':']
colors = ['blue', 'red', 'green']

# Create a plot
plt.figure(figsize=(12, 8))

# Loop over each file and plot
for i, (file, label) in enumerate(zip(files, labels)):
    # Read data
    data = pd.read_csv(file, sep='\t')
    
    # Extract columns
    timesteps = data['Timestep'].to_numpy()
    atom_counts = data['Number of Atoms'].to_numpy()
    
    # Plot the data
    plt.plot(
        timesteps,
        atom_counts,
        linestyle=line_styles[i % len(line_styles)],
        color=colors[i % len(colors)],
        label=label
    )

# Axis labels
plt.xlabel('$t$')
plt.ylabel('$n$')
plt.xlim(0, 5e6)

# Scientific notation on x-axis
ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.xaxis.get_offset_text().set_size(SMALL_SIZE)
ax.xaxis.get_offset_text().set_fontfamily("Arial")

# Legend
plt.legend(loc='upper left', ncol=3, handlelength=2.5, handletextpad=0.7, fontsize=20, labelspacing=0.1, columnspacing=1.0)

plt.grid(False)

# Save and show
output_file = 'number_of_atoms_vs_timestep_comparison_10curves.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.show()

print(f"Plot saved as {output_file}")


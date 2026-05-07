import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SMALL_SIZE = 30
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rcParams["font.family"] = "Arial"

# List of input files and their labels
files = [
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na50_modified.dump',
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na100_modified.dump',
   # 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_Na200_modified.dump'
]

# Corresponding labels for each file
labels = [
    r'$N_L = 50$',
    r'$N_L = 100$'
]

# Plot styles
line_styles = ['-', '--']
colors = ['blue', 'red']
line_widths = [2.5, 2.5]  # Optional: consistent line thickness

def parse_custom_file(file_path):
    timesteps = []
    atom_counts = []

    with open(file_path, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            if lines[i].strip() == "ITEM: TIMESTEP":
                try:
                    timestep = int(float(lines[i + 1].strip()))
                    i += 2
                except ValueError:
                    i += 1
                    continue
            elif lines[i].strip() == "ITEM: NUMBER OF ATOMS":
                try:
                    atom_count = int(float(lines[i + 1].strip()))
                    timesteps.append(timestep)
                    atom_counts.append(atom_count)
                    i += 2
                except ValueError:
                    i += 1
                    continue
            else:
                i += 1

    return timesteps, atom_counts

# Create a plot
plt.figure(figsize=(10, 7))

# Loop over each file and plot
for i, (file, label) in enumerate(zip(files, labels)):
    timesteps, atom_counts = parse_custom_file(file)
    
    plt.plot(
        timesteps,
        atom_counts,
        linestyle=line_styles[i % len(line_styles)],
        color=colors[i % len(colors)],
        linewidth=line_widths[i],
        label=label
    )

# Add labels, title, and legend
plt.xlabel(r'$\mathit{t}$', fontsize=40)
plt.ylabel(r'$\mathit{n}$', fontsize=40)
plt.xlim(0, 3e6)
plt.ylim(0, 110)

# Format x-axis labels as scientific notation
ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.xaxis.get_offset_text().set_size(SMALL_SIZE)
ax.xaxis.get_offset_text().set_fontfamily("Arial")

# Customize legend
plt.legend(loc='lower right', fontsize=18, framealpha=1, edgecolor='black', bbox_to_anchor=(1, 0))

plt.grid(False)

# Save the plot
output_file = 'number_of_atoms_vs_timestep_comparison_NA_mod.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.show()

print(f"Plot saved as {output_file}")


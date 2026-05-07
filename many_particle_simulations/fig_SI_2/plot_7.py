import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SMALL_SIZE = 30
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rcParams["font.family"] = "Arial"

files = [
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-08_1000_modified_timestep_vs_atom_counts.txt',
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-07_1000_modified_timestep_vs_atom_counts.txt',
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-06_1000_modified_timestep_vs_atom_counts.txt',
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_modified_timestep_vs_atom_counts.txt',
    'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-04_1000_modified_timestep_vs_atom_counts.txt'
]

labels = [
    r'$\rho_{\mathrm{id}} = 10^{-8}$',
    r'$\rho_{\mathrm{id}} = 10^{-7}$',
    r'$\rho_{\mathrm{id}} = 10^{-6}$',
    r'$\rho_{\mathrm{id}} = 10^{-5}$',
    r'$\rho_{\mathrm{id}} = 10^{-4}$'
]

markers = ['o', 's', '^', 'x', '*']
line_styles = ['-', '--', '-.', ':', '--']
colors = ['blue', 'red', 'green', 'purple', 'orange']

plt.figure(figsize=(12, 8))

for i, (file, label) in enumerate(zip(files, labels)):
    data = pd.read_csv(file, sep='\t')
    timesteps = data['Timestep'].to_numpy()
    atom_counts = data['Number of Atoms'].to_numpy()
    
    plt.plot(
        timesteps,
        atom_counts,
        linestyle=line_styles[i % len(line_styles)],
        color=colors[i % len(colors)],
        label=label
    )

# Add labels, title, and legend
plt.xlabel(r'$\mathit{t}$', fontsize=40)
plt.ylabel(r'$\mathit{n}$', fontsize=40)

plt.xlim(0, 5e6)

ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.xaxis.get_offset_text().set_size(SMALL_SIZE)
ax.xaxis.get_offset_text().set_fontfamily("Arial")

# --- ONLY legend size change ---
plt.legend(
    loc='upper right',
    ncol=2,
    handlelength=2.5,
    handletextpad=0.7,
    fontsize=18,   # changed from 20 → 18
    labelspacing=0.1,
    columnspacing=1.0,
    framealpha=1,
    edgecolor='black'
)
# ------------------------------

plt.grid(False)

output_file = 'number_of_atoms_vs_timestep_comparison_3.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')

plt.show()

print(f"Plot saved as {output_file}")

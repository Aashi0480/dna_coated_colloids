import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# =======================
# Plot styling
# =======================
SMALL_SIZE = 30
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rcParams["font.family"] = "Arial"

# =======================
# tau0 values and files
# =======================
tau0_values = [1, 10, 100]

files = {
    1:   'long_Gillespie_BD_tau0_1.dump',
    10:  'long_Gillespie_BD_tau0_10.dump',
    100: 'long_Gillespie_BD_tau0_100.dump'
}

labels = {
    1:   r'$\tau_0 = 1$',
    10:  r'$\tau_0 = 10$',
    100: r'$\tau_0 = 100$'
}

# =======================
# Plot styles
# =======================
line_styles = ['-', '--', '-.']
colors = ['blue', 'red', 'green']
line_width = 2.5

# =======================
# Dump file parser
# =======================
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

            elif lines[i].strip() == "ITEM: NUMBER OF ATOMS":
                try:
                    atom_count = int(float(lines[i + 1].strip()))
                    timesteps.append(timestep)
                    atom_counts.append(atom_count)
                    i += 2
                except ValueError:
                    i += 1

            else:
                i += 1

    return timesteps, atom_counts

# =======================
# Plot
# =======================
plt.figure(figsize=(10, 7))

for i, tau0 in enumerate(tau0_values):
    timesteps, atom_counts = parse_custom_file(files[tau0])

    plt.plot(
        timesteps,
        atom_counts,
        linestyle=line_styles[i],
        color=colors[i],
        linewidth=line_width,
        label=labels[tau0]
    )

# =======================
# Axes formatting
# =======================
plt.xlabel(r'$\mathit{t}$', fontsize=40)
plt.ylabel(r'$\mathit{n(t)}$', fontsize=40)

plt.xlim(0, 3e6)
plt.ylim(0, 160)

ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.xaxis.get_offset_text().set_size(SMALL_SIZE)
ax.xaxis.get_offset_text().set_fontfamily("Arial")

plt.legend(loc='lower right', fontsize=20, framealpha=1, edgecolor='black')

plt.grid(False)

# =======================
# Save figure
# =======================
output_file = 'n_vs_t_tau0_1_10_100.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.show()

print(f"Plot saved as {output_file}")


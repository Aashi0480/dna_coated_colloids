import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Plot settings
SMALL_SIZE = 25
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rcParams["font.family"] = "Arial"

# File
file_path = 'averages_output.txt'

# Initialize figure
fig, ax = plt.subplots(figsize=(8, 6))

# Load data
data = np.loadtxt(file_path, skiprows=1)
x = data[:, 0] / 1e6  # rescale timestep to millions
y1 = data[:, 1]       # Average_sur_2sbrg
y2 = data[:, 2]       # Average_sur_3sbrg

# Plot both series
ax.plot(x, y1, marker='None', markersize=6, linewidth=1.5, label=r'$n_{ps}^{(2)}$',alpha=1)
ax.plot(x, y2, marker='None', markersize=6, linewidth=1.5, label=r'$n_{ps}^{(3)}$',alpha=0.5)

# Customize main plot
ax.set_xlabel(r'$t$')
ax.set_ylabel('$n_\mathrm{ps}$')
ax.set_xlim(0, 5.0)   # explicitly set axis limit to 0–5

# Annotation for units
ax.text(0.63, -0.15, r'$(\times 10^6)$', fontsize=20,
        color='black', ha='center', va='center', transform=ax.transAxes)

# Inset: zoom in on early timesteps
inset_ax = inset_axes(ax, width="30%", height="30%", loc='upper center')
inset_ax.plot(x, y1, linewidth=1.5, alpha=1)
inset_ax.plot(x, y2, linewidth=1.5, alpha=0.7)
inset_ax.set_xlim(0, 0.1)   # zoom on first 0.1 million timesteps
inset_ax.set_ylim(0, 20)    # fixed y-limits

# Remove inset labels
inset_ax.set_xlabel(None)
inset_ax.set_ylabel(None)

# Axis styling
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.legend(fontsize=20, loc='best')

# Final layout and save
plt.tight_layout()
plt.savefig('plot_output_with_inset.png', dpi=300, bbox_inches='tight')
plt.show()


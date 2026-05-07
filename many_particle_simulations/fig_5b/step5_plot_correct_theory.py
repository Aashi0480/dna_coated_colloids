import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Plot settings
LARGE_SIZE = 25
MEDIUM_SIZE = 22
SMALL_SIZE = 20
plt.rc('font', size=LARGE_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rcParams["font.family"] = "Arial"

# Function to read the first dataset
def plot_first_data(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.dropna(subset=['avg_z', 'err_z', 'avg_a1b1_a2b2_sur_2sbrg', 'err_a1b1_a2b2_sur_2sbrg'])

    if isinstance(df.columns, pd.MultiIndex):
        df.columns = ['_'.join(col).strip() for col in df.columns.values]

    avg_z1 = df['avg_z'].to_numpy()
    err_z1 = df['err_z'].to_numpy()
    avg_l1 = df['avg_a1b1_a2b2_sur_2sbrg'].to_numpy()
    err_l1 = df['err_a1b1_a2b2_sur_2sbrg'].to_numpy()

    return avg_z1, err_z1, avg_l1, err_l1


# Function to read the second dataset
def plot_second_data(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.dropna(subset=['avg_z', 'err_z', 'avg_a1a2_a2a1', 'err_a1a2_a2a1'])

    if isinstance(df.columns, pd.MultiIndex):
        df.columns = ['_'.join(col).strip() for col in df.columns.values]

    avg_z2 = df['avg_z'].to_numpy()
    err_z2 = df['err_z'].to_numpy()
    avg_l2 = df['avg_a1a2_a2a1'].to_numpy()
    err_l2 = df['err_a1a2_a2a1'].to_numpy()

    return avg_z2, err_z2, avg_l2, err_l2


# Function to read theory data
def read_theory_data(theory_file):
    with open(theory_file, 'r') as file:
        code = file.read()
    local_vars = {}
    exec(code, {'np': np}, local_vars)
    return local_vars['n_bridges_3']


# Combined plotting function
def read_and_plot_data(file_path1, file_path2, theory_file):

    avg_z1, err_z1, avg_l1, err_l1 = plot_first_data(file_path1)
    avg_z2, err_z2, avg_l2, err_l2 = plot_second_data(file_path2)
    theory_values = read_theory_data(theory_file)

    theory_x = avg_z1[:3]

    # First dataset
    plt.errorbar(avg_z1, avg_l1,
                 xerr=err_z1, yerr=err_l1,
                 fmt='o',
                 ecolor='blue',
                 capsize=5,
                 markersize=12,
                 elinewidth=2,
                 label=r'$\langle n_{a_1;b_1}+n_{a_2;b_2}\rangle$',
                 color='blue')

    # Second dataset
    plt.errorbar(avg_z2, avg_l2,
                 xerr=err_z2, yerr=err_l2,
                 fmt='s',
                 ecolor='magenta',
                 capsize=5,
                 markersize=12,
                 elinewidth=2,
                 label=r'$\langle n_{a_1;a_2}+n_{a_2;a_1}\rangle$',
                 color='magenta')

    # Theory
    plt.plot(theory_x, theory_values,
             'D',
             color='green',
             markersize=12,
             label='Theory')

    plt.xlabel(r'$r_S$')
    plt.ylabel(r'$n^{(2)}_\mathrm{pp+ps}$')

    plt.xlim(0, 60)
    plt.ylim(0, 40)

    plt.legend(loc='upper right',
               bbox_to_anchor=(1, 1),
               fontsize=18,
               frameon=True)

    plt.tight_layout()

    plt.savefig(
        'na1a2_na2a1_na1b1_na2b2_with_theory.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


# run
read_and_plot_data(
    "na1b1_na2b2_sum_averages.txt",
    "na1a2_na2a1_sum_averages.txt",
    "theory.txt"
)

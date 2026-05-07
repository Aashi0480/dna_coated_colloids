import pandas as pd
import matplotlib.pyplot as plt

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
    df = df.dropna(subset=['avg_z', 'err_z', 'avg_l', 'err_l'])

    if isinstance(df.columns, pd.MultiIndex):
        df.columns = ['_'.join(col).strip() for col in df.columns.values]

    avg_z1 = df['avg_z'].to_numpy()
    err_z1 = df['err_z'].to_numpy()
    avg_l1 = df['avg_l'].to_numpy()
    err_l1 = df['err_l'].to_numpy()

    return avg_z1, err_z1, avg_l1, err_l1


# Function to read the second dataset
def plot_second_data(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.dropna(subset=['avg_z', 'err_z',
                           'avg_na1_na2_nb1_nb2',
                           'err_na1_na2_nb1_nb2'])

    if isinstance(df.columns, pd.MultiIndex):
        df.columns = ['_'.join(col).strip() for col in df.columns.values]

    avg_z2 = df['avg_z'].to_numpy()
    err_z2 = df['err_z'].to_numpy()
    avg_l2 = df['avg_na1_na2_nb1_nb2'].to_numpy()
    err_l2 = df['err_na1_na2_nb1_nb2'].to_numpy()

    return avg_z2, err_z2, avg_l2, err_l2


# Combined plotting function
def read_and_plot_data(file_path1, file_path2):

    avg_z1, err_z1, avg_l1, err_l1 = plot_first_data(file_path1)
    avg_z2, err_z2, avg_l2, err_l2 = plot_second_data(file_path2)

    # First dataset
    plt.errorbar(avg_z1, avg_l1,
                 xerr=err_z1,
                 yerr=err_l1,
                 fmt='o',
                 ecolor='blue',
                 capsize=6,
                 markersize=12,
                 elinewidth=2,
                 label=r'$\langle n^{(l)}\rangle$',
                 color='blue')

    # Second dataset
    plt.errorbar(avg_z2, avg_l2,
                 xerr=err_z2,
                 yerr=err_l2,
                 fmt='s',
                 ecolor='magenta',
                 capsize=6,
                 markersize=12,
                 elinewidth=2,
                 label=r'$\langle n^{(f)}\rangle$',
                 color='magenta')

    plt.xlabel(r'$r_S$')
    plt.ylabel(r'$n^{(l)}, n^{(f)}$')

    plt.xticks(range(0, 90, 10))
    plt.yticks(range(0, 90, 10))

    plt.xlim(0, 60)
    plt.ylim(0, 80)

    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig(
        'combined_plot_correct_DGT_-4.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


# run
read_and_plot_data(
    "l_sum_averages.txt",
    "na1_na2_nb1_nb2_sum_averages.txt"
)

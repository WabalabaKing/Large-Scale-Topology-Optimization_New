import pandas as pd
import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties
#plt.style.use('seaborn-v0_8-deep')
#plt.style.use('seaborn-v0_8-dark-palette')
plt.style.use('petroff10')
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"

SCATTER_STYLE = dict(
    s=28,
    facecolors="black",
    edgecolors="gray",
    linewidths=1.2,
    marker="P"
)

# -------------------------------------------------------------------
# Shared line styling taken from plot_optimality and reused everywhere
# -------------------------------------------------------------------
LINE_COLORS = ['C0', 'C2', 'C1']
LINESTYLES = ["-", "dashdot"]

"""
Generate publication-quality convergence plots from IPOPT optimization logs.

This script reads one or more CSV files containing per-iteration IPOPT
optimization history and produces a consistent set of convergence figures
using a shared plotting style. The generated figures include:

    - Objective function versus iteration
    - Dual residual (optimality) versus iteration
    - Primal residual (feasibility) versus iteration
    - Newton step perturbation norm versus iteration
    - (Optional) Bar chart showing the number of free-mu iterations

Each curve is drawn using a common color and linestyle convention to ensure
consistent presentation across figures. Iterations for which the optimizer
is not operating in free-mu mode are highlighted using scatter markers,
allowing transitions between globalization modes to be identified visually.

The script is intended for producing publication-quality figures for
performance comparisons of different IPOPT barrier parameter update
strategies and KKT-error globalization settings.
"""


def get_plot_styles(n_curves):
    """
    Return line colors and linestyles for n_curves, reusing the same
    styling pattern defined for the optimality plot.
    """
    colors = [LINE_COLORS[i % len(LINE_COLORS)] for i in range(n_curves)]
    linestyles = [LINESTYLES[i % len(LINESTYLES)] for i in range(n_curves)]
    return colors, linestyles


def scatter_false_points(ax, x, y, mask):
    """Overlay markers at points where mask is True (False rows in CSV)."""
    if mask.any():
        ax.scatter(x[mask], y[mask], **SCATTER_STYLE)


def last_col_bool_mask(df):
    col = df.columns[-1]
    s = df[col]
    if s.dtype == bool:
        return s.values
    # Normalize common truthy strings/numbers
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "t", "yes"}).values


def apply_plot_style(ax, xlabel, ylabel, show_legend=True):
    """Apply the consistent perturbation-style formatting + legend style."""

    plt.xlabel(xlabel, fontsize=22, fontname="Times New Roman")
    plt.ylabel(ylabel, fontsize=22, fontname="Times New Roman")

    # Grid: major + minor
    ax.grid(which='major', color='black', linestyle=':', linewidth=0.01)
    ax.minorticks_on()
    ax.grid(which='minor', color='black', linestyle=':', linewidth=0.01)

    # Spine thickness
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)

    # Tick fonts
    plt.xticks(fontname="Times New Roman", fontsize=20)
    plt.yticks(fontname="Times New Roman", fontsize=20)

    # Tick params
    ax.tick_params(bottom=True, top=True, left=True, right=True)
    ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax.tick_params(which='major', length=5, width=1.2, direction='in')
    ax.tick_params(which='minor', length=5, width=1.2, direction='in')

    plt.xlim(-10, 350)

    # Shared legend formatting
    if show_legend:
        #legend_font = FontProperties(family='Times New Roman', size=22)
        #ax.legend(frameon=False, ncol=1, prop=legend_font)

        # ⬇ Use monospace font so columns align vertically
        ax.legend(
            prop={'family': 'monospace', 'size': 20},
            labelspacing=2.0,
            handlelength=3.5,
            frameon=False, ncol=1
        )


def plot_objective(input_files):
    # labels = [r"$k_{hist}:2,  \mathcal{C}: 12.249$", r"$k_{hist}:6,  \mathcal{C}: 12.546$", r"$k_{hist}:12, \mathcal{C}: 12.507 $"]

    labels = [
        "l:$2$",
        "l:$4$"
    ]
    output_file = "Plots/PROBING/objective_vs_iter_multi.png"

    colors, linestyles = get_plot_styles(len(input_files))

    fig, ax = plt.subplots()
    for file, label, ls, color in zip(input_files, labels, linestyles, colors):
        df = pd.read_csv(file)
        x, y = df["iter"], df["objective"]
        ax.plot(x, y, linewidth=2.0, label=label, linestyle=ls, color=color)
        mask = (last_col_bool_mask(df) == False)
        scatter_false_points(ax, x, y, mask)

    apply_plot_style(ax, "Iteration", r"$\mathcal{C}$ (J)")

    plt.ylim(-40, 800)
    

    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.0, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Combined plot saved to {output_file}")


def plot_optimality(input_files):
    #labels = [r"$k_{hist}:2$", r"$k_{hist}:6$", r"$k_{hist}:12$"]
    output_file = "Plots/PROBING/optimality_vs_iter_multi.png"

    colors, linestyles = get_plot_styles(len(input_files))

    fig, ax = plt.subplots()
    for file, ls, color in zip(input_files, linestyles, colors):
        df = pd.read_csv(file)
        x, y = df["iter"], df["inf_du"]
        ax.semilogy(x, y, linewidth=2.0, linestyle=ls, color=color)
        mask = (last_col_bool_mask(df) == False)
        scatter_false_points(ax, x, y, mask)

    apply_plot_style(ax, "Iteration", r"$\mathcal{R}_{du}$")

    plt.ylim(1E-08, 1E3)

    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.0, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Combined plot saved to {output_file}")


def plot_feasibility(input_files):
   # labels = [r"$k_{hist}:2$", r"$k_{hist}:6$", r"$k_{hist}:12$"]
    output_file = "Plots/PROBING/feasibility_vs_iter_multi.png"

    colors, linestyles = get_plot_styles(len(input_files))

    fig, ax = plt.subplots()
    for file, ls, color in zip(input_files, linestyles, colors):
        df = pd.read_csv(file)
        x, y = df["iter"], df["inf_pr"]
        ax.semilogy(x, y, linewidth=2.0, linestyle=ls, color=color)
        mask = (last_col_bool_mask(df) == False)

        scatter_false_points(ax, x, y, mask)

    apply_plot_style(ax, "Iteration", r"$\mathcal{R}_{pr}$")

    plt.ylim(1E-8, 1E3)

    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.0, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Combined plot saved to {output_file}")


def plot_perturbation(input_files):
   # labels = [r"$k_{hist}:2$", r"$k_{hist}:6$", r"$k_{hist}:12$"]
    output_file = "Plots/PROBING/perturbation_vs_iter_multi.png"

    colors, linestyles = get_plot_styles(len(input_files))

    fig, ax = plt.subplots()
    for file, ls, color in zip(input_files, linestyles, colors):
        df = pd.read_csv(file)
        x, y = df["iter"], df["||d||"]
        ax.semilogy(x, y, linewidth=2.0, linestyle=ls, color=color)
        mask = (last_col_bool_mask(df) == False)

        scatter_false_points(ax, x, y, mask)

    apply_plot_style(ax, "Iteration", r"$||d||$")

    plt.ylim(1E-08, 1E2)

    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.0, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Combined plot saved to {output_file}")


def plot_free_mode_true_counts(input_files, labels=None,
                               output_file="Plots/AUTO_SCALING_0.001/free_mode_true_counts.png",
                               bar_width=0.20, bar_spacing=0.25):
    """
    Reads 'free_mu_mode' column from each CSV and plots only the count of True entries
    across configurations, allowing manual bar width and spacing.
    """
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if labels is None:
        labels = ["l:$2$", "$l:$6$", r"$k_{hist}:12$"]

    # Match your shared line plot colors
    bar_colors = [LINE_COLORS[i % len(LINE_COLORS)] for i in range(len(labels))]

    # Compute True counts
    true_counts = []
    for file, label in zip(input_files, labels):
        df = pd.read_csv(file)
        if "free_mu_mode" not in df.columns:
            raise KeyError(f"'free_mu_mode' not found in {file}")
        true_counts.append((df["free_mu_mode"] == True).sum())

    print("True counts:", dict(zip(labels, true_counts)))

    # === Bar positioning with spacing control ===
    n = len(labels)
    total_bar_and_space = bar_width + bar_spacing
    x = [i * total_bar_and_space for i in range(n)]

    # --- Plot ---
    fig, ax = plt.subplots()
    bars = ax.bar(x, true_counts, width=bar_width,
                  color=bar_colors[:n], edgecolor="black", linewidth=1.2)

    # Set x-tick positions and labels
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontname="Times New Roman", fontsize=20)

    ax.set_ylabel("Count of True (Free Mode)", fontsize=22, fontname="Times New Roman")
    ax.set_xlabel("Configuration", fontsize=22, fontname="Times New Roman")

    # Grid and style (consistent with your other plots)
    ax.grid(which='major', axis='y', color='black', linestyle=':', linewidth=0.01)
    ax.minorticks_on()
    ax.grid(which='minor', axis='y', color='black', linestyle=':', linewidth=0.01)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params(which='major', length=10, width=1.2, direction='in')
    ax.tick_params(which='minor', length=10, width=1.2, direction='in')
    plt.yticks(fontname="Times New Roman", fontsize=20)

    # Annotate bar tops
    for bar, val in zip(bars, true_counts):
        ax.annotate(f"{val}", xy=(bar.get_x() + bar.get_width() / 2, val),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=14, fontname="Times New Roman")

    # Final formatting
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.0, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Free-mode True count bar chart saved to {output_file}")


# === Example usage ===
input_files = [
    "/media/prateek/MDO/calTop/TestCases/SCB/NEW/100MM/R_019/CGx_8/PROBING/KKT_ERR_2/CGx_8_PROBING_KKT_RED_2.csv",
#    "/media/prateek/MDO/calTop/TestCases/SCB/NEW/100MM/R_019/CGx_8/LOQO/KKT_ERR_6/CGx_8_LOQO_KKT_RED_6.csv",
    "/media/prateek/MDO/calTop/TestCases/SCB/NEW/100MM/R_019/CGx_8/PROBING/KKT_ERR_12/CGx_8_PROBING_KKT_RED_12.csv"
]

plot_perturbation(input_files)
plot_objective(input_files)
plot_optimality(input_files)
plot_feasibility(input_files)

#plot_free_mode_true_counts(input_files, labels=None,
#                               output_file="Plots/LOQO/free_mode_true_counts.png",
#                               bar_width=0.2, bar_spacing=0.20)
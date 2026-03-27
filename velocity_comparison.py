import os
import numpy as np
import matplotlib.pyplot as plt

def plot_avg_velocity_all_alpha():
    """
    Plot <velocity outside wound>(t) vs time for multiple alpha values
    on a single publication-quality figure.

    Expected data files (2 columns: time[min], <v_outside_wound>(t)):
        avg_velocity_vs_time_alpha_10.dat
        avg_velocity_vs_time_alpha_35.dat
        avg_velocity_vs_time_alpha_60.dat
        avg_velocity_vs_time_alpha_120.dat
        avg_velocity_vs_time_alpha_200.dat

    Only every 'stride'-th data point is plotted to avoid overcrowding.
    Output:
        avg_velocity_vs_time_allalpha.eps
    """

    # ------------------------ USER PARAMETERS ------------------------
    alpha_list = [10, 35, 60, 120, 200]

    stride = 10   # plot every 10th point (change to 5 if needed)

    lw_line = 2.5
    ms_mark = 6

    colors = [
        (0.0,   0.0,   0.0),    # black      α=10
        (0.0,   0.447, 0.741),  # blue       α=35
        (0.85,  0.325, 0.098),  # orange     α=60
        (0.466, 0.674, 0.188),  # green      α=120
        (0.494, 0.184, 0.556),  # purple     α=200
    ]
    linestyles = ['-', '--', '-.', ':', '-']

    # ------------------------ PLOT STYLING ------------------------
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 16,
        "axes.labelsize": 22,
        "axes.titlesize": 22,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 10,
        "axes.linewidth": 1.6,
        "xtick.major.width": 1.6,
        "ytick.major.width": 1.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.visible": False,
        "ytick.minor.visible": False,
        "figure.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    # ------------------------ LOAD & PLOT ------------------------
    for alpha, color, ls in zip(alpha_list, colors, linestyles):
        fname = f"avg_velocity_vs_time_alpha_{alpha}.dat"
        if not os.path.exists(fname):
            print(f"Warning: file not found: {fname}")
            continue

        data = np.loadtxt(fname)
        if data.ndim == 1:
            data = data[np.newaxis, :]

        if data.shape[1] < 2:
            print(f"Warning: {fname} does not have 2 columns; skipped.")
            continue

        t = data[:, 0]
        v = data[:, 1]

        mask = np.isfinite(t) & np.isfinite(v)
        t = t[mask]
        v = v[mask]

        if t.size == 0:
            print(f"No valid data in {fname}; skipped.")
            continue

        t_plot = t[::stride]
        v_plot = v[::stride]

        ax.plot(
            t_plot, v_plot,
            linestyle=ls,
            linewidth=lw_line,
            color=color,
            marker='o',
            markersize=ms_mark,
            markerfacecolor=color,
            markeredgecolor='k',
            label=fr"$\alpha = {alpha}$"
        )

    # ------------------------ AXES & LABELS ------------------------
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(r"$\langle v_{\mathrm{outside\ wound}}(t) \rangle$  [$\mu$m/min]")

    # === NEW: Fixed y-axis ===
    ax.set_ylim(0, 1.1)

    # === NEW: Remove all gridlines ===
    ax.grid(False)

    ax.tick_params(which="both", direction="in")
    ax.legend(frameon=False, loc="best")

    fig.tight_layout()

    # ------------------------ SAVE EPS ------------------------
    out_eps = "avg_velocity_vs_time_allalpha.eps"
    fig.savefig(out_eps, format="eps", bbox_inches="tight")
    print(f"Saved {out_eps}")

    plt.close(fig)

if __name__ == "__main__":
    plot_avg_velocity_all_alpha()


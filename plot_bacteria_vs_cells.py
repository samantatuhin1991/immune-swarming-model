import os
import numpy as np
import matplotlib.pyplot as plt

def plot_bacteria_vs_rate_immune_allalpha():
    """
    Plot bacteria vs max rate and bacteria vs immune cells
    for three different alpha values (10, 35, 80).

    Data files for each alpha (XX = 10, 35, 80):
        bacteria_vs_max_rate_sim_alpha_XX.dat    (sim:  [|m|, B, log10(B)])
        bacteria_vs_max_rate_anal_alpha_XX.dat   (anal: [|m|, B, log10(B)])
        bacteria_vs_immune_cells_sim_alpha_XX.dat  (sim:  [immune, B, log10(B)])
        bacteria_vs_immune_cells_anal_alpha_XX.dat (anal: [immune, B, log10(B)])

    Generates 4 publication-quality figures:
        1) bacteria_vs_maxrate_linear_allalpha.(png|pdf)
            y = B, x = |m|
        2) bacteria_vs_maxrate_log10_allalpha.(png|pdf)
            y = log10(B), x = |m|
        3) bacteria_vs_immune_linear_allalpha.(png|pdf)
            y = B, x = immune cells n_w at t=90
        4) bacteria_vs_immune_log10_allalpha.(png|pdf)
            y = log10(B), x = immune cells n_w at t=90
    """

    # ---------- Styling (publication-level) ----------
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times', 'DejaVu Serif', 'Computer Modern Roman']
    fs_ax   = 18   # axis tick font size
    fs_lab  = 22   # axis label font size
    fs_leg  = 16   # legend font size
    lw_ax   = 1.5  # axis line width
    lw_line = 2.5  # line width
    ms_sim  = 8    # marker size (simulation)
    ms_anal = 9    # marker size (analytical)
    dpi_out = 600  # PNG resolution

    # Alphas and associated colors (RGB)
    alpha_list = [10, 35, 60]
    colors = np.array([
        [0.0000, 0.4470, 0.7410],  # blue  (alpha=10)
        [0.8500, 0.3250, 0.0980],  # orange (alpha=35)
        [0.4660, 0.6740, 0.1880],  # green (alpha=80)
    ])

    # Utility: safe load that always returns 2D array
    def load_ascii_2d(filename):
        """Load ASCII file and ensure result is 2D: shape (n, m)."""
        data = np.loadtxt(filename)
        data = np.atleast_2d(data)
        return data

    # Utility: deduplicate legend entries
    def dedup_legend(ax):
        handles, labels = ax.get_legend_handles_labels()
        seen = {}
        new_handles = []
        new_labels = []
        for h, lab in zip(handles, labels):
            if lab not in seen:
                seen[lab] = True
                new_handles.append(h)
                new_labels.append(lab)
        ax.legend(new_handles, new_labels, loc='best',
                  frameon=False, fontsize=fs_leg)

    # =========================================================================
    # 1) BACTERIA VS MAX RATE  (linear y, y = B, x = |m|)
    # =========================================================================
    fig1, ax1 = plt.subplots(figsize=(9, 7))

    for k, alpha_val in enumerate(alpha_list):
        col = colors[k, :]

        # ----- Simulation data: bacteria vs max rate -----
        sim_file = f"bacteria_vs_max_rate_sim_alpha_{alpha_val}.dat"
        if os.path.exists(sim_file):
            sim_data = load_ascii_2d(sim_file)  # [|m|, B, log10(B)]
            if sim_data.shape[1] >= 2:
                m_sim  = sim_data[:, 0]
                B_sim  = sim_data[:, 1]
                good   = np.isfinite(m_sim) & np.isfinite(B_sim)
                m_sim  = m_sim[good]
                B_sim  = B_sim[good]

                if m_sim.size > 0:
                    ax1.plot(m_sim, B_sim, 'o',
                             markersize=ms_sim,
                             markerfacecolor=col,
                             markeredgecolor='k',
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Simulation, α = {alpha_val}")
        else:
            print(f"Warning: Simulation file not found: {sim_file}")

        # ----- Analytical data: bacteria vs max rate -----
        anal_file = f"bacteria_vs_max_rate_anal_alpha_{alpha_val}.dat"
        if os.path.exists(anal_file):
            anal_data = load_ascii_2d(anal_file)  # [|m|, B, log10(B)]
            if anal_data.shape[1] >= 2:
                m_anal = anal_data[:, 0]
                B_anal = anal_data[:, 1]
                goodA  = np.isfinite(m_anal) & np.isfinite(B_anal)
                m_anal = m_anal[goodA]
                B_anal = B_anal[goodA]

                if m_anal.size > 0:
                    ax1.plot(m_anal, B_anal, 's',
                             markersize=ms_anal,
                             markerfacecolor='w',
                             markeredgecolor=col,
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Analytic, α = {alpha_val}")
        else:
            print(f"Warning: Analytical file not found: {anal_file}")

    ax1.set_xlabel('Max. rate of recruitment', fontsize=fs_lab)
    ax1.set_ylabel('Bacteria amount', fontsize=fs_lab)
    ax1.tick_params(axis='both', labelsize=fs_ax)
    for spine in ax1.spines.values():
        spine.set_linewidth(lw_ax)
    ax1.tick_params(direction='out')
    ax1.grid(False)

    dedup_legend(ax1)

    fig1.tight_layout()
    fig1.savefig('bacteria_vs_maxrate_linear_allalpha.png', dpi=dpi_out, bbox_inches='tight')
    fig1.savefig('bacteria_vs_maxrate_linear_allalpha.pdf', dpi=dpi_out, bbox_inches='tight')
    # plt.close(fig1)  # keep open if you want to inspect

    # =========================================================================
    # 2) LOG10(Bacteria) VS MAX RATE  (data-converted log, y = log10 B)
    # =========================================================================
    fig2, ax2 = plt.subplots(figsize=(9, 7))

    for k, alpha_val in enumerate(alpha_list):
        col = colors[k, :]

        # ----- Simulation data -----
        sim_file = f"bacteria_vs_max_rate_sim_alpha_{alpha_val}.dat"
        if os.path.exists(sim_file):
            sim_data = load_ascii_2d(sim_file)  # [|m|, B, log10(B)]
            if sim_data.shape[1] >= 3:
                m_sim    = sim_data[:, 0]
                logB_sim = sim_data[:, 2]
                good     = np.isfinite(m_sim) & np.isfinite(logB_sim)
                m_sim    = m_sim[good]
                logB_sim = logB_sim[good]

                if m_sim.size > 0:
                    ax2.plot(m_sim, logB_sim, 'o',
                             markersize=ms_sim,
                             markerfacecolor=col,
                             markeredgecolor='k',
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Simulation, α = {alpha_val}")
        else:
            print(f"Warning: Simulation file not found: {sim_file}")

        # ----- Analytical data -----
        anal_file = f"bacteria_vs_max_rate_anal_alpha_{alpha_val}.dat"
        if os.path.exists(anal_file):
            anal_data = load_ascii_2d(anal_file)  # [|m|, B, log10(B)]
            if anal_data.shape[1] >= 3:
                m_anal    = anal_data[:, 0]
                logB_anal = anal_data[:, 2]
                goodA     = np.isfinite(m_anal) & np.isfinite(logB_anal)
                m_anal    = m_anal[goodA]
                logB_anal = logB_anal[goodA]

                if m_anal.size > 0:
                    ax2.plot(m_anal, logB_anal, 's',
                             markersize=ms_anal,
                             markerfacecolor='w',
                             markeredgecolor=col,
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Analytic, α = {alpha_val}")
        else:
            print(f"Warning: Analytical file not found: {anal_file}")

    ax2.set_xlabel('Max. rate of recruitment', fontsize=fs_lab)
    ax2.set_ylabel('log(Bacteria amount)', fontsize=fs_lab)
    ax2.tick_params(axis='both', labelsize=fs_ax)
    for spine in ax2.spines.values():
        spine.set_linewidth(lw_ax)
    ax2.tick_params(direction='out')
    ax2.grid(False)

    dedup_legend(ax2)

    fig2.tight_layout()
    fig2.savefig('bacteria_vs_maxrate_log10_allalpha.png', dpi=dpi_out, bbox_inches='tight')
    fig2.savefig('bacteria_vs_maxrate_log10_allalpha.pdf', dpi=dpi_out, bbox_inches='tight')
    # plt.close(fig2)

    # =========================================================================
    # 3) BACTERIA VS IMMUNE CELLS (n_w(90))  (linear y, y = B)
    # =========================================================================
    fig3, ax3 = plt.subplots(figsize=(9, 7))

    for k, alpha_val in enumerate(alpha_list):
        col = colors[k, :]

        # ----- Simulation data: bacteria vs immune -----
        sim_file = f"bacteria_vs_immune_cells_sim_alpha_{alpha_val}.dat"
        if os.path.exists(sim_file):
            sim_data = load_ascii_2d(sim_file)  # [immune, B, log10(B)]
            if sim_data.shape[1] >= 2:
                Imm_sim = sim_data[:, 0]
                B_sim   = sim_data[:, 1]
                good    = np.isfinite(Imm_sim) & np.isfinite(B_sim)
                Imm_sim = Imm_sim[good]
                B_sim   = B_sim[good]

                if Imm_sim.size > 0:
                    ax3.plot(Imm_sim, B_sim, 'o',
                             markersize=ms_sim,
                             markerfacecolor=col,
                             markeredgecolor='k',
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Simulation, α = {alpha_val}")
        else:
            print(f"Warning: Simulation immune file not found: {sim_file}")

        # ----- Analytical data: bacteria vs immune -----
        anal_file = f"bacteria_vs_immune_cells_anal_alpha_{alpha_val}.dat"
        if os.path.exists(anal_file):
            anal_data = load_ascii_2d(anal_file)  # [immune, B, log10(B)]
            if anal_data.shape[1] >= 2:
                Imm_anal = anal_data[:, 0]
                B_anal   = anal_data[:, 1]
                goodA    = np.isfinite(Imm_anal) & np.isfinite(B_anal)
                Imm_anal = Imm_anal[goodA]
                B_anal   = B_anal[goodA]

                if Imm_anal.size > 0:
                    ax3.plot(Imm_anal, B_anal, 's',
                             markersize=ms_anal,
                             markerfacecolor='w',
                             markeredgecolor=col,
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Analytic, α = {alpha_val}")
        else:
            print(f"Warning: Analytical immune file not found: {anal_file}")

    ax3.set_xlabel('Immune cells n_w', fontsize=fs_lab)
    ax3.set_ylabel('Bacteria amount', fontsize=fs_lab)
    ax3.tick_params(axis='both', labelsize=fs_ax)
    for spine in ax3.spines.values():
        spine.set_linewidth(lw_ax)
    ax3.tick_params(direction='out')
    ax3.grid(False)

    dedup_legend(ax3)

    fig3.tight_layout()
    fig3.savefig('bacteria_vs_immune_linear_allalpha.png', dpi=dpi_out, bbox_inches='tight')
    fig3.savefig('bacteria_vs_immune_linear_allalpha.pdf', dpi=dpi_out, bbox_inches='tight')
    # plt.close(fig3)

    # =========================================================================
    # 4) LOG10(Bacteria) VS IMMUNE CELLS (n_w(90)) (data log)
    # =========================================================================
    fig4, ax4 = plt.subplots(figsize=(9, 7))

    for k, alpha_val in enumerate(alpha_list):
        col = colors[k, :]

        # ----- Simulation data -----
        sim_file = f"bacteria_vs_immune_cells_sim_alpha_{alpha_val}.dat"
        if os.path.exists(sim_file):
            sim_data = load_ascii_2d(sim_file)  # [immune, B, log10(B)]
            if sim_data.shape[1] >= 3:
                Imm_sim  = sim_data[:, 0]
                logB_sim = sim_data[:, 2]
                good     = np.isfinite(Imm_sim) & np.isfinite(logB_sim)
                Imm_sim  = Imm_sim[good]
                logB_sim = logB_sim[good]

                if Imm_sim.size > 0:
                    ax4.plot(Imm_sim, logB_sim, 'o',
                             markersize=ms_sim,
                             markerfacecolor=col,
                             markeredgecolor='k',
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Simulation, α = {alpha_val}")
        else:
            print(f"Warning: Simulation immune file not found: {sim_file}")

        # ----- Analytical data -----
        anal_file = f"bacteria_vs_immune_cells_anal_alpha_{alpha_val}.dat"
        if os.path.exists(anal_file):
            anal_data = load_ascii_2d(anal_file)  # [immune, B, log10(B)]
            if anal_data.shape[1] >= 3:
                Imm_anal  = anal_data[:, 0]
                logB_anal = anal_data[:, 2]
                goodA     = np.isfinite(Imm_anal) & np.isfinite(logB_anal)
                Imm_anal  = Imm_anal[goodA]
                logB_anal = logB_anal[goodA]

                if Imm_anal.size > 0:
                    ax4.plot(Imm_anal, logB_anal, 's',
                             markersize=ms_anal,
                             markerfacecolor='w',
                             markeredgecolor=col,
                             linestyle='none',
                             linewidth=lw_line,
                             label=f"Analytic, α = {alpha_val}")
        else:
            print(f"Warning: Analytical immune file not found: {anal_file}")

    ax4.set_xlabel('Immune cells n_w', fontsize=fs_lab)
    ax4.set_ylabel('log(Bacteria amount)', fontsize=fs_lab)
    ax4.tick_params(axis='both', labelsize=fs_ax)
    for spine in ax4.spines.values():
        spine.set_linewidth(lw_ax)
    ax4.tick_params(direction='out')
    ax4.grid(False)

    dedup_legend(ax4)

    fig4.tight_layout()
    fig4.savefig('bacteria_vs_immune_log10_allalpha.png', dpi=dpi_out, bbox_inches='tight')
    fig4.savefig('bacteria_vs_immune_log10_allalpha.pdf', dpi=dpi_out, bbox_inches='tight')
    # plt.close(fig4)

    print('Generated 4 combined plots for alpha = 10,35,80 (simulation + analytical).')

if __name__ == "__main__":
    plot_bacteria_vs_rate_immune_allalpha()


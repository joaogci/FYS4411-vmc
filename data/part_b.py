from process_data import read_data_file
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

N_vals = [1, 10, 100, 500]
filenames = ["./data_files/data_N1_d3.csv", "./data_files/data_N10_d3.csv", "./data_files/data_N100_d3.csv", "./data_files/data_N500_d3.csv"]

for i, filename in enumerate(filenames):
    N = N_vals[i]
    
    (alphas, 
    runtime, 
    acceptance_ratio, 
    sampler_param, 
    steps, 
    mean_E, 
    std_E, 
    std_E_blocking) = read_data_file(filename)

    plt.figure(1, figsize=(3, 2.25))
    plt.plot(alphas, mean_E[-1, :] / N, ".-", label=rf"$N=${N}")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\langle E \rangle$ / N")
    plt.legend()

    plt.figure(2, figsize=(3, 2.25))
    plt.plot(alphas, std_E[-1, :], ".-", label=rf"$N=${N}")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\sigma$")
    plt.legend()
    
    if N == 100:
        for i, alpha in enumerate(alphas):
            print(f"{alpha} & {mean_E[-1, i]:.6f} & {std_E[-1, i]:.6f} & {std_E_blocking[-1, i]:.6f} & {runtime[i]:.3f} \\\\")
    
plt.figure(1, figsize=(3, 2.25))
plt.savefig("E_vs_alphas.eps", bbox_inches="tight")

plt.figure(2, figsize=(3, 2.25))
plt.savefig("std_E_vs_alphas.eps", bbox_inches="tight")

(alphas, 
runtime, 
acceptance_ratio, 
sampler_param, 
steps, 
mean_E, 
std_E, 
std_E_blocking) = read_data_file("./data_files/N100_d3_2.csv")
E_exact = 150

for a in [0.3, 0.42, 0.62, 0.7]:
    idx = np.where(alphas == a)[0][0]
    
    plt.figure(3, figsize=(3, 2.25))
    plt.plot(np.log2(steps[2:, idx]), np.abs(mean_E[2:, idx] - E_exact), label=rf"$\alpha=${alphas[idx]}")
    plt.xlabel("MC cycles")
    plt.ylabel(r"$|\langle E \rangle - E_{0}|$")
    # plt.legend(loc='lower right')
    # plt.xticks(np.log2(steps[::2, idx]), [r"$2^{11}$", r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])
    plt.legend()
    plt.xticks(np.log2(steps[2::2, idx]), [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])

    plt.figure(4, figsize=(3, 2.25))
    plt.plot(np.log2(steps[2:, idx]), std_E[2:, idx], label=rf"$\alpha=${alphas[idx]}")
    plt.xlabel("MC cycles")
    plt.ylabel(r"$\sigma$")
    # plt.legend(loc='upper right')
    # plt.xticks(np.log2(steps[::2, idx]), [r"$2^{11}$", r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])
    plt.legend()
    plt.xticks(np.log2(steps[2::2, idx]), [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])

plt.figure(3, figsize=(3, 2.25))
plt.savefig("E_vs_mc_cycles.eps", bbox_inches="tight")

plt.figure(4, figsize=(3, 2.25))
plt.savefig("std_E_vs_mc_cycles.eps", bbox_inches="tight")


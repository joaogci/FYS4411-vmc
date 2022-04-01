from process_data import read_data_file
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

E_exact = 150
N = 100
a = 0.66
sl_vals = [0.05, 0.1, 0.5, 1.0]
filenames_sl = ["./data_files/data_sl0.05.csv", "./data_files/data_sl0.1.csv", "./data_files/data_sl0.5.csv", "./data_files/data_sl1.0.csv"]

dt_vals = [0.001, 0.01, 0.1, 1.0, 10.0]
filenames_dt = ["./data_files/data_dt0.001.csv", "./data_files/data_dt0.01.csv", "./data_files/data_dt0.1.csv", "./data_files/data_dt1.0.csv", "./data_files/data_dt10.0.csv"]

sl = 0.5
filename_sl = "./data_files/data_sl0.5.csv"
dt = 0.1
filename_dt = "./data_files/data_dt0.1.csv"

(alphas, 
runtime, 
acceptance_ratio, 
sampler_param, 
steps, 
mean_E, 
std_E, 
std_E_blocking) = read_data_file(filename_sl)


plt.figure(1, figsize=(3, 2.25))
plt.plot(alphas, mean_E[-1, :] / N, ".-", label=rf"Metropolis $sl=${sl}")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\langle E \rangle$")

plt.figure(2, figsize=(3, 2.25))
plt.plot(alphas, std_E[-1, :], ".-", label=rf"Metropolis $sl=${sl}")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\sigma$")

plt.figure(3, figsize=(3, 2.25))
idx = np.where(alphas == a)[0][0]
plt.plot(np.log2(steps[:, idx]), np.abs(mean_E[:, idx] - E_exact), label=rf"Metropolis $sl=${sl}")
plt.xticks(np.log2(steps[::2, idx]), [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle E \rangle$")

plt.figure(4, figsize=(3, 2.25))
idx = np.where(alphas == a)[0][0]
plt.plot(np.log2(steps[:, idx]), std_E[:, idx], label=rf"Metropolis $sl=${sl}")
plt.xticks(np.log2(steps[::2, idx]), [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])
plt.xlabel("MC cycles")
plt.ylabel(r"$\sigma$")


(alphas, 
runtime, 
acceptance_ratio, 
sampler_param, 
steps, 
mean_E, 
std_E, 
std_E_blocking) = read_data_file(filename_dt)


plt.figure(1, figsize=(3, 2.25))
plt.plot(alphas, mean_E[-1, :] / N, ".-", label=rf"Importance $\Delta t=${dt}")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\langle E \rangle$")
plt.legend()
plt.savefig("E_vs_alphas_comp.eps", bbox_inches="tight")

plt.figure(2, figsize=(3, 2.25))
plt.plot(alphas, std_E[-1, :], ".-", label=rf"Importance $\Delta t=${dt}")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\sigma$")
plt.legend()
plt.savefig("std_E_vs_alphas_comp.eps", bbox_inches="tight")

plt.figure(3, figsize=(3, 2.25))
idx = np.where(alphas == a)[0][0]
plt.plot(np.log2(steps[:, idx]), np.abs(mean_E[:, idx] - E_exact), label=rf"Importance $\Delta t=${dt}")
plt.xticks(np.log2(steps[::2, idx]), [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle E \rangle$")
plt.legend()
plt.savefig("E_vs_mc_cycles_comp.eps", bbox_inches="tight")

plt.figure(4, figsize=(3, 2.25))
idx = np.where(alphas == a)[0][0]
plt.plot(np.log2(steps[:, idx]), std_E[:, idx], label=rf"Importance $\Delta t=${dt}")
plt.xticks(np.log2(steps[::2, idx]), [r"$2^{13}$", r"$2^{15}$", r"$2^{17}$", r"$2^{19}$", r"$2^{21}$", r"$2^{23}$"])
plt.xlabel("MC cycles")
plt.ylabel(r"$\sigma$")
plt.legend()
plt.savefig("std_E_vs_mc_cycles_comp.eps", bbox_inches="tight")

for i, filename in enumerate(filenames_sl):
    sl = sl_vals[i]
    
    (alphas, 
    runtime, 
    acceptance_ratio, 
    sampler_param, 
    steps, 
    mean_E, 
    std_E, 
    std_E_blocking) = read_data_file(filename)

    idx = np.where(alphas == a)[0][0]
    print(f"{sl} & {mean_E[-1, idx]:.6f} & {std_E[-1, idx]:.6f} & {std_E_blocking[-1, idx]:.6f} & {acceptance_ratio[idx]:.3f} & {runtime[idx]:.3f} \\\\")
        
print()

for i, filename in enumerate(filenames_dt):
    dt = dt_vals[i]
    
    (alphas, 
    runtime, 
    acceptance_ratio, 
    sampler_param, 
    steps, 
    mean_E, 
    std_E, 
    std_E_blocking) = read_data_file(filename)

    idx = np.where(alphas == a)[0][0]
    print(f"{dt} & {mean_E[-1, idx]:.6f} & {std_E[-1, idx]:.6f} & {std_E_blocking[-1, idx]:.6f} & {acceptance_ratio[idx]:.3f} & {runtime[idx]:.3f} \\\\")
    
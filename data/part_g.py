from cProfile import run
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib

from process_data import block

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

def f(filename):
    with open(directory_name + filename, "r") as file:
        line = file.readline()
        fread = file.readline().strip().split(",")
        mc_cycles = int(fread[0])
        measure_after = int(fread[1])

    energies = np.zeros((mc_cycles, 1))
    alphas = np.zeros(1)
    acceptance_ratio = np.zeros(1)
    runtime = np.zeros(1)
    sampling_param = np.zeros(1)

    with open(directory_name + filename, "r") as file:
        file.readline()
        fread = file.readline().strip().split(",")
        file.readline()
        data = np.loadtxt(file, delimiter=",")
        
        alphas = float(fread[2])
        runtime = float(fread[3])
        acceptance_ratio = float(fread[4])
        sampling_param = float(fread[5])
        energies[:, 0] = data[:, 0]
        
    steps = np.power(2, np.arange(int(np.log2(mc_cycles - measure_after)) + 1, int(np.log2(mc_cycles)) + 1))
    mean_E = np.zeros((len(steps), 1))
    std_E = np.zeros((len(steps), 1))
    std_E_blocking = np.zeros((len(steps), 1))

    for i, step in enumerate(steps):
        mean_E[i] = np.mean(energies[2**(int(np.log2(mc_cycles - measure_after))):step])
        std_E[i] = np.std(energies[2**(int(np.log2(mc_cycles - measure_after))):step])
        tmp, _ = block(energies[2**(int(np.log2(mc_cycles - measure_after)) ):step])
        std_E_blocking[i] = np.sqrt(tmp)

    return alphas, runtime, acceptance_ratio, sampling_param, steps, mean_E, std_E, std_E_blocking

N_vals = [10, 50, 100]

directory_name = "./data_files/part_g_interaction/"
filenames = ["data_interaction_N10_alpha0.505400.csv", "data_interaction_N50_alpha0.505400.csv", "data_interaction_N100_alpha0.505400.csv"]

alphas = list()
runtime = list()
acceptance_ratio = list()
sampling_param = list()
steps = list()
mean_E = list()
std_E = list() 
std_E_blocking = list()

for i, filename in enumerate(filenames):
    alphas_t, runtime_t, acceptance_ratio_t, sampling_param_t, steps_t, mean_E_t, std_E_t, std_E_blocking_t = f(filename)
    alphas.append(alphas_t)
    runtime.append(runtime_t)
    acceptance_ratio.append(acceptance_ratio_t)
    sampling_param.append(sampling_param_t)
    steps.append(steps_t)
    mean_E.append(mean_E_t)
    std_E.append(std_E_t)
    std_E_blocking.append(std_E_blocking_t)

for i in range(len(N_vals)):
    print(f"{N_vals[i]} & {mean_E[i][-1, 0]:.6f} & {mean_E[i][-1, 0] / N_vals[i]:.6f} & {std_E[i][-1, 0]:.6f} & {std_E_blocking[i][-1, 0]:.6f} & {runtime[i]:.3f} \\\\")
    
N_vals = [10, 50, 100]

directory_name = "./data_files/part_g_interaction/"
filenames = ["data_non_interaction_N10_alpha0.505400.csv", "data_non_interaction_N50_alpha0.505400.csv", "data_non_interaction_N100_alpha0.505400.csv"]

alphas = list()
runtime = list()
acceptance_ratio = list()
sampling_param = list()
steps = list()
mean_E = list()
std_E = list() 
std_E_blocking = list()

for i, filename in enumerate(filenames):
    alphas_t, runtime_t, acceptance_ratio_t, sampling_param_t, steps_t, mean_E_t, std_E_t, std_E_blocking_t = f(filename)
    alphas.append(alphas_t)
    runtime.append(runtime_t)
    acceptance_ratio.append(acceptance_ratio_t)
    sampling_param.append(sampling_param_t)
    steps.append(steps_t)
    mean_E.append(mean_E_t)
    std_E.append(std_E_t)
    std_E_blocking.append(std_E_blocking_t)

for i in range(len(N_vals)):
    print(f"{N_vals[i]} & {mean_E[i][-1, 0]:.6f} & {mean_E[i][-1, 0] / N_vals[i]:.6f} & {std_E[i][-1, 0]:.6f} & {std_E_blocking[i][-1, 0]:.6f} & {runtime[i]:.3f} \\\\")
    


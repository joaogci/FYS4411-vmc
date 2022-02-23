from cProfile import label
from statistics import variance
import numpy as np
import matplotlib.pyplot as plt

N = 100
n_alpha = 20
omega_ho = 1
cycles_exp = 6

E_exact = 0.5 * N * omega_ho

alpha = np.zeros(n_alpha)
energy = np.zeros(n_alpha)
variance = np.zeros(n_alpha)

file_name = "results/vmc_ho_N_" + str(N) +  "_omega_" + str(omega_ho) + "_MC_1E" + str(cycles_exp) + "_na_" + str(n_alpha) + ".txt"
with open(file_name, "r") as file:
    i = 0
    for line in file:
        line = line.strip().split(",")
        alpha[i] = float(line[0])
        energy[i] = float(line[1])
        variance[i] = float(line[2])
        i += 1

var_min = np.where(variance == np.min(variance))[0][0]
e_min = np.where(energy == np.min(energy))[0][0]
        
plt.figure(1)

plt.subplot(2, 1, 1)
plt.title("Energies vs Alpha")
plt.plot(alpha, energy, "-.b")
plt.plot(alpha[var_min], energy[var_min], "or", label="Min Variance")
plt.plot(alpha[e_min], energy[e_min], "ok", label="Min Energy")
plt.plot(alpha, E_exact * np.ones(n_alpha), '--g', label="Exact Energy")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$E$")
plt.legend()

plt.subplot(2, 1, 2)
plt.title("Variance vs Alpha")
plt.plot(alpha, variance, "-.b")
plt.plot(alpha[var_min], variance[var_min], "or", label="Min Variance")
plt.plot(alpha[e_min], variance[e_min], "ok", label="Min Energy")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\sigma$")
plt.legend()

plt.show()

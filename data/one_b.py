import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

N = 10

rho_int = np.loadtxt("one_body_density.csv")
r_int = np.loadtxt("r_one_body_density.csv")

rho_non_int = np.loadtxt("one_body_density_non_int.csv")
r_non_int = np.loadtxt("r_one_body_density_non_int.csv")

rho_int = rho_int * N / (np.sum(rho_int))
rho_non_int = rho_non_int * N / (np.sum(rho_non_int))

rho_int, r_int = zip(*sorted(zip(rho_int, r_int)))
rho_non_int, r_non_int = zip(*sorted(zip(rho_non_int, r_non_int)))

plt.figure(1, figsize=(3, 2.25))
plt.plot(r_int, rho_int, '-', label="Interacting")
plt.plot(r_non_int, rho_non_int, '-', label="Non Interacting")
plt.xlabel(r"$r$")
plt.ylabel(r"$\rho (r)$")
plt.xlim([0, np.max(r_int)])
plt.legend()

plt.savefig("./figs/one_body_densities.eps", bbox_inches="tight")

plt.show()

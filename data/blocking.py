import sys
import numpy as np

def autocov(arr, N):
    s1 = 0
    s2 = 0
    meanarr = np.mean(arr)

    for k in range(1, N):
        s1 += arr[k] - meanarr

    for k in range(N-1):
        s2 += arr[k] - meanarr

    return s1*s2/N

def block(x):
    n = len(x)
    d = int(np.log2(n))
    sample_var = np.zeros(d)
    gamma = np.zeros(d)

    for i in range(d):
        n = len(x)
        gamma[i] = autocov(x, n)
        sample_var[i] = np.var(x)

        x_1 = x[0::2]
        x_2 = x[1::2]

        if (len(x_1) > len(x_2)):
            x_1 = x_1[:-1]
        elif (len(x_2) > len(x_1)):
            x_2 = x_2[:-1]

        x = 1/2 * (x_1 + x_2)

    # Generate M_k values by using magic
    M = (np.cumsum( ( (gamma/sample_var)**2 * 2**np.arange(1, d + 1)[::-1] )[::-1] ) )[::-1]

    # Chi squared numbers
    q = np.array([6.634897, 9.210340, 11.344867, 13.276704, 15.086272, 16.811894, \
            18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, \
            27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, \
            36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, \
            44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    for k in range(d):
        if (M[k] < q[k]):
            res = sample_var[k]/(2**(d - k))
            break
        if (k >= d - 1):
            print("More data needed")

    return res, k


if __name__ == "__main__":

    path = sys.argv[1]      #"../data/N100_d3/data_alpha0.300000.csv"
    infile = open(path)

    infile.readline()
    mc_cycles, measure_cycles, alpha, runtime, a_ratio = infile.readline().split(",")
    mc_cycles = int(mc_cycles)
    measure_cycles = int(measure_cycles)
    alpha = float(alpha)
    runtime = float(runtime)
    a_ratio = float(a_ratio)

    infile.readline()

    data = np.loadtxt(path, delimiter=',', skiprows=3)
    print("Read successful")
    start = int(mc_cycles - measure_cycles)
    print(f"mc cycles: {mc_cycles}, log2 mc cycles: {np.log2(mc_cycles)}")
    print(f"start idx: {start}, log2 start idx: {np.log2(start)}")

    energy, energy2 = data[:, 0][start:], data[:, 1][start:]

    sigma_sqrd, blocking_iterations = block(energy)
    sigma = np.sqrt(sigma_sqrd)

    print("variance after blocking:")
    print(sigma_sqrd)

    print("std after blocking:")
    print(sigma)

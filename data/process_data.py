import os
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
    np.seterr(divide='ignore', invalid='ignore')
    M = (np.cumsum( ( (gamma/sample_var)**2 * 2**np.arange(1, d + 1)[::-1] )[::-1] ) )[::-1]

    # Chi squared numbers
    q = np.array([6.634897, 9.210340, 11.344867, 13.276704, 15.086272, 16.811894, \
            18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, \
            27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, \
            36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, \
            44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    for k in range(d):
        if (M[k] < q[k]):
            break
    
    res = sample_var[k]/(2**(d - k))
    return res, k

def process_data(directory_name, save_file_name, before=0):
    files = os.listdir(directory_name)
    try:
        files.remove(".DS_Store")
    except:
        ...
    files.sort()

    with open(directory_name + files[0], "r") as file:
        line = file.readline()
        fread = file.readline().strip().split(",")
        mc_cycles = int(fread[0])
        measure_after = int(fread[1])

    energies = np.zeros((mc_cycles, len(files)))
    alphas = np.zeros(len(files))
    acceptance_ratio = np.zeros(len(files))
    runtime = np.zeros(len(files))
    sampling_param = np.zeros(len(files))

    for j, filename in enumerate(files):
        with open(directory_name + filename, "r") as file:
            file.readline()
            fread = file.readline().strip().split(",")
            file.readline()
            data = np.loadtxt(file, delimiter=",")
            
            alphas[j] = float(fread[2])
            runtime[j] = float(fread[3])
            acceptance_ratio[j] = float(fread[4])
            sampling_param[j] = float(fread[5])
            energies[:, j] = data[:, 0]
            
        print("read alpha =", alphas[j], end="\r")
    print()
    print("Read files with success!")

    steps = np.power(2, np.arange(int(np.log2(mc_cycles - measure_after)) + 1 - before, int(np.log2(mc_cycles)) + 1))
    mean_E = np.zeros((len(steps), len(alphas)))
    std_E = np.zeros((len(steps), len(alphas)))
    std_E_blocking = np.zeros((len(steps), len(alphas)))

    for j in range(len(alphas)):
        for i, step in enumerate(steps):
            mean_E[i, j] = np.mean(energies[mc_cycles - measure_after - np.power(2, before):step, j])
            std_E[i, j] = np.std(energies[mc_cycles - measure_after - np.power(2, before):step, j])
            tmp, _ = block(energies[mc_cycles - measure_after - np.power(2, before):step, j])
            std_E_blocking[i, j] = np.sqrt(tmp)
        
        print("processed alpha =", alphas[j], end="\r")
    print()

    print("Writing output file.")
    with open(save_file_name, "w") as file:
        file.write("n_steps,n_alphas\n")
        file.write(f"{len(steps)},{len(alphas)}\n")
        
        for j in range(len(alphas)):
            file.write("\nalpha,runtime,acceptance_ratio,sampling_param\n")
            file.write(",".join((str(alphas[j]), str(runtime[j]), str(acceptance_ratio[j]), str(sampling_param[j]))) + "\n")
            
            file.write("step,mean_E,std_E,std_E_blocking\n")
            for i, step in enumerate(steps):
                file.write(",".join((str(step), str(mean_E[i, j]), str(std_E[i, j]), str(std_E_blocking[i, j]))) + "\n")

def read_data_file(filename):
    with open(filename, "r") as file:
        file.readline()
        n_steps, n_alphas = [int(x) for x in file.readline().strip().split(",")]

    alphas = np.zeros(n_alphas)
    runtime = np.zeros(n_alphas)
    acceptance_ratio = np.zeros(n_alphas)
    sampler_param = np.zeros(n_alphas)

    steps = np.zeros((n_steps, n_alphas))
    mean_E = np.zeros((n_steps, n_alphas))
    std_E = np.zeros((n_steps, n_alphas))
    std_E_blocking = np.zeros((n_steps, n_alphas))

    with open(filename, "r") as file:
        file.readline()
        file.readline()
        file.readline()
        
        for j in range(n_alphas):
            file.readline()
            alphas[j], runtime[j], acceptance_ratio[j], sampler_param[j] = [float(x) for x in file.readline().strip().split(",")]
            file.readline()
            
            for i in range(n_steps):
                steps[i, j], mean_E[i, j], std_E[i, j], std_E_blocking[i, j] = [float(x) for x in file.readline().strip().split(",")]
            
            file.readline()

    return alphas, runtime, acceptance_ratio, sampler_param, steps, mean_E, std_E, std_E_blocking


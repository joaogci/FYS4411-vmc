import numpy as np
import matplotlib.pyplot as plt

def wave_function(r, alpha, N, dim, a=0.00433):
    r2_sum = 0
    for i in range(N):
        for d in range(dim):
            if d != 3:
                r2_sum += r[i, d] * r[i, d] 
            elif d == 3:
                r2_sum += np.sqrt(8) * r[i, d] * r[i, d]
    single_part = np.exp(- alpha * r2_sum)

    jastrow = 1.0
    for i in range(N):
        for j in range(i):
            dist = np.linalg.norm(r[j, :] - r[i, :])

            if dist > a:
                jastrow *= 1.0 - a / dist
            else:
                jastrow *= 0

    return single_part * jastrow

def one_body_density(M, m, lim, alpha, N, dim, a=0.00433):
    r1 = np.linspace(-lim, lim, m)

    one_body = np.zeros(m**3)
    r = np.zeros((N, dim))
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                for i in range(M):
                    r[0, :] = np.array([r1[ix], r1[iy], r1[iz]])
                    r[1:, :] = np.random.rand(N - 1, dim)
                    wf = wave_function(r, alpha, N, dim, a)
                    one_body[ix * m**2 + iy * m + iz] += wf * wf
                one_body[ix * m**2 + iy * m + iz] /= M
                print(f"{ix=}; {iy=}; {iz=}", end='\r')
    
    r = np.zeros(m**3)
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                r[ix * m**2 + iy * m + iz] = np.sqrt(r1[ix]**2 + r1[iy]**2 + r1[iz]**2)

    return r, one_body


if __name__ == '__main__':
    N = 10
    d = 3
    alpha = 0.5054
    
    r, one_body = one_body_density(100000, 10, 3, 0.5, 10, 3)

    np.savetxt("one_body_density.csv", one_body)
    np.savetxt("r_one_body_density.csv", one_body)

    plt.plot(r, one_body, '.')
    plt.savefig("one_body_test.eps")


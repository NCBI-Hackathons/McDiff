import numpy as np
import sims
import matplotlib.pyplot as plt

def proposal(x, sigma, xmin, xmax):
    x_new = x + np.random.normal(0, sigma)
    breakage = False
    cnt = 0
    while (x_new > xmax) | (x_new < xmin):
        x_new = x + np.random.normal(0, sigma)
        cnt += 1
        if cnt > 100:
            breakage = True
            break
    if breakage == True:
        return x, breakage
    else:
        return x_new, breakage

def MCMC(D0, f_mobile0, f_bleached, nuc, roi, N, T, sigma1, sigma2, fmin, fmax, dmin, dmax, sim_len, data, data_pre, data_norm):
    s2 = 2*T**2. #temperature
    all_params = np.zeros((2, N+1)) #store parameters
    all_params[0,0] = D0
    all_params[1,0] = f_mobile0
    errores = np.zeros(N+1) #store error
    stuck, roi_pre = sims.simulate(D0, f_mobile0, 0.5, nuc, roi, sim_len) #do a simulation
    stuck_norm = stuck / roi_pre
    stuck_time = np.arange(sim_len+1) * 0.1
    error = sims.compute_error(data, data_norm, stuck_time, stuck_norm)
    chi2 = np.sum(error)
    errores[0] = chi2
    old_params = [D0, f_mobile0]
    new_params = [0, 0]
    for i in range(1, N+1):
        new_params[0], b1 = proposal(old_params[0], sigma1, dmin, dmax)
        new_params[1], b2 = proposal(old_params[1], sigma2, fmin, fmax)
        if (b1 == True) | (b2 == True):
            return old_params, errores, all_params, b1, b2, i
        else:
            stuck, roi_pre = sims.simulate(new_params[0], new_params[1], 0.5, nuc, roi, sim_len)
            stuck_norm = stuck / roi_pre
            stuck_time = np.arange(sim_len+1) * 0.1
            error = sims.compute_error(data, data_norm, stuck_time, stuck_norm)
            chi2_new = np.sum(error)
            if chi2_new  < chi2:
                old_params = new_params
                chi2 = chi2_new
                all_params[0,i] = new_params[0]
                all_params[1,i] = new_params[1]
                errores[i] = chi2
                print("Updated Parameters")
            else:
                coin = np.random.uniform()
                r = np.exp(chi2 - chi2_new)/s2 #likelihood function
                if coin < r:
                    old_params = new_params
                    chi2 = chi2_new
                    all_params[0,i] = new_params[0]
                    all_params[1,i] = new_params[1]
                    errores[i] = chi2
                    print("Updated Parameters")
                else:
                    all_params[0,i] = old_params[0]
                    all_params[1,i] = old_params[1]
                    errores[i] = chi2
                    print("Keeping Old Parameters")
    return old_params, errores, all_params, b1, b2, i

def rand_sam(f_bleached, nuc, roi, N, fmin, fmax, dmin, dmax):
    D = np.random.uniform(dmin, dmax, N)
    F = np.random.uniform(fmin, fmax, N)
    E = np.zeros(N)
    for i in range(N):
        stuck, roi_pre = simulate(D[i], F[i], f_bleached, nuc, roi, sim_len)
        stuck_norm = stuck / roi_pre
        stuck_time = np.arange(sim_len+1) * 0.1
        error = compute_error(data, data_norm, stuck_time, stuck_norm)
        E[i] = np.sum(error)
    return D, F, E
def course_fine(f_bleached, nuc, roi, N, fmin, fmax, dmin, dmax, r1, r2, N1, N2):
    D,F,E = rand_sam(f_bleached, nuc, roi, N1, fmin, fmax, dmin, dmax)
    x = E.argmin()
    D2, F2, E2 = rand_sam(f_bleached, nuc, roi, N2, F[x] - r1, F[x] + r1, D[x] - r2, D[x] + r2)
    return D,F,E, D2, F2, E2

def CF(f_bleached, nuc, roi, fmin, fmax, dmin, dmax, s1, s2, N, L):
    Params = np.zeros((3, N*(L+1)))
    D,F,E = rand_sam(f_bleached, nuc, roi, N, fmin, fmax, dmin, dmax)
    Params[0,0:N] = D
    Params[1,0:N] = F
    Params[2,0:N] = E
    x = E.argmin()
    for i in range(L):
        D2, F2, E2 = rand_sam(f_bleached, nuc, roi, N, F[x] - (F[x] - fmin)*s1, F[x] + (F[x] + fmax)*s1, D[x] - (D[x] - dmin)*s2, D[x] + (D[x] + dmax)*s2)
        Params[0, (N*(i+1)):(N*(i+1)) + N] = D2
        Params[1, (N*(i+1)):(N*(i+1)) + N] = F2
        Params[2, (N*(i+1)):(N*(i+1)) + N] = E2
    return Params

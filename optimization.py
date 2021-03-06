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

def MCMC(D0, f_mobile0, f_bleached, nuc, roi, N, T, sigma1, sigma2, fmin, fmax, dmin, dmax, sim_len, data, data_pre, data_norm, x0, y0):
    s2 = 2*T**2. #temperature
    all_params = np.zeros((2, N+1)) #store parameters
    all_params[0,0] = D0
    all_params[1,0] = f_mobile0
    errores = np.zeros(N+1) #store error
    stuck_norm = sims.simulate(D0, f_mobile0, 0.5, nuc, roi, sim_len, x0, y0) #do a simulation
    stuck_time = np.arange(sim_len+1) * 0.18
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
            stuck_norm = sims.simulate(new_params[0], new_params[1], 0.5, nuc, roi, sim_len, x0, y0)
            stuck_time = np.arange(sim_len+1) * 0.18
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
                #chi2_new > chi2, thus exponent is negative, numerator is between 0 and 1
                #the worse the proposed parameters are, the less likely you are to accept them
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

def rand_sam(f_bleached, nuc, roi, N, fmin, fmax, dmin, dmax, x0, y0, sim_len, data, data_norm):
    D = np.random.uniform(dmin, dmax, N)
    F = np.random.uniform(fmin, fmax, N)
    E = np.zeros(N)
    for i in range(N):
        stuck_norm = sims.simulate(D[i], F[i], f_bleached, nuc, roi, sim_len, x0, y0)
        stuck_time = np.arange(sim_len+1) * 0.18
        error = sims.compute_error(data, data_norm, stuck_time, stuck_norm)
        E[i] = np.sum(error)
    return D, F, E

def CF(f_bleached, nuc, roi, fmin, fmax, dmin, dmax, s1, s2, N, L, x0, y0, sim_len, data, data_norm):
    #random search, pick the best point, scale the parameter bounds inwards by some factor, repeat over L layers
    Params = np.zeros((3, N*(L+1)))
    D,F,E = rand_sam(f_bleached, nuc, roi, N, fmin, fmax, dmin, dmax, x0, y0, sim_len, data, data_norm)
    Params[0,0:N] = D
    Params[1,0:N] = F
    Params[2,0:N] = E
    x = E.argmin()
    bounds = np.zeros(4)
    bounds[0] = F[x] - (F[x] - fmin)*s1
    bounds[1] = F[x] + (fmax - F[x])*s1
    bounds[2] = D[x] - (D[x] - dmin)*s2
    bounds[3] = D[x] + (dmax - D[x])*s2
    print("initial", bounds)
    for i in range(L):
        D2, F2, E2 = rand_sam(f_bleached, nuc, roi, N, bounds[0], bounds[1], bounds[2], bounds[3], x0, y0, sim_len, data, data_norm)
        Params[0, (N*(i+1)):(N*(i+1)) + N] = D2
        Params[1, (N*(i+1)):(N*(i+1)) + N] = F2
        Params[2, (N*(i+1)):(N*(i+1)) + N] = E2
        x = E2.argmin()
        #print(F2[x], " is my best f coordinate")
        #print(bounds[0], "is my previous lower bound")
        #print(F2[x] - bounds[0], "is the difference from previous bound")
        #print(s1*(F2[x] - bounds[0]), "is the difference from previous bound times s")
        #print(F2[x] - (F2[x] - bounds[0])*s1, "is proposed new bound")
        bounds[0] = F2[x] - (F2[x] - bounds[0])*s1
        bounds[1] = F2[x] + (bounds[1]- F2[x])*s1
        bounds[2] = D2[x] - (D2[x] - bounds[2])*s2
        bounds[3] = D2[x] + (bounds[3] - D2[x])*s2
        print(bounds)
    return Params

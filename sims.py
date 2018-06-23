import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon
from scipy.interpolate import interp1d
from descartes import PolygonPatch
from shapely import geometry
import shapely.vectorized

# exec(open("initialize_points_w.py", 'r').read())

def parse_mask(x):
    with open(x) as f:
        a = f.read()
    b = a.split('\n')
    c = [i.split(', ') for i in b]
    c.pop()
    d = np.zeros((2, len(c)))
    for i in range(len(c)):
        d[0,i] = float(c[i][0])
        d[1,i] = float(c[i][1])
    return d


def parse_roi(x):
    with open(x) as f:
        a = f.read()
    b = a.split(', ')
    d = np.zeros(len(b))
    for i in range(len(b)):
        d[i] = float(b[i])
    return d


def parse_data(data, offset):
    with open(data, 'r') as f:
        a = f.read()
    b = a.split('\n')
    d = np.zeros((2, len(b)-5))
    for i in range(len(b)-5):
        c = b[i+4].split(',')
        d[0,i] = float(c[0])
        d[1,i] = float(c[11])
    x1 = d[:,:6]
    x2 = d[:,6:]
    x1[0,:] -= offset
    x2[0,:] -= offset
    return x1, x2


def generate_random_points(N, poly):
	list_of_points = np.zeros((2, N))
	minx,miny,maxx,maxy = poly.bounds
	counter = 0
	x = np.random.uniform(minx, maxx, N*3)
	y = np.random.uniform(miny, maxy, N*3)
	innie = shapely.vectorized.contains(poly, x, y)
	return x[innie][:N], y[innie][:N]

def D_2_x(D, h):
    x = np.sqrt((2*D*h))
    return x

def update_positions(x_cord, y_cord, mu, sigma, nucleus, roi):
    l = len(x_cord)
    x = np.random.normal(mu, sigma, l)
    y = np.random.normal(mu, sigma, l)
    fx = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    fy = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    x *= fx
    y *= fy
    x_new = x_cord + x
    y_new = y_cord + y
    #if you kicked a protein outside the nucleus, restore to initial position
    # out_nuc = check_inside(points_new, nucleus) == False
    out_nuc = shapely.vectorized.contains(nucleus, x_new, y_new) == False
    x_new[out_nuc] = x_cord[out_nuc]
    y_new[out_nuc] = y_cord[out_nuc]
    #now check if you have particles stuck in the middle
    in_roi = shapely.vectorized.contains(roi, x_new, y_new)
    out_roi = in_roi == False
    N_stuck = np.sum(in_roi)
    # return x_new[out_roi], y_new[out_roi], x_new[in_roi], y_new[in_roi]
    return x_new[out_roi], y_new[out_roi], N_stuck


def simulate(D, f_mobile, f_bleached, nuc, roi, runtime):
    h = 0.18 #ms 0.19 or 0.16 in other code
    microns_2_pixels = .08677
    N = 12000
    dx = D_2_x(D, h) / microns_2_pixels # dx is in pixels
    x, y = generate_random_points(12000, nuc) #positions of all points
    in_roi = shapely.vectorized.contains(roi, x, y) #return a boolean mask
    out_roi = in_roi == False
    N0_roi = np.sum(in_roi)
    stuck = int(N0_roi * (1 - f_bleached))
    all_stuck = np.zeros(runtime+1)
    all_stuck[0] = stuck
    N_sim = int((N - N0_roi) * f_mobile) #N = (N*f_mobile) - (in_roi * f_mobile)
    x = x[out_roi][:N_sim]
    y = y[out_roi][:N_sim]
    for i in range(1,runtime+1):
        x, y, N_stuck = update_positions(x, y, dx, 0.001, nuc, roi)
        stuck += N_stuck
        all_stuck[i] = stuck
    return all_stuck, N0_roi#, x_stuck, y_stuck

def compute_error(data, data_norm, stuck_time, stuck_norm):
    times = np.zeros(len(data[0,:]))
    y_ = np.zeros(len(times))
    error = np.zeros(len(times))
    for i in range(len(times)):
        dx = (np.abs(stuck_time - data[0,i])).argmin()
        times[i] = stuck_time[dx]
        y_[i] = stuck_norm[dx]
        error[i] = (stuck_norm[dx] - data_norm[i])**2.
    return error #np.sqrt(error)


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


def MCMC(D0, f_mobile0, f_bleached, nuc, roi, N, T, sigma1, sigma2, fmin, fmax, dmin, dmax):
    s2 = 2*T**2. #temperature
    all_params = np.zeros((2, N+1)) #store parameters
    all_params[0,0] = D0
    all_params[1,0] = f_mobile0
    errores = np.zeros(N+1) #store error
    stuck, roi_pre = simulate(D0, f_mobile0, 0.5, nuc, roi, sim_len) #do a simulation
    stuck_norm = stuck / roi_pre
    stuck_time = np.arange(sim_len+1) * 0.1
    error = compute_error(data, data_norm, stuck_time, stuck_norm)
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
            stuck, roi_pre = simulate(new_params[0], new_params[1], 0.5, nuc, roi, sim_len)
            stuck_norm = stuck / roi_pre
            stuck_time = np.arange(sim_len+1) * 0.1
            error = compute_error(data, data_norm, stuck_time, stuck_norm)
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

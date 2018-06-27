import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon
from scipy.interpolate import interp1d
from descartes import PolygonPatch
from shapely import geometry
import shapely.vectorized

def parse_mask(x):
    ''' Parses a file containing coordinate points for a nucleus outline.
        Input file is a txt file with coordinate points defining the boundaries
        as ordered pairs (i.e., each row has x coord, y coord)
        Output is a 2-d numpy array, where row 1 has x, row 2 had y'''
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
    '''Reads in a .txt file with points defining the region of interest
       within the nucleus. Returns a one-D numpy array.'''
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


def init_sim(N, nuc):
    x, y = generate_random_points(N, nuc) #positions of all points
    return x, y

def simulate(D, f_mobile, f_bleached, nuc, roi, runtime, x0, y0):
    h = 0.18 #ms 0.19 or 0.16 in other code
    microns_2_pixels = .08677
    N = 12000
    dx = D_2_x(D, h) / microns_2_pixels # dx is in pixels
    in_roi = shapely.vectorized.contains(roi, x0, y0) #return a boolean mask
    out_roi = in_roi == False
    N0_roi = np.sum(in_roi)
    stuck = int(N0_roi * (1 - f_bleached))
    all_stuck = np.zeros(runtime+1)
    all_stuck[0] = stuck
    N_sim = int((N - N0_roi) * f_mobile) #N = (N*f_mobile) - (in_roi * f_mobile)
    x = x0[out_roi][:N_sim]
    y = y0[out_roi][:N_sim]
    # expected_N0 = (nuc.intersection(roi).area / nuc.area)
    for i in range(1,runtime+1):
        x, y, N_stuck = update_positions(x, y, dx, 0.001, nuc, roi)
        stuck += N_stuck
        all_stuck[i] = stuck
    all_stuck = all_stuck / (N * (nuc.intersection(roi).area / nuc.area))
    return all_stuck #, x_stuck, y_stuck


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

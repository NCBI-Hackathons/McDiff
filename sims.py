import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon
from shapely.vectorized import contains
from descartes import PolygonPatch

exec(open("initialize_points_w.py", 'r').read())

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

def D_2_x(D, h):
    x = np.sqrt((2*D*h))
    return x

def check_inside(points, poly):
    a = np.zeros(len(points[0,:]), dtype = bool)
    ps = np.array([Point(points[0,i], points[1,i]) for i in range(points.size//2)], dtype=object)
    for i in range(len(a)):
        a[i] = poly.contains(ps[i])
    return a
"""
def check_recursion(points, poly):
    if(len(points) != 0):

        check_recursion(points, poly)
    else:
        return 0

    return points
    """

def check_fast(points, poly):
    return shapely.vectorized.contains(poly, points[0,:], points[1,:])


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
    out_nuc = shapely.vectorized.contains(nuc, x_new, y_new) == False
    x_new = x_new[out_nuc]
    y_new = y_new[out_nuc]
    #now check if you have particles stuck in the middle
    in_roi = shapely.vectorized.contains(roi, x_new, y_new)
    out_roi = in_roi == False
    N_stuck = np.sum(in_roi)
    return x_new, y_new, N_stuck


def simulate(D, f_mobile, f_bleached, nuc, roi, runtime):
    h = 0.1 #ms 0.19 or 0.16 in other code
    microns_2_pixels = .08677
    N = 12000
    x = D_2_x(D, h)
    points = generate_random_points(12000, nuc) #positions of all points
    in_roi = check_inside(points, roi)#return a boolean mask
    out_roi = in_roi == False
    stuck = np.sum(in_roi)
    stuck *= f_bleached
    N_sim = int((N - stuck) * f_mobile) #N = (N*f_mobile) - (in_roi * f_mobile)
    points = points[:,out_roi][:,0:N_sim]
    x = points[0,:]
    y = points[1,:]
    for i in range(runtime):
        x, y, points_stuck = update_positions(x, y, 10, 0.01, nuc, roi)
        stuck += points_stuck
        print(i)
    return points, stuck

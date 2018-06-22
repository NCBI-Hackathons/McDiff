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

def update_positions(points, mu, sigma, nucleus, roi):
    l = len(points[0,:])
    x = np.random.normal(mu, sigma, l)
    y = np.random.normal(mu, sigma, l)
    fx = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    fy = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    x *= fx
    y *= fy
    points_new = np.zeros((2, l))
    points_new[0,:] += points[0,:] + x
    points_new[1,:] += points[1,:] + y
    #if you kicked a protein outside the nucleus, restore to initial position
    out_nuc = check_inside(points_new, nucleus) == False
    points_new[0,out_nuc] = points[0,out_nuc]
    points_new[1,out_nuc] = points[1,out_nuc]
    #now check if you have particles stuck in the middle
    in_roi = check_inside(points_new, roi)
    out_roi = in_roi == False
    N_stuck = np.sum(in_roi)
    return points_new[:, out_roi], N_stuck


def simulate(D, f_mobile, f_bleached, nuc, roi):
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
    for i in range(50):
        points, points_stuck = update_positions(points, 10, 0.01, nuc, roi)
        stuck += points_stuck
        print(i)
    return points, stuck

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon

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

f = "1.31.18_GFPP1_Hela_1min_002NuclMask.txt"
mask = parse_mask(f)

a = LineString(list(zip(mask[0,:], mask[1,:])))
x, y = a.xy
plt.plot(x, y, "g.")
plt.plot(mask[0,:], mask[1,:], "b.", alpha = .1)
a = Polygon(list(zip(mask[0,:], mask[1,:])))
#sanity check
a.contains(Point(250, 250))
a.contains(Point(250, 2000))
#


#check that we get the same shape

def parse_roi(x):
    with open(x) as f:
        a = f.read()
    b = a.split(', ')
    d = np.zeros(len(b))
    for i in range(len(b)):
        d[i] = float(b[i])
    return d

roi = "1.31.18_GFPP1_Hela_1min_002ROI.txt"
d = parse_roi(roi)

plt.plot(mask[0,:], mask[1,:], ".")
plt.plot(d[0], d[1], ".")
plt.plot(d[0] + d[2], d[1], ".")
plt.plot(d[0], d[1] + d[3], ".")
plt.plot(d[0] + d[2], d[1] + d[3], ".")

#first 6 frames nothing happens
#

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

# def parse_data(data):
#     with open(data, 'r') as f:
#         a = f.read()
#     b = a.split('\n')
#     d = np.zeros((2, len(b)-5))
#     for i in range(len(b)-5):
#         c = b[i+4].split(',')
#         d[0,i] = float(c[0])
#         d[1,i] = float(c[11])
#     return d

data_pre, data = parse_data("1.31.18_GFPP1_Hela_1min_002.csv")
plt.plot(data[0,:], data[1,:], ".")
plt.show()

#1. initialize 12,000
#
#
def D_2_x(D, h):
    x = np.sqrt((2*D*h))
    return x


x = np.random.normal(mu, sigma, l) * (((np.random.uniform(0, 1, l) > 0.5)*2) - 1)


def update_positions(points, mu, sigma, nucleus):
    l = len(points)
    x = np.random.normal(mu, sigma, l)
    y = np.random.normal(mu, sigma, l)
    fx = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    fy = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
    x *= fx
    y *= fy
    points_new = np.zeros((2, l))
    points_new[0,:] += points[0,:] + x
    points_new[1,:] += points[1,:] + y
    out_nuc = nucleus.contains(points_new) == False
    points_new[0,out_nuc] = points[0,out_nuc]
    points_new[1,out_nuc] = points[1,out_nuc]









def simualte(D, f_mobile, f_bleached, nucleus, roi):
    dt = 0.1 #ms 0.19 or 0.16 in other code
    microns_2_pixels = .08677
    N = 12000
    x = D_2_x(D, h)
    points = initialize_points #positions of all points
    in_roi = roi.contains(points) #return a boolean mask
    out_roi = in_roi == False
    N_sim = (N - np.sum(in_roi)) * f_mobile #N = (N*f_mobile) - (in_roi * f_mobile)
    points = points[out_roi][:N_sim]

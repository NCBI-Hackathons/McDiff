import numpy as np
import matplotlib.pyplot as plt

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

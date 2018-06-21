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
    return d[:2], d[2:]



roi = "1.31.18_GFPP1_Hela_1min_002ROI.txt"
roi = parse_roi(roi)



def parse_data(x):
    with open(x) as f:
        a = f.read()
    b = a.split('\n')



def sim_rand_walk(N, runtime, h, mask, roi, f, D, color):
    ntotal = int(runtime/h)

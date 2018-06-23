#actaully do some simulations
exec(open('sims.py').read())
exec(open('optimization.py').read())

f = "./test_files/1.31.18_GFPP1_Hela_1min_002NuclMask.txt"
mask = parse_mask(f)

# IV = generate_random_points(12000, a)

roi = "./test_files/1.31.18_GFPP1_Hela_1min_002ROI.txt"
d = parse_roi(roi)

# roi = Polygon([(d[0], d[1]), (d[0], d[1]+d[3]), (d[0]+d[2], d[1]+d[3]), (d[0]+d[2], d[1])])
roi = Polygon([(d[1], d[0]), (d[1], d[0]+d[2]), (d[1]+d[3], d[0]+d[2]), (d[1]+d[3], d[0])])
nuc = Polygon(list(zip(mask[0,:], mask[1,:])))

# fig, ax = plt.subplots(1)
# ax.plot(mask[0,:], mask[1,:], "b.", alpha = .1)
# ax.plot(IV[0,:], IV[1,:], "r.", ms = 1.)
# ring = PolygonPatch(roi_points)
# ax.add_patch(ring)

data_pre, data = parse_data("./test_files/1.31.18_GFPP1_Hela_1min_002.csv", 10.5)
#check_inside(IV, a)
data_norm = data[1,:] / np.mean(data_pre[1,:])

sim_len = 650
np.random.seed(1200)
x0, y0 = init_sim(12000, nuc)

#########################MCMC
#
# sigmaD = 2.
# sigmaF = .1
# N = 150
# OP, Error, AP, bool_flag_1, bool_flag_2, Iterate_ended = MCMC(5, .5, .5, nuc, roi, N, 1, sigmaD, sigmaF, 0, 1, 0, 20,sim_len, data, data_pre, data_norm, x0, y0)
#
# lo_mejor = Error.argmin() #index of parameter optimal
# los_mejores = AP[:, lo_mejor] #best parameters
#
# stuck_norm = simulate(los_mejores[0], los_mejores[1], 0.5, nuc, roi, sim_len, x0, y0)
#
# stuck_time = np.arange(sim_len+1) * 0.18 #converts array indices into seconds
#
# fig, ax = plt.subplots(4)
# ax[0].plot(stuck_time, stuck_norm, ".", label = "Simulation")
# ax[0].plot(data[0,:], data_norm, ".", label = "Data")
# ax[0].legend()
# ax[0].set_xlabel("Time (s)")
# ax[1].hist(AP[0,:], 20)
# ax[2].hist(AP[1,:], 20)
# ax[3].plot(Error)
# plt.show()

f_bleached = .54
N = 50
L = 3
s1 = .1
s2 = .1
fmin = 0
fmax = 1
dmin = 0
dmax = 20
results = CF(f_bleached, nuc, roi, fmin, fmax, dmin, dmax, s1, s2, N, L, x0, y0)

x = results[2,:].argmin()
results[2,x]
stuck_norm = simulate(results[0,x], results[1,x], 0.5, nuc, roi, sim_len, x0, y0)
stuck_time = np.arange(sim_len+1) * 0.18 #converts array indices into seconds
error = compute_error(data, data_norm, stuck_time, stuck_norm)


fig, ax = plt.subplots(2)
ax[0].scatter(results[0,:], results[1,:], 10., results[2,:], ".")
ax[1].plot(stuck_time, stuck_norm, ".", label = "Simulation")
ax[1].plot(data[0,:], data_norm, ".", label = "Data")
ax[1].legend()
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Fraction of Proteins Bound/Baseline")
plt.show()




    # newfile = open("MCMC_results.txt", 'w')
    # for i in range(len(E)):
    #     text = "{0} {1} {2}\n".format(E[i], AP[0,i], AP[1,i])
    #     newfile.write(text)
    #
    #

# fig, ax = plt.subplots(1)
# ax.plot(mask[0,:], mask[1,:], ".")
# ring = PolygonPatch(roi)
# ax.add_patch(ring)

# plt.plot(x, y, "r.", ms = 1.)
# plt.plot(x_stuck, y_stuck, "r.", ms = 1.)
# plt.show()



############################# Random Sampling

# N = 50
# f_bleached = .46
# fmin = 0
# fmax = 1
# dmin = 0
# dmax = 20
# D, F, E = rand_sam(f_bleached, nuc, roi, N, fmin, fmax, dmin, dmax)

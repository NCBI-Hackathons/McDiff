#actaully do some simulations
exec(open("sims.py").read())

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

#########################MCMC

sigmaD = 2.
sigmaF = .05
N = 200
OP, Error, AP, bool_flag_1, bool_flag_2, Iterate_ended = MCMC(20, .1, .5, nuc, roi, N, 1, sigmaD, sigmaF, 0, 1, 0, 20)

lo_mejor = Error.argmin() #index of parameter optimal
los_mejores = AP[:, lo_mejor] #best parameters

stuck, roi_pre = simulate(los_mejores[0], los_mejores[1], 0.5, nuc, roi, sim_len)

stuck_norm = stuck / roi_pre
stuck_time = np.arange(sim_len+1) * 0.18 #converts array indices into seconds

fig, ax = plt.subplots(1)
ax.plot(stuck_time, stuck_norm, ".", label = "Simulation")
ax.plot(data[0,:], data_norm, ".", label = "Data")
ax.legend()
ax.set_xlabel("Time (s)")
ax.set_ylabel("Fraction of Protiens Bound/Baseline")
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

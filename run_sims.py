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

x, y, stuck, roi_pre, x_stuck, y_stuck = simulate(6.69, 0.43 , 0.5, nuc, roi, sim_len)

stuck_norm = stuck / roi_pre
stuck_time = np.arange(sim_len+1) * 0.1

# fig, ax = plt.subplots(1)
# ax.plot(stuck_time, stuck_norm, ".", label = "Simulation")
# ax.plot(data[0,:], data_norm, ".", label = "Data")
# ax.legend()
# ax.set_xlabel("Time (s)")
# ax.set_ylabel("Fraction of Protiens Bound/Baseline")
# plt.show()

# fig, ax = plt.subplots(1)
# ax.plot(mask[0,:], mask[1,:], ".")
# ring = PolygonPatch(roi)
# ax.add_patch(ring)
# plt.plot(x, y, "r.", ms = 1.)
# plt.plot(x_stuck, y_stuck, "r.", ms = 1.)
# plt.show()

#testing pieces of the simulation
# f_mobile = .5
# dt = 0.1 #ms 0.19 or 0.16 in other code
# microns_2_pixels = .08677
# N = 12000
# # x = D_2_x(D, h)
# points = generate_random_points(12000, a) #positions of all points
# in_roi = check_inside(points, roi) #return a boolean mask
# out_roi = in_roi == False
# N_sim = int((N - np.sum(in_roi)) * f_mobile) #N = (N*f_mobile) - (in_roi * f_mobile)
# pointz = points[:,out_roi][:,0:N_sim]
# pointz_new, points_stuck = update_positions(pointz, 10, 0.01, nuc, roi)



#testing just the update position function
# mu = 10
# sigma = .01
# l = len(pointz[0,:])
# x = np.random.normal(mu, sigma, l)
# y = np.random.normal(mu, sigma, l)
# fx = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
# fy = ((np.random.uniform(0, 1, l) > 0.5)*2) - 1
# x *= fx
# y *= fy
# points_new = np.zeros((2, l))
# points_new[0,:] += pointz[0,:] + x
# points_new[1,:] += pointz[1,:] + y
#




#space to play with shapely.vectorized

#
# import shapely.vectorized
# minx,miny,maxx,maxy = nuc.bounds
# x = np.random.uniform(minx, maxx, 100)
# y = np.random.uniform(miny, maxy, 100)
# pnts = np.zeros((2, 100))
# pnts[0,:] = x
# pnts[1,:] = y
# # shapely.vectorized.contains(nuc, x, y)
# c1 = check_inside(pnts, nuc)
# c2 = check_fast(pnts, nuc)
# c1 == c2
#

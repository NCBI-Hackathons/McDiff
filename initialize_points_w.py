from shapely import geometry
# import random
# import numpy as np

# numParticles = 120

# point_list = [[0,0],[0,1],[1,1],[1,0]]
# poly = geometry.Polygon(point_list)

def generate_random_points(N, poly):
	list_of_points = np.zeros((2, N))
	minx,miny,maxx,maxy = poly.bounds
	counter = 0
	while counter < N:
		punto = (np.random.uniform(minx, maxx), np.random.uniform(miny,maxy))
		p = geometry.Point(punto)
		if poly.contains(p):
			list_of_points[0,counter] = punto[0]
			list_of_points[1,counter] = punto[1]
			counter += 1
	return list_of_points

# print generate_random_points(numParticles, poly)

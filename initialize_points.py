from shapely import geometry
import random
import numpy as np

numParticles = 120

point_list = [[0,0],[0,1],[1,1],[1,0]]
poly = geometry.Polygon(point_list)

def generate_random_points(number, poly):
	list_of_points=[]
	minx,miny,maxx,maxy = poly.bounds
	counter = 0
	while counter < number:
		p = geometry.Point(random.uniform(minx, maxx), random.uniform(miny,maxy))
		if poly.contains(p):
			list_of_points.append(p)
			counter += 1
	return list_of_points

print generate_random_points(numParticles, poly)

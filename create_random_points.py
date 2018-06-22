from shapely import geometry
import random
import numpy as np
import sys


def generate_random_points(x_array, y_array, total_points):
	point_list = zip(numpy_x, numpy_y)
	nucleus = geometry.Polygon([point for point in point_list])

	list_of_points=[]
	minx,miny,maxx,maxy = nucleus.bounds
	current_points = 0
	while current_points < total_points:
		current_point = geometry.Point(random.uniform(minx, maxx), random.uniform(miny,maxy))
		if nucleus.contains(current_point):
			list_of_points.append(current_point)
			current_points += 1

	# print(list_of_points)

	x_list = []
	y_list = []
	for point in list_of_points:
		x_list.append(point.x)
		y_list.append(point.y)

	numpy_out_x = np.array(x_list)
	numpy_out_y = np.array(y_list)

	return numpy_out_x,numpy_out_y

if __name__ == "__main__":
	numpy_x = np.array([-1,2,3,6])
	numpy_y = np.array([0,1,1,0])
	total_points = 4

	print(generate_random_points(numpy_x, numpy_y, total_points))
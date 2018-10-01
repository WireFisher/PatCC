#!/usr/bin/env python2

import random
from math import sin, cos, tan, sqrt

precision = 10000000

random.seed(2333)

def lonlat2xyz(lon, lat):
    x = cos(lon)
    y = sin(lon)
    z = sin(lat)
    return x, y, z

num_points = 10000
# random on lat-lon global grid
min_lon = 0
max_lon = 360
min_lat = -90
max_lat = 90

fp = open("lonlat_global_random.txt", "w")
for _ in xrange(num_points):
    lon = random.randint(min_lon * precision, max_lon * precision - 1) / precision
    lat = random.randint(min_lat * precision, max_lat * precision) / precision
    x, y, z = lonlat2xyz(lon, lat)
    fp.write("%.10lf, %.10lf, %.10lf\n" % (x, y, z))
fp.close()


# random on lat-lon regional grid
min_lon = 0
max_lon = 180
min_lat = -45
max_lat = 45

fp = open("lonlat_regional_random.txt", "w")
for _ in xrange(num_points):
    lon = random.randint(min_lon * precision, max_lon * precision - 1) / precision
    lat = random.randint(min_lat * precision, max_lat * precision) / precision
    x, y, z = lonlat2xyz(lon, lat)
    fp.write("%.10lf, %.10lf, %.10lf\n" % (x, y, z))
fp.close()


# random on lat-lon global grid
min_lon = 0
max_lon = 360
min_lat = -90
max_lat = 90

fp = open("lonlat_global_uniform.txt", "w")
edge_points = sqrt(num_points)
for i in xrange(edge_points):
    for j in xrange(edge_points):
    lon = min_lon + (max_lon - min_lon) * i / edge_points
    lat = min_lat + (max_lat - min_lat) * i / (edge_points-1)
    x, y, z = lonlat2xyz(lon, lat)
    fp.write("%.10lf, %.10lf, %.10lf\n" % (x, y, z))
fp.close()

#!/bin/bash

clear
mpigxx -g -DDEBUG -DUNITTEST -fopenmp -isystem /opt/intel/impi/3.2.0.011/include64 -isystem /home/yanghy/opt/opencv/include -isystem dependency/googletest/include -isystem dependency/googlemock/include -isystem /opt/netCDF-gcc4.4.7/include/ -L/opt/netCDF-gcc4.4.7/lib  -L/home/yanghy/opt/opencv/lib64 -lopencv_core -lopencv_imgproc -lopencv_imgcodecs -lnetcdf -pthread processing_unit_mgt.cxx delaunay_grid_decomposition_mgt.cxx delaunay_voronoi_2D.cxx ccpl_utils.cxx opencv_utils.cxx component.cxx unittest/*.cxx dependency/googlemock/libgmock.a dependency/googletest/libgtest.a -o run_all_test

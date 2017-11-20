#!/bin/bash

clear
mpigxx -g -DDEBUG -Wall -fopenmp -isystem /opt/intel/impi/3.2.0.011/include64 -isystem /home/yanghy/opt/opencv/include -isystem dependency/googletest/include -isystem dependency/googlemock/include -L/home/yanghy/opt/opencv/lib64 -lopencv_core -lopencv_imgproc -lopencv_imgcodecs -pthread processing_unit_mgt.cxx delaunay_grid_decomposition_mgt.cxx unittest/*.cxx dependency/googlemock/libgmock.a dependency/googletest/libgtest.a -o run_all_test

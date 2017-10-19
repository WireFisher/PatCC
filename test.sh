#!/bin/bash

clear
mpigxx -g -DDEBUG -fopenmp -isystem dependency/googletest/include -isystem dependency/googlemock/include -pthread processing_unit_mgt.cxx delaunay_grid_decomposition_mgt.cxx unittest/*.cxx dependency/googlemock/libgmock.a dependency/googletest/libgtest.a -o run_all_test

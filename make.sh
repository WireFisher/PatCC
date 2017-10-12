#!/bin/bash

clear
mpigxx -fopenmp -Idependency/googletest/include -Idependency/googlemock/include processing_unit_mgt.cxx delaunay_grid_decomposition_mgt.cxx testcase.cxx -o pDelaunay

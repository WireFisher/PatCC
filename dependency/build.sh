#!/bin/bash

GXX=g++
GTEST_DIR=$(pwd)/googletest/
GMOCK_DIR=$(pwd)/googlemock/
CPP=/opt/gcc-5.4.0/bin/g++
cd ${GTEST_DIR}
${GXX} -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv libgtest.a gtest-all.o

cd ${GMOCK_DIR}
${GXX} -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
${GXX} -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} -pthread -c ${GMOCK_DIR}/src/gmock-all.cc
ar -rv libgmock.a gtest-all.o gmock-all.o

cd ${GMOCK_DIR}/../

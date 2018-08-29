#!/bin/bash

GTEST_DIR=$(pwd)/googletest/
GMOCK_DIR=$(pwd)/googlemock/

cd ${GTEST_DIR}
g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv libgtest.a gtest-all.o

cd ${GMOCK_DIR}
g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} -pthread -c ${GMOCK_DIR}/src/gmock-all.cc
ar -rv libgmock.a gtest-all.o gmock-all.o

cd ${GMOCK_DIR}/../

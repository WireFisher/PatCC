#!/bin/bash

GXX=mpiicpc
GTEST_DIR=$(pwd)/googletest/
GMOCK_DIR=$(pwd)/googlemock/
cd ${GTEST_DIR}

rm googlemock/*.a
rm googlemock/*.o
rm googletest/*.o
rm googletest/*.a

${GXX} -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv libgtest.a gtest-all.o

cd ${GMOCK_DIR}
${GXX} -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
${GXX} -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} -pthread -c ${GMOCK_DIR}/src/gmock-all.cc
ar -rv libgmock.a gtest-all.o gmock-all.o

cd ${GMOCK_DIR}/../

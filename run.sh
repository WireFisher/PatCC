#!/bin/bash
set -x
./clean.sh
export OMP_STACKSIZE="128M"
export PDLN_DEBUG=true 
export PDLN_USE_OPENCV=true 
export PDLN_USE_NETCDF=true 
echo `/opt/gcc-5.4.0/bin/g++ --version`
rm test_log/*
optimize_option=("-O0" "-O1" "-O2" "-O3")
process_num=(1 3 4 5 8 11 16 32)
thread_num=(32 5 8 6 3 2 2 1)
for((i=0;i<4;i++));
do
	make clean
	sed -i '/^ADDED_FLAGS :=/c ADDED_FLAGS := '${optimize_option[$i]} Makefile
	make -j8 test 1> test_log/makelog$i.out 2> test_log/makelog$i.err
	for((j=0;j<8;j++));
	do
	
		echo ${process_num[$j]} ${thread_num[$j]}
		echo `/opt/gcc-5.4.0/bin/g++ --version` >> test_log/runlog$i.${process_num[$j]}.out
		OMP_STACKSIZE="128M" OMP_NUM_THREADS=${thread_num[$j]} mpiexec -n ${process_num[$j]} ./run_all_test --gtest_filter=FullProcess.* 1>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out 2>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
		echo "===================================================="
	done
done

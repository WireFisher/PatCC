#!/bin/bash

./clean.sh
export OMP_STACKSIZE="128M"
export PDLN_DEBUG=true 
export PDLN_USE_OPENCV=true 
export PDLN_USE_NETCDF=true 
source /opt/intel/composer_xe_2015.3.187/bin/compilervars_arch.sh
source /opt/intel/composer_xe_2015.3.187/bin/compilervars.sh intel64
echo `icpc --version`
rm test_log/*

optimize_option=("-O0" "-O1" "-O2" "-O3")
for((i=0;i<4;i++));
do
	make clean
	sed -i '/^ADDED_FLAGS :=/c ADDED_FLAGS := '${optimize_option[$i]} Makefile
	make -j8 test 1> test_log/makelog$i.out 2> test_log/makelog$i.err
	for process_num in 1 4 5 8 16 32
	do
		echo `pwd`
		thread_num=`expr 32 / $process_num`
		echo `mpiicpc --version` >> test_log/runlog$i.$process_num.out
		OMP_STACKSIZE="128M" OMP_NUM_THREADS=$thread_num mpiexec -n $process_num ./run_all_test --gtest_filter=FullProcess.* 1>>test_log/runlog$i.$process_num.$thread_num.out 2>>test_log/runlog$i.$process_num.$thread_num.err
		mv log log.$i.$process_num
		mkdir log
		echo "===================================================="
	done
done

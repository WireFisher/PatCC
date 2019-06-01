#!/bin/bash
set -x
export PDLN_DEBUG=true 
export PDLN_USE_OPENCV=true 
export PDLN_USE_NETCDF=true 
#source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64  #icpc15
source /opt/gcc-4.4.7/bin/gcc-env
optimize_option=( "-O3" )
process_num=(1  1 3 4 5 8 11 16 32 64 64 64 64 64 )
thread_num=(10 32 5 8 6 3  2  2  1  1  2  4  8 13 )
for((i=0;i<${#optimize_option[@]};i++));
do
	for((j=0;j<${#process_num[@]};j++));
	do
		total=$(expr ${thread_num[$j]} \* ${process_num[$j]})
		OMP_STACKSIZE="128M" OMP_NUM_THREADS=${thread_num[$j]} mpiexec -n ${process_num[$j]} ./run_all_test MonteCarlo_1000000.dat --gtest_filter=FullProcess.P*
	    echo $file $total `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> all_grid_checksum
	    echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
	done
done

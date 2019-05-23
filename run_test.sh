#!/bin/bash
set -x
source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64  #icpc15
#source /opt/intel/composer_xe_2013_sp1.2.144/bin/compilervars.sh intel64 #icpc14
#source /opt/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64 #icpc13
echo `mpiicpc --version`
echo $I_MPI_ROOT

rm nohup.out
rm -rf test_log
rm -rf log
mkdir test_log
mkdir -p log/check
rm all_grid_checksum

optimize_option=( "-O3" )
process_num=(1  1 3 4 5 8 11 16 32 64 64 64 64 64 )
thread_num=(10 32 5 8 6 3  2  2  1  1  2  4  8 13 )
old_test_list=( "OLD_Basic" "OLD_LatLonGrid" "OLD_LatLonSinglePolar" "OLD_LatLonMutiPolars" "OLD_ThreePolar" "OLD_ThreePolarBig" )
if [ "${#process_num[@]}" -ne "${#thread_num[@]}" ]; then
	echo "error the length of process_num and thread_num is different"
	exit 1
fi

for((i=0;i<${#optimize_option[@]};i++));
do
	make clean
	sed -i '/^ADDED_FLAGS :=/c ADDED_FLAGS := '${optimize_option[$i]} Makefile
	make -j8 test 1> test_log/makelog$i.out 2> test_log/makelog$i.err
	for((j=0;j<${#process_num[@]};j++));
	do
		./clean.sh
		mkdir log
		mkdir log/check
		echo `mpiicpc --version`
		total=$(expr ${thread_num[$j]} \* ${process_num[$j]})

		for((k=0;k<${#old_test_list[@]};k++));
		do
            if [ -f log/global_triangles_$total ]; then
	            rm log/global_triangles_$total
	        fi
			echo ${old_test_list[$k]} >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
            echo ${old_test_list[$k]} >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
	        OMP_NUM_THREADS=${thread_num[$j]} mpiexec -n ${process_num[$j]} ./run_all_test --gtest_filter=FullProcess.${old_test_list[$k]} 1>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out 2>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
	        echo ${old_test_list[$k]} `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> all_grid_checksum
	        echo ${old_test_list[$k]} `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
	        mv log/global_triangles_$total log/check/global_triangles_${old_test_list[$k]}
		done

		for file in `cat many_grid_list`
		do
		    if [ -f log/global_triangles_$total ]; then
		         rm log/global_triangles_$total
	        fi
			echo $file >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
			echo $file >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
			OMP_NUM_THREADS=${thread_num[$j]} mpiexec -n ${process_num[$j]} ./run_all_test ${file} --gtest_filter=FullProcess.M* 1>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out 2>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
		    echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> all_grid_checksum
		    echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
		    mv log/global_triangles_$total log/check/global_triangles_$file
		done

		for file in `cat performance_grid_list`
		do
            if [ -f log/global_triangles_$total ]; then
	            rm log/global_triangles_$total
	        fi
			echo $file >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
			echo $file >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
			OMP_NUM_THREADS=${thread_num[$j]} mpiexec -n ${process_num[$j]} ./run_all_test ${file} --gtest_filter=FullProcess.P* 1>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out 2>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
    		echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> all_grid_checksum
		    echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
		    mv log/global_triangles_$total log/check/global_triangles_$file
		done

		for file in `cat coupler_grid_list`
		do
            if [ -f log/global_triangles_$total ]; then
                rm log/global_triangles_$total
	        fi
			echo $file >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
			echo $file >>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
			OMP_NUM_THREADS=${thread_num[$j]} mpiexec -n ${process_num[$j]} ./run_all_test ${file} --gtest_filter=FullProcess.R* 1>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out 2>>test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.err
		    echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> all_grid_checksum
		    echo $file `md5sum log/global_triangles_$total | awk -F" " '{print $1}'` >> test_log/runlog$i.${process_num[$j]}.${thread_num[$j]}.out
		    mv log/global_triangles_$total log/check/global_triangles_$file
		done
		mv log test_log/log${i}_${process_num[$j]}_${thread_num[$j]}
		mv all_grid_checksum test_log/all_grid_checksum${i}_${process_num[$j]}_${thread_num[$j]}
		echo "===================================================="
	done
done

#!/usr/bin/bash

make
echo "#Strong scaling"
thr=1
export OMP_NUM_THREADS=1
echo "#NUM PROCS $thr"
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
let thr=thr*24
while [ $thr -lt 49 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpirun --bind-to-core -np $thr	./numatest_mpi 4294967296 &> numa_ss_np4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpirun -oversubscribe --bind-to-core -np $thr	./numatest_mpi 4294967296 &> numa_ss_np4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done
let thr=1
mpirun --bind-to-core -np $thr ./numatest_mpi 4294967296 &> numa_ss_np4_mm_$thr

echo "#Weak scaling"
thr=1
size=44739243
size1=44739243
export OMP_NUM_THREADS=1
echo "#NUM PROCS $thr"
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
let thr=thr*24
while [ $thr -lt 49 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpirun --bind-to-core -np $thr	./numatest_mpi $size &> numa_ws_np4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpirun -oversubscribe --bind-to-core -np $thr	./numatest_mpi $size &> numa_ws_np4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done
let size=44739243
let thr=1
mpirun --bind-to-core -np $thr ./numatest_mpi $size &> numa_ws_np4_mm_$thr

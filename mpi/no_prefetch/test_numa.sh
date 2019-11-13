#!/usr/bin/bash

make
echo "#Strong scaling"
thr=1
export OMP_NUM_THREADS=1
mpirun -np $thr ./numatest_mpi 17179869184
#cp sicm_numa_config sicm_numa_config_o1
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
mpirun -np $thr	./numatest_mpi 17179869184
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done

echo "#Weak scaling"
thr=1
size=178956971
size1=178956971
export OMP_NUM_THREADS=1
mpirun -np $thr ./numatest_mpi $size
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
mpirun -np $thr	./numatest_mpi $size
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done

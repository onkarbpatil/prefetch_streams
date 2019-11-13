#!/usr/bin/bash

make
echo "#Strong scaling"
thr=1
export OMP_NUM_THREADS=1
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=2147483648 stream_mpi.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to core -np $thr ./stream_mpi 
#cp sicm_numa_config sicm_numa_config_o1
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=2147483648 stream_mpi.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to core -np $thr	./stream_mpi 
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done

echo "#Weak scaling"
thr=1
size=22369622
size1=22369622
export OMP_NUM_THREADS=1
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=$size stream_mpi.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to core -np $thr ./stream_mpi $size
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=$size stream_mpi.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to core -np $thr	./stream_mpi 
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done

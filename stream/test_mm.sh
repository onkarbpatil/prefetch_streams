#!/usr/bin/bash

make
echo "#Strong scaling"
thr=1
export OMP_NUM_THREADS=1
echo "#NUM PROCS $thr"
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=536870912 stream_mpi.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr ./stream_mpi &> stream_ss_pre4_mm_$thr  
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
let thr=thr*24
while [ $thr -lt 49 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=536870912 stream_mpi.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr ./stream_mpi &> stream_ss_pre4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=536870912 stream_mpi.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to-core -np $thr	./stream_mpi &> stream_ss_pre4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done

echo "#Weak scaling"
thr=1
size=5592406
size1=5592406
export OMP_NUM_THREADS=1
echo "#NUM PROCS $thr"
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=$size stream_mpi.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr ./stream_mpi $size &> stream_ws_pre4_mm_$thr
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
let thr=thr*24
while [ $thr -lt 49 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=$size stream_mpi.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr	./stream_mpi &> stream_ws_pre4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi -DSTREAM_ARRAY_SIZE=$size stream_mpi.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to-core -np $thr	./stream_mpi &> stream_ws_pre4_mm_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done

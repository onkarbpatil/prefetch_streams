#!/usr/bin/bash

make
echo "#Strong scaling"
thr=1
export OMP_NUM_THREADS=1
echo "#NUM PROCS $thr"
mpicc -o stream_mpi_np -DSTREAM_ARRAY_SIZE=2147483648 stream_mpi_np.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr ./stream_mpi_np >> stream_ss_np_$thr 
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
let thr=thr*24
while [ $thr -lt 49 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi_np -DSTREAM_ARRAY_SIZE=2147483648 stream_mpi_np.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr ./stream_mpi_np >> stream_ss_np_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi_np -DSTREAM_ARRAY_SIZE=2147483648 stream_mpi_np.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to-core -np $thr	./stream_mpi_np >> stream_ss_np_$thr 
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let thr=thr+24
done

echo "#Weak scaling"
thr=1
size=22369622
size1=22369622
export OMP_NUM_THREADS=1
echo "#NUM PROCS $thr"
mpicc -o stream_mpi_np -DSTREAM_ARRAY_SIZE=$size stream_mpi_np.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr ./stream_mpi_np $size >> stream_ws_np_$thr
#cp sicm_numa_config sicm_numa_config_o1
let size=size*24
let thr=thr*24
while [ $thr -lt 49 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi_np -DSTREAM_ARRAY_SIZE=$size stream_mpi_np.c -O1 -lnuma -lm
mpirun --bind-to-core -np $thr	./stream_mpi_np >> stream_ws_np_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done
while [ $thr -lt 97 ]; do
#	export OMP_NUM_THREADS=$thr
echo "#NUM PROCS $thr"
mpicc -o stream_mpi_np -DSTREAM_ARRAY_SIZE=$size stream_mpi_np.c -O1 -lnuma -lm
mpirun -oversubscribe  --bind-to-core -np $thr	./stream_mpi_np >> stream_ws_np_$thr
#	cp sicm_numa_config "sicm_numa_config_o$thr"
	let size=size+24*size1
	let thr=thr+24
done

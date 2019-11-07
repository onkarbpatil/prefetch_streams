#!/usr/bin/bash

thr=8
node_range=0-23
export OMP_NUM_THREADS=1
echo "0-23 Upto 24 threads\n" >> memory_characterization
echo "1:\n" >> memory_characterization
numactl -C $node_range ./numatest >> memory_characterization
while [ $thr -lt 25 ]; do
	export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
	numactl -C $node_range ./numatest >> memory_characterization
	let thr=thr+8
done
thr=8
node_range=24-47 
export OMP_NUM_THREADS=1
echo "24-47 Upto 24 threads\n" >> memory_characterization
echo "1:\n" >> memory_characterization
numactl -C $node_range ./numatest >> memory_characterization
while [ $thr -lt 25 ]; do
        export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=8
node_range=48-71 
export OMP_NUM_THREADS=1
echo "48-71 Upto 24 threads\n" >> memory_characterization
echo "1:\n" >> memory_characterization
numactl -C $node_range ./numatest >> memory_characterization
while [ $thr -lt 25 ]; do
        export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=8
node_range=72-95 
export OMP_NUM_THREADS=1
echo "72-95 Upto 24 threads\n" >> memory_characterization
echo "1:\n" >> memory_characterization
numactl -C $node_range ./numatest >> memory_characterization
while [ $thr -lt 25 ]; do
        export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done


thr=32
node_range="0-23,48-71"
echo "0-23 48-71 Upto 48 threads\n" >> memory_characterization
while [ $thr -lt 49 ]; do
        export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=32
node_range="24-47,72-95"
echo "24-47 72-95 Upto 48 threads\n" >> memory_characterization
while [ $thr -lt 49 ]; do
        export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done

thr=32
node_range=0-47
echo "0-47 Upto 48 threads\n" >> memory_characterization
while [ $thr -lt 49 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=32
node_range=24-71
echo "24-71 Upto 48 threads\n" >> memory_characterization
while [ $thr -lt 49 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=32
node_range=48-95
echo "48-95 Upto 48 threads\n" >> memory_characterization
while [ $thr -lt 49 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done



node_range=0-71
echo "0-71 Upto 72 threads\n" >> memory_characterization
while [ $thr -lt 73 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=56
node_range=24-95
echo "24-95 Upto 72 threads\n" >> memory_characterization
while [ $thr -lt 73 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=56
node_range="0-23,48-95"
echo "0-23 48-95 Upto 72 threads\n" >> memory_characterization
while [ $thr -lt 73 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done
thr=56
node_range="0-47,72-95"
echo "0-47 72-95 Upto 72 threads\n" >> memory_characterization
while [ $thr -lt 73 ]; do
        export OMP_NUM_THREADS=$thr
        echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done


echo "0-96 Upto 96 threads\n" >> memory_characterization
node_range=0-95
while [ $thr -lt 97 ]; do
        export OMP_NUM_THREADS=$thr
	echo "$thr:\n" >> memory_characterization
        numactl -C $node_range ./numatest >> memory_characterization
        let thr=thr+8
done


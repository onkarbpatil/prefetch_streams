#! /usr/bin/bash

#sh mpi/test_numa.sh
#cd mpi/no_prefetch
#sh test_numa.sh
cd stream
sh test_numa.sh
cd no_prefetch
sh test_numa.sh

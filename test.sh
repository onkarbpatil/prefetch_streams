#! /usr/bin/bash

cd mpi
sh test_numa.sh
cd no_prefetch
sh test_numa.sh
cd ../../stream
sh test_numa.sh
cd no_prefetch
sh test_numa.sh

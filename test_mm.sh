#! /usr/bin/bash

cd mpi
sh test_mm.sh
cd no_prefetch
sh test_mm.sh
cd ../../stream
sh test_mm.sh
cd no_prefetch
sh test_mm.sh

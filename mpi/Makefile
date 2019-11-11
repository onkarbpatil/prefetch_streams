CC = mpicc
LDFLAGS = -O1 -fopenmp -lnuma -lm

numatest_mpi: test3.c numatest_mpi.c numatest_mpi.h
		$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
	clean:
		rm -f $(obj) numatest_mpi

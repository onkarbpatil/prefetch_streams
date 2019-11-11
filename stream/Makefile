CC = mpicc
LDFLAGS = -O1 -fopenmp -lnuma -lm

stream_mpi: stream_mpi.c
		$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
	clean:
		rm -f $(obj) stream_mpi

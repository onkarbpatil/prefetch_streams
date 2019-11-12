#include "numatest_mpi.h"

int main(int argc, char ** argv){
		int rank, size;
		unsigned long bytes;
		char * ptr;
	MPI_Init(&argc, &argv);
		bytes = strtoul(argv[1], &ptr, 10);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	char *labels[] = {"NVME","DRAM"};
	numatest(2,labels, rank, size, bytes);
	MPI_Finalize();
}

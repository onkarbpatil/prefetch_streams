#include "numatest_omp.h"

int main(int argc, char ** argv){
	char *labels[] = {"NVME","DRAM"};
	numatest(2,labels);
}

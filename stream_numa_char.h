#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <sched.h>
#include <sys/time.h>
#include "numa.h"
#include <mpi.h>

#define BILLION  1000000000L;

struct numa_node_bw{
	int numa_id;
	char * mem_type;
	long double wr_only_avg;
	long double owor_avg;
        long double owtr_avg;
        long double owthr_avg;
        long double owfr_avg;
        long double twor_avg;
        long double twtr_avg;
        long double twthr_avg;
        long double twfr_avg;
        long double thwor_avg;
        long double thwtr_avg;
        long double thwthr_avg;
        long double thwfr_avg;
        long double fwor_avg;
        long double fwtr_avg;
        long double fwthr_avg;
        long double fwfr_avg;
        long double str_avg;
        long double rand_avg;
        long double diff_avg;
        long double row_avg;
        long double col_avg;
        long double rc_avg;
        long double l2cache_avg;
        long double t_sten_avg;
        long double f_sten_avg;
        long double n_sten_avg;
	struct numa_node_bw * next;
};

struct stream_ranking{
	int * rank;
	long double * bw;
};

struct stream_ranking stream_priority[27];
struct numa_node_bw * numa_node_list;
struct numa_node_bw * numa_list_head;
int mem_types;
int max_node;
int numt;
int total_numa_nodes;
int * numa_node_ids;
struct bitmask * numa_nodes;
char ** mem_tech;
long double * means;
int * cluster_sizes;

void calculate_distances();
void calculate_mean();
void classify();
void sort_list(struct numa_node_bw * new_node);
void write_config_file();
void sort_stream_priority();
void stream_numa_char_omp(int argc, char ** argv);
void stream_numa_char_mpi(int argc, char ** argv);

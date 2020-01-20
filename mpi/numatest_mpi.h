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
		long double s_sten_avg;
		long double t7_sten_avg;
		long double wr_only_t;
		long double l2cache_t;
		long double rand_t;
		long double str_t;
		long double owfr_t;
		long double t_sten_t;
		long double f_sten_t;
		long double n_sten_t;
		long double s_sten_t;
		long double t7_sten_t;
		long double wr_only_tmin;
		long double l2cache_tmin;
		long double rand_tmin;
		long double str_tmin;
		long double owfr_tmin;
		long double t_sten_tmin;
		long double f_sten_tmin;
		long double n_sten_tmin;
		long double s_sten_tmin;
		long double t7_sten_tmin;
		long double wr_only_min;
		long double l2cache_min;
		long double rand_min;
		long double str_min;
		long double owfr_min;
		long double t_sten_min;
		long double f_sten_min;
		long double n_sten_min;
		long double s_sten_min;
		long double t7_sten_min;
		long double wr_only_o;
		long double l2cache_o;
		long double rand_o;
		long double str_o;
		long double owfr_o;
		long double t_sten_o;
		long double f_sten_o;
		long double n_sten_o;
		long double s_sten_o;
		long double t7_sten_o;
		long double wr_only_omin;
		long double l2cache_omin;
		long double rand_omin;
		long double str_omin;
		long double owfr_omin;
		long double t_sten_omin;
		long double f_sten_omin;
		long double n_sten_omin;
		long double s_sten_omin;
		long double t7_sten_omin;
	struct numa_node_bw * next;
};

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
void numatest(int argc, char ** argv, int rank, int procs, unsigned long bytes);

#include "numatest_mpi.h"

struct numa_node_bw * numa_node_list = NULL;
struct numa_node_bw * numa_list_head = NULL;
int mem_types;
int max_node;
int numt;
int total_numa_nodes = 0;
int * numa_node_ids;
struct bitmask * numa_nodes;
char ** mem_tech;
long double * means;
int * cluster_sizes;

void calculate_distances(){
	int i;
        struct numa_node_bw * bw_it = numa_list_head;
	while(bw_it != NULL){
		i = 0;
		long double delta = 9999999.9999;
		while(i < mem_types){
			long double dist = abs(sqrt(abs((means[i] - bw_it->owfr_avg))*abs((means[i] - bw_it->owfr_avg))));
			if(dist < delta){
				delta = dist;
				if(strcmp(bw_it->mem_type, mem_tech[i])!=0){
					if(((i-1)>=0)&&(strcmp(bw_it->mem_type, mem_tech[i-1])==0)){
						cluster_sizes[i-1]--;
						bw_it->mem_type = mem_tech[i];
						cluster_sizes[i]++;
					}
					else if(((i+1)<mem_types)&&(strcmp(bw_it->mem_type, mem_tech[i+1])==0)){
                                                cluster_sizes[i+1]--;
						bw_it->mem_type = mem_tech[i];
						cluster_sizes[i]++;
                    			}
				}
			}
			i++;
		}
		bw_it = bw_it->next;
	}
}

void calculate_mean(){
	int i = 0;
	struct numa_node_bw * bw_it = numa_list_head;
	while(i < mem_types){
		int j = 0;
		means[i] = 0.0;
		while(j < cluster_sizes[i]){
			means[i] += bw_it->owfr_avg;
			j++;
			bw_it = bw_it->next;
		}
		means[i] /= cluster_sizes[i];
		i++;
	}
	calculate_distances();
}

void classify(){
	int cluster_size;
	int last_cluster_size;
	cluster_size = total_numa_nodes/mem_types;
	last_cluster_size = cluster_size + (total_numa_nodes%mem_types);
	cluster_sizes = (int *)malloc(sizeof(int)*mem_types);
	means = (long double *)(malloc(mem_types*sizeof(long double)));
	struct numa_node_bw * bw_it = numa_list_head;
	int i = 0;
	while(i < mem_types){
		if(i == (mem_types - 1)){
			cluster_sizes[i] = last_cluster_size;
		}
		else{
			cluster_sizes[i] = cluster_size;
		}
		i++;
	}
	i = 0;
	int j = 1;
	while(bw_it != NULL){
		bw_it->mem_type = mem_tech[i];
		if(j < cluster_sizes[i]){
			j++;
		}
		else{
			j = 1;
			i++;
		}
		bw_it = bw_it->next;
	}
	bw_it = numa_list_head;
	i = 0;
	while(i < 10){
		calculate_mean();
		i++;
	}

}

void sort_list(struct numa_node_bw * new_node){
	struct numa_node_bw * bw_it = numa_list_head;
	struct numa_node_bw * prev_bw_it = NULL;
	while(bw_it != NULL){
		if((bw_it->owtr_avg > new_node->owtr_avg)){
			if(prev_bw_it == NULL){
				new_node->next = bw_it;
				numa_list_head = new_node;
			}else{
				prev_bw_it->next = new_node;
				new_node->next = bw_it;
			}
			return;
		}
		prev_bw_it = bw_it;
		bw_it = bw_it->next;
	}
	prev_bw_it->next = new_node;
	return;

}

void write_config_file(){
	FILE * conf;
	char fname[50];
	char thr[10];
	snprintf(thr, 10, "%d", numt);
	strcpy(fname, "sicm_numa_config");
	strcat(fname, thr);
	conf = fopen(fname, "w");
	struct numa_node_bw * bw_it = numa_list_head;
		printf("#NUMA id WR-only_avg_bw 1W4R_avg_bw Str_avg_bw Rand_avg_bw LLC_avg_bw WR-only_pk_bw 1W4R_pk_bw Str_pk_bw Rand_pk_bw LLC_pk_bw WR-only_avg_lat 1W4R_avg_lat Str_avg_lat Rand_avg_lat LLC_avg_lat WR-only_min 1W4R_min Str_min Rand_min LLC_min WR-only_avg_ovr 1W4R_avg_ovr Str_avg_ovr Rand_avg_ovr LLC_avg_ovr WR-only_omin 1W4R_omin Str_omin Rand_omin LLC_omin\n");
	while(bw_it != NULL){	
		fprintf(conf, "%d %s %Lf %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF\n", bw_it->numa_id, bw_it->mem_type, bw_it->wr_only_avg, bw_it->owor_avg, bw_it->owtr_avg, bw_it->owthr_avg, bw_it->owfr_avg, bw_it->twor_avg, bw_it->twtr_avg, bw_it->twthr_avg, bw_it->twfr_avg, bw_it->thwor_avg, bw_it->thwtr_avg, bw_it->thwthr_avg, bw_it->thwfr_avg, bw_it->fwor_avg, bw_it->fwtr_avg, bw_it->fwthr_avg, bw_it->fwfr_avg, bw_it->str_avg, bw_it->rand_avg, bw_it->diff_avg, bw_it->row_avg, bw_it->col_avg, bw_it->rc_avg, bw_it->t_sten_avg, bw_it->f_sten_avg, bw_it->n_sten_avg, bw_it->l2cache_avg);
		printf("%d %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF\n", bw_it->numa_id, bw_it->wr_only_avg, bw_it->owfr_avg, bw_it->str_avg, bw_it->rand_avg, bw_it->l2cache_avg, bw_it->wr_only_min, bw_it->owfr_min, bw_it->str_min, bw_it->rand_min, bw_it->l2cache_min, bw_it->wr_only_t, bw_it->owfr_t, bw_it->str_t, bw_it->rand_t, bw_it->l2cache_t, bw_it->wr_only_tmin, bw_it->owfr_tmin, bw_it->str_tmin, bw_it->rand_tmin, bw_it->l2cache_tmin, bw_it->wr_only_o, bw_it->owfr_o, bw_it->str_o, bw_it->rand_o, bw_it->l2cache_o, bw_it->wr_only_omin, bw_it->owfr_omin, bw_it->str_omin, bw_it->rand_omin, bw_it->l2cache_omin);
		//printf("%d %s %Lf %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF %LF\n", bw_it->numa_id, bw_it->mem_type, bw_it->wr_only_avg, bw_it->owor_avg, bw_it->owtr_avg, bw_it->owthr_avg, bw_it->owfr_avg, bw_it->twor_avg, bw_it->twtr_avg, bw_it->twthr_avg, bw_it->twfr_avg, bw_it->thwor_avg, bw_it->thwtr_avg, bw_it->thwthr_avg, bw_it->thwfr_avg, bw_it->fwor_avg, bw_it->fwtr_avg, bw_it->fwthr_avg, bw_it->fwfr_avg, bw_it->str_avg, bw_it->rand_avg, bw_it->diff_avg, bw_it->row_avg, bw_it->col_avg, bw_it->rc_avg, bw_it->t_sten_avg, bw_it->f_sten_avg, bw_it->n_sten_avg, bw_it->l2cache_avg);
		bw_it = bw_it->next;
	}
	fclose(conf);
}

void numatest(int argc, char ** argv, int rank, int procs, unsigned long bytes){
	max_node = numa_max_node() + 1;
	int cpu_count = numa_num_possible_cpus();
	numa_node_ids = (int*)malloc(sizeof(int)*max_node);
	struct bitmask * numa_nodes = numa_get_membind();
	int i = 0;
        while(i < numa_nodes->size){
                if(numa_bitmask_isbitset(numa_nodes, i)){
                        numa_node_ids[total_numa_nodes] = i;
                        total_numa_nodes++;
                }
                i++;
        }
		total_numa_nodes++;
		MPI_Request reqs[32768]; 
		MPI_Status stat[32768];
	unsigned long size = bytes/procs;
	int mbs = size/sizeof(double);
	int r_size = 32768;
	int c_size = 32768;
	if(procs != 1){
		r_size = (32768)/(12) + 2*sizeof(double);
		c_size = (32768)/(procs/12) + 2*sizeof(double);
	}
	int rs = 0;
	int z = 0;
	int bdist = 0;
	int sdist = 0;
	int rd_dist, wr_dist;
	int * rand_tab;
	int inner = 0;
	rand_tab = (int*)malloc(mbs*sizeof(int));
	double *a, *b, *c, *d, *e, *f, *g, *h;
	double **aa, **bb, **cc;
	clock_t start, end;
	struct timespec begin, stop, obegin, ostop;
	srand(10725);
	//sleep(10);
	if(argc == 0){
		printf("Enter memory technologies available in ascending order of speed. eg: GPU NVRAM DRAM HBM\n");
		return;
	}
	else{
		mem_types = argc;
		mem_tech = (char**)malloc(argc*sizeof(char*));
		int a;
		for(a = 0; a < argc; a++){
			mem_tech[a] = argv[a];
		}

	}
	for(i=0; i< size/sizeof(double); i++)
		rand_tab[i]=rand()%(size/sizeof(double));

//#ifdef _OPENMP
//#pragma omp parallel private(numt)
    {
    numt = omp_get_num_threads();
    }
//#endif
  	i = 0;
	while(i < total_numa_nodes){
			if(i <= 1){
				rd_dist = 8192/sizeof(double);
				wr_dist = 2048/sizeof(double);
				bdist = rd_dist;
				sdist = wr_dist;
				if(procs > 48){
					rd_dist = 0;
					wr_dist = 0;
					bdist = 0;
					sdist = 0;
				}
			}
			if((i > 1)&&(i < (total_numa_nodes-1))){
				rd_dist = 8192/sizeof(double);
				wr_dist = 32768/sizeof(double);
				bdist = wr_dist;
				sdist = rd_dist;
				if(procs > 48){
					rd_dist = 0;
					wr_dist = 32768/sizeof(double);
					bdist = wr_dist;
					sdist = rd_dist;
				}
			}
			if(i == (total_numa_nodes - 1)){
				rd_dist = 8192/sizeof(double);
				wr_dist = 2048/sizeof(double);
				bdist = rd_dist;
				sdist = wr_dist;
				if(procs > 48){
					rd_dist = 0;
					wr_dist = 0;
					bdist = 0;
					sdist = 0;
				}
			}
	// Dynamically allocate the three arrays using "posix_memalign()"
		int iters = 0;
		int stride;
		long double wr_only_avg=0.0;
		long double owor_avg=0.0;
		long double owtr_avg=0.0;
        long double owthr_avg=0.0;
		long double owfr_avg=0.0;
		long double twor_avg=0.0;
		long double twtr_avg=0.0;
		long double twthr_avg=0.0;
		long double twfr_avg=0.0;
		long double thwor_avg=0.0;
		long double thwtr_avg=0.0;
		long double thwthr_avg=0.0;
		long double thwfr_avg=0.0;
		long double fwor_avg=0.0;
		long double fwtr_avg=0.0;
		long double fwthr_avg=0.0;
		long double fwfr_avg=0.0;
		long double str_avg=0.0;
		long double rand_avg = 0.0;
		long double diff_avg = 0.0;
		long double row_avg = 0.0;
		long double col_avg = 0.0;
		long double rc_avg = 0.0;
		long double l2cache_avg = 0.0;
		long double t_sten_avg = 0.0;
		long double f_sten_avg = 0.0;
		long double n_sten_avg = 0.0;
		long double wr_only_t = 0.0;
		long double owfr_t = 0.0;
		long double l2cache_t = 0.0;
		long double str_t = 0.0;
		long double rand_t = 0.0;
		long double wr_only_tmin = 999999.9999;
		long double owfr_tmin = 999999.9999;
		long double l2cache_tmin = 999999.9999;
		long double str_tmin = 999999.9999;
		long double rand_tmin = 999999.9999;
		long double wr_only_min = 0.0;
		long double owfr_min = 0.0;
		long double l2cache_min = 0.0;
		long double str_min = 0.0;
		long double rand_min = 0.0;
		long double wr_only_o = 0.0;
		long double owfr_o = 0.0;
		long double l2cache_o = 0.0;
		long double str_o = 0.0;
		long double rand_o = 0.0;
		long double wr_only_omin = 999999.9999;
		long double owfr_omin = 999999.9999;
		long double l2cache_omin = 999999.9999;
		long double str_omin = 999999.9999;
		long double rand_omin = 999999.9999;
		long double accum;
		for( iters = 0; iters < 10; iters++)
		{
			int j = 0;
			int k = 0;
			if(i == (total_numa_nodes-1)){
				a = (double*)numa_alloc_onnode(size, numa_node_ids[0]);
				b = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				c = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				d = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				e = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				//f = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				//g = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				//h = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
			/*	aa = (double**)numa_alloc_onnode(r_size, numa_node_ids[0]);
				for(j = 0; j < r_size/sizeof(double*); j++){
					aa[j] = (double*)numa_alloc_onnode(c_size, numa_node_ids[0]);
				}
				bb = (double**)numa_alloc_onnode(r_size, numa_node_ids[2]);
				for(j = 0; j < r_size/sizeof(double*); j++){
					bb[j] = (double*)numa_alloc_onnode(c_size, numa_node_ids[2]);
				}
				cc = (double**)numa_alloc_onnode(r_size, numa_node_ids[2]);
				for(j = 0; j < r_size/sizeof(double*); j++){
					cc[j] = (double*)numa_alloc_onnode(c_size, numa_node_ids[2]);
				}*/
			}else{
				a = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				b = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				c = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				d = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				e = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				//f = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				//g = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				//h = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				/*aa = (double**)numa_alloc_onnode(r_size, numa_node_ids[i]);
				for(j = 0; j < r_size/sizeof(double*); j++){
					aa[j] = (double*)numa_alloc_onnode(c_size, numa_node_ids[i]);
				}
				bb = (double**)numa_alloc_onnode(r_size, numa_node_ids[i]);
				for(j = 0; j < r_size/sizeof(double*); j++){
					bb[j] = (double*)numa_alloc_onnode(c_size, numa_node_ids[i]);
				}
				cc = (double**)numa_alloc_onnode(r_size, numa_node_ids[i]);
				for(j = 0; j < r_size/sizeof(double*); j++){
					cc[j] = (double*)numa_alloc_onnode(c_size, numa_node_ids[i]);
				}*/
			}
			long double empty=0.0;
			long double empty2=0.0;

/*			for(j =0; j < (r_size/sizeof(double*)); j++){
         		for(k = 0; k < (c_size/sizeof(double)); k++)
					aa[j][k] = (double)rand();
			}*/
redo1:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
			for(j = 0;j < wr_dist;j++){
					__builtin_prefetch (&a[j], 1, 0);
					__builtin_prefetch (&b[j], 1, 0);
					__builtin_prefetch (&c[j], 1, 0);
					__builtin_prefetch (&d[j], 1, 0);
					__builtin_prefetch (&e[j], 1, 0);
				//	__builtin_prefetch (&f[j], 1, 0);
				//	__builtin_prefetch (&g[j], 1, 0);
				//	__builtin_prefetch (&h[j], 1, 0);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			if(rank == 0){
			accum = ( ostop.tv_sec - obegin.tv_sec ) + (long double)( ostop.tv_nsec - obegin.tv_nsec ) / (long double)BILLION;
			wr_only_o += accum;
			if(accum < wr_only_omin)
					wr_only_omin = accum;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(j = 0;j < ((size/sizeof(double)) - wr_dist);j+=6){
					__builtin_prefetch (&a[j+wr_dist], 1, 0);
					__builtin_prefetch (&b[j+wr_dist], 1, 0);
					__builtin_prefetch (&c[j+wr_dist], 1, 0);
					__builtin_prefetch (&d[j+wr_dist], 1, 0);
					__builtin_prefetch (&e[j+wr_dist], 1, 0);
					//__builtin_prefetch (&f[j+wr_dist], 1, 0);
					//__builtin_prefetch (&g[j+wr_dist], 1, 0);
					//__builtin_prefetch (&h[j+wr_dist], 1, 0);
					for(inner = 0;inner<6;inner++){
						a[j+inner] = 1.0;
						b[j+inner] = 2.0;
						c[j+inner] = 3.0;
						d[j+inner] = 4.0;
						e[j+inner] = 5.0;
					}
			//	f[j] = 6.0;
			//	g[j] = 7.0;
			//	h[j] = 8.0;
			}
			for(j = ((size/sizeof(double)) - wr_dist); j < (size/sizeof(double)); j++){
				a[j] = 1.0;
				b[j] = 2.0;
				c[j] = 3.0;
				d[j] = 4.0;
				e[j] = 5.0;
			//	f[j] = 6.0;
			//	g[j] = 7.0;
			//	h[j] = 8.0;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			wr_only_t += accum;
			if(accum < wr_only_tmin)
					wr_only_tmin = accum;
			if(accum <= empty){
				goto redo1;
			}
			wr_only_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
/*redo2:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = 20 + b[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo2;
			}
			owor_avg += ((2*size*1.0E-06)/(long double)(accum - empty));
			}
redo3:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] + b[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo3;
			}
			owtr_avg += ((3*size*1.0E-06)/(long double)(accum - empty));
			}
redo4:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] + d[j] + b[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo4;
			}
			owthr_avg += ((4*size*1.0E-06)/(long double)(accum - empty));
			}*/
redo5:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
			if(wr_dist >= rd_dist){
				for(j =0; j < wr_dist - rd_dist; j++){
					__builtin_prefetch (&a[j], 1, 0);
				}
				for(j =(wr_dist - rd_dist); j < wr_dist; j++){
					__builtin_prefetch (&a[j], 1, 0);
					__builtin_prefetch (&b[j-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&c[j-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&d[j-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&e[j-(wr_dist - rd_dist)], 0, 0);
				}
			}else{
				for(j =0; j < rd_dist - wr_dist; j++){
					__builtin_prefetch (&b[j], 0, 0);
					__builtin_prefetch (&c[j], 0, 0);
					__builtin_prefetch (&d[j], 0, 0);
					__builtin_prefetch (&e[j], 0, 0);
				}
				for(j =(rd_dist - wr_dist); j < rd_dist; j++){
					__builtin_prefetch (&a[j-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[j], 0, 0);
					__builtin_prefetch (&c[j], 0, 0);
					__builtin_prefetch (&d[j], 0, 0);
					__builtin_prefetch (&e[j], 0, 0);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			if(rank == 0){
			accum = ( ostop.tv_sec - obegin.tv_sec ) + (long double)( ostop.tv_nsec - obegin.tv_nsec ) / (long double)BILLION;
			owfr_o += accum;
			if(accum < owfr_omin)
					owfr_omin = accum;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(j =0; j < ((size/sizeof(double))-bdist); j+=4){
					__builtin_prefetch (&a[j+wr_dist], 1, 0);
					__builtin_prefetch (&b[j+rd_dist], 0, 0);
					__builtin_prefetch (&c[j+rd_dist], 0, 0);
					__builtin_prefetch (&d[j+rd_dist], 0, 0);
					__builtin_prefetch (&e[j+rd_dist], 0, 0);
						for(inner=0;inner<4;inner++)
                            a[j+inner] = c[j+inner] + d[j+inner] + e[j+inner] + b[j+inner];
            }
			if(wr_dist >= rd_dist){
				for(j = ((size/sizeof(double))-wr_dist);j <((size/sizeof(double))-rd_dist);j+=4){
						__builtin_prefetch (&b[j+rd_dist], 0, 0);
						__builtin_prefetch (&c[j+rd_dist], 0, 0);
						__builtin_prefetch (&d[j+rd_dist], 0, 0);
						__builtin_prefetch (&e[j+rd_dist], 0, 0);
						for(inner=0;inner<4;inner++)
						a[j+inner] = c[j+inner] + d[j+inner] + e[j+inner] + b[j+inner];
				}
				for(j = ((size/sizeof(double))-rd_dist);j <(size/sizeof(double));j++){
						a[j] = c[j] + d[j] + e[j] + b[j];
				}
			}else{
				for(j = ((size/sizeof(double))-rd_dist);j <((size/sizeof(double))-wr_dist);j+=4){
						__builtin_prefetch (&a[j+wr_dist], 1, 0);
						for(inner=0;inner<4;inner++)
						a[j+inner] = c[j+inner] + d[j+inner] + e[j+inner] + b[j+inner];
				}
				for(j = ((size/sizeof(double))-wr_dist);j <(size/sizeof(double));j++){
						a[j] = c[j] + d[j] + e[j] + b[j];
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			owfr_t += accum;
			if(accum < owfr_tmin)
					owfr_tmin = accum;
			if(accum <= empty){
				goto redo5;
			}
			owfr_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
/*redo6:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = b[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo6;
			}
			twor_avg += ((3*size*1.0E-06)/(long double)(accum - empty));
			}
redo7:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] + b[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo7;
			}
			twtr_avg += ((4*size*1.0E-06)/(long double)(accum - empty));
			}
redo8:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] + b[j]+ e[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo8;
			}
			twthr_avg += ((5*size*1.0E-06)/(long double)(accum - empty));
			}
redo9:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] + b[j]+ e[j] + f[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo9;
			}
			twfr_avg += ((6*size*1.0E-06)/(long double)(accum - empty));
			}
redo10:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = e[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo10;
			}
			thwor_avg += ((4*size*1.0E-06)/(long double)(accum - empty));
			}
redo11:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j]+ e[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo11;
			}
			thwtr_avg += ((5*size*1.0E-06)/(long double)(accum - empty));
			}
redo12:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j] + e[j] +f[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo12;
			}
			thwthr_avg += ((6*size*1.0E-06)/(long double)(accum - empty));
			}
redo13:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j] + e[j] + f[j] +g[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo13;
			}
			thwfr_avg += ((7*size*1.0E-06)/(long double)(accum - empty));
			}
redo14:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j] = e[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo14;
			}
			fwor_avg += ((5*size*1.0E-06)/(long double)(accum - empty));
			}
redo15:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j] = e[j] + f[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo15;
			}
			fwtr_avg += ((6*size*1.0E-06)/(long double)(accum - empty));
			}
redo16:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j] = e[j] + f[j] +g[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo16;
			}
			fwthr_avg += ((7*size*1.0E-06)/(long double)(accum - empty));
			}
redo17:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] = d[j] = b[j] = e[j] + f[j] +g[j] + h[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo17;
			}
			fwfr_avg += ((8*size*1.0E-06)/(long double)(accum - empty));
			}*/
redo18:
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
					if(wr_dist >= rd_dist){
                        for(j =0; j < (wr_dist - rd_dist); j++){
							__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
							stride +=3;
						}
						for(j =(wr_dist - rd_dist); j < wr_dist; j++){
							__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
							__builtin_prefetch (&b[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							__builtin_prefetch (&c[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							__builtin_prefetch (&d[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							__builtin_prefetch (&e[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							stride +=3;
						}
					}else{
                        for(j =0; j < (rd_dist - wr_dist); j++){
							__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
							stride +=3;
						}
						for(j =(rd_dist - wr_dist); j < rd_dist; j++){
							__builtin_prefetch (&a[stride%(size/sizeof(double))-(rd_dist - wr_dist)], 1, 0);
							__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
							stride +=3;
						}
					}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			if(rank == 0){
			accum = ( ostop.tv_sec - obegin.tv_sec ) + (long double)( ostop.tv_nsec - obegin.tv_nsec ) / (long double)BILLION;
			str_o += accum;
			if(accum < str_omin)
					str_omin = accum;
			}
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
                        for(j =0; j < ((size/sizeof(double)) - bdist); j+=8){
					__builtin_prefetch (&a[stride%(size/sizeof(double))+wr_dist], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double))+rd_dist], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))+rd_dist], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))+rd_dist], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))+rd_dist], 0, 0);
						for(inner=0;inner<8;inner++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    stride +=3;
						}
                        }
						if(wr_dist >= rd_dist){
                        	for(j = ((size/sizeof(double)) - wr_dist); j < ((size/sizeof(double)) -rd_dist); j+=8){
								__builtin_prefetch (&b[stride%(size/sizeof(double))+rd_dist], 0, 0);
								__builtin_prefetch (&c[stride%(size/sizeof(double))+rd_dist], 0, 0);
								__builtin_prefetch (&d[stride%(size/sizeof(double))+rd_dist], 0, 0);
								__builtin_prefetch (&e[stride%(size/sizeof(double))+rd_dist], 0, 0);
							for(inner=0;inner<8;inner++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
							}
							for(j = ((size/sizeof(double)) - rd_dist); j < (size/sizeof(double)); j++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
						}else{
                        	for(j = ((size/sizeof(double)) - rd_dist); j < ((size/sizeof(double)) -wr_dist); j+=8){
								__builtin_prefetch (&a[stride%(size/sizeof(double))+wr_dist], 1, 0);
						for(inner=0;inner<8;inner++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
							}
							for(j = ((size/sizeof(double)) - rd_dist); j < (size/sizeof(double)); j++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			str_t += accum;
			if(accum < str_tmin)
					str_tmin = accum;
			if(accum <= empty){
				goto redo18;
			}
			str_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
redo19:
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
                        for(j =0; j < rd_dist; j++){
					__builtin_prefetch (&rand_tab[j], 0, 0);
			}
                       /* for(j =0; j < rd_dist; j++){
					__builtin_prefetch (&rand_tab[j+rd_dist], 0, 0);
					__builtin_prefetch (&a[rand_tab[j]], 1, 0);
					__builtin_prefetch (&b[rand_tab[j]], 0, 0);
					__builtin_prefetch (&c[rand_tab[j]], 0, 0);
					__builtin_prefetch (&d[rand_tab[j]], 0, 0);
					__builtin_prefetch (&e[rand_tab[j]], 0, 0);
						}
			*/
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			if(rank == 0){
			accum = ( ostop.tv_sec - obegin.tv_sec ) + (long double)( ostop.tv_nsec - obegin.tv_nsec ) / (long double)BILLION;
			rand_o += accum;
			if(accum < rand_omin)
					rand_omin = accum;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
                        for(j =0; j < ((size/sizeof(double)) - (rd_dist)); j+=8){
					__builtin_prefetch (&rand_tab[j+rd_dist], 0, 0);
				//	__builtin_prefetch (&a[rand_tab[j+rd_dist]], 1, 0);
				//	__builtin_prefetch (&b[rand_tab[j+rd_dist]], 0, 0);
				//	__builtin_prefetch (&c[rand_tab[j+rd_dist]], 0, 0);
				//	__builtin_prefetch (&d[rand_tab[j+rd_dist]], 0, 0);
				//	__builtin_prefetch (&e[rand_tab[j+rd_dist]], 0, 0);
						for(inner=0;inner<8;inner++)
			  a[rand_tab[j+inner]] = b[rand_tab[j+inner]] + c[rand_tab[j+inner]] + d[rand_tab[j+inner]] + e[rand_tab[j+inner]];
                        }
                        /*for(j = ((size/sizeof(double)) - (2*rd_dist));j < ((size/sizeof(double)) - rd_dist); j+=32){
					__builtin_prefetch (&a[rand_tab[j+rd_dist]], 1, 0);
					__builtin_prefetch (&b[rand_tab[j+rd_dist]], 0, 0);
					__builtin_prefetch (&c[rand_tab[j+rd_dist]], 0, 0);
					__builtin_prefetch (&d[rand_tab[j+rd_dist]], 0, 0);
					__builtin_prefetch (&e[rand_tab[j+rd_dist]], 0, 0);
						for(inner=0;inner<32;inner++)
			  a[rand_tab[j+inner]] = b[rand_tab[j+inner]] + c[rand_tab[j+inner]] + d[rand_tab[j+inner]] + e[rand_tab[j+inner]];
                        }*/
                        for(j = ((size/sizeof(double)) - rd_dist);j < (size/sizeof(double)); j++){
			  a[rand_tab[j]] = b[rand_tab[j]] + c[rand_tab[j]] + d[rand_tab[j]] + e[rand_tab[j]];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			rand_t += accum;
			if(accum < rand_tmin)
					rand_tmin = accum;
			if(accum <= empty){
				goto redo19;
			}
			rand_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
/*redo20:
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = e[j] + f[j] +g[j] + h[j];
			    c[j] = e[j] + f[j] +g[j] + h[j];
			    d[j] = e[j] + f[j] +g[j] + h[j];
			    b[j] = e[j] + f[j] +g[j] + h[j];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty){
				goto redo20;
			}
			diff_avg += ((8*size*1.0E-06)/(long double)(accum - empty));
			}
redo21:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
            for(j = 1; j < (r_size/sizeof(double*)) - 1; j++){
				for(k = 1; k < dist; k++){
					__builtin_prefetch (&aa[j][k], 1, 0);
					__builtin_prefetch (&bb[j][k], 0, 0);
					__builtin_prefetch (&cc[j][k], 0, 0);
				}
				for(k = 1; k < (c_size/sizeof(double)) - dist - 1; k++){
					__builtin_prefetch (&aa[j][k+wr_dist], 1, 0);
					__builtin_prefetch (&bb[j][k+rd_dist], 0, 0);
					__builtin_prefetch (&cc[j][k+rd_dist], 0, 0);
                            		aa[j][k] = bb[j][k]*cc[j][k];
				}
				for(k = (c_size/sizeof(double)) - dist - 1; k < (c_size/sizeof(double)) - 1; k++){
						aa[j][k] = bb[j][k]*cc[j][k];
				}
            }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty2){
				goto redo21;
			}
			row_avg += ((long double)(3*(long)(r_size-2)*(c_size-2)*1.0E-06)/(long double)(accum - empty2));
			}	
			
redo22:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
            for(j = 1; j < (c_size/sizeof(double)) - 1; j++){
				for(k = 1; k < dist; k++){
					__builtin_prefetch (&aa[k][j], 1, 0);
					__builtin_prefetch (&bb[k][j], 0, 0);
					__builtin_prefetch (&cc[k][j], 0, 0);
				}
				for(k = 1; k < ((r_size/sizeof(double*)) - dist - 1); k++){
					__builtin_prefetch (&aa[k+wr_dist][j], 1, 0);
					__builtin_prefetch (&bb[k+rd_dist][j], 0, 0);
					__builtin_prefetch (&cc[k+rd_dist][j], 0, 0);
                    aa[k][j] = bb[k][j]*cc[k][j];
                }
				for(k = (r_size/sizeof(double*)) - dist - 1; k < (r_size/sizeof(double*)) - 1; k++){
						aa[k][j] = bb[k][j]*cc[k][j];
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty2){
				goto redo22;
			}
			col_avg += ((long double)(3*(long)(r_size-2)*(c_size-2)*1.0E-06)/(long double)(accum - empty2));
			}
redo23:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
            for(j = 1; j < (r_size/sizeof(double*)) - 1; j++){
				for(k = 1; k < dist; k++){
					__builtin_prefetch (&aa[j][k], 1, 0);
					__builtin_prefetch (&bb[j][k], 0, 0);
					__builtin_prefetch (&cc[k][j], 0, 0);
				}
				for(k = 1; k < ((c_size/sizeof(double)) - dist - 1); k++){
					__builtin_prefetch (&aa[j][k+wr_dist], 1, 0);
					__builtin_prefetch (&bb[j][k+rd_dist], 0, 0);
					__builtin_prefetch (&cc[k+rd_dist][j], 0, 0);
                            		aa[j][k] = bb[j][k]*cc[k][j];
                        }
		printf("HEre: %d %d %d\n", rank, k, dist);
		fflush(NULL);
				for(k = ((c_size/sizeof(double)) - dist - 1); k < (c_size/sizeof(double)) - 1; k++){
						aa[j][k] = bb[j][k]*cc[k][j];
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty2){
				goto redo23;
			}
			rc_avg += ((long double)(3*(long)(r_size-2)*(c_size-2)*1.0E-06)/(long double)(accum - empty2));
			}
redo24:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (r_size/sizeof(double*)); j++){
				for(k = 0; k < (c_size/sizeof(double)); k++)
					if((k!=0)&&(k!=((c_size/sizeof(double))-1))&&(j!=0)&&(j!=((r_size/sizeof(double*))-1)))
                            			aa[j][k] = aa[j][k-1]+aa[j][k+1] + aa[j-1][k] + aa[j+1][k];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty2){
				goto redo24;
			}
			f_sten_avg += 5*((long double)(((long)r_size-1)*((long)c_size-1)*1.0E-06)/(long double)(accum - empty2));
			}
redo25:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (r_size/sizeof(double*)); j++){
                                for(k = 0; k < (c_size/sizeof(double)); k++)
                                        if((k!=0)&&(k!=((c_size/sizeof(double))-1)))
                                                aa[j][k] = aa[j][k-1]+aa[j][k+1];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty2){
				goto redo25;
			}
			t_sten_avg += 3*((long double)(((long)r_size)*(c_size)*1.0E-06)/(long double)(accum - empty2));
			}
redo26:
			rs = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			if(rank >= 12){
				MPI_Isend(&aa[1], (c_size/sizeof(double)), MPI_DOUBLE, (rank - 12), rs, MPI_COMM_WORLD, &reqs[rs]);
				rs++;
			}
			if(rank < (procs - 12)){
				MPI_Irecv(&aa[(r_size/sizeof(double*))-1], (c_size/sizeof(double)), MPI_DOUBLE, (rank + 12), MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[rs]);
				rs++;
			}
			if(rank < (procs - 12)){
				MPI_Isend(&aa[(r_size/sizeof(double*))-2], (c_size/sizeof(double)), MPI_DOUBLE, (rank + 12), rs, MPI_COMM_WORLD, &reqs[rs]);
				rs++;
			}
			if(rank >= 12){
				MPI_Irecv(&aa[0], (c_size/sizeof(double)), MPI_DOUBLE, (rank - 12), MPI_ANY_TAG , MPI_COMM_WORLD, &reqs[rs]);
				rs++;
			}
			if((rank%12 < 11)){
					for(z=1; z<(r_size/sizeof(double*))-2; z++) {
						MPI_Isend(&aa[z][(c_size/sizeof(double)) - 2], 1, MPI_DOUBLE, (rank + 1), rs, MPI_COMM_WORLD, &reqs[rs]);
						rs++;
					}
			}
			if((rank%12 != 0)){
					for(z=1; z<(r_size/sizeof(double*))-2; z++) {
						MPI_Irecv(&aa[z][0], 1, MPI_DOUBLE, (rank - 1), MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[rs]);
						rs++;
					}
			}
			if((rank%12 != 0)){
					for(z=1; z<(r_size/sizeof(double*))-2; z++) {
						MPI_Isend(&aa[z][1], 1, MPI_DOUBLE, (rank - 1), rs, MPI_COMM_WORLD, &reqs[rs]);
						rs++;
					}
			}
			if((rank%12 < 11)){
					for(z=1; z<(r_size/sizeof(double*))-2; z++) {
						MPI_Irecv(&aa[z][(c_size/sizeof(double)) - 1], 1, MPI_DOUBLE, (rank + 1), MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[rs]);
						rs++;
					}
			}
			MPI_Waitall (rs , reqs , stat );
			//printf("Here %d\n", rank);
			//fflush(NULL);

//#pragma omp parallel for
                        for(j =1; j < (r_size/sizeof(double*)) -1; j++){
					__builtin_prefetch (&aa[j-1][k-1], 0, 2);
					__builtin_prefetch (&aa[j-1][k], 0, 2);
					__builtin_prefetch (&aa[j][k-1], 1, 2);
					__builtin_prefetch (&aa[j][k], 1, 2);
					__builtin_prefetch (&aa[j+1][k-1], 1, 2);
					__builtin_prefetch (&aa[j+1][k], 1, 2);
                                for(k = 1; k < dist; k++){
					__builtin_prefetch (&aa[j-1][k+1], 0, 2);
					__builtin_prefetch (&aa[j][k+1], 1, 2);
					__builtin_prefetch (&aa[j+1][k+1], 1, 2);
								}
                                for(k = 1; k < (c_size/sizeof(double))- dist - 1; k++){
                                      //  if((k!=0)&&(k!=((c_size/sizeof(double))-1))&&(j!=0)&&(j!=((r_size/sizeof(double*))-1)))
					__builtin_prefetch (&aa[j-1][k+1+rd_dist], 0, 2);
					__builtin_prefetch (&aa[j][k+1+wr_dist], 1, 2);
					__builtin_prefetch (&aa[j+1][k+1+wr_dist], 1, 2);
                                                aa[j][k] = aa[j][k-1]+aa[j][k+1] + aa[j-1][k] + aa[j+1][k] + aa[j-1][k-1] + aa[j-1][k+1] + aa[j+1][k-1] + aa[j+1][k+1];
                        }
                                for(k = ((c_size/sizeof(double))- dist - 1); k < (c_size/sizeof(double))- 1; k++){
                                                aa[j][k] = aa[j][k-1]+aa[j][k+1] + aa[j-1][k] + aa[j+1][k] + aa[j-1][k-1] + aa[j-1][k+1] + aa[j+1][k-1] + aa[j+1][k+1];
								}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			if(accum <= empty2){
				goto redo26;
			}
			n_sten_avg += 9*((long double)(((long)r_size)*((long)c_size)*1.0E-06)/(long double)(accum - empty2));
			}*/
redo27:
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
                        	for(j =0; j < wr_dist - rd_dist; j++){
					__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
                        	for(j =(wr_dist - rd_dist); j < wr_dist; j++){
					__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}else{
                        	for(j =0; j < rd_dist - wr_dist; j++){
					__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
                        	for(j =(rd_dist - wr_dist); j < rd_dist; j++){
					__builtin_prefetch (&a[stride%(size/sizeof(double))-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			if(rank == 0){
			accum = ( ostop.tv_sec - obegin.tv_sec ) + (long double)( ostop.tv_nsec - obegin.tv_nsec ) / (long double)BILLION;
			l2cache_o += accum;
			if(accum < l2cache_omin)
					l2cache_omin = accum;
			}
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
                        for(j =0; j < (size/sizeof(double)) - bdist; j+=32){
					__builtin_prefetch (&a[stride%(size/sizeof(double)) + wr_dist], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double)) + rd_dist], 0, 0);
						for(inner=0;inner<32;inner++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if(((j+inner)%8 == 0)&&((j+inner) != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
                        }
				}
						if(wr_dist >= rd_dist){
                        for(j = ((size/sizeof(double)) - wr_dist); j < (size/sizeof(double)) - rd_dist; j+=32){
					__builtin_prefetch (&b[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double)) + rd_dist], 0, 0);
						for(inner=0;inner<32;inner++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if(((j+inner)%8 == 0)&&((j+inner) != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
						 for(j = (size/sizeof(double))-rd_dist; j < (size/sizeof(double)); j++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}else{
                        for(j = ((size/sizeof(double)) - rd_dist); j < (size/sizeof(double)) - wr_dist; j+=32){
					__builtin_prefetch (&a[stride%(size/sizeof(double)) + wr_dist], 1, 0);
						for(inner=0;inner<32;inner++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if(((j+inner)%8 == 0)&&((j+inner) != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
						 for(j = (size/sizeof(double))-wr_dist; j < (size/sizeof(double)); j++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			l2cache_t += accum;
			if(accum < l2cache_tmin)
					l2cache_tmin = accum;
                        if(accum <= empty){
                                goto redo27;
                        }
                        l2cache_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
			numa_free(a, size);
			numa_free(b, size);
			numa_free(c, size);
			numa_free(d, size);
			numa_free(e, size);
			//numa_free(f, size);
			//numa_free(g, size);
			//numa_free(h, size);
			/*for(j = 0; j < r_size/sizeof(double*); j++){
				numa_free(aa[j], c_size);
			}
			for(j = 0; j < r_size/sizeof(double*); j++){
				numa_free(bb[j], c_size);
			}
			for(j = 0; j < r_size/sizeof(double*); j++){
				numa_free(cc[j], c_size);
			}
			numa_free(aa, r_size);
			numa_free(bb, r_size);
			numa_free(cc, r_size);*/
		}
		if(rank == 0){
		struct numa_node_bw * node_bw = (struct numa_node_bw *)malloc(sizeof(struct numa_node_bw));
		node_bw->numa_id = numa_node_ids[i];
		node_bw->wr_only_avg = wr_only_avg/10;
		node_bw->owor_avg = owor_avg/10;
		node_bw->owtr_avg = owtr_avg/10;
		node_bw->owthr_avg = owthr_avg/10;
		node_bw->owfr_avg = owfr_avg/10;
		node_bw->twor_avg = twor_avg/10;
		node_bw->twtr_avg = twtr_avg/10;
		node_bw->twthr_avg = twthr_avg/10;
		node_bw->twfr_avg = twfr_avg/10;
		node_bw->thwor_avg = thwor_avg/10;
		node_bw->thwtr_avg = thwtr_avg/10;
		node_bw->thwthr_avg = thwthr_avg/10;
		node_bw->thwfr_avg = thwfr_avg/10;
		node_bw->fwor_avg = fwor_avg/10;
		node_bw->fwtr_avg = fwtr_avg/10;
		node_bw->fwthr_avg = fwthr_avg/10;
		node_bw->fwfr_avg = fwfr_avg/10;
		node_bw->str_avg = str_avg/10;
		node_bw->rand_avg = rand_avg/10;
		node_bw->diff_avg = diff_avg/10;
		node_bw->row_avg = row_avg/10;
		node_bw->col_avg = col_avg/10;
		node_bw->rc_avg = rc_avg/10;
		node_bw->l2cache_avg = l2cache_avg/10;
		node_bw->t_sten_avg = t_sten_avg/10;
		node_bw->f_sten_avg = f_sten_avg/10;
		node_bw->n_sten_avg = n_sten_avg/10;
		node_bw->wr_only_t = wr_only_t/10;
		node_bw->l2cache_t = l2cache_t/10;
		node_bw->rand_t = rand_t/10;
		node_bw->str_t = str_t/10;
		node_bw->owfr_t = owfr_t/10;
		node_bw->wr_only_tmin = wr_only_tmin;
		node_bw->l2cache_tmin = l2cache_tmin;
		node_bw->rand_tmin = rand_tmin;
		node_bw->str_tmin = str_tmin;
		node_bw->owfr_tmin = owfr_tmin;
		node_bw->wr_only_min = ((5*size*procs*1.0E-06)/(long double)wr_only_tmin);
		node_bw->l2cache_min = ((5*size*procs*1.0E-06)/(long double)l2cache_tmin);
		node_bw->rand_min = ((5*size*procs*1.0E-06)/(long double)rand_tmin);
		node_bw->str_min = ((5*size*procs*1.0E-06)/(long double)str_tmin);
		node_bw->owfr_min = ((5*size*procs*1.0E-06)/(long double)owfr_tmin);
		node_bw->wr_only_o = wr_only_o/10;
		node_bw->l2cache_o = l2cache_o/10;
		node_bw->rand_o = rand_o/10;
		node_bw->str_o = str_o/10;
		node_bw->owfr_o = owfr_o/10;
		node_bw->wr_only_omin = wr_only_omin;
		node_bw->l2cache_omin = l2cache_omin;
		node_bw->rand_omin = rand_omin;
		node_bw->str_omin = str_omin;
		node_bw->owfr_omin = owfr_omin;
		node_bw->next = NULL;
		if(numa_node_list == NULL){
			numa_node_list = node_bw;
			numa_list_head = numa_node_list;
		}
		else{
			sort_list(node_bw);
		}
		}
		i++;
	}
	if(rank == 0){
	classify();
	write_config_file();
	}
	free(rand_tab);
}

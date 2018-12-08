#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#define set_bit(buffer,value) (*(buffer + (value / BITS_PER_LL)) |= (ONE << (value % BITS_PER_LL)))
#define  ull unsigned long long
#define  WORDS_PER_LL  sizeof(ull)
#define  BITS_PER_LL  (WORDS_PER_LL * 8)
ull ONE = 1;

/**
 * merges removing duplicates
 * returns number of duplicates that were removed while merging
 */
void merge_bit_lists(ull *a, ull *b, ull a_offset, ull b_offset, ull n_longs) {

	ull *ao = a + a_offset;
	ull *bo = b + b_offset;

	for (ull i = 0; i < n_longs; i++) {
		*(ao + i) |= *(bo + i);
	}

}

/**
 * allocates memory for a bit list of size n_longs
 */
ull *allocate_and_clear(ull n_longs) {
	//ull *array = malloc(sizeof(ull) * n_longs);
	//memset(array, 0, sizeof(ull) * n_longs);
	return calloc(WORDS_PER_LL, n_longs);
}

// macro is faster due to no stack usage, no jumps,etc, will be inplace instead
// for n=100k , saves 6 seconds on 4 cores
/*void set_bit( ull * buffer, ull value) {

 ull l = value / BITS_PER_LL;
 ull b = value % BITS_PER_LL;

 ull mask = 1;

 *(buffer + l) |= (mask << b);

 }*/

/**
 * generates the bit list for this processor with the specified range
 */
void get_bit_list(ull register start, ull register end, ull register *buffer) {

	ull register col, row;

	col = 1;
	row = 1;

	// go to starting position of this process' range

	ull ind = start + 1;
	col = floor(sqrt(ind * 2) + 0.5);
	row = col * col - col;
	row /= 2;
	row = ind - row;

	for (ull register index = start; index < end; index++) {

		set_bit(buffer, row * col);
		row++;
		if (row > col) {
			row = 1;
			col++;

		}

	}

}
/**
 * counts the number of high bits in the buffer
 */
ull count_bits(ull *buffer, ull offset, ull n_longs) {

	ull count = 0;
	ull *bo = buffer + offset;

	for (ull i = 0; i < n_longs; i++) {

		ull ll = *(bo + i);

		while (ll > 0) {
			count += (ll & 1);
			ll >>= 1;

		}

	}
	return count;

}

/**
 * gets time elapsed from a start and end time
 */
double get_seconds(struct timeval start_time, struct timeval end_time) {

	double us = end_time.tv_usec - start_time.tv_usec;
	double s = end_time.tv_sec - start_time.tv_sec;

	return s + (us / 1.0e6);

}

/**
 * nicely prints time elapsed to the screen from number of seconds
 */
void print_elapsed_sec(double took, char * job) {

	if (took < 1) {
		took *= 1000;

		if (took < 1) {

			took *= 1000;
			printf("%s took %f microseconds\n", job, took);
		} else {

			printf("%s took %f milliseconds\n", job, took);
		}
	} else if (took > 60) {
		took /= 60;

		if (took > 60) {
			took /= 60;
			printf("%s took %f hours\n", job, took);
		} else {
			printf("%s took %f minutes\n", job, took);
		}
	} else {
		printf("%s took %f seconds\n", job, took);
	}

}

/**
 * nicely prints time elapsed to the screen from a start and end time
 */
void print_elapsed_range(struct timeval start_time, struct timeval end_time, char * job) {

	double took = get_seconds(start_time, end_time);

	print_elapsed_sec(took, job);

}

int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	MPI_Status status;
	struct timeval stop_time, start_time, merge_start_time, generate_end_time;
	ull start, end, partition_size, num_upper_tri;
	ull *a, *b;
	int id, p;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//MPI_Barrier(MPI_COMM_WORLD);

	gettimeofday(&start_time, NULL);

	ull n = atol(argv[1]);

	num_upper_tri = (n * n - n) / 2 + n;

	// for each process, calculate their start, range
	partition_size = num_upper_tri / p;
	start = id * partition_size;

	if (id == (p - 1)) {
		end = num_upper_tri;
	} else {
		end = start + partition_size;
	}

	ull n_bits = n * n;
	n_bits += 1;

	// calculate how many long longs are needed to hold a bit for each unique element in the nxn table

	ull n_longs = n_bits / BITS_PER_LL;

	if ((n_bits % BITS_PER_LL) > 0) {
		n_longs++;
	}

	// allocate memory for this process' bit list
	a = allocate_and_clear(n_longs);
	b = NULL;

	// generate list of all unique elements in this process' range
	get_bit_list(start, end, a);

	gettimeofday(&generate_end_time, NULL);

	//MPI_Barrier(MPI_COMM_WORLD);

	// now we scatter results, proc 0 gets slice 0 from all, proc i gets slice i from all

	gettimeofday(&merge_start_time, NULL);

	ull longs_per_core = n_longs / p;
	ull extra_longs = n_longs % p;
	ull offset = 0;
	ull offset_for_this_core = longs_per_core * id;

	ull longs_for_this_core = longs_per_core;
	if (id < extra_longs) {
		longs_for_this_core++;
		offset_for_this_core += id;
	} else {
		offset_for_this_core += extra_longs;
	}

	//allocate bit list to store recieved data before merge
	b = allocate_and_clear(longs_for_this_core * p);

	MPI_Request r_jobs[p];
	MPI_Request s_jobs[p];

	// send and recieve slices of bit list and merge them with slice of this process' bit list

	for (int i = 0; i < p; i++) {
		ull longs_for_core_i = longs_per_core;
		if (i < extra_longs) {
			longs_for_core_i++;
		}
		if (i != id) {	// for all but this core
			MPI_Isend(a + offset, longs_for_core_i, MPI_UNSIGNED_LONG_LONG, i, 1, MPI_COMM_WORLD, s_jobs + i);

		} else {

			for (int j = 0; j < p; j++) {

				if (j != id) {

					MPI_Irecv(b + (longs_for_this_core * j), longs_for_this_core, MPI_UNSIGNED_LONG_LONG, j, 1,
							MPI_COMM_WORLD, r_jobs + j);

					//merge_bit_lists(a, b, offset_for_this_core,longs_for_this_core);
				}
			}
		}

		offset += longs_for_core_i;
	}

	for (int i = 0; i < p; i++) {
		if (i != id) {

			MPI_Wait(r_jobs + i, &status);
			MPI_Wait(s_jobs + i, &status);

		}

	}

	struct timeval count_start_time;
	gettimeofday(&count_start_time, NULL);

	for (int i = 0; i < p; i++) {
		if (i != id) {
			merge_bit_lists(a, b, offset_for_this_core, longs_for_this_core * i, longs_for_this_core);
		}
	}

	free(b);

	// count high bits in bit list in the slice this process is responsible for

	ull count = count_bits(a, offset_for_this_core, longs_for_this_core);

	free(a);

	// send all counts to process 0 to get master count

	if (id == 0) { // receive counts
		ull i_count = 0;
		for (int i = 1; i < p; i++) {
			MPI_Recv(&i_count, 1, MPI_UNSIGNED_LONG_LONG, i, 2, MPI_COMM_WORLD, &status);
			count += i_count;
		}
	} else {
		MPI_Send(&count, 1, MPI_UNSIGNED_LONG_LONG, 0, 2, MPI_COMM_WORLD);
	}

	//old merging code , 1/2 cores working 1/2 of those and then 1/2 of those and so on  7 or passes for 256 cores
	// merge bit lists from all processes
	// on first pass , half the processes will send their lists to the other half to be merged
	// on subsequent passes half the processes that recieved will send to the other half and so on
	/*int height, process, proc_recv;

	 height = log2(p);
	 for (int i = 0; i < height; i++) {
	 MPI_Barrier(MPI_COMM_WORLD);

	 process = p / pow(2, i + 1);
	 for (int j = 0; j < process; j++) {
	 if (id == j) {

	 proc_recv = process + j;

	 if (b == NULL) {
	 b = allocate_and_clear(n_longs);
	 }
	 MPI_Recv(b, n_longs, MPI_UNSIGNED_LONG_LONG, proc_recv, 3, MPI_COMM_WORLD, &status);

	 merge_bit_lists(a, b, n_longs);

	 } else if (id == (process + j)) {

	 MPI_Send(a, n_longs, MPI_UNSIGNED_LONG_LONG, j, 3, MPI_COMM_WORLD);
	 free(a);
	 if (b != NULL) {
	 free(b);
	 }
	 }
	 }
	 }*/

	//MPI_Barrier(MPI_COMM_WORLD);
	double generating_time = get_seconds(start_time, generate_end_time);

	// process 0 now has all count info
	// so now it's time to print output
	if (id == 0) {

		printf("Using %d cores\n", p);

		// print answer
		printf("M(%lli) = %lli\n", n, count);

		// print timing info
		gettimeofday(&stop_time, NULL);
		print_elapsed_range(start_time, merge_start_time, "Generating");

		// print generation time for each process
		print_elapsed_sec(generating_time, "  Generating for Process 0");
		for (int proc_recv = 1; proc_recv < p; proc_recv++) {
			char c[40];
			sprintf(c, "  Generating for Process %d", proc_recv);
			MPI_Recv(&generating_time, 1, MPI_DOUBLE, proc_recv, 4, MPI_COMM_WORLD, &status);
			print_elapsed_sec(generating_time, c);
		}

		print_elapsed_range(merge_start_time, count_start_time, "Sending slices");
		print_elapsed_range(count_start_time, stop_time, "Merging & Counting");

		print_elapsed_range(start_time, stop_time, "The whole job");

	} else {
		// send timing info to process 0
		MPI_Send(&generating_time, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);

	}

	MPI_Finalize();
	return 0;
}

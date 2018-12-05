#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

int ID = -1;

int BITS_PER_LL = sizeof(long long) * 8;


/**
* merges removing duplicates
* returns number of duplicates that were removed while merging
*/
void merge_bit_lists(long long *a, long long *b,long long n_longs ) {
    
	for (long long i = 0; i < n_longs; i++){
	
		* (a+i) |= *(b+i);
		
	}
	

}


void set_bit( long long * buffer, long long value){
	
	long long l = value / BITS_PER_LL;
	long long b = value % BITS_PER_LL;
	
	long long mask = 1;
	
	mask <<= b;
	
	* ( buffer+l) |= mask;
	
	
}


/*
i = floor ( sqrt(index*2)+1/2)
j = index - (i*i-i)/2
*/
long long get_bit_list(long long start, long long end, long long* buffer){
    
    
  
    long long size = end - start;
    long long newsize = 0;
    long long col,row;
    long long *temp1 ;
    long long *temp2 ;
    
	
    long long count = 0;

    
    long long ind=start+1;
    col = floor(sqrt(ind * 2) + 0.5);
    row = ind - (col * col - col) / 2;
	
	
	fprintf(stderr,"ID:%d, col:%lli , row%lli \n", ID,col, row);
    
  
    
    for (long long index = start; index < end; index++){
        long long value = row* col;
		count ++;
		
					
		
		
		set_bit(buffer, row * col);
        row++;
        if (row > col){
            row = 1;
            col ++;
            
           
            
            
        }
        
    }
	
	fprintf(stderr,"ID:%d, col:%lli , row%lli \n", ID,col, row);
    
	
	return count;
    
}


long long count_bits (long long *buffer, long long n_longs){
	
	long long count = 0;
	
	for (long long i = 0; i < n_longs; i++){
		long long l = *(buffer+i);
		for (long long i2 = 0; i2 < BITS_PER_LL; i2 ++){
			count += (l & 1);
			l >>=1;
			
		}
		
	}
	return count;
	
	
	
}



int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Status status;

    long long start, end, partition_size, num_upper_tri;
    long long *a;
    int id, p;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    struct timeval stop_time, start_time;
    gettimeofday(&start_time, NULL);

  
    
    
    ID = id;
   
    long long n = atol(argv[1]);

    num_upper_tri = (n*n-n)/2+n;

    // for each process, calculate their start, range, sort and remove duplicates
    partition_size = num_upper_tri/p;
    start = id*partition_size;
    
   

    
    if(id == (p-1)) {
        end = num_upper_tri;
    } else {
        end = start + partition_size;
    }
    
    long long newsize = end-start;
	
	
	
	long long n_bits = n*n+1;
	long long n_longs = n_bits / BITS_PER_LL;
	
	if ((n_bits  %BITS_PER_LL ) >0){
		n_longs++;
	}
	fprintf( stderr,"ID:%d, generating size:%lli (%lli bits)\n", ID, n_longs, n_bits);
	
	a = malloc ( sizeof(long long) * n_longs);
	
    memset(a, 0, sizeof(long long)*n_longs );
	
	
	long long *b;

    b = malloc ( sizeof(long long) * n_longs);
	
    memset(b, 0, sizeof(long long)*n_longs );
    
	
	fflush(stderr);
	
	
	 
	
    long long count = get_bit_list(start, end, a);

   
	//fprintf( stderr,"ID:%d, generated, count:%lli, start:%lli, end%lli\n", ID,count_bits(a,n_longs), i_start, i_end);
	fprintf( stderr,"ID:%d, generated, count:%lli, start:%lli, end:%lli\n", ID,count, start, end);
    fflush(stderr);
    
   // printf("Process %d has array, start:%d, end:%d, size:%d, newsize:%d, n:%d \n", id, start, end,partition_size,newsize, n );
    
    
    
    int height, process, proc_recv;
	
	
    

    height = log2(p);

	
	
    for (int i = 0; i < height; i++){
         MPI_Barrier(MPI_COMM_WORLD);

         process = p / pow(2,i+1);
        for (int j = 0; j < process; j++){
            //printf("Process %d merges into process %d at height %d of height %d \n", process+j, j, i, height);
            if(id == j){
                
              
                
                proc_recv = process+j;
				
			     
                MPI_Recv(b, n_longs, MPI_LONG_LONG, proc_recv, 3, MPI_COMM_WORLD, &status);
               fprintf(stderr, "Process %d recieved %d longs, tag 2 \n", ID, n_longs);
                fflush(stderr);
                
                /*
                // merge the 2 arrays, sort, remove duplicates
                newa = malloc(sizeof(long) * size2);
                serial_merge(&a[0], &a2[0], &newa, newsize, mergesize);
                newsize=find_dups(size2,newa);
                */
                merge_bit_lists(a, b, n_longs);
                
				   
                
                /*a = malloc(sizeof(long) * size2);
                for(int k = 0; k < newsize; k++){
                    a[k] = newa[k];
                }*/

            }else if(id == (process+j)) {
                fprintf(stderr,"Process %d sends to process %d at height %d of height %d \n", ID, j, i, height);
                fflush(stderr);
                
				
				merge_bit_lists(a, b, n_longs);
                
                MPI_Send(a, n_longs, MPI_LONG_LONG, j, 3, MPI_COMM_WORLD);
            }
        }
    }

     MPI_Barrier(MPI_COMM_WORLD);
  
    
    if(id == 0){
		
		
		long long count = count_bits (a, n_longs);
        
        printf("M(%lli) = %lli\n", n, count);
        gettimeofday(&stop_time, NULL);
        long tl = stop_time.tv_usec - start_time.tv_usec;
        double tls = stop_time.tv_sec - start_time.tv_sec;
        double took = tl;
        
        if (tls >0){
			took /= 1e6;
            took += tls;
            
            
                    if (took > 60){
                        took /= 60;
                        
                        if (took >60){
                            took /= 60;
                            printf("took %f hours\n", took);
                        }else{
                             printf("took %f minutes\n", took);
                        }
                    }else{
                         printf("took %f seconds\n", took);
                    }
               
            
        }else{
            
        
        
            if (took > 1000){
                took /= 1000;
                if (took > 1000){
                    took /= 1000;
                    if (took > 60){
                        took /= 60;
                        
                        if (took >60){
                            took /= 60;
                            printf("took %f hours\n", took);
                        }else{
                             printf("took %f minutes\n", took);
                        }
                    }else{
                         printf("took %f seconds\n", took);
                    }
                }else{
                     printf("took %f milliseconds\n", took);
                }
            }else{
                 printf("took %f microseconds\n", took);
            }
        }
        
                   
               
            
           
        
        
    }
	
	free (a);
	free (b);

   
    
    
   
    MPI_Finalize();
    return 0;
}

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>


#define  ull unsigned long long

int ID = -1;

int BITS_PER_LL = sizeof(ull) * 8;

ull MX_BIT = 0;
ull MX_LONG = 0;



/**
* merges removing duplicates
* returns number of duplicates that were removed while merging
*/
void merge_bit_lists(ull *a, ull *b,ull n_longs ) {

    for (ull i = 0; i < n_longs; i++){

        * (a+i) |= *(b+i);

    }


}


void set_bit( ull * buffer, ull value){

  //  if ( value >= MX_BIT){
//	fprintf(stderr, "***value:%lli mx_bit:%lli \n", value, MX_BIT);
     //   fflush(stderr);
   // }

    ull l = value / BITS_PER_LL;
    ull b = value % BITS_PER_LL;

    if (l >= MX_LONG){
        fprintf(stderr, "***l:%lli max:%lli \n", l, MX_LONG);
        fflush(stderr);
    }
    //assert (l < MX_LONG);

    ull mask = 1;

    mask <<= b;

    * ( buffer+l) |= mask;


}


/*
i = floor ( sqrt(index*2)+1/2)
j = index - (i*i-i)/2
*/
ull get_bit_list(ull start, ull end, ull* buffer){



    ull size = end - start;
    ull newsize = 0;
    ull col,row;
    ull *temp1 ;
    ull *temp2 ;


    ull count = 0;


    ull ind=start+1;
    col = floor(sqrt(ind * 2) + 0.5);
    //row = ind - (col * col - col) / 2;
   row=col*col-col;
	fprintf(stderr, "**row1:%lli\n",row);
	row=row/2;
fprintf(stderr, "**row2:%lli\n",row);
	row=ind-row;
fprintf(stderr, "**row3:%lli\n",row);

if ( (row*col) >= MX_BIT){
                fprintf(stderr, "***rowxcol:%lli mx_bit:%lli row:%llu col:%llu start:%llu end:%llu \n", row*col, MX_BIT, row, col,start,end);
                fflush(stderr);

        }
    //fprintf(stderr,"ID:%d, starting at col:%lli , row:%lli \n", ID,col, row);



    for (ull index = start; index < end; index++){
        ull value = row* col;
        count ++;


    	if ( value >= MX_BIT){
        	fprintf(stderr, "***value:%lli mx_bit:%lli row:%llu col:%llu \n", value, MX_BIT, row, col);
        	fflush(stderr);

	}

        set_bit(buffer, row * col);
        row++;
        if (row > col){
            row = 1;
            col ++;




        }

    }

    //fprintf(stderr,"ID:%d, finished at 1 before col:%lli , row:%lli \n", ID,col, row);


    return count;

}


ull count_bits (ull *buffer, ull n_longs){

    ull count = 0;

    for (ull i = 0; i < n_longs; i++){
        ull ll = *(buffer+i);
        while (ll >0){
        //for (ull i2 = 0; i2 < BITS_PER_LL; i2 ++){
            count += (ll & 1);
            ll >>=1;

        }

    }
    return count;



}



int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Status status;

    ull start, end, partition_size, num_upper_tri;
    ull *a;
    int id, p;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);


    MPI_Barrier(MPI_COMM_WORLD);
    struct timeval stop_time, start_time;
    gettimeofday(&start_time, NULL);




    ID = id;

    ull n = atol(argv[1]);

    num_upper_tri = (n*n-n)/2+n;

    // for each process, calculate their start, range, sort and remove duplicates
    partition_size = num_upper_tri/p;
    start = id*partition_size;




    if(id == (p-1)) {
        end = num_upper_tri;
    } else {
        end = start + partition_size;
    }

    ull newsize = end-start;



    ull n_bits = n*n;
    n_bits+=1;
    ull n_longs = n_bits / BITS_PER_LL;


    MX_BIT = n_bits;


    if ((n_bits  %BITS_PER_LL ) >0){
        n_longs++;
    }

    MX_LONG = n_longs;

    fprintf( stderr,"ID:%d, generating size:%lli (%lli bits)\n", ID, n_longs, n_bits);

    a = malloc ( sizeof(ull) * n_longs);

    memset(a, 0, sizeof(ull)*n_longs );


    ull *b;

    b = malloc ( sizeof(ull) * n_longs);

    memset(b, 0, sizeof(ull)*n_longs );


    fflush(stderr);




    ull count = get_bit_list(start, end, a);


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


                MPI_Recv(b, n_longs, MPI_UNSIGNED_LONG_LONG, proc_recv, 3, MPI_COMM_WORLD, &status);
                fprintf(stderr, "Process %d recieved %d longs, tag 2 \n", ID, n_longs);
                fflush(stderr);


                merge_bit_lists(a, b, n_longs);





            }else if(id == (process+j)) {
                fprintf(stderr,"Process %d sends to process %d at height %d of height %d \n", ID, j, i, height);
                fflush(stderr);



                MPI_Send(a, n_longs, MPI_UNSIGNED_LONG_LONG, j, 3, MPI_COMM_WORLD);
            }
        }
    }

     MPI_Barrier(MPI_COMM_WORLD);


    if(id == 0){


        ull count = count_bits (a, n_longs);

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

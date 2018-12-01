#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>


int ID = -1;







/**
* merges removing duplicates
* returns number of duplicates that were removed while merging
*/
int serial_merge_no_dup(long *a, long *b, long *c, long size_a, long size_b) {
    long i = 0;
    long j = 0;
    long k = 0;
    long b_j;
    long a_i;
    long size_c = size_a + size_b;
    int duplicates = 0;

    
    while (k < size_c && i < size_a && j < size_b) {
        a_i = a[i];
        b_j = b[j];
        if (a_i <= b_j) {
            c[k] = a_i;
            i++;
        } else if (a_i > b_j){
            c[k] = b_j;
            j++;
        }
        
        
        if (k >0){
            if (c[k] == c[k-1]){
                k--;
                size_c--;
                duplicates++;
            }
        }
        k++;
        
        
    }
    if (i >= size_a) {
        for (; j < size_b; j++) {
            b_j = b[j];
            c[k] = b_j;
            if (k >0){
                if (c[k] == c[k-1]){
                    k--;
                    size_c--;
                    duplicates++;
                }
            }
            k++;
        }
    } else if (j >= size_b) {
        for (; i < size_a; i++) {
            a_i = a[i];
            c[k] = a_i;
            if (k >0){
                if (c[k] == c[k-1]){
                    k--;
                    size_c--;
                    duplicates++;
                }
            }
            k++;
        }
    }
    return duplicates;

}



/*
i = floor ( sqrt(index*2)+1/2)
j = index - (i*i-i)/2
*/
int get_array(long start, long end, long* buffer){
    
    
  
    long size = end - start;
    int newsize = 0;
    int col,row;
    long *temp1 ;
    long *temp2 ;
    
    temp1 = malloc(sizeof(long)*size);
    temp2 = malloc(sizeof(long)*size);
    
    
    long *temp = temp1;
    long *otemp = temp2;
    
    int prev_end = -1;
    
    int ind=start+1;
    col = floor(sqrt(ind * 2) + 0.5);
    row = ind - (col * col - col) / 2;
    
    
    
    for (int index = start; index < end; index++){
        temp[index-start]= row * col;
        row++;
        if (row > col){
            row = 1;
            col ++;
            
            if (prev_end > 1){ // merge this column with the previous merge
                
                int dup= serial_merge_no_dup(temp, temp+prev_end,otemp,prev_end, index - start-prev_end+ 1);
                end -= dup;
                index -= dup;
                
                long *tt = otemp;
                otemp = temp;
                temp = tt;
                
            }
            prev_end = index-start +1 ;
            
            
        }
        
    }
    
    if (prev_end > 1 && end  - start - prev_end > 0 ){  // merge this partial column with the previous merge
        int index = end-1;
      
        
        int dup= serial_merge_no_dup(temp, temp+prev_end,otemp,prev_end,index - start-prev_end+ 1);
        end -= dup;
        index -= dup;
        
        long *tt = otemp;
        otemp = temp;
        temp = tt;

    }
    

   
    newsize = end - start;
    
    

    for (int i=0; i<newsize;i++){
        buffer[i]=temp[i];
    }
    
    
    free (temp1);
    free (temp2);
    
    
    return newsize;
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Status status;

    long start, end, partition_size, num_upper_tri;
    long *a;
    int id, p;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    ID = id;

    long n = atol(argv[1]);

    num_upper_tri = (n*n-n)/2+n;

    // for each process, calculate their start, range, sort and remove duplicates
    partition_size = num_upper_tri/p;
    start = id*partition_size;
    
   

    
    if(id == (p-1)) {
        end = num_upper_tri;
    } else {
        end = start + partition_size;
    }
    
    long newsize = end-start;
    
    a = malloc(sizeof(long) * newsize);

    
    newsize = get_array(start, end, a);

    
    
    
   // printf("Process %d has array, start:%d, end:%d, size:%d, newsize:%d, n:%d \n", id, start, end,partition_size,newsize, n );
    
    
    char filename[100] ;
    
    sprintf(filename, "log2_%d.txt", id);
    FILE *f = fopen(filename, "w");
    
    for (int c = 0; c < newsize; c++){
        fprintf(f,"%d ", a[c]);
        
    }
    fclose(f); 
   
    
    
    

    int height, process, proc_recv;
    long mergesize, size2;
    long *a2;
    long *newa;

    height = log2(p);

   
    for (int i = 0; i < height; i++){
         MPI_Barrier(MPI_COMM_WORLD);

         process = p / pow(2,i+1);
        for (int j = 0; j < process; j++){
            //printf("Process %d merges into process %d at height %d of height %d \n", process+j, j, i, height);
            if(id == j){
                
                proc_recv = process+j;
                //printf("Process %d recieves from process %d at height %d of height %d \n", ID, proc_recv, i, height);
                MPI_Recv( &mergesize, 1, MPI_LONG, proc_recv, 1, MPI_COMM_WORLD, &status);
                //printf("Process %d recieved %d, tag 1 \n", ID, mergesize);
                
                
                size2 = newsize + mergesize;
                a2 = malloc(sizeof(long) * mergesize);
                MPI_Recv(a2, mergesize, MPI_LONG, proc_recv, 2, MPI_COMM_WORLD, &status);
               // printf("Process %d recieved %d longs, tag 1 \n", ID, mergesize);
                
                
                /*
                // merge the 2 arrays, sort, remove duplicates
                newa = malloc(sizeof(long) * size2);
                serial_merge(&a[0], &a2[0], &newa, newsize, mergesize);
                newsize=find_dups(size2,newa);
                */
                newa = malloc(sizeof(long) * size2);
                size2 -= serial_merge_no_dup(a,a2,newa, newsize, mergesize);
                
                free(a);
                free(a2);
                
                
                a= newa;
                newsize = size2;
                /*a = malloc(sizeof(long) * size2);
                for(int k = 0; k < newsize; k++){
                    a[k] = newa[k];
                }*/

            }else if(id == (process+j)) {
               // printf("Process %d sends to process %d at height %d of height %d \n", ID, j, i, height);
                
                MPI_Send(&newsize, 1, MPI_LONG, j, 1, MPI_COMM_WORLD);
                MPI_Send(a, newsize, MPI_LONG, j, 2, MPI_COMM_WORLD);
            }
        }
    }

    if(id == 0){
        printf("M(%d) = %d", n, newsize);
    }

    MPI_Finalize();
    return 0;
}

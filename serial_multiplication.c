#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

//used https://stackoverflow.com/questions/15821123/removing-elements-from-an-array-in-c
//as reference for
void remove_element(long *array, int index, int array_length)
{
    int i;
    for(i = index; i < array_length - 1; i++) array[i] = array[i + 1];
}

// Radix Sort
//https://austingwalters.com/radix-sort-in-c/
long findLargestNum(long * array, int size){

    int i;
    long largestNum = -1;

    for(i = 0; i < size; i++){
        if(array[i] > largestNum)
            largestNum = array[i];
    }

    return largestNum;
}

void radixSort(long * array, int size){

    // Base 10 is used
    int i;
    long semiSorted[size];
    int significantDigit = 1;
    long largestNum = findLargestNum(array, size);

    // Loop until we reach the largest significant digit
    while (largestNum / significantDigit > 0){


        long bucket[10] = { 0 };

        // Counts the number of "keys" or digits that will go into each bucket
        for (i = 0; i < size; i++)
            bucket[(array[i] / significantDigit) % 10]++;

        /**
         * Add the count of the previous buckets,
         * Acquires the indexes after the end of each bucket location in the array
             * Works similar to the count sort algorithm
         **/
        for (i = 1; i < 10; i++)
            bucket[i] += bucket[i - 1];

        // Use the bucket to fill a "semiSorted" array
        for (i = size - 1; i >= 0; i--)
            semiSorted[--bucket[(array[i] / significantDigit) % 10]] = array[i];


        for (i = 0; i < size; i++)
            array[i] = semiSorted[i];

        // Move to next significant digit
        significantDigit *= 10;

    }
}

int find_dups(int oldsize, long *a) {
    int newsize=oldsize;
    long current = a[0];
    long previous = 0;
    for (int i=0; i<newsize; i++){
        current = a[i];
        if (current == previous){
            remove_element(a,i,newsize);
            newsize--;
            i--;
        }
        previous = current;
    }
    return newsize;
}

int sort_array(long* a, int oldsize ){
    radixSort(a, oldsize);
    
//    for (int i=0;i<oldsize;i++){
//        printf("%d ",(int) a[i]);
//    }
//    printf("\n old size %d\n",oldsize);
    
    //sorted
    int newsize=oldsize;
    newsize=find_dups(oldsize,a);

    return newsize;
}

void serial_merge(long *a, long *b, long **c, long size_a, long size_b) {
    long i = 0;
    long j = 0;
    long k = 0;
    long b_j;
    long a_i;
    long n2 = size_a + size_b;

    *c=malloc(sizeof(long)*n2);
    while (k < n2 && i < size_a && j < size_b) {
        a_i = a[i];
        b_j = b[j];
        if (a_i <= b_j) {
            (*c)[k] = a_i;
            i++;
        } else {
            (*c)[k] = b_j;
            j++;
        }
        k++;
    }
    if (i >= size_a) {
        for (; j < size_b; j++) {
            b_j = b[j];
            (*c)[k] = b_j;
//j++;
            k++;
        }
    } else if (j >= size_b) {
        for (; i < size_a; i++) {
            a_i = a[i];
            (*c)[k] = a_i;
            //i++;
            k++;
        }
    }

}

/*
i = floor ( sqrt(index*2)+1/2)
j = index - (i*i-i)/2
*/
int get_array(long start, long end, int size, long* buffer){
    int newsize = 0;
    int i,j;
    long temp[size];
    for (int index=start; index<end; index++){
        int ind=index+1;
        i = floor(sqrt(ind * 2) + 0.5);
        j = ind - (i * i - i) / 2;
        temp[index-start]= i * j;
    }
    newsize = sort_array(temp, size);

    for (i=0; i<newsize;i++){
        buffer[i]=temp[i];
    }
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

    long n = atol(argv[1]);

    num_upper_tri = (n*n-n)/2+n;

    // for each process, calculate their start, range, sort and remove duplicates
    partition_size = num_upper_tri/p;
    start = id*partition_size;

    a = malloc(sizeof(long) * partition_size);

    if(id == (p-1)) {
        end = num_upper_tri;
    } else {
        end = start + partition_size;
    }

    int newsize = get_array(start, end, partition_size, &a[0]);

    printf("Process %d has array\n", id);
    // tree merge starting at the leaf
    // process array a(process+j) merge into process a(j)
    // if a(j) recv(a(j))
        // receive array, sizeof(array)
        // serial merge
        // remove duplicates
    // if a(process+j) send(a(process+j))
        // send array, newsize

    int height, process, proc_recv;
    long mergesize, size2;
    long *a2;
    long *newa;

    height = log2(p);

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < height; i++){
         process = p / pow(2,i+1);
        for (int j = 0; j < process; j++){
            printf("Process %d merges into process %d\n", process+j, j);
            if(id == j){
                proc_recv = process+j;
                MPI_Recv(&mergesize, 1, MPI_LONG, proc_recv, 1, MPI_COMM_WORLD, &status);
                size2 = newsize + mergesize;
                a2 = malloc(sizeof(long) * mergesize);
                MPI_Recv(a2, mergesize, MPI_LONG, proc_recv, 2, MPI_COMM_WORLD, &status);

                // merge the 2 arrays, sort, remove duplicates
                newa = malloc(sizeof(long) * size2);
                serial_merge(&a[0], &a2[0], &newa, newsize, mergesize);
                newsize=find_dups(size2,newa);

                free(a);
                a = malloc(sizeof(long) * size2);
                for(int k = 0; k < newsize; k++){
                    a[k] = newa[k];
                }

            }else if(id == (process+j)) {
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

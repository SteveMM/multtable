#include <stdio.h>
#include <math.h>
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

int sort_array(long* a, int oldsize ){
    radixSort(a, oldsize);
    for (int i=0;i<oldsize;i++){
        printf("%d ",(int) a[i]);
    }
    printf("\n old size %d\n",oldsize);
    //sorted
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

void serial_merge(long *a, long *b, long **c, long size_a, long size_b) {
    long i = 0;
    long j = 0;
    long k = 0;
    long b_j;
    long a_i;
    long n2 = size_a + size_b;

    printf("doing serial merge\n");

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
            j++;
            k++;
        }
    } else if (j >= size_b) {
        for (; i < size_a; i++) {
            a_i = a[i];
            (*c)[k] = a_i;
            i++;
            k++;
        }
    }

}


/*
i = floor ( sqrt(index*2)+1/2)
j = index - (i*i-i)/2
*/
int get_array(int start, int size, long* buffer){
    int newsize = 0;
    int end = start + size;
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

int main(void) {
    int size = 50;
    long buffer[50];
    int newsize = get_array(0, size, &buffer[0]);
    // printf("the number of unique numbers in first %d elements of upper triangular of multiplication table is %d\n", size, newsize);
    for (int i=0;i<newsize;i++){
        printf("%d ",(int) buffer[i]);
    }

    printf("\n new size %d\n",newsize);

    long buffer1[50];
    int newsize1 = get_array(52, size, &buffer1[0]);
    // printf("the number of unique numbers in first %d elements of upper triangular of multiplication table is %d\n", size, newsize);
    for (int i=0;i<newsize1;i++){
        printf("%d ",(int) buffer1[i]);
    }

    printf("\n new size %d\n",newsize1);

    printf("\n printing combined\n");

    //long *buffer2;
    //*buffer2 = malloc(sizeof(long)*(newsize+newsize1));
    long buffer2[86];
    serial_merge(&buffer[0],&buffer1[0],&buffer2[0],newsize,newsize1);

    for (int i=0;i<(newsize+newsize1);i++){
        printf("%d\n ",(int) buffer2[i]);
    }

    return 0;
}

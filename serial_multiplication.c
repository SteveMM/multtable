#include <stdio.h>
#include <math.h>


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
        }
        previous = current;
    }
    return newsize;
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
  return 0;
}

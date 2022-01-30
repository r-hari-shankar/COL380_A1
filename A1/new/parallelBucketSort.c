#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define max_threads 80

struct bucket {
    int n_elem;
    int index; // [start : n_elem)
    int start; //starting point in B array
};

int cmpfunc (const void * a, const void * b)
{
 return ( *(int*)a - *(int*)b );
}

int main(int argc, char *argv[]) {

    int *A, *B, *temp;
    int n, p, i, w, limit, num_threads, workload, b_index;
    struct bucket* buckets; //array of buckets
    
    double t1; // Timing variable
    float total; //total time

    printf("Give length of array to sort \n");
    if (scanf("%d", &n) != 1){
        printf("error\n");
        return -1;
    }
    printf("Give number of buckets \n");
    if (scanf("%d", &p) != 1){
        printf("error\n");
        return -1;
    }
	int *tmp;
    //global buckets
    int global_n_elem[p]; //number of elements in each bucket
    int global_starting_position[p]; //starting position in A for each bucket
    memset(global_n_elem, 0, sizeof(int)*p);
    memset(global_starting_position, 0, sizeof(int)*p);

    num_threads = p;
    omp_set_num_threads(num_threads);

    limit = 100000;
    w = limit/p;
    A = (int *) malloc(sizeof(int)*n);
    B = (int *) malloc(sizeof(int)*n);

    for(i=0;i<n;i++) {
	    A[i] = random() % limit;
	}
	
	if (n <= 40) {
        printf("Unsorted data \n");
        for(i=0;i<n;i++) {
            printf("%d ",A[i]);
        }
        printf("\n");
    }

    //local buckets, p for each thread
    buckets = (struct bucket *) calloc(p*num_threads, sizeof(struct bucket));

// ****************************
// Starting the main algorithm
// ****************************

    t1 = omp_get_wtime();
#pragma omp parallel
{
    num_threads = omp_get_num_threads();

    int j,k;
    int local_index; // [0 : p)
    int real_bucket_index; // [0 : p * num_threads)
    int tid = omp_get_thread_num();
    workload = n/num_threads;
	int pindex;

   #pragma omp for private(i,local_index)
    for (i=0; i< n;i++){
        local_index = A[i]/w;
        if (local_index > p-1)
                local_index = p-1;
        real_bucket_index = local_index + tid*p;
        buckets[real_bucket_index].n_elem++;
	}   
   
	int local_sum=0;
	for (j=tid; j< p*num_threads; j=j+num_threads){
		local_sum += buckets[j].n_elem;
	}
	global_n_elem[tid]=local_sum;
	
	#pragma omp barrier
	
	#pragma omp master
    {
	 for (j=1; j<p; j++){
		        global_starting_position[j] = global_starting_position[j-1] + global_n_elem[j-1];	 
				buckets[j].start = buckets[j-1].start + global_n_elem[j-1];
				buckets[j].index = buckets[j-1].index + global_n_elem[j-1];
		}
	}

	#pragma omp barrier
	for (j=tid+p; j< p*num_threads; j=j+num_threads){
		int pindex = j-p;
		buckets[j].start = buckets[pindex].start + buckets[pindex].n_elem;
		buckets[j].index = buckets[pindex].index + buckets[pindex].n_elem;	
	}
	#pragma omp barrier

    #pragma omp for private(i, b_index)
    for (i=0; i< n ;i++){
        j = A[i]/w;
        if (j > p -1)
                j = p-1;
        k = j + tid*p;
        b_index = buckets[k].index++;
        B[b_index] = A[i];
    }

#pragma omp for private(i)
    for(i=0; i<p; i++)
        qsort(B+global_starting_position[i], global_n_elem[i], sizeof(int), cmpfunc);
}
   total = omp_get_wtime() - t1;
	tmp = A;
	A = B;
	B = tmp;

   if (n <= 40) {
        printf("A \n");
        for(i=0;i<n;i++) {
            printf("%d ",A[i]);
        }
        printf("\n");
    }
    printf("Sorting %d elements took %f seconds\n", n,total);

    int sorted = 1;
    for(i=0;i<n-1;i++) {
        if (A[i] > A[i+1])
            sorted = 0;
    }
    if (!sorted)
        printf("The data is not sorted!!!\n");
}
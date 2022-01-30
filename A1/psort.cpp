#include "psort.h"
#include <omp.h>
#include <iostream>
#include <cstring>
struct bin{
    int n;
    int index;
    int start;
};
int cmpfunc (const void * a, const void * b)
{
 return ( *(int*)a - *(int*)b );
}
void ParallelSort(uint32_t *data, uint32_t n, int p)
{
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.
    uint32_t *B,*temp,workload;
    std::cout<<data[0]<<" "<<n<<" "<<p<<"\n";
    uint32_t bucket[p];
    //std::cout<<"Hi"<<"\n";
    struct bin* buckets;
    //std::cout<<"Hi"<<"\n";
    uint32_t start_bucket[p];
    uint32_t w=100000/p;
    B=(uint32_t *)malloc(sizeof(uint32_t)*n);
    //std::cout<<"Hi"<<"\n";
    memset(bucket,0,sizeof(uint32_t)*p);
    memset(start_bucket,0,sizeof(uint32_t)*p);    
    omp_set_num_threads(p);
    uint32_t bindex,i;
    //std::cout<<"Hi"<<"\n";
    buckets = (struct bin *) calloc(p*p,sizeof(struct bin));
    int numthreads;
    //std::cout<<"Hi"<<"\n";
    #pragma omp parallel
    {
        numthreads=p;
        int j,k,local_index,real_bucket_index,pindex,tid=omp_get_thread_num();
        workload = n/numthreads;
        #pragma omp for private(i,local_index)
            for(i=0;i<n;i++){
                local_index=data[i]/w;
                if(local_index>p-1)
                    local_index=p-1;
                real_bucket_index =local_index+tid*p;
                buckets[real_bucket_index].n++;
            }
        //std::cout<<"Hi"<<"\n";
        int local_sum=0;
        for(j=tid;j<p*numthreads;j=j+numthreads){
            local_sum+=buckets[j].n;
        }
        bucket[tid]=local_sum;
        #pragma omp barrier
        //std::cout<<"Hi"<<"\n";
        #pragma omp master
        {
            for(j=1;j<p;j++){
                start_bucket[j]=start_bucket[j-1]+bucket[j-1];
                buckets[j].start=buckets[j-1].start+bucket[j-1];
                buckets[j].index = buckets[j-1].index+bucket[j-1];
            }
        }
        //std::cout<<"Hi"<<"\n";
        #pragma omp barrier
        for(j=tid+p;j<p*numthreads;j+=numthreads){
            int pindex=j-p;
            buckets[j].start = buckets[pindex].start+buckets[pindex].n;
            buckets[j].index = buckets[pindex].index+buckets[pindex].n;
        }
        
        #pragma omp barrier
        std::cout<<"Hi"<<"\n";

        #pragma omp for private(i,bindex)
        for(i=0;i<n;i++){
                j=data[i]/w;
                if(j>p-1)
                    j = p-1;
                k=j+tid*p;
                bindex = buckets[k].index++;
                B[bindex]=data[i];
        }    
        
        std::cout<<"Hi"<<"\n";
    #pragma omp for private(i)
        for(i=0;i<p;i++)
            qsort(B+start_bucket[i],bucket[i],sizeof(int),cmpfunc);
    std::cout<<"Hi"<<"\n";
    }
    data = B;

    for(int i=0;i<10;i++){
        std::cout<<B[i]<<" ";
    }


}
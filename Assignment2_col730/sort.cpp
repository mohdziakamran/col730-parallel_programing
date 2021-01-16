#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include "sort.h"
#include "omp.h"

using namespace std;

#define dataType pSort::dataType 

#define base (1<<8)
#define numproc omp_get_num_procs()

/*-------------------definer-------------------*/
void seq_mergesort(dataType * data,int ndata);
void seq_merge(dataType * data, int ndata);
void omp_mergesort(dataType * data,int ndata,int num_procs);
void omp_merge(dataType * data, int ndata);

void seq_quicksort(dataType * data, int ndata);
void omp_quicksort(dataType * data, int ndata);

int radix_partition( dataType * data, int ndata);
void omp_radixsort(dataType * data, int ndata,int i,int shift, dataType * buf);
/*---------------------------------------------*/

void pSort::init(){
    
}

void pSort::close(){

}

void pSort::sort(dataType *data, int ndata, SortType sorter){

    int num_procs;int index0;

    switch (sorter)
    {
    case QUICK:
        /* code */
        // seq_quicksort(data,ndata);
        
        #pragma omp parallel num_threads(numproc)
        {
            #pragma omp single nowait
            {
                omp_quicksort(data,ndata);
            }
        }

        // for(int i=0;i<ndata;i++){
        //     cout<<data[i].key<<"\n";
        // }

        break;
    case MERGE:
        // seq_mergesort(data,ndata);
        #pragma omp parallel num_threads(numproc)
        {
            num_procs=omp_get_num_threads();
            #pragma omp single nowait
            {
                omp_mergesort(data,ndata,num_procs);
            }
        }
        break;
    case RADIX:
        //code
        index0=radix_partition(data,ndata);
        if(index0>0){
            
            dataType * buf=(dataType *)malloc(index0*sizeof(dataType));
            int total_digits=32;
            int i;
            for(int shift=0;shift<total_digits;shift+=8){
                omp_radixsort(data,index0,i,shift,buf);
            }
            free(buf);
            dataType * buff=(dataType *)malloc((ndata-index0)*sizeof(dataType));
            for(int shift=0;shift<total_digits;shift+=8){
                omp_radixsort(data+index0,ndata-index0,i,shift,buff);
            }
            free(buff);

            // cout<<"********************\n";
            // for(int i=0;i<ndata;i++){
            //     cout<<data[i].key<<"\n";
            // }
        }else{
            dataType * buf=(dataType *)malloc(ndata*sizeof(dataType));
            int total_digits=32;
            int i;
            for(int shift=0;shift<total_digits;shift+=8){
                omp_radixsort(data,ndata,i,shift,buf);
            }
            free(buf);
            // for(int i=0;i<ndata;i++){
            //     cout<<data[i].key<<"\n";
            // }
        }
        break;
    case BEST:
        sort(data,ndata,RADIX);
    default:
        //best
        break;
    }
}
void seq_mergesort(dataType * data,int ndata){
    if(ndata<2){
        return;
    }
    seq_mergesort(data,ndata/2);
    seq_mergesort(data+(ndata/2),ndata-(ndata/2));
    seq_merge(data,ndata);
}
void seq_merge(dataType * data,int ndata){
    int i=0;
    int j=ndata/2;
    int k=0;
    dataType * tmp=(dataType *)malloc((ndata)*sizeof(dataType));
    while(i<ndata/2&&j<ndata){
        if(data[i].key<=data[j].key){
            tmp[k]=data[i];
            k++;i++;
        }else{
            tmp[k]=data[j];
            j++;k++;
        }
    }
    while(i<ndata/2){
        tmp[k]=data[i];
        k++;i++;
    }
    while(j<ndata){
        tmp[k]=data[j];
        j++;k++;
    }
    memcpy (data,tmp,(ndata)*sizeof(dataType));
    free(tmp);
}
void omp_mergesort(dataType * data,int ndata,int num_procs){

    if(num_procs==1){seq_mergesort(data,ndata);return;}

    #pragma omp task 
    {
        omp_mergesort(data,ndata/2,num_procs/2);
    }
    #pragma omp task 
    {
        omp_mergesort(data+(ndata/2),ndata-(ndata/2),num_procs-(num_procs/2));
    }
    #pragma omp taskwait
    {
        omp_merge(data,ndata);
    }

}
void omp_merge(dataType * data, int ndata){
    int i=0;
    int j=ndata/2;
    int k=0;
    dataType * tmp=(dataType *)malloc(ndata*sizeof(dataType));
    while(i<ndata/2&&j<ndata){
        if(data[i].key<=data[j].key){
            tmp[k]=data[i];
            k++;i++;
        }else{
            tmp[k]=data[j];
            k++;j++;
        }
    }
    while(i<ndata/2){
        tmp[k]=data[i];
        k++;i++;
    }
    while(j<ndata){
        tmp[k]=data[j];
        k++;j++;
    }

    memcpy(data,tmp,ndata*sizeof(dataType));
    free(tmp);
}

void seq_quicksort(dataType * data, int ndata){
    if(ndata<2)return;
    int i=1;
    int j=ndata-1;
    int pivot=0;
    while(i<=j){
        while(i<ndata&&data[i].key<data[pivot].key){
            i++;
        }
        while(data[j].key>data[pivot].key){
            j--;
        }
        if(i<=j){
            dataType tmp=data[i];
            data[i]=data[j];
            data[j]=tmp;
            i++;j--;
        }
    }
    pivot=j;
    dataType tmp=data[0];
    data[0]=data[pivot];
    data[pivot]=tmp;

    seq_quicksort(data,pivot);
    seq_quicksort(data+pivot+1,ndata-pivot-1);
}
void omp_quicksort(dataType * data, int ndata){
    if(ndata<10000){seq_quicksort(data,ndata);return;}
    int i=1;
    int j=ndata-1;
    int pivot=0;
    while(i<=j){
        while( i<ndata && data[i].key<data[pivot].key){
            i++;
        }
        while(data[j].key>data[pivot].key){
            j--;
        }
        if(i<j){
            dataType tmp=data[i];
            data[i]=data[j];
            data[j]=tmp;
            i++;j--;
        }
    }
    pivot=j;
    dataType tmp=data[0];
    data[0]=data[pivot];
    data[pivot]=tmp;

    #pragma omp task
    {
        omp_quicksort(data,pivot);
    }
    #pragma omp task
    {
        omp_quicksort(data+pivot+1,ndata-pivot-1);
    }
    #pragma omp taskwait
    
}

int radix_partition( dataType * data, int ndata){
    int i=0;
    int j=ndata-1;
    int pivot=0;
    while(i<j){
        while(data[i].key<pivot){
            i++;
        }
        while(data[j].key>pivot){
            j--;
        }
        if(i<j){
            dataType tmp=data[i];
            data[i]=data[j];
            data[j]=tmp;
            if(data[i].key!=pivot)i++;
            if(data[j].key!=pivot)j--;
        }
    }
    return i;
}
void omp_radixsort(dataType * data, int ndata,int i,int shift, dataType * buf){
    // int base=1<<8;
    int mask=base-1;

    int bucket[base]={0};
    int local_bucket[base]={0};
    #pragma omp parallel num_threads(numproc) firstprivate(local_bucket)
    {
        // cout<<omp_get_num_threads()<<"..........\n";//////////////////////////////
        #pragma omp for schedule(static) nowait
        for(i=0;i<ndata;i++){
            local_bucket[((data[i].key>>shift) & mask)]++;
        }
        #pragma omp critical
        for(i=0;i<base;i++){
            bucket[i] +=local_bucket[i];
        }
        #pragma omp barrier
        #pragma omp single
        for(i=1;i<base;i++){
            bucket[i]+=bucket[i-1];
        }
        int num_threads=omp_get_num_threads();
        int tid=omp_get_thread_num();
        for(int j=num_threads-1;j>=0;j--){
            if(j==tid){
                for(i=0;i<base;i++){
                    bucket[i]-=local_bucket[i];
                    local_bucket[i]=bucket[i];
                }
            }else{
                #pragma omp barrier
            }
        }
        #pragma omp for schedule(static)
        for(i=0;i<ndata;i++){
            buf[local_bucket[((data[i].key>>shift) & mask)]++]=data[i];
        }
        // #pragma omp barrier
        // #pragma omp for schedule(static)
        // for(i=0;i<ndata;i++){
        //     data[i]=buf[i];
        // }
    }
    // memcpy(data,buf,ndata);
    #pragma omp parallel
    #pragma omp for schedule(static)
    for(i=0;i<ndata;i++){
        data[i]=buf[i];
    }

    // dataType * tmp=data;
    // data=buf;
    // buf=tmp;
    
}
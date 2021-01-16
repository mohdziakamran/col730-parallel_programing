#include "sort.h"
#include "mpi.h"
#include <bits/stdc++.h>

using namespace std;

// //---Global Variables-----
// int size, rank;

void msort(pSort::dataType *, int);
void _qsort(pSort::dataType *, int );
void rsort(pSort::dataType *, int );
void bsort(pSort::dataType *, int );

//-----------quicksort
void qsort_swap(pSort::dataType *,int ,int );
int qsort_partition(pSort::dataType *,int ,int );
void qsort_recc(pSort::dataType *,int ,int );
int pqsort_partition(pSort::dataType *,pSort::dataType *,pSort::dataType *,int , int,int *,int *);
pSort::dataType * pqsort_recc(pSort::dataType *data,MPI_Datatype mystruct,int *newdata_newlength,MPI_Comm mycomm);
//-----------mergesort func
void seq_merge(pSort::dataType *data,int l,int m,int r);
void p_mergelow(pSort::dataType *odata,int lenodata,pSort::dataType *rbuff,int lenrbuff);
void p_mergehi(pSort::dataType *odata,int lenodata,pSort::dataType *rbuff,int lenrbuff);
void seq_msort(pSort::dataType *data,int i,int j);
void pmsort(pSort::dataType *data,int ndata,int endrank,int height,MPI_Datatype mystruct,int rankshift);

void pSort::init()
{
    MPI_Init(NULL,NULL);
}

void pSort::close()
{
    MPI_Finalize();
}

void pSort::sort(dataType *data, int ndata, SortType sorter)
{
    // cout<<"hello world\n";
    switch (sorter)
    {
    case QUICK:
        /* code */
        _qsort(data,ndata);
        break;
    case MERGE:
        //code
        msort(data,ndata);
        break;
    case RADIX:
        //code
        rsort(data,ndata);
        break;
    default:
        //best
        bsort(data,ndata);
        break;
    }
}

void msort(pSort::dataType *data, int ndata){
    int size,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    seq_msort(data,0,ndata-1);
    
    //-----------mpi_struct------
    MPI_Datatype mystruct;
    int          blocklens[2];
    MPI_Aint     indices[2];
    MPI_Datatype old_types[2];

    blocklens[0] = 1;
    blocklens[1] = LOADSIZE;
    old_types[0] = MPI_INT;
    old_types[1] = MPI_CHAR;
    indices[1] = 4;
    indices[0] = 0;
    MPI_Type_struct( 2, blocklens, indices, old_types, &mystruct );
    MPI_Type_commit( &mystruct );
    //---------mpi_struct-END--------------
    
    // pSort::dataType *rbuff;
    // pSort::dataType *bbuff;
    // MPI_Status status;
    // int recv_count;

    // for(int i=0;i<size;i++){
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if(rank>=i){
    //         if(rank==i){
    //             bbuff=(pSort::dataType *)malloc(ndata*sizeof(pSort::dataType));
    //             for(int i=0;i<ndata;i++){
    //                 bbuff[i]=data[i];
    //             }
    //             recv_count=ndata;
    //         }
    //             MPI_Bcast(&recv_count,1,MPI_INT,i,MPI_COMM_WORLD);
    //         if(rank!=i){
    //             bbuff=(pSort::dataType *)malloc(recv_count*sizeof(pSort::dataType));
    //         }
    //             MPI_Bcast(bbuff,recv_count,mystruct,i,MPI_COMM_WORLD);
    //         if(rank==i){
    //             MPI_Probe(rank+1,0,MPI_COMM_WORLD,&status);
    //             MPI_Get_count(&status,mystruct,&recv_count);
    //             rbuff=(pSort::dataType *)malloc(recv_count*sizeof(pSort::dataType));
    //             MPI_Recv(rbuff,recv_count,mystruct,rank+1,0,MPI_COMM_WORLD,&status);
    //             p_mergelow(data,ndata,rbuff,recv_count);
                
    //             delete[] rbuff;

    //         }else if(rank==i+1){
    //             MPI_Send(data,ndata,mystruct,rank-1,0,MPI_COMM_WORLD);
    //             p_mergehi(data,ndata,bbuff,recv_count);
    //             delete[] bbuff;
    //         }else if(rank>i){
                
    //             p_mergehi(data,ndata,bbuff,recv_count);
    //             delete[] bbuff;
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    for(int i=0;i<size-1;i++){
        int height=ceil(log2(size-i));
        if(rank>=i){
            pmsort(data,ndata,size-1,height,mystruct,i);
        }

    }

//   // chek seq-merge-sort
//     if(rank==2){
//         for(int i=0;i<ndata;i++){
//             cout<<"rank="<<rank<<" "<<data[i].key<<"->"<<data[i].payload[0]<<data[i].payload[1]<<data[i].payload[2]<<data[i].payload[3]<<"\n";
//         }cout<<endl;
//     }

}
void _qsort(pSort::dataType *data, int ndata){
    if(true){rsort(data,ndata);return;}
    int MASTER=0;
    int rank,size,localsize;
    pSort::dataType *localdata;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    //-----------mpi_struct------
    MPI_Datatype mystruct;
    int          blocklens[2];
    MPI_Aint     indices[2];
    MPI_Datatype old_types[2];

    blocklens[0] = 1;
    blocklens[1] = LOADSIZE;
    old_types[0] = MPI_INT;
    old_types[1] = MPI_CHAR;
    indices[1] = 4;
    indices[0] = 0;
    MPI_Type_struct( 2, blocklens, indices, old_types, &mystruct );
    MPI_Type_commit( &mystruct );

    //---------mpi_struct-END--------------
    int newdata=ndata;
    pSort::dataType *sorteddata;
    sorteddata=pqsort_recc(data,mystruct,&newdata,MPI_COMM_WORLD);
    // cout<<newdata<<"---"<<ndata<<"---------\n";///////////////////////
    // for(int i=0;i<newdata;i++){
    //     cout<<data[i].key<<" - "<<(data+i)->payload[0]<<(data+i)->payload[1]<<(data+i)->payload[2]<<(data+i)->payload[3]<<"\n";
    // }

    // pSort::dataType *rbuff;
    // pSort::dataType *sbuff=sorteddata;

    // int *rcount=(int *)malloc(size*sizeof(int));
    // MPI_Gather(&newdata,1,MPI_INT,rcount,1,MPI_INT,0,MPI_COMM_WORLD);

    // int *displs=(int *)malloc(size*sizeof(int));
    // int n;
    // if (rank==0)
    // {
    //     displs[0]=0;
    //     for (int i = 0; i < size; i++)
    //     {
    //         displs[i]=displs[i-1]+rcount[i-1];
    //     }
    //     n=displs[size-1]+rcount[size-1];
    //     rbuff=(pSort::dataType *)malloc(n*sizeof(pSort::dataType));
    // }
    
    // MPI_Gatherv(sbuff,newdata,mystruct,rbuff,rcount,displs,mystruct,0,MPI_COMM_WORLD);

    // // if(rank==0){
    // //     cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~before~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    // //     // for(int i=0;i<n;i++){
    // //     //     cout<<rbuff[i].key<<" - "<<rbuff[i].payload[0]<<rbuff[i].payload[1]<<rbuff[i].payload[2]<<rbuff[i].payload[3]<<"\n";
    // //     // }
    // //     cout<<rcount[0]<<" "<<rcount[1]<<"\n";
    // // }

    // int *scount=(int *)malloc(size*sizeof(int));
    // MPI_Gather(&ndata,1,MPI_INT,scount,1,MPI_INT,0,MPI_COMM_WORLD);
    // if (rank==0)
    // {
    //     displs[0]=0;
    //     for (int i = 1; i < size; i++)
    //     {
    //         displs[i]=displs[i-1]+scount[i-1];
    //     }
    // }
    // sbuff=(pSort::dataType *)malloc(ndata*sizeof(pSort::dataType));
    // MPI_Scatterv(rbuff,scount,displs,mystruct,sbuff,ndata,mystruct,0,MPI_COMM_WORLD);
    // // for(int i=0;i<ndata;i++){
    // //     cout<<"@-"<<rank<<"-@"<<sbuff[i].key<<" - "<<sbuff[i].payload[0]<<sbuff[i].payload[1]<<sbuff[i].payload[2]<<sbuff[i].payload[3]<<"\n";
    // // }
    // for(int i=0;i<ndata;i++){
    //     *(data+i)=*(sbuff+i);
    // }

//  for(int i=0;i<ndata;i++){
//         cout<<"@-"<<rank<<"-@"<<data[i].key<<" - "<<data[i].payload[0]<<data[i].payload[1]<<data[i].payload[2]<<data[i].payload[3]<<"\n";
//     }

 for(int i=0;i<ndata;i++){
        cout<<"@-"<<rank<<"-@"<<sorteddata[i].key<<" - "<<sorteddata[i].payload[0]<<sorteddata[i].payload[1]<<sorteddata[i].payload[2]<<sorteddata[i].payload[3]<<"\n";
    }


}
void rsort(pSort::dataType *data, int ndata){

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int M=16,g=4;
    int *myhist=(int *)malloc(M*sizeof(int));
    for(int h=0;h<M;h++){myhist[h]=0;}
    vector<pSort::dataType> bucketlist[M];//=(vector<pSort::dataType> *)malloc(M*sizeof(vector<pSort::dataType>));
    
    /* distributing in buckets according to MSD */
    for(int i=0;i<ndata;i++){
        int v=data[i].key>>28 &15;
        if(data[i].key>>31 & 1){
            // cout<<"is neg v="<<v<<endl;//////////////////////////////////////////
            bucketlist[15-v].push_back(data[i]);
            myhist[15-v]++;
        }else{
            // cout<<"is pos v="<<v<<endl;///////////////////////////////////////////////
            bucketlist[8+v].push_back(data[i]);
            myhist[8+v]++;
        }
    }
    /* gather ing list of ndata from all*/
    int *ndatalist=(int *)malloc(size*sizeof(int));
    MPI_Allgather(&ndata,1,MPI_INT,ndatalist,1,MPI_INT,MPI_COMM_WORLD);

    /* gathering all myhist and making Global hist */
    int *histrbuff=(int *)malloc(M*size*sizeof(int));
    MPI_Allgather(myhist,M,MPI_INT,histrbuff,M,MPI_INT,MPI_COMM_WORLD);
    long *ghist=(long *)malloc(M*sizeof(long));
    for(int h=0;h<M;h++){ghist[h]=0;}
    for(int i=0;i<M;i++){
        for(int j=0;j<M*size;j+=M){
            ghist[i]+=histrbuff[i+j];
        }
    }
    
    // delete[] histrbuff;//-----------------------------------------------------------------------------

    /*count array to be sent or receive by every processor using Global histogram*/
    pair<int,int> *sendrecv_indexcount=(pair<int,int> *)malloc(size*sizeof(pair<int,int>));
    int j=0;
    for(int i=0;i<size;i++){
        if(j<M){
            sendrecv_indexcount[i].first=j;
            int val=ghist[j];
            j++;
            while(val<ndatalist[i]&&j<M){
                val+=ghist[j];
                j++;
            }
            sendrecv_indexcount[i].second=val;
        }else{
            sendrecv_indexcount[i].first=-1;
            sendrecv_indexcount[i].second=-1;
        }
    }


// //debuger
// for(int l=0;l<M;l++){/////////////////////////////////////
// cout<<"~"<<ghist[l];
// }cout<<endl;

//debuger
// for(int l=0;l<size;l++){/////////////////////////////////////
// cout<<"~("<<sendrecv_indexcount[l].first<<","<<sendrecv_indexcount[l].second<<")";
// }cout<<endl;


    //-----------mpi_struct------
    MPI_Datatype mystruct;
    int          blocklens[2];
    MPI_Aint     indices[2];
    MPI_Datatype old_types[2];

    blocklens[0] = 1;
    blocklens[1] = LOADSIZE;
    old_types[0] = MPI_INT;
    old_types[1] = MPI_CHAR;
    indices[1] = 4;
    indices[0] = 0;
    MPI_Type_struct( 2, blocklens, indices, old_types, &mystruct );
    MPI_Type_commit( &mystruct );

    //---------mpi_struct-END--------------



    /*do the distribution in forloop */
    pSort::dataType *sbuff;
    pSort::dataType *bcasrbuffrecv;
    int rbuff_count;
    int send_count;
    MPI_Status status;
    int *recvcountforgatherv;
    int totalcount_receivedinbcast;
    bool receivedbcast=false;
    int *dispsl;
    for(int i=0;i<size;i++){
        if(sendrecv_indexcount[i].first!=-1&&sendrecv_indexcount[i].second>0){

            int startindex_of_ghist=sendrecv_indexcount[i].first;
            int endindex_of_ghist=(i+1<size&&sendrecv_indexcount[i+1].first>0) ? (sendrecv_indexcount[i+1].first): M;
            //preparing sendbuff for bcast
            send_count=0;
            for(int k=startindex_of_ghist;k<endindex_of_ghist;k++){
                send_count+=myhist[k];
            }
            sbuff=(pSort::dataType *)malloc(send_count*sizeof(pSort::dataType));
            int p=0;
            for(int k=startindex_of_ghist;k<endindex_of_ghist;k++){
                for(int l=0;l<bucketlist[k].size();l++){
                    sbuff[p]=bucketlist[k][l];p++;
                }
            }

            if(rank==i){
                totalcount_receivedinbcast=sendrecv_indexcount[i].second;
                bcasrbuffrecv=(pSort::dataType *)malloc(totalcount_receivedinbcast*sizeof(pSort::dataType));
                recvcountforgatherv=(int *)malloc(size*sizeof(int));
            }
                MPI_Gather(&send_count,1,MPI_INT,recvcountforgatherv,1,MPI_INT,i,MPI_COMM_WORLD);
            if(rank==i){
                dispsl=(int *)malloc(size*sizeof(int));
                dispsl[0]=0;
                for(int k=1;k<size;k++){
                    dispsl[k]=recvcountforgatherv[k-1]+dispsl[k-1];
                }
            }

                MPI_Gatherv(sbuff,send_count,mystruct,bcasrbuffrecv,recvcountforgatherv,dispsl,mystruct,i,MPI_COMM_WORLD);
            if(rank==i){receivedbcast=true;}

        }
    }

    // delete[] sendrecv_indexcount;//---------------------------------------------------------------------
    // delete[] myhist;//-----------------------------------------------------------------------------
// cout<<"**rank "<<rank<<"bool "<<receivedbcast<<endl;///////////////////////////
/*Now every processor is having bcasr buffer to sort except some last one*/
pSort::dataType *newdata=bcasrbuffrecv;
int newdatalen=totalcount_receivedinbcast;
    pSort::dataType *fsendbuff;
    pSort::dataType *frecvbuff;
    int finalscount;
    int finalrcount;
    if(receivedbcast){
        seq_msort(newdata,0,newdatalen-1);
        // for(int i=0;i<newdatalen;i++){
        //         cout<<"~"<<newdata[i].key;
        //     }cout<<endl;
        if(rank==0){
            /*update localdata data*/
            for(int i=0;i<ndata;i++){data[i]=newdata[i];}
            /*make send buff*/
            finalscount=newdatalen-ndata;
            fsendbuff=(pSort::dataType *)malloc(finalscount*sizeof(pSort::dataType));
            for(int i=0;i<finalscount;i++){fsendbuff[i]=newdata[i+ndata];} 
            // for(int i=0;i<finalscount;i++){
            //     cout<<"~"<<fsendbuff[i].key;
            // }cout<<endl;
        }else{
            /*receve data from left*/
            MPI_Probe(rank-1,0,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status,mystruct,&finalrcount);
            frecvbuff=(pSort::dataType *)malloc(finalrcount*sizeof(pSort::dataType));
            MPI_Recv(frecvbuff,finalrcount,mystruct,rank-1,0,MPI_COMM_WORLD,&status);
            
            /*update localdata data and make senf buff*/
            if(finalrcount>ndata){
                /*updating data*/
                for(int i=0;i<ndata;i++){data[i]=frecvbuff[i];}
                /*makin send buff*/
                finalscount=newdatalen+finalrcount-ndata;
                fsendbuff=(pSort::dataType *)malloc(finalscount*sizeof(pSort::dataType));
                for(int i=ndata;i<finalrcount;i++){fsendbuff[i-ndata]=frecvbuff[i];}
                for(int i=0;i<newdatalen;i++){fsendbuff[i+finalrcount-ndata]=newdata[i];}

            }else if(finalrcount<ndata){
                /*updating localdata data*/
                for(int i=0;i<finalrcount;i++){data[i]=frecvbuff[i];}
                for(int i=finalrcount;i<ndata;i++){data[i]=newdata[i-finalrcount];}
                /*making send buffer*/
                finalscount=finalrcount+newdatalen-ndata;
                fsendbuff=(pSort::dataType *)malloc(finalscount*sizeof(pSort::dataType));
                for(int i=0;i<finalscount;i++){fsendbuff[i]=newdata[i+ndata-finalrcount];}

            }else{
                /*updating localdata data*/
                for(int i=0;i<ndata;i++){data[i]=frecvbuff[i];}
                /*making sendBuff*/
                finalscount=newdatalen;
                fsendbuff=(pSort::dataType *)malloc(finalscount*sizeof(pSort::dataType));
                for(int i=0;i<newdatalen;i++){fsendbuff[i]=newdata[i];}
            }

        }
        /*perform send action*/
        if(rank!=(size-1)){
            MPI_Send(fsendbuff,finalscount,mystruct,rank+1,0,MPI_COMM_WORLD);
        }
    }else{
        
        /*Receive data */
        MPI_Probe(rank-1,0,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status,mystruct,&finalrcount);
        frecvbuff=(pSort::dataType *)malloc(finalrcount*sizeof(pSort::dataType));
        MPI_Recv(frecvbuff,finalrcount,mystruct,rank-1,0,MPI_COMM_WORLD,&status);
        
        /*update localdata*/
        for(int i=0;i<ndata;i++){data[i]=frecvbuff[i];}

        /*make send bbuffer*/
        finalscount=finalrcount-ndata;
        fsendbuff=(pSort::dataType *)malloc(finalscount*sizeof(pSort::dataType));
        for(int i=0;i<finalscount;i++){fsendbuff[i]=frecvbuff[i+ndata];}
        if(rank<(size-1)){
            MPI_Send(fsendbuff,finalscount,mystruct,rank+1,0,MPI_COMM_WORLD);
            // for(int i=0;i<finalscount;i++){
            //     cout<<"~"<<fsendbuff[i].key;
            // }cout<<endl;
        }

    }



    /*Now every processor is having bcasr buffer to sort except some last one*/
    // if(receivedbcast){
    //     seq_msort(bcasrbuffrecv,0,totalcount_receivedinbcast-1);

    //     /*send excess*/
    //     pSort::dataType *final_sbuff;
    //     pSort::dataType *final_rbuff;
    //     int final_sendcount;
    //     int final_recvcount;
    //     if(rank!=0){
    //         MPI_Probe(rank-1,0,MPI_COMM_WORLD,&status);
    //         MPI_Get_count(&status,mystruct,&final_recvcount);
    //         final_rbuff=(pSort::dataType *)malloc(final_recvcount*sizeof(pSort::dataType));
    //         MPI_Recv(final_rbuff,final_recvcount,mystruct,rank-1,0,MPI_COMM_WORLD,&status);

    //         if(final_recvcount<ndata){
                
    //             for(int i=0;i<final_recvcount;i++){
    //                 data[i]=final_rbuff[i];
    //             }
    //             for(int i=final_recvcount,ii=0;i<ndata;i++){
    //                 data[i]=bcasrbuffrecv[ii];
    //                 ii++;
    //             }
    //             final_sendcount=final_recvcount+totalcount_receivedinbcast-ndata;
    //             final_sbuff=(pSort::dataType *)malloc(final_sendcount*sizeof(pSort::dataType));
    //             for(int i=totalcount_receivedinbcast-final_sendcount,ii=0;i<totalcount_receivedinbcast;i++){
    //                 final_sbuff[ii]=bcasrbuffrecv[i];
    //                 ii++;
    //             }

    //         }else if(final_recvcount==ndata){
    //             for(int i=0;i<final_recvcount;i++){
    //                 data[i]=final_rbuff[i];
    //             }
    //             final_sendcount=totalcount_receivedinbcast;
    //             final_sbuff=(pSort::dataType *)malloc(final_sendcount*sizeof(pSort::dataType));
    //             for(int i=0;i<totalcount_receivedinbcast;i++){
    //                 final_sbuff[i]=bcasrbuffrecv[i];
    //             }

    //         }else{
    //             for(int i=0;i<ndata;i++){
    //                 data[i]=final_rbuff[i];
    //             }
    //             final_sendcount=final_recvcount+totalcount_receivedinbcast-ndata;
    //             final_sbuff=(pSort::dataType *)malloc(final_sendcount*sizeof(pSort::dataType));
    //             for(int i=ndata,ii=0;i<final_recvcount;i++){
    //                 final_sbuff[ii]=final_rbuff[i];
    //                 ii++;
    //             }
    //             for(int i=final_recvcount-ndata,ii=0;ii<totalcount_receivedinbcast;i++){
    //                 final_sbuff[i]=bcasrbuffrecv[ii];
    //                 ii++;
    //             }

    //         }
            
    //     }else {
    //         for(int i=0;i<ndata;i++){
    //             data[i]=bcasrbuffrecv[i];
    //         }
    //         final_sendcount=totalcount_receivedinbcast-ndata;
    //         final_sbuff=(pSort::dataType *)malloc(final_sendcount*sizeof(pSort::dataType));
    //         for(int i=ndata,ii=0;i<totalcount_receivedinbcast;i++){
    //             final_sbuff[ii]=bcasrbuffrecv[i];
    //             ii++;
    //         }
    //     }
    //     if(rank!=size-1){
    //         MPI_Send(final_sbuff,final_sendcount,mystruct,rank+1,0,MPI_COMM_WORLD);
    //     }
        
    // }else{
    //     /*receive excess*/
    //     pSort::dataType *finaldata;
    //     int finalcount;
        
        
        
    //     MPI_Probe(rank-1,0,MPI_COMM_WORLD,&status);
    //     MPI_Get_count(&status,mystruct,&finalcount);
    //     finaldata=(pSort::dataType *)malloc(finalcount*sizeof(pSort::dataType));
    //     MPI_Recv(finaldata,finalcount,mystruct,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     for(int i=0;i<ndata;i++){
    //         data[i]=finaldata[i];
    //     }
    //     if(rank!=(size-1)){
    //         pSort::dataType *finalsend;
    //         int finalsendcount;
    //         finalsendcount=finalcount-ndata;
    //         finalsend=(pSort::dataType *)malloc(send_count*sizeof(pSort::dataType));
    //         for(int i=ndata,ii=0;i<finalcount;i++){
    //             finalsend[ii]=finaldata[i];
    //             ii++;
    //         }
    //         MPI_Send(finalsend,finalsendcount,mystruct,rank+1,0,MPI_COMM_WORLD);

    //     }

    // }



    
//             if(rank==i){
// //reinting bcastrecvcount[]//////////////////////////////////////////////////
// cout<<"bcastrcvcount\n";
// for(int m=0;m<size;m++){
// cout<<recvcountforgatherv[m]<<",";
// }cout<<endl;
// //reinting bcastrecvcount[]
// cout<<"^displs\n";
// for(int m=0;m<size;m++){
// cout<<dispsl[m]<<",";
// }cout<<endl;
//             }

// if(rank==i){
//     for(int l=0;l<totalcount_receivedinbcast;l++){/////////////////////////////////////
// cout<<"########"<<bcasrbuffrecv[l].key<<"~"<<bcasrbuffrecv[l].payload[0]<<bcasrbuffrecv[l].payload[1]<<bcasrbuffrecv[l].payload[2]<<bcasrbuffrecv[l].payload[3]<<endl;
// }cout<<endl;
// }

//         //debuger
// for(int l=0;l<totalcount_receivedinbcast;l++){/////////////////////////////////////
// cout<<"@@@@@@@@"<<bcasrbuffrecv[l].key<<"~"<<bcasrbuffrecv[l].payload[0]<<bcasrbuffrecv[l].payload[1]<<bcasrbuffrecv[l].payload[2]<<bcasrbuffrecv[l].payload[3]<<endl;
// }cout<<endl;

// debuger
// for(int l=0;l<ndata;l++){/////////////////////////////////////
// cout<<"@@@@@@@@"<<data[l].key<<"~"<<data[l].payload[0]<<data[l].payload[1]<<data[l].payload[2]<<data[l].payload[3]<<endl;
// }cout<<endl;
// cout<<"just atest----------------";















}
void bsort(pSort::dataType *data, int ndata){
    msort(data,ndata);
    return;
}

/*-------------------------------qsort funcs-----------------------------------------*/
void qsort_swap(pSort::dataType *data,int m,int n){
    pSort::dataType temp;
    temp=*(data+m);
    *(data+m)=*(data+n);
    *(data+n)=temp;
}
int qsort_partition(pSort::dataType *data,int first,int last){
    int pivot;
    int i,j;
    pivot=last;
    i=first-1;
    j=first;
    for(j=first;j<=last-1;j++){
        if((data+j)->key<(data+pivot)->key){
            i++;
            qsort_swap(data,i,j);
        }

    }
    qsort_swap(data,i+1,last);
    return i+1;    
}
void qsort_recc(pSort::dataType *data,int first,int last){
    int pivot;
    if(first>last){return;}
    pivot=qsort_partition(data,first,last);
    qsort_recc(data,first,pivot-1);
    qsort_recc(data,pivot+1,last);
}
int pqsort_partition(pSort::dataType *data,pSort::dataType *low,pSort::dataType *high,int pivot, int len,int *lowlen,int *highlen){
    // int pivot;
    int l=0,h=0;
    for(int i=0;i<len;i++){
        if((data+i)->key<=pivot){
            *(low+l)=(*(data+i));
            l++;
        }else
        {
            *(high+h)=(*(data+i));
            h++;
        }
    }
    *lowlen=l;
    *highlen=h;    
}
pSort::dataType * pqsort_recc(pSort::dataType *data,MPI_Datatype mystruct,int *newdata_newlength,MPI_Comm mycomm){
    int rank;
    int group_size;
    int len=*newdata_newlength;
    int mid;//=(mygend+mygstart)/2; 
    MPI_Comm_rank(mycomm,&rank);
    MPI_Comm_size(mycomm,&group_size);

    if(group_size==1){qsort_recc(data,0,len-1);return data;}

/*select pivot from 1st and bcast it to other member in group*/
    int pivotkey;
    if(rank==0){
        pivotkey=data[len-1].key;
    }
    MPI_Bcast(&pivotkey,1,MPI_INT,0,mycomm);

/* with respect to pivot divide localdata into low & high for every proceddor in group*/
    pSort::dataType *low;
    pSort::dataType *high;
    low=(pSort::dataType *)malloc(len*sizeof(pSort::dataType));
    high=(pSort::dataType *)malloc(len*sizeof(pSort::dataType));
    int lowlen=0,highlen=0;
    pqsort_partition(data,low,high,pivotkey,len,&lowlen,&highlen); 
    /*till here we partitioned every local data of involed processor*/

/*calculate |S| and |L| for global in group*/
    int S,L;
    MPI_Allreduce(&lowlen,&S,1,MPI_INT,MPI_SUM,mycomm);
    MPI_Allreduce(&highlen,&L,1,MPI_INT,MPI_SUM,mycomm);

/* define  Mid number of processor for left part involve*/
    mid=ceil(S*group_size/(S+L) +0.5);

// /*making low histogram and high histofram*/
//     int *lowhist=(int *)malloc(group_size*sizeof(int));
//     int *highhist=(int *)malloc(group_size*sizeof(int));

//             /*gather length in histogram*/
//             MPI_Allgather(&lowlen,1,MPI_INT,lowhist,1,MPI_INT,mycomm);
//             MPI_Allgather(&highlen,1,MPI_INT,highhist,1,MPI_INT,mycomm);

//             /*basis of divissin of length to be equal*/
//             int bod=S/(mid);
        
// /*assign pairs for the processor to receive (simille to rsort)*/
//     pair<int,int> *lowsendrecv_inexlength=(pair<int,int> *)malloc(mid*sizeof(pair<int,int>));
//     pair<int,int> *lowsendrecv_inexlength=(pair<int,int> *)malloc((group_size-mid)*sizeof(pair<int,int>));

//         /*for  low*/
//         int j=0;
//         for(int i=0;i<mid;i++){
//             lowsendrecv_inexlength[i].first=i;
//             int v=lowhist[j];
//         }

/*gather all S to rank 0 and all right to rank mid*/
    pSort::dataType *recvlow;
    pSort::dataType *recvhigh;
    int *recvcounts_arrlow=(int *)malloc(group_size*sizeof(int));
    int *recvcounts_arrhigh=(int *)malloc(group_size*sizeof(int));
    int *displslow=(int *)malloc(group_size*sizeof(int));
    int *displshigh=(int *)malloc(group_size*sizeof(int));

    /*recv_count[] for low*/
        MPI_Allgather(&lowlen,1,MPI_INT,recvcounts_arrlow,1,MPI_INT,mycomm);
    /*recv_count[] for high*/
        MPI_Allgather(&lowlen,1,MPI_INT,recvcounts_arrhigh,1,MPI_INT,mycomm);
    /*displs[] for low*/
        displslow[0]=0;
        for(int i=1;i<group_size;i++){
            displslow[i]=displslow[i-1]+recvcounts_arrlow[i-1];
        }
    /*displs[] for low*/
        displshigh[0]=0;
        for(int i=1;i<group_size;i++){
            displshigh[i]=displshigh[i-1]+recvcounts_arrhigh[i-1];
        }

    if(rank==0){
        /*make receive buffer for gather*/
        recvlow=(pSort::dataType *)malloc(S*sizeof(pSort::dataType));
    }
    if(rank==mid){
        /*make receive buffer for gather*/
        recvlow=(pSort::dataType *)malloc(S*sizeof(pSort::dataType));
    }
    /*gathering low in rank 0*/
        MPI_Gatherv(low,lowlen,mystruct,recvlow,recvcounts_arrlow,displslow,mystruct,0,mycomm);
    /*gathering high in rank mid*/
        MPI_Gatherv(high,highlen,mystruct,recvhigh,recvcounts_arrhigh,displshigh,mystruct,mid,mycomm);

/*split group*/
    int g_colour;
    MPI_Comm newcomm;
    // MPI_Comm newcommhigh;
    if(rank<mid){g_colour=0;}else{g_colour=1;}
    MPI_Comm_split(mycomm,g_colour,rank,&newcomm);
    // MPI_Comm_split(mycomm,g_colour,rank,&newcommhigh);

/*scatter in group*/
    int newrank;
    int newg_size;
    int avgs=S/newg_size;
    int avgl=L/newg_size;

    MPI_Comm_rank(newcomm,&newrank);
    MPI_Comm_size(newcomm,&newg_size);
    pSort::dataType *recvscatteredbuff;
    int *sendcounts_arrlow=(int *)malloc(newg_size*sizeof(int));
    int *sendcounts_arrhigh=(int *)malloc(newg_size*sizeof(int));
    delete[] displslow;
    delete[] displshigh;
    displslow=(int *)malloc(newg_size*sizeof(int));
    displshigh=(int *)malloc(newg_size*sizeof(int));

    if(rank==0){
            int v=S;
          for(int i=0;i<newg_size;i++){
              if(v>=avgs){sendcounts_arrlow[i]=avgs;v-=avgs;}else{sendcounts_arrlow[i]=v;v=0;}
          }
          displslow[0]=0;
          for(int i=1;i<newg_size;i++){
                displslow[i]=displslow[i-1]+sendcounts_arrlow[i-1];
            }
    }
    
    if(rank==mid){
        int v=S;
          for(int i=0;i<newg_size;i++){
              if(v>=avgs){sendcounts_arrhigh[i]=avgl;v-=avgl;}else{sendcounts_arrhigh[i]=v;v=0;}
          }
          displshigh[0]=0;
          for(int i=1;i<newg_size;i++){
                displshigh[i]=displshigh[i-1]+sendcounts_arrhigh[i-1];
            }
    }

    recvscatteredbuff=(pSort::dataType *)malloc(sendcounts_arrlow[newrank]*sizeof(pSort::dataType));

    if(rank<mid){
        MPI_Scatterv(recvlow,sendcounts_arrlow,displslow,mystruct,recvscatteredbuff,sendcounts_arrlow[newrank],mystruct,0,newcomm);
    }
    if(rank>=mid){
        MPI_Scatterv(recvhigh,sendcounts_arrhigh,displshigh,mystruct,recvscatteredbuff,sendcounts_arrhigh[newrank],mystruct,0,newcomm);
    }

/*run seperately recusrivly*/
    *newdata_newlength=sendcounts_arrlow[newrank];
    data=pqsort_recc(recvscatteredbuff,mystruct,newdata_newlength,newcomm);
    return data;




////////////////////////////////////////////////////////////////////////////////////////


    // pSort::dataType *recvlow;
    // pSort::dataType *recvhigh;
    // MPI_Status status;if(rank==mygstart){
    //     pivotkey=data[len-1].key;
    // }
    // int recv_count;

    // int newlength;
    // pSort::dataType *newdata;
    
    // if (rank<=mid)
    // {
    //     MPI_Type_commit( &mystruct );
        
    //     // if(rank==mid&&n%2==0){}///////////oddcase
    //     MPI_Send(high,highlen,mystruct,rank+(n/2)+1,0,MPI_COMM_WORLD);
    //     MPI_Probe(rank+(n/2)+1,0,MPI_COMM_WORLD,&status);
    //     MPI_Get_count(&status,mystruct,&recv_count);
    //     recvlow=(pSort::dataType *)malloc(recv_count*sizeof(pSort::dataType));
    //     MPI_Recv(recvlow,recv_count,mystruct,rank+(n/2)+1,0,MPI_COMM_WORLD,&status);
        
    //     newlength=lowlen+ recv_count;
    //     newdata=(pSort::dataType *)malloc(newlength*sizeof(pSort::dataType));
    //     int i=0;
    //     for(int j=0;j<lowlen;j++){
    //         newdata[i]=low[j];
    //         i++;
    //     }
    //     for(int j=0;j<recv_count;j++){
    //         newdata[i]=recvlow[j];
    //         i++;
    //     }
    //     delete recvlow;

    // }else
    // {
    //     MPI_Type_commit( &mystruct );
        
    //     MPI_Probe(rank-(n/2)-1,0,MPI_COMM_WORLD,&status);
    //     MPI_Get_count(&status,mystruct,&recv_count);
    //     recvhigh=(pSort::dataType *)malloc(recv_count*sizeof(pSort::dataType));
    //     MPI_Recv(recvhigh,recv_count,mystruct,rank-(n/2)-1,0,MPI_COMM_WORLD,&status);
    //     MPI_Send(low,lowlen,mystruct,rank-(n/2)-1,0,MPI_COMM_WORLD);

    //     newlength=highlen+ recv_count;
    //     newdata=(pSort::dataType *)malloc(newlength*sizeof(pSort::dataType));
    //     int i=0;
    //     for(int j=0;j<highlen;j++){
    //         newdata[i]=high[j];
    //         i++;
    //     }
    //     for(int j=0;j<recv_count;j++){
    //         newdata[i]=recvhigh[j];
    //         i++;
    //     }
    //     delete recvhigh;
    // }
    //     delete low;                        // free(low);
    //     // delete recvlow;                        // free(recvlow);
    //     delete high;                        // free(high);
    //     // delete(recvhigh);                        // free(recvhigh);
    
    // *newdata_newlength=newlength;
    // data=pqsort_recc(newdata,mygstart,mid,mystruct,newdata_newlength);
    // data=pqsort_recc(newdata,mid+1,mygend,mystruct,newdata_newlength);
    // return data;

}

/*-----------------------------msort funcs--------------------------------------------*/
void seq_merge(pSort::dataType *data,int l,int m,int r){
    int n1=m-l+1;
    int n2=r-m;
    pSort::dataType *L=(pSort::dataType *)malloc(n1*sizeof(pSort::dataType));
    pSort::dataType *R=(pSort::dataType *)malloc(n2*sizeof(pSort::dataType));
    for(int i=0;i<n1;i++){
        L[i]=data[l+i];
    }
    for(int i=0;i<n2;i++){
        R[i]=data[m+1+i];
    }
    int i=0,j=0,k=l;
    while (i<n1&&j<n2)
    {
        if (L[i].key<=R[j].key)
        {
            data[k]=L[i];
            i++;
        }else
        {
            data
            [k]=R[j];
            j++;
        }
        k++;
    }
    while (i<n1)
    {
        data[k]=L[i];
        i++;k++;
    }
    while (j<n2)   
    {
        data[k]=R[j];
        j++;k++;
    }

    delete[] L;
    delete[] R;
}
void p_mergelow(pSort::dataType *odata,int lenodata,pSort::dataType *rbuff,int lenrbuff){
    if(odata[lenodata-1].key<=rbuff[0].key){return;}
    pSort::dataType *odatacopy;
    odatacopy=(pSort::dataType *)malloc(lenodata*sizeof(pSort::dataType));
    for(int a=0;a<lenodata;a++){
        odatacopy[a]=odata[a];
    }
    int i=0;
    int j=0;
    int k=0;
    while(i<lenodata && j<lenrbuff && k<lenodata){
        if(odatacopy[i].key<=rbuff[j].key){
            odata[k]=odatacopy[i];
            i++;
        }else{
            odata[k]=rbuff[j];
            j++;
        }
        k++;
    }
    while(i<lenodata && k<lenodata){
        odata[k]=odatacopy[i];
        i++;k++;
    }
    while(j<lenrbuff && k<lenodata){
        odata[k]=rbuff[j];
        j++;k++;
    }
    delete[] odatacopy;
}
void p_mergehi(pSort::dataType *odata,int lenodata,pSort::dataType *rbuff,int lenrbuff){
    // for(int f=0;f<lenrbuff;f++){
    //     cout<<"#"<<rbuff[f].key<<",";
    // }
    // cout<<(odata[0].key>=rbuff[lenrbuff-1].key)<<endl;/////////
    if(odata[0].key>=rbuff[lenrbuff-1].key){return;}
    pSort::dataType *odatacopy;
    odatacopy=(pSort::dataType *)malloc(lenodata*sizeof(pSort::dataType));
    for(int a=0;a<lenodata;a++){
        odatacopy[a]=odata[a];
    }
    int i=lenodata-1;
    int j=lenrbuff-1;
    int k=lenodata-1;
    while(i>=0 && j>=0 && k>=0){
        if(odatacopy[i].key>=rbuff[j].key){
            odata[k]=odatacopy[i];
            i--;
        }else{
            odata[k]=rbuff[j];
            j--;
        }
        k--;
    }
    while(i>=0 && k>=0){
        odata[k]=odatacopy[i];
        i--;k--;
    }
    while(j>=0 && k>=0){
        odata[k]=rbuff[j];
        j--;k--;
    }
    delete[] odatacopy;
    
}
void seq_msort(pSort::dataType *data,int i,int j){
    if(i<j){
        int mid=(i+j)/2;
        seq_msort(data,i,mid);
        seq_msort(data,mid+1,j);
        seq_merge(data,i,mid,j);
        
    }
}
void pmsort(pSort::dataType *data,int ndata,int endrank,int height,MPI_Datatype mystruct,int rankshift){//height =ceil(log2(numprocs))

    int rank;
    int count;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int id=rank-rankshift;
    int parent,rightchild,myheight;
    pSort::dataType *half2,mergeresult;

    myheight=0;
    while(myheight<height){
        parent=(id & (~(1<<myheight)));

        if(parent==id){//left child
        
            rightchild=(id | (1<<myheight));
            if((rightchild+rankshift)<=endrank){
                MPI_Probe((rightchild+rankshift),0,MPI_COMM_WORLD,&status);
                MPI_Get_count(&status,mystruct,&count);
                half2=(pSort::dataType *)malloc(count*sizeof(pSort::dataType));
                MPI_Recv(half2,count,mystruct,(rightchild+rankshift),0,MPI_COMM_WORLD,&status);
                MPI_Send(data,ndata,mystruct,rightchild+rankshift,0,MPI_COMM_WORLD);
                p_mergelow(data,ndata,half2,count);

                free(half2);
            }
            myheight++;
        }else{//rightchild
            MPI_Send(data,ndata,mystruct,parent+rankshift,0,MPI_COMM_WORLD);
            MPI_Probe((parent+rankshift),0,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status,mystruct,&count);
            half2=(pSort::dataType *)malloc(count*sizeof(pSort::dataType));
            MPI_Recv(half2,count,mystruct,parent+rankshift,0,MPI_COMM_WORLD,&status);
            p_mergehi(data,ndata,half2,count);
            myheight=height;
            free(half2);
        }

    }
}

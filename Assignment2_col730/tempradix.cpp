private static final int CUTOFF = 15;
private static final int BITS_PER_INT = 32;
private static final int BITS_PER_BYTE = 8;
private static final int R = 256;

public void sort(int[] a){
    int firstPositiveIndex = partition(0, a, 0, a.length-1);
    int[] aux =new int[a.length];
    if(firstPositiveIndex>0){
        recSort(a, firstPositiveIndex, a.length-1, 0,aux);
        recSort(a, 0, firstPositiveIndex-1, 0,aux);
    }else{//all positive
        recSort(a, 0, a.length-1, 0, aux);
    }
}

private void recSort(int[] a, int lo, int hi, int d, int[] aux){
    if(d>4)return;
    if(hi-lo<CUTOFF){
        insertionSort(a,lo, hi);
        return;
    }

    int[] count = new int[R+1];

    //compute counts
    int bitsToShift = BITS_PER_INT-BITS_PER_BYTE*d-BITS_PER_BYTE;
    int mask = 0b1111_1111;
    for(int i = lo; i<=hi; i++){
        int c = (a[i]>>bitsToShift) & mask;
        count[c+1]++;
    }

    //compute indices
    for(int i = 0; i<R; i++){
        count[i+1]=count[i]+count[i+1];
    }

    //distribute
    for(int i = lo; i<=hi; i++){
        int c = (a[i]>>bitsToShift) & mask;
        aux[count[c]+lo] = a[i];
        count[c]++;
    }
    //copy back
    for(int i = lo; i<=hi; i++){
        a[i]=aux[i];
    }

    if(count[0]>0)
        recSort(a, lo, lo+count[0]-1, d+1, aux);
    for(int i = 1; i<R; i++){
        if(count[i]>0)
            recSort(a, lo+count[i-1], lo+count[i]-1, d+1, aux);
    }
}

// insertion sort a[lo..hi], starting at dth character
private void insertionSort(int[] a, int lo, int hi) {
    for (int i = lo; i <= hi; i++)
        for (int j = i; j > lo && a[j] < a[j-1]; j--)
            swap(a, j, j-1);
}


//returns the index of the partition or to the right of where it should be if the pivot is not in the array 
public int partition(int pivot, int[] a, int lo, int hi){
    int curLo = lo;
    int curHi = hi;
    while(curLo<curHi){
        while(a[curLo]<pivot){
            if((curLo+1)>hi)return hi+1;
            curLo++;
        }

        while(a[curHi]>pivot){
            if((curHi-1)<lo)return lo-1;
            curHi--;
        }
        if(curLo<curHi){
            swap(a, curLo, curHi);
            if(a[curLo]!=pivot)curLo++;
            if(a[curHi]!=pivot)curHi--;             
        }
    }
    return curLo;
}


private void swap(int[] a, int i1, int i2){
    int t = a[i1];
    a[i1]=a[i2];
    a[i2]=t;
}
/*********************************************************************************************/
void lsdradixsort(int* a, size_t n)
{
    // isolate integer byte by index.
    auto bmask = [](int x, size_t i)
    {
        return (static_cast<unsigned int>(x) >> i*8) & 0xFF;
    };

    // allocate temporary buffer.
    auto m = std::make_unique<int[]>(n);
    int* b = m.get();

    // for each byte in integer (assuming 4-byte int).
    for ( size_t i, j = 0; j < 4; j++ ) {
        // initialize counter to zero;
        size_t h[256] = {}, start;

        // histogram.
        // count each occurrence of indexed-byte value.
        for ( i = 0; i < n; i++ )
            h[bmask(a[i], j)]++;

        // accumulate.
        // generate positional offsets. adjust starting point
        // if most significant digit.
        start = (j != 3) ? 0 : 128;
        for ( i = 1+start; i < 256+start; i++ )
            h[i % 256] += h[(i-1) % 256];

        // distribute.
        // stable reordering of elements. backward to avoid shifting
        // the counter array.
        for ( i = n; i > 0; i-- )
            b[--h[bmask(a[i-1], j)]] = a[i-1];
        std::swap(a, b);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

#define BASE_BITS 8
#define BASE (1 <> shift) & MASK)
 
void omp_lsd_radix_sort(size_t n, unsigned data[n]) {
    unsigned * buffer = malloc(n*sizeof(unsigned));
    int total_digits = sizeof(unsigned)*8;
 
    //Each thread use local_bucket to move data
    size_t i;
    for(int shift = 0; shift < total_digits; shift+=BASE_BITS) {
        size_t bucket[BASE] = {0};
 
        size_t local_bucket[BASE] = {0}; // size needed in each bucket/thread
        //1st pass, scan whole and check the count
        #pragma omp parallel firstprivate(local_bucket)
        {
            #pragma omp for schedule(static) nowait
            for(i = 0; i < n; i++){
                local_bucket[DIGITS(data[i], shift)]++;
            }
            #pragma omp critical
            for(i = 0; i < BASE; i++) {
                bucket[i] += local_bucket[i];
            }
            #pragma omp barrier
            #pragma omp single
            for (i = 1; i = 0; cur_t--) {
                if(cur_t == tid) {
                    for(i = 0; i < BASE; i++) {
                        bucket[i] -= local_bucket[i];
                        local_bucket[i] = bucket[i];
                    }
                } else { //just do barrier
                    #pragma omp barrier
                }
 
            }
            #pragma omp for schedule(static)
            for(i = 0; i < n; i++) { //note here the end condition
                buffer[local_bucket[DIGITS(data[i], shift)]++] = data[i];
            }
        }
        //now move data
        unsigned* tmp = data;
        data = buffer;
        buffer = tmp;
    }
    free(buffer);
}
/*
nvevweoib wv ifjjcnjcnisddc  dsnnvnjksnvjsnckchabvufjfivsnvuueuinejnv
bhubvuhvhbvhbhsbdhshcjsdcj dhcbbsjajjviajoqwjnfuwfwqifuvfhbvuhsdihbcfvhs
ibqbbahhbidqiwbciqwncqwic
cbiuqvqjrq

*/

/*
Copyright (c) 2014, Haichuan Wang
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Code for the blog https://haichuanwang.wordpress.com/2014/05/26/a-faster-openmp-radix-sort-implementation/ */

#define BASE_BITS 8
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
 
void omp_lsd_radix_sort(size_t n, unsigned data[n]) {
    unsigned * buffer = malloc(n*sizeof(unsigned));
    int total_digits = sizeof(unsigned)*8;
 
    //Each thread use local_bucket to move data
    size_t i;
    for(int shift = 0; shift < total_digits; shift+=BASE_BITS) {
        size_t bucket[BASE] = {0};
 
        size_t local_bucket[BASE] = {0}; // size needed in each bucket/thread
        //1st pass, scan whole and check the count
        #pragma omp parallel firstprivate(local_bucket)
        {
            #pragma omp for schedule(static) nowait
            for(i = 0; i < n; i++){
                local_bucket[DIGITS(data[i], shift)]++;
            }
            #pragma omp critical
            for(i = 0; i < BASE; i++) {
                bucket[i] += local_bucket[i];
            }
            #pragma omp barrier
            #pragma omp single
            for (i = 1; i < BASE; i++) {
                bucket[i] += bucket[i - 1];
            }
            int nthreads = omp_get_num_threads();
            int tid = omp_get_thread_num();
            for(int cur_t = nthreads - 1; cur_t >= 0; cur_t--) {
                if(cur_t == tid) {
                    for(i = 0; i < BASE; i++) {
                        bucket[i] -= local_bucket[i];
                        local_bucket[i] = bucket[i];
                    }
                } else { //just do barrier
                    #pragma omp barrier
                }
 
            }
            #pragma omp for schedule(static)
            for(i = 0; i < n; i++) { //note here the end condition
                buffer[local_bucket[DIGITS(data[i], shift)]++] = data[i];
            }
        }
        //now move data
        unsigned* tmp = data;
        data = buffer;
        buffer = tmp;
    }
    free(buffer);
}



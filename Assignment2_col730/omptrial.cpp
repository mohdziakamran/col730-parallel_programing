#include "omp.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    int n, m;
    int c=0;
    n=omp_get_num_procs();
    #pragma omp parallel 
    {
         m=omp_get_num_threads();
         #pragma omp critical
         {
             c++;
         }
         
    }
    cout<<n<<"\n";
    cout<<m<<"\n";
    cout<<c<<"\n";
    return 0;
}

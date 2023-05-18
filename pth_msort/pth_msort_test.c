//C header files
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "timer.h"
#include "pth_msort.h"

//Dynamically allocate an array of size N in heap
void fillArray ( int** values, unsigned int N )
{
    srand(time(0));
	unsigned int i;
    for (i=0; i<N; i++)
        (*values)[i] = rand() - (RAND_MAX/2);
}

//Check to see if the sorted array is really sorted!
int checkArray (const int* values, unsigned int N, const int* sorted)
{
    unsigned int i;
	for (i=0; i <= N-2; i++)
        if (sorted[i]>sorted[i+1])
        {
		printf("ERROR i=%d %d %d", i, sorted[i], sorted[i+1]);
		return 0;
        }
	return 1;
}

int main ( int argc, char* argv[] )
{
    if ( argc != 2 ){
        printf("Please enter: \n ./a.out  M\n");
        exit( 1 );
    }

	//input array
	unsigned int N;

	//Fill in the input array with random data
	N = pow( 2, atoi(argv[1]) );
	int* values = (int*) malloc ( sizeof(int) * N );
	fillArray ( &values, N );

	//Sorted array
	int* sorted = (int*) malloc ( sizeof(int) * N );

	//Sort the input array and report its runtime
	double t1,t2;
	GET_TIME(t1);
	mergeSortParallel (values, N, sorted);
	GET_TIME(t2);
	

    //Check for correctness
	if (checkArray ( values, N, sorted ))
		printf("Runtime= %g s N=%d\n", t2-t1,N);
	else{
		printf("Runtime= %g s N=%d\n", t2-t1,N);
		printf("Error");
	}

	return 0;
}

//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "bmm.h"

#define tx threadIdx.x
#define ty threadIdx.y

#define bx blockIdx.x
#define by blockIdx.y

// TILEX and TILEY are used to set the number of threads in a CUDA block 
#define TILEX 32
#define TILEY 16

// you may define other parameters here!
// it's lower than tilex and tiley or bigger than both of them
#define TILE 128
// you may define other macros here!
// you may define other functions here!

dim3 getDimGrid(const int m, const int n) {
	dim3 dimGrid(n / TILEX, n / TILEY);
	return dimGrid;
}
dim3 getDimBlock(const int m, const int n) {
	dim3 dimBlock(TILEX, TILEY);
	return dimBlock;
}

__global__ void kernelFunc(float* ad, float* bd, float* cd, const int m, const int n) {

	__shared__ float mad[TILEY][TILE];
	__shared__ float mbd[TILE][TILEX];

	int Row = by * TILEY + ty;
	int Col = bx * TILEX + tx;
	float Pvalue = 0;

	for (int i = 0; i < n / TILE; ++i) {

		if (TILE <= TILEX){
			if (tx < TILE)
				mad[ty][tx] = ad[Row * n + i * TILE + tx];
		}
		else
			for (int j = 0; j < TILE/TILEX; ++j)
				mad[ty][tx + j*TILEX] = ad[Row * n + ((TILE/TILEX)*i + j) * TILEX + tx];

		if (TILE <= TILEY){
			if (ty < TILE)
				mbd[ty][tx] = bd[(i * TILE + ty) * n + Col];
		}	
		else	
			for (int j = 0; j < TILE / TILEY; ++j) 
				mbd[ty + j * TILEY][tx] = bd[(((TILE / TILEY) * i + j) * TILEY + ty) * n + Col];
				

		__syncthreads();

		for (int k = 0; k < TILE; ++k) {
			Pvalue += mad[ty][k] * mbd[k][tx];
		}
		__syncthreads();
	}
	cd[Row * n + Col] = Pvalue;
}

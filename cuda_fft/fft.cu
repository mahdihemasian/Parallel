//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "fft.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

//-----------------------------------------------------------------------------

__global__ void sort_rad4(float* input, float* input1, unsigned int M, int k) 
{

    int i = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;

    unsigned int j = i;
    j = ((j & 0xcccccccc) >> 2) | ((j & 0x33333333) << 2);
    j = ((j & 0xf0f0f0f0) >> 4) | ((j & 0x0f0f0f0f) << 4);
    j = ((j & 0xff00ff00) >> 8) | ((j & 0x00ff00ff) << 8);
    j = (j >> 16) | (j << 16);
    j >>= 32-M;

    input[j]=input1[i + k];
    
}


__global__ void fft_rad2 (float* x_r_d, float* x_i_d ,const unsigned int N, int M){

    int i = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;

    float angel = -1 * 2 * PI * ((i*(N/(2*M)))-(i/M)*(N/2)) / N;
	float Wr = cos(angel);                   
	float Wi = sin(angel);

	float x1_r, x2_r, x1_i, x2_i;
	
	x1_r = x_r_d[i+(i/M)*M];         
	x2_r = x_r_d[i+(i/M)*M+(M)];
	
	x1_i = x_i_d[i+(i/M)*M];
	x2_i = x_i_d[i+(i/M)*M+(M)];
   
	x_r_d[i+(i/M)*M] = x1_r + Wr * x2_r - Wi * x2_i;
	x_i_d[i+(i/M)*M] = x1_i + Wr * x2_i + Wi * x2_r;
	
	x_r_d[i+(i/M)*M+(M)] = x1_r - Wr * x2_r + Wi * x2_i;
	x_i_d[i+(i/M)*M+(M)] = x1_i - Wr * x2_i - Wi * x2_r;		
		
}


__global__ void fft_rad4 (float* x_r_d, float* x_i_d, const unsigned int N, int M,unsigned int k)
{

    int i = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x + tx;

	float x1_r, x2_r, x1_i, x2_i, x4_r, x3_r, x4_i, x3_i;	
	float y2_i, y3_i, y4_i, y2_r, y3_r, y4_r;
	float teta  = -2*PI*(i%M) / (M*4);	
	
	
	x1_r = x_r_d[(i/M)*(4*M)+(i%M) + k];      
	x2_r = x_r_d[(i/M)*(4*M)+(i%M) + M + k];
	x3_r = x_r_d[(i/M)*(4*M)+(i%M) + 2*M + k];
	x4_r = x_r_d[(i/M)*(4*M)+(i%M) + 3*M + k];
	
	x1_i = x_i_d[(i/M)*(4*M)+(i%M) + k];
	x2_i = x_i_d[(i/M)*(4*M)+(i%M) + M + k];
	x3_i = x_i_d[(i/M)*(4*M)+(i%M) + 2*M + k];
	x4_i = x_i_d[(i/M)*(4*M)+(i%M) + 3*M + k];	
	
	float aaa = cos(teta);
	float bbb = sin(teta);

	y2_r = x2_r * aaa - x2_i * bbb;
	y2_i = x2_r * bbb + x2_i * aaa;

    aaa = cos(2*teta);
	bbb = sin(2*teta);

	y3_r = x3_r * aaa - x3_i * bbb;
	y3_i = x3_r * bbb + x3_i * aaa;

    aaa = cos(3*teta);
	bbb = sin(3*teta);

	y4_r = x4_r * aaa - x4_i * bbb;
	y4_i = x4_r * bbb + x4_i * aaa;	
	
	
	x_r_d[(i/M)*(M*4)+(i%M) + k] = x1_r + y2_r + y3_r + y4_r;
	x_i_d[(i/M)*(M*4)+(i%M) + k] = x1_i + y2_i + y3_i + y4_i;
	
	x_r_d[(i/M)*(M*4)+(i%M) + M + k] = x1_r + y2_i - y3_r - y4_i;
	x_i_d[(i/M)*(M*4)+(i%M) + M + k] = x1_i - y2_r - y3_i + y4_r;
	
	x_r_d[(i/M)*(M*4)+(i%M) + 2*M + k] = x1_r - y2_r + y3_r - y4_r;
	x_i_d[(i/M)*(M*4)+(i%M) + 2*M + k] = x1_i - y2_i + y3_i - y4_i;
	
	x_r_d[(i/M)*(M*4)+(i%M) + 3*M + k] = x1_r - y2_i - y3_r + y4_i;
	x_i_d[(i/M)*(M*4)+(i%M) + 3*M + k] = x1_i + y2_r - y3_i - y4_r;
	
}


__global__ void transfer(float* x, float* temp, float* x1, float* temp1) {

	int i = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;

	x[i] = temp[i];
	x1[i] = temp1[i];

}


__global__ void transfer1(float* x, float* temp, int k) {

	int i = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;

	x[i + k] = temp[i];
}


__global__ void transpose(float* x, float* tmp, float* x1, float* tmp1, const unsigned int N) 
{
    int i = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;

    if(i % 2 == 0)
    {
        tmp[i/2] = x[i];
        tmp1[i/2] = x1[i];
    }
    else
    {
        tmp[i/2+N/2] = x[i];
        tmp1[i/2+N/2] = x1[i];
    }

}


void sort_even_number(float* x_r_d, float* x_i_d ,const unsigned int N,const unsigned int M, int k)
{
    float* tmp;

    cudaMalloc((void**)&tmp, sizeof(float) * N);

    dim3 dimGrid((N / (512*512)), 32, 32);
	dim3 dimBlock(256, 1, 1);

    sort_rad4 <<< dimGrid, dimBlock >>>(tmp, x_r_d, M, k);
    transfer1 <<< dimGrid, dimBlock >>>(x_r_d,tmp, k);
    sort_rad4 <<< dimGrid, dimBlock >>>(tmp, x_i_d, M, k);
    transfer1 <<< dimGrid, dimBlock >>>(x_i_d,tmp, k);

    cudaFree(tmp);

}


void sort_odd_number(float* x_r_d, float* x_i_d ,const unsigned int N,const unsigned int M)
{
    float* tmp_r;
    float* tmp_i;

    cudaMalloc((void**)&tmp_r, sizeof(float) * N);
    cudaMalloc((void**)&tmp_i, sizeof(float) * N);

    dim3 dimGrid(N/1024, 1, 1);
	dim3 dimBlock(1024, 1, 1);
    transpose<<<dimGrid , dimBlock>>>(x_r_d, tmp_r, x_i_d, tmp_i, N);
    transfer<<<dimGrid , dimBlock>>>(x_r_d, tmp_r, x_i_d, tmp_i);

    cudaFree(tmp_r);
    cudaFree(tmp_i);

    sort_even_number(x_r_d, x_i_d, N/2, M-1, 0);
    sort_even_number(x_r_d, x_i_d, N/2, M-1, N/2);

}



void gpuKernel(float* x_r_d, float* x_i_d, const unsigned int N, const unsigned int M)
{
 	
	if(M%2 == 0)
    {
        sort_even_number(x_r_d, x_i_d, N, M, 0);

        dim3 dimGrid((N / (16*256)), 8, 1);
	    dim3 dimBlock(128, 1, 1);

	    for (int i=1; i<N; i*=4)  
	    {
	        fft_rad4 <<< dimGrid, dimBlock >>>(x_r_d, x_i_d, N, i, 0);
		}

    }

    else
    {
	    sort_odd_number(x_r_d, x_i_d, N, M);

        dim3 dimGrid((N / (32*256)), 8, 1);
	    dim3 dimBlock(128, 1, 1);

        for (int i=1; i<N/2; i*=4)  
	    {
	        fft_rad4 <<< dimGrid, dimBlock >>>(x_r_d, x_i_d, N, i, 0);
	        fft_rad4 <<< dimGrid, dimBlock >>>(x_r_d, x_i_d, N, i, N/2);
		}

        dim3 dimGrid1((N / (1024*256*2)), 32, 32);
	    dim3 dimBlock1(256, 1, 1);
        fft_rad2 <<< dimGrid1, dimBlock1 >>>(x_r_d, x_i_d, N, N/2);
        
    }

	
}
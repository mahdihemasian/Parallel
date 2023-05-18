//Do NOT MODIFY THIS FILE

#include "fft.h"

// ===========================> Functions Prototype <===============================
void fill(float* data, int size);
double calc_mse(float* data1_r, float* data1_i, float* data2_r, float* data2_i, int size);
void get_inputs(int argc, char *argv[], unsigned int& N, unsigned int& M);
void cpuKernel(float* X_serial_r, float* X_serial_i, int n, float* tmp_r, float* tmp_i);
void gpuKernels(float* x_r, float* x_i, float* X_r, float* X_i, unsigned int N, unsigned int M, double* gpu_kernel_time);
// =================================================================================

int main(int argc, char *argv[]) {

    struct cudaDeviceProp p;
    cudaGetDeviceProperties(&p, 0);
    printf("Device Name: %s\n", p.name);

    // get parameters from command line
    unsigned int N, M;
    get_inputs(argc, argv, N, M);

    // allocate memory in CPU for calculation
    float* x_r; // real part
    float* x_i; // imaginary part
    float* X_serial_r;
    float* X_serial_i;
    float* X_r;
    float* X_i;
    x_r = (float*) malloc(N * sizeof(float));
    x_i = (float*) malloc(N * sizeof(float));
    X_serial_r = (float*) malloc(N * sizeof(float));
    X_serial_i = (float*) malloc(N * sizeof(float));
    X_r = (float*) malloc(N * sizeof(float));
    X_i = (float*) malloc(N * sizeof(float));

    // fill x_r and x_i arrays with random values between -8.0f and 8.0f
    srand(0);
    fill(x_r, N);
    fill(x_i, N);
	int i; for (i = 0; i < N; i++) {
		X_serial_r[i] = x_r[i];
		X_serial_i[i] = x_i[i];
	}

    // time measurement for CPU calculation
	float *tmp_r, *tmp_i;
	tmp_r = (float*) malloc(N * sizeof(float));
    tmp_i = (float*) malloc(N * sizeof(float));
    clock_t t0 = clock();
    cpuKernel(X_serial_r, X_serial_i, N, tmp_r, tmp_i);
    clock_t t1 = clock();
	free(tmp_r); free(tmp_i);

    // time measurement for GPU calculation
	double gpu_kernel_time = 0.0;
    clock_t t2 = clock();
	gpuKernels(x_r, x_i, X_r, X_i, N, M, &gpu_kernel_time);
    clock_t t3 = clock();

    // check correctness of calculation
    double mse = calc_mse(X_serial_r, X_serial_i, X_r, X_i, N);
	printf("m=%d n=%d CPU=%g ms GPU=%g ms GPU-Kernels=%g ms mse=%g\n",
	M, N, (t1-t0)/1000.0, (t3-t2)/1000.0, gpu_kernel_time, mse);
	
	/*
	for (i = 0; i<N; i++) {
		printf("%f\t%f\n", x_r[i], x_i[i]);
	}
	printf("\n");
	for (i = 0; i<N; i++) {
		printf("%f\t%f\n", X_serial_r[i], X_serial_i[i]);
	}
	*/
	
    // free allocated memory for later use
    free(x_r);
    free(x_i);
    free(X_serial_r);
    free(X_serial_i);
    free(X_r);
    free(X_i);

    return 0;
}

//-----------------------------------------------------------------------------
void gpuKernels(float* x_r, float* x_i, float* X_r, float* X_i, unsigned int N, unsigned int M, double* gpu_kernel_time) {
    float* x_r_d;
    float* x_i_d;
    //float* X_r_d;
    //float* X_i_d;

    HANDLE_ERROR(cudaMalloc((void**)&x_r_d, N * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&x_i_d, N * sizeof(float)));
    //HANDLE_ERROR(cudaMalloc((void**)&X_r_d, N * sizeof(float)));
    //HANDLE_ERROR(cudaMalloc((void**)&X_i_d, N * sizeof(float)));

    HANDLE_ERROR(cudaMemcpy(x_r_d, x_r, N * sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(x_i_d, x_i, N * sizeof(float), cudaMemcpyHostToDevice));

	GpuTimer timer;
    timer.Start();
	gpuKernel(x_r_d, x_i_d, N, M);/*<<<dim3(32,1,1),dim3(32,1,1)>>>(x_r_d, x_i_d, N, M);*/
	timer.Stop();
	*gpu_kernel_time = timer.Elapsed();
	
    //HANDLE_ERROR(cudaMemcpy(X_r, X_r_d, N * sizeof(float), cudaMemcpyDeviceToHost));
    //HANDLE_ERROR(cudaMemcpy(X_i, X_i_d, N * sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(X_r, x_r_d, N * sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(X_i, x_i_d, N * sizeof(float), cudaMemcpyDeviceToHost));

    HANDLE_ERROR(cudaFree(x_r_d));
    HANDLE_ERROR(cudaFree(x_i_d));
    //HANDLE_ERROR(cudaFree(X_r_d));
    //HANDLE_ERROR(cudaFree(X_i_d));
}
//-----------------------------------------------------------------------------
void cpuKernel(float* X_serial_r, float* X_serial_i, int n, float* tmp_r, float* tmp_i) {
	if(n > 1) {	// otherwise, do nothing and return
		int k, m;
		float z_r, z_i, w_r, w_i;
		float *vo_r, *vo_i, *ve_r, *ve_i;
		ve_r = tmp_r; ve_i = tmp_i;
		vo_r = tmp_r + n/2; vo_i = tmp_i + n/2;
		
		for(k=0; k<n/2; k++) {
			ve_r[k] = X_serial_r[2*k]; ve_i[k] = X_serial_i[2*k];
			vo_r[k] = X_serial_r[2*k+1]; vo_i[k] = X_serial_i[2*k+1];
		}
		cpuKernel(ve_r, ve_i, n/2, X_serial_r, X_serial_i);	// FFT on even-indexed elements of v[]
		cpuKernel(vo_r, vo_i, n/2, X_serial_r, X_serial_i);	// FFT on odd-indexed elements of v[]
		
		for(m=0; m<n/2; m++) {
			w_r =  cos((2*PI*m)/n);
			w_i = -sin((2*PI*m)/n);
			z_r = w_r*vo_r[m] - w_i*vo_i[m];	// Re(w*vo[m])
			z_i = w_r*vo_i[m] + w_i*vo_r[m];	// Im(w*vo[m])
			X_serial_r[  m  ] = ve_r[m] + z_r;
			X_serial_i[  m  ] = ve_i[m] + z_i;
			X_serial_r[m+n/2] = ve_r[m] - z_r;
			X_serial_i[m+n/2] = ve_i[m] - z_i;
		}
	}
	return;
}
//-----------------------------------------------------------------------------
void get_inputs(int argc, char *argv[], unsigned int& N, unsigned int& M)
{
    if (
	argc != 2 || 
	atoi(argv[1]) < 0 || atoi(argv[1]) > 26 
	) {
        printf("<< Error >>\n");
        printf("Enter the following command:\n");
        printf("\t./a.out  M\n");
        printf("\t\tM must be between 0 and 26\n");
		exit(-1);
    }
	M = atoi(argv[1]);
    N = (1 << M);
}
//-----------------------------------------------------------------------------
void fill(float* data, int size) {
    for (int i = 0; i < size; i++)
        data[i] = (float)(rand() % 17 - 8);
}
double calc_mse(float* data1_r, float* data1_i, float* data2_r, float* data2_i, int size) {
    double mse = 0.0;
    int i;
    for (i = 0; i < size; i++) {
        double e_r = data1_r[i] - data2_r[i];
        double e_i = data1_i[i] - data2_i[i];
        double e = e_r * e_r + e_i * e_i;
        mse += e;
    }
    return mse/size;
}

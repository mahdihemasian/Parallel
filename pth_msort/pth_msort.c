// Include your C header files here

#include "pth_msort.h"
#include <stdlib.h>
#include <semaphore.h>
//pthread_mutex_t* mutex;



/*global variable*/
sem_t barrier_sem;
sem_t count_sem;
int  count=0;
long thread_count;
int* x_val;
int* x_sort;
int N_val;
int left[7];
int right[7];







void  swap(int* arr1,int* arr2){
	int temp;
	temp=*arr1;
	*arr1=*arr2;
	*arr2=temp;
}
int partition (int arr[], int low, int high)
{
    int pivot = arr[high]; // pivot
    int i = (low - 1),j; // Index of smaller element and indicates the right position of pivot found so far
 
    for (j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (arr[j] < pivot)
        {
            i++; // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
void quickSort(int arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
        at right place */
        int pi = partition(arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void merge (int *x,int* temp,int first,int middle,int last){
	int i,j,k;
	int n1=middle-first+1;
	int n2=last-middle;
	int* L=x+first;
        int* R=x+middle+1;	
    
    i = 0; 
    j = 0; 
    k = first; 
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            temp[k] = L[i];
            i++;
        }
        else {
            temp[k] = R[j];
            j++;
        }
        k++;
    }
	
    while (i < n1) {
        temp[k] = L[i];
        i++;
        k++;
    }
  
   
    while (j < n2) {
        temp[k] = R[j];
        j++;
        k++;
    }
	
	int z;
	for(z=first;z<=last;z++)
		x[z]=temp[z];
	
		
}
void merge1 (int *x,int* temp,int first,int middle,int last){
	int i,j,k;
	int n1=middle-first+1;
	int n2=last-middle;
	int* L=x+first;
        int* R=x+middle+1;	
    
    i = 0; 
    j = 0; 
    k = first; 
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            temp[k] = L[i];
            i++;
        }
        else {
            temp[k] = R[j];
            j++;
        }
        k++;
    }
	
    while (i < n1) {
        temp[k] = L[i];
        i++;
        k++;
    }
  
   
    while (j < n2) {
        temp[k] = R[j];
        j++;
        k++;
    }	
}
	
void mergesort(int* x,int* temp,int i,int j){
	if(j>i) {
	int k;
	k=i+(j-i)/2;
	mergesort(x,temp,i,k);
	mergesort(x,temp,k+1,j);
	merge(x,temp,i,k,j);
	}
}

int binarySearch(int arr[], int l, int r, int x){
    if (r >= l) {
        int mid = l + (r - l) / 2;

        if ((x <= arr[mid]) && (x >= arr[mid-1]) ) 
            return mid;

   
        if (arr[mid] > x)
            return binarySearch(arr, l, mid - 1, x);

       
        return binarySearch(arr, mid + 1, r, x);
    }


    return -1;
}
void index_finder(int*arr,int *a,int* b,int Length){
	int idx1,idx2,idx3,i;
	int x=Length/4;
	idx1=binarySearch(a,0,Length-1,b[x]);
	idx2=binarySearch(a,0,Length-1,b[x*2]);
	idx3=binarySearch(a,0,Length-1,b[3*x]);
	
	int temp[7];
	arr[0]=0;
	arr[1]=x;
	arr[2]=2*x;
	arr[3]=3*x;
	arr[4]=idx1>=0?idx1:Length;
	arr[5]=idx2>=0?idx2:Length;
	arr[6]=idx3>=0?idx3:Length;
	
	mergesort(arr,temp,0,6);

}
void copy(int* a,int* b,int length){
	
	int i=0;
	int j=0;
	int* arr1 = a;
        int* arr2 = &a[length/2];
	for(i=left[0];i<left[1];i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[0];i<right[1];i++)
	{
		b[j]=arr2[i];
		j++;
	}
	for(i=left[1];i<left[2];i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[1];i<right[2];i++)
	{
		b[j]=arr2[i];
		j++;
	}
	for(i=left[2];i<left[3];i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[2];i<right[3];i++)
	{
		b[j]=arr2[i];
		j++;
	}	
	for(i=left[3];i<left[4];i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[3];i<right[4];i++)
	{
		b[j]=arr2[i];
		j++;
	}
	for(i=left[4];i<left[5];i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[4];i<right[5];i++)
	{
		b[j]=arr2[i];
		j++;
	}
	for(i=left[5];i<left[6];i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[5];i<right[6];i++)
	{
		b[j]=arr2[i];
		j++;
	}
	for(i=left[6];i<length/2;i++)
	{
		b[j]=arr1[i];
		j++;
	}
	for(i=right[6];i<length/2;i++)
	{
		b[j]=arr2[i];
		j++;
	}	
}

//****************************************************

void* mergelvl1(void* Rank){
	long R=(long )Rank;
	mergesort(x_val,x_sort,R*(N_val/thread_count),(R+1)*(N_val/thread_count)-1);
	//quickSort(x_val,R*(N_val/thread_count),(R+1)*(N_val/thread_count)-1);

}

void* mergelvl2(void* Rank){
	long R=(long )Rank;
	//merge1(x_val,x_sort,R*N_val/2,(int)(((2*R+1)*(N_val/2)-1)/2),(R+1)*(N_val/2)-1);
	merge1(x_sort,x_val,R*N_val/2,(int)(((2*R+1)*(N_val/2)-1)/2),(R+1)*(N_val/2)-1);
}

void* mergelvl3(void* Rank){
	long R=(long )Rank;
	switch (R){
		case 0:
			merge1(x_sort,x_val,left[0]+right[0],left[1]+right[0]-1,left[1]+right[1]-1);
			merge1(x_sort,x_val,left[1]+right[1],left[2]+right[1]-1,left[2]+right[2]-1);	
			break;              
		case 1:              
			merge1(x_sort,x_val,left[2]+right[2],left[3]+right[2]-1,left[3]+right[3]-1);
			merge1(x_sort,x_val,left[3]+right[3],left[4]+right[3]-1,left[4]+right[4]-1);
			                    
			break;              
		case 2:              
			merge1(x_sort,x_val,left[4]+right[4],left[5]+right[4]-1,left[5]+right[5]-1);
			merge1(x_sort,x_val,left[5]+right[5],left[6]+right[5]-1,left[6]+right[6]-1);
			break;              
		case 3:	             
			merge1(x_sort,x_val,left[6]+right[6],right[6]+N_val/2-1,N_val-1);
			break;	
	}
	    
}



void mergeSortParallel (const int* values, unsigned int N, int* sorted) {
	int i=0;
	long 	thread;
	pthread_t* thread_handles;
	x_val= (int*)values;
	x_sort=sorted;
	N_val=N;
	thread_count=4;
	
	thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t));
	//*********** first level of MERGE*****************
	for(thread=0;thread<thread_count;thread++){
	pthread_create(&thread_handles[thread], (pthread_attr_t*) NULL,
	mergelvl1, (void*) thread);
	}
	for(thread = 0; thread < thread_count; thread++) {
      pthread_join(thread_handles[thread], NULL);
   	 }			
	//*********** second level of MERGE*****************
	thread_count=2;
	for(thread=0;thread<thread_count;thread++){
	pthread_create(&thread_handles[thread], (pthread_attr_t*) NULL,
	mergelvl2, (void*) thread);
	}
	for(thread = 0; thread < thread_count; thread++) {
      pthread_join(thread_handles[thread], NULL);
   	 }
	
	//*********** third level of MERGE*****************
	thread_count=4;
    index_finder(left,x_val,x_val+N/2,N/2);
    index_finder(right,x_val+N/2,x_val,N/2);
    copy(x_val,x_sort,N);
    
    
	for(thread=0;thread<thread_count;thread++){
	pthread_create(&thread_handles[thread], (pthread_attr_t*) NULL,
	mergelvl3, (void*) thread);
	}
	for(thread = 0; thread < thread_count; thread++) {
      pthread_join(thread_handles[thread], NULL);
   	 }
	
 	 for(i=0;i<N;i++){
   	 	sorted[i]=x_val[i];
		}
		
}

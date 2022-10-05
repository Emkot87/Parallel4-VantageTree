extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

#include "utilFunctions.h"

struct timeval startwtime, endwtime;
double seq_time;

int LIMIT = 80000 ;
int LIMIT2 = 800000 ;
int THREADS = 128 ;
int BSIZE = 4 ;


__global__ void distCalc(float* points, int* localIndex,int d, int threads, float* distance,int numPoints, int bSize){
      int index = threadIdx.x + blockIdx.x * threads;
      float dist;
      for(int j = index*bSize ; j <(index+1)*bSize ; j++){
        distance[j] = 0;
        
        for(int i = 0 ; i < d ; i++){
            dist = points[localIndex[numPoints-1]*d + i] - points[localIndex[j]*d + i]  ;
            distance[j] += dist*dist;
        }

        
        
    }
      
}

void calculateDistancesCuda(float* cValues, int* localIndex,int numPoints, int dimensPoints, float* distances){
    
    int threadsBlock = THREADS;
    float* cDistances;
    int* clocalIndex;

    // allocate space 
    cudaMallocManaged(&cDistances, numPoints*sizeof(float));
    cudaMallocManaged(&clocalIndex, numPoints*sizeof(int));

    // copy data
    cudaMemcpy(clocalIndex, localIndex, numPoints*sizeof(int),cudaMemcpyHostToDevice);
    
// find out how many blocks and threads we need given how many threads per block and how many points a thread calculates
    int s = ceil(numPoints/(float)BSIZE);
    int blocks = ceil(s/(float)threadsBlock);


    int threads =  ceil(s/(float)blocks);

    distCalc<<<blocks,threads>>>(cValues, clocalIndex, dimensPoints,threads,cDistances,numPoints,BSIZE);
    

    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n",
               cudaGetErrorString(cudaerr));
    
    // send back the distances
    cudaMemcpy(distances,cDistances, numPoints*sizeof(float),cudaMemcpyDeviceToHost);

    //and free the space
    cudaFree(cDistances);
    cudaFree(clocalIndex);
}

void calculateDistances(float* pointOfInterest, float* points, int* localIndex, int numPoints, int dimensPoints, float* distances){

    if((numPoints)*dimensPoints>LIMIT){
        #pragma omp parallel for
        for(int i = 0 ; i < numPoints ; i++){
            distances[i] = 0;
            for(int j = 0 ; j < dimensPoints ; j++){
                distances[i] += pow(*(points +(*(localIndex+i))*dimensPoints + j) - *(pointOfInterest + j),2); 
            }
        }
    }

    else{
        for(int i = 0 ; i < numPoints ; i++){
            distances[i] = 0;
            for(int j = 0 ; j < dimensPoints ; j++){
                distances[i] += pow(*(points +(*(localIndex+i))*dimensPoints + j) - *(pointOfInterest + j),2); 
            }
        }
    }
}


// splits the points by picking the last point and calculating its distance to the rest 
// uses quickselect to find the median
VtreePoint * makeTree(float *values, int N, int d, int *localIndex, VtreePoint *TlocalPar, float* cValues){
    
    if(N == 0)
        return NULL; 

    int index = localIndex[N-1];
    float * Vpointer = (values + d*index);

    int halfSize = N/2 + 1;

    float *distances = (float *)malloc(N * sizeof(float));

	int *innerPoints = (int*)malloc((halfSize) * sizeof(int));
	int *outterPoints = (int*)malloc((halfSize) * sizeof(int));

    if((N*d) > LIMIT2){
        calculateDistancesCuda(cValues,localIndex, N, d, distances);
    }
    
    else{
        calculateDistances(Vpointer, values, localIndex, N-1, d, distances);
    }

    float median = findMedian(distances, N-1);

    int count_inner = 0;
    int count_outter = 0;

    for(int i = 0 ; i < (N-1) ; i++ ){
        if( distances[i] >= median){
            outterPoints[count_outter] = localIndex[i];
            count_outter++;
        }
        else{
            innerPoints[count_inner] = localIndex[i];
            count_inner++;
        }
    }

    free(distances);

    VtreePoint *node = (VtreePoint*)malloc(sizeof(VtreePoint));

    node->parent = TlocalPar;
    node->median = median;
    node->index = index;

    if( N*d >= LIMIT){
        float *c = cValues;
        float *c2 = cValues;
        //parallel section
	  	#pragma omp parallel 
		{
		   //2 sections with thread creation for inner and outer function call
		    #pragma omp sections nowait
		    {
                #pragma omp section
                node->innerPoint = makeTree(values, count_inner, d,innerPoints,node, c);
                #pragma omp section
                node->outterPoint = makeTree(values, count_outter, d,outterPoints,node, c2);

            }
        }
    }
    else{
        node->innerPoint = makeTree(values, count_inner, d,innerPoints,node, cValues);
        node->outterPoint = makeTree(values, count_outter, d,outterPoints,node, cValues);
    }

    free(innerPoints);
    free(outterPoints);

	return node;

}



int main(int argc, char* argv[]){

    if(argc < 4){
        printf("Need 3 arguments, number of points, dimensions, and if to search nearest neighbours (0/1)\n");
        return;
    }
	// it reads the arguments 
	int N, d, searchNeis;
    N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);
    searchNeis = (int) strtol(argv[3],NULL,10);

    // read the limits if they were given
    if(argc > 7){
        LIMIT = strtol(argv[4],NULL,10);
        LIMIT2 = strtol(argv[5],NULL,10);
        THREADS = strtol(argv[6],NULL,10);
        BSIZE = strtol(argv[7],NULL,10);
    }
    printf("running ./vtCudaOmp.out with %d %d %d %d %d %d %d\n",N,d,searchNeis,LIMIT,LIMIT2,THREADS,BSIZE);

    omp_set_num_threads(32);

    // allocate space for the points on the gpu
    float* values = (float*)malloc(sizeof(float)*N*d);
    float* cValues;

    cudaMallocManaged(&cValues, N*d*sizeof(float));

    // read the points from the file vtSerial.out created
    FILE *fp;
    fp = fopen("points.bin","rb");

    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    fread(values, sizeof(float), N*d, fp);

    fclose(fp);
    // copy the points to the gpu
    cudaMemcpy(cValues, values, N*d*sizeof(float), cudaMemcpyHostToDevice);

    // create the default indexes array
    int* idxs = (int*)malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++)
        idxs[i] = i;


    gettimeofday (&startwtime, NULL);

    // and finally create the tree
    VtreePoint* root = makeTree(values, N, d, idxs,NULL,cValues);

    gettimeofday (&endwtime, NULL);
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("\n\n-=-=-=-=-=-=-=-+++total time to create tree was %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
    
    int* arrayOfIndex = (int*)malloc(N*sizeof(int));
    int i = 0;
    treeToArray(root,arrayOfIndex,&i);

    int* arraySerial = (int*)malloc(N*sizeof(int));


    FILE *fp2;
    fp2 = fopen("treeCheck.bin","rb");
    fread(arraySerial,sizeof(int),N,fp2);
    fclose(fp2);
    
    int check = 0;
    for(int k = 0 ; k < N ; k++){
        if(arrayOfIndex[k] != arraySerial[k])
            check++;
    }

    // check back if the results are correct
    if(check > 1000)
        printf("%d indexes were wrong\n",check);
    else
        printf("everything is correct\n");
    
    gettimeofday(&startwtime, NULL);

    // search up to searchNeis nearest neighbours of every point and save them in a file
    if(searchNeis){
    
        FILE *nearestNeis;

        nearestNeis = fopen("nearestNeis.txt","w");

        gettimeofday(&startwtime, NULL);

        for(int j = 0 ; j < N ; j++){
            for(int i = 2 ; i <= searchNeis ; i = i*2 ){
                struct neighbours * head = NULL;
                float thress = __FLT_MAX__;
                knnSearch(root, values, j, d, &head, &thress, i);
                print(head,values,d,nearestNeis);
                fprintf(nearestNeis,"\n");
                free_queue(&head);
            }
        }
        

        gettimeofday (&endwtime, NULL);
        seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
        printf("\n\n-=-=-=-=-=-=-=-+++total time to search everyones knn %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
        fclose(nearestNeis);
    }


    cudaFree(cValues);
    
    free(arrayOfIndex);
    free(arraySerial);
    free(root);
}
}
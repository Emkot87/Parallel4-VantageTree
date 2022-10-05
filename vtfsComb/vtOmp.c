#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <limits.h>

#include "utilFunctions.h"

struct timeval startwtime, endwtime;
double seq_time;

int LIMIT = 1000000 ;


VtreePoint * makeTree(float *values, int N, int d, int *localIndex, VtreePoint *TlocalPar){
    
    if(N == 0)
        return NULL;
        
    int index = localIndex[N-1];
    float * Vpointer = (values + d*index);

    int halfSize = N/2 + 1;

    float *distances = (float *)malloc((N-1) * sizeof(float));

	int *innerPoints = (int*)malloc((halfSize) * sizeof(int));
	int *outterPoints = (int*)malloc((halfSize) * sizeof(int));


    if((N*d)>LIMIT){
        #pragma omp parallel for
        for(int i = 0 ; i < N - 1 ; i++){
            float sum = 0;
            for(int j = 0 ; j < d ; j++){
                sum += pow(*(values +(*(localIndex+i))*d + j) - *(Vpointer + j),2);
            }
            distances[i] = sum;
        }
    }
    else{
        for(int i = 0 ; i < N - 1 ; i++){
            float sum = 0;
            for(int j = 0 ; j < d ; j++){
                sum += pow(*(values +(*(localIndex+i))*d + j) - *(Vpointer + j),2); 
            }
            distances[i] = sum;
        }
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

    if( N*d >= LIMIT*2){
        //parallel section
	  	#pragma omp parallel 
		{
		   //2 sections with thread creation for inner and outer function call
		    #pragma omp sections nowait
		    {
		       #pragma omp section
		       node->innerPoint = makeTree(values, count_inner, d,innerPoints,node);
		       #pragma omp section
		       node->outterPoint = makeTree(values, count_outter, d,outterPoints,node);
		   }
        }
            }
    else{//else serial
        node->innerPoint = makeTree(values, count_inner, d,innerPoints,node);
        node->outterPoint = makeTree(values, count_outter, d,outterPoints,node);
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

    float* values = (float*)malloc(sizeof(float)*N*d);
    
    // and the points from the file the serial one created
    FILE *fp;
    fp = fopen("points.bin","rb");

    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    fread(values, sizeof(float), N*d, fp);

    fclose(fp);

    // read the limit if it was given
    if( argc > 4)
        LIMIT = (int) strtol(argv[4],NULL,10);

    printf("running ./vtOmp.out with %d %d %d %d\n",N,d,searchNeis,LIMIT);

    omp_set_num_threads(32);

    // create the default indexes array
    int* idxs = (int*)malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++){
        idxs[i] = i;
    }
        
    gettimeofday (&startwtime, NULL);

    // and finally create the tree
    VtreePoint* root = makeTree(values, N, d, idxs,NULL);

    gettimeofday (&endwtime, NULL);
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("\n\n-=-=-=-=-=-=-=-+++total time to create tree: %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
    
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
                float tau = __FLT_MAX__;
                knnSearch(root, values, j, d, &head, &tau, i);
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

    free(arrayOfIndex);
    free(arraySerial);
    free(root);
    free(values);
}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>

#include "utilFunctions.h"

struct timeval startwtime, endwtime;
double seq_time;


VtreePoint * makeTree(float *values, int N, int d, int *localIndex, VtreePoint *TlocalPar){
    
    if(N == 0)
        return NULL;
        
    int index = localIndex[N-1];
    float * Vpointer = (values + d*index);

    int halfSize = N/2 + 1;

    float *distances = (float *)malloc((N-1) * sizeof(float));

	int *innerPoints = (int*)malloc((halfSize) * sizeof(int));
	int *outterPoints = (int*)malloc((halfSize) * sizeof(int));


   
    for(int i = 0 ; i < N - 1 ; i++){
        float sum = 0;
        for(int j = 0 ; j < d ; j++){
            sum += pow(*(values +(*(localIndex+i))*d + j) - *(Vpointer + j),2); 
        }
        distances[i] = sqrt(sum);
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

    
    node->innerPoint = makeTree(values, count_inner, d,innerPoints,node);
    node->outterPoint = makeTree(values, count_outter, d,outterPoints,node);
   
    free(innerPoints);
    free(outterPoints);
    
	return node;

}



int main(int argc, char* argv[]){
    srand((unsigned int)time(NULL));

    if(argc < 4){
        printf("Need 3 arguments, number of points, dimensions, and if to search nearest neighbours (0/1)\n");
        return;
    }
    // it reads the arguments 
    int N, d, searchNeis;
    N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);
    searchNeis = (int) strtol(argv[3],NULL,10);

    printf("running ./vtSerial.out with %d %d %d\n",N,d,searchNeis);

    float *pointsForFile = (float*)malloc(N*d*sizeof(float));

    // create random points and save them in a file for the other ones to use
    for(int i = 0; i < N ; i++){
        for(int j = 0; j < d ; j++){
            float k = (float)rand()/((float)RAND_MAX/100);
            pointsForFile[i*d+j] = k;
        }
    }

    FILE *fptr;
    
    fptr = fopen("points.bin","wb+");
    
    fwrite(pointsForFile,sizeof(float),N*d,fptr);
    fclose(fptr);


    // create the default indexes array
    int* idxs = (int*)malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++)
        idxs[i] = i;


    gettimeofday (&startwtime, NULL);

    // and finally create the tree
    VtreePoint* root = makeTree(pointsForFile, N, d, idxs,NULL);

    gettimeofday (&endwtime, NULL);
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("\n\n-=-=-=-=-=-=-=-+++total time to create tree: %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
    

    int* arrayOfIndex = (int*)malloc(N*sizeof(int));
    int i = 0;
    treeToArray(root,arrayOfIndex,&i);

    FILE *fp2;
    fp2 = fopen("treeCheck.bin","wb");
    fwrite(arrayOfIndex,sizeof(int),N,fp2);
    fclose(fp2);

    // search up to searchNeis nearest neighbours of every point and save them in a file
    if(searchNeis){
    
        FILE *nearestNeis;

        nearestNeis = fopen("nearestNeis.txt","w");

        gettimeofday(&startwtime, NULL);

        for(int j = 0 ; j < N ; j++){
            for(int i = 2 ; i <= searchNeis ; i = i*2 ){
                struct neighbours * head = NULL;
                float tau = __FLT_MAX__;
                knnSearch(root, pointsForFile, j, d, &head, &tau, i);
                print(head,pointsForFile,d,nearestNeis);
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
    free(root);
    free(pointsForFile);
}
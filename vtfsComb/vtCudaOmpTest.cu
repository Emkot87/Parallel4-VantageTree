#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>
#include <omp.h>


struct timeval startwtime, endwtime;
double seq_time;

int LIMIT = 80000 ;
int LIMIT2 = 800000 ;
int THREADS = 128 ;
int BSIZE = 1 ;

typedef struct Vtree {

    struct Vtree *innerPoint;
    struct Vtree *outterPoint;
    struct Vtree *parent;
    float *values, median;
    int index;
    
    //auto den kserw giati to exei alla parto
    float *ptr_vp;
} VtreePoint;

struct neighbours {
    //float *point;
    int index;
    float distance;
    struct neighbours* nextNei;
};


void popNei(struct neighbours** head)
{
    struct neighbours* previous = NULL;
    struct neighbours* temp = *head;


    while(temp->nextNei != NULL){

        previous = temp;
        temp = temp->nextNei;
    }
    
    previous->nextNei = NULL;
    free(temp);
}


struct neighbours* peek(struct neighbours** head)
{
    struct neighbours* temp = *head;
    while(temp->nextNei != NULL)
    {
        temp = temp->nextNei;
    }
    return temp;
}

int isEmpty(struct neighbours** head)
{
    return (*head) == NULL;
}

struct neighbours* newNode(float d,int index)
{
    struct neighbours* temp = (struct neighbours*)malloc(sizeof(struct neighbours));
    temp->distance = d;
    temp->nextNei = NULL;
    temp->index = index;
 
    return temp;
}

int getCount(struct neighbours** head)
{   
    if(isEmpty(head))
        return 0;

      // Initialize count

    struct neighbours* current = *head;  // Initialize current

    int count = 0;

    while (current != NULL)
    {
        count++;
        //printf("to count einai e? %d",count);
        current = current->nextNei;
    }
    return count;
}

void pushNei(struct neighbours** head, float d,int index)
{
    struct neighbours* start = (*head);
 
    // Create new Node
    struct neighbours* temp = newNode(d,index); 
 
    // Special Case: The head of list has lesser
    // priority than new node. So insert new
    // node before head node and change head node.
    
    if(isEmpty(head)){

        *head = temp;
        return;

    }

    
    if ((*head)->distance > d) {
 
        // Insert New Node before head
        temp->nextNei = *head;
        (*head) = temp;
    }
    else {
 
        // Traverse the list and find a
        // position to insert new node
        while (start->nextNei != NULL && start->nextNei->distance < d) {
            start = start->nextNei;
        }
 
        // Either at the ends of the list
        // or at required position
        temp->nextNei = start->nextNei;
        start->nextNei = temp;
    }
    
}

void queue_to_arr(struct neighbours* head, int* array, int k)
{
    struct neighbours* temp = head;
    for(int i = 0; i < k; i++)
    {
        if(temp != NULL)
        {
            array[i] = temp->index;
            temp = temp->nextNei;
        }
    }
}

void free_queue(struct neighbours** head)
{
    struct neighbours* temp = (*head);
    while(temp != NULL)
    {
        struct neighbours* next = temp->nextNei;
        free(temp);
        temp = next;
    }
}

void treeToArray(VtreePoint* point, int* arrayOfIndex,int* i){
    if(point != NULL){
        arrayOfIndex[*i] = point->index;
        (*i)++;
        treeToArray(point->innerPoint, arrayOfIndex, i);
        treeToArray(point->outterPoint, arrayOfIndex, i);
    }
}


#define swap(x, y) { float temp = x; x = y; y = temp; }

int partition(float arr[], int l, int r)
{
    float x = arr[r];
	int	i = l;
    for (int j = l; j <= r - 1; j++) {
        if (arr[j] <= x) {
            swap(arr[i], arr[j]);
            i++;
        }
    }
    swap(arr[i], arr[r]);
    return i;
}

__global__ void distCalc( int* localIndex,int d, int threads, float* distance,int numPoints, int bSize){
      int index = threadIdx.x + blockIdx.x * threads;
      float dist;
      //printf("index : %d\n",index);
      for(int j = index*bSize ; j <(index+1)*bSize ; j++){
        distance[j] = 0;
        
        for(int i = 0 ; i < d ; i++){
            dist = localIndex[numPoints-1]*d + i - localIndex[j]*d + i  ;
            distance[j] += dist*dist;
        }

        
        
    }
      
}

void calculateDistancesCuda( int* localIndex,int numPoints, int dimensPoints, float* distances){
    
    int threadsBlock = THREADS;
    float* cDistances;
    int* clocalIndex;

    cudaMalloc(&cDistances, numPoints*sizeof(float));
    cudaMalloc(&clocalIndex, numPoints*sizeof(int));

    //cudaMemcpy(cPointsCut, points, numPoints*dimensPoints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(clocalIndex, localIndex, numPoints*sizeof(int),cudaMemcpyHostToDevice);
    
    int s = ceil(numPoints/(float)BSIZE);
    int blocks = ceil(s/(float)threadsBlock);

    //printf("blocks %d\n", blocks);
    int threads =  ceil(s/(float)blocks);
    //printf("threads %d\n", threads);
    //printf("ola popa ws edw ?\n");
    distCalc<<<blocks,threads>>>( clocalIndex, dimensPoints,threads,cDistances,numPoints,BSIZE);
    

    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n",
               cudaGetErrorString(cudaerr));
    

    cudaMemcpy(distances,cDistances, numPoints*sizeof(float),cudaMemcpyDeviceToHost);

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


float kthSmallest(float arr[], int l, int r, int k)
{
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= r - l + 1) {
 
        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, l, r);
 
        // If position is same as k
        if (index - l == k - 1)
            return arr[index];
 
        // If position is more, recur
        // for left subarray
        if (index - l > k - 1)
            return kthSmallest(arr, l, index - 1, k);
 
        // Else recur for right subarray
        return kthSmallest(arr, index + 1, r,
                            k - index + l - 1);
    }
 
    // If k is more than number of
    // elements in array
    return 12;
}

int RandRange(int Min, int Max)
{	
    srand((unsigned int)time(NULL));
    int diff = Max-Min;
    return (int) (((double)(diff+1)/RAND_MAX) * rand() + Min);
}

float findMedian(float *distances, int N){
    int n1, n2;
    float median;
    float* tempDist = (float*)malloc(N*sizeof(float));
    // Find median by creating a copy of distances and using it to QuickSelect
	for(int i=0;i<N;i++){
		tempDist[i] = *(distances + i);
	}
	if(N % 2 == 0) {
        n1 = (N + 2) / 2;
        n2 = N / 2;
        median = (kthSmallest(tempDist, 0, N - 1, n1) + kthSmallest(distances, 0, N - 1, n2)) / 2.0;
    }
    else {
        median = kthSmallest(tempDist, 0, N - 1, (N+1)/2 );
    }

    return median;
}

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
        calculateDistancesCuda(localIndex, N, d, distances);
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
        //parallel section
	  	#pragma omp parallel 
		{
		   //2 sections with thread creation for inner and outer function call
		    #pragma omp sections nowait
		    {
                #pragma omp section
                node->innerPoint = makeTree(values, count_inner, d,innerPoints,node, cValues);
                #pragma omp section
                node->outterPoint = makeTree(values, count_outter, d,outterPoints,node, cValues);

            }
        }
    }
    else{
        node->innerPoint = makeTree(values, count_inner, d,innerPoints,node, cValues);
        node->outterPoint = makeTree(values, count_outter, d,outterPoints,node, cValues);
    }
	return node;

}

// Oxi re mia xara leitourgei
void knnSearch(VtreePoint *TreeNode,float* points, int index, int d, struct neighbours **neiList, float* tau, int k){
    if(TreeNode == NULL)
        return;

    float distance = 0;


    for(int j = 0 ; j < d ; j++){
        //printf("points[index*d + j] : %f TreeNode->index*d + j : %d \n",points[index*d + j], TreeNode->index*d + j);
        distance += pow(points[index*d + j] - points[TreeNode->index*d + j],2);
    }


    if(distance < *tau){
        //printf("eftasa edw mia fora\n");
        if (getCount(neiList) == k){
            popNei(neiList);
        }
        pushNei(neiList,distance,TreeNode->index);

        if(getCount(neiList) == k){
            *tau = peek(neiList)->distance;
        }
    }



   if(distance < TreeNode->median){
       if(distance - *tau <= TreeNode->median){
           knnSearch(TreeNode->innerPoint,points,index,d,neiList,tau,k);
       }
       if(distance + *tau >= TreeNode->median){
           knnSearch(TreeNode->outterPoint,points,index,d,neiList,tau,k);
       }
    }
   else{
       if(distance + *tau >= TreeNode->median){
           knnSearch(TreeNode->outterPoint,points,index,d,neiList,tau,k);
       }
       if(distance - *tau <= TreeNode->median){
           knnSearch(TreeNode->innerPoint,points,index,d,neiList,tau,k);
       }
    }
    

}


void print(struct neighbours *head, float* points, int d) {
    struct neighbours *current_node = head;
   	while ( current_node != NULL) {
        for(int i = 0 ; i < d ; i++){
        printf("%f ", points[current_node->index*d + i]);
        }
        current_node = current_node->nextNei;
        printf("\n");
    }
}



int main(int argc, char* argv[]){
    int N, d;

    if(argc < 3){
		printf("2 arguments, number of points and dimension\n");
		return 0;
	}
	
	N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);

    if(argc > 6){
        LIMIT = strtol(argv[3],NULL,10);
        LIMIT2 = strtol(argv[4],NULL,10);
        THREADS = strtol(argv[5],NULL,10);
        BSIZE = strtol(argv[6],NULL,10);
    }

    //omp_set_num_threads(8);

    float* values = (float*)malloc(sizeof(float)*N*d);
    float* cValues;

    cudaMalloc(&cValues, N*d*sizeof(float));

    FILE *fp;
    fp = fopen("points.bin","rb");

    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    fread(values, sizeof(float), N*d, fp);

    fclose(fp);

    cudaMemcpy(cValues, values, N*d*sizeof(float), cudaMemcpyHostToDevice);

    
    int* idxs = (int*)malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++)
        idxs[i] = i;


    gettimeofday (&startwtime, NULL);


    VtreePoint* root = makeTree(values, N, d, idxs,NULL,cValues);

    gettimeofday (&endwtime, NULL);
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("\n\n-=-=-=-=-=-=-=-+++total time to create tree was %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
    
    int* arrayOfIndex = (int*)malloc(N*sizeof(int));
    int i = 0;
    treeToArray(root,arrayOfIndex,&i);

    for(int k = 0; k < N/10000 ; k++){
        printf("%d\n",arrayOfIndex[k]);
    }

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
    
    if(check == 0)
        printf("everything is correct\n");
    else
        printf("%d indexes were wrong\n",check);
    
    gettimeofday(&startwtime, NULL);

    for(int j = 0 ; j < N/1000 ; j++){
        for(int i = 2 ; i <= 8 ; i = i*2 ){
            struct neighbours * head = NULL;
            float tau = __FLT_MAX__;
            knnSearch(root, values, j, d, &head, &tau, i);
            //print(head,values,d);
            free_queue(&head);
        }
    }
    

    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("\n\n-=-=-=-=-=-=-=-+++total time to search everyones knn was %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
	
    cudaFree(cValues);
    free(arrayOfIndex);
    free(arraySerial);
    free(root);
}
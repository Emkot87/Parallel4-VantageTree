#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <limits.h>

struct timeval startwtime, endwtime;
double seq_time;

int LIMIT = 1000000 ;
int LIMIT2 = 2000 ;
int THREADS = 256 ;
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

__global__ void distCalc(float* point, float* points,int d, int threads, float* distance, int bSize){
      int index = threadIdx.x + blockIdx.x * threads;
      //printf("index : %d\n",index);
      for(int j = index*bSize ; j <(index+1)*bSize ; j++){
        distance[j] = 0;
        
        for(int i = 0 ; i < d ; i++){
            float dist = point[i] - points[j*d + i];
            distance[j] += dist*dist;
        }
        
    }
      
}

void calculateDistancesCuda(float* pointOfInterest, float* points, int numPoints, int dimensPoints, float* distances){
    
    int threadsBlock = THREADS;
    

    float* cDistances;
    float* cPointsCut;
    float* cPoint;

    cudaMalloc(&cDistances, numPoints*sizeof(float));
    cudaMalloc(&cPointsCut, numPoints*dimensPoints*sizeof(float));
    cudaMalloc(&cPoint, dimensPoints*sizeof(float));

    cudaMemcpy(cPointsCut, points, numPoints*dimensPoints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(cPoint, pointOfInterest, dimensPoints*sizeof(float),cudaMemcpyHostToDevice);

    int s = ceil(numPoints/(float)BSIZE);
    int blocks = ceil(s/(float)threadsBlock);

    //printf("blocks %d\n", blocks);
    int threads =  ceil(s/(float)blocks);
    //printf("threads %d\n", threads);

    distCalc<<<blocks,threads>>>(cPoint, cPointsCut, dimensPoints,threads,cDistances,BSIZE);
    
    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n",
               cudaGetErrorString(cudaerr));
    

    cudaDeviceSynchronize();

    cudaMemcpy(distances,cDistances, numPoints*sizeof(float),cudaMemcpyDeviceToHost);

}

void calculateDistances(float* pointOfInterest, float* points,int* localIndex, int numPoints, int dimensPoints, float* distances){

    //float* distances = (float*)malloc(sizeof(float)*numPoints);
    //float* pointsCut = (float*)malloc(sizeof(float)*numPoints*dimensPoints);
    //float* distancesTest = (float*)malloc(sizeof(float)*numPoints);



    if((numPoints)*dimensPoints>LIMIT){
        #pragma omp parallel for
        for(int i = 0 ; i < numPoints - 1 ; i++){
            float sum = 0;
            for(int j = 0 ; j < dimensPoints ; j++){
                sum += pow(*(points +(*(localIndex+i))*dimensPoints + j) - *(pointOfInterest + j),2); 
            }
            distances[i] = sqrt(sum);
        }
    }
    else{
        for(int i = 0 ; i < numPoints - 1 ; i++){
            float sum = 0;
            for(int j = 0 ; j < dimensPoints ; j++){
                sum += pow(*(points +(*(localIndex+i))*dimensPoints + j) - *(pointOfInterest + j),2); 
            }
            distances[i] = sqrt(sum);
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

VtreePoint * makeTree(float *values, int N, int d, int *localIndex, VtreePoint *TlocalPar){
    
    if(N == 0)
        return NULL; 

    int index = localIndex[N-1];
    float * Vpointer = (values + d*index);

    int halfSize = N/2 + 1;

    float *distances = (float *)malloc((N-1) * sizeof(float));

	int *innerPoints = (int*)malloc((halfSize) * sizeof(int));
	int *outterPoints = (int*)malloc((halfSize) * sizeof(int));
    

    if(N*d > LIMIT2){
        float *cutPoints = (float*)malloc((N-1)*d*sizeof(float));


        for(int i = 0 ; i < N-1 ; i++){
            for(int j = 0 ; j < d ; j++){
                cutPoints[i*d + j] = values[localIndex[i]*d + j];
            }
        }

        calculateDistancesCuda(Vpointer, cutPoints, N-1, d, distances);
        
    }
    else{
        calculateDistances(Vpointer, values,localIndex, N, d , distances);
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

   if(argc > 2){
    LIMIT = strtol(argv[3],NULL,10);
    LIMIT2 = strtol(argv[4],NULL,10);
    THREADS = strtol(argv[5],NULL,10);
    BSIZE = strtol(argv[6],NULL,10);
   }

    omp_set_num_threads(8);

    float* values = (float*)malloc(sizeof(float)*N*d);
    

    for(int i = 0; i < N ; i++){
        for(int j = 0; j < d ; j++){
            values[i*d+j] = (float)rand()/((float)RAND_MAX/100);
        }
    }

    
    int* idxs = (int*)malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++)
        idxs[i] = i;


    gettimeofday (&startwtime, NULL);


    VtreePoint* root = makeTree(values, N, d, idxs,NULL);

    gettimeofday (&endwtime, NULL);
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("\n\n-=-=-=-=-=-=-=-+++total time to create tree was %f+++-=-=-=-=-=-=-=-=-=-\n\n",seq_time);
    
    
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
	
}
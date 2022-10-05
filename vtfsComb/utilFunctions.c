#include "utilFunctions.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define swap(x, y) { float temp = x; x = y; y = temp; }

void print(struct neighbours *head, float* points, int d, FILE* neis) {
    struct neighbours *current_node = head;
   	while ( current_node != NULL) {
        for(int i = 0 ; i < d ; i++){
        fprintf(neis, "%f ", points[current_node->index*d + i]);
        }
        current_node = current_node->nextNei;
        fprintf(neis,"\n");
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

void knnSearch(VtreePoint *TreeNode,float* points, int index, int d, struct neighbours **neiList, float* tau, int k){
    if(TreeNode == NULL)
        return;

    float distance = 0;


    for(int j = 0 ; j < d ; j++){
        
        distance += pow(points[index*d + j] - points[TreeNode->index*d + j],2);
    }


    if(distance < *tau){

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
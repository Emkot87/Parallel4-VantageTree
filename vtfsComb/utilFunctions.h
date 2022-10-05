#ifndef UTILFUNTIONS_H_
#define UTILFUNTIONS_H_

#include <stdio.h>
#include <limits.h>

// struct for the nodes of the tree
typedef struct Vtree {

    struct Vtree *innerPoint;
    struct Vtree *outterPoint;
    struct Vtree *parent;
    float *values, median;
    int index;
    float *ptr_vp;

} VtreePoint;

// struct for the priority queue
struct neighbours {

    int index;
    float distance;
    struct neighbours* nextNei;

};

// print the nearest neighbours
void print(struct neighbours *head, float* points, int d, FILE* neis);

float kthSmallest(float arr[], int l, int r, int k);

int RandRange(int Min, int Max);

// quickselect
float findMedian(float *distances, int N);

// remove one from the queue
void popNei(struct neighbours** head);

// find the last 
struct neighbours* peek(struct neighbours** head);

// check if the queue is empty
int isEmpty(struct neighbours** head);

// create and allocate memory for a new neighbour 
// basically a constructor
struct neighbours* newNode(float d,int index);

// find how many members we have
int getCount(struct neighbours** head);

// add a member to the queue
void pushNei(struct neighbours** head, float d,int index);

// turn the queue into an array
void queue_to_arr(struct neighbours* head, int* array, int k);

// free the memory for every member of the queue
void free_queue(struct neighbours** head);

// read the tree and save it to an array
void treeToArray(VtreePoint* point, int* arrayOfIndex,int* i);

// partition points for quickselect
int partition(float arr[], int l, int r);

// search  k nearest neighbours of indexth point
void knnSearch(VtreePoint *TreeNode,float* points, int index, int d, struct neighbours **neiList, float* tau, int k);

#endif
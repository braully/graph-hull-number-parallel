//#include "UndirectedSparseGraph.h"


__host__ __device__ void printCombination(int *currentCombination,
        int sizeComb);

__host__ __device__
void initialCombination(int n, int k, int* combinationArray, int idx);

void parallelFindHullNumber(UndirectedCSRGraph *graph);

void serialFindHullNumber(UndirectedCSRGraph *graph);

int checkConvexityP3(UndirectedCSRGraph *graph,
        unsigned char *aux,
        int auxSize,
        int *currentCombination,
        int sizeComb, long idx);


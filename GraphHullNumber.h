//#include "UndirectedSparseGraph.h"


void parallelFindHullNumber(UndirectedCSRGraph *graph);

void serialFindHullNumber(UndirectedCSRGraph *graph);

int checkConvexityP3(UndirectedCSRGraph *graph,
        unsigned char *aux,
        int auxSize,
        int *currentCombination,
        int sizeComb);


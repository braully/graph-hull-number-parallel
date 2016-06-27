#include "UndirectedSparseGraph.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

#define DEFAULT_THREAD_PER_BLOCK 512
#define MAX_DEFAULT_SIZE_QUEUE 50
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))



//References and Examples:
//https://msdn.microsoft.com/en-us/library/aa289166(v=vs.71).aspx
//D. Knuth. The Art of Computer Programming: Generating All Combinations and Partitions. Number v. 3-4 in Art of Computer Programming. Addison-Wesley, 2005.

__host__ __device__
int maxCombinations(int n, int k) {
    if (n == 0 || k == 0) {
        return 0;
    }
    if (n < k) {
        return 0;
    }
    if (n == k) {
        return 1;
    }
    long delta, idxMax;
    if (k < n - k) {
        delta = n - k;
        idxMax = k;
    } else {
        delta = k;
        idxMax = n - k;
    }

    long ans = delta + 1;
    for (long i = 2; i <= idxMax; ++i) {
        ans = (ans * (delta + i)) / i;
    }
    return ans;
}

__host__ __device__
int maxCombinations(int n, int k, unsigned int *cacheMaxCombination) {
    if (cacheMaxCombination[k]) {
        return cacheMaxCombination[k];
    }
    return cacheMaxCombination[k] = maxCombinations(n, k);
}

__host__ __device__
void initialCombination(int n, int k, int* combinationArray, int idx, unsigned int *cacheMaxCombination) {
    int a = n;
    int b = k;
    long x = (maxCombinations(n, k) - 1) - idx;
    for (int i = 0; i < k; ++i) {
        combinationArray[i] = a - 1;
        while (maxCombinations(combinationArray[i], b) > x) {
            --combinationArray[i];
        }
        x = x - maxCombinations(combinationArray[i], b);
        a = combinationArray[i];
        b = b - 1;
    }

    for (int i = 0; i < k; ++i) {
        combinationArray[i] = (n - 1) - combinationArray[i];
    }
}

void initialCombination(int n, int k, int* combinationArray) {
    for (int i = 0; i < k; i++) {
        combinationArray[i] = i;
    }
}

__host__ __device__
void nextCombination(int n,
        int k,
        int* currentCombination) {
    if (currentCombination[0] == n - k) {
        return;
    }
    int i;
    for (i = k - 1; i > 0 && currentCombination[i] == n - k + i; --i);
    ++currentCombination[i];
    for (int j = i; j < k - 1; ++j) {
        currentCombination[j + 1] = currentCombination[j] + 1;
    }
}

__host__ __device__ void printCombination(int *currentCombination,
        int sizeComb) {
    printf("Combination: {");
    for (int i = 0; i < sizeComb; i++) {
        printf("%d", currentCombination[i]);
        if (i < sizeComb - 1) {
            printf(", ");
        }
    }
    printf("}");
}

bool verboseSerial = false;
bool verboseParallel = false;

__host__ __device__
void printQueue(int *queue, int headQueue, int tailQueue) {
    printf("\nQueue(%d):{", tailQueue - headQueue);
    for (int i = headQueue; i <= tailQueue; i++) {
        printf("%d", queue[i]);
        if (i < tailQueue) {
            printf(", ");
        }
    }
    printf("}\n");
}

__host__ __device__
int checkConvexityP3CSR(int *csrColIdxs, int nvertices,
        int *csrRowOffset, int sizeRowOffset,
        unsigned char *aux,
        int auxSize,
        int *currentCombination,
        int sizeComb) {
    //clean aux vector            
    for (int i = 0; i < auxSize; i++) {
        aux[i] = 0;
    }
    int closeCount = 0;
    int maxSizeQueue = MAX(auxSize / 2, MAX_DEFAULT_SIZE_QUEUE);
    int *queue = (int *) malloc(maxSizeQueue * sizeof (int));
    int headQueue = 0;
    int tailQueue = -1;

    for (int i = 0; i < sizeComb; i++) {
        tailQueue = (tailQueue + 1) % maxSizeQueue;
        queue[tailQueue] = currentCombination[i];
    }

    //    if (verboseSerial) {
    //        printCombination(currentCombination, sizeComb);
    //    }

    while (headQueue <= tailQueue) {
        //                if (verboseSerial) {
        //                    printQueue(queue, headQueue, tailQueue);
        //                }
        int verti = queue[headQueue];
        headQueue = (headQueue + 1) % maxSizeQueue;
        //        if (verboseSerial) {
        //            printf("vi: %d", verti);
        //        }

        if (aux[verti] < PROCESSED && verti < nvertices) {
            closeCount++;
            int end = csrColIdxs[verti + 1];
            for (int i = csrColIdxs[verti]; i < end; i++) {
                int vertn = csrRowOffset[i];
                if (vertn != verti) {
                    unsigned char previousValue = aux[vertn];
                    if (previousValue < INCLUDED) {
                        aux[vertn] = aux[vertn] + NEIGHBOOR_COUNT_INCLUDED;
                    }
                    if (previousValue < INCLUDED && aux[vertn] >= INCLUDED) {
                        tailQueue = (tailQueue + 1) % maxSizeQueue;
                        queue[tailQueue] = vertn;
                        //                        if (verboseSerial) printf("\n\t vert-includ: %d, ", vertn);
                    }
                }
            }
            aux[verti] = PROCESSED;
        }
    }
    free(queue);
    return closeCount;
}

int checkConvexityP3(UndirectedCSRGraph *graph,
        unsigned char *aux,
        int auxSize,
        int *currentCombination,
        int sizeComb) {
    return checkConvexityP3CSR(graph->getCsrColIdxs(), graph->getVerticesCount(),
            graph->getCsrRowOffset(), graph->getSizeRowOffset(),
            aux, auxSize, currentCombination, sizeComb);

    //    int checkConvexityP3CSR(int *csrColIdxs, int nvertices,
    //        int *csrRowOffset, int sizeRowOffset,
    //        unsigned char *aux,
    //        int auxSize,
    //        int *currentCombination,
    //        int sizeComb)
}

void serialFindHullNumber(UndirectedCSRGraph *graph) {
    int nvs = graph->getVerticesCount();
    int k;
    unsigned char *aux = new unsigned char [nvs];
    int *currentCombination;

    int currentSize = 0;
    int maxSize = nvs;
    int sizeCurrentHcp3 = 0;

    bool found = false;

    graph->begin_serial_time = clock();
    while (currentSize < maxSize && !found) {
        currentSize++;
        k = currentSize;
        int maxCombination = maxCombinations(nvs, k);
        currentCombination = (int *) malloc(k * sizeof (int));
        initialCombination(nvs, k, currentCombination);
        if (verboseSerial)
            printf("\nComb(%d,%d)=%d", nvs, k, maxCombination);
        for (int i = 0; i < maxCombination && !found; i++) {
            sizeCurrentHcp3 = checkConvexityP3(graph, aux, nvs, currentCombination, k);
            //                        printCombination(currentCombination, k);
            found = (sizeCurrentHcp3 == nvs);
            if (!found)
                nextCombination(nvs, k, currentCombination);
        }
        if (found) {
            printf("Result Serial: ");
            printCombination(currentCombination, currentSize);
            printf(" |S| = %d |hcp3(S)| = |V(g)| = %d\n", k, sizeCurrentHcp3);
        }
        free(currentCombination);
    }
    graph->end_serial_time = clock();
    free(aux);
}

__global__ void kernelFindHullNumber(int *csrColIdxs, int nvertices,
        int *csrRowOffset, int sizeRowOffset, int maxCombination,
        int k, int offset, int *result, unsigned char *aux,
        unsigned int *cacheMaxCombination) {
    //    printf("\n%d-Comb(%d,%d)=%d", k, nvertices, k, maxCombination);
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //    printf("\nT(%2d)= G(%d,%d),C(%d,%d)=%d", idx, nvertices, sizeRowOffset, nvertices, k, maxCombination);
    //    int idx = 1;
    bool found = false;
    int *currentCombination = (int *) malloc(k * sizeof (int));
    int auxoffset = idx*nvertices;
    //     = (unsigned char *) malloc(nvertices * sizeof (unsigned char));
    int sizeCurrentHcp3 = 0;
    int limmit = (idx + 1) * offset;
    if (limmit > maxCombination) {
        limmit = maxCombination;
    }
    int i = idx * offset;

    int maxSizeQueue = MAX(nvertices / 2, MAX_DEFAULT_SIZE_QUEUE);
    int *queue = (int *) malloc(maxSizeQueue * sizeof (int));

    //    bool debugar = false;

    //    printf("\nT(%2d:|%d|-|%d|)", idx, i, limmit);

    //    printf("\nN. vertices: %d", nvertices);
    //    printf("\nSize RowOffset: %d", sizeRowOffset);
    //    printf("\nR: {");
    //    for (int i = 0; i < nvertices + 1; i++) {
    //        printf("%d", csrColIdxs[i]);
    //        if (i < nvertices) {
    //            printf(", ");
    //        }
    //    }
    //    printf("} \nC: {");
    //    for (int i = 0; i < sizeRowOffset; i++) {
    //        printf("%d", csrRowOffset[i]);
    //        if (i < sizeRowOffset - 1) {
    //            printf(", ");
    //        }
    //    }
    //    printf("}");

    initialCombination(nvertices, k, currentCombination, i, cacheMaxCombination);
    //    printf("\nT(%2d:|%d|-|%d|): ", idx, i, limmit);
    //    printCombination(currentCombination, k);
    while (i < limmit && !found && !result[0]) {
        //        debugar = k == 14 && currentCombination[0] == 0
        //                && currentCombination[1] == 7 && currentCombination[2] == 8
        //                && currentCombination[3] == 10 && currentCombination[4] == 12
        //                && currentCombination[5] == 13 && currentCombination[6] == 15
        //                && currentCombination[7] == 17 && currentCombination[8] == 19
        //                && currentCombination[9] == 20 && currentCombination[10] == 21
        //                && currentCombination[11] == 22 && currentCombination[12] == 23
        //                && currentCombination[13] == 24
        //        ;
        //        sizeCurrentHcp3 = checkConvexityP3CSR(csrColIdxs, nvertices,
        //                csrRowOffset, sizeRowOffset,
        //                aux, nvertices, currentCombination, k);

        int headQueue = 0;
        int tailQueue = -1;
        sizeCurrentHcp3 = 0;

        for (int y = 0; y < nvertices; y++) {
            aux[auxoffset + y] = 0;
        }

        for (int j = 0; j < k; j++) {
            tailQueue = (tailQueue + 1) % maxSizeQueue;
            queue[tailQueue] = currentCombination[j];
        }
        //        if (debugar) {
        //            printCombination(currentCombination, k);
        //            printf("\n");
        //        }
        while (headQueue <= tailQueue) {
            //            if (debugar) {
            //                printQueue(queue, headQueue, tailQueue);
            //            }

            int verti = queue[headQueue];
            headQueue = (headQueue + 1) % maxSizeQueue;

            //            if (debugar) {
            //                printf("vi: %d", verti);
            //            }

            if (aux[auxoffset + verti] < PROCESSED) {
                sizeCurrentHcp3++;
                int end = csrColIdxs[verti + 1];
                int x = csrColIdxs[verti];
                //                if (debugar) {
                //                    printf("-- processando d(vi[%d-%d])=%d", x, end, end - csrColIdxs[verti]);
                //                }
                for (; x < end; x++) {
                    int vertn = csrRowOffset[x];
                    if (vertn != verti) {
                        unsigned char previousValue = aux[auxoffset + vertn];
                        if (previousValue < INCLUDED) {
                            aux[auxoffset + vertn] = aux[auxoffset + vertn] + NEIGHBOOR_COUNT_INCLUDED;
                        }
                        if (previousValue < INCLUDED && aux[auxoffset + vertn] >= INCLUDED) {
                            tailQueue = (tailQueue + 1) % maxSizeQueue;
                            queue[tailQueue] = vertn;
                        }
                    }
                }
                aux[auxoffset + verti] = PROCESSED;
            }
            //            else {
            //                printf("vi: %d -- ja processado", verti);
            //            }
        }
        //        if (debugar) {
        //            printf("\nHcp3(idx-%d,il-%d,k-%d)=%d ", idx, i, k, sizeCurrentHcp3);
        //        }
        found = (sizeCurrentHcp3 == nvertices);
        if (!found) {
            nextCombination(nvertices, k, currentCombination);
            i++;
        }
    }
    if (found) {
        result[0] = sizeCurrentHcp3;
        result[1] = (i - 1);
        //        printf("\n\nFind Parallel\n");
    }
    free(queue);
    free(currentCombination);
    //    free(aux);
}

void parallelFindHullNumber(UndirectedCSRGraph *graph) {
    int nvs = graph->getVerticesCount();
    int k;
    int currentSize = 0;
    int maxSize = nvs;
    int result[2];
    result[0] = result[1] = 0;
    int* csrColIdxs = graph->getCsrColIdxs();
    int verticesCount = graph->getVerticesCount();
    int* csrRowOffset = graph->getCsrRowOffset();
    int sizeRowOffset = graph->getSizeRowOffset();

    int* csrColIdxsGpu;
    int* csrRowOffsetGpu;
    int *resultGpu;
    unsigned char *auxGpu;
    unsigned int *cacheMaxCombination;

    //    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 50 * 1024 * 1024);

    int numBytesClsIdx = sizeof (int)*(verticesCount + 1);
    cudaMalloc((void**) &csrColIdxsGpu, numBytesClsIdx);

    int numBytesRowOff = sizeof (int)*sizeRowOffset;
    cudaMalloc((void**) &csrRowOffsetGpu, numBytesRowOff);

    cudaMalloc((void**) &auxGpu, sizeof (unsigned char)*DEFAULT_THREAD_PER_BLOCK * nvs);
    cudaMalloc((void**) &cacheMaxCombination, sizeof (unsigned int)*nvs);
    cudaMemset(cacheMaxCombination, 0, sizeof (unsigned int)*nvs);

    int numBytesResult = sizeof (int)*2;
    cudaMalloc((void**) &resultGpu, numBytesResult);

    if (resultGpu == NULL || csrRowOffsetGpu == NULL || csrColIdxsGpu == NULL) {
        perror("Failed allocate memory in GPU");
    }

    cudaError_t r = cudaMemcpy(csrColIdxsGpu, csrColIdxs, numBytesClsIdx, cudaMemcpyHostToDevice);
    if (r != cudaSuccess) {
        perror("Failed to copy memory");
    }
    r = cudaMemcpy(csrRowOffsetGpu, csrRowOffset, numBytesRowOff, cudaMemcpyHostToDevice);
    if (r != cudaSuccess) {
        perror("Failed to copy memory");
    }
    r = cudaMemcpy(resultGpu, result, numBytesResult, cudaMemcpyHostToDevice);
    if (r != cudaSuccess) {
        perror("Failed to copy memory");
    }

    bool found = false;

    graph->begin_parallel_time = clock();
    while (currentSize < maxSize - 1 && !found) {
        //    while (currentSize < 14 && !found) {
        currentSize++;
        k = currentSize;
        int maxCombination = maxCombinations(nvs, k);

        int threadsPerBlock = DEFAULT_THREAD_PER_BLOCK;
        if (DEFAULT_THREAD_PER_BLOCK > maxCombination) {
            threadsPerBlock = maxCombination / 3;
        }
        //        int threadsPerBlock = 1;

        int offset = maxCombination / threadsPerBlock;
        //        int offset = maxCombination;
        //        clock_t t = clock();
        kernelFindHullNumber << < 1, threadsPerBlock >>> (csrColIdxsGpu, verticesCount,
                csrRowOffsetGpu, sizeRowOffset, maxCombination, k, offset, resultGpu, auxGpu,
                cacheMaxCombination);
        //        t = clock() - t;
        //        printf("\n%d-Comb(%d,%d)=%d --> tempo: %d", k, verticesCount, k, maxCombination, t);
        cudaMemcpy(result, resultGpu, numBytesResult, cudaMemcpyDeviceToHost);
        found = (result[0] == nvs);
    }
    if (found) {
        printf("Result Parallel: S=%d-Comb(%d,%d) |S| = %d |Hcp3(S)| = |V(g)| = %d\n", result[1], nvs, k, k, result[0]);
    }
    //    else {
    //        printf("Result Parallel: Not found!");
    //    }
    graph->end_parallel_time = clock();
    cudaFree(resultGpu);
    cudaFree(csrRowOffsetGpu);
    cudaFree(csrColIdxsGpu);
    cudaFree(auxGpu);
    cudaFree(cacheMaxCombination);
}

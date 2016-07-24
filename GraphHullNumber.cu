#include "UndirectedSparseGraph.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

#define DEFAULT_THREAD_PER_BLOCK 512
#define MAX_DEFAULT_SIZE_QUEUE 200
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))



//References and Examples:
//https://msdn.microsoft.com/en-us/library/aa289166(v=vs.71).aspx
//D. Knuth. The Art of Computer Programming: Generating All Combinations and Partitions. Number v. 3-4 in Art of Computer Programming. Addison-Wesley, 2005.


#define verboseSerial false
#define verboseParallel false

__host__ __device__
long maxCombinations(long n, long k) {
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
long maxCombinations(long n, long k, unsigned long *cacheMaxCombination) {
    if (cacheMaxCombination[k]) {
        return cacheMaxCombination[k];
    }
    return cacheMaxCombination[k] = maxCombinations(n, k);
}

__host__ __device__
void initialCombination(long n, long k, long* combinationArray, long idx, unsigned long *cacheMaxCombination) {
    long a = n;
    long b = k;
    long x = (maxCombinations(n, k) - 1) - idx;
    for (long i = 0; i < k; ++i) {
        combinationArray[i] = a - 1;
        while (maxCombinations(combinationArray[i], b) > x) {
            --combinationArray[i];
        }
        x = x - maxCombinations(combinationArray[i], b);
        a = combinationArray[i];
        b = b - 1;
    }

    for (long i = 0; i < k; ++i) {
        combinationArray[i] = (n - 1) - combinationArray[i];
    }
}

void initialCombination(long n, long k, long* combinationArray) {
    for (long i = 0; i < k; i++) {
        combinationArray[i] = i;
    }
}

__host__ __device__
void nextCombination(long n,
        long k,
        long* currentCombination) {
    if (currentCombination[0] == n - k) {
        return;
    }
    long i;
    for (i = k - 1; i > 0 && currentCombination[i] == n - k + i; --i);
    ++currentCombination[i];
    for (long j = i; j < k - 1; ++j) {
        currentCombination[j + 1] = currentCombination[j] + 1;
    }
}

__host__ __device__ void printCombination(long *currentCombination,
        long sizeComb) {
    printf("Combination: {");
    for (long i = 0; i < sizeComb; i++) {
        printf("%2d", currentCombination[i]);
        if (i < sizeComb - 1) {
            printf(", ");
        }
    }
    printf("}");
}

__host__ __device__
void prlongQueue(long *queue, long headQueue, long tailQueue) {
    printf("Queue(%d):{", tailQueue - headQueue);
    for (long i = headQueue; i <= tailQueue; i++) {
        printf("%2d", queue[i]);
        if (i < tailQueue) {
            printf(", ");
        }
    }
    printf("}\n");
}

__host__ __device__
long checkConvexityP3CSR(long *csrColIdxs, long nvertices,
        long *csrRowOffset, long sizeRowOffset,
        unsigned char *aux,
        long auxSize,
        long *currentCombination,
        long sizeComb, long idx) {
    //clean aux vector            
    for (long i = 0; i < auxSize; i++) {
        aux[i] = 0;
    }
    long closeCount = 0;
    long maxSizeQueue = MAX((auxSize / 2), MAX_DEFAULT_SIZE_QUEUE);
    long *queue = (long *) malloc(maxSizeQueue * sizeof (long));
    long headQueue = 0;
    long tailQueue = -1;

    for (long i = 0; i < sizeComb; i++) {
        tailQueue = (tailQueue + 1) % maxSizeQueue;
        queue[tailQueue] = currentCombination[i];
    }

    long countExec = 1;

    while (headQueue <= tailQueue) {
        if (verboseSerial) {
            printf("\nP3(k=%2d,c=%ld)-%ld: ", sizeComb, idx, countExec++);
            prlongQueue(queue, headQueue, tailQueue);
        }
        long verti = queue[headQueue];
        headQueue = (headQueue + 1) % maxSizeQueue;
        if (verboseSerial) {
            printf("\tv-rm: %d", verti);
        }

        if (aux[verti] < PROCESSED && verti < nvertices) {
            closeCount++;
            long end = csrColIdxs[verti + 1];
            for (long i = csrColIdxs[verti]; i < end; i++) {
                long vertn = csrRowOffset[i];
                if (vertn != verti && vertn < nvertices) {
                    unsigned char previousValue = aux[vertn];
                    if (previousValue < INCLUDED) {
                        aux[vertn] = aux[vertn] + NEIGHBOOR_COUNT_INCLUDED;
                    }
                    if (previousValue < INCLUDED && aux[vertn] >= INCLUDED) {
                        tailQueue = (tailQueue + 1) % maxSizeQueue;
                        queue[tailQueue] = vertn;
                        if (verboseSerial)
                            printf("\t v-inc: %d,", vertn);
                    }
                }
            }
            aux[verti] = PROCESSED;
        }
    }
    free(queue);
    return closeCount;
}

long checkConvexityP3(UndirectedCSRGraph *graph,
        unsigned char *aux,
        long auxSize,
        long *currentCombination,
        long sizeComb, long idx) {
    return checkConvexityP3CSR(graph->getCsrColIdxs(), graph->getVerticesCount(),
            graph->getCsrRowOffset(), graph->getSizeRowOffset(),
            aux, auxSize, currentCombination, sizeComb, idx);
}

void serialFindHullNumber(UndirectedCSRGraph *graph) {
    graph->begin_serial_time = clock();

    long nvs = graph->getVerticesCount();
    long k;
    unsigned char *aux = new unsigned char [nvs];
    long *currentCombination;

    long currentSize = 0;
    long maxSize = nvs;
    long sizeCurrentHcp3 = 0;

    bool found = false;

    while (currentSize < maxSize && !found) {
        currentSize++;
        k = currentSize;
        long maxCombination = maxCombinations(nvs, k);
        currentCombination = (long *) malloc(k * sizeof (long));
        initialCombination(nvs, k, currentCombination);

        for (long i = 0; i < maxCombination && !found; i++) {
            sizeCurrentHcp3 = checkConvexityP3(graph, aux, nvs, currentCombination, k, i);
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
    free(aux);
    graph->end_serial_time = clock();
}

__global__ void kernelFindHullNumber(long *csrColIdxs, long nvertices,
        long *csrRowOffset, long sizeRowOffset, long maxCombination,
        long k, long offset, long *result, unsigned char *aux,
        unsigned long *cacheMaxCombination) {
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    bool found = false;
    long *currentCombination = (long *) malloc(k * sizeof (long));
    long auxoffset = idx*nvertices;

    long sizeCurrentHcp3 = 0;
    long limmit = (idx + 1) * offset;
    if (limmit > maxCombination) {
        limmit = maxCombination;
    }
    long i = idx * offset;

    long maxSizeQueue = MAX(nvertices / 2, MAX_DEFAULT_SIZE_QUEUE);
    long *queue = (long *) malloc(maxSizeQueue * sizeof (long));

    initialCombination(nvertices, k, currentCombination, i, cacheMaxCombination);

    while (i < limmit && !found && !result[0]) {
        long headQueue = 0;
        long tailQueue = -1;
        sizeCurrentHcp3 = 0;

        for (long y = 0; y < nvertices; y++) {
            aux[auxoffset + y] = 0;
        }

        for (long j = 0; j < k; j++) {
            tailQueue = (tailQueue + 1) % maxSizeQueue;
            queue[tailQueue] = currentCombination[j];
        }

        while (headQueue <= tailQueue) {
            long verti = queue[headQueue];
            headQueue = (headQueue + 1) % maxSizeQueue;

            if (aux[auxoffset + verti] < PROCESSED) {
                sizeCurrentHcp3++;
                long end = csrColIdxs[verti + 1];
                long x = csrColIdxs[verti];

                for (; x < end; x++) {
                    long vertn = csrRowOffset[x];
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
        }
        found = (sizeCurrentHcp3 == nvertices);
        if (!found) {
            nextCombination(nvertices, k, currentCombination);
            i++;
        }
    }
    if (found) {
        result[0] = sizeCurrentHcp3;
        result[1] = (i - 1);
    }
    free(queue);
    free(currentCombination);
}

void parallelFindHullNumber(UndirectedCSRGraph *graph) {
    graph->begin_parallel_time = clock();
    clock_t parcial = clock();

    long nvs = graph->getVerticesCount();
    long k;
    long currentSize = 0;
    long maxSize = nvs;
    long result[2];
    result[0] = result[1] = 0;
    long* csrColIdxs = graph->getCsrColIdxs();
    long verticesCount = graph->getVerticesCount();
    long* csrRowOffset = graph->getCsrRowOffset();
    long sizeRowOffset = graph->getSizeRowOffset();

    long* csrColIdxsGpu;
    long* csrRowOffsetGpu;
    long *resultGpu;
    unsigned char *auxGpu;
    unsigned long *cacheMaxCombination;

    long numBytesClsIdx = sizeof (long)*(verticesCount + 1);
    cudaMalloc((void**) &csrColIdxsGpu, numBytesClsIdx);

    long numBytesRowOff = sizeof (long)*sizeRowOffset;
    cudaMalloc((void**) &csrRowOffsetGpu, numBytesRowOff);

    cudaMalloc((void**) &auxGpu, sizeof (unsigned char)*DEFAULT_THREAD_PER_BLOCK * nvs);
    cudaMalloc((void**) &cacheMaxCombination, sizeof (unsigned long)*nvs);
    cudaMemset(cacheMaxCombination, 0, sizeof (unsigned long)*nvs);

    long numBytesResult = sizeof (long)*2;
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

    while (currentSize < maxSize - 1 && !found) {
        currentSize++;
        k = currentSize;
        long maxCombination = maxCombinations(nvs, k);

        long threadsPerBlock = DEFAULT_THREAD_PER_BLOCK;
        if (DEFAULT_THREAD_PER_BLOCK > maxCombination) {
            threadsPerBlock = maxCombination / 3;
        }


        long offset = maxCombination / threadsPerBlock;

        kernelFindHullNumber << < 1, threadsPerBlock >>> (csrColIdxsGpu, verticesCount,
                csrRowOffsetGpu, sizeRowOffset, maxCombination, k, offset, resultGpu, auxGpu,
                cacheMaxCombination);

        cudaMemcpy(result, resultGpu, numBytesResult, cudaMemcpyDeviceToHost);
        found = (result[0] == nvs);
    }
    if (found) {
        printf("Result Parallel: S=%d-Comb(%d,%d) |S| = %d |Hcp3(S)| = |V(g)| = %d\n", result[1], nvs, k, k, result[0]);
    }
    cudaFree(resultGpu);
    cudaFree(csrRowOffsetGpu);
    cudaFree(csrColIdxsGpu);
    cudaFree(auxGpu);
    cudaFree(cacheMaxCombination);
    graph->end_parallel_time = clock();
}

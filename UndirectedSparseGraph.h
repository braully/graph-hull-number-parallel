/* 
 * File:   UndirectedSparseGraph.h
 * Author: braully
 *
 * Created on 10 de Junho de 2016, 09:52
 */
#define PROCESSED 3
#define INCLUDED 2
#define NEIGHBOOR_COUNT_INCLUDED 1
//
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <time.h>

class UndirectedSparseGraph {
public:
    clock_t begin_serial_time, end_serial_time;
    clock_t begin_parallel_time, end_parallel_time;

    UndirectedSparseGraph(const int size) {
        vertices = (int *) malloc(size * sizeof (int));
        adj_list = new std::list < int > [size];
        nVertices = size;
        for (int i = 0; i < size; i++) {
            vertices[i] = i;
        }
    }

    void addEdge(int v1, int v2) {
        if (v1 >= nVertices || v2 >= nVertices) {
            perror("Invalid vertices!");
        }
        adj_list[vertices[v1]].push_back(vertices[v2]);
        adj_list[vertices[v2]].push_back(vertices[v1]);
        //        vertices[v1].
    }

    void printGraph() {
        printf("\nVertices: {");
        for (int i = 0; i < nVertices; i++) {
            printf("v%d(d%ld)", vertices[i], adj_list[i].size());
            if (i < nVertices - 1) {
                printf(", ");
            }
        }
        printf("}\n");
    }

    int degree(int vert) {
        if (vert >= nVertices) {
            perror("Invalid vertices!");
        }
        return adj_list[vert].size();
    }

    int getVerticesCount() {
        return nVertices;
    }

    std::list<int> getAdjacency(int vert) {
        if (vert > nVertices) {
            perror("Invalid vertice");
        }
        return adj_list[vert];
    }

    void clear() {

    }
private:
    int *vertices;
    int nVertices;
    std::list < int >* adj_list;
};

class UndirectedCSRGraph {
public:
    clock_t begin_serial_time, end_serial_time;
    clock_t begin_parallel_time, end_parallel_time;

    clock_t getTotalTimeSerial() {
        return end_serial_time - begin_serial_time;
    }

    clock_t getTotalTimeParallel() {
        return end_parallel_time - begin_parallel_time;
    }

    UndirectedCSRGraph(int *csrClIdxs, int nVerts,
            int *csrRowOffst, int sizeRowOffst) {
        nVertices = nVerts;
        sizeRowOffset = sizeRowOffst;
        csrColIdxs = csrClIdxs;
        csrRowOffset = csrRowOffst;
    }

    UndirectedCSRGraph(UndirectedSparseGraph* graph) {
        nVertices = graph->getVerticesCount();
        sizeRowOffset = 0;
        csrColIdxs = new int[nVertices + 1];
        for (int i = 0; i < nVertices; i++) {
            csrColIdxs[i] = 0;
            sizeRowOffset = sizeRowOffset + graph->degree(i);
        }

        csrRowOffset = new int[sizeRowOffset];
        int idx = 0;
        for (int i = 0; i < nVertices; i++) {
            csrColIdxs[i] = idx;
            std::list<int>::iterator it;
            std::list<int> tmp = graph->getAdjacency(i);
            for (std::list<int>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
                //                if (verboseCsr) {
                //                    printf("csrRow[%d] = %d\n", idx, *it);
                //                }
                csrRowOffset[idx++] = *it;
            }
        }
        csrColIdxs[nVertices] = idx;
    }

    int *getCsrColIdxs() {
        return csrColIdxs;
    }

    int *getCsrRowOffset() {
        return csrRowOffset;
    }

    int getOffsetBeginNeighbors(int vert) {
        if (vert >= nVertices) {
            perror("Invalid vertices!");
        }
        return csrColIdxs[vert];
    }

    int getIdxAdjList(int idx) {
        return csrRowOffset[idx];
    }

    int getOffsetEndNeighbors(int vert) {
        if (vert >= nVertices) {
            perror("Invalid vertices!");
        }
        return csrColIdxs[vert + 1];
    }

    int degree(int vert) {
        if (vert >= nVertices) {
            perror("Invalid vertices!");
        }
        return csrColIdxs[vert + 1] - csrColIdxs[vert];
    }

    int getVerticesCount() {
        return nVertices;
    }

    int getSizeRowOffset() {
        return sizeRowOffset;
    }

    void printGraph() {
        printf("\nN. vertices: %d", nVertices);
        printf("\nSize RowOffset: %d", sizeRowOffset);
        printf("\nR: {");
        for (int i = 0; i < nVertices + 1; i++) {
            printf("%d", csrColIdxs[i]);
            if (i < nVertices) {
                printf(", ");
            }
        }
        printf("} \nC: {");
        for (int i = 0; i < sizeRowOffset; i++) {
            printf("%d", csrRowOffset[i]);
            if (i < sizeRowOffset - 1) {
                printf(", ");
            }
        }
        printf("}");
    }
    //    virtual ~UndirectedCSRGraph();
private:
    int nVertices;
    int sizeRowOffset;
    int *csrColIdxs;
    int *csrRowOffset;
};


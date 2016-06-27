#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <string>
#include <sstream>
#include "UndirectedSparseGraph.h"
#include "GraphHullNumber.h"

#define CHARACTER_INIT_COMMENT '#'

//References and Examples:
//http://www.boost.org/doc/libs/1_39_0/libs/graph/example/undirected.cpp
//Example: http://devblogs.nvidia.com/parallelforall/wp-content/uploads/2014/07/CSR.png
//Example: https://devblogs.nvidia.com/parallelforall/accelerating-graph-betweenness-centrality-cuda/


bool verboseGraph = false;
//bool verboseCsr = false;

void printHelp() {
    printf("\n\tp: Parallel execution");
    printf("\n\ts: Serial execution");
    printf("\n\th: Print this help message");
    printf("\n\tinput: Graph file in format CSR, see graph-test.txt");
}

int main(int argc, char** argv) {
    int opt = 0;
    //    char* strFile = "graph-test.txt";
    char* strFile = "graph-csr-2124643542179835849.txt";
    bool serial = false;
    bool parallel = false;
    bool verbose = false;

    if ((argc <= 1) || (argv[argc - 1] == NULL) || (argv[argc - 1][0] == '-')) {
        serial = true;
        parallel = true;
        //        printf("\nusage: graph-hull-number [-ps] [input]");
        //        printHelp();
        //        return 1;
    } else {
        strFile = argv[argc - 1];
    }

    while ((opt = getopt(argc, argv, "psv")) != -1) {
        switch (opt) {
            case 'p':
                parallel = true;
                break;
            case 's':
                serial = true;
                break;
            case 'v':
                verbose = true;
                break;
            case '?':
                printf("Unknow option: %c", char(opt));
                break;
        }
    }

    std::string line, strCArray, strRArray;
    std::ifstream infile(strFile);

    if (infile) {
        while (getline(infile, line)) {
            if (line.at(0) != CHARACTER_INIT_COMMENT) {
                if (strCArray.empty()) {
                    strCArray = line;
                } else if (strRArray.empty()) {
                    strRArray = line;
                }
            }
        }
    } else {
        printf("file '%s' not found!", strFile);
        return 1;
    }

    infile.close();

    if (strCArray.empty() || strRArray.empty()) {
        printf("Invalid file format");
        return 1;
    }

    std::stringstream stream(strCArray);
    std::vector<int> values;
    int n;
    while (stream >> n) {
        values.push_back(n);
    }
    strCArray.clear();

    int numVertices = values.size() - 1;
    int *colIdx = new int[numVertices + 1];
    std::copy(values.begin(), values.end(), colIdx);
    values.clear();
    stream.str("");

    std::stringstream stream2(strRArray);
    while (stream2 >> n) {
        values.push_back(n);
    }
    stream2.str("");
    strRArray.clear();

    int sizeRowOffset = values.size();
    int *rowOffset = new int[sizeRowOffset];
    std::copy(values.begin(), values.end(), rowOffset);
    values.clear();

    UndirectedCSRGraph csr(colIdx, numVertices, rowOffset, sizeRowOffset);
    if (verbose) {
        printf("Graph detail: ");
        csr.printGraph();
        printf("\n");
    }

    if (serial)
        serialFindHullNumber(&csr);
    if (parallel)
        parallelFindHullNumber(&csr);

    //    int n1 = 25, k = 14, idx = 1000;
    //    int* combinationArray = new int[25];
    //    initialCombination(n1, k, combinationArray, idx);
    //    printCombination(combinationArray, k);

    printf("\nTotal time serial: %ld\nTotal time parallel: %ld\n",
            csr.getTotalTimeSerial(),
            csr.getTotalTimeParallel());
    return 0;
}


//    UndirectedCSRGraph(int *csrClIdxs, int nVerts, int *csrRowOffst,
//            int sizeRowOffst)

//    int numVertices = 0;
//    int numCols = 0;
//    int numAdjsOffset = 0;
//
//    if (!strRArray.empty()) {
//        unsigned int lastI = 0;
//        for (unsigned int i = 0; i < strRArray.length(); i++)
//            if (strRArray.at(i) == ' ') {
//                if (i > lastI + 1) {
//                    numCols++;
//                }
//                lastI = i;
//            }
//    }
//
//    if (!strCArray.empty()) {
//        unsigned int lastI = 0;
//        for (unsigned int i = 0; i < strCArray.length(); i++)
//            if (strCArray.at(i) == ' ') {
//                if (i > lastI + 1) {
//                    numAdjsOffset++;
//                }
//                lastI = i;
//            }
//    }

//    int *colIdx = new int[numCols];
//    int *rowOffset = new int[numAdjsOffset];

//    UndirectedSparseGraph graph(0);
//    UndirectedCSRGraph csr(&graph);
//    csr.printGraph();

//    UndirectedSparseGraph graph(9);
//
//    graph.addEdge(0, 1);
//    graph.addEdge(0, 2);
//    graph.addEdge(0, 3);
//    graph.addEdge(1, 2);
//    graph.addEdge(2, 3);
//
//    graph.addEdge(3, 4);
//    graph.addEdge(3, 5);
//
//    graph.addEdge(4, 5);
//    graph.addEdge(4, 7);
//    graph.addEdge(4, 6);
//
//
//    graph.addEdge(5, 6);
//    graph.addEdge(5, 7);
//
//    graph.addEdge(6, 8);
//
//    graph.addEdge(7, 6);
//
//    if (verboseGraph)
//        graph.printGraph();


//    int nVertices = graph.getVerticesCount();
//    unsigned char* aux = new unsigned char[nVertices];
//    int currentCombination[3] = {0, 2, 3};
//    int sizeComb = 3;
//
//    printf("\n\nGrafo 1\n");
//    int sizeHcp3 = checkConvexityP3(&csr, aux, nVertices, currentCombination, sizeComb);
//    printf(": hcp3(G) = %d", sizeHcp3);
//
//    printf("\n\nGrafo 2\n");
//    //    currentCombination[0] = 2;
//    //    currentCombination[1] = 5;
//    //    currentCombination[2] = 8;
//    currentCombination[0] = 0;
//    currentCombination[1] = 5;
//    currentCombination[2] = 6;
//    sizeHcp3 = checkConvexityP3(&csr, aux, nVertices, currentCombination, sizeComb - 1);
//    printf(": hcp3(G) = %d", sizeHcp3);
//
//    printf("\n\nGrafo 3\n");
//    sizeHcp3 = checkConvexityP3(&csr, aux, nVertices, currentCombination, sizeComb);
//    printf(": hcp3(G) = %d", sizeHcp3);

//    int n = 9, k = 3;
//    int ncombs = maxCombinations(n, k);
//    printf("\ncomb(%d,%d)=%d", n, k, ncombs);
//    int *combs = new int[10];
//    initialCombination(n, k, combs);
//    for (int i = 0; i < ncombs; i++) {
//        printf("%3ld - ", i);
//        printCombination(combs, k);
//        nextCombination(n, k, combs);
//    }

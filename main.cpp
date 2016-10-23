//C++
#include <iostream>
//C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
//MPI
#include "mpi.h"

//Defines
#define ArrToMatrix(ROW, COL, NODES) (ROW) * (NODES) + (COL)
#define ROOT 0

int **theMatrix;
int *theMatrixData;
int **smallMatrix1;
int **smallMatrix2;
int **smallMatrix3;
int **smallMatrix4;

int myRank, cartRank;
int rowMatrix, columnMatrix;
int Q, nProc, nNodes;
int myRow, myCol;
MPI_Comm cartComm, rowComm, colComm;

int checkIfPossible( int p, int n){
    printf("Está a verificar se é possivel....\n");
    int tempQ = sqrt(p);
    if(tempQ * tempQ == p){
        if(n%tempQ == 0){
            return tempQ;
        }
    }
    printf("Não é possível aplicar o algoritmo Fox\n");
    return 1;
}

void dealWithInput(){
    printf("Inserir valores da matriz:\n");

    for(int i =0; i<nNodes; i++){
        for(int j=0; j<nNodes; j++){
            scanf("%d", &theMatrix[i][j] );
            if(theMatrix[i][j] == 0 && i != j){
                theMatrix[i][j] = -1;
            }
        }
    }
}


void prepareMatrixes(){
    printf("Entrou no prepareMatrixes\n");
    int i;
    int dims[2];
    int period[2] = {1,1};
    int cartCoords[2];

    MPI_Bcast(&nNodes, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    int* auxMatrix = new int[Q*Q];
    smallMatrix1 = new int*[Q];
    for(i=0; i<Q; i++){
        smallMatrix1[i] = &auxMatrix[i*Q];
    }

    auxMatrix = new int[Q*Q];
    smallMatrix2 = new int*[Q];
    for(i=0; i<Q;i++){
        smallMatrix2[i] = &auxMatrix[i*Q];
    }

    auxMatrix = new int[Q*Q];
    smallMatrix3 = new int*[Q];
    for(i=0; i<Q; i++){
        smallMatrix3[i] = &auxMatrix[i*Q];
    }

    auxMatrix = new int[Q*Q];
    smallMatrix4 = new int*[Q];
    for(i=0; i<Q; i++){
        smallMatrix4[i] = &auxMatrix[i*Q];
    }


    MPI_Dims_create(nProc, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &cartComm);
    MPI_Comm_rank(cartComm, &cartRank);

    MPI_Cart_coords(cartComm, cartRank, 2, cartCoords);
    myRow = cartCoords[0];
    myCol = cartCoords[1];



    printf("Saiu do prepareMatrixes\n");
}

void printMatrix(){
    for (int i = 0; i < nNodes; i++){
        for (int j = 0; j < nNodes; j++){
            printf("%d%c", theMatrix[i][j], j == nNodes - 1 ? '\n' : '\t');
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if(myRank == ROOT){
        //printf("nProc = %d\n", nProc);
        scanf("%d", &nNodes);
        //printf("Nodes = %d\n", nNodes);

        int Q = checkIfPossible(nProc, nNodes);
        if(Q == 1){
            MPI_Abort(MPI_COMM_WORLD, 0);
            return 1;
        }

        theMatrix = new int*[nNodes];
        theMatrixData = new int[nNodes*nNodes];

        for (int i = 0; i < nNodes ; ++i) {
            theMatrix[i] = &theMatrixData[i * nNodes];
        }

        dealWithInput();
        //printf("\nINFOWARS.COM\n");

    }

    MPI_Barrier(MPI_COMM_WORLD);
    double startTiem = MPI_Wtime();

    prepareMatrixes();


    double finish = MPI_Wtime();

    MPI_Finalize();
    return 0;
}

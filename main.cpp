//C++
#include <iostream>
#include <algorithm>
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
int **localMatrix;
int **smallMatrixA;
int **smallMatrixB;
int **smallMatrixC;

int myRank, cartRank;
int rowMatrix, columnMatrix;
int Q, nProc, nNodes;
int myRow, myCol;
MPI_Comm cartComm, rowComm, colComm;
MPI_Status status;


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
    int dimsSub[2] = {};

    MPI_Bcast(&nNodes, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    int* auxMatrix = new int[Q*Q];
    localMatrix = new int*[Q];
    for(i=0; i<Q; i++){
        localMatrix[i] = &auxMatrix[i*Q];
    }

    auxMatrix = new int[Q*Q];
    smallMatrixA = new int*[Q];
    for(i=0; i<Q;i++){
        smallMatrixA[i] = &auxMatrix[i*Q];
    }

    auxMatrix = new int[Q*Q];
    smallMatrixB = new int*[Q];
    for(i=0; i<Q; i++){
        smallMatrixB[i] = &auxMatrix[i*Q];
    }

    auxMatrix = new int[Q*Q];
    smallMatrixC = new int*[Q];
    for(i=0; i<Q; i++){
        smallMatrixC[i] = &auxMatrix[i*Q];
    }

    MPI_Dims_create(nProc, 2, dims);;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &cartComm);
    MPI_Comm_rank(cartComm, &cartRank);

    MPI_Cart_coords(cartComm, cartRank, 2, cartCoords);
    myRow = cartCoords[0];
    myCol = cartCoords[1];

    dimsSub[0] = 0;
    dimsSub[1] = 1;
    MPI_Cart_sub(cartComm, dimsSub, &rowComm);

    dimsSub[0] = 1;
    dimsSub[1] = 0;
    MPI_Cart_sub(cartComm, dimsSub, &colComm);

    printf("Saiu do prepareMatrixes\n");
}

void writeToMatrixes(){
    int rankToSend;
    MPI_Datatype qqMatrix;
    MPI_Type_vector(Q,Q,nNodes,MPI_INT, &qqMatrix);
    MPI_Type_commit(&qqMatrix);

    if(myRank == ROOT){
        for(int i =0; i<=Q; i+=Q){
            for(int j=0; j<=Q; j+=Q){
                int coordsToRank[2] = {i/Q, j/Q};

                MPI_Cart_rank(cartComm, coordsToRank, &rankToSend);

                if(rankToSend == ROOT){
                    for(int ii = i; ii< i+Q; ii++){
                        for(int jj = j; jj<j+Q; jj++){
                            localMatrix[i][j] = theMatrix[i][j];
                        }
                    }
                }
                else {
                    MPI_Send(&theMatrix[i][j], 1, qqMatrix, rankToSend, 1, cartComm);
                }

            }
        }
    }
    else{
        MPI_Recv(localMatrix, 1, qqMatrix, ROOT, 1, cartComm, &status);
    }

}



void multiplySomething(int** matrixA, int** matrixB, int** matrixC, int size){
    int i, j, k;
    int sum = 0;
    for(i=0; i<size; i++ ){
        for(j=0; j<size; j++){
            for(k=0; k<size; k++){
                if(matrixA[i][k] != -1 && matrixB[k][j] != -1 && matrixC[i][j] != -1) {
                    matrixC[i][j] = std::min(matrixC[i][j], matrixA[i][k] + matrixB[k][j]);
                }
                else if(matrixA[i][k] != -1 && matrixB[k][j] != -1 ){
                    matrixC[i][j] = matrixA[i][k] + matrixB[k][j];
                }
            }
        }
    }
}


void printMatrix(){
    for (int i = 0; i < nNodes; i++){
        for (int j = 0; j < nNodes; j++){
            printf("%d%c", theMatrix[i][j], j == nNodes - 1 ? '\n' : '\t');
        }
    }
}

int main(int argc, char *argv[]) {

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

        MPI_Bcast(&Q, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

        theMatrix = new int*[nNodes];
        theMatrixData = new int[nNodes*nNodes];

        for (int i = 0; i < nNodes ; ++i) {
            theMatrix[i] = &theMatrixData[i * nNodes];
        }

        dealWithInput();
        //printf("\nINFOWARS.COM\n");

    }

    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = MPI_Wtime();

    prepareMatrixes();
    writeToMatrixes();



    double finish = MPI_Wtime();

    MPI_Finalize();
    return 0;
}

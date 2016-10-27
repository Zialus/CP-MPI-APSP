//C++
#include <iostream>
#include <algorithm>
//C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//MPI
#include "mpi.h"

//Defines
#define ROOT 0

int **theMatrix;
int *theMatrixData;
int **localMatrix;
int **smallMatrixA;
int **smallMatrixB;
int **smallMatrixC;

int myRank, cartRank;
int rowMatrix, columnMatrix;
int Q, nProc, nNodes, nByQ;
int myRow, myCol;
MPI_Comm cartComm, rowComm, colComm;
MPI_Status status;

void fox();

int checkIfPossible(int p, int nNodes){
    printf("Está a verificar se é possivel....\n");
    int tempQ = sqrt(p);
    if(tempQ * tempQ == p){
        if(nNodes%tempQ == 0){
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

void APSP(){
    for (int d = 2; d <= 2*nNodes; d=d*2) {
        fox();
    }
}

void fox() {

    MPI_Datatype littleMatrix;
    MPI_Type_vector(nByQ*nByQ,1,1,MPI_INT, &littleMatrix);
    MPI_Type_commit(&littleMatrix);




    // Calculate indices of matrices above and below (on the same column)
    int source = (myRow + 1) % Q;
    int dest = (myRow + Q - 1) % Q;

    // Save their ranks on rankUP and rankDOWN
    int rankUP, rankDOWN;
    int coords[1];
    coords[0]= source;
    MPI_Cart_rank(colComm,coords,&rankUP);
    coords[0]=dest;
    MPI_Cart_rank(colComm,coords,&rankDOWN);


    for (int stage = 0; stage < Q; stage++) {

        int bcastROOT = ( myRow + stage) % Q;
        if (bcastROOT == myCol) {
            MPI_Bcast(localMatrix[0],1,littleMatrix,bcastROOT,rowComm);
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

void writeToMatrices(){

    // rank of process that will receive a particular subMatrix inside the CartGrid
    int rankToSend;

    MPI_Datatype subMatrix;
    MPI_Type_vector(nByQ,nByQ,nNodes,MPI_INT, &subMatrix);
    MPI_Type_commit(&subMatrix);

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
                    MPI_Send(&theMatrix[i][j], 1, subMatrix, rankToSend, 1, cartComm);
                }

            }
        }
    }
    else{
        MPI_Recv(localMatrix, 1, subMatrix, ROOT, 1, cartComm, &status);
    }

    MPI_Type_free(&subMatrix);
}



void localMultiply(int** matrixA, int** matrixB, int** matrixC, int size){
    int i, j, k;
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

        nByQ = nNodes / Q ;

        MPI_Bcast(&Q, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

        theMatrix = new int*[nNodes];
        theMatrixData = new int[nNodes*nNodes];

        for (int i = 0; i < nNodes ; ++i) {
            theMatrix[i] = &theMatrixData[i * nNodes];
        }

        dealWithInput();
        //printf("\n!!!!!INFOWARS DOT COM!!!!!!\n");
        //printf("\n!!!!!MY LIBIDO IS THRU THE ROOF!!!!!!\n");


    }

    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = MPI_Wtime();

    prepareMatrixes();
    writeToMatrices();


    double finishTime = MPI_Wtime();

    double elapsedTime = finishTime - startTime;
    std::cout << elapsedTime << std::endl;

    MPI_Finalize();
    return 0;
}

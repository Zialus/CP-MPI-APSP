//C++
#include <iostream>
//C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//MPI
#include "mpi.h"

//Defines
#define ArrToMatrix(ROW, COL, NODES) (ROW) * (NODES) + (COL)
#define ROOT 0

int *theMatrix;
int *localMatrix;
int nNodes;


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
  printf("Entrou no dealWithInput---------\n");
  theMatrix = (int*) malloc((nNodes*nNodes) * sizeof(int));
  int aux;
  for(int i =0; i<nNodes; i++){
    for(int j=0; j<nNodes; j++){
      scanf("%d", &theMatrix[ArrToMatrix(i,j,nNodes)] );
      if(theMatrix[ArrToMatrix(i,j,nNodes)] == 0 && i != j){
        theMatrix[ArrToMatrix(i,j,nNodes)] = -1;
      }
    }
  }
}


void printMatrix(){
  for (int i = 0; i < nNodes; i++){
    for (int j = 0; j < nNodes; j++){
      printf("%d%c", theMatrix[ArrToMatrix(i,j,nNodes)], j == nNodes - 1 ? '\n' : ' ');
    }
  }
}

int main(int argc, char *argv[]) {
  int myrank;
  int rowMatrix, columnMatrix;
  int Q, nProc;
  int myRow, myCol;

  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank == ROOT){
    printf("nProc = %d\n", nProc);
    scanf("%d", &nNodes);
    printf("Nodes = %d\n", nNodes);

    int Q = checkIfPossible(nProc, nNodes);
    if(Q == 1){
      MPI_Abort(MPI_COMM_WORLD, 0);
      return 1;
    }

    dealWithInput();
    printf("INFOWARS.COM\n");
    printMatrix();

  }



  MPI_Finalize();
  return 0;
}

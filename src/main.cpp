//C++
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
//C
#include <cmath>
#include <cstring>
//MPI
#include "mpi.h"
//Defines
#define ROOT 0

int** theMatrix;
int* theMatrixData;

int** localMatrix;
int** smallMatrixA;
int** smallMatrixB;
int** smallMatrixC;

int myRank, cartRank;
int Q, nProcs, nNodes, N_By_Q;
int myRow, myCol;
MPI_Comm cartComm, rowComm, colComm;
MPI_Status status;

template<typename InputIterator1, typename InputIterator2>
bool range_equal(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2) {
    while (first1 != last1 && first2 != last2) {
        if (*first1 != *first2) { return false; }
        ++first1;
        ++first2;
    }
    return (first1 == last1) && (first2 == last2);
}

bool compare_files(const std::string& filename1, const std::string& filename2) {
    std::ifstream file1(filename1);
    std::ifstream file2(filename2);

    std::istreambuf_iterator<char> begin1(file1);
    std::istreambuf_iterator<char> begin2(file2);

    std::istreambuf_iterator<char> end;

    return range_equal(begin1, end, begin2, end);
}

int checkIfPossible(int nProcs, int nNodes) {
    std::cout << "Checking if is possible to apply Fox algorithm..." << std::endl;

    double doubleQ = sqrt(nProcs);
    int tempQ = (int) doubleQ;

    if (tempQ != doubleQ) {
        perror("Can't apply Fox algorithm");
        perror("Number of processors is not a perfect square");
        return -1;
    }

    if (nNodes % tempQ != 0) {
        perror("Can't apply Fox algorithm");
        perror("Number of nodes is not divisible by the square root of the number of processors");
        return -1;
    }

    std::cout << "...Fox algorithm can be applied!" << std::endl;
    return tempQ;
}


void localMultiply(int** matrixA, int** matrixB, int** matrixC, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                if (matrixA[i][k] != -1 && matrixB[k][j] != -1 && matrixC[i][j] != -1) {
                    matrixC[i][j] = std::min(matrixC[i][j], matrixA[i][k] + matrixB[k][j]);
                } else if (matrixA[i][k] != -1 && matrixB[k][j] != -1) {
                    matrixC[i][j] = matrixA[i][k] + matrixB[k][j];
                }
            }
        }
    }
}

void freeMemory() {
    delete[] theMatrix;
    delete[] theMatrixData;

    delete[] localMatrix;
    delete[] smallMatrixA;
    delete[] smallMatrixB;
    delete[] smallMatrixC;
}

void prepareMatrices() {
    int i;
    int dimsCart[2] = {Q, Q};
    int period[2] = {1, 1};
    int cartCoords[2];
    int dimsSub[2] = {};

    int* auxMatrix;

    auxMatrix = new int[N_By_Q * N_By_Q];
    localMatrix = new int* [N_By_Q];
    for (i = 0; i < N_By_Q; i++) {
        localMatrix[i] = &auxMatrix[i * N_By_Q];
    }

    auxMatrix = new int[N_By_Q * N_By_Q];
    smallMatrixA = new int* [N_By_Q];
    for (i = 0; i < N_By_Q; i++) {
        smallMatrixA[i] = &auxMatrix[i * N_By_Q];
    }

    auxMatrix = new int[N_By_Q * N_By_Q];
    smallMatrixB = new int* [N_By_Q];
    for (i = 0; i < N_By_Q; i++) {
        smallMatrixB[i] = &auxMatrix[i * N_By_Q];
    }

    auxMatrix = new int[N_By_Q * N_By_Q];
    smallMatrixC = new int* [N_By_Q];
    for (i = 0; i < N_By_Q; i++) {
        smallMatrixC[i] = &auxMatrix[i * N_By_Q];
    }

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimsCart, period, 1, &cartComm);
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
}


void divideTheMatrix() {
    // rank of process that will receive a particular subMatrix inside the CartGrid
    int rankToSend;

    MPI_Datatype subMatrix;
    MPI_Type_vector(N_By_Q, N_By_Q, nNodes, MPI_INT, &subMatrix);
    MPI_Type_commit(&subMatrix);

    if (myRank == ROOT) {
        for (int i = 0; i < nNodes; i += N_By_Q) {
            for (int j = 0; j < nNodes; j += N_By_Q) {
                int coordsToRank[2] = {i / N_By_Q, j / N_By_Q};

                MPI_Cart_rank(cartComm, coordsToRank, &rankToSend);

                if (rankToSend == ROOT) {
                    for (int ii = i; ii < i + N_By_Q; ii++) {
                        for (int jj = j; jj < j + N_By_Q; jj++) {
                            localMatrix[ii][jj] = theMatrix[ii][jj];
                        }
                    }
                } else {
                    MPI_Send(&theMatrix[i][j], 1, subMatrix, rankToSend, 1, cartComm);
                }

            }
        }
    } else {
        MPI_Recv(localMatrix[0], N_By_Q * N_By_Q, MPI_INT, ROOT, 1, cartComm, &status);
    }

    MPI_Type_free(&subMatrix);
}


void conquerTheMatrix() {
    int rankToRecv;
    MPI_Datatype subMatrix;
    MPI_Type_vector(N_By_Q, N_By_Q, nNodes, MPI_INT, &subMatrix);
    MPI_Type_commit(&subMatrix);

    if (myRank == ROOT) {
        for (int i = 0; i < nNodes; i += N_By_Q) {
            for (int j = 0; j < nNodes; j += N_By_Q) {
                int coordsToRank[2] = {i / N_By_Q, j / N_By_Q};
                MPI_Cart_rank(cartComm, coordsToRank, &rankToRecv);

                if (rankToRecv == ROOT) {
                    for (int ii = i; ii < i + N_By_Q; ii++) {
                        for (int jj = j; jj < j + N_By_Q; jj++) {
                            theMatrix[ii][jj] = localMatrix[ii][jj];
                        }
                    }
                } else {
                    MPI_Recv(&theMatrix[i][j], 1, subMatrix, rankToRecv, 1, cartComm, &status);
                }

            }
        }
    } else {
        MPI_Send(localMatrix[0], N_By_Q * N_By_Q, MPI_INT, ROOT, 1, cartComm);
    }
}


void fox() {
    int i, j;
    MPI_Datatype littleMatrix;
    MPI_Type_vector(N_By_Q * N_By_Q, 1, 1, MPI_INT, &littleMatrix);
    MPI_Type_commit(&littleMatrix);

    // Calculate indices of matrices above and below (on the same column)
    int source = (myRow + 1) % Q;
    int dest = (myRow + Q - 1) % Q;

    // Save their ranks on rankUP and rankDOWN
    int rankUP, rankDOWN;
    int coords[1];
    coords[0] = source;
    MPI_Cart_rank(colComm, coords, &rankUP);
    coords[0] = dest;
    MPI_Cart_rank(colComm, coords, &rankDOWN);


    for (i = 0; i < N_By_Q; i++) {
        for (j = 0; j < N_By_Q; j++) {
            smallMatrixB[i][j] = localMatrix[i][j];
            smallMatrixC[i][j] = localMatrix[i][j];
        }
    }

    for (int stage = 0; stage < Q; stage++) {

        int bcastROOT = (myRow + stage) % Q;
        coords[0] = bcastROOT;
        int bcastROOTrank;
        MPI_Cart_rank(rowComm, coords, &bcastROOTrank);
        if (bcastROOT == myCol) {
            MPI_Bcast(localMatrix[0], 1, littleMatrix, bcastROOTrank, rowComm);
            localMultiply(localMatrix, smallMatrixB, smallMatrixC, N_By_Q);
        } else {
            MPI_Bcast(smallMatrixA[0], 1, littleMatrix, bcastROOTrank, rowComm);
            localMultiply(smallMatrixA, smallMatrixB, smallMatrixC, N_By_Q);
        }

        MPI_Sendrecv_replace(smallMatrixB[0], 1, littleMatrix, dest, 1, source, 1, colComm, &status);

    }

    MPI_Type_free(&littleMatrix);

    for (i = 0; i < N_By_Q; i++) {
        for (j = 0; j < N_By_Q; j++) {
            localMatrix[i][j] = smallMatrixC[i][j];
        }
    }

}

void APSP() {
    for (int d = 2; d <= 2 * nNodes; d = d * 2) {
        fox();
    }
}


void printMatrix(FILE* file) {
    std::cout << "Printing solution to file... ";
    for (int i = 0; i < nNodes; i++) {
        for (int j = 0; j < nNodes; j++) {
            if (theMatrix[i][j] == -1) {
                fprintf(file, "0%c", j == nNodes - 1 ? '\n' : ' ');
            } else {
                fprintf(file, "%d%c", theMatrix[i][j], j == nNodes - 1 ? '\n' : ' ');
            }
        }
    }
    std::cout << "DONE!" << std::endl;
}

void dealWithInput(int argc, char* const* argv) {
    if (argc == 3) {

        if (freopen(argv[1], "r", stdin) == nullptr) {
            perror("freopen() failed");
            MPI_Abort(MPI_COMM_WORLD, 1);
        };

        std::cout << "Insert number of nodes:" << std::endl;
        std::cin >> nNodes;

        Q = checkIfPossible(nProcs, nNodes);
        if (Q == -1) {
            MPI_Abort(MPI_COMM_WORLD, 0);
        }

        N_By_Q = nNodes / Q;

        theMatrix = new int* [nNodes];
        theMatrixData = new int[nNodes * nNodes];
        for (int i = 0; i < nNodes; ++i) {
            theMatrix[i] = &theMatrixData[i * nNodes];
        }

        std::cout << "Insert " << nNodes << " by " << nNodes << " values for the matrix:" << std::endl;

        for (int i = 0; i < nNodes; i++) {
            for (int j = 0; j < nNodes; j++) {
                std::cin >> theMatrix[i][j];
                if (theMatrix[i][j] == 0 && i != j) {
                    theMatrix[i][j] = -1;
                }
            }
        }

    } else {
        std::cout << "Wrong number of arguments" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void dealWithOutput(char* const* argv, double elapsedTime) {
    const char* outfile = "tmp.txt";
    std::string expected_file = argv[2];
    FILE* file = fopen(outfile, "w");
    printMatrix(file);
    fclose(file);
    std::cout << "Comparing the output... ";
    std::string result = compare_files(outfile, expected_file) ? "CORRECT OUTPUT" : "WRONG OUTPUT";
    std::cout << result << "!" << std::endl;
    std::cout << "The execution time was: " << elapsedTime << std::endl;
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == ROOT) {
        dealWithInput(argc, argv);
    }

    MPI_Bcast(&nNodes, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&N_By_Q, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&Q, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = MPI_Wtime();

    prepareMatrices();
    divideTheMatrix();
    APSP();
    conquerTheMatrix();

    MPI_Barrier(MPI_COMM_WORLD);
    double finishTime = MPI_Wtime();

    double elapsedTime = finishTime - startTime;

    if (myRank == ROOT) {
        dealWithOutput(argv, elapsedTime);
    }

    freeMemory();
    MPI_Finalize();
    return 0;
}

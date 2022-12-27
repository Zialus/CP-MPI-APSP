#pragma once

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

bool compare_files(const std::string&, const std::string&);

int checkIfPossible(int, int);

void dealWithOutput(char* const*, double);

void localMultiply(int**, int**, int**, int);

void freeMemory();

void prepareMatrices();

void divideTheMatrix();

void conquerTheMatrix();

void fox();

void APSP();

void printMatrix(FILE*);

void dealWithInput(int, char* const*);

#pragma once


#include "CRS.h"
#include <mpi.h>

int calcPartSize(const int procRank, const int procNum, const  int worldSize, int *start,
	 int *finish);
void sendData(const int procRank, const int procNum, const  int N, crsMatrix &sourceA, crsMatrix &sourceB, crsMatrix *receiveA,
	crsMatrix *receiveB,  int *rowSizes);
void mergeResult(const int procRank, const int procNum,  int *rowSizes, crsMatrix &C, crsMatrix *buffC);
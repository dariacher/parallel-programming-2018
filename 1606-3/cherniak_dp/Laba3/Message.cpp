#include "Message.h"

int _min(const int a, const int b) { if (a < b) return a; return b; }

int calcPartSize(const int procRank, const int procNum, const  int Size, int *start,	int *finish)
{
	*start = procRank * (Size / procNum) + _min(procRank, Size % procNum);
	*finish = (procRank + 1) * (Size / procNum) + _min(procRank + 1, Size % procNum) + _min(1, (procNum - 1 - procRank));

	return *finish - *start;
}

void sendData(const int rProc, const int nProc, const int N, crsMatrix &sA, crsMatrix &sB, crsMatrix *rA, crsMatrix *rB, int *rS)
{
	int *rowStarts, *rowFinishs;
	int *rowSizes;
	int *regularColSizes;

	MPI_Datatype RowIndexs;
	MPI_Datatype regularColIndexs;
	MPI_Datatype regularValIndexs;
	int *indices;

	rowStarts = (int *)malloc(sizeof(int)* nProc);
	rowFinishs = (int *)malloc(sizeof(int)* nProc);
	rowSizes = (int *)malloc(sizeof(int)* nProc);

	regularColSizes = (int *)malloc(sizeof(int)* nProc);

	indices = (int *)malloc(sizeof(int)* nProc);

	for (int rank = 0; rank < nProc; rank++)
		rowSizes[rank] = calcPartSize(rank, nProc, N + 1, &rowStarts[rank], &rowFinishs[rank]);

	if (rProc == 0)
	for (int rank = 0; rank < nProc; rank++)
		regularColSizes[rank] = sA.rowIndexs[rowFinishs[rank] - 1] - sA.rowIndexs[rowStarts[rank]];

	MPI_Bcast(regularColSizes, nProc, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&sB.N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&sB.NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);

	initialize_matrix(rowSizes[rProc] - 1, regularColSizes[rProc], rA);
	initialize_matrix(sB.N, sB.NZ, rB);

	if (rProc == 0) {
		memcpy(rB->rowIndexs, sB.rowIndexs, sizeof(int)* (sB.N + 1));
		memcpy(rB->columns, sB.columns, sizeof(int)* (sB.NZ));
		memcpy(rB->els, sB.els, sizeof(double)* (sB.NZ));
		memcpy(rS, rowSizes, sizeof(unsigned int)* nProc);
	}

	indices[0] = 0;
	for (int rank = 1; rank < nProc; rank++)
		indices[rank] = indices[rank - 1] + rowSizes[rank - 1] - 1;
	MPI_Type_indexed(nProc, (int *)rowSizes, indices, MPI_INT, &RowIndexs);
	MPI_Type_commit(&RowIndexs);//регистрация нового произвольного типа данных

	MPI_Scatterv(sA.rowIndexs, (int *)rowSizes, indices, MPI_INT, rA->rowIndexs, 1, RowIndexs, 0, MPI_COMM_WORLD);

	MPI_Bcast(rB->rowIndexs, sB.N + 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(rB->columns, sB.NZ, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(rB->els, sB.NZ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int rank = 1; rank < nProc; rank++)
		indices[rank] = indices[rank - 1] + regularColSizes[rank - 1];
	MPI_Type_indexed(nProc, (int *)regularColSizes, indices, MPI_INT, &regularColIndexs);
	MPI_Type_commit(&regularColIndexs);

	MPI_Type_indexed(nProc, (int *)regularColSizes, indices, MPI_DOUBLE, &regularValIndexs);
	MPI_Type_commit(&regularValIndexs);

	MPI_Scatterv(sA.columns, (int *)regularColSizes, indices, MPI_INT, rA->columns, 1, regularColIndexs, 0, MPI_COMM_WORLD);
	MPI_Scatterv(sA.els, (int *)regularColSizes, indices, MPI_DOUBLE, rA->els, 1, regularValIndexs, 0, MPI_COMM_WORLD);


	if (rProc != 0) {
		int *temp = new int[rA->N];
		temp[0] = 0;
		for (int i = 1; i < rA->N; ++i)
			temp[i] = rA->rowIndexs[i] - rA->rowIndexs[i - 1] + temp[i - 1];
		memcpy(rA->rowIndexs, temp, sizeof(int)* rA->N);
		rA->rowIndexs[rA->N] = rA->NZ;

		delete[] temp;
	}

	free(rowStarts), free(rowFinishs), free(rowSizes);
	free(indices); free(regularColSizes);
	MPI_Type_free(&RowIndexs); MPI_Type_free(&regularColIndexs);
	MPI_Type_free(&regularValIndexs);
}

void mergeResult(const int rProc, const int nProc, int *rowSizes, crsMatrix &C, crsMatrix *buffC)
{
	int N = 0, NZ = 0;
	int *sendbuf, scount;
	int *nzSizes = new int[nProc];
	int *rowDispls = new int[nProc];
	int *colDispls = new int[nProc];

	MPI_Gather(&C.NZ, 1, MPI_INT, nzSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rProc == 0) {
		for (int rank = 1; rank < nProc; rank++) rowSizes[rank]--;
		for (int rank = 0; rank < nProc; rank++) N += rowSizes[rank];
		for (int rank = 0; rank < nProc; rank++) NZ += nzSizes[rank];
		colDispls[0] = 0;
		for (int rank = 1; rank < nProc; rank++)
			colDispls[rank] = colDispls[rank - 1] + nzSizes[rank - 1];
		rowDispls[0] = 0;
		for (int rank = 1; rank < nProc; ++rank)
			rowDispls[rank] = rowDispls[rank - 1] + rowSizes[rank - 1];

		initialize_matrix(N - 1, NZ, buffC);
		sendbuf = C.rowIndexs;
		scount = C.N + 1;
	}
	else {
		sendbuf = &C.rowIndexs[1];
		scount = C.N;
	}
	MPI_Gatherv(sendbuf, scount, MPI_INT, buffC->rowIndexs, (int *)rowSizes, rowDispls, MPI_INT, 0, MPI_COMM_WORLD);
	if (rProc == 0) {
		for (int rank = 1; rank < nProc; rank++) rowSizes[rank] += rowSizes[rank - 1];
		for (int rank = 0; rank < nProc - 1; rank++)
		for (int i = rowSizes[rank]; i < rowSizes[rank + 1]; i++)
			buffC->rowIndexs[i] += buffC->rowIndexs[rowSizes[rank] - 1];
	}
	MPI_Gatherv(C.columns, C.NZ, MPI_INT, buffC->columns, nzSizes, colDispls, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(C.els, C.NZ, MPI_DOUBLE, buffC->els, nzSizes, colDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	delete[] nzSizes; delete[] colDispls; delete[] rowDispls;
}
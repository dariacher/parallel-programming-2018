#include "mpi.h"
#include <cstdlib>
#include <time.h>
#include <Windows.h>
#include <iostream>
#include "CRS.h"
#include "Message.h"


using namespace std;

int main(int argc, char **argv) {
	int nProc, rProc;
	int N, NZ;

	crsMatrix A, BT, C;
	crsMatrix buffA, buffB, buffBT, buffC;
	int *rowSizes = nullptr;

	double startTime, minStartTime;
	double finishTime, maxFinishTime;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rProc);
	N = atoi(argv[1]);
	NZ = atoi(argv[2]);

	if (rProc == 0) {
		rowSizes = (int *)malloc(sizeof(int)* nProc);
		generateRegularCrs(N, NZ, &buffA);
		generateSpecialCrs(N, NZ, &buffB);
		initialize_matrix(buffB.N, buffB.NZ, &buffBT);
		matr_transpose(&buffB, &buffBT);
		free_matrix(&buffB);
	}

	startTime = MPI_Wtime();

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);

	sendData(rProc, nProc, N, buffA, buffBT, &A, &BT, rowSizes);

	if (rProc == 0) { free_matrix(&buffA), free_matrix(&buffBT); }

	symbolicMult(N, &A, &BT, &C);
	numericMult(N, &A, &BT, &C);

	mergeResult(rProc, nProc, rowSizes, C, &buffC);

	finishTime = MPI_Wtime();

	MPI_Reduce(&startTime, &minStartTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&finishTime, &maxFinishTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rProc == 0) {
		if (N <= 10) {
			cout << "A: elemets = ";
			for (int i = 0; i < NZ; i++) {
				cout << A.els[i] << " ";
			}
			cout << endl;
			cout <<endl<< "A: columns = ";
			for (int i = 0; i < NZ; i++) {
				cout << A.columns[i] << " ";
			}
			cout << endl;

			cout << "B: elemets = ";
			for (int i = 0; i < NZ; i++) {
				cout << BT.els[i] << " ";
			}
			cout << endl;
			cout << endl << "B: columns = ";
			for (int i = 0; i < NZ; i++) {
				cout << BT.columns[i] << " ";
			}
			cout << endl;
			cout <<" C: elemets = ";
			for (int i = 0; i < NZ; i++) {
				cout << C.els[i] << " ";
			}
			cout << endl;
			cout << endl << "C: columns = ";
			for (int i = 0; i < NZ; i++) {
				cout << C.columns[i] << " ";
			}
		}
		cout << endl;
		cout << "Num of process = " << nProc << " N = " << N << " NZ = " << NZ;
		cout << " Time: " << maxFinishTime - minStartTime;
		free_matrix(&buffC);
	}

	MPI_Finalize();
	free_matrix(&A), free_matrix(&BT); free_matrix(&C);
	return 0;
}

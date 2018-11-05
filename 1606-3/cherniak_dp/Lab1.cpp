#include <iostream>
#include <mpi.h>
#include <string>
using namespace std;

void printMatrix(int *data, int row, int col) {
	if (row < 15 && col < 15) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				cout << data[i * col + j] << " ";
			}
			cout << endl;
		}
	}
}

int* createMatrix(int row, int col) {
	int *matrix;
	matrix = new int[row*col];
	return matrix;
}

void fullMatrix(int *matrix, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			matrix[i*col + j] = rand() % 100;
		}
	}
}
int maxSearch(int a, int b) {
	if (a >= b)
		return a;
	else return b;
}
int main(int argc, char **argv) {
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const int rows = stoi(string(argv[1]));
	const int cols = stoi(string(argv[2]));
	
	int *matrix = nullptr;
	if (rank == 0) {
	
		matrix = new int[rows * cols];
		fullMatrix(matrix, rows, cols);
		printMatrix(matrix, rows, cols);
	}

	int partSize = rows / size;
	int *vec = new int[cols * partSize];
	MPI_Scatter(matrix, cols * partSize, MPI_INT, vec, cols*partSize, MPI_INT, 0, MPI_COMM_WORLD);

	int *localMax = new int[partSize];

	for (int i = 0; i < partSize; i++) {
		int max = INT_MIN;
		for (int j = 0; j < cols; j++) {
			if (vec[i * cols + j] > max) {
				max = vec[i * cols + j];
			}
		}
		localMax[i] = max;
	}
	int *totalMax = nullptr;
	if (rank == 0) {
		totalMax = new int[rows];
	}
	MPI_Gather(localMax,  partSize, MPI_INT, totalMax, partSize, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		int tail = rows - size*partSize;
		for (int i = tail + 1; i < rows; i++) {
			int max = INT_MIN;
			for (int j = 0; j < cols; j++) {
				max = maxSearch(matrix[i*cols + j], max);
			}
			totalMax[i] = max;
		}

		for (int i = 0; i < rows; i++) {
			cout << totalMax[i] << endl;
		}
	}

	MPI_Finalize();
	return 0;
}

#include<iostream>
#include <mpi.h>
#include <math.h>
#include <random>
#include <string>
#include<ctime>
using namespace std;

void fillMatrix(int* matrix, int size) {
	for (int i = 0; i < size*size; i++) {
		matrix[i] = rand() % 50;
	}
}
void fillVector(int* vec, int size) {
	for (int i = 0; i < size; i++) {
		vec[i] = rand() % 50;
	}
}
void printVector(int *vec, int size) {
	if (size <= 15) {
		for (int i = 0; i < size; i++)
			cout << vec[i] << " ";
	}
}
void printMatrix(int *data, int size) {
	if (size < 15) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				cout << data[i * size + j] << " ";
			}
			cout << endl;
		}
	}
}
int main(int argc, char **argv) {
	int ProcRank, ProcNum;
	int *matrix = nullptr, *vec = nullptr, *result = nullptr, *linearRes = nullptr;
	int *block = nullptr, *displs = nullptr, *partRes = nullptr, *displRes = nullptr;
	double start, finish, linStart, linFinish;
	int size;
	int RestRows;//���������� ����� �������, ������� ��� �� ������������
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	size = atoi(argv[1]);
	int size1 = size;
	// ��������� ������ � ������������� �������� ������
	matrix = new int[size * size];
	vec = new int[size];
	result = new int[size];
	linearRes = new int[size];
	block = new int[ProcNum];
	displs = new int[ProcNum];
	partRes = new int[ProcNum];
	displRes = new int[ProcNum];
	if (ProcRank == 0) {
		//���������� ������� � �������
		fillMatrix(matrix, size);
		fillVector(vec, size);
		//������������ ������
		//������ - ���������� ���������, ������������ ��������
		int remainder = size1 % ProcNum;
		size1 -= remainder;						//����� �������� ��� �������
		int count = size1 / ProcNum;			// ���������� �����, ������������� ��������
		for (int i = 0; i < ProcNum; i++) {
			block[i] = count * size;			//���������� ������ ��� ������ ����
		}
		for (int i = 0; i < remainder; i++) {
			block[i] += size;					//������� �� �������
		}
		//������������ ��������
		int sm = 0;
		for (int i = 0; i < ProcNum; i++) {
			displs[i] = sm;
			sm += block[i];
		}
		//��������� ����� � �������� ���������� ��� ������� ����� 
		for (int i = 0; i < ProcNum; i++) {
			partRes[i] = block[i] / size;
		}//���������� ����� ��� ������� �����
			int Count = 0;
		for (int i = 0; i < ProcNum; i++) {
			displRes[i] = Count; //�������� ����������
			Count += partRes[i];
		}
	}
	start = MPI_Wtime();
	int partSize = (size1 / ProcNum + 1)*size;
	int *part = new int[partSize];
	MPI_Bcast(vec, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(block, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(partRes, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displRes, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	int *countRes = new int[block[ProcRank] / size];
	for (int i = 0; i < block[ProcRank] / size; i++) {
		countRes[i] = 0;
	}
	MPI_Scatterv(matrix, block, displs, MPI_INT, part, partSize, MPI_INT, 0, MPI_COMM_WORLD);
	//����������
	int d = 0, q = 0;
	for (int i = 0; i < block[ProcRank]; i++) {
		countRes[q] += part[i] * vec[d];
		d++;
		if (d == size) {
			d = 0;
			q++;
		}
	}
	//���� ����������
	MPI_Gatherv(countRes, q, MPI_INT, result, partRes, displRes, MPI_INT, 0, MPI_COMM_WORLD);
	finish = MPI_Wtime();
	//������ ����������
	if (ProcRank == 0) {
		cout << "Source matrix" << endl;
		printMatrix(matrix, size);
		cout << endl<< endl;
		cout << "Source vector" << endl;
		printVector(vec, size);
		cout <<endl<< "-------------------------------------"<< endl;
		cout << "result" << endl;
		printVector(result, size);
		cout << endl << "time: " << finish - start;
	}
	MPI_Finalize();
	return 0;
}
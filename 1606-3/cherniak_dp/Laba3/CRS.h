#pragma once

#include <math.h>
#include <vector>
#include <mpi.h>
using namespace std;

#define MAX_VAL 100
#define SUCCESS 0 
#define FAILURE 1
#define EMPTY -1
#define ZERO 0.000001
#define EPSILON 0.000001

struct crsMatrix {
	int N, NZ;

	//������ �������� ��������� ���������
	double *els;

	//������ ������ �������� ��� ������� ��������
	int *columns;

	//�������� ������ ����� �� ������ ������ ������ ������. ���������� ��������� ����� N+1
	int *rowIndexs;
};

void initialize_matrix(const int N, const int NZ, crsMatrix *mtx);
void free_matrix(crsMatrix *mtx);


//���������� ���������� �������
//� ������ ������ cntInRow ��������� ���������
void generateRegularCrs(const int N, const int cntInRow, crsMatrix *mtx);

//���������� ���������� �������
//����� ��������� ��������� � ������� ������ �� 1 �� cntInRow
//����� ����� - ���������� ��������
void generateSpecialCrs(const int N, const int cntInRow, crsMatrix *mtx);


//� ������ ���������� �����, ������� ����������� ����������, ������������ ���������� �������� ����������
//��������� ����� ������� � ���������� �������� ����������
//������ ������������� ���� - ��������� ������� ������� �. ���, ����� ���������������� ������� � ���������� ��� ������ ���������
//�������� ����������, ���� �� ���� �� ���� ���� ���������, ����������� � ����� � ��� �� �������.
//���� �������������� ������� ������� �, ����� ������ �������� �� ��� � ��������� ������ �� ������� ��� ������. ������������ ������ � � ������ ��
int symbolicMult(const int N, const crsMatrix *A, const crsMatrix *B, crsMatrix *C);
int numericMult(const int N, const crsMatrix *A, const crsMatrix *B, crsMatrix *C);

void matr_transpose(const crsMatrix *A, crsMatrix *AT);


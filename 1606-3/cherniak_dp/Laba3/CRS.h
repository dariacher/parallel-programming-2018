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

	//хранит значения элементов построчно
	double *els;

	//хранит номера столбцов для каждого элемента
	int *columns;

	//заменяет номера строк на индекс начала каждой строки. Количество элементов равно N+1
	int *rowIndexs;
};

void initialize_matrix(const int N, const int NZ, crsMatrix *mtx);
void free_matrix(crsMatrix *mtx);


//генерирует квадратную матрицу
//в каждой строке cntInRow ненулевых элементов
void generateRegularCrs(const int N, const int cntInRow, crsMatrix *mtx);

//генерирует квадратную матрицу
//число ненулевых элементов в строках растет от 1 до cntInRow
//закон роста - кубическая парабола
void generateSpecialCrs(const int N, const int cntInRow, crsMatrix *mtx);


//В рамках символьной части, которая выполняется однократно, производится построение портрета результата
//Численная часть состоит в заполнении портрета значениями
//Задача символической фазы - построить портрет матрицы С. Так, после транспонирования матрицы В необходимо при каждом умножении
//векторов определить, есть ли хотя бы одна пара элементов, находящихся в одном и том же столбце.
//Имея сформированный портрет матрицы С, можно просто пройтись по ней и вычислить каждый ее элемент как скаляр. произведение строки А и строки ВТ
int symbolicMult(const int N, const crsMatrix *A, const crsMatrix *B, crsMatrix *C);
int numericMult(const int N, const crsMatrix *A, const crsMatrix *B, crsMatrix *C);

void matr_transpose(const crsMatrix *A, crsMatrix *AT);


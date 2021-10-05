#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#define SIZE 3
#define EPS 0.2

//Создание двумерного массива
double** СreatingArray(int sizeLine, int sizeColumn)
{
	int i;
	double **Array;
	Array = (double **)malloc(sizeLine * sizeof(double*));
	for(i = 0; i < sizeLine; i++) 
		Array[i] = (double *)malloc(sizeColumn * sizeof(double));
	return Array;
}
//Печать двумерного массива
void PrintDoubleArray(double **Array, int siseLine, int sizeColumn)
{
	int i, j;
	for(i = 0; i < siseLine; i++)
	{
		for(j = 0; j < sizeColumn; j++)
			printf("%e   ", Array[i][j]);
		printf("\n");
	}
}
//Освобождение двумерного массива
void FreeDoubleArray(double **Array, int size)
{
	int i;
	for(i = 0; i < size; i++) 
		 free(Array[i]);
	free(Array);
}

//рандомное число
double Random(double a, double b)
{
	double x;
	int number = 0;
	number = rand();
	x = a + (b - a) * number / ( number + 1);
	return x;
}

//создание диагональной матрицы с рандомными числами
double** CreatingD(size)
{
	double **D;
	int i;
	D = СreatingArray(size, size);
	for (i = 0; i < SIZE ; i++)
		D[i][i] = Random(1 - EPS, 1 + EPS);
	return D;
}

//создание столбца рандомных чисел
double** CreatingW(int size)
{
	double **W;
	int i;
	W = СreatingArray(size, 1);
	for (i = 0; i < size ; i++)
		W[i][0] = Random(1 - EPS, 1 + EPS);
	return W;
}

//транспонированние столбца
double** CreatingWT(double **W, int size)
{
	double **WT;
	int i;
	WT = СreatingArray(1, size);
	for (i = 0; i < size ; i++)
		WT[0][i] = W[i][0];
	return W;
}
//Создание матрицы А
double** CreatingA(int size)
{
	double **A, **D, **W, **WT;
	A = СreatingArray(size, size);
	D = CreatingD(size);
	W = CreatingW(size);
	WT = CreatingWT(W, size);
	return A;
}
////////////////////////////////////////////////////////////////
double* Decision(double **U, double **L, double **B, int size)
{
	double *Y, *X, sum;
	int i, k;
	Y = (double*)malloc(size * sizeof(double));
	X = (double*)malloc(size * sizeof(double));
	Y[0] = B[0][0];
	for(i = 1; i < size; i++)
	{
		sum = 0;
		for (k = 0; k < i; k++)
                sum += L[i][k] * Y[k];
		Y[i] = B[i][0] - sum;
	}
	X[size-1] = Y[size-1]/U[size-1][size-1];
	for(i = 1; i < size; i++)
	{
		sum = 0;
		for (k = 0; k < i; k++)
                sum += U[size-1-i][size-1-k] * X[size-1-k];
		X[size-1-i] = (Y[size-1-i] - sum)/U[size-1-i][size-1-i];
	}
	printf("\nВектор B:\n");
	for(i = 0; i < size; i++)
		printf("%e   ", B[i][0]);
	/*printf("\nY\n");
	for(i = 0; i < size; i++)
		printf("%e   ", Y[i]);*/
	printf("\nРешение X:\n");
	for(i = 0; i < size; i++)
		printf("%e   ", X[i]);
	printf("\n");
	free(Y);
	return X;
}
void LU(double **A, double **L, double **U, int size)
{
	double sum;
	int i, j, k;
	printf("Матрица A:\n");
	PrintDoubleArray(A, size, size);
	for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            U[0][i] = A[0][i];
            L[i][0] = A[i][0] / U[0][0];
            sum = 0;
            for (k = 0; k < i; k++)
                sum += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - sum;
            if (i > j)
                L[j][i] = 0;
            else
            {
                sum = 0;
                for (k = 0; k < i; k++)
                    sum += L[j][k] * U[k][i];
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }
	printf("Матрица L:\n");
	PrintDoubleArray(L, size, size);
	printf("Маnрица U:\n");
	PrintDoubleArray(U, size, size);
}
double** Vozmushcheniye(double** B, int size)
{
	double** B1 = B;
	int i, j;
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			B[i][j] = B[i][j] + Random(-B[i][j]/100, B[i][j]/100);
	return B1;
}
void Nevyazka(double* X, double* X1, int size)
{
	int i;
	printf("Вектор невязки:\n");
	for(i = 0; i < size; i++)
		printf("%e   ", X[i]-X1[i]);
	printf("\n");

}
////////////////////////////////////////////////////////////////

main()
{
	double **A, **A1, **B, **B1, **U, **L, *X, *X1;
	setlocale(LC_ALL, "Rus");
	A = СreatingArray(SIZE, SIZE);
	A1 = СreatingArray(SIZE, SIZE);
	B = СreatingArray(SIZE, 1);
	B1 = СreatingArray(SIZE, 1);
	L = СreatingArray(SIZE, SIZE);
	U = СreatingArray(SIZE, SIZE);
	//Матрица А
    A[0][0]=7; A[0][1]=1; A[0][2]= 3;
    A[1][0]=2; A[1][1]=4; A[1][2]=1;
    A[2][0]=5; A[2][1]=8; A[2][2]=10;
	//Столбец В
	B[0][0]=1;
	B[1][0]=5;
	B[2][0]=2;

	LU(A, L, U, 3);
	X = Decision(U, L, B, SIZE);
	B1 = Vozmushcheniye(B, SIZE);
	printf("\nВнесем возмущение в вектор В, получим новый вектор В1 и решим систему для нового вектора В");
	X1 = Decision(U, L, B1, SIZE);
	Nevyazka(X, X1, SIZE);

	printf("\nВнесем возмущение в матрицу А, получим новый матрицу А1 и решим систему для новой матрицы А1");
	A1 = Vozmushcheniye(A, SIZE);
	LU(A, L, U, 3);
	X1 = Decision(U, L, B, SIZE);
	Nevyazka(X, X1, SIZE);
	return 0;
}
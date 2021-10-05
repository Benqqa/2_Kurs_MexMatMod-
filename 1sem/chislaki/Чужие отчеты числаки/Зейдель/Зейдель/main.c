#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#pragma warning(disable: 4996)
#define SIZE 3
#define size 3
#define M 3
#define EPS 0.000001

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
			printf("%e ", Array[i][j]);
		printf("\n");
	}
}
////Освобождение двумерного массива
//void FreeDoubleArray(double **Array, int size)
//{
//	int i;
//	for(i = 0; i < size; i++) 
//		 free(Array[i]);
//	free(Array);
//}


void Nevyazka(double **A, double *B, double* X)
{
	double *C, max = 0;
	int i, j;
	C = (double *)malloc(size * sizeof(double));
	printf("\n\nNevyazka A B:\n");
	for(i = 0; i < size; i++)
	{
		C[i] = 0;
		for(j = 0; j < size; j++)
			C[i] += A[i][j] * X[j];
		C[i] = C[i] - B[i];
		if(fabs(C[i]) > max)
			max = fabs(C[i]);
	}
	printf("%e\n ", max);
	free(C);
}

double Norma(double** A)
{
	int i, j;
	double norm = 0,normMax = 0;
	for (i = 0; i < size; i++)
	{
		norm = 0;
		for (j = 0; j < size; j++)
			norm += fabs(A[i][j]);
		if (norm > normMax)
			normMax = norm;
	}
	return normMax;
}

//double Norma(double** A)
//{
//	int i, j;
//	double norm = 0;
//	for (i = 0; i < size; i++)
//	{
//		for (j = 0; j < size; j++)
//			norm += A[i][j] * A[i][j];
//	}
//	norm = sqrt(norm);
//	return norm;
//}

void CreateAlfBet(double **A, double *B, double **alfa, double *betta)
{
	int i, j;
	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
			alfa[i][j] = (-1)*A[i][j] / A[i][i];
		alfa[i][i] = 0;
		betta[i] = B[i] / A[i][i];
	}
}

int MethodZeidel(double** alfa, double *betta, double *x1, int param) // 0 - норма < 1, 1 - At
{
	int i, j;
	int tmp = 0, count = 0;
	double X, q, norma;
	double *x2;
	x2 = (double *)malloc(SIZE * sizeof(double));
	for (i = 0; i < SIZE; i++)
		x2[i] = x1[i] = betta[i];
	norma = Norma(alfa);
	if (param = 1)
		norma = 0.5 ;
	q = fabs(norma / (1 - norma)) * EPS;
	//printf("\nq = %lf\n", q);
	while (tmp != SIZE)
	{
		tmp = 0;
		for (i = 0; i < SIZE; i++)
		{
			X = 0;
			for (j = 0; j < SIZE; j++)
				X += alfa[i][j] * x1[j];
			x1[i] = X + betta[i];
			if (fabs(x1[i] - x2[i]) < q)
				tmp++;
		}
		for (i = 0; i < SIZE; i++)
			x2[i] = x1[i];
		count++;
	}
	return count;
}

void Transposition(double **A, double **B)
{
	int i, j;
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			B[i][j] = A[j][i];
}

void MartixVectorMultiplication(double **A, double *B, double *C)
{
	int i, j;
	for (i = 0; i < SIZE; i++)
	{
		C[i] = 0;
		for (j = 0; j < SIZE; j++)
		{

			C[i] += A[i][j] * B[j];
		}
	}
}
  
void MartixMultiplication(double **A, double **B, double **C)
{
	int i, j, k;
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
		{
			C[i][j] = 0;
			for (k = 0; k < SIZE; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
}

int main(void)
{
	int i, j, count, k, par;
	double norma;
	double **A, **alfa, **At, **AA;
	double *B, *betta, *x1, *BB, *d;

	A = СreatingArray(SIZE, SIZE);
	alfa = СreatingArray(SIZE, SIZE);
	At = СreatingArray(SIZE, SIZE);
	AA = СreatingArray(SIZE, SIZE);
	B = (double *)malloc(SIZE * sizeof(double));
	betta = (double *)malloc(SIZE * sizeof(double));
	x1 = (double *)malloc(SIZE * sizeof(double));
	BB = (double *)malloc(SIZE * sizeof(double));
	d = (double *)malloc(SIZE * sizeof(double));

	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			scanf("%lf", &A[i][j]); 
	for (i = 0; i < SIZE; i++)
	{
		scanf("%lf", &B[i]);
		x1[i] = B[i];
	}
	CreateAlfBet(A, B, alfa, betta);
	par = 0;
	if (Norma(alfa) > 1)
	{
		Transposition(A, At);
		MartixVectorMultiplication(At, B, BB);
		MartixMultiplication(At, A, AA);
		CreateAlfBet(AA, BB, alfa, betta);
		par = 1;
	}
	for (i = 0; i < SIZE; i++)
	{
		k = 0;
		printf("X(%d) = ", i);
		for (j = 0; j < SIZE; j++)
		{
			printf("%lf*x(%d) ", alfa[i][j],k );
			k++;
		}
		printf("+%lf \n",betta[i]);
	}
	norma = Norma(alfa);
	printf("\nX : ");
	count = MethodZeidel(alfa, betta, x1,par);
	for (i = 0; i < SIZE; i++)
		printf("%lf ", x1[i]);
	printf("\n\nCount of iterations : %d", count);
	Nevyazka(A, B, x1);
	return 0;
}
/*
1 4 7
8 2 3
6 7 10
1 2 3
*/

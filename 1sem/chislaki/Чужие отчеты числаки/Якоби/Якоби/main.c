#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#pragma warning(disable: 4996)
#define SIZE 10
#define size 10
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




//double Norma(double** A)
//{
//	int i, j;
//	double norm = 0,normMax = 0;
//	for (i = 0; i < size; i++)
//	{
//		norm = 0;
//		for (j = 0; j < size; j++)
//			norm += fabs(A[i][j]);
//		if (norm > normMax)
//			normMax = norm;
//	}
//	return normMax;
//}

double Norma(double** A)
{
	int i, j;
	double norm = 0;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
			norm += A[i][j] * A[i][j];
	}
	norm = sqrt(norm);
	return norm;
}

void Transposition(double **A, double **B)
{
	int i, j;
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			B[i][j] = A[j][i];
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

//Поиск максимального элемента 
double MaxFind(double **A)
{
	int i, j;
	double max = A[0][1];
	for(i = 0; i < SIZE; i++)
		for(j = 0; j < SIZE; j++)
			if(A[i][j] > max && i != j)
				max = A[i][j];
	return max;
}
//Невязка
void Diff(double **A, double **B) 
{
	double max=0;
	double *res, *res1;
	int s = 0, i, j, k;
	res = (double*)malloc(SIZE * sizeof(double));
	res1 = (double*)malloc(SIZE * sizeof(double));
	while (s < size)
	{
		for (i=0; i<size; i++)
		{
			for (j=0; j<1; j++)
			{
				res[i] = 0;
				for (k=0; k<size; k++)
					res[i]+= A[i][k]*B[k][s];
			}
		}
		for (i=0; i<size; i++) 
			res1[i] = B[i][s]*A[s][s];
		for (i=0; i<size; i++) 
			if(res[i]-res1[i]>max) 
				max = res[i]-res1[i];
		s++;
	}
	printf("\n\n||Ax-lamb*x|| :%e\n\n", max);
}
//Матрица поворота
void MatrixRotation(double **G, double **GT, double **A)
{
	int i, j, ii=0, jj=0;
	double sum, max = 0, tan, sin, cos;
	for (i = 0; i < SIZE; i++)
	{
		sum = 0;
		for (j = 0; j < SIZE; j++)
		{
			if (i != j) 
				sum += (A[i][j])*(A[i][j]);
		}
		if (max < sum)
		{
			max = sum;
			ii = i;
		}
	}
	max = 0;
	for (j = 0; j < size; j++)
	{
		if (j != ii)
		{
			if (fabs(A[ii][j])>fabs(max))
			{
				max = A[ii][j];
				jj = j;
			}
		}
	}
	for (i=0; i<size; i++)
	{
		for (j=0; j<size; j++)
		{
			if (i!=j) 
				G[i][j] = 0;
			else 
				G[i][j] = 1;
		}
	}
	if ((A[jj][jj] - A[ii][ii])!=0)
	{
		tan = 2*A[ii][jj]/(A[jj][jj] - A[ii][ii]);
		//printf("fi = %g\n", 0.5*atan(tan));
		sin = sqrt(0.5*(1-1/(sqrt(1+tan*tan))));
		cos = sqrt(0.5*(1+1/(sqrt(1+tan*tan))));
		if (tan < 0) sin=-sin;
	}
	else
	{
		sin = sqrt(0.5);
		cos = sqrt(0.5);
	}
	G[ii][ii] = cos;
	G[jj][ii] = -sin;
	G[ii][jj] = sin;
	G[jj][jj] = cos;
	Transposition(G, GT);
}

//Решение Якоби
double **Jacobi(double **A, double **G, double **GT, double **B)
{
	double sum = 0, max = 0, **C, **CB, **RES, **HELP;
	int i, j, ii, jj, count = 0;
	C = СreatingArray(SIZE, SIZE);
	CB = СreatingArray(SIZE, SIZE);
	RES = СreatingArray(SIZE, SIZE);
	for(i = 0; i < SIZE; i++)
			for(j = 0; j < SIZE; j++)
				RES[i][j] = 1;
	HELP = СreatingArray(SIZE, SIZE);
	while (MaxFind(A) > EPS)
	{
		for (i = 0; i < SIZE; i++)
		{
			sum = 0;
			for (j = 0; j < SIZE; j++)
			{
				if (i != j) 
					sum += (A[i][j])*(A[i][j]);
			}
			if (max < sum)
			{
				max = sum;
				ii = i;
			}
		}
		max = 0;
		for (j = 0; j < size; j++)
		{
			if (j != ii)
			{
				if (fabs(A[ii][j])>fabs(max))
				{
					max = A[ii][j];
					jj = j;
				}
			}
		}
		MatrixRotation(G, GT, A);
		if(count == 0)
			for(i = 0; i < SIZE; i++)
				for(j = 0; j < SIZE; j++)
					RES[i][j] = G[i][j];
		else
		{
			MartixMultiplication (RES, G, HELP);
			for(i = 0; i < SIZE; i++)
				for(j = 0; j < SIZE; j++)
					RES[i][j] = HELP[i][j];
		}
		MartixMultiplication (B, G, CB);
		for(i = 0; i < SIZE; i++)
			for(j = 0; j < SIZE; j++)
				B[i][j] = CB[i][j];
		MartixMultiplication (GT, A, C);
		MartixMultiplication (C, G, A);
		count++; 
	}
	printf("\nCount of iteration: %i\n", count);
	return RES;
}


int main(void)
{
	int i, j;
	double **A, **GT, **G, **B, **RES;

	A = СreatingArray(SIZE, SIZE);
	G = СreatingArray(SIZE, SIZE);
	GT = СreatingArray(SIZE, SIZE);
	B = СreatingArray(SIZE, SIZE);

	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			scanf("%lf", &A[i][j]); 

	RES = Jacobi(A, G, GT, B);
	for(i = 0; i < SIZE; i++)
		printf("\nlamb%i : %e", i, A[i][i]);
	Diff(A, G);
	printf("X1		X2		 X3\n");
	PrintDoubleArray(RES, SIZE, SIZE);
	return 0;
}
/*
1 8 1
8 5 2
1 2 13
*/
/* 
4 2 7
2 -3 5
7 5 1
*/
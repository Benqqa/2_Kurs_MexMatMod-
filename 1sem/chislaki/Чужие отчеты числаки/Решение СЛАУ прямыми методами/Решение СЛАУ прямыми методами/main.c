#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#pragma warning(disable: 4996)
#define SIZE 10


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
//рандомное число
double Random(double a, double b)
{
	double x;
	int number = 0;
	number = rand();
	x = a + (b - a) * number / ( number + 1);
	return x;
}

////////////////////////////////////////////////////////////////
double* Decision(double **L, double **U, double *B, int size)
{
	double *Y, *X, sum;
	int i, k;
	Y = (double*)malloc(size * sizeof(double));
	X = (double*)malloc(size * sizeof(double));
	Y[0] = B[0];
	for(i = 1; i < size; i++)
	{
		sum = 0;
		for (k = 0; k < i; k++)
                sum += L[i][k] * Y[k];
		Y[i] = B[i] - sum;
	}
	X[size-1] = Y[size-1]/U[size-1][size-1];
	for(i = 1; i < size; i++)
	{
		sum = 0;
		for (k = 0; k < i; k++)
                sum += U[size-1-i][size-1-k] * X[size-1-k];
		X[size-1-i] = (Y[size-1-i] - sum)/U[size-1-i][size-1-i];
	}
	/*printf("\nB\n");
	for(i = 0; i < size; i++)
		printf("%e ", B[i]);*/
	printf("\nX:\n");
	for(i = 0; i < size; i++)
		printf("%e\n ", X[i]);
	printf("\n");
	free(Y);
	return X;
}
void LU(double **A, double **L, double **U, int size)
{
	double sum;
	int i, j, k;
	
	/*for (i = 0; i < size; i++)
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
    }*/
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			U[i][j] = A[i][j];
	for(k = 1; k < size; k++)
	{
		for(i = k-1; i < size; i++)
			for(j = i; j < size; j++)
				L[j][i]=U[j][i]/U[i][i];
		for(i = k; i < size; i++)
			for(j = k-1; j < size; j++)
				U[i][j]=U[i][j]-L[i][k-1]*U[k-1][j];
	}
	printf("\nmatrica L:\n");
	PrintDoubleArray(L, size, size);
	printf("\nmatrica U:\n");
	PrintDoubleArray(U, size, size);
}
double* Nevyazka(double **A, double *B, double* X, int size)
{
	double* C;
	int i, j;
	C = (double *)malloc(size * sizeof(double));
	for(i = 0; i < size; i++)
	{
		C[i] = 0;
		for(j = 0; j < size; j++)
			C[i] += A[i][j] * X[j];
		C[i] = C[i] - B[i];
	}
	return C;
}

double** VozmushcheniyeA(double** A, int size)
{
	double** A1;
	int i, j;
	A1 = СreatingArray(SIZE, SIZE);
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			A1[i][j] = A[i][j] + Random(-A[i][j]/100, A[i][j]/100);
	return A1;
}

double* VozmushcheniyeB(double* B, int size)
{
	double* B1;
	int i;
	B1 = (double *)malloc(SIZE * sizeof(double));
	for(i = 0; i < size; i++)
		B1[i] = B[i] + Random(-B[i]/100, B[i]/100);
	return B1;
}

double* Sum(double* X, double* X1, int size)
{
	int i;
	double* C;
	C = (double *)malloc(SIZE * sizeof(double));
	for(i = 0; i < size; i++)
		C[i] = X[i]+X1[i];
	return C;
}

double* Diff(double* X, double* X1, int size)
{
	int i;
	double* C;
	C = (double *)malloc(SIZE * sizeof(double));
	for(i = 0; i < size; i++)
		C[i] = X[i]-X1[i];
	return C;
}

double** DiffMatr(double** X, double** X1, int size)
{
	int i, j;
	double** C;
	C = СreatingArray(size, size);
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			C[i][j] = (X[i][j]-X1[i][j]);
	return C;
}

double** DiffMatrA1(double** X, double** X1, int size)
{
	int i, j;
	double** C;
	C = СreatingArray(size, size);
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			C[i][j] = (X[i][j]-X1[i][j])/1000;
	return C;
}

//double Norma(double** A, int size)
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

double Norma(double** A, int size)
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

double NormaStolb(double* A, int size)
{
	int i;
	double norm = 0;
	for (i = 0; i < size; i++)
		norm += A[i] * A[i];
	norm = sqrt(norm);
	return norm;
}

double SearchK1(double* B, double* B1, double* X, double* X1, int size)
{
	double normB, normDB, normX, normDX;
	normB = NormaStolb(B, size);
	normDB = NormaStolb(Diff(B, B1, size), size);
	normX = NormaStolb(X, size);
	normDX = NormaStolb(Diff(X, X1, size), size);
	printf("\n\nK1: %e\n",  normB * normDX / (normDB * normX));
	return  normB * normDX / (normDB * normX);
}

double SearchK2(double** A, double** A1, double* X, double* X1, int size)
{
	double normA, normDA, normDX, normXDX;
	normA = Norma(A, size);
	normDA = Norma(DiffMatr(A, A1, size), size);
	normDX = NormaStolb(Diff(X, X1, size), size); 
	normXDX = NormaStolb(Sum(X, Diff(X, X1, size), size), size);
	Diff(X, X1, size);
	printf("\n\nK2: %e\n",  normA * normDX / (normDA * normXDX));
	return  normA * normDX / (normDA * normXDX);
}

double SearchNK2(double** A, double** A1, double* X, double* X1, int size)
{
	double normA, normDA, normDX, normXDX;
	normA = 1.8320635;
	normDA = Norma(DiffMatrA1(A, A1, size), size);
	normDX = NormaStolb(Diff(X, X1, size), size); 
	normXDX = NormaStolb(Sum(X, Diff(X, X1, size), size), size);
	Diff(X, X1, size);
	printf("\n\nK2: %e\n",  normA * normDX / (normDA * normXDX));
	return  normA * normDX / (normDA * normXDX);
}

void Write(FILE *f, double** A, int size)
{
	int i, j;
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
			fprintf(f, "%lf ", A[i][j]);

}

void WriteStolb(FILE *f, double* A, int size)
{
	int i, j;
	for(i = 0; i < size; i++)
		fprintf(f, "%lf ", A[i]);
	fprintf(f, "\n");

}
////////////////////////////////////////////////////////////////

main()
{
	int i, j;
	double **A, **A1, *B, *B1, **U, **L, *X, *X1, *X2, *Nev, **C, norma;
	FILE *f, *fMatlab;

	A = СreatingArray(SIZE, SIZE);
	B = (double *)malloc(SIZE * sizeof(double));
	L = СreatingArray(SIZE, SIZE);
	U = СreatingArray(SIZE, SIZE);
	f = fopen("matrix.txt", "r");
	if(f == NULL)
	{
		printf("ERROR");
		return 0;
	}
	fMatlab = fopen("newMatrix.txt", "w");
	if(fMatlab == NULL)
	{
		printf("ERROR");
		return 0;
	}
	//Хорошо обусловленная матрица
	for(i = 0; i < SIZE; i++)
		for(j = 0; j < SIZE; j++)
			fscanf(f, "%lf", &A[i][j]);
	for(i = 0; i < SIZE; i++)
			fscanf(f, "%lf", &B[i]);

	printf("A:\n");
	PrintDoubleArray(A, SIZE, SIZE);
	printf("\nB:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", B[i]);

	LU(A, L, U, 10);

	X = Decision(L, U, B, SIZE);
	Nev = Nevyazka(A, B, X, SIZE);
	printf("\nNevyazka A B:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", Nev[i]);

	B1 = VozmushcheniyeB(B, SIZE);
	printf("\nB1:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", B1[i]);
	X1 = Decision(L, U, B1, SIZE);
	Nev = Nevyazka(A, B1, X1, SIZE);
	printf("\nNevyazka A B1:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", Nev[i]);
	norma = SearchK1(B, B1, X, X1, SIZE);

	A1 = VozmushcheniyeA(A, SIZE);
	printf("\nA1:\n");
	PrintDoubleArray(A1, SIZE, SIZE);
	LU(A1, L, U, 10);
	X2 = Decision(L, U, B, SIZE);
	Nev = Nevyazka(A1, B, X2, SIZE);
	printf("\nNevyazka A1 B:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", Nev[i]);
	norma = SearchK2(A, A1, X, X2, SIZE);

	Write(fMatlab, A, SIZE);
	Write(fMatlab, DiffMatr(A, A1, SIZE), SIZE);
	WriteStolb(fMatlab, B, SIZE);
	WriteStolb(fMatlab, Diff(B, B1, SIZE), SIZE);
	WriteStolb(fMatlab, X, SIZE);
	WriteStolb(fMatlab, Diff(X, X1, SIZE), SIZE);
	WriteStolb(fMatlab, Diff(X, X2, SIZE), SIZE);
	WriteStolb(fMatlab, Sum(X, Diff(X, X2, SIZE), SIZE), SIZE);
	

	//Плохо обусловленная матрица
	printf("\n\n\n/////////////////////////\n\n");
	for(i = 0; i < SIZE; i++)
		for(j = 0; j < SIZE; j++) 
			fscanf(f, "%lf", &A[i][j]);
	for(i = 0; i < SIZE; i++)
			fscanf(f, "%lf", &B[i]);

	printf("\nA:\n");
	PrintDoubleArray(A, SIZE, SIZE);
	printf("\nB:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", B[i]);

	LU(A, L, U, 10);

	X = Decision(L, U, B, SIZE);
	Nev = Nevyazka(A, B, X, SIZE);
	printf("\nNevyazka A B:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", Nev[i]);
	
	B1 = VozmushcheniyeB(B, SIZE);
	printf("\nB1:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", B1[i]);
	X1 = Decision(L, U, B1, SIZE);
	Nev = Nevyazka(A, B1, X1, SIZE);
	printf("\nNevyazka A B1:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", Nev[i]);
	norma = SearchK1(B, B1, X, X1, SIZE);

	A1 = VozmushcheniyeA(A, SIZE);
	printf("\nA1:\n");
	PrintDoubleArray(A1, SIZE, SIZE);
	LU(A1, L, U, 10);
	X2 = Decision(L, U, B, SIZE);
	Nev = Nevyazka(A1, B, X2, SIZE);
	printf("\nNevyazka A1 B:\n");
	for(i = 0; i < SIZE; i++)
			printf("%e\n ", Nev[i]);
	norma = SearchNK2(A, A1, X, X2, SIZE);
	
	Write(fMatlab, A, SIZE);
	Write(fMatlab, DiffMatrA1(A, A1, SIZE), SIZE);
	WriteStolb(fMatlab, B, SIZE);
	WriteStolb(fMatlab, Diff(B, B1, SIZE), SIZE);
	WriteStolb(fMatlab, X, SIZE);
	WriteStolb(fMatlab, Diff(X, X1, SIZE), SIZE);
	WriteStolb(fMatlab, Diff(X, X2, SIZE), SIZE);
	WriteStolb(fMatlab, Sum(X, Diff(X, X2, SIZE), SIZE), SIZE);

	FreeDoubleArray(A, SIZE);
	FreeDoubleArray(L, SIZE);
	FreeDoubleArray(U, SIZE);
	free(B);
	fclose (f);
	fclose (fMatlab);
	return 0;
}
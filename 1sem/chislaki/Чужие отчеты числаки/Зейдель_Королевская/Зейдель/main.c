#include <stdio.h>
#include <math.h>
#pragma warning (disable:4996)
#define M 3

double VectorNevjazki(double A[M][M], double x[M], double b[M], double d[M])
{
	int i, j;
	double sum = 0, max = 0;
	for (i = 0; i < M; i++)
	{
		sum = 0;
		for (j = 0; j < M; j++)
		{
			sum = sum + A[i][j] * x[j];
		}
		d[i] = sum - b[i];
		if (max < fabs(d[i]))
			max = fabs(d[i]);
	}
	return max;
}

double Norma(double alfa[M][M])
{
	int i, j;
	double norm=0,normmax=0;
	for (i = 0; i < M; i++)
	{
		norm = 0;
		for (j = 0; j < M; j++)
			norm = norm + fabs(alfa[i][j]);
		norm = norm;
		if (norm > normmax)
			normmax = norm;
	}
	return norm;
}

void Create(double A[M][M], double b[M], double alfa[M][M], double betta[M])
{
	int i, j;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
			alfa[i][j] = (-1)*A[i][j] / A[i][i];
		alfa[i][i] = 0;
		betta[i] = b[i] / A[i][i];
	}
}

int Zeidel (double alfa[M][M], double betta[M], double x1[M], int param) // 0 - норма < 1, 1 - At
{
	int i, j;
	int pogr = 0, itr=0;
	double X, q,norma;
	double x2[M];
	for (i = 0; i < M; i++)
		x2[i] = x1[i] = betta[i];
	norma = Norma(alfa);
	if (param = 1)
		norma = 0.5 ;
	q = fabs(norma / (1 - norma)) * 0.000001;
	printf("\nq = %lf\n", q);
	while (pogr != M)
	{
		pogr = 0;
		for (i = 0; i < M; i++)
		{
			X = 0;
			for (j = 0; j < M; j++)
			{
				X = X + alfa[i][j] * x1[j];
			}
			x1[i] = X + betta[i];
			if (fabs(x1[i] - x2[i]) < q)
				pogr++;
		}
		for (i = 0; i < M; i++)
			x2[i] = x1[i];
		itr++;
		for (i = 0; i < M; i++)
			printf("%lf ", x2[i]);
	}
	return itr;
}

void Trans(double A[M][M], double B[M][M])
{
	int i, j;
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			B[i][j] = A[j][i];
}

void Mulmatvec(double A[M][M], double B[M], double C[M])
{
	int i, j;
	for (i = 0; i < M; i++)
	{
		C[i] = 0;
		for (j = 0; j < M; j++)
		{

			C[i] += A[i][j] * B[j];
		}
	}
}
  
void Mulmatrix(double A[M][M], double B[M][M], double C[M][M])
{
	int i, j, k;
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
		{
			C[i][j] = 0;
			for (k = 0; k < M; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
}

double Det(double B[M][M])
{
	double A[M][M], det, p = 1;
	int i, j, k = 0, null = 1;
	//scanf("%d", &n);
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			A[i][j] = B[i][j];
	i = 0;
	/*	if(A[i][i] == 0)
		{
			null=1;
			while (null + i < n && A[i+null][i] == 0)
				null++;
			if((null+i)>=n)
			{
				printf("0");
				return 0;
			}
			for(j=i; j < n; j++)
			{
				det = A[i][j];
				A[i][j]=A[i+null][j];
				A[i+null][j]=det;
			}
			if (p==1)
				p=-1;
			else
				p=1;
		}
	for(j = 0; j < n - 1; j++)
	{
		for(k = 1; k < n; k++)
		{
			A[k][j+1] = A[k][j+1]*A[i][i] - A[k][i]*A[i][j+1];
		}
	}*/
	for (i = 0; i < M - 1; i++)
	{
		if (A[i][i] == 0)
		{
			null = 1;
			while (null + i < M && A[M + null][i] == 0)
				null++;
			if ((null + i) >= M)
			{
				printf("0");
				return 0;
			}
			for (j = i; j < M; j++)
			{
				det = A[i][j];
				A[i][j] = A[i + null][j];
				A[i + null][j] = det;
			}
			if (p == 1)
				p = -1;
			else
				p = 1;
		}
		for (j = i; j < M - 1; j++)
		{
			for (k = i + 1; k < M; k++)
			{
				A[k][j + 1] = (A[k][j + 1] * A[i][i] - A[k][i] * A[i][j + 1]);
				if (i)
					A[k][j + 1] /= A[i - 1][i - 1];
			}
		}
	}
	printf("%lf\n", p*A[M - 1][M - 1]);
	return p * A[M - 1][M - 1];
}



void ReverseMatrix(double A[M][M], double AA[M][M])
{
	int i, j;
	double det;
	double Atr[M][M];
	det = Det(A);
	Trans(A, Atr);
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			AA[i][j] = Atr[i][j] / det;
		}
	}
}

double Cond(double A[M][M])
{
	double AA[M][M];
	double nA, nAA;
	//int i, j;
	ReverseMatrix(A, AA);
	nA = Norma(A);
	nAA = Norma(AA);
	//for (i = 0; i < M; i++)
	//	for (j = 0; j < M; j++)
	//		printf("%lf ", AA[i][j]);
	//printf("Cond = %lf\n", (nA*nAA));
	return (nA*nAA);
}

int main(void)
{
	int i, j, itr,k,par;
	double norma;
	double A[M][M], alfa[M][M],At[M][M], AA[M][M];
	double b[M], betta[M], x1[M], bb[M], d[M];
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			scanf("%lf", &A[i][j]); 
	for (i = 0; i < M; i++)
	{
		scanf("%lf", &b[i]);
		x1[i] = b[i];
	}
	Create(A, b, alfa, betta);
	par = 0;
	printf("Cond(A) = %lf\n", Cond(A));
	if (Norma(alfa) > 1)
	{
		printf("Norma(alfa) > 1 \n");
		Trans(A, At);
		Mulmatvec(At, b, bb);
		Mulmatrix(At, A, AA);
		Create(AA, bb, alfa, betta);
		par = 1;
		printf("Cond(A*At) = %lf\n", Cond(AA));
	}
	printf("Cond(alfa) = %lf\n", Cond(alfa));
	for (i = 0; i < M; i++)
	{
		k = 0;
		printf("x(%d) = ", i);
		for (j = 0; j < M; j++)
		{
			printf("%lf*x(%d) ", alfa[i][j],k );
			k++;
		}
		printf("+%lf \n",betta[i]);
	}
	norma = Norma(alfa);
	printf("\n||alfa|| = %lf\n x = { ", norma);
	itr = Zeidel(alfa, betta, x1,par);
	for (i = 0; i < M; i++)
		printf("%lf ", x1[i]);
	printf("}\nITTERATION = %d", itr);
	printf("\n  Ax* - b = %e\n", VectorNevjazki(A, x1, b, d));
	scanf("%d", &i);
	return 0;
}
/*
10 1 1
2 10 1
2 2 10
12 13 14
*/

/* 0.0010 0.0010 -0.0009 -0.0009 -0.0002 -0.0002 -0.0012 0.0001 -0.0009 -0.0002 
0 0.0010 -0.0003 -0.0003 -0.0007 -0.0017 -0.0007 -0.0003 -0.0013 -0.0007 
0 0 0.0006 0.0004 -0.0006 0.0002 -0.0008 0.0004 -0.0007 -0.0004 
0 0 0.0002 0.0010 0.0002 0.0000 -0.0000 -0.0000 -0.0001 0.0004 
0 0 0.0003 0.0001 0.0014 0.0002 0.0001 0.0001 -0.0000 0.0003 
0 0 -0.0001 -0.0004 -0.0002 0.0006 0.0005 -0.0004 0.0006 -0.0003 
0 0 0.0005 0.0003 0.0007 0.0005 0.0014 0.0003 0.0002 -0.0004 
0 0 0.0004 0.0002 0.0005 0.0004 0.0003 0.0012 0.0002 -0.0003 
0 0 0.0004 0.0002 0.0005 0.0004 0.0003 0.0002 0.0012 -0.0003 
0 0 0.0005 0.0003 0.0007 0.0005 0.0004 0.0003 0.0002 0.0006 
*/
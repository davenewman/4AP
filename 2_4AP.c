

#include <stdio.h>
#include <stdlib.h>

#define MAXCHAR 1000

void solveTriDiag(double *a, double *b, double *c, double *f, double *u, int n);

void printFunction(int n, double *u, double *x);

void readConfig(int argc, char *argv[], int *numPoints, double *k, double *A, double *leftEnd, double *rightEnd, double *bc1, double *bc2);

int main(int argc, char *argv[])
{
//	double a[];
//	double b[];
//	double c[];
//	double f[];
//	double u[3];
	int numPoints;
	double k, A, leftEnd, rightEnd, bc1, bc2;
	readConfig(argc, argv, &numPoints, &k, &A, &leftEnd, &rightEnd, &bc1, &bc2);
	printf("numPoints = %d\n",numPoints);
	printf("k = %lf\n",k);
	printf("A = %lf\n",A);
	printf("leftEnd = %lf\n",leftEnd);
	printf("rightEnd = %lf\n",rightEnd);
	printf("bc1 = %lf\n",bc1);
	printf("bc2 = %lf\n",bc2);


	return 0;
}



void solveTriDiag(double *a, double *b, double *c, double *f, double*u, int n)
{
	//solves the tri-diagonal system. "a" is the main diagonal, "b" and "c" are the lower and upper diagonals, respectively.
	int i,k;

	for (i = 1; i < n; i++)
	{
		a[i] = a[i] - (1/a[i-1])*b[i-1]*c[i-1];
		f[i] = f[i] - (1/a[i-1])*b[i-1]*f[i-1];

	}

	u[n-1] = f[n-1]/a[n-1];

	for (k = n - 2; k > -1; k--)
	{
		u[k] = (1/a[k])*(f[k] - c[k]*u[k+1]);
	}
}


void printFunction(int n, double * u, double * x)
{
	int i;
	for (i = 0; i < n; i++)
	{
		printf("u[%d] = %.10lf\ta[%d] = %.10lf\n", i, u[i], i, x[i]);
	}

}

void readConfig(int argc, char *argv[], int *numPoints, double *k, double *A, double *leftEnd, double *rightEnd, double *bc1, double *bc2)

{
	*numPoints = atoi(argv[1]);
	*k = atof(argv[2]);
	*A = atof(argv[3]);
	*leftEnd = atof(argv[4]);
	*rightEnd = atof(argv[5]);
	*bc1 = atof(argv[6]);
	*bc2 = atof(argv[7]);
}

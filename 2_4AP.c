

#include <stdio.h>
void solveTriDiag(double *a, double *b, double *c, double *f, double*u, int n);

void printFunction(int n, double * u, double * x);

int main()
{
	double a[] = {2, 5, 1};
	double b[] = {2, 1};
	double c[] = {1, 2};
	double f[] = {1, 2, 3};
	double u[3];
	int n = 3;
	solveTriDiag(a, b, c, f, u, n);
	printFunction(n, u, a);


	return 0;
}



void solveTriDiag(double *a, double *b, double *c, double *f, double*u, int n)
{

	int i,k;

	for (i = 1; i < n; i++)
	{
		a[i] = a[i] - (1/a[i-1])*b[i-1]*c[i-1];//ok
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
		printf("u[%d] = %lf\ta[%d] = %lf\n", i, u[i], i, x[i]);
	}

}

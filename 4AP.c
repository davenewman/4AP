#include <stdio.h>
#include <math.h>

void createAlphaAndG(int nplus1, double * alpha, double * g, double A, double a, double b, double c, double u_0);

void printFunction(int nplus1, double * u, double * x);

void createU(int nplus1, double * alpha, double * u, double * g, double c);

void createX(int nplus1, double * x, double h);

int main()
{
	//the value for L will be inputted in the future using both endpoints of the domain
	double L = 1;

	//r = k^2
	double r = 100.0;
	double A = 1.0;

	//will cast this as a double and an int so that we can do math with it
	double n_f = 5.0;
	int n = 5;

	//BCs
	double leftBC = 1.0;
	double rightBC = 0.0;

	//calculate the step size
	double h = L/(n_f-1);

	//vector a has constant values
	double a = h*h*r - 2;

	//vectors b and c have constant values
	double b = 1.0;
	double c = 1.0;

	//initialize three vectors to hold alpha, g, and u
	double alpha[n];
	double g[n];
	double u[n + 1];

	//initialize vector to hold x (only for printing)
	double x[n + 1];

	u[0] = leftBC;
	u[nplus1] = rightBC;

	createAlphaAndG(n, alpha, g, A, a, b, c, u[0]);
	createU(nplus1, u, alpha, g, c);
	createX(nplus1, x, h);
	printFunction(n, u, x);


	return 0;
}

void createAlphaAndG(int n, double * alpha, double * g, double A, double a, double b, double c, double u_0)
{
	int j;
	alpha[0] = a;
	g[0] = A - u_0;

	for (j = 1; j < n - 1; j++)
	{
		alpha[j] = a - (b/alpha[j-1])*c;
		g[j] = A - (b/alpha[j-1])*g[j-1];
	}

}


void printFunction(int n, double * u, double * x)
{

	int i;

	for (i = 0; i < n; i++)
	{
		printf("u[%d] = %lf\tx[%d] = %lf\n", i, u[i], i, x[i]);
	}



}

void createU(int n, double * u, double * alpha, double * g, double c)
{

	int m;

	u[n-1] = g[n-1]/alpha[n-1];

	for (m = n - 2; m > 0; m--)
	{
		u[m] = (1/alpha[m])*(g[m] - c*u[m+1]);
	}


}

void createX(int n, double * x, double h)
{
	int p;

	for (p = 0; p < n; p++)
	{
		x[p] = (double)p*h;
	}

}

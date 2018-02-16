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
	double n_plus_1 = 100.0;
	int nplus1 = 100;

	//BCs
	double leftBC = 1.0;
	double rightBC = 0.0;

	//calculate the step size
	double h = L/n_plus_1;

	//vector a has constant values
	double a = h*h*r - 2;

	//vectors b and c have constant values
	double b = 1.0;
	double c = 1.0;

	//initialize three vectors to hold alpha, g, and u
	double alpha[nplus1];
	double g[nplus1];
	double u[nplus1 + 1];

	//initialize vector to hold x (only for printing)
	double x[nplus1 + 1];

	u[0] = leftBC;
	u[nplus1] = rightBC;

	createAlphaAndG(nplus1, alpha, g, A, a, b, c, u[0]);
	createU(nplus1, u, alpha, g, c);
	createX(nplus1, x, h);
	printFunction(nplus1, u, x);


	return 0;
}

void createAlphaAndG(int nplus1, double * alpha, double * g, double A, double a, double b, double c, double u_0)
{
	int j,k;
	alpha[0] = a;
	g[0] = A - u_0;

	for (j = 1; j < nplus1 - 1; j++)
	{
		alpha[j] = a - (b/alpha[j-1])*c;
		g[j] = A - (b/alpha[j-1])*g[j-1];
	}

}


void printFunction(int nplus1, double * u, double * x)
{

	int i;

	for (i = 0; i < nplus1; i++)
	{
		printf("u[%d] = %lf\tx[%d] = %lf\n", i, u[i], i, x[i]);
	}



}

void createU(int nplus1, double * u, double * alpha, double * g, double c)
{

	int m;

	u[nplus1-1] = g[nplus1-1]/alpha[nplus1-1];

	for (m = nplus1 - 2; m > 0; m--)
	{
		u[m] = (1/alpha[m])*(g[m] - c*u[m+1]);
	}


}

void createX(int nplus1, double * x, double h)
{
	int p;

	for (p = 0; p < nplus1; p++)
	{
		x[p] = (double)p*h;
	}

}

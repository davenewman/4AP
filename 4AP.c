#include <stdio.h>
#include <math.h>

void createAlphaAndG(int n, double * alpha, double * g, double A, double a, double b, double c, double u_0);

void printFunction(int n, double * u, double * x);

void createU(int n, double * alpha, double * u, double * g, double c);

void createX(int n, double * x, double h);

int main()
{
	//the value for L will be inputted in the future using both endpoints of the domain
	double L = 1;

	//r = k^2
	double r = 100.0;


	//will cast this as a double and an int so that we can do math with it
	double n_f = 200.0;
	int n = 200;

	//BCs
	double leftBC = 1.0;
	double rightBC = 0.0;

	//calculate the step size
	double h = L/(n_f-1);
	double A = h*h*1.0;

	//vector a has constant values
	double a = h*h*r - 2;

	//vectors b and c have constant values
	double b = 1.0;
	double c = 1.0;

	//initialize three vectors to hold alpha, g, and u
	double alpha[n-1];
	double g[n-1];
	double u[n];

	//initialize vector to hold x (only for printing)
	double x[n];

	u[0] = leftBC;
	u[n-1] = rightBC;

	//testing
	printf("u0 = %lf\nu_end = %lf\n\n",u[0],u[n-1]);

	createAlphaAndG(n, alpha, g, A, a, b, c, u[0]);
	createU(n, u, alpha, g, c);
	createX(n, x, h);
	printFunction(n, u, x);
	//testing
//	printf("u[99] = %lf ok\n",u[99]);

	return 0;
}

void createAlphaAndG(int n, double * alpha, double * g, double A, double a, double b, double c, double u_0)
{
	int j;
	alpha[0] = a;
	g[0] = A - u_0;
	alpha[1] = a - (b/alpha[0])*c;
	g[1] = g[0] - (b/a)*g[0];

	for (j = 2; j < n - 1; j++)
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

	u[n-2] = g[n-2]/alpha[n-2];

	for (m = n - 3; m > 0; m--)
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

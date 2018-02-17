

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCHAR 1000

void solveTriDiag(double *a, double *b, double *c, double *f, double *u, int n);

void createVectors(double *a, double *b, double *c, double *f, int n, double A);

void printFunction(int n, double * u, int type, double bc1, double bc2, double leftEnd, double rightEnd, double h);

int readConfig(int argc, char *argv[], int *numPoints, double *k, double *A, double *leftEnd, double *rightEnd, double *bc1, double *bc2);

int main(int argc, char *argv[])
{
	
	int numPoints,n;
	double k, A, leftEnd, rightEnd, bc1, bc2, firstC, lastB, firstF, lastF, h, diagElement;
	int type = readConfig(argc, argv, &numPoints, &k, &A, &leftEnd, &rightEnd, &bc1, &bc2);
	
	h = (rightEnd - leftEnd)/( (double)numPoints - 1.0);						//step size h
	//check the type of boundary conditions, will affect what we put in to the build vectors function
	//n is the number of elements in the main diagonal
	//this could probably be reformatted to be shorter, but for my sake I need the readability
	 
	if (type == 0) {															// Both Dirichlet
		n = numPoints - 2;
		firstC = 1.0;
		lastB = 1.0;
		firstF = A*h*h - bc1;
		lastF = A*h*h - bc2; }
	else if (type ==1) {														// Neumann left, Dirichlet right
		n = numPoints - 1;
		firstC = 2.0;
		lastB = 1.0;
		firstF = A*h*h + 2*h*bc1;
		lastF = A*h*h - bc2; }
	else if (type == 2) {														// Dirichlet left, Neumann right
		n = numPoints - 1;
		firstC = 1.0;
		lastB = 2.0;
		firstF = A*h*h - bc1;
		lastF = A*h*h - 2*h*bc2; }//double check?
	else {																		// Both Neumann 
		n = numPoints;
		firstC = 2.0;
		lastB = 2.0;
		firstF = A*h*h + 2*h*bc1;
		lastF = A*h*h - 2*h*bc2; }
		
	double a[n], b[n-1], c[n-1], f[n], u[n];
	
	//now we need to fill these vectors with some function
	a[0] = -1.0*k*h*h - 2;
	a[n-1] = a[0];
	b[n-2] = lastB;
	c[0] = firstC;
	f[0] = firstF;
	f[n-1] = lastF;
	A = A*h*h;
	createVectors(a, b, c, f, n, A);
	for (int y = 0; y < sizeof(a)/sizeof(double); y++)
		printf("a[%d] = %.10lf\tf[%d] = %.10lf\n",y,a[y],y,f[y]);
	
	printf("\n");
	
	for (int z = 0; z < sizeof(b)/sizeof(double); z++)
		printf("b[%d] = %.10lf\tc[%d] = %.10lf\n",z,b[z],z,c[z]);
	
	printf("\nnumpoints = %d\n", numPoints);
	printf("h = %lf\n", h);
	printf("type = %d\n",type);
	
	
	solveTriDiag(a, b, c, f, u, n);
	
	printf("number of elements in a = %lu\n",sizeof(a)/sizeof(double));
	printf("number of elements in b = %lu\n",sizeof(b)/sizeof(double));
	printf("number of elements in c = %lu\n",sizeof(c)/sizeof(double));
	printf("number of elements in f = %lu\n",sizeof(f)/sizeof(double));
	printf("number of elements in u = %lu\n",sizeof(u)/sizeof(double));
	printf("\nNew A values:\n");
	for (int s = 0; s < sizeof(a)/sizeof(double); s++)
		printf("a[%d] = %lf\tf[%d] = %lf\n",s,a[s],s,f[s]);
		
	printf("\n");

	printFunction(n, u, type, bc1, bc2, leftEnd, rightEnd, h);
	return 0;
}

void solveTriDiag(double *a, double *b, double *c, double *f, double*u, int n)	//solves the tri-diagonal system. 
																				//"a" is the main diagonal, 																				//"b" and "c" are the lower and upper diagonals, respectively.
{
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

void printFunction(int n, double * u, int type, double bc1, double bc2, double leftEnd, double rightEnd, double h)
{
	int numtoPrint, i;
	if (type == 0)
	{
		numtoPrint = n + 2;
		printf("%.10lf\t%.10lf\n",leftEnd,bc1);
		for (i = 1; i < numtoPrint-1; i++)
			printf("%.10lf\t%.10lf\n",leftEnd+((double)i)*h,u[i-1]);
		printf("%.10lf\t%.10lf\n",rightEnd,bc2);
	}
	else if (type == 1)
	{
		numtoPrint = n + 1;
		for (i = 0; i < numtoPrint-1; i++)
			printf("%.10lf\t%.10lf\n",leftEnd+((double)i)*h,u[i]);
		printf("%.10lf\t%.10lf\n",rightEnd,bc2);
	}
	else if (type == 2)
	{
		numtoPrint = n + 1;
		printf("%.10lf\t%.10lf\n",leftEnd,bc1);
		for (i = 1; i < numtoPrint; i++)
			printf("%.10lf\t%.10lf\n",leftEnd+((double)i)*h,u[i]);
	}
	else
	{
		numtoPrint = n;
		for(i = 0; i < numtoPrint; i++)
			printf("%.10lf\t%.10lf\n",leftEnd+((double)i)*h,u[i]);
	}
}

int readConfig(int argc, char *argv[], int *numPoints, double *k, double *A, double *leftEnd, double *rightEnd, double *bc1, double *bc2)

{
	//returns 0 for two Dirichlet conditions, 1 for Neumann left and Dirichlet right, 2 for Dirichlet left and Neumann right
	//and 3 for two Neumann conditions
	*numPoints = atoi(argv[1]);
	*k = atof(argv[2]);
	*A = atof(argv[3]);
	*leftEnd = atof(argv[4]);
	*rightEnd = atof(argv[5]);
	
	char *temp1 = argv[6];
	char *temp2 = argv[7];
	if (temp1[0] == '*' & temp2[0] == '*')
	{
		memmove(temp1, temp1+1,strlen(temp1));
		memmove(temp2, temp2+1,strlen(temp2));
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 3;
	}
	else if (temp1[0] == '*' & temp2[0] != '*')
	{
		memmove(temp1, temp1+1, strlen(temp1));
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 1;
	}
	else if (temp1[0] != '*' & temp2[0] == '*')
	{
		memmove(temp2, temp2+1, strlen(temp2));
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 2;
	}
	else 
	{
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 0;
	}
}


void createVectors(double *a, double *b, double *c, double *f, int n, double A)
{
	//creates all vectors, waste of memory since most values will be 1
	int p, k;
	for (p = 1; p < n; p++)
	{
		a[p] = a[0];
		f[p] = A;
	}
	
	for (k = 0; k < n - 2; k++)
	{
		b[k] = 1.0;
		c[k+1] = 1.0;
	}
}

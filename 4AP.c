

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void solveTriDiag(double *a, double *b, double *c, double *f, double *u, int n);
void createVectors(double *a, double *b, double *c, double *f, int n, double A);
void printFunction(int n, int numPoints, double * u, int type, double bc1, double bc2, double leftEnd, double rightEnd, double h, double k, double A);
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
	 
	if (type == 0) {				// Both Dirichlet
		n = numPoints - 2;
		firstC = 1.0;
		lastB = 1.0;
		firstF = A*h*h - bc1;
		lastF = A*h*h - bc2; }
	else if (type ==1) {			// Neumann left, Dirichlet right
		n = numPoints - 1;
		firstC = 2.0;
		lastB = 1.0;
		firstF = A*h*h + 2*h*bc1;
		lastF = A*h*h - bc2; }
	else if (type == 2) {			// Dirichlet left, Neumann right
		n = numPoints - 1;
		firstC = 1.0;
		lastB = 2.0;
		firstF = A*h*h - bc1;
		lastF = A*h*h - 2*h*bc2; }
	else {							// Both Neumann 
		n = numPoints;
		firstC = 2.0;
		lastB = 2.0;
		firstF = A*h*h + 2*h*bc1;
		lastF = A*h*h - 2*h*bc2; }
		
	double a[n], b[n-1], c[n-1], f[n], u[n];
	
	a[0] = -1.0*k*h*h - 2;
	a[n-1] = a[0];
	b[n-2] = lastB;
	c[0] = firstC;
	f[0] = firstF;
	f[n-1] = lastF;
	createVectors(a, b, c, f, n, A*h*h);
	solveTriDiag(a, b, c, f, u, n);
	printFunction(n, numPoints, u, type, bc1, bc2, leftEnd, rightEnd, h, k, A);
	return 0;
}

void solveTriDiag(double *a, double *b, double *c, double *f, double*u, int n)
{
	// "a": main diagonal | "b": lower diagonal | "c": upper diagonal
	int i,k;
	for (i = 1; i < n; i++) {
		a[i] = a[i] - (1/a[i-1])*b[i-1]*c[i-1];
		f[i] = f[i] - (1/a[i-1])*b[i-1]*f[i-1];
	}
	u[n-1] = f[n-1]/a[n-1];
	for (k = n - 2; k > -1; k--) {
		u[k] = (1/a[k])*(f[k] - c[k]*u[k+1]);
	}
}

void printFunction(int n, int numPoints, double * u, int type, double bc1, double bc2, double leftEnd, double rightEnd, double h, double k, double A)
{
	// Dynamically print to filename based on values of command line arguments
	const char* s1 = "t_";
	const char* s2 = "_n_";
	const char* s3 = "_k2_";
	const char* s4 = "_A_";
	const char* s5 = "_L_";
	const char* s6 = "_R_";
	const char* s7 = ".txt";
	char filename[128];
	FILE* fp = NULL;
	sprintf(filename,"%s%d%s%d%s%.1f%s%.1f%s%.1f%s%.1f%s",s1,type,s2,numPoints,s3,k,s4,A,s5,bc1,s6,bc2,s7);
	printf("Saved to %s\n",filename);
	fp = fopen(filename,"w");
	
	int numtoPrint, i;
	if (type == 0) {
		numtoPrint = n + 2;
		fprintf(fp,"%.15lf\t%.15lf\n",leftEnd,bc1);
		for (i = 1; i < numtoPrint-1; i++)
			fprintf(fp,"%.15lf\t%.15lf\n",leftEnd+((double)i)*h,u[i-1]);
		fprintf(fp,"%.15lf\t%.15lf\n",rightEnd,bc2);
	}
	else if (type == 1) {
		numtoPrint = n + 1;
		for (i = 0; i < numtoPrint-1; i++)
			fprintf(fp,"%.15lf\t%.15lf\n",leftEnd+((double)i)*h,u[i]);
		fprintf(fp,"%.15lf\t%.15lf\n",rightEnd,bc2);
	}
	else if (type == 2) {
		numtoPrint = n + 1;
		fprintf(fp,"%.15lf\t%.15lf\n",leftEnd,bc1);
		for (i = 1; i < numtoPrint; i++)
			fprintf(fp,"%.15lf\t%.15lf\n",leftEnd+((double)i)*h,u[i]);
	}
	else {
		numtoPrint = n;
		for(i = 0; i < numtoPrint; i++)
			fprintf(fp,"%.15lf\t%.15lf\n",leftEnd+((double)i)*h,u[i]);
	}
}

int readConfig(int argc, char *argv[], int *numPoints, double *k, double *A, double *leftEnd, double *rightEnd, double *bc1, double *bc2)
{
	// parses the command line arguments passed to the program
	// returns 0: both Dirichlet | 1: Neumann left Dirichlet right | 2: Dirichlet left Neumann right | 3: both Neumann
	*numPoints = atoi(argv[1]);
	*k = atof(argv[2]);
	*A = atof(argv[3]);
	*leftEnd = atof(argv[4]);
	*rightEnd = atof(argv[5]);
	
	char *temp1 = argv[6];
	char *temp2 = argv[7];
	if (temp1[0] == '*' & temp2[0] == '*') {
		memmove(temp1, temp1+1,strlen(temp1));
		memmove(temp2, temp2+1,strlen(temp2));
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 3;
	}
	else if (temp1[0] == '*' & temp2[0] != '*') {
		memmove(temp1, temp1+1, strlen(temp1));
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 1;
	}
	else if (temp1[0] != '*' & temp2[0] == '*') {
		memmove(temp2, temp2+1, strlen(temp2));
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 2;
	}
	else {
		*bc1 = atof(temp1);
		*bc2 = atof(temp2);
		return 0;
	}
}


void createVectors(double *a, double *b, double *c, double *f, int n, double A)
{
	//creates all vectors, waste of memory since most values will be 1
	//however the solve function expects vectors
	int p, k;
	for (p = 1; p < n; p++) {
		a[p] = a[0];
		f[p] = A;
	}
	
	for (k = 0; k < n - 2; k++) {
		b[k] = 1.0;
		c[k+1] = 1.0;
	}
}

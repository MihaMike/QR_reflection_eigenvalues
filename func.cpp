#include "func.h"

#define A(i,j) a[(i)*n+(j)]

#define ds sizeof(real)

const real eps = 1e-50;

void PrintVector(real * arr, int n);

real RowNorm(real * a, int n)
{
    real N = 0; for (int i = 0; i < n;i++)
    { real xx = 0; for (int j = 0; j < n; j++) xx += A(i,j)*A(i,j); xx = sqrt(xx); if (xx > N) N = xx; }
    return N;
}

bool QRReflect(real * a, int n, int e, real * q)
{
    real x0,x1;
    real Nx, Na; // Norm of vector
    real s, dp; // k-th el. square and dot product

    for (int i = 0; i < e; i++)
	for (int j = 0; j < e; j++) q[i*n+j] = (i == j);
    
    for (int k = 0; k < e-1; k++) 
    {
	x0 = A(k,k);
	x1 = A(k+1,k);
	
	s = x1*x1;
	Na = sqrt(x0*x0+s);
	x0 -= Na;
	Nx = sqrt(x0*x0+s);
	if (Nx < eps) return false;
	Nx = 1.0 / Nx;
	    
	x0 *= Nx;
	x1 *= Nx;
	
        for (int i = k+1; i < e; i++) // not optimal yet!!!
	{
	    real dp = 2.0*(A(k,i)*x0+A(k+1,i)*x1);
    	    A(k,i) -= x0*dp;
    	    A(k+1,i) -= x1*dp;
	}
	A(k,k) = Na;
	A(k+1,k)=0.0;
	
	for (int i = 0; i < e; i++) // row cycle
	{
	    dp = 2.0*(x0*q[i*n+k] + x1*q[i*n+k+1]);
	    q[i*n+k] -= x0*dp;
	    q[i*n+k+1] -= x1*dp;
	}
    }
    return true;
}

void Multiply(real * r, real * q, real * res, int n, int e)
{
    for (int i = 0; i < e; i++)
	for (int j = 0; j < e; j++)
	{
	    real s = 0;
	    for (int k = 0; k < e; k++) s += r[i*n+k]*q[k*n+j];
	    res[i*n+j] = s;
	}
}

// =============================================
bool S_Reflect(real b, real * a, int n, real * q, real * ev)
{
    real Row_norm = RowNorm(a,n);
    real * res = new real[n*n];
    q[0] = ev[0]=0;
    real * x = new real[n];   // vector X
    real s, Nx, Na; // Norm of vector

    for (int k = 0; k < n-2; k++) // steps. let's begin...
    {
	// calc norm of first column (firsk means k-th !!!)
	s = 0.0;
	// Get vector a_1 
	for (int i=k+2;i<n;i++) { s+=A(i,k)*A(i,k); x[i] = A(i,k); }
	x[k+1] = A(k+1,k);

	bool flag = false;
	for (int i=k+2; i < n; i++) if (fabs(x[i]) > eps) flag = true;
	if (!flag) continue;
	Na = sqrt(s + A(k+1,k)*A(k+1,k));
	x[k+1] -= Na; // subtract (e_k * norm) -- in fact, subtr one coordinate

	// normalize x: get its norm and divide by it
	Nx = sqrt(s + x[k+1]*x[k+1]);
	if (Nx < eps) return false;
	Nx = 1.0 / Nx;
	    
	for (int i=k+1;i<n;i++) x[i] *= Nx;
	
	for (int col = k + 1; col < n; col++) // cycle through cols of A
	{
    	    // calc dot prod. of x and y
	    real dp = 0.0;
	    for (int i = k+1; i < n; i++) dp += x[i]*A(i,col);
	    dp *= 2.0;
	    for (int i = k+1; i < n; i++) A(i,col) -= x[i]*dp;
	}
	A(k+1,k) = Na;
	for (int i = k+2; i < n; i++) A(i,k) = 0.0;
	for (int row = k + 1; row < n; row++) // cycle through cols of A
	{
	    // calc dot prod. of x and y
	    real dp = 0.0;
	    for (int i = k+1; i < n; i++) dp += x[i]*A(row,i);
	    dp *= 2.0;
	    for (int i = k+1; i < n; i++) A(row,i) -= x[i]*dp;
	}
	A(k,k+1) = Na;
	for (int j = k+2; j < n; j++) A(k,j) = 0.0;
    }
    delete x;
    for (int k = n; k > 2; k--)
    {
//	int it = 0;
	printf("Step: k = %d\n", k);
	while (fabs(A(k-1,k-2)) > b * Row_norm)
	{
    	    QRReflect(a, n, k, q);
//	    printf("a (R), k = %d:\n", k); PrintMatrix(a, 0, n);

//	    printf("q, k = %d:\n", k); PrintMatrix(q, 0, n);
	    
	    Multiply(a, q, res, n, k); 
	    
	    for (int i = 0; i < n*n; i++) a[i] = res[i];
//	    printf("a: \n"); PrintMatrix(a, 0, n);
//	    it++; if (it > 10) { printf("Error: too many iterations.\n"); exit(1); }
	}
    }
    
    if ( n > 1)
    {
	real t = a[0]+a[n+1];
	real d = a[0]*a[n+1]-a[1]*a[n];
	d = t*t-4*d;
	if (d < 0) printf("suxx: inconsistent equation!!!\n");
	else
	{
	    d= sqrt(d);
	    a[0] = (t+d)/2.0;
	    a[n+1] = (t-d)/2.0;
	}	
    }
//    PrintMatrix(a,0,n);
    for (int i = 0; i < n; i++)ev[i] = A(i,i);
    delete res;
    return true;
}

void PrintMatrix(real * arr, int k, int n)
{
    for (int i = k; i < n; i++)
    {
	for (int j = k; j < n; j++) printf("%1.7Lf ", arr[i*n+j]);
	printf("\n");
    }
}

void PrintVector(real * arr, int n)
{
    printf("arr = ( "); for (int i = 0; i < n; i++) printf("%1.3Lf ", arr[i]); printf(")\n");
}
    
bool ReadMatrix(char * fn, real * array, int m, int n)
{    
    FILE * F = fopen(fn, "r");
    if (!F) return false;
    int c = 0;
    while (!feof(F))
    {
	real x;
	int r = fscanf(F, "%Lf", &x);
	if (r != 1) { fclose(F); return false; }
	array[c] = x;
	c++;
    }
    if (c < n*m) return false;
    fclose(F);
    return true;
}

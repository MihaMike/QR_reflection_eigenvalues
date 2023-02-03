#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "func.h"

real f1(int k, int l)
{
    int i = k+1;
    int j = l+1;
    real A = 10.0;
    real Pi = 4.0*atan(1.0);
    
    real d = real(i-j);
    real s = real(i+j);
    real p = 2*Pi*i;
    real q = Pi*i/A;
    if (i!=j) return 2*A*A/(Pi*Pi) * (1.0/(d*d)-1.0/(s*s));
    else return A*A/3.0-2.0*A*A/(p*p)+(q*q);
}

real f2(int k, int l)
{
    int i = k+1;
    int j = l+1;
    real A = 7.0;
    real Pi = 4.0*atan(1.0);
    
    real d = real(i-j);
    real s = real(i+j);
    real p = Pi*(2*i-1);
    real q = (k-0.5)*Pi/A;
    if (i!=j) return 2*A*A/(Pi*Pi) * (1.0/(d*d)-1.0/(s*s-1));
    else return A*A/3.0-2.0*A*A/(p*p)+(q*q);
}

real f3(int k, int l)
{
    int i = k+1;
    int j = l+1;
    real A = 7.0;
    real Pi = 4.0*atan(1.0);
    
    real d = real(i-j);
    real s = real(i+j);
    real p = Pi*(2*i-1);
    real q = Pi*i/A;
    real c = pow(-1, k+l);
    if (i!=j) return c*2*A*A/(Pi*Pi) * (1.0/(d*d)-1.0/(s*s-1));
    else return A*A/3.0-2.0*A*A/(p*p)+(q*q);
}

real f4(int k, int l)
{
    int i = k+1;
    int j = l+1;
    real A = 7.0;
    real Pi = 4.0*atan(1.0);
    
    real d = real(i-j);
    real s = real(i+j);
    real p = Pi*(2*i-1);
    real q = (k-0.5)*Pi/A;
    real c= pow(-1,k+l);
    if (i!=j) return c*2*A*A/(Pi*Pi) * (1.0/(d*d)-1.0/(s*s-1));
    else return A*A/3.0-2.0*A*A/(p*p)+(q*q);
}

real CalcError(real * a, real * ev, int n)
{
    real tr = 0;
    real trd = 0;
    for (int i = 0; i < n; i++)
    {
	tr += a[i*n+i];
	trd += ev[i];
    }
    return fabs(tr -trd);
}

void FillMatrix(real * a, int n)
{
    for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++) a[i*n+j] = f1(i,j);
}

int main()
{
    try
    {
	int n;
        char fn[256];
        printf("Input dimension (n): ");
        scanf("%d", &n);
	if (n <= 0) throw -10;
//	printf("Input source file name (or 'func' to use function): ");
//	scanf("%255s", fn);
	real * a = new real[n*n];
	real * ac = new real[n*n];
	real * q = new real[n*n];
	real * ev = new real[n]; // eigenvalues
	//if (fn[0] == 'f' && fn[1] == 0)
	FillMatrix(a, n);
	//else ReadMatrix(fn, a, n, n);
	memcpy(ac, a, sizeof(real)*n*n);
//	PrintMatrix(a,0,n);
	printf("Input accuracy: ");
	real DEPbMO;
	scanf("%Lf", &DEPbMO);
	printf("let's begin!..");
	clock_t ts = clock(); 
	if (!S_Reflect(DEPbMO, a, n, q, ev)) { printf("Bad matrix!\n"); throw -100; }
	clock_t te = clock();
	printf("done.\n Eigenvalues of matrix are:\n");
	for (int i = 0; i < n; i++)
	{
	    printf("%1.10Lf ", ev[i]);    
	}
	printf("\n");
	
	printf("Elapsed time: %.3Lf\n", real(te-ts)/real(CLOCKS_PER_SEC));
	
	printf("Error: %1.15Lf\n", CalcError(ac, ev, n));
	delete a;
	delete ac;
        delete q;
        throw 0;
    }
    catch(int err)
    {
	switch (err)
	{
	    case 0: return 0;
	    default: printf("Error in main program. err = %d\n", err); return err; 
	}
    }
}

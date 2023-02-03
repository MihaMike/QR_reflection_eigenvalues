#ifndef ___XPEH_BAM_BCEM___
#define ___XPEH_BAM_BCEM___

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

typedef long double real;
// read matrix from file FileName
bool ReadMatrix(char * FileName, real * a, int m, int n);
// debug method: display matrix
void PrintMatrix(real * arr, int k, int n);

// Solve system using Reflection method
bool S_Reflect(real XPEHOTEHb, real * arr, int n, real * q, real * ev);

#endif


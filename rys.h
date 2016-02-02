#ifndef _RYS_H_
#define _RYS_H_
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "util.h"
#define MAXROOTS 11
#define MAXRP1 12
#define MAXRP2 13
#define CSIZE 144
#define TSIZE 288
#define RSIZE 121
#define FSIZE 23
#define THRESHOLD 1.e-14

/* rys.c */
extern void Rys_setB00 ( double x );
extern void Rys_setB1 ( double x );
extern void Rys_setB1p ( double x );
extern void Rys_setC ( double x );
extern void Rys_setCp ( double x );
extern const double *RysWeights ();
extern const double *RysRoots ();
extern void Rys_findRoots ( int n, double X );
extern void RysInit ();
extern void RysRecur ( double G[][11], int l12, int l34 );
extern double RysShift ( double abx, double cdx, double G[][11], int l12, int l2,
                         int l34, int l4 );
#endif

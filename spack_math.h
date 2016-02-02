#ifndef _SPACKMATH_H_
#define _SPACKMATH_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

/* spack_math.h */
extern void copy_identity ( const int nr, const int nc,
                            double *z, const int ldz );

extern void copy_trans ( const int n, double *z );
extern void spack_diag ( const int n, double *a, double *z,
                         double *evals, double *tmp );

extern void spack_trans ( int n, double *a, double *z, double *tmp );

extern double spack_trace_prod ( const int n, const double *a,
                                 const double *b );
#endif

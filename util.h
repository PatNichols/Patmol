#ifndef _UTIL_H_
#define _UTIL_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* putil.c */
extern void exit_program();
extern void fatalError ( const char *message );
extern void *pmalloc ( size_t nbytes );
extern void *pcalloc ( size_t nmem, size_t memsize );
extern void *prealloc ( void *old_ptr, size_t nbytes );
extern FILE *openFile ( const char *fname, const char *mode );
extern double *dalloc ( size_t n );
extern double **dmatrix ( size_t n1, size_t n2 );
extern double ***dtensor ( size_t n1, size_t n2, size_t n3 );
extern void free_dmatrix ( double **m, size_t n1 );
extern void free_dtensor ( double ***t, size_t n1, size_t n2 );
extern int *ialloc ( size_t n );
extern int **imatrix ( size_t n1, size_t n2 );
extern int ***itensor ( size_t n1, size_t n2, size_t n3 );
extern void free_imatrix ( int **m, size_t n1 );
extern void free_itensor ( int ***t, size_t n1, size_t n2 );
extern int feqZero ( const double x );
extern double distSqr ( const double *a, const double *b );
extern double normSqr ( const double *a );
extern int fcmp ( const double x, const double y );
#endif

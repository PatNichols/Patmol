
#include "util.h"

#ifdef USE_MPI
#include <mpi.h>

inline void exit_program()
{
    MPI_Abort ( MPI_COMM_WORLD, EXIT_FAILURE );
}

#else

inline void exit_program()
{
    exit ( EXIT_FAILURE );
}

#endif

inline void
fatalError ( const char *message )
{
    fprintf ( stderr, "FATAL ERROR: %s\n", message );
    exit ( EXIT_FAILURE );
}

inline void *
pmalloc ( size_t nbytes )
{
    void *p = malloc ( nbytes );
    if ( !p ) {
        fprintf ( stderr, "Could not allocate %s bytes\n", nbytes );
        exit_program();
    }
    return p;
}

inline void *
pcalloc ( size_t nmem, size_t memsize )
{
    void *p = calloc ( nmem, memsize );
    if ( !p ) {
        fprintf ( stderr, "Could not allocate %d objects of size %d bytes\n", nmem, memsize );
        exit_program();
    }
    return p;
}

inline void *
prealloc ( void *old_ptr, size_t nbytes )
{
    void *p = realloc ( old_ptr, nbytes );
    if ( !p ) {
        fprintf ( stderr, "Could not reallocate %s bytes\n", nbytes );
        exit_program();
    }
    return p;
}

inline FILE * openFile ( const char *fname, const char *mode )
{
    FILE *fp = fopen ( fname, mode );
    if ( !fp ) {
        fprintf ( stderr, "Could not open the file %s in mode %s!\n", fname, mode );
        exit_program();
    }
    return fp;
}

inline double *
dalloc ( size_t n )
{
    return ( double * ) pmalloc ( sizeof ( double ) * n );
}

inline double **
dmatrix ( size_t n1, size_t n2 )
{
    int i;
    double **m = ( double ** ) pmalloc ( sizeof ( double * ) * n1 );
    for ( i = 0; i < n1; ++i ) {
        m[i] = ( double * ) pmalloc ( sizeof ( double ) * n2 );
    }
    return m;
}

inline double ***
dtensor ( size_t n1, size_t n2, size_t n3 )
{
    int i, j;
    double ***t = ( double *** ) pmalloc ( sizeof ( double ** ) * n1 );
    for ( i = 0; i < n1; ++i ) {
        t[i] = ( double ** ) pmalloc ( sizeof ( double * ) * n2 );
        for ( j = 0; j < n2; ++j ) {
            t[i][j] = ( double * ) pmalloc ( sizeof ( double ) * n3 );
        }
    }
    return t;
}

inline void
free_dmatrix ( double **m, size_t n1 )
{
    register int i = n1;
    for ( ; i--; )
        free ( m[i] );
    free ( m );
}

inline void
free_dtensor ( double ***t, size_t n1, size_t n2 )
{
    register int i, j;
    for ( i = n1; i--; ) {
        for ( j = n2; j--; ) {
            free ( t[i][j] );
        }
        free ( t[i] );
    }
    free ( t );
}

inline int *
ialloc ( size_t n )
{
    return ( int * ) pmalloc ( sizeof ( int ) * n );
}

inline int **
imatrix ( size_t n1, size_t n2 )
{
    int i;
    int **m = ( int ** ) pmalloc ( sizeof ( int * ) * n1 );
    for ( i = 0; i < n1; ++i ) {
        m[i] = ( int * ) pmalloc ( sizeof ( int ) * n2 );
    }
    return m;
}

inline int ***
itensor ( size_t n1, size_t n2, size_t n3 )
{
    int i, j;
    int ***t = ( int *** ) pmalloc ( sizeof ( int ** ) * n1 );
    for ( i = 0; i < n1; ++i ) {
        t[i] = ( int ** ) pmalloc ( sizeof ( int * ) * n2 );
        for ( j = 0; j < n2; ++j ) {
            t[i][j] = ( int * ) pmalloc ( sizeof ( int ) * n3 );
        }
    }
    return t;
}

inline void
free_imatrix ( int **m, size_t n1 )
{
    register int i = n1;
    for ( ; i--; )
        free ( m[i] );
    free ( m );
}

inline void
free_itensor ( int ***t, size_t n1, size_t n2 )
{
    register int i, j;
    for ( i = n1; i--; ) {
        for ( j = n2; j--; ) {
            free ( t[i][j] );
        }
        free ( t[i] );
    }
    free ( t );
}

/** extra math functions **/
inline int
feqZero ( const double x )
{
    if ( fabs ( x ) <= DBL_EPSILON )
        return 1;
    return 0;
}

inline double
distSqr ( const double *a, const double *b )
{
    register double x, y, z;
    x = a[0] - b[0];
    y = a[1] - b[1];
    z = a[2] - b[2];
    return ( x * x + y * y + z * z );
}

inline double
normSqr ( const double *a )
{
    return ( a[0] * a[0] + a[1] * a[1] + a[2] * a[2] );
}

inline int
fcmp ( const double x, const double y )
{
    register double z = x - y;
    if ( fabs ( z ) <= DBL_EPSILON )
        return 0;
    if ( z > 0. )
        return 1;
    return -1;
}



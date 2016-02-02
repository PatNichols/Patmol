#ifndef _ERIFS_H_
#define _ERIFS_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "util.h"
#include "pconfig.h"

typedef struct {
    int rank, nproc, nfiles;
    int *numints;
    char *basename;
} ERIFS;

typedef struct {
    double val;
    unsigned short i, j, k, l;
} ERINT __attribute__((aligned(16)));

/* erifs.c */
extern FILE *openERIFile ( int rank, int file_no, const char *bname );
extern FILE *createERIFile ( int rank, int file_no, const char *bname );
extern ERIFS *init_ERIFS ( int nst, int nend, const char *sname );
extern void clean_ERIFS ( ERIFS * e );
#endif

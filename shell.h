#ifndef _SHELL_H_
#define _SHELL_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "util.h"

typedef struct {
    int npr, lsh, cen, off, nst;
    double *al;
    double *co;
} Shell;

/* shell.c */
extern void read_shell ( FILE * infile, Shell * sh, int *off );
extern void normalize_shell ( Shell * sh );
extern void write_shell ( FILE * outfile, Shell * sh );
extern int shell_eq ( const Shell * s1, const Shell * s2 );

#endif

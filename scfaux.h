#ifndef _SCFAUX_H_
#define _SCFAUX_H_
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "basis.h"
#include "spack_math.h"
#include "erifs.h"
/* scfaux.c */
extern void form_cmatrix ( int n, double *cm, const double *xm, double *wm );
extern void form_pmatrix ( int no, int noc, const double *cm, double *pm );
extern void form_pmatrix2 ( int no, int *noc, const double *cm, double *pm );
extern void scfAccelerator ( int iter, int no2,
                             double *p0, double *p1, double *p2 );
extern double calc_energy ( int no, const double *pm,
                            const double *hm, const double *gm );
extern double pdiffNorm ( int no, const double *p0, const double *p1 );
extern void rhf_analyzeMoments ( int no2, const double *pmat, double *cmat );
extern void uhf_analyzeMoments ( int no2, const double *pmatA,
                                 const double *pmatB, double *cmatA,
                                 double *cmatB );
extern void rhf_form_gmatrix ( const ERIFS * tfiles, const double *Pmat,
                               double *Gmat );
extern void uhf_form_gmatrix ( const ERIFS * tfiles, const double *PmatA,
                               const double *PmatB, double *GmatA,
                               double *GmatB );
#endif

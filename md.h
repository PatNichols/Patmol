#ifndef _MD_H_
#define _MD_H_
#include "util.h"
#include "pconfig.h"
#include <math.h>
/* mdr.c */
extern void formMDR ( double ***r, double sr, double t, double w,
                      const double *pq, int ltot );
/* mdd.c */
extern void formMDD ( double ***d, double abi, double ax,
                      double bx, int l1, int l2 );
#endif

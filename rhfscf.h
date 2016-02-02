#ifndef _RHFSCF_H_
#define _RHFSCF_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "config.h"
#include "util.h"
#include "basis.h"
#include "shell.h"
#include "center.h"
#include "erifs.h"
#include "oneints.h"
#include "twoints.h"
#include "momints.h"
#include "dpmints.h"
#include "scfaux.h"
#include "spack_math.h"
#include "stopwatch.h"
/* rhfscf.c */
extern int rhf_scf_update ();
extern void rhf_energy ();
extern void rhf_init ();
extern void rhf_finite_field_calculation ();
extern void rhf_dpm_calculation ();
extern void rhf_MullPopAnalysis ( FILE * outfile );
extern void rhf_output ();
extern void rhf_form_xmatrix ();
extern void rhf_form_guess ();
extern void rhf_getDipoleHam ( const double *fld );
extern void rhf_clean();
#endif

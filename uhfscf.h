#ifndef _UHFSCF_H_
#define _UHFSCF_H_
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
/* uhfscf.c */
extern void uhf_clean();
extern int uhf_scf_update ();
extern void uhf_energy ();
extern void uhf_init ();
extern void uhf_finite_field_calculation ();
extern void uhf_dpm_calculation ();
extern void uhf_MullPopAnalysis ( FILE * outfile );
extern void uhf_output ();
extern void uhf_form_xmatrix ();
extern void uhf_form_guess ();
extern void uhf_getDipoleHam ( const double *fld );
#endif

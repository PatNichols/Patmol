#ifndef _RHFSCF_MPI_H_
#define _RHFSCF_MPI_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <mpi.h>
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
/* rhfscf_mpi.h */
extern void rhf_mpi_free();
extern int rhf_mpi_scf_update ();
extern void rhf_mpi_energy ();
extern void rhf_mpi_init ();
extern void rhf_mpi_finite_field_calculation ();
extern void rhf_mpi_dpm_calculation ();
extern void rhf_mpi_MullPopAnalysis ( FILE * outfile );
extern void rhf_mpi_output ();
extern void rhf_mpi_form_xmatrix ();
extern void rhf_mpi_form_guess ();
extern void rhf_mpi_getDipoleHam ( const double *fld );
//extern void rhf_mpi_clean();
#endif

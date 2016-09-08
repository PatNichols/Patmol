#ifndef _UHFSCF_MPI_H_
#define _UHFSCF_MPI_H_
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

/* uhfscf_mpi.h */
extern void uhf_mpi_free();
extern int uhf_mpi_scf_update ();
extern void uhf_mpi_energy ();
extern void uhf_mpi_init ();
extern void uhf_mpi_finite_field_calculation ();
extern void uhf_mpi_dpm_calculation ();
extern void uhf_mpi_MullPopAnalysis ( FILE * outfile );
extern void uhf_mpi_output ();
extern void uhf_mpi_form_xmatrix ();
extern void uhf_mpi_form_guess ();
extern void uhf_mpi_getDipoleHam ( const double *fld );
#endif

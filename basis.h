#ifndef BASIS_H
#define BASIS_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pconfig.h"
#include "util.h"
#include "shell.h"
#include "center.h"
#define DPM_NSHELL 7
#define DPM_MAXL 1
#define DPM_NORB 13

/* basis.c */
extern int number_of_orbitals ();
extern int number_of_shells ();
extern int number_of_centers ();
extern int total_number_of_orbitals ();
extern int total_number_of_shells ();
extern int total_number_of_centers ();
extern int number_of_electrons ();
extern int spin_multiplicity ();
extern int max_scf_iterations ();
extern int skip_center ();
extern int max_lvalue ();
extern int max_lstates ();
extern int do_dpm_calculation ();
extern int do_dft_calculation ();
extern int dft_correlation ();
extern int dft_exchange ();
extern int do_finite_field_calculation ();
extern int scf_acceleration ();
extern int scf_guess ();
extern int print_eigenvectors ();
extern int print_density ();
extern int print_population ();
extern double finite_field_electric_field ();
extern double scf_convergence ();
extern const Shell *shells ();
extern const Center *centers ();
extern int ***lvector_components ();
extern double **normalizationFactors ();
extern void set_center_position ( int icen, const double *Rnew );
extern void set_center_charge ( int icen, double qnew );
extern double nuclear_repulsion ();
extern int row ( const Center * cen );
extern void assign_dpm_shells ( Shell * sh, int poscen, int off );
extern void init_Basis ();
extern void dpm_augment ();
extern void basis_output ( double *, double *, double * );
#endif

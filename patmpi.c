#include <stdio.h>
#include <stdlib.h>
#define _USING_MPI
#include <mpi.h>
#include "pconfig.h"
#include "util.h"
#include "basis.h"
#include "rhfscf_mpi.h"
#include "uhfscf_mpi.h"

int
main ( int argc, char **argv )
{
    int ne, sm, id, np;
    double st, et;
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &id);
    MPI_Comm_size( MPI_COMM_WORLD, &np);
    if (id==0) 
        fprintf(stderr,"there are %d processes\n",np);
    
    st = MPI_Wtime();
    init_Basis ();
    ne = number_of_electrons ();
    sm = spin_multiplicity ();
    if ( ! ( ne % 2 ) && ( sm == 1 ) ) {
        rhf_mpi_init ();
        rhf_mpi_energy ();
        if ( do_finite_field_calculation () )
            rhf_mpi_finite_field_calculation ();
        if ( do_dpm_calculation () )
            rhf_mpi_dpm_calculation ();
        rhf_mpi_free();
    } else {
        uhf_mpi_init ();
        uhf_mpi_energy ();
        if ( do_finite_field_calculation () )
            uhf_mpi_finite_field_calculation ();
        if ( do_dpm_calculation () )
            uhf_mpi_dpm_calculation ();
        uhf_mpi_free();
    }
    et = MPI_Wtime() - st;
    fprintf ( stderr, "WALL TIME = %lg seconds\n", et );
    MPI_Finalize ();
    return EXIT_SUCCESS;
}

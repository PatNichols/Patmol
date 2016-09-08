#include <stdio.h>
#include <stdlib.h>
#include "pconfig.h"
#include "util.h"
#include "basis.h"
#include "rhfscf.h"
#include "uhfscf.h"

int
main ( int argc, char **argv )
{
    int ne, sm;
    init_Basis ();
    ne = number_of_electrons ();
    sm = spin_multiplicity ();
//  fprintf (stderr, "Ne = %d %d\n", ne, sm);
//  fprintf (stderr, "%d %d\n", (ne % 2), (sm == 1));
    if ( ! ( ne % 2 ) && ( sm == 1 ) ) {
        rhf_init ();
        rhf_energy ();
        if ( do_finite_field_calculation () )
            rhf_finite_field_calculation ();
        if ( do_dpm_calculation () )
            rhf_dpm_calculation ();
        rhf_clean();
    } else {
        uhf_init ();
        uhf_energy ();
        if ( do_finite_field_calculation () )
            uhf_finite_field_calculation ();
        if ( do_dpm_calculation () )
            uhf_dpm_calculation ();
        uhf_clean();
    }
}

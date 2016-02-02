#include "rhfscf.h"

stopwatch_t gtimer;
stopwatch_t dtimer;
static double scf_eps, ediff, pdiff, eold, enrg, nucrep, enrg_gs, enrg0;
static double *fock, *hmat, *tmat, *smat, *wmat, *xmat, *work;
static double *gmat, *pmat, *cmat;
static double *pmat1, *pmat2, *evals, *pmat_gs, *hmat_gs;
static int no2, tno2, its, nocc, norb, tnorb, scf_accel, iteration;
static int max_scf_its;
static ERIFS *mints;
static ERIFS *xints = 0;
static const double EPS10 = 10 * DBL_EPSILON;

void rhf_clean()
{
    clean_ERIFS ( mints );
    if ( xints ) clean_ERIFS ( xints );
}

int
rhf_scf_update ()
{
    int i;
    ++iteration;
    for ( i = 0; i < no2; ++i ) {
        gmat[i] = 0.;
    }
    stopwatch_start(&gtimer);
    rhf_form_gmatrix ( mints, pmat, gmat );
    if ( xints ) {
        rhf_form_gmatrix ( xints, pmat, gmat );
    }
    stopwatch_stop(&gtimer);
    for ( i = 0; i < no2; ++i ) {
        fock[i] = hmat[i] + gmat[i];
    }
    stopwatch_start(&dtimer);
    spack_trans ( norb, fock, xmat, work );
    spack_diag ( norb, fock, wmat, evals, work );
    stopwatch_stop(&dtimer);
    form_cmatrix ( norb, cmat, xmat, wmat );
    for ( i = 0; i < no2; ++i ) {
        pmat2[i] = pmat1[i];
        pmat1[i] = pmat[i];
    }
    form_pmatrix ( norb, nocc, cmat, pmat );
    pdiff = pdiffNorm ( norb, pmat, pmat1 );
    enrg = 2. * calc_energy ( norb, pmat1, hmat, gmat );
    ediff = enrg - eold;
    eold = enrg;
    if ( iteration == 1 )
        return 0;
    if ( fabs ( ediff ) < EPS10 ) {
        return 1;
    }
    if ( pdiff < scf_eps ) {
        return 1;
    }
    if ( ediff > 0. ) {
        if ( scf_accel == 1 ) {
            for ( i = 0; i < no2; ++i ) {
                pmat[i] = 0.5 * ( pmat1[i] + pmat2[i] );
            }
        }
        if ( scf_accel == 2 ) {
            scfAccelerator ( iteration, no2, pmat, pmat1, pmat2 );
        }
    }
    return 0;
}

void
rhf_energy ()
{
    stopwatch_t timer;
    int done;
    
    stopwatch_clear(&timer);
    stopwatch_clear(&gtimer);
    stopwatch_clear(&dtimer);
    stopwatch_start(&timer);
    eold = 0.0;
    nucrep = nuclear_repulsion ();
    rhf_scf_update ();
    enrg0 = enrg + nucrep;
    done = 0;
    iteration = 0;
    fprintf ( stderr, "ITERATION = %12d EDIFF= %lg\n", iteration, ediff );
    fprintf ( stderr, "  energy = %lg\n", ( enrg + nucrep ) );
    while ( !done && iteration < max_scf_its ) {
        done = rhf_scf_update ();
        fprintf ( stderr, "ITERATION = %12d EDIFF= %lg\n", iteration, ediff );
        fprintf ( stderr, "  energy = %lg\n", ( enrg + nucrep ) );
    }
    if ( iteration == max_scf_its ) {
        fprintf ( stderr, "WARNING CONVERGENCE NOT REACHED IN SCF\n" );
    } else {
        fprintf ( stderr, "Convergence reached!\n" );
    }
    stopwatch_stop(&timer);
    fprintf(stderr,"scf time = %lg\n",stopwatch_elapsed_time(&timer));
    fprintf(stderr,"gmat time= %lg\n",stopwatch_elapsed_time(&gtimer));
    fprintf(stderr,"diag time= %lg\n",stopwatch_elapsed_time(&dtimer));
    memset ( pmat_gs, 0, sizeof ( double ) * no2 );
    memcpy ( pmat_gs, pmat, sizeof ( double ) * no2 );
    memcpy ( hmat_gs, hmat, sizeof ( double ) * no2 );
    enrg_gs = enrg + nucrep;
    fprintf ( stderr, "Doing output\n" );
    rhf_output ();
    stopwatch_clear(&timer);
    stopwatch_start(&timer);
    fprintf ( stderr, "Analyzing moments\n" );
    rhf_analyzeMoments ( norb, pmat, cmat );
    stopwatch_stop(&timer);
    fprintf(stderr,"moments time = %lg\n",stopwatch_elapsed_time(&timer));    
    fprintf ( stderr, "Done Here!\n" );
}

void
rhf_init ()
{
    int tnsq, ne;
    ne = number_of_electrons ();
    if ( ( ne % 2 ) || ( spin_multiplicity () != 1 ) ) {
        fatalError ( "Odd number of electrons in RHF" );
    }
    nocc = ne / 2;
    norb = number_of_orbitals ();
    no2 = ( norb * ( norb + 1 ) ) / 2;
    tnorb = total_number_of_orbitals ();
    tno2 = ( tnorb * ( tnorb + 1 ) ) / 2;
    tnsq = tnorb * tnorb;
    fock = dalloc ( tno2 );
    gmat = dalloc ( tno2 );
    hmat = dalloc ( tno2 );
    pmat = dalloc ( tno2 );
    pmat1 = dalloc ( tno2 );
    pmat2 = dalloc ( tno2 );
    smat = dalloc ( tno2 );
    tmat = dalloc ( tno2 );
    pmat_gs = dalloc ( tno2 );
    hmat_gs = dalloc ( tno2 );
    cmat = dalloc ( tnsq );
    xmat = dalloc ( tnsq );
    wmat = dalloc ( tnsq );
    work = dalloc ( tnsq );
    evals = dalloc ( tnorb );
    mints = init_ERIFS ( 0, norb, "MINTS.DAT" );
    xints = 0;
    one_electron_ints ( smat, tmat, hmat );
    basis_output ( smat, tmat, hmat );
    fprintf ( stderr, "One electron ints done!\n" );
    two_electron_ints ( mints, 0 );
    fprintf ( stderr, "Two Electron ints done!\n" );
    momints ();
    fprintf ( stderr, "Moments done\n" );
    max_scf_its = max_scf_iterations ();
    scf_accel = scf_acceleration ();
    scf_eps = scf_convergence ();
    if ( scf_eps < ( EPS10 * norb * norb ) ) {
        scf_eps = EPS10 * norb * norb;
        fprintf ( stderr, "EPS was set to low. Now set to %lg\n", scf_eps );
    }
    fprintf ( stderr, "norb = %d nshell = %d ncen = %d  skipcen= %d\n",
              number_of_orbitals (),
              number_of_shells (), number_of_centers (), skip_center () );
    fprintf ( stderr, "tnorb = %d tnshell = %d tncen = %d  skipcen= %d\n",
              total_number_of_orbitals (),
              total_number_of_shells (),
              total_number_of_centers (), skip_center () );
    rhf_form_xmatrix ();
    rhf_form_guess ();
}

void
rhf_finite_field_calculation ()
{
    int idir, done;
    double efx, alfa0, alfa2;
    double field[3];
    double alfa[3];
    FILE *fp;
    efx = finite_field_electric_field ();
    for ( idir = 0; idir < 3; ++idir ) {
        field[0] = field[1] = field[2] = 0.0;
        field[idir] = efx;
        memcpy ( hmat, hmat_gs, sizeof ( double ) * no2 );
        memcpy ( pmat, pmat_gs, sizeof ( double ) * no2 );
        rhf_getDipoleHam ( field );
        rhf_scf_update ();
        enrg0 = enrg + nucrep;
        done = 0;
        iteration = 0;
        while ( !done && iteration < max_scf_its ) {
            done = rhf_scf_update ();
        }
        if ( iteration == max_scf_its ) {
            fprintf ( stderr, "WARNING CONVERGENCE NOT REACHED IN SCF\n" );
        }
        ediff = enrg + nucrep - enrg0;
        alfa[idir] = -2. * ediff / efx / efx;
    }
    fp = openFile ( "FiniteField.out", "w" );
    fprintf ( fp, "  FINITE FIELD ANALYSIS\n" );
    fprintf ( fp, "  FIELD = %15.6le\n", efx );
    fprintf ( fp, "  Alfa xx = %15.6le\n", alfa[0] );
    fprintf ( fp, "  Alfa yy = %15.6le\n", alfa[1] );
    fprintf ( fp, "  Alfa zz = %15.6le\n", alfa[2] );
    alfa0 = ( alfa[0] + alfa[1] + alfa[2] ) / 3.;
    alfa2 = ( 2.0 * alfa[2] - alfa[1] - alfa[2] ) / 3.;
    fprintf ( fp, "  Alfa 0  = %15.6le\n", alfa0 );
    fprintf ( fp, "  Alfa 2  = %15.6le\n", alfa2 );
    fclose ( fp );
}

void
rhf_dpm_calculation ()
{
    int ipt, npts, done, ono, onsh, poscen;
    double rp[3], vpol, vstat, energ0, r2, alfa;
    FILE *infile, *outfile, *potfile;
    poscen = number_of_centers ();
    ono = number_of_orbitals ();
    onsh = number_of_shells ();
    xints = init_ERIFS ( norb, tnorb, "XINTS.DAT" );
    dpm_augment ();
    norb = tnorb;
    no2 = ( tnorb * ( tnorb + 1 ) ) / 2;
    infile = openFile ( "posgrid.dat", "r" );
    fscanf ( infile, "%d", &npts );
    outfile = openFile ( "dpmvpol.dat", "w" );
    for ( ipt = 0; ipt < npts; ++ipt ) {
        fscanf ( infile, "%lg %lg %lg", rp, rp + 1, rp + 2 );
        set_center_position ( poscen, rp );
        one_electron_ints ( smat, tmat, hmat );
        dpm_ints ( hmat );
        two_electron_ints ( xints, onsh );
        nucrep = nuclear_repulsion ();
        memcpy ( pmat, pmat_gs, sizeof ( double ) * tno2 );
        rhf_form_xmatrix ();
        rhf_scf_update ();
        enrg0 = enrg + nucrep;
        done = 0;
        iteration = 0;
        while ( !done && iteration < max_scf_its ) {
            done = rhf_scf_update ();
        }
        if ( iteration == max_scf_its ) {
            fprintf ( stderr, "WARNING CONVERGENCE NOT REACHED IN SCF\n" );
        }
        ediff = enrg + nucrep - enrg0;
        vpol = ediff;
        vstat = energ0 - enrg_gs;
        r2 = normSqr ( rp );
        alfa = -2. * vpol * r2 * r2;
        fprintf ( outfile, "%15.6le %15.6le %15.6le %20.10le %15.6le %20.10le\n",
                  rp[0], rp[1], rp[2], vpol, alfa, vstat );
        fflush ( outfile );
    }
    fclose ( outfile );
    fclose ( infile );
}

void
rhf_MullPopAnalysis ( FILE * outfile )
{
    int i, j, k, ii, ir, ik, jk, lv, cn, nls;
    double sum, qc;
    int *icen;
    double *Overlap, *NetChg;
    FILE *poutfile;
    const int ncen = number_of_centers ();
    const int nshell = number_of_shells ();
    const Shell *shell = shells ();
    const Center *center = centers ();
    Overlap = dalloc ( norb * norb );
    NetChg = dalloc ( ncen );
    icen = ialloc ( norb );
    for ( i = 0; i < nshell; ++i ) {
        lv = ( shell + i )->lsh;
        cn = ( shell + i )->cen;
        ii = ( shell + i )->off;
        nls = ( shell + i )->nst;
        for ( j = 0; j < nls; ++j ) {
            ir = ii + j;
            icen[ir] = cn;
        }
    }
    for ( i = 0; i < norb; ++i ) {
        ii = i * ( i + 1 ) / 2;
        for ( j = 0; j <= i; ++j ) {
            ik = ii;
            jk = j * ( j + 1 ) / 2;
            sum = 0.0;
            for ( k = 0; k <= j; ++ik, ++jk, ++k ) {
                sum += pmat[ik] * smat[jk];
            }
            --jk;
            for ( k = j + 1; k <= i; ++ik, ++k ) {
                jk = jk + k;
                sum += pmat[ik] * smat[jk];
            }
            --ik;
            for ( k = i + 1; k < norb; ++k ) {
                jk = jk + k;
                ik = ik + k;
                sum += pmat[ik] * smat[jk];
            }
            Overlap[j * norb + i] = Overlap[i * norb + j] = ( sum + sum );
        }
    }
    for ( i = 0; i < ncen; ++i ) {
        NetChg[i] = ( center + i )->q;
    }
    for ( i = 0; i < norb; ++i ) {
        cn = icen[i];
        NetChg[cn] -= Overlap[i * norb + i];
    }
    fprintf ( outfile, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    fprintf ( outfile, "           Mulliken Populations \n" );
    fprintf ( outfile, "Orbital  Net Population\n" );
    for ( i = 0; i < norb; ++i ) {
        fprintf ( outfile, "%7u %25.16le\n", i + 1, Overlap[i * norb + i] );
    }
    fprintf ( outfile, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    fprintf ( outfile, "Mulliken Atomic Charges \n" );
    fprintf ( outfile, "Orbital   Nuclear Charge   Net Charge \n" );
    for ( i = 0; i < ncen; ++i ) {
        qc = ( center + i )->q;
        fprintf ( outfile, "%7u %15.10lf %15.10lf \n", i + 1, qc, NetChg[i] );
    }
    if ( print_population () ) {
        poutfile = openFile ( "overlaps.dat", "w" );
        fprintf ( poutfile, "   MULLIKEN OVERLAP POPULATIONS\n" );
        fprintf ( poutfile, "  i      j                   value \n" );
        for ( i = 0; i < norb; ++i ) {
            for ( j = 0; j <= i; ++j ) {
                fprintf ( poutfile, "%5d %5d  %24.16le\n", i, j,
                          Overlap[i * norb + j] );
            }
        }
        fclose ( poutfile );
    }
    free ( icen );
    free ( NetChg );
    free ( Overlap );
}

void
rhf_output ()
{
    int i, j, nex, occup, end;
    double trace_t, virial;
    const double *Cm;
    FILE *out;
    trace_t = 2.0 * spack_trace_prod ( norb, pmat, tmat );
    virial = fabs ( ( enrg_gs - trace_t ) * 0.5 / ( trace_t ) );
    out = openFile ( "short.gs.dat", "w" );
    fprintf ( out, "%25.16le\n", enrg0 );
    fprintf ( out, "%25.16le\n", enrg_gs );
    fprintf ( out, "%25.16le\n", ediff );
    fclose ( out );
    out = openFile ( "scfout.gs.dat", "w" );
    if ( max_scf_its <= iteration ) {
        fprintf ( out, "WARNING CONVERGENCE _NOT_ REACHED in %d iterations!\n",
                  its );
    }
    fprintf ( out, "Final Iteration              = %12d\n", iteration );
    fprintf ( out, "Hartree Fock Energy          = %25.15le Hartree\n",
              ( enrg + nucrep ) );
    fprintf ( out, "Electronic Energy            = %25.15le Hartree\n", enrg );
    fprintf ( out, "Nuclear Rep. Energy          = %25.15le Hartree\n", nucrep );
    fprintf ( out, "Kinetic Energy               = %25.15le Hartree\n", trace_t );
    fprintf ( out, "Difference in Final Energies = %25.15le Hartree\n", ediff );
    fprintf ( out, "Pmatrix diff norm            = %25.15le\n", pdiff );
    fprintf ( out, "virial                       = %25.15le\n", virial );
    fprintf ( out, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    fprintf ( out, "               Orbital Energies\n" );
    fprintf ( out, "Orbital Energy          Occupancy\n" );
    for ( i = 0; i < nocc; i++ ) {
        fprintf ( out, "%7u %25.16le %12u\n", i + 1, evals[i], 2 );
    }
    for ( i = nocc; i < norb; i++ ) {
        fprintf ( out, "%7u %25.16le %12u\n", i + 1, evals[i], 0 );
    }
    fprintf ( out, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    rhf_MullPopAnalysis ( out );
    if ( print_eigenvectors () ) {
        fprintf ( out, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
        fprintf ( out, "              Fock Eigenvectors\n" );
        fprintf ( out, "\n" );
        for ( i = 0; i < norb; ++i ) {
            occup = ( i < nocc ) ? 2 : 0;
            fprintf ( out, "Column %5d Occupancy= %2d EigenEnergy= %25.15le\n",
                      i, occup, evals[i] );
            nex = norb - 4 * norb / 4;
            end = 4 * ( norb / 4 );
            Cm = cmat + i;
            if ( norb > 4 ) {
                for ( j = 0; j < norb; j += 4 ) {
                    fprintf ( out, "%18.10le %18.10le %18.10le %18.10le\n",
                              ( *Cm ), ( * ( Cm + norb ) ), ( * ( Cm + 2 * norb ) ),
                              ( * ( Cm + 3 * norb ) ) );
                    Cm += 4 * norb;
                }
            }
            for ( j = end; j < norb; ++j ) {
                fprintf ( out, "%18.10le ", *Cm );
                Cm += norb;
            }
            fprintf ( out, "\n" );
        }
    }
    fclose ( out );
    out = openFile ( "ORBS.DAT", "w" );
    fwrite ( evals, sizeof ( double ), norb, out );
    fclose ( out );
    if ( print_density () ) {
        out = openFile ( "PMATRIX.OUT", "w" );
        fwrite ( pmat, sizeof ( double ), no2, out );
        fclose ( out );
    }
    if ( print_eigenvectors () ) {
        out = openFile ( "CMATRIX.OUT", "w" );
        fwrite ( cmat, sizeof ( double ), norb * norb, out );
        fclose ( out );
    }
}


void
rhf_form_xmatrix ()
{
    int i, j;
    double tmp;
    memcpy ( fock, smat, sizeof ( double ) * no2 );
    spack_diag ( norb, fock, wmat, evals, work );
    for ( i = 0; i < norb; ++i ) {
        tmp = evals[i];
        if ( tmp <= DBL_EPSILON ) {
            fprintf ( stderr, "Change your basis set!\n" );
            fatalError ( "Negative or Zero Eigenvalues in Overlap Matrix!" );
        }
        evals[i] = 1. / sqrt ( tmp );
    }
    for ( i = 0; i < norb; ++i ) {
        for ( j = 0; j < norb; ++j ) {
            xmat[i * norb + j] = wmat[i * norb + j] * evals[j];
        }
    }
}

void
rhf_form_guess ()
{
    FILE *fp;
    if ( scf_guess () == 1 ) {
        fprintf ( stderr, "Reading Pmatrix guess\n" );
        fp = openFile ( "PMATRIX.GUESS", "r" );
        fread ( pmat, sizeof ( double ), no2, fp );
        fclose ( fp );
        return;
    }
    fprintf ( stderr, "Generating Pmatrix Guess!\n" );
    memcpy ( fock, hmat, sizeof ( double ) * no2 );
    spack_trans ( norb, fock, xmat, work );
    spack_diag ( norb, fock, wmat, evals, work );
    form_cmatrix ( norb, cmat, xmat, wmat );
    form_pmatrix ( norb, nocc, cmat, pmat );
}

void
rhf_getDipoleHam ( const double *fld )
{
    int i, ir, jr, ij;
    double dm[3], qm[6];
    FILE *infile = openFile ( "MOMINTS.DAT", "r" );
    for ( i = 0; i < no2; ++i ) {
        fread ( &ir, sizeof ( int ), 1, infile );
        fread ( &jr, sizeof ( int ), 1, infile );
        fread ( dm, sizeof ( double ), 3, infile );
        fread ( qm, sizeof ( double ), 6, infile );
        ij = ( ir * ( ir + 1 ) ) / 2 + jr;
        hmat[ij] -= ( fld[0] * dm[0] + fld[1] * dm[1] + fld[2] * dm[2] );
    }
    fclose ( infile );
}

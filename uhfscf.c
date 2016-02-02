
#include "uhfscf.h"

static double scf_eps, ediff, pdiff, eold, enrg, nucrep, enrg_gs, enrg0;
static double *fock, *hmat, *tmat, *smat, *wmat, *xmat, *work;
static double *gmatA, *pmatA, *cmatA, *pmat1A, *pmat2A;
static double *gmatB, *pmatB, *cmatB;
static double *pmat1B, *pmat2B, *evalsA, *evalsB, *pmatA_gs,
       *pmatB_gs, *hmat_gs;
static int no2, tno2, its, noccA, noccB, norb, tnorb, scf_accel, iteration;
static int max_scf_its;
static ERIFS *mints;
static ERIFS *xints = 0;
static const double EPS10 = 10 * DBL_EPSILON;

void uhf_clean()
{
    clean_ERIFS ( mints );
    if ( xints ) clean_ERIFS ( xints );
}

int
uhf_scf_update ()
{
    int i;
    ++iteration;
    __builtin_memset(gmatA,0,sizeof(double)*no2);
    __builtin_memset(gmatB,0,sizeof(double)*no2);
    uhf_form_gmatrix ( mints, pmatA, pmatB, gmatA, gmatB );
    if ( xints ) {
        uhf_form_gmatrix ( xints, pmatA, pmatB, gmatA, gmatB );
    }
    for ( i = 0; i < no2; ++i ) {
        fock[i] = hmat[i] + gmatA[i];
    }
    spack_trans ( norb, fock, xmat, work );
    spack_diag ( norb, fock, wmat, evalsA, work );
    form_cmatrix ( norb, cmatA, xmat, wmat );
    __builtin_memcpy(pmat2A,pmat1A,sizeof(double)*no2);
    __builtin_memcpy(pmat1A,pmatA,sizeof(double)*no2);
    form_pmatrix ( norb, noccA, cmatA, pmatA );
    pdiff = pdiffNorm ( norb, pmatA, pmat1A );
    for ( i = 0; i < no2; ++i ) {
        fock[i] = hmat[i] + gmatB[i];
    }
    spack_trans ( norb, fock, xmat, work );
    spack_diag ( norb, fock, wmat, evalsB, work );
    form_cmatrix ( norb, cmatB, xmat, wmat );
    __builtin_memcpy(pmat2B,pmat1B,sizeof(double)*no2);
    __builtin_memcpy(pmat1B,pmatB,sizeof(double)*no2);
    form_pmatrix ( norb, noccB, cmatB, pmatB );
    pdiff += pdiffNorm ( norb, pmatB, pmat1B );
    enrg = calc_energy ( norb, pmat1A, hmat, gmatA ) +
           calc_energy ( norb, pmat1B, hmat, gmatB );
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
                pmatA[i] = 0.5 * ( pmat1A[i] + pmat2A[i] );
                pmatB[i] = 0.5 * ( pmat1B[i] + pmat2B[i] );
            }
        }
        if ( scf_accel == 2 ) {
            scfAccelerator ( iteration, no2, pmatA, pmat1A, pmat2A );
            scfAccelerator ( iteration, no2, pmatB, pmat1B, pmat2B );
        }
    }
    return 0;
}

void
uhf_energy ()
{
    int done;
    eold = 0.0;
    nucrep = nuclear_repulsion ();
    uhf_scf_update ();
    enrg0 = enrg + nucrep;
    done = 0;
    iteration = 0;
    fprintf ( stderr, "ITERATION = %12d EDIFF= %lg\n", iteration, ediff );
    fprintf ( stderr, "  energy = %lg\n", ( enrg + nucrep ) );
    while ( !done && iteration < max_scf_its ) {
        done = uhf_scf_update ();
        fprintf ( stderr, "ITERATION = %12d EDIFF= %lg\n", iteration, ediff );
        fprintf ( stderr, "  energy = %lg\n", ( enrg + nucrep ) );
    }
    if ( iteration == max_scf_its ) {
        fprintf ( stderr, "WARNING CONVERGENCE NOT REACHED IN SCF\n" );
    } else {
        fprintf ( stderr, "Convergence reached!\n" );
    }
    memset ( pmatB_gs, 0, sizeof ( double ) * no2 );
    memcpy ( pmatB_gs, pmatB, sizeof ( double ) * no2 );
    memset ( pmatA_gs, 0, sizeof ( double ) * no2 );
    memcpy ( pmatA_gs, pmatA, sizeof ( double ) * no2 );
    memcpy ( hmat_gs, hmat, sizeof ( double ) * no2 );
    enrg_gs = enrg + nucrep;
    fprintf ( stderr, "Doing output\n" );
    uhf_output ();
    fprintf ( stderr, "Analyzing moments\n" );
    uhf_analyzeMoments ( norb, pmatA, pmatB, cmatA, cmatB );
    fprintf ( stderr, "Done Here!\n" );
}

void
uhf_init ()
{
    int tnsq, ne, sm;
    ne = number_of_electrons ();
    sm = spin_multiplicity ();
    noccB = ( ne - ( sm - 1 ) ) / 2;
    noccA = ne - noccB;
    norb = number_of_orbitals ();
    no2 = ( norb * ( norb + 1 ) ) / 2;
    tnorb = total_number_of_orbitals ();
    tno2 = ( tnorb * ( tnorb + 1 ) ) / 2;
    tnsq = tnorb * tnorb;
    fock = dalloc ( tno2 );
    hmat = dalloc ( tno2 );
    gmatA = dalloc ( tno2 );
    pmatA = dalloc ( tno2 );
    pmat1A = dalloc ( tno2 );
    pmat2A = dalloc ( tno2 );
    pmatA_gs = dalloc ( tno2 );
    cmatA = dalloc ( tnsq );
    gmatB = dalloc ( tno2 );
    pmatB = dalloc ( tno2 );
    pmat1B = dalloc ( tno2 );
    pmat2B = dalloc ( tno2 );
    pmatB_gs = dalloc ( tno2 );
    cmatB = dalloc ( tnsq );
    smat = dalloc ( tno2 );
    tmat = dalloc ( tno2 );
    hmat_gs = dalloc ( tno2 );
    xmat = dalloc ( tnsq );
    wmat = dalloc ( tnsq );
    work = dalloc ( tnsq );
    evalsA = dalloc ( tnorb );
    evalsB = dalloc ( tnorb );
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
    uhf_form_xmatrix ();
    uhf_form_guess ();
}

void
uhf_finite_field_calculation ()
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
        memcpy ( pmatA, pmatA_gs, sizeof ( double ) * no2 );
        memcpy ( pmatB, pmatB_gs, sizeof ( double ) * no2 );
        uhf_getDipoleHam ( field );
        uhf_scf_update ();
        enrg0 = enrg + nucrep;
        done = 0;
        iteration = 0;
        while ( !done && iteration < max_scf_its ) {
            done = uhf_scf_update ();
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
uhf_dpm_calculation ()
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
    no2 = tno2;
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
        memcpy ( pmatA, pmatA_gs, sizeof ( double ) * tno2 );
        memcpy ( pmatB, pmatB_gs, sizeof ( double ) * tno2 );
        uhf_form_xmatrix ();
        uhf_scf_update ();
        enrg0 = enrg + nucrep;
        done = 0;
        iteration = 0;
        while ( !done && iteration < max_scf_its ) {
            done = uhf_scf_update ();
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
uhf_MullPopAnalysis ( FILE * outfile )
{
    int i, j, k, ii, ir, ik, jk, lv, cn, nls;
    double sum, qc;
    int *icen;
    double *OverlapA, *OverlapB, *NetChg;
    FILE *poutfile;
    const int ncen = number_of_centers ();
    const int nshell = number_of_shells ();
    const Shell *shell = shells ();
    const Center *center = centers ();
    OverlapA = dalloc ( norb * norb );
    OverlapB = dalloc ( norb * norb );
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
                sum += pmatA[ik] * smat[jk];
            }
            --jk;
            for ( k = j + 1; k <= i; ++ik, ++k ) {
                jk = jk + k;
                sum += pmatA[ik] * smat[jk];
            }
            --ik;
            for ( k = i + 1; k < norb; ++k ) {
                jk = jk + k;
                ik = ik + k;
                sum += pmatA[ik] * smat[jk];
            }
            OverlapA[j * norb + i] = OverlapA[i * norb + j] = sum;
        }
    }
    for ( i = 0; i < norb; ++i ) {
        ii = i * ( i + 1 ) / 2;
        for ( j = 0; j <= i; ++j ) {
            ik = ii;
            jk = j * ( j + 1 ) / 2;
            sum = 0.0;
            for ( k = 0; k <= j; ++ik, ++jk, ++k ) {
                sum += pmatB[ik] * smat[jk];
            }
            --jk;
            for ( k = j + 1; k <= i; ++ik, ++k ) {
                jk = jk + k;
                sum += pmatB[ik] * smat[jk];
            }
            --ik;
            for ( k = i + 1; k < norb; ++k ) {
                jk = jk + k;
                ik = ik + k;
                sum += pmatB[ik] * smat[jk];
            }
            OverlapB[j * norb + i] = OverlapB[i * norb + j] = sum;
        }
    }
    for ( i = 0; i < ncen; ++i ) {
        NetChg[i] = ( center + i )->q;
    }
    for ( i = 0; i < norb; ++i ) {
        cn = icen[i];
        NetChg[cn] -= OverlapA[i * norb + i];
    }
    fprintf ( outfile, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    fprintf ( outfile, "           Mulliken Populations \n" );
    fprintf ( outfile, "Orbital  Net Population\n" );
    for ( i = 0; i < norb; ++i ) {
        fprintf ( outfile, "%7u %25.16le %25.16le\n", i + 1,
                  OverlapA[i * norb + i], OverlapB[i * norb + i] );
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
                fprintf ( poutfile, "%5d %5d  %24.16le %24.16le\n", i, j,
                          OverlapA[i * norb + j], OverlapB[i * norb + j] );
            }
        }
        fclose ( poutfile );
    }
    free ( icen );
    free ( NetChg );
    free ( OverlapA );
    free ( OverlapB );
}

void
uhf_output ()
{
    int i, j, nex, occup, end, fa, fb;
    double trace_t, virial;
    const double *Cm;
    FILE *out;
    trace_t = spack_trace_prod ( norb, pmatA, tmat ) +
              spack_trace_prod ( norb, pmatB, tmat );
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
    for ( i = 0; i < norb; i++ ) {
        fa = fb = 0;
        if ( i < noccA )
            fa = 1;
        if ( i < noccB )
            fb = 1;
        fprintf ( out, "%7u %25.16le %6d %25.16le %6d \n", (i + 1), evalsA[i], fa,
                  evalsB[i], fb );
    }
    fprintf ( out, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
    uhf_MullPopAnalysis ( out );
    if ( print_eigenvectors () ) {
        fprintf ( out, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
        fprintf ( out, "       Alpha  Fock Eigenvectors\n" );
        fprintf ( out, "\n" );
        for ( i = 0; i < norb; ++i ) {
            occup = ( i < noccA ) ? 2 : 0;
            fprintf ( out, "Column %5d Occupancy= %2d EigenEnergy= %25.15le\n",
                      i, occup, evalsA[i] );
            nex = norb - 4 * norb / 4;
            end = 4 * ( norb / 4 );
            Cm = cmatA + i;
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
        fprintf ( out, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
        fprintf ( out, "       Beta  Fock Eigenvectors\n" );
        fprintf ( out, "\n" );
        for ( i = 0; i < norb; ++i ) {
            occup = ( i < noccB ) ? 2 : 0;
            fprintf ( out, "Column %5d Occupancy= %2d EigenEnergy= %25.15le\n",
                      i, occup, evalsB[i] );
            nex = norb - 4 * norb / 4;
            end = 4 * ( norb / 4 );
            Cm = cmatB + i;
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
    fwrite ( evalsA, sizeof ( double ), norb, out );
    fwrite ( evalsB, sizeof ( double ), norb, out );
    fclose ( out );
    if ( print_density () ) {
        out = openFile ( "PMATRIX.OUT", "w" );
        fwrite ( pmatA, sizeof ( double ), no2, out );
        fwrite ( pmatB, sizeof ( double ), no2, out );
        fclose ( out );
    }
    if ( print_eigenvectors () ) {
        out = openFile ( "CMATRIX.OUT", "w" );
        fwrite ( cmatA, sizeof ( double ), norb * norb, out );
        fwrite ( cmatB, sizeof ( double ), norb * norb, out );
        fclose ( out );
    }
}


void
uhf_form_xmatrix ()
{
    int i, j;
    double tmp;
    __builtin_memcpy ( fock, smat, sizeof ( double ) * no2 );
    spack_diag ( norb, fock, wmat, evalsA, work );
    for ( i = 0; i < norb; ++i ) {
        tmp = evalsA[i];
        if ( tmp <= DBL_EPSILON ) {
            fprintf ( stderr, "Change your basis set!\n" );
            fatalError ( "Negative or Zero Eigenvalues in Overlap Matrix!" );
        }
        evalsA[i] = 1. / sqrt ( tmp );
    }
    for ( i = 0; i < norb; ++i ) {
        for ( j = 0; j < norb; ++j ) {
            xmat[i * norb + j] = wmat[i * norb + j] * evalsA[j];
        }
    }
}

void
uhf_form_guess ()
{
    FILE *fp;
    if ( scf_guess () == 1 ) {
        fprintf ( stderr, "Reading Pmatrix guess\n" );
        fp = openFile ( "PMATRIX.GUESS", "r" );
        fread ( pmatA, sizeof ( double ), no2, fp );
        fread ( pmatB, sizeof ( double ), no2, fp );
        fclose ( fp );
        return;
    }
    fprintf ( stderr, "Generating Pmatrix Guess!\n" );
    memcpy ( fock, hmat, sizeof ( double ) * no2 );
    spack_trans ( norb, fock, xmat, work );
    spack_diag ( norb, fock, wmat, evalsA, work );
    form_cmatrix ( norb, cmatA, xmat, wmat );
    __builtin_memcpy ( cmatB, cmatA, sizeof ( double ) *norb * norb );
    form_pmatrix ( norb, noccA, cmatA, pmatA );
    form_pmatrix ( norb, noccB, cmatB, pmatB );
}

void
uhf_getDipoleHam ( const double *fld )
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

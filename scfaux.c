#include "scfaux.h"
#define BINSIZE 4000

void
form_cmatrix ( int n, double *cm, const double *xm, double *wm )
{
    register int i, j, k;
    register double sum;
    double *wr, *wc;
    const double *xi, *wj;
    for ( i = 0; i < n; ++i ) {
        xi = xm + i * n;
        for ( j = 0; j < n; ++j ) {
            sum = 0.0;
            for ( k = 0; k < n; ++k ) {
                sum += xi[k] * wm[k * n + j];
            }
            cm[i * n + j] = sum;
        }
    }
}

void
form_pmatrix ( int no, int noc, const double *cm, double *pm )
{
    register int k, ij, j, i;
    register double sum;
    const double *ci, *cj;
    ij = 0;
    for ( i = 0; i < no; ++i ) {
        ci = cm + i * no;
        for ( j = 0; j <= i; ++j, ++ij ) {
            cj = cm + j * no;
            sum = 0.0;
            for ( k = 0; k < noc; ++k ) {
                sum += ci[k] * cj[k];
            }
            pm[ij] = sum;
        }
    }
}

void
form_pmatrix2 ( int no, int *noc, const double *cm, double *pm )
{
    register int k, ij, j, i;
    register double sum;
    const double *ci, *cj;
    ij = 0;
    for ( i = 0; i < no; ++i ) {
        ci = cm + i * no;
        for ( j = 0; j <= i; ++j, ++ij ) {
            cj = cm + j * no;
            sum = 0.0;
            for ( k = 0; k < no; ++k )
                sum += noc[k] * ci[k] * cj[k];
            pm[ij] = sum;
        }
    }
}


void
scfAccelerator ( int iter, int no2, double *p0, double *p1, double *p2 )
{
    register int i;
    register double p00, p11, p01, beta, s0, s1, betam1;
    if ( iter > 2 ) {
        p00 = p01 = p11 = 0.0;
        for ( i = 0; i < no2; ++i ) {
            s0 = p0[i] - p1[i];
            s1 = p1[i] - p2[i];
            p00 += s0 * s0;
            p01 += s0 * s1;
            p11 += s1 * s1;
        }
        beta = ( p00 - p01 ) / ( p00 - 2. * p01 + p11 );
        if ( beta < 1.e-3 )
            beta = 1.e-3;
        if ( beta > 1.5 )
            beta = 1.5;
        betam1 = 1. - beta;
        for ( i = 0; i < no2; ++i ) {
            p2[i] = p1[i];
            p1[i] = p0[i];
            p0[i] = p1[i] * beta + p2[i] * betam1;
        }
        return;
    }
    for ( i = 0; i < no2; ++i ) {
        p2[i] = p1[i];
        p1[i] = p0[i];
        p0[i] = 0.5 * ( p1[i] + p2[i] );
    }
}



double
calc_energy ( int no, const double *pm, const double *hm, const double *gm )
{
    register int i, j, ij;
    register double sum1, sum2, pij;
    sum1 = sum2 = 0.0;
    ij = 0;
    for ( i = 0; i < no; ++i ) {
        for ( j = 0; j < i; ++j, ++ij ) {
            pij = pm[ij];
            sum1 += hm[ij] * pij;
            sum2 += gm[ij] * pij;
        }
        pij = pm[ij];
        sum1 += pij * hm[ij] * 0.5;
        sum2 += pij * gm[ij] * 0.5;
        ++ij;
    }
    return ( sum1 + sum1 + sum2 );
}

double
pdiffNorm ( int no, const double *p0, const double *p1 )
{
    int i, j, ij;
    double t, sum;
    sum = 0.0;
    ij = 0;
    for ( i = 0; i < no; ++i, ++ij ) {
        for ( j = 0; j < i; ++j, ++ij ) {
            t = p0[ij] - p1[ij];
            sum += t * t;
        }
        t = p0[ij] - p1[ij];
        sum += t * t * 0.5;
    }
    return ( sqrt ( sum * 2. ) / no );
}


void
rhf_analyzeMoments ( int norb, const double *pmat, double *cmat )
{
    int ij, ir, jr, ijr, icen;
    double pij, f, qi;
    double edx, edy, edz, eqxx, eqxy, eqxz, eqyy, eqyz, eqzz;
    double tdx, tdy, tdz, tqxx, tqxy, tqxz, tqyy, tqyz, tqzz;
    double **edints, **eqints, *work;
    const double *ri;
    double ndints[3], nqints[6];
    FILE *infile, *outfile;
    int ncen, no2;
    const Center *center;
    ncen = number_of_centers ();
    center = centers ();
    no2 = ( ( norb + 1 ) * norb ) / 2;
    edints = dmatrix ( 3, no2 );
    eqints = dmatrix ( 6, no2 );
    infile = openFile ( "MOMINTS.DAT", "r" );
    edx = edy = edz = eqxx = eqxy = eqxz = eqyy = eqyz = eqzz = 0.0;
    for ( ij = 0; ij < no2; ++ij ) {
        fread ( &ir, sizeof ( int ), 1, infile );
        fread ( &jr, sizeof ( int ), 1, infile );
        fread ( ndints, sizeof ( double ), 3, infile );
        fread ( nqints, sizeof ( double ), 6, infile );
        ijr = ir * ( ir + 1 ) / 2 + jr;
        pij = pmat[ijr];
        edints[0][ijr] = pij * ndints[0];
        edints[1][ijr] = pij * ndints[1];
        edints[2][ijr] = pij * ndints[2];
        eqints[0][ijr] = pij * nqints[0];
        eqints[1][ijr] = pij * nqints[1];
        eqints[2][ijr] = pij * nqints[2];
        eqints[3][ijr] = pij * nqints[3];
        eqints[4][ijr] = pij * nqints[4];
        eqints[5][ijr] = pij * nqints[5];
        f = 4.;
        if ( ir == jr )
            f = 2.;
        pij *= f;
        edx += pij * ndints[0];
        edy += pij * ndints[1];
        edz += pij * ndints[2];
        eqxx += pij * nqints[0];
        eqxy += pij * nqints[1];
        eqxz += pij * nqints[2];
        eqyy += pij * nqints[3];
        eqyz += pij * nqints[4];
        eqzz += pij * nqints[5];
    }
    fclose ( infile );
    ndints[0] = 0.0;
    ndints[1] = 0.0;
    ndints[2] = 0.0;
    nqints[0] = 0.0;
    nqints[1] = 0.0;
    nqints[2] = 0.0;
    nqints[3] = 0.0;
    nqints[4] = 0.0;
    nqints[5] = 0.0;
    for ( icen = 0; icen < ncen; ++icen ) {
        qi = ( center + icen )->q;
        ri = ( center + icen )->r;
        ndints[0] += qi * ri[0];
        ndints[1] += qi * ri[1];
        ndints[2] += qi * ri[2];
        nqints[0] += qi * ri[0] * ri[0];
        nqints[1] += qi * ri[0] * ri[1];
        nqints[2] += qi * ri[0] * ri[2];
        nqints[3] += qi * ri[1] * ri[1];
        nqints[4] += qi * ri[1] * ri[2];
        nqints[5] += qi * ri[2] * ri[2];
    }
    tdx = ndints[0] - edx;
    tdy = ndints[1] - edy;
    tdz = ndints[2] - edz;
    tqxx = nqints[0] - eqxx;
    tqxy = nqints[1] - eqxy;
    tqxz = nqints[2] - eqxz;
    tqyy = nqints[3] - eqyy;
    tqyz = nqints[4] - eqyz;
    tqzz = nqints[5] - eqzz;
    outfile = openFile ( "moments.out", "w" );
    fprintf ( outfile, "\t\t\tPATMOL MULTIPOLE MOMENTS ANALYSIS\n" );
    fprintf ( outfile, "----------------------------------------------------\n" );
    fprintf ( outfile, "\t\t\tTotal Cartesian Moments\n" );
    fprintf ( outfile, "   Electronic      Nuclear         Total\n" );
    fprintf ( outfile, "x  %15.6le %15.6le %15.6le\n", edx, ndints[0], tdx );
    fprintf ( outfile, "y  %15.6le %15.6le %15.6le\n", edy, ndints[1], tdy );
    fprintf ( outfile, "z  %15.6le %15.6le %15.6le\n", edz, ndints[2], tdz );
    fprintf ( outfile, "xx %15.6le %15.6le %15.6le\n", eqxx, nqints[0], tqxx );
    fprintf ( outfile, "xy %15.6le %15.6le %15.6le\n", eqxy, nqints[1], tqxy );
    fprintf ( outfile, "xy %15.6le %15.6le %15.6le\n", eqxz, nqints[2], tqxz );
    fprintf ( outfile, "yy %15.6le %15.6le %15.6le\n", eqyy, nqints[3], tqyy );
    fprintf ( outfile, "yz %15.6le %15.6le %15.6le\n", eqyz, nqints[4], tqyz );
    fprintf ( outfile, "zz %15.6le %15.6le %15.6le\n", eqzz, nqints[5], tqzz );
    fprintf ( outfile, "\n" );
    fflush ( outfile );
    fclose ( outfile );
    /*
      work = dalloc (norb * norb);
      spack_trans (norb, &(edints[0][0]), cmat, work);
      spack_trans (norb, &(edints[1][0]), cmat, work);
      spack_trans (norb, &(edints[2][0]), cmat, work);
      spack_trans (norb, &(eqints[0][0]), cmat, work);
      spack_trans (norb, &(eqints[1][0]), cmat, work);
      spack_trans (norb, &(eqints[2][0]), cmat, work);
      spack_trans (norb, &(eqints[3][0]), cmat, work);
      spack_trans (norb, &(eqints[4][0]), cmat, work);
      spack_trans (norb, &(eqints[5][0]), cmat, work);
      outfile = createFile ("TMOM.DAT");
      fwrite (edints[0], sizeof (double), no2, outfile);
      fwrite (edints[1], sizeof (double), no2, outfile);
      fwrite (edints[2], sizeof (double), no2, outfile);
      fwrite (eqints[0], sizeof (double), no2, outfile);
      fwrite (eqints[1], sizeof (double), no2, outfile);
      fwrite (eqints[2], sizeof (double), no2, outfile);
      fwrite (eqints[3], sizeof (double), no2, outfile);
      fwrite (eqints[4], sizeof (double), no2, outfile);
      fwrite (eqints[5], sizeof (double), no2, outfile);
      fclose (outfile);
      free (work);
    */
    free_dmatrix ( eqints, 6 );
    free_dmatrix ( edints, 3 );
}

void
uhf_analyzeMoments ( int norb,
                     const double *pmatA, const double *pmatB,
                     double *cmatA, double *cmatB )
{
    int ij, ir, jr, ijr, icen, ncen, n2;
    double pij, pijA, pijB, f, qi;
    double edx, edy, edz, eqxx, eqxy, eqxz, eqyy, eqyz, eqzz;
    double tdx, tdy, tdz, tqxx, tqxy, tqxz, tqyy, tqyz, tqzz;
    double **edintsA, **eqintsA, *work;
    double **edintsB, **eqintsB;
    const double *ri;
    double ndints[3], nqints[6];
    FILE *infile, *outfile;
    const Center *center;
    n2 = ( norb * ( norb + 1 ) ) / 2;
    ncen = number_of_centers ();
    center = centers ();
    norb = number_of_orbitals ();
    edintsA = dmatrix ( 3, n2 );
    eqintsA = dmatrix ( 6, n2 );
    edintsB = dmatrix ( 3, n2 );
    eqintsB = dmatrix ( 6, n2 );
    infile = openFile ( "MOMINTS.DAT", "r" );
    edx = edy = edz = eqxx = eqxy = eqxz = eqyy = eqyz = eqzz = 0.0;
    for ( ij = 0; ij < n2; ++ij ) {
        fread ( &ir, sizeof ( int ), 1, infile );
        fread ( &jr, sizeof ( int ), 1, infile );
        fread ( ndints, sizeof ( double ), 3, infile );
        fread ( nqints, sizeof ( double ), 6, infile );
        ijr = ir * ( ir + 1 ) / 2 + jr;
        pijA = pmatA[ijr];
        pijB = pmatB[ijr];
        edintsA[0][ijr] = pijA * ndints[0];
        edintsA[1][ijr] = pijA * ndints[1];
        edintsA[2][ijr] = pijA * ndints[2];
        eqintsA[0][ijr] = pijA * nqints[0];
        eqintsA[1][ijr] = pijA * nqints[1];
        eqintsA[2][ijr] = pijA * nqints[2];
        eqintsA[3][ijr] = pijA * nqints[3];
        eqintsA[4][ijr] = pijA * nqints[4];
        eqintsA[5][ijr] = pijA * nqints[5];
        edintsB[0][ijr] = pijB * ndints[0];
        edintsB[1][ijr] = pijB * ndints[1];
        edintsB[2][ijr] = pijB * ndints[2];
        eqintsB[0][ijr] = pijB * nqints[0];
        eqintsB[1][ijr] = pijB * nqints[1];
        eqintsB[2][ijr] = pijB * nqints[2];
        eqintsB[3][ijr] = pijB * nqints[3];
        eqintsB[4][ijr] = pijB * nqints[4];
        eqintsB[5][ijr] = pijB * nqints[5];
        f = 2.;
        if ( ir == jr )
            f = 1.;
        pij = f * ( pijA + pijB );
        edx += pij * ndints[0];
        edy += pij * ndints[1];
        edz += pij * ndints[2];
        eqxx += pij * nqints[0];
        eqxy += pij * nqints[1];
        eqxz += pij * nqints[2];
        eqyy += pij * nqints[3];
        eqyz += pij * nqints[4];
        eqzz += pij * nqints[5];
    }
    fclose ( infile );
    ndints[0] = 0.0;
    ndints[1] = 0.0;
    ndints[2] = 0.0;
    nqints[0] = 0.0;
    nqints[1] = 0.0;
    nqints[2] = 0.0;
    nqints[3] = 0.0;
    nqints[4] = 0.0;
    nqints[5] = 0.0;
    for ( icen = 0; icen < ncen; ++icen ) {
        qi = ( center + icen )->q;
        ri = ( center + icen )->r;
        ndints[0] += qi * ri[0];
        ndints[1] += qi * ri[1];
        ndints[2] += qi * ri[2];
        nqints[0] += qi * ri[0] * ri[0];
        nqints[1] += qi * ri[0] * ri[1];
        nqints[2] += qi * ri[0] * ri[2];
        nqints[3] += qi * ri[1] * ri[1];
        nqints[4] += qi * ri[1] * ri[2];
        nqints[5] += qi * ri[2] * ri[2];
    }
    tdx = ndints[0] - edx;
    tdy = ndints[1] - edy;
    tdz = ndints[2] - edz;
    tqxx = nqints[0] - eqxx;
    tqxy = nqints[1] - eqxy;
    tqxz = nqints[2] - eqxz;
    tqyy = nqints[3] - eqyy;
    tqyz = nqints[4] - eqyz;
    tqzz = nqints[5] - eqzz;
    outfile = openFile ( "moments.out", "w" );
    fprintf ( outfile, "\t\t\tPATMOL MULTIPOLE MOMENTS ANALYSIS\n" );
    fprintf ( outfile,
              "--------------------------------------------------------\n" );
    fprintf ( outfile, "\t\t\tTotal Cartesian Moments\n" );
    fprintf ( outfile, "   Electronic      Nuclear         Total\n" );
    fprintf ( outfile, "x  %15.6le %15.6le %15.6le\n", edx, ndints[0], tdx );
    fprintf ( outfile, "y  %15.6le %15.6le %15.6le\n", edy, ndints[1], tdy );
    fprintf ( outfile, "z  %15.6le %15.6le %15.6le\n", edz, ndints[2], tdz );
    fprintf ( outfile, "xx %15.6le %15.6le %15.6le\n", eqxx, nqints[0], tqxx );
    fprintf ( outfile, "xy %15.6le %15.6le %15.6le\n", eqxy, nqints[1], tqxy );
    fprintf ( outfile, "xy %15.6le %15.6le %15.6le\n", eqxz, nqints[2], tqxz );
    fprintf ( outfile, "yy %15.6le %15.6le %15.6le\n", eqyy, nqints[3], tqyy );
    fprintf ( outfile, "yz %15.6le %15.6le %15.6le\n", eqyz, nqints[4], tqyz );
    fprintf ( outfile, "zz %15.6le %15.6le %15.6le\n", eqzz, nqints[5], tqzz );
    fprintf ( outfile, "\n" );
    fclose ( outfile );
    work = dalloc ( norb * norb );
    spack_trans ( norb, & ( edintsA[0][0] ), cmatA, work );
    spack_trans ( norb, & ( edintsA[1][0] ), cmatA, work );
    spack_trans ( norb, & ( edintsA[2][0] ), cmatA, work );
    spack_trans ( norb, & ( eqintsA[0][0] ), cmatA, work );
    spack_trans ( norb, & ( eqintsA[1][0] ), cmatA, work );
    spack_trans ( norb, & ( eqintsA[2][0] ), cmatA, work );
    spack_trans ( norb, & ( eqintsA[3][0] ), cmatA, work );
    spack_trans ( norb, & ( eqintsA[4][0] ), cmatA, work );
    spack_trans ( norb, & ( eqintsA[5][0] ), cmatA, work );
    spack_trans ( norb, & ( edintsB[0][0] ), cmatB, work );
    spack_trans ( norb, & ( edintsB[1][0] ), cmatB, work );
    spack_trans ( norb, & ( edintsB[2][0] ), cmatB, work );
    spack_trans ( norb, & ( eqintsB[0][0] ), cmatB, work );
    spack_trans ( norb, & ( eqintsB[1][0] ), cmatB, work );
    spack_trans ( norb, & ( eqintsB[2][0] ), cmatB, work );
    spack_trans ( norb, & ( eqintsB[3][0] ), cmatB, work );
    spack_trans ( norb, & ( eqintsB[4][0] ), cmatB, work );
    spack_trans ( norb, & ( eqintsB[5][0] ), cmatB, work );
    outfile = openFile ( "TMOM.DAT", "w" );
    fwrite ( edintsA[0], sizeof ( double ), n2, outfile );
    fwrite ( edintsA[1], sizeof ( double ), n2, outfile );
    fwrite ( edintsA[2], sizeof ( double ), n2, outfile );
    fwrite ( eqintsA[0], sizeof ( double ), n2, outfile );
    fwrite ( eqintsA[1], sizeof ( double ), n2, outfile );
    fwrite ( eqintsA[2], sizeof ( double ), n2, outfile );
    fwrite ( eqintsA[3], sizeof ( double ), n2, outfile );
    fwrite ( eqintsA[4], sizeof ( double ), n2, outfile );
    fwrite ( eqintsA[5], sizeof ( double ), n2, outfile );
    fwrite ( edintsB[0], sizeof ( double ), n2, outfile );
    fwrite ( edintsB[1], sizeof ( double ), n2, outfile );
    fwrite ( edintsB[2], sizeof ( double ), n2, outfile );
    fwrite ( eqintsB[0], sizeof ( double ), n2, outfile );
    fwrite ( eqintsB[1], sizeof ( double ), n2, outfile );
    fwrite ( eqintsB[2], sizeof ( double ), n2, outfile );
    fwrite ( eqintsB[3], sizeof ( double ), n2, outfile );
    fwrite ( eqintsB[4], sizeof ( double ), n2, outfile );
    fwrite ( eqintsB[5], sizeof ( double ), n2, outfile );
    fclose ( outfile );
    free ( work );
    free_dmatrix ( eqintsB, 6 );
    free_dmatrix ( edintsB, 3 );
    free_dmatrix ( eqintsA, 6 );
    free_dmatrix ( edintsA, 3 );
}



void
rhf_form_gmatrix ( const ERIFS * tfiles, const double *Pmat, double *Gmat )
{
    register int ix, i, j, k, l, ij, ik, il, jk, jl, kl;
    register double val, da, db, sik, sil, sjk, sjl;
    int ifile, nints, ibin, nbin, nex, ii;
    const int nfiles = tfiles->nfiles;
    const int rank = tfiles->rank;
    const char *basename = tfiles->basename;
    ERINT sints[BINSIZE];
    FILE *in;
    for ( ifile = 0; ifile < nfiles; ++ifile ) {
        nints = ( tfiles->numints ) [ifile];
        nbin = nints / BINSIZE;
        nex = nints % BINSIZE;
        in = openERIFile ( rank, ifile, basename );
        if ( nex ) {
            fread ( sints, sizeof ( ERINT ), nex, in );
            for ( ix = 0; ix < nex; ++ix ) {
                val = ( sints + ix )->val;
                i = ( sints + ix )->i;
                j = ( sints + ix )->j;
                k = ( sints + ix )->k;
                l = ( sints + ix )->l;
                ii = ( i * ( i + 1 ) ) / 2;
                ij = ii + j;
                ik = ii + k;
                il = ii + l;
                kl = ( k * ( k + 1 ) ) / 2 + l;
                if ( j > k ) {
                    jk = ( j * ( j + 1 ) ) / 2 + k;
                } else {
                    jk = ( k * ( k + 1 ) ) / 2 + j;
                }
                if ( j > l ) {
                    jl = ( j * ( j + 1 ) ) / 2 + l;
                } else {
                    jl = ( l * ( l + 1 ) ) / 2 + j;
                }
                da = val * 2.0 * Pmat[ij];
                db = val * 2.0 * Pmat[kl];
                sjl = val * Pmat[ik];
                sjk = val * Pmat[il];
                sik = val * Pmat[jl];
                sil = val * Pmat[jk];
                if ( k != l ) {
                    db = db + db;
                    Gmat[ik] -= sik;
                    if ( i != j && j >= k )
                        Gmat[jk] -= sjk;
                }
                Gmat[il] -= sil;
                Gmat[ij] += db;
                if ( i != j && j >= l )
                    Gmat[jl] -= sjl;
                if ( ij != kl ) {
                    if ( i != j )
                        da = da + da;
                    if ( j <= k ) {
                        Gmat[jk] -= sjk;
                        if ( i != j && i <= k )
                            Gmat[ik] -= sik;
                        if ( k != l && j <= l )
                            Gmat[jl] -= sjl;
                    }
                    Gmat[kl] += da;
                }
            }
        }
        for ( ibin = 0; ibin < nbin; ++ibin ) {
            fread ( sints, sizeof ( ERINT ), BINSIZE, in );
            for ( ix = 0; ix < BINSIZE; ++ix ) {
                val = ( sints + ix )->val;
                i = ( sints + ix )->i;
                j = ( sints + ix )->j;
                k = ( sints + ix )->k;
                l = ( sints + ix )->l;
                ii = ( i * ( i + 1 ) ) / 2;
                ij = ii + j;
                ik = ii + k;
                il = ii + l;
                kl = ( k * ( k + 1 ) ) / 2 + l;
                if ( j > k ) {
                    jk = ( j * ( j + 1 ) ) / 2 + k;
                } else {
                    jk = ( k * ( k + 1 ) ) / 2 + j;
                }
                if ( j > l ) {
                    jl = ( j * ( j + 1 ) ) / 2 + l;
                } else {
                    jl = ( l * ( l + 1 ) ) / 2 + j;
                }
                da = val * 2.0 * Pmat[ij];
                db = val * 2.0 * Pmat[kl];
                sjl = val * Pmat[ik];
                sjk = val * Pmat[il];
                sik = val * Pmat[jl];
                sil = val * Pmat[jk];
                if ( k != l ) {
                    db = db + db;
                    Gmat[ik] -= sik;
                    if ( i != j && j >= k )
                        Gmat[jk] -= sjk;
                }
                Gmat[il] -= sil;
                Gmat[ij] += db;
                if ( i != j && j >= l )
                    Gmat[jl] -= sjl;
                if ( ij != kl ) {
                    if ( i != j )
                        da = da + da;
                    if ( j <= k ) {
                        Gmat[jk] -= sjk;
                        if ( i != j && i <= k )
                            Gmat[ik] -= sik;
                        if ( k != l && j <= l )
                            Gmat[jl] -= sjl;
                    }
                    Gmat[kl] += da;
                }
            }
        }
        fclose ( in );
    }
}

void
uhf_form_gmatrix ( const ERIFS * tfiles,
                   const double *PmatA, const double *PmatB,
                   double *GmatA, double *GmatB )
{
    register int ix, i, j, k, l, ij, ik, il, jk, jl, kl;
    register double val, da, db, sikA, silA, sjkA, sjlA, sikB, silB, sjkB, sjlB;
    int ifile, nints, ibin, nbin, nex, ii;
    const int nfiles = tfiles->nfiles;
    const int rank = tfiles->rank;
    const char *basename = tfiles->basename;
    ERINT sints[BINSIZE];
    FILE *in;
    for ( ifile = 0; ifile < nfiles; ++ifile ) {
        nints = ( tfiles->numints ) [ifile];
        nbin = nints / BINSIZE;
        nex = nints - BINSIZE * nbin;
        in = openERIFile ( rank, ifile, basename );
        if ( nex ) {
            fread ( sints, sizeof ( ERINT ), nex, in );
            for ( ix = 0; ix < nex; ++ix ) {
                val = ( sints + ix )->val;
                i = ( sints + ix )->i;
                j = ( sints + ix )->j;
                k = ( sints + ix )->k;
                l = ( sints + ix )->l;
                ii = ( i * ( i + 1 ) ) / 2;
                ij = ii + j;
                ik = ii + k;
                il = ii + l;
                kl = ( k * ( k + 1 ) ) / 2 + l;
                if ( j > k ) {
                    jk = ( j * ( j + 1 ) ) / 2 + k;
                } else {
                    jk = ( k * ( k + 1 ) ) / 2 + j;
                }
                if ( j > l ) {
                    jl = ( j * ( j + 1 ) ) / 2 + l;
                } else {
                    jl = ( l * ( l + 1 ) ) / 2 + j;
                }
                da = val * ( PmatA[ij] + PmatB[ij] );
                db = val * ( PmatA[kl] + PmatB[kl] );
                sjlA = val * PmatA[ik];
                sjkA = val * PmatA[il];
                sikA = val * PmatA[jl];
                silA = val * PmatA[jk];
                sjlB = val * PmatB[ik];
                sjkB = val * PmatB[il];
                sikB = val * PmatB[jl];
                silB = val * PmatB[jk];
                if ( k != l ) {
                    db = db + db;
                    GmatA[ik] -= sikA;
                    GmatB[ik] -= sikB;
                    if ( i != j && j >= k ) {
                        GmatA[jk] -= sjkA;
                        GmatB[jk] -= sjkB;
                    }
                }
                GmatA[il] -= silA;
                GmatA[ij] += db;
                GmatB[il] -= silB;
                GmatB[ij] += db;
                if ( i != j && j >= l ) {
                    GmatA[jl] -= sjlA;
                    GmatB[jl] -= sjlB;
                }
                if ( ij != kl ) {
                    if ( i != j )
                        da = da + da;
                    if ( j <= k ) {
                        GmatA[jk] -= sjkA;
                        GmatB[jk] -= sjkB;
                        if ( i != j && i <= k ) {
                            GmatA[ik] -= sikA;
                            GmatB[ik] -= sikB;
                        }
                        if ( k != l && j <= l ) {
                            GmatA[jl] -= sjlA;
                            GmatB[jl] -= sjlB;
                        }
                    }
                    GmatA[kl] += da;
                    GmatB[kl] += da;
                }
            }
        }
        for ( ibin = 0; ibin < nbin; ++ibin ) {
            fread ( sints, sizeof ( ERINT ), BINSIZE, in );
            for ( ix = 0; ix < BINSIZE; ++ix ) {
                val = ( sints + ix )->val;
                i = ( sints + ix )->i;
                j = ( sints + ix )->j;
                k = ( sints + ix )->k;
                l = ( sints + ix )->l;
                ii = ( i * ( i + 1 ) ) / 2;
                ij = ii + j;
                ik = ii + k;
                il = ii + l;
                kl = ( k * ( k + 1 ) ) / 2 + l;
                if ( j > k ) {
                    jk = ( j * ( j + 1 ) ) / 2 + k;
                } else {
                    jk = ( k * ( k + 1 ) ) / 2 + j;
                }
                if ( j > l ) {
                    jl = ( j * ( j + 1 ) ) / 2 + l;
                } else {
                    jl = ( l * ( l + 1 ) ) / 2 + j;
                }
                da = val * ( PmatA[ij] + PmatB[ij] );
                db = val * ( PmatA[kl] + PmatB[kl] );
                sjlA = val * PmatA[ik];
                sjkA = val * PmatA[il];
                sikA = val * PmatA[jl];
                silA = val * PmatA[jk];
                sjlB = val * PmatB[ik];
                sjkB = val * PmatB[il];
                sikB = val * PmatB[jl];
                silB = val * PmatB[jk];
                if ( k != l ) {
                    db = db + db;
                    GmatA[ik] -= sikA;
                    GmatB[ik] -= sikB;
                    if ( i != j && j >= k ) {
                        GmatA[jk] -= sjkA;
                        GmatB[jk] -= sjkB;
                    }
                }
                GmatA[il] -= silA;
                GmatA[ij] += db;
                GmatB[il] -= silB;
                GmatB[ij] += db;
                if ( i != j && j >= l ) {
                    GmatA[jl] -= sjlA;
                    GmatB[jl] -= sjlB;
                }
                if ( ij != kl ) {
                    if ( i != j )
                        da = da + da;
                    if ( j <= k ) {
                        GmatA[jk] -= sjkA;
                        GmatB[jk] -= sjkB;
                        if ( i != j && i <= k ) {
                            GmatA[ik] -= sikA;
                            GmatB[ik] -= sikB;
                        }
                        if ( k != l && j <= l ) {
                            GmatA[jl] -= sjlA;
                            GmatB[jl] -= sjlB;
                        }
                    }
                    GmatA[kl] += da;
                    GmatB[kl] += da;
                }
            }
        }
        fclose ( in );
    }
}



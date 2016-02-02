#include "basis.h"
#define DPM_NSHELL 7
#define DPM_MAXL 1
#define DPM_NORB 13

static int norb, nshell, ncen;
static int tnorb, tnshell, tncen;
static int dpm, dft, dft_x, dft_c, ffc;
static int scf_c, scf_g, scf_a;
static int prt_e, prt_p, prt_d;
static int nelec, spinm;
static int max_scf_its;
static double scf_eps, ffc_efld;
static int skipcen;
static Shell *shell;
static Center *center;
static int ***lxyz;
static double **normFact;
static int maxlsh, maxlst;

int
number_of_orbitals ()
{
    return norb;
};
int
number_of_shells ()
{
    return nshell;
};
int
number_of_centers ()
{
    return ncen;
};
int
total_number_of_orbitals ()
{
    return norb;
};
int
total_number_of_shells ()
{
    return nshell;
};
int
total_number_of_centers ()
{
    return ncen;
};
int
number_of_electrons ()
{
    return nelec;
};
int
spin_multiplicity ()
{
    return spinm;
};
int
max_scf_iterations ()
{
    return max_scf_its;
};
int
skip_center ()
{
    return skipcen;
};
int
max_lvalue ()
{
    return maxlsh;
};
int
max_lstates ()
{
    return maxlst;
};
int
do_dpm_calculation ()
{
    return dpm;
};
int
do_dft_calculation ()
{
    return dft;
};
int
dft_correlation ()
{
    return dft_c;
};
int
dft_exchange ()
{
    return dft_x;
};
int
do_finite_field_calculation ()
{
    return ffc;
};
int
scf_acceleration ()
{
    return scf_a;
};
int
scf_guess ()
{
    return scf_g;
};
int
print_eigenvectors ()
{
    return prt_e;
};
int
print_density ()
{
    return prt_d;
};
int
print_population ()
{
    return prt_p;
};
double
finite_field_electric_field ()
{
    return 5.e-3;
};
double
scf_convergence ()
{
    return scf_eps;
};
const Shell *
shells ()
{
    return shell;
};
const Center *
centers ()
{
    return center;
};
int ***
lvector_components ()
{
    return lxyz;
};
double **
normalizationFactors ()
{
    return normFact;
};

void
set_center_position ( int icen, const double *Rnew )
{
    memcpy ( ( center + icen )->r, Rnew, sizeof ( double ) * 3 );
}

void
set_center_charge ( int icen, double qnew )
{
    ( center + icen )->q = qnew;
}

double
nuclear_repulsion ()
{
    int i, j;
    double sum, qi, qj, x, y, z, r2;
    const double *ri, *rj;
    sum = 0.0;
    for ( i = 0; i < ncen; ++i ) {
        qi = ( center + i )->q;
        ri = ( center + i )->r;
        for ( j = i + 1; j < ncen; ++j ) {
            qj = ( center + j )->q;
            rj = ( center + j )->r;
            x = ri[0] - rj[0];
            y = ri[1] - rj[1];
            z = ri[2] - rj[2];
            r2 = x * x + y * y + z * z;
            sum += qi * qj / sqrt ( r2 );
        }
    }
    return sum;
}

int
row ( const Center * cen )
{
    int iz;
    double dz;
    dz = cen->q;
    iz = ( int ) rint ( dz );
    if ( iz <= 2 )
        return 1;
    if ( iz <= 10 )
        return 2;
    if ( iz <= 18 )
        return 3;
    if ( iz <= 36 )
        return 4;
    if ( iz <= 54 )
        return 5;
    if ( iz <= 86 )
        return 6;
    return 7;
}


void
assign_dpm_shells ( Shell * sh, int poscen, int off )
{
    int ish;
    const double alf[] =
    { 13.3615, 2.0133, 0.453757, 1.2331, 2.5719, 0.59976, 0.139365 };
    const size_t lsh[] = { 0, 0, 0, 0, 1, 1, 1 };
    const size_t nls[] = { 1, 1, 1, 1, 3, 3, 3 };
    int ix;
    for ( ish = 0; ish < DPM_NSHELL; ++ish ) {
        ix = ish + nshell;
        ( shell + ix )->npr = 1;
        ( shell + ix )->lsh = lsh[ish];
        ( shell + ix )->cen = poscen;
        ( shell + ix )->off = off;
        off += nls[ish];
        ( shell + ix )->al = dalloc ( 1 );
        ( shell + ix )->co = dalloc ( 1 );
        ( shell + ix )->al[0] = alf[ish];
        ( shell + ix )->co[0] = 1.0;
    }
}

void
init_Basis ()
{
    int icen, ish, l, lx, ly, off_, lst, maxlp1;
    double xt, yt, zt;
    const double dfact[] = { 1., 1., 3., 15., 105., 945., 10395. };
    int *nlstates;
    FILE *infile;
    infile = openFile ( "patin.dat", "r" );
    fscanf ( infile, "%d %d %d", &nshell, &norb, &ncen );
    fscanf ( infile, "%d %le", &max_scf_its, &scf_eps );
    fscanf ( infile, "%d %d", &nelec, &spinm );
    fscanf ( infile, "%d %d %d %d", &dpm, &dft, &dft_x, &dft_x );
    fscanf ( infile, "%d %d %d", &ffc, &scf_a, &scf_g );
    fscanf ( infile, "%d %d %d", &prt_e, &prt_p, &prt_d );
    tncen = ncen;
    tnshell = nshell;
    tnorb = norb;
    if ( dpm ) {
        ++tncen;
        tnshell += DPM_NSHELL;
        tnorb += DPM_NORB;
    }
    skipcen = ncen;
    center = ( Center * ) pmalloc ( sizeof ( Center ) * tncen );
    shell = ( Shell * ) pmalloc ( sizeof ( Shell ) * tnshell );
    for ( icen = 0; icen < ncen; ++icen ) {
        read_center ( infile, center + icen );
    }
    off_ = 0;
    maxlsh = 0;
    for ( ish = 0; ish < nshell; ++ish ) {
        read_shell ( infile, shell + ish, &off_ );
        if ( ( ( shell + ish )->lsh ) > maxlsh ) {
            maxlsh = ( shell + ish )->lsh;
        }
    }
    if ( dpm ) {
        maxlsh = ( maxlsh < DPM_MAXL ) ? DPM_MAXL : maxlsh;
        assign_dpm_shells ( shell + nshell, ncen, off_ );
        ( center + ncen )->q = 1.0;
    }
    maxlst = ( maxlsh + 1 ) * ( maxlsh + 2 ) / 2;
    /*  for (ish = 0; ish < tnshell; ++ish)
        normalize_shell (shell + ish);
    */
    maxlp1 = maxlsh + 1;
    nlstates = ialloc ( maxlp1 );
    for ( l = 0; l <= maxlsh; ++l ) {
        nlstates[l] = ( l + 1 ) * ( l + 2 ) / 2;
    }
    lxyz = ( int *** ) pmalloc ( sizeof ( int ** ) * maxlp1 );
    for ( l = 0; l <= maxlsh; ++l ) {
        lxyz[l] = imatrix ( nlstates[l], 3 );
        lst = 0;
        lx = l + 1;
        for ( ; lx--; ) {
            ly = l - lx + 1;
            for ( ; ly--; ) {
                lxyz[l][lst][0] = lx;
                lxyz[l][lst][1] = ly;
                lxyz[l][lst][2] = l - lx - ly;
                ++lst;
            }
        }
    }
    normFact = ( double ** ) pmalloc ( sizeof ( double * ) * maxlp1 );
    for ( l = 0; l <= maxlsh; ++l ) {
        normFact[l] = dalloc ( nlstates[l] );
        lx = nlstates[l];
        for ( lst = 0; lst < lx; ++lst ) {
            xt = dfact[lxyz[l][lst][0]];
            yt = dfact[lxyz[l][lst][1]];
            zt = dfact[lxyz[l][lst][2]];
            normFact[l][lst] = 1. / sqrt ( xt * yt * zt );
        }
    }
    free ( nlstates );
}

void
dpm_augment ()
{
    norb = tnorb;
    nshell = tnshell;
    ncen = tncen;
}

void
basis_output ( double *smat, double *tmat, double *hmat )
{
    int i, j, ij;
    FILE *out = openFile ( "intsout.dat", "w" );
    fprintf ( out, "  patmol alpha version 2.5 \n" );
    fprintf ( out, "\n" );
    fprintf ( stderr, "Using " );
    switch ( scf_a ) {
    case 0:
        fprintf ( stderr, " No SCF Acceleration\n" );
        break;
    case 1:
        fprintf ( stderr, "Simple Averaging scf\n" );
        break;
    default:
        fprintf ( stderr, "Anderson Averaging scf\n" );
        break;
    }
    fprintf ( out, "Number of Shells   = %12u\n", nshell );
    fprintf ( out, "Number of Orbitals = %12u\n", norb );
    fprintf ( out, "Number of Centers  = %12u\n", ncen );
    fprintf ( out, "Max L value        = %12u\n", maxlsh );
    fprintf ( out, "Max SCF Iterations = %12u\n", max_scf_its );
    fprintf ( out, "Number of Electrons= %12u\n", nelec );
    fprintf ( out, " Centers \n" );
    fprintf ( out, "Charge     Position(x,y,z  in bohr)\n" );
    for ( i = 0; i < ncen; ++i ) {
        write_center ( out, center + i );
    }
    fprintf ( out, "\n" );
    fprintf ( out, " Shells\n" );
    for ( i = 0; i < nshell; ++i ) {
        write_shell ( out, shell + i );
    }
    fprintf ( out, "\n" );
    fprintf ( out, "Overlap Matrix \n" );
    ij = 0;
    for ( i = 0; i < norb; ++i ) {
        for ( j = 0; j <= i; ++j ) {
            fprintf ( out, "%7u %7u %24.16le\n", i, j, smat[ij] );
            ++ij;
        }
    }
    fprintf ( out, "\n" );
    fprintf ( out, "Kinetic Energy Matrix \n" );
    ij = 0;
    for ( i = 0; i < norb; ++i ) {
        for ( j = 0; j <= i; ++j ) {
            fprintf ( out, "%7u %7u %24.16le\n", i, j, tmat[ij] );
            ++ij;
        }
    }
    fprintf ( out, "\n" );
    fprintf ( out, "Core Hamiltonian Matrix \n" );
    ij = 0;
    for ( i = 0; i < norb; ++i ) {
        for ( j = 0; j <= i; ++j ) {
            fprintf ( out, "%7u %7u %24.16le\n", i, j, hmat[ij] );
            ++ij;
        }
    }
    fclose ( out );
}

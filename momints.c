#include "momints.h"

static double absq;
static double *nstates;
static int *lstates;
static int *ostates;
static double ***dx, ***dy, ***dz;
static const double *alf1, *cof1, *alf2, *cof2;
static const double *a, *b;
static int lv1, lv2, npr1, npr2;
static int ***lxyz;

inline void
calc_pmom ( int nlst, double **dints, double **qints )
{
    int kc, ipr, jpr, ls1, ls2;
    int l1, m1, n1, l2, m2, n2, lvt, lv1p2, lv2p2;
    double axp, c1, c12, bxp, pxp, abi, s12, p[3];
    double ovl, nf, sov, pvx, pvy, pvz;
    const int *lst;
    const double *dxp, *dyp, *dzp, *nfc;
    const double PITERM = 5.5683279968317079;
    lvt = lv1 + lv2;
    lv1p2 = lv1 + 2;
    lv2p2 = lv2 + 2;
    for ( ipr = 0; ipr < npr1; ++ipr ) {
        axp = alf1[ipr];
        c1 = cof1[ipr];
        for ( jpr = 0; jpr < npr2; ++jpr ) {
            c12 = c1 * cof2[jpr];
            bxp = alf2[jpr];
            pxp = axp + bxp;
            abi = 1. / pxp;
            p[0] = ( axp * a[0] + bxp * b[0] ) * abi;
            p[1] = ( axp * a[1] + bxp * b[1] ) * abi;
            p[2] = ( axp * a[2] + bxp * b[2] ) * abi;
            s12 = PITERM * exp ( -abi * axp * bxp * absq ) * abi * sqrt ( abi );
            abi *= 0.5;
            formMDD ( dx, abi, ( p[0] - a[0] ), ( p[0] - b[0] ), lv1p2, lv2p2 );
            formMDD ( dy, abi, ( p[1] - a[1] ), ( p[1] - b[1] ), lv1p2, lv2p2 );
            formMDD ( dz, abi, ( p[2] - a[2] ), ( p[2] - b[2] ), lv1p2, lv2p2 );
            lst = lstates;
            for ( kc = 0; kc < nlst; ++kc ) {
                ovl = s12 * c12 * nstates[kc];
                ls1 = ( *lst++ );
                ls2 = ( *lst++ );
                l1 = lxyz[lv1][ls1][0];
                m1 = lxyz[lv1][ls1][1];
                n1 = lxyz[lv1][ls1][2];
                l2 = lxyz[lv2][ls2][0];
                m2 = lxyz[lv2][ls2][1];
                n2 = lxyz[lv2][ls2][2];
                dxp = dx[l1][l2];
                dyp = dy[m1][m2];
                dzp = dz[n1][n2];
                sov = ovl * dxp[0] * dyp[0] * dzp[0];
                pvx = ovl * dyp[0] * dzp[0];
                pvy = ovl * dxp[0] * dzp[0];
                pvz = ovl * dxp[0] * dyp[0];
                /* dipole moments */
                dints[kc][0] += sov * p[0] + dxp[1] * pvx;
                dints[kc][1] += sov * p[1] + dyp[1] * pvy;
                dints[kc][2] += sov * p[2] + dzp[1] * pvz;
                /* quadrupole moments */
                qints[kc][0] += sov * ( p[0] * p[0] + abi ) +
                                2. * pvx * ( dxp[2] + p[0] * dxp[1] );
                qints[kc][1] += dzp[0] * ovl * ( dxp[1] * dyp[1] +
                                                 p[0] * dxp[0] * dyp[1] +
                                                 p[1] * dxp[1] * dyp[0] +
                                                 p[0] * p[1] * dxp[0] * dyp[0] );
                qints[kc][2] += dyp[0] * ovl * ( dxp[1] * dzp[1] +
                                                 p[0] * dxp[0] * dzp[1] +
                                                 p[2] * dxp[1] * dzp[0] +
                                                 p[0] * p[2] * dxp[0] * dzp[0] );
                qints[kc][3] += sov * ( p[1] * p[1] + abi ) +
                                2. * pvy * ( dyp[2] + p[1] * dyp[1] );
                qints[kc][4] += dxp[0] * ovl * ( dyp[1] * dzp[1] +
                                                 p[1] * dyp[0] * dzp[1] +
                                                 p[2] * dyp[1] * dzp[0] +
                                                 p[1] * p[2] * dyp[0] * dzp[0] );
                qints[kc][5] += sov * ( p[2] * p[2] + abi ) +
                                2. * pvz * ( dzp[2] + p[2] * dzp[1] );
            }
        }
    }
}

void
momints ()
{
    int ils, jls, ish, jsh, ir, jr, knt;
    int cen1, off1, cen2, off2, nls1, nls2;
    int *ost, *lst;
    double *nst, **dints, **qints;
    int dsize, asize;
    const Shell *sh1;
    const Shell *sh2;
    const Shell *shell = shells ();
    const Center *center = centers ();
    const int nshell = number_of_shells ();
    const int maxlsh = max_lvalue ();
    double **normFactors;
    FILE *outfile;
    lxyz = lvector_components ();
    normFactors = normalizationFactors ();
    dsize = maxlsh + 3;
    asize = maxlsh + maxlsh + 7;
    dx = dtensor ( dsize, dsize, asize );
    dy = dtensor ( dsize, dsize, asize );
    dz = dtensor ( dsize, dsize, asize );
    asize = max_lstates ();
    asize = asize * asize;
    dints = dmatrix ( asize, 3 );
    qints = dmatrix ( asize, 6 );
    nstates = dalloc ( asize );
    lstates = ialloc ( asize * 2 );
    ostates = ialloc ( asize * 2 );
    outfile = openFile ( "MOMINTS.DAT", "w" );
    for ( ish = 0; ish < nshell; ++ish ) {
        sh1 = ( shell + ish );
        npr1 = sh1->npr;
        lv1 = sh1->lsh;
        cen1 = sh1->cen;
        off1 = sh1->off;
        nls1 = sh1->nst;
        alf1 = sh1->al;
        cof1 = sh1->co;
        a = ( center + cen1 )->r;
        for ( jsh = 0; jsh <= ish; ++jsh ) {
            sh2 = ( shell + jsh );
            npr2 = sh2->npr;
            lv2 = sh2->lsh;
            cen2 = sh2->cen;
            off2 = sh2->off;
            nls2 = sh2->nst;
            alf2 = sh2->al;
            cof2 = sh2->co;
            b = ( center + cen2 )->r;
            absq = distSqr ( a, b );
            knt = 0;
            ost = ostates;
            lst = lstates;
            for ( ils = 0; ils < nls1; ++ils ) {
                ir = off1 + ils;
                for ( jls = 0; jls < nls2; ++jls ) {
                    jr = off2 + jls;
                    if ( jr > ir )
                        break;
                    ( *ost++ ) = ir;
                    ( *ost++ ) = jr;
                    ( *lst++ ) = ils;
                    ( *lst++ ) = jls;
                    nstates[knt] = normFactors[lv1][ils] *
                                   normFactors[lv2][jls];
                    dints[knt][0] = 0.0;
                    dints[knt][1] = 0.0;
                    dints[knt][2] = 0.0;
                    qints[knt][0] = 0.0;
                    qints[knt][1] = 0.0;
                    qints[knt][2] = 0.0;
                    qints[knt][3] = 0.0;
                    qints[knt][4] = 0.0;
                    qints[knt][5] = 0.0;
                    ++knt;
                }
            }
            calc_pmom ( knt, dints, qints );
            ost = ostates;
            for ( ils = 0; ils < knt; ++ils ) {
                ir = ( *ost++ );
                jr = ( *ost++ );
                fwrite ( &ir, sizeof ( int ), 1, outfile );
                fwrite ( &jr, sizeof ( int ), 1, outfile );
                fwrite ( dints[ils], sizeof ( double ), 3, outfile );
                fwrite ( qints[ils], sizeof ( double ), 6, outfile );
            }
        }
    }
    fclose ( outfile );
    free ( ostates );
    free ( lstates );
    free ( nstates );
    free_dmatrix ( qints, asize );
    free_dmatrix ( dints, asize );
    free_dtensor ( dz, dsize, dsize );
    free_dtensor ( dy, dsize, dsize );
    free_dtensor ( dx, dsize, dsize );
}


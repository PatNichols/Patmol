#include "oneints.h"
static double absq;
static double *nrmstates;
static int *lstates;
static double ***dx, ***dy, ***dz, ***rsum, ***rterm;
static const double *alf1, *cof1, *alf2, *cof2, *a, *b;
static int lv1, npr1, lv2, npr2;
static int ***lxyz;
static int skipcen, ncen;
static const Center *center;

static inline void
calc_poeints ( int nlst, double *svals, double *tvals, double *vvals )
{
    register int ix, iy, iz, kc, ls1, ls2;
    register double psum, sum;
    int l1, m1, n1, l2, m2, n2, l12, m12, n12;
    int ip, jp, ic;
    double axp, c1, c12, bxp, pxp, abi, s12, pc2, tx, ty, tz, qc, sr, nfct;
    const double *rc, *rp, *dxp, *dyp, *dzp, *nrm;
    const int *lst;
    double *rs;
    const double piterm = 5.5683279968317079;
    double p[3];
    double pc[3];
    const int lvt = lv1 + lv2;
    for ( ip = 0; ip < npr1; ip++ ) {
        axp = alf1[ip];
        c1 = cof1[ip];
        for ( jp = 0; jp < npr2; jp++ ) {
            c12 = c1 * cof2[jp];
            bxp = alf2[jp];
            pxp = axp + bxp;
            abi = 1.0 / pxp;
            s12 = piterm * exp ( -axp * bxp * abi * absq ) / ( pxp * sqrt ( pxp ) );
            p[0] = ( axp * a[0] + bxp * b[0] ) * abi;
            p[1] = ( axp * a[1] + bxp * b[1] ) * abi;
            p[2] = ( axp * a[2] + bxp * b[2] ) * abi;
            abi = abi * 0.5;
            formMDD ( dx, abi, p[0] - a[0], p[0] - b[0], lv1 + 1, lv2 + 1 );
            formMDD ( dy, abi, p[1] - a[1], p[1] - b[1], lv1 + 1, lv2 + 1 );
            formMDD ( dz, abi, p[2] - a[2], p[2] - b[2], lv1 + 1, lv2 + 1 );
            sr = 2.0 * s12 * sqrt ( pxp / M_PI );
            for ( ix = 0; ix <= lvt; ix++ ) {
                for ( iy = 0; iy <= lvt; iy++ ) {
                    for ( iz = 0; iz <= lvt; iz++ ) {
                        rsum[ix][iy][iz] = 0.0;
                    }
                }
            }
            for ( ic = 0; ic < ncen; ic++ ) {
                if ( ic == skipcen )
                    continue;
                qc = ( center + ic )->q;
                rc = ( center + ic )->r;
                pc[0] = p[0] - rc[0];
                pc[1] = p[1] - rc[1];
                pc[2] = p[2] - rc[2];
                pc2 = pc[0] * pc[0] + pc[1] * pc[1] + pc[2] * pc[2];
                formMDR ( rterm, sr, ( pxp * pc2 ), pxp, pc, lvt );
                for ( ix = 0; ix <= lvt; ix++ ) {
                    for ( iy = 0; iy <= lvt; iy++ ) {
                        rs = rsum[ix][iy];
                        rp = rterm[ix][iy];
                        for ( iz = 0; iz <= lvt; iz++ ) {
                            rs[iz] -= qc * rp[iz];
                        }
                    }
                }
            }
            for ( kc = 0; kc < nlst; kc++ ) {
                nfct = c12 * nrmstates[kc];
                ls1 = lstates[kc];
                ls2 = ( ls1 & 255 );
                ls1 >>= 8;
                l1 = lxyz[lv1][ls1][0];
                m1 = lxyz[lv1][ls1][1];
                n1 = lxyz[lv1][ls1][2];
                l2 = lxyz[lv2][ls2][0];
                m2 = lxyz[lv2][ls2][1];
                n2 = lxyz[lv2][ls2][2];
                l12 = l1 + l2;
                m12 = m1 + m2;
                n12 = n1 + n2;
                dxp = dx[l1][l2];
                dyp = dy[m1][m2];
                dzp = dz[n1][n2];
                svals[kc] += nfct * s12 * dxp[0] * dyp[0] * dzp[0];
                tx = 2.0 * axp * bxp * dx[l1 + 1][l2 + 1][0];
                if ( l1 != 0 ) {
                    tx -= l1 * bxp * dx[l1 - 1][l2 + 1][0];
                }
                if ( l2 != 0 ) {
                    tx -= l2 * axp * dx[l1 + 1][l2 - 1][0];
                    if ( l1 != 0 ) {
                        tx += 0.5 * l1 * l2 * dx[l1 - 1][l2 - 1][0];
                    }
                }
                tx *= dyp[0] * dzp[0];
                ty = 2.0 * axp * bxp * dy[m1 + 1][m2 + 1][0];
                if ( m1 != 0 ) {
                    ty -= m1 * bxp * dy[m1 - 1][m2 + 1][0];
                }
                if ( m2 != 0 ) {
                    ty -= m2 * axp * dy[m1 + 1][m2 - 1][0];
                    if ( m1 != 0 ) {
                        ty += 0.5 * m1 * m2 * dy[m1 - 1][m2 - 1][0];
                    }
                }
                ty *= dxp[0] * dzp[0];
                tz = 2.0 * axp * bxp * dz[n1 + 1][n2 + 1][0];
                if ( n1 != 0 ) {
                    tz -= n1 * bxp * dz[n1 - 1][n2 + 1][0];
                }
                if ( n2 != 0 ) {
                    tz -= n2 * axp * dz[n1 + 1][n2 - 1][0];
                    if ( n1 != 0 ) {
                        tz += 0.5 * n1 * n2 * dz[n1 - 1][n2 - 1][0];
                    }
                }
                tz *= dxp[0] * dyp[0];
                tvals[kc] += nfct * s12 * ( tx + ty + tz );
                sum = 0.0;
                for ( ix = 0; ix <= l12; ix++ ) {
                    for ( iy = 0; iy <= m12; iy++ ) {
                        psum = dxp[ix] * dyp[iy];
                        rp = rsum[ix][iy];
                        for ( iz = 0; iz <= n12; iz++ ) {
                            sum += psum * dzp[iz] * rp[iz];
                        }
                    }
                }
                vvals[kc] += nfct * sum;
            }
        }
    }
}

void
one_electron_ints ( double *Smat, double *Tmat, double *Hmat )
{
    int ils, jls, knt, ir, jr, kc, ijr, ir2;
    int ishell, jshell, cn1, off1, nls1, cn2, off2, nls2, dsize, rsize;
    int maxlsh, asize;
    int *ostates;
    double *svals, *tvals, *vvals;
    const Shell *sh1, *sh2, *shell;
    double **normFactors;
    const int nshell = number_of_shells ();
    normFactors = normalizationFactors ();
    lxyz = lvector_components ();
    skipcen = skip_center ();
    shell = shells ();
    center = centers ();
    ncen = number_of_centers ();
    maxlsh = max_lvalue ();
    dsize = maxlsh + 2;
    rsize = 2 * maxlsh + 5;
    dx = dtensor ( dsize, dsize, rsize );
    dy = dtensor ( dsize, dsize, rsize );
    dz = dtensor ( dsize, dsize, rsize );
    rsize = 2 * maxlsh + 1;
    rsum = dtensor ( rsize, rsize, rsize );
    rterm = dtensor ( rsize, rsize, rsize );
    asize = max_lstates ();
    asize = asize * asize;
    ostates = ialloc ( asize );
    svals = dalloc ( asize );
    tvals = dalloc ( asize );
    vvals = dalloc ( asize );
    nrmstates = dalloc ( asize );
    lstates = ialloc ( asize );
    for ( ishell = 0; ishell < nshell; ishell++ ) {
        sh1 = ( shell + ishell );
        npr1 = sh1->npr;
        lv1 = sh1->lsh;
        cn1 = sh1->cen;
        off1 = sh1->off;
        nls1 = sh1->nst;
        alf1 = sh1->al;
        cof1 = sh1->co;
        a = ( center + cn1 )->r;
        for ( jshell = 0; jshell <= ishell; jshell++ ) {
            sh2 = ( shell + jshell );
            npr2 = sh2->npr;
            lv2 = sh2->lsh;
            cn2 = sh2->cen;
            off2 = sh2->off;
            nls2 = sh2->nst;
            alf2 = sh2->al;
            cof2 = sh2->co;
            b = ( center + cn2 )->r;
            absq = distSqr ( a, b );
            knt = 0;
            for ( ils = 0; ils < nls1; ils++ ) {
                ir = off1 + ils;
                ir2 = ( ir * ( ir + 1 ) / 2 );
                for ( jls = 0; jls < nls2; ++jls ) {
                    jr = off2 + jls;
                    if ( jr > ir )
                        break;
                    lstates[knt] = ( ils << 8 ) + jls;
                    ostates[knt] = ir2 + jr;
                    nrmstates[knt] = normFactors[lv1][ils] *
                                     normFactors[lv2][jls];
                    svals[knt] = 0.0;
                    tvals[knt] = 0.0;
                    vvals[knt] = 0.0;
                    ++knt;
                }
            }
            if ( !knt )
                continue;
            calc_poeints ( knt, svals, tvals, vvals );
            for ( kc = 0; kc < knt; kc++ ) {
                ijr = ostates[kc];
                Smat[ijr] = svals[kc];
                Tmat[ijr] = tvals[kc];
                Hmat[ijr] = tvals[kc] + vvals[kc];
            }
        }
    }
    free ( lstates );
    free ( ostates );
    free ( nrmstates );
    free ( vvals );
    free ( tvals );
    free ( svals );
    free_dtensor ( rterm, rsize, rsize );
    free_dtensor ( rsum, rsize, rsize );
    free_dtensor ( dz, dsize, dsize );
    free_dtensor ( dy, dsize, dsize );
    free_dtensor ( dx, dsize, dsize );
}

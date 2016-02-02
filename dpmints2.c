#include "dpmints.h"
#define SRTERM 34.9868366552497250

static double ***dx, ***dy, ***dz, ***rs;
static double absq;
static const double *alf1, *cof1, *alf2, *cof2;
static const double *a, *b;
static int lv1, lv2, npr1, npr2;
static int *lstates;
static int *ostates;
static double *nstates;
static int ***lxyz;
static double **normFactors;
static const double *q;

inline void
calcPDPM ( double *vints, int nlst )
{
    int iroot, kc, lpr, jpr, je, ipr;
    int ix, iy, iz;
    int lvt, order;
    int l1, m1, n1, l2, m2, n2, ls1, ls2;
    int l12, m12, n12;
    double abx, aby, abz, b00;
    double axp, c1, f12, bxp, c12, pxp, abi, s12, pa[3], p[3];
    double c34, qxp, cdi, sr, pq[3], pq2;
    double txp, w, t, rt, xt, yt, zt, nf, sum, nfc, fff;
    const int npalf = 6;
    const double palf[] = { 6.8505018000, 4.0491646000, 3.5941062989,
                            1.2478274000, 0.7927690989, 0.3377107978
                          };
    const double pcof[] = { 0.0766926784, 0.1483475912, 0.0462334641,
                            0.0717376427, 0.0447149792, 0.0069678529
                          };
    double pa[3],pb[3];
    const double *rz, *dzp;
    double cx, cy;
    abx = a[0] - b[0];
    aby = a[1] - b[1];
    abz = a[2] - b[2];
    lvt = lv1 + lv2;
    for ( ipr = 0; ipr < npr1; ++ipr ) {
        axp = alf1[ipr];
        c1 = cof1[ipr];
        je = npr2;
        f12 = 1.0;
        if ( alf1 == alf2 ) {
            je = ipr + 1;
            f12 = 2.0;
        }
        for ( jpr = 0; jpr < je; ++jpr ) {
            if ( jpr == ipr )
                f12 = 1.0;
            c12 = c1 * f12 * cof2[jpr];
            bxp = alf2[jpr];
            pxp = axp + bxp;
            abi = 1. / pxp;
            p[0] = ( axp * a[0] + bxp * b[0] ) * abi;
            p[1] = ( axp * a[1] + bxp * b[1] ) * abi;
            p[2] = ( axp * a[2] + bxp * b[2] ) * abi;
            pa[0] = p[0] - a[0];
            pa[1] = p[1] - a[1];
            pa[2] = p[2] - a[2];
            pb[0] = p[0] - b[0];
            pb[1] = p[1] - b[1];
            pb[2] = p[2] - b[2];
            s12 = exp ( -abi * axp * bxp * absq );
            abi = 0.5 * abi;
            formMDR ( dx, abi, pa[0], pb[0], lvt );
            formMDR ( dy, abi, pa[1], pb[1], lvt );
            formMDR ( dz, abi, pa[2], pb[2], lvt );
            for ( lpr = 0; lpr < npalf; ++lpr ) {
                c34 = pcof[lpr];
                qxp = palf[lpr];
                cdi = 1. / qxp;
                pq[0] = p[0] - q[0];
                pq[1] = p[1] - q[1];
                pq[2] = p[2] - q[2];
                pq2 = normSqr ( pq );
                txp = pxp + qxp;
                sr = c12 * c34 * SRTERM * s12 * abi * cdi / sqrt ( txp );
                w = pxp * qxp / txp;
                t = w * pq2;
                formMDR ( rs, sr, t, w, pq, lvt );
                for ( kc = 0; kc < nlst; ++kc ) {
                    nf = nstates[kc] * sr;
                    ls1 = lstates[kc];
                    ls2 = ls1 & 255;
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
                    nfc = nf * sr;
                    sum = 0.0;
                    for ( ix = 0; ix <= l12; ++ix ) {
                        cx = dx[l1][l2][ix];
                        for ( iy = 0; iy <= m12; ++iy ) {
                            cy = dy[m1][m2][iy] * cx;
                            dzp = dz[n1][n2];
                            rz = rs[ix][iy];
                            for ( iz = 0; iz <= n12; ++iz ) {
                                sum += dz[iz] * rz[iz] * cy;
                            }
                        }
                    }
                    vints[kc] += sum * nf;
                }
            }
        }
    }
}

void
dpm_ints ( double *hmat )
{
    int ils, jls, ish, jsh, ir, ir2, jr, ijr;
    int cen1, off1, cen2, off2, nls1, nls2, nlst, asize;
    double *vints;
    const Shell *sh1;
    const Shell *sh2;
    const Shell *shell = shells ();
    const int nshell = number_of_shells ();
    const Center *center = centers ();
    const int skipcen = skip_center ();
    q = ( center + skipcen )->r;
    normFactors = normalizationFactors ();
    lxyz = lvector_components ();
    asize = max_lstates ();
    asize = asize * asize;
    ostates = ialloc ( asize );
    lstates = ialloc ( asize );
    nstates = dalloc ( asize );
    vints = dalloc ( asize );
    int dsize = 2 * max_lvalue() + 1;
    dx = dtensor ( dsize, dsize, 2 * dsize + 1 );
    dy = dtensor ( dsize, dsize, 2 * dsize + 1 );
    dz = dtensor ( dsize, dsize, 2 * dsize + 1 );
    rs = dtensor ( dsize, dsize, dsize );
    RysInit ();
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
            nlst = 0;
            for ( ils = 0; ils < nls1; ++ils ) {
                ir = off1 + ils;
                ir2 = ( ir * ( ir + 1 ) ) / 2;
                for ( jls = 0; jls < nls2; ++jls ) {
                    jr = off2 + jls;
                    if ( jr > ir )
                        break;
                    ostates[nlst] = ir2 + jr;
                    lstates[nlst] = ( ils << 8 ) + jls;
                    nstates[nlst] = normFactors[lv1][ils] *
                                    normFactors[lv2][jls];
                    vints[nlst] = 0.;
                    ++nlst;
                }
            }
            calcPDPM ( vints, nlst );
            for ( ils = 0; ils < nlst; ++ils ) {
                ijr = ostates[ils];
                hmat[ijr] -= vints[ils];
            }
        }
    }
    free_dtensor ( rs, dsize, dsize );
    free_dtensor ( dz, dsize, dsize );
    free_dtensor ( dy, dsize, dsize );
    free_dtensor ( dx, dsize, dsize );
    free ( ostates );
    free ( lstates );
    free ( nstates );
    free ( vints );
}


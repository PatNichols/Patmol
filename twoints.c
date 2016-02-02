#ifndef _TWOELEC_CC_
#define _TWOELEC_CC_

#include "twoints.h"
#define L_SHIFT 4
#define L_MASK 0xF
#define L_SHIFT2 8
#define L_TSHIFT1 12
#define L_TSHIFT2 24
#define L_TSHIFT3 36

typedef u_int64_t l_state_t;

static double Gx[MAXROOTS][MAXROOTS][MAXROOTS],
       Gy[MAXROOTS][MAXROOTS][MAXROOTS], Gz[MAXROOTS][MAXROOTS][MAXROOTS];
static double absq, cdsq;
static int npr1, lv1, npr2, lv2, npr3, lv3, npr4, lv4;
static const double *al1, *co1, *al2, *co2,
       *al3, *co3, *al4, *co4, *a, *b, *c, *d;
static double *nrm_states;
static l_state_t *lstates;

inline static void
calc_pteints ( ERINT * sints, int nlst )
{
    register int iroot, kc;
    register double sum, nfact;
    l_state_t ils;
    int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, l12, m12, n12, l34, m34,
        n34;
    int i, j, k, l, jend, lend;
    double axp, c1, bxp, c12, cxp, c3, dxp, c34, s12, s34, f12, f34;
    double pxp, abi, qxp, cdi, sr, txp, t, w, rr, drt, fff;
    double p[3], q[3], pq[3], pa[3], qc[3], pq2, b00;
    const double SRterm = 34.9868366552497250;
    const double threshold = 1.e-15;
    const double abx = a[0] - b[0];
    const double aby = a[1] - b[1];
    const double abz = a[2] - b[2];
    const double cdx = c[0] - d[0];
    const double cdy = c[1] - d[1];
    const double cdz = c[2] - d[2];
    const int lvt12 = lv1 + lv2;
    const int lvt34 = lv3 + lv4;
    const int lvt = lvt12 + lvt34;
    const int nroots = ( lvt / 2 ) + 1;
    const double *roots, *wghts;
    if ( a == b ) {
        s12 = 1.;
        __builtin_memcpy ( p, a, sizeof ( double ) * 3 );
        __builtin_memset ( pa, 0, sizeof ( double ) * 3 );
    }
    if ( c == d ) {
        s34 = 1.;
        __builtin_memcpy ( q, c, sizeof ( double ) * 3 );
        __builtin_memset ( qc, 0, sizeof ( double ) * 3 );
    }
    for ( i = 0; i < npr1; ++i ) {
        axp = al1[i];
        c1 = co1[i];
        f12 = 1.0;
        jend = npr2;
        if ( al1 == al2 ) {
            f12 = 2.0;
            jend = i + 1;
        }
        for ( j = 0; j < jend; ++j ) {
            if ( i == j )
                f12 = 1.0;
            c12 = c1 * f12 * co2[j];
            bxp = al2[j];
            pxp = axp + bxp;
            abi = 1.0 / pxp;
            if ( a != b ) {
                s12 = exp ( -axp * bxp * absq * abi );
                p[0] = ( axp * a[0] + bxp * b[0] ) * abi;
                p[1] = ( axp * a[1] + bxp * b[1] ) * abi;
                p[2] = ( axp * a[2] + bxp * b[2] ) * abi;
                pa[0] = p[0] - a[0];
                pa[1] = p[1] - a[1];
                pa[2] = p[2] - a[2];
            }
            for ( k = 0; k < npr3; ++k ) {
                cxp = al3[k];
                c3 = co3[k];
                f34 = 1.0;
                lend = npr4;
                if ( al3 == al4 ) {
                    f34 = 2.0;
                    lend = k + 1;
                }
                for ( l = 0; l < lend; ++l ) {
                    if ( k == l )
                        f34 = 1.0;
                    c34 = c3 * f34 * co4[l];
                    dxp = al4[l];
                    qxp = cxp + dxp;
                    cdi = 1.0 / qxp;
                    if ( c != d ) {
                        s34 = exp ( -cxp * dxp * cdsq * cdi );
                        q[0] = ( cxp * c[0] + dxp * d[0] ) * cdi;
                        q[1] = ( cxp * c[1] + dxp * d[1] ) * cdi;
                        q[2] = ( cxp * c[2] + dxp * d[2] ) * cdi;
                        qc[0] = q[0] - c[0];
                        qc[1] = q[1] - c[1];
                        qc[2] = q[2] - c[2];
                    }
                    txp = pxp + qxp;
                    sr = SRterm * s12 * s34 * abi * cdi / sqrt ( txp );
                    if ( sr < threshold )
                        continue;
                    sr *= c12 * c34;
                    pq[0] = p[0] - q[0];
                    pq[1] = p[1] - q[1];
                    pq[2] = p[2] - q[2];
                    pq2 = pq[0] * pq[0] + pq[1] * pq[1] + pq[2] * pq[2];
                    t = pq2 * pxp * qxp / txp;
                    Rys_findRoots ( nroots, t );
                    roots = RysRoots ();
                    for ( iroot = 0; iroot < nroots; ++iroot ) {
                        rr = roots[iroot];
                        drt = rr / ( 1.0 + rr );
                        fff = drt / txp;
                        b00 = 0.5 * fff;
                        Rys_setB00 ( b00 );
                        Rys_setB1 ( ( ( 0.5 - b00 * qxp ) / pxp ) );
                        Rys_setB1p ( ( ( 0.5 - b00 * pxp ) / qxp ) );
                        Rys_setC ( ( pa[0] - qxp * pq[0] * fff ) );
                        Rys_setCp ( ( qc[0] + pxp * pq[0] * fff ) );
                        RysRecur ( Gx[iroot], lvt12, lvt34 );
                        Rys_setC ( ( pa[1] - qxp * pq[1] * fff ) );
                        Rys_setCp ( ( qc[1] + pxp * pq[1] * fff ) );
                        RysRecur ( Gy[iroot], lvt12, lvt34 );
                        Rys_setC ( ( pa[2] - qxp * pq[2] * fff ) );
                        Rys_setCp ( ( qc[2] + pxp * pq[2] * fff ) );
                        RysRecur ( Gz[iroot], lvt12, lvt34 );
                    }
                    wghts = RysWeights ();
                    for ( kc = 0; kc < nlst; ++kc ) {
                        ils = lstates[kc];
                        n4 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        m4 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        l4 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        n3 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        m3 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        l3 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        n2 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        m2 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        l2 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        n1 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        m1 = ( int ) ( ils & L_MASK );
                        ils >>= L_SHIFT;
                        l1 = ( int ) ( ils & L_MASK );
                        l12 = l1 + l2;
                        m12 = m1 + m2;
                        n12 = n1 + n2;
                        l34 = l3 + l4;
                        m34 = m3 + m4;
                        n34 = n3 + n4;
                        sum = 0.0;
                        for ( iroot = 0; iroot < nroots; ++iroot ) {
                            sum +=
                                RysShift ( abx, cdx, Gx[iroot], l12, l2, l34,
                                           l4 ) *
                                RysShift ( aby, cdy, Gy[iroot], m12, m2, m34,
                                           m4 ) *
                                RysShift ( abz, cdz, Gz[iroot], n12, n2, n34,
                                           n4 ) * wghts[iroot];
                        }
                        ( sints + kc )->val += sum * nrm_states[kc] * sr;
                    }
                }
            }
        }
    }
}

void
two_electron_ints ( ERIFS * tfiles, int start )
{
    int kc, lls, kls, jls, ils, ir, jr, kr, lr, l1, l2, l3, l4, knt;
    int cen1, off1, nls1, cen2, off2, nls2, cen3, off3, nls3, cen4, off4, nls4,
        it;
    int switch12, switch34, nfiles, asize, fsize;
    int ish, jsh, ksh, lsh, rank, nproc, maxfsize;
    double nf12, nf123, etime;
    int *numints, *lst, pknt;
    const double *dp;
    const double threshold = 1.e-14;
    const char *basename = tfiles->basename;
    ERINT *sints;
    const Shell *sh1, *sh2, *sh3, *sh4;
    const Shell *shell = shells ();
    const Center *center = centers ();
    const int nshell = number_of_shells ();
    FILE *out;
    stopwatch_t timer;
    int ***lxyz = lvector_components();
    double **normFactors = normalizationFactors();
    l_state_t lk1, lk2, lk3, lk4;
    maxfsize = MAXFILESIZE / sizeof ( ERINT );
    rank = tfiles->rank;
    nproc = tfiles->nproc;
    numints = tfiles->numints;
    asize = max_lstates ();
    asize = asize * asize;
    asize = asize * asize;
    nrm_states = dalloc ( asize );
    sints = ( ERINT * ) pmalloc ( sizeof ( ERINT ) * asize );
    lstates = ( l_state_t* ) pmalloc ( sizeof ( l_state_t ) * asize );
    RysInit ();
    nfiles = 0;
    fsize = 0;
    out = createERIFile ( rank, nfiles, basename );
    stopwatch_clear ( &timer );
    stopwatch_start ( &timer );
    pknt = 0;
    for ( ish = start; ish < nshell; ++ish ) {
        sh1 = ( shell + ish );
        npr1 = sh1->npr;
        l1 = lv1 = sh1->lsh;
        cen1 = sh1->cen;
        off1 = sh1->off;
        nls1 = sh1->nst;
        al1 = sh1->al;
        co1 = sh1->co;
        a = ( center + cen1 )->r;
        for ( jsh = 0; jsh <= ish; ++jsh ) {
            sh2 = ( shell + jsh );
            npr2 = sh2->npr;
            l2 = lv2 = sh2->lsh;
            cen2 = sh2->cen;
            off2 = sh2->off;
            nls2 = sh2->nst;
            al2 = sh2->al;
            co2 = sh2->co;
            b = ( center + cen2 )->r;
            absq = distSqr ( a, b );
            switch12 = l1 < l2;
            if ( switch12 ) {
                it = npr1;
                npr1 = npr2;
                npr2 = it;
                it = lv1;
                lv1 = lv2;
                lv2 = it;
                dp = al1;
                al1 = al2;
                al2 = dp;
                dp = co1;
                co1 = co2;
                co2 = dp;
                dp = a;
                a = b;
                b = dp;
            }
            for ( ksh = 0; ksh <= ish; ++ksh ) {
                sh3 = ( shell + ksh );
                npr3 = sh3->npr;
                l3 = lv3 = sh3->lsh;
                cen3 = sh3->cen;
                off3 = sh3->off;
                nls3 = sh3->nst;
                al3 = sh3->al;
                co3 = sh3->co;
                c = ( center + cen3 )->r;
                for ( lsh = 0; lsh <= ksh; ++lsh ) {
#ifdef USE_MPI
                	/* pass out the integrals like cards with no communication */
                    ++pknt;
                    if ( ( pknt % nproc ) != rank ) {
                        continue;
                    }
                    /* reset pknt so it does not get very large */
                    pknt = rank;
#endif
                    sh4 = ( shell + lsh );
                    npr4 = sh4->npr;
                    l4 = lv4 = sh4->lsh;
                    cen4 = sh4->cen;
                    off4 = sh4->off;
                    nls4 = sh4->nst;
                    al4 = sh4->al;
                    co4 = sh4->co;
                    d = ( center + cen4 )->r;
                    cdsq = distSqr ( c, d );
                    switch34 = l3 < l4;
                    if ( switch34 ) {
                        it = npr3;
                        npr3 = npr4;
                        npr4 = it;
                        it = lv3;
                        lv3 = lv4;
                        lv4 = it;
                        dp = al3;
                        al3 = al4;
                        al4 = dp;
                        dp = co3;
                        co3 = co4;
                        co4 = dp;
                        dp = c;
                        c = d;
                        d = dp;
                    }
                    knt = 0;
                    for ( ils = 0; ils < nls1; ++ils ) {
                        ir = off1 + ils;
                        for ( jls = 0; jls < nls2; ++jls ) {
                            jr = off2 + jls;
                            if ( jr > ir )
                                break;
                            lk1 = (lxyz[l1][ils][0]<< L_SHIFT2) | (lxyz[l1][ils][1]<<L_SHIFT) | lxyz[l1][ils][2] ; 
                            lk2 = (lxyz[l2][jls][0]<< L_SHIFT2) | (lxyz[l2][jls][1]<<L_SHIFT) | lxyz[l2][jls][2] ; 
                            nf12 = normFactors[l1][ils] * normFactors[l2][jls];
                            for ( kls = 0; kls < nls3; ++kls ) {
                                kr = off3 + kls;
                                if ( kr > ir )
                                    break;
                                lk3 = (lxyz[l3][kls][0]<<L_SHIFT2) | (lxyz[l3][kls][1]<<L_SHIFT) | lxyz[l3][kls][2] ; 
                                nf123 = nf12 * normFactors[l3][kls];
                                for ( lls = 0; lls < nls4; ++lls ) {
                                    lr = off4 + lls;
                                    if ( lr > kr )
                                        break;
                                    if ( ir == kr && lr > jr )
                                        break;
                                    lk4 = (lxyz[l4][lls][0]<<L_SHIFT2) | (lxyz[l4][lls][1]<<L_SHIFT) | lxyz[l4][lls][2] ; 
                                    ( sints + knt )->val = 0.0;
                                    ( sints + knt )->i = ( unsigned short ) ir;
                                    ( sints + knt )->j = ( unsigned short ) jr;
                                    ( sints + knt )->k = ( unsigned short ) kr;
                                    ( sints + knt )->l = ( unsigned short ) lr;
                                    if ( !switch12 ) {
                                        if ( !switch34 ) {
                                            lstates[knt] = (lk1 << L_TSHIFT3 ) | ( lk2 << L_TSHIFT2 ) | (lk3 << L_TSHIFT1 ) | lk4;
                                        } else {
                                            lstates[knt] =( lk1 << L_TSHIFT3 ) | ( lk2 <<  L_TSHIFT2 ) | ( lk4 << L_TSHIFT1 ) | lk3;
                                        }
                                    } else {
                                        if ( !switch34 ) {
                                            lstates[knt] = (lk2 << L_TSHIFT3) | (lk1 << L_TSHIFT2) | (lk3<< L_TSHIFT1) | lk4;
                                        } else {
                                            lstates[knt] = (lk2 << L_TSHIFT3) | (lk1 << L_TSHIFT2) | (lk4<< L_TSHIFT1) | lk3;
                                        }
                                    }
                                    nrm_states[knt] = nf123 *
                                                      normFactors[l4][lls];
                                    ++knt;
                                }
                            }
                        }
                    }
                    if ( !knt ) {
                        if ( switch34 ) {
                            npr3 = npr4;
                            lv3 = lv4;
                            al3 = al4;
                            co3 = co4;
                            c = d;
                        }
                        continue;
                    }
                    calc_pteints ( sints, knt );
                    for ( kc = 0; kc < knt; ++kc ) {
                        if ( fabs ( ( sints + kc )->val ) > threshold ) {
                            fwrite ( ( sints + kc ), sizeof ( ERINT ), 1, out );
                            ++fsize;
                        }
                    }
                    if ( fsize > maxfsize ) {
                        fclose ( out );
                        numints[nfiles] = fsize;
                        ++nfiles;
                        out = createERIFile ( rank, nfiles, basename );
                        fsize = 0;
                    }
                    if ( switch34 ) {
                        npr3 = npr4;
                        lv3 = lv4;
                        al3 = al4;
                        co3 = co4;
                        c = d;
                    }
                }
            }
            if ( switch12 ) {
                npr1 = npr2;
                lv1 = lv2;
                al1 = al2;
                co1 = co2;
                a = b;
            }
        }
    }
    stopwatch_stop ( &timer );
    etime = stopwatch_elapsed_time ( &timer );
    fprintf ( stderr, "Time for Two Electrons Integrals = %10.2le seconds\n",
              etime );
    if ( fsize ) {
        fclose ( out );
        numints[nfiles] = fsize;
        ++nfiles;
    }
    tfiles->nfiles = nfiles;
    free ( lstates );
    free ( sints );
    free ( nrm_states );
}
#endif


#include "shell.h"

void
read_shell ( FILE * infile, Shell * sh, int *off )
{
    int npr, lsh, cen, nls, ipr;
    fscanf ( infile, "%d %d %d", &npr, &lsh, &cen );
    sh->npr = npr;
    sh->lsh = lsh;
    sh->cen = cen;
    sh->off = ( *off );
    nls = ( ( lsh + 1 ) * ( lsh + 2 ) ) / 2;
    sh->nst = nls;
    ( *off ) += nls;
    sh->al = dalloc ( npr );
    sh->co = dalloc ( npr );
    for ( ipr = 0; ipr < npr; ++ipr ) {
        fscanf ( infile, "%lg %lg", ( sh->al ) + ipr, ( sh->co ) + ipr );
    }
    normalize_shell ( sh );
}



void
normalize_shell ( Shell * sh )
{
    int ipr, jpr, npr;
    double lpow, sum, c1, c2, a1, a2;
    double *alf, *cof;
// twofact = 2 * sqrt(2)
    const double twofact = 2.8284271247461903;
// PITERM =  (pi)^(1.5) = pi * sqrt(pi)
    const double PITERM = 5.568327996831707;
    npr = ( sh->npr );
    lpow = 1.5 + ( sh->lsh );
    alf = ( sh->al );
    cof = ( sh->co );
    sum = 0.0;
    for ( ipr = 0; ipr < npr; ++ipr ) {
        a1 = alf[ipr];
        c1 = cof[ipr];
        for ( jpr = 0; jpr < npr; ++jpr ) {
            a2 = alf[jpr];
            c2 = cof[jpr];
            sum += c1 * c2 * pow ( sqrt ( a1 * a2 ) / ( a1 + a2 ), lpow );
        }
    }
    sum *= twofact;
    sum = 1. / sqrt ( sum );
    for ( ipr = 0; ipr < npr; ++ipr ) {
        a1 = sum * sqrt ( pow ( 2. * alf[ipr], lpow ) / PITERM );
        cof[ipr] *= a1;
    }
}

void
write_shell ( FILE * outfile, Shell * sh )
{
    int ipr;
    fprintf ( outfile, "%3d %3d %3d %3d\n", sh->npr, sh->lsh, sh->cen, sh->off );
    for ( ipr = 0; ipr < sh->npr; ++ipr ) {
        fprintf ( outfile, "%20.8le %20.8le\n", ( sh->al ) [ipr], ( sh->co ) [ipr] );
    }
}

int
shell_eq ( const Shell * s1, const Shell * s2 )
{
    int np, i;
    const double *a1, *c1, *a2, *c2;
    if ( s1 == s2 )
        return 1;
    if ( s1->npr != s2->npr || s1->lsh != s2->lsh )
        return 0;
    np = s1->npr;
    a1 = s1->al;
    c1 = s1->co;
    a2 = s2->al;
    c2 = s2->co;
    for ( i = 0; i < np; ++i ) {
        if ( fabs ( a1[i] - a2[i] ) > DBL_EPSILON ||
                fabs ( c1[i] - c2[i] ) > DBL_EPSILON ) {
            return 0;
        }
    }
    return 1;
}

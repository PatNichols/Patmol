#include "md.h"

void
formMDD ( double ***d, double abi, double ax, double bx, int l1, int l2 )
{
    register int n, ipj;
    int i, j;
    double *dx;
    const double *dm1;
    const int ltot = l1 + l2;
    if ( ltot == 0 ) {
        d[0][0][0] = 1.0;
        return;
    }
    for ( i = 0; i <= l1; i++ ) {
        for ( j = 0; j <= l2; j++ ) {
            for ( n = 0; n <= ltot; n++ )
                d[i][j][n] = 0.0;
        }
    }
    d[0][0][0] = 1.0;
    for ( i = 1; i <= l2; i++ ) {
        dm1 = d[0][i - 1];
        dx = d[0][i];
        dx[0] = bx * dm1[0] + dm1[1];
        for ( n = 1; n < i; n++ ) {
            dx[n] = abi * dm1[n - 1] + bx * dm1[n] + ( n + 1 ) * dm1[n + 1];
        }
        dx[i] = abi * dm1[i - 1];
    }
    for ( j = 1; j <= l1; j++ ) {
        for ( i = 0; i <= l2; i++ ) {
            dm1 = d[j - 1][i];
            dx = d[j][i];
            ipj = i + j;
            dx[0] = ax * dm1[0] + dm1[1];
            for ( n = 1; n < ipj; n++ ) {
                dx[n] = abi * dm1[n - 1] + ax * dm1[n] + ( n + 1 ) * dm1[n + 1];
            }
            dx[ipj] = abi * dm1[ipj - 1];
        }
    }
}


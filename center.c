
#include "center.h"

void
read_center ( FILE * infile, Center * cen )
{
    fscanf ( infile, "%lg %lg %lg %lg\n",
             & ( cen->q ), cen->r, cen->r + 1, cen->r + 2 );
}

void
write_center ( FILE * outfile, const Center * cen )
{
    fprintf ( outfile, "%15.6le %15.6le %15.6le %15.6le\n",
              ( cen->q ), ( cen->r ) [0], ( cen->r ) [1], ( cen->r ) [2] );
}

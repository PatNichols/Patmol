#include "erifs.h"

ERIFS *
init_ERIFS ( int nst, int nend, const char *sname )
{
    size_t slen;
    int64_t sn, fn;
    ERIFS *e = pmalloc ( sizeof ( ERIFS ) );
#ifdef USE_MPI
    MPI_Comm_size ( MPI_COMM_WORLD, & ( e->nproc ) );
    MPI_Comm_rank ( MPI_COMM_WORLD, & ( e->rank ) );
#else
    e->nproc = 1;
    e->rank = 0;
#endif
    e->nfiles = 0;
    sn = ( nst * ( nst + 1 ) ) / 2;
    sn = ( sn * ( sn + 1 ) ) / 2;
    fn = ( nend * ( nend + 1 ) ) / 2;
    fn = ( fn * ( fn + 1 ) ) / 2;
    fn -= sn;
    fn *= sizeof ( ERINT );
    fn /= MAXFILESIZE;
    ++fn;
    e->numints = ialloc ( ( size_t ) fn );
    fprintf ( stderr, "MAX NUMBER OF INT FILES = %d\n", ( int ) fn );
    slen = strlen ( sname ) + 1;
    e->basename = ( char * ) pmalloc ( slen );
    strncpy ( e->basename, sname, slen );
    return e;
}

inline void
create_filename ( int rank, int file_no, const char *bname, char *fname )
{
    const char totext[] = "0123456789";
    char buffer[6];
    if ( rank > 99999 ) {
        fatalError ( "rank is too large for ERIFS" );
    }
    if ( file_no > 99999 ) {
        fatalError ( "file number is too large for ERIFS" );
    }
    strncpy ( fname, TMPDIR_NAME, TMPDIR_LEN );
    sprintf(buffer,"%d",rank);
    strncat(fname, buffer, 6);
    sprintf(buffer,"%d",file_no);
    strncat(fname,buffer,6);
    strncat ( fname, bname, strlen ( bname ) + 1 );
}

FILE *
openERIFile ( int rank, int file_no, const char *bname )
{
    char fname[128];
    create_filename ( rank, file_no, bname, fname );
    return openFile ( fname, "r" );
}

FILE *
createERIFile ( int rank, int file_no, const char *bname )
{
    char fname[128];
    create_filename ( rank, file_no, bname, fname );
    return openFile ( fname, "w" );
}

void
clean_ERIFS ( ERIFS * e )
{
    int i;
    FILE *fp;
    char fname[128];
    for ( i = 0; i < ( e->nfiles ); ++i ) {
        create_filename ( e->rank, i, e->basename, fname );
        remove ( fname );
    }
    free ( e->basename );
    free ( e->numints );
    free ( e );
}



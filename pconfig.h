#ifndef _PCONFIG_H_
#define _PCONFIG_H_
#include "config.h"

#define TMPDIR_NAME "/glade/scratch/pnichols/"
#define TMPDIR_LEN 25
#define MAXFILESIZE 1048576*1024
//#define BINSIZE 1000
#define MAXL 5
#define MAXORB (1<<15)

#ifndef inline
#ifndef __inline__
#define inline
#else
#define inline __inline__
#endif
#endif

#endif

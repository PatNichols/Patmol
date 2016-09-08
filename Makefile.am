
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = patmpi patmol

patmpi_SOURCES = shell.c center.c basis.c dpmints.c erifs.c mdd.c mdr.c momints.c oneints.c twoints.c util.c spack_math.c \
rys.c scfaux.c patmpi.c rhfscf_mpi.c uhfscf_mpi.c 
 
patmpi_CFLAGS  =  -DUSE_MPI -O2 -g -funroll-all-loops -ffast-math -mtune=native -mavx
 
patmol_SOURCES = shell.c center.c basis.c dpmints.c erifs.c mdd.c mdr.c  momints.c oneints.c twoints.c util.c spack_math.c \
rys.c scfaux.c patmol.c rhfscf.c uhfscf.c 

patmol_CFLAGS  =  -O2 -funroll-all-loops -ffast-math -mtune=native -mavx
BUILT_SOURCES  = basis.c center.c dpmints.c erifs.c mdd.c mdr.c momints.c oneints.c patmol.c patmpi.c \
rhfscf.c rhfscf_mpi.c rys.c scfaux.c shell.c spack_math.c twoints.c uhfscf.c uhfscf_mpi.c util.c



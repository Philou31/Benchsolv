PROJECT=/scratch/algo/leleux/Benchsolv/
#####################################
## LIBRARIES
#####################################
LIB_HOME=/scratch/algo/leleux/lib
#############
# MPI
#############
OMPIDIR=$(LIB_HOME)/openmpi-2.0.1/build
IOMPI2= -I$(OMPIDIR)/include
LOMPI2= -L$(OMPIDIR)/lib -lmpi -lmpi_mpifh
#############
# QR_MUMPS
#############
QRM_HOME=$(LIB_HOME)/qr_mumps-2.0/build
IQRMUMPS2= -I$(QRM_HOME)/include -Dhave_colamd -Dhave_metis -Dhave_scotch -Dhave_starpu -Ddprec
LQRMUMPS2= -L${QRM_HOME}/lib -ldqrm -lsqrm -lcqrm -lzqrm -lgomp -lqrm_common
#############
#STARPU
#############
STARPU_HOME=$(LIB_HOME)/starpu-1.2.0/build
ISTARPU2= -I$(STARPU_HOME)/include/starpu/1.2
LSTARPU2= -L$(STARPU_HOME)/lib -lstarpu-1.2 -lstarpumpi-1.2
#############
# MUMPS
#############
MUMPS_HOME=$(LIB_HOME)/MUMPS_5.0.2
IMUMPS2= -I$(MUMPS_HOME)/include
LMUMPS2= -L$(MUMPS_HOME)/lib -lcmumps -lsmumps -ldmumps -lzmumps -lmumps_common -lpord
#############
# METIS
#############
METIS_HOME=$(LIB_HOME)/metis-5.1.0/build
IMETIS2= -I$(METIS_HOME)/include
LMETIS2= -L$(METIS_HOME)/lib -lmetis
PARMETIS_HOME=$(LIB_HOME)/parmetis-4.0.3/build
IPARMETIS2= -I$(PARMETIS_HOME)/include
LPARMETIS2= -L$(PARMETIS_HOME)/lib -lparmetis
#############
# SCOTCH
#############
SCOTCH_HOME=${LIB_HOME}/scotch_6.0.4
ISCOTCH2= -I$(SCOTCH_HOME)/include
LSCOTCH2= -L$(SCOTCH_HOME)/lib -lesmumps -lscotch -lscotcherr -lscotcherrexit -lscotchmetis
LPTSCOTCH2= -lptesmumps -lptscotch -lptscotcherr -lptscotcherrexit -lptscotchparmetis
#############
# SCALAPACK
#############
SCALAPACK_HOME=$(LIB_HOME)/scalapack-2.0.2/
LSCALAPACK2= -L$(SCALAPACK_HOME) -lscalapack
#############
# BLAS
#############
BLAS_HOME=$(LIB_HOME)/OpenBLAS-0.2.19
LBLAS2= -L$(BLAS_HOME) -lopenblas
#############
# SUITESPARSE
#############
SUITESPARSE_HOME=$(LIB_HOME)/SuiteSparse
ISUITESPARSE2= -I$(SUITESPARSE_HOME)/include
LSUITESPARSE2= -L$(SUITESPARSE_HOME)/lib -lsuitesparseconfig
LCOLAMD2= -L$(SUITESPARSE_HOME)/lib -lcolamd

#####################################
## COMPILER
#####################################
CC=$(OMPIDIR)/bin/mpic++
CFLAGS=-c -O2 -Wall -std=c++11 $(ISTARPU2) $(IMETIS2) $(ISCOTCH2) $(IMUMPS2) $(IOMPI2) $(IPARMETIS2) $(ICOLAMD2) $(ISUITESPARSE2) $(IQRMUMPS2) -MMD -MP -fopenmp
LDFLAGS=$(LQRMUMPS2) $(LSTARPU2) $(LMUMPS2) $(LSCOTCH2) $(LPTSCOTCH2) $(LPARMETIS2) $(LMETIS2) $(LCOLAMD2) $(LSUITESPARSE2) $(LSCALAPACK2) $(LBLAS2) $(LOMPI2) -std=c++11 -fopenmp -lgfortran

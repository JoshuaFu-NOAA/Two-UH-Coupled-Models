export

SHELL = /bin/sh
 
ARCH  = LINUX

srcdir = .
top_srcdir = .

prefix = .
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include
oldincludedir = /usr/include
infodir = ${prefix}/info
mandir = ${prefix}/man

sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var

program_transform_name = s,x,x,

MPIROOT        = /usr/mpi/pgi/mvapich2-1.7
MPI_LIB        = -L$(MPIROOT)/lib 
MPI_INCLUDE    = -I$(MPIROOT)/include

NETCDFROOT     = /home/xfu/bin/netcdf-413
NETCDF_LIB     = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff
NETCDF_INCLUDE = -I$(NETCDFROOT)/include

LIB      = -L../lib -lsupport -llapack -lblas
LIBS     = $(LIB) $(NETCDF_LIB) $(MPI_LIB)

MODOPT   = -I
MODULES  = ../modules

INCLUDE  = ../include
INCLUDES = $(MODOPT)$(MODULES) -I$(INCLUDE) $(NETCDF_INCLUDE) $(MPI_INCLUDE)

F90      = /usr/mpi/pgi/mvapich2-1.7/bin/mpif90
FC       = /usr/mpi/pgi/mvapich2-1.7/bin/mpif90
CC       = gcc
CPP      = gcc -E
AR       = ar
AS       = as 

DEFS     = -DHAVE_CONFIG_H

CFLAGS   = -I../config -O -DpgiFortran
FFLAGS   = -O
F90FLAGS = $(INCLUDES) -O -DpgiFortran -r8 -Mmpi -Mpreprocess -Mbyteswapio 
CPPFLAGS = 
ARFLAGS  = crv
LDFLAGS  = -O -DpgiFortran -r8 -Mmpi -Mpreprocess -Mbyteswapio

SRCDIRS = lapack blas support modules src

all:
	@for DIR in $(SRCDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done 

clean:
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f config.cache
	-rm -f lib/*.a bin/echam5
	-rm -f html/[a-z]*

tar:
	@tarfile=../echam5.f90.`date +%y%m%d`.taz ; gtar zcvf $$tarfile \
	`ls */*.f90 */*.[fhc] */*inc */Makefile Makefile.in Makefile run/hjob*`

index:
	-rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html support modules src include

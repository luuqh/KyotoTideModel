#!/bin/sh
#
#  $Id: Makefile 13 2008-12-12 10:53:49Z ishikawa $
#
SDIR		= .
XDIR		= ..
#EXEC		= nwnpac2007.bin
#EXEC		= jp68m.bin
EXEC		= rokka2007.bin
#EXEC		= a.out
#
#OMP_NUM_THREADS=8
OMP_NUM_THREADS=4
#OMP_NUM_THREADS=2


export OMP_NUM_THREADS
#
# Fortran compiler
F90		= mpifrt
FFLAGS		= -Qt,p,i,a,d,x -Kparallel,instance=$(OMP_NUM_THREADS) -X9 -Am

CPPFLAGS	= -Cpp -DNESTED -DROKKA -DDEBUG -DBISMFRIC -DNCIO -DNCEPSF  -DEVAPO -DSGOUT -DXTIDE -DEXP2007 -DMASKLVIS \
			-DBNDTIDE -DT3OUT -DEXPTEST -DEXPTIDE

OFLAGS		= -Kparallel,array_subscript,hardbarrier

MODULES = param.o date_manage.o mod_time.o mod_mpi.o
OGCM = main.o init.o
FWD = set.o tracli.o brtrop.o write.o mmchk.o \
  calpp.o ioutil.o ref_month.o surf_bulk.o \
  noh_mix.o gm_dif.o ice.o my_mpi.o tidelib.o
NEST = nest.o
NCDF = ncutils.o
TIDE = mod_nao.o mod_dump.o
OBJS = $(MODULES) $(OGCM) $(FWD) $(NEST) $(NCDF) calender.o $(TIDE)

LIB = -lnetcdf
NCDIR = /thin/local/NETCDF/netcdf-3.6.1
LIBDIR = -L$(NCDIR)/lib
INCDIR = -I$(NCDIR)/include

# Object name
FSRC   = $(OBJS:.o=.F)

all : $(XDIR)/$(EXEC)

FSRC   = $(OBJS:.o=.F90)
LIST = ${addsuffix .lst,$(FSRC)}
ERROR = ${addsuffix .err,$(FSRC)}
$(XDIR)/$(EXEC) :	$(OBJS)
	$(F90) $(OFLAGS) $(OBJS) -o $@ $(LIBDIR) $(LIB)

main.o		: param.o mod_mpi.o common.h mod_time.o
my_mpi.o	: param.o mod_mpi.o
ioutil.o	: param.o mod_mpi.o
set.o		: param.o mod_mpi.o common.h mod_time.o tidelib.o mod_nao.o mod_dump.o
init.o		: param.o mod_mpi.o common.h mod_time.o
ref_month.o	: param.o common.h
surf_bulk.o	: param.o mod_mpi.o common.h
tracli.o	: param.o mod_mpi.o common.h
noh_mix.o	: param.o mod_mpi.o common.h
gm_dif.o	: param.o mod_mpi.o common.h
ice.o		: param.o mod_mpi.o common.h
brtrop.o	: param.o mod_mpi.o common.h
write.o		: param.o mod_mpi.o common.h
mmchk.o		: param.o mod_mpi.o common.h
calpp.o		: param.o common.h
nest.o		: param.o mod_mpi.o common.h mod_time.o
ncutils.o	: param.o mod_mpi.o
calender.o	: param.o mod_mpi.o common.h
param.o		:
mod_mpi.o	: param.o
mod_time.o	: date_manage.o
tidelib.o   : param.o mod_mpi.o mod_nao.o common.h
date_manage.o : mod_mpi.o
		$(F90) -c $(CPPFLAGS) $(FFLAGS) $(INCDIR) -Qodate_manage.F90.lst -Zdate_manage.F90.err -o $@ date_manage.F90
mod_nao.o	: param.o mod_mpi.o
mod_dump.o	: param.o mod_mpi.o

# Fortran compilation
.SUFFIXES	:
.SUFFIXES	: .o .F90
.F90.o:
#%.o : %.F
	$(F90) -c $(CPPFLAGS) $(FFLAGS) $(INCDIR) -Fixed -Qo$<.lst -Z$<.err -o $@ $<
#	$(F90) -c $(CPPFLAGS) $(FFLAGS) -Qo$<.lst -Z$<.err -o $@ $<
#####	$(F77) -c $(CPPFLAGS) $(FFLAGS) -o $@ $<
clean:
	rm -f $(OBJS)
cleanall:
	rm -f $(OBJS) $(LIST) $(ERROR) *.mod
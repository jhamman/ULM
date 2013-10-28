# Makefile for UW implementation of Unified Noah SAC models
# Author: Ben Livneh, blivneh@hydro.washington.edu
# RCS Id string: $Id: makefile,v 1.2 2007/10/04 20:43:51 vicadmin Exp $

OBJS =  obj/driverMod.o \
	obj/MAIN_DRIVER.o \
	obj/READCNTL.o \
	obj/DATE_TIME.o \
        obj/GET_GRID.o \
        obj/GET_GRID_ASC.o \
        obj/GET_GRID_NETCDF.o \
	obj/ALLOCATE_ARRAYS.o \
	obj/READ_LSC.o \
        obj/READ_INT_ASC.o \
        obj/READ_REAL_ASC.o \
        obj/READ_INT_SPATIAL.o \
        obj/READ_REAL_SPATIAL.o \
        obj/READ_LSC_NETCDF.o \
        obj/READ_SNOWBANDS.o \
	obj/READ_INITIAL.o \
	obj/CHECK_INITIAL.o \
	obj/OPEN_FORCING.o \
	obj/OPEN_RESTART.o \
	obj/OPEN_OUTPUT.o \
	obj/INTERP_MONTHLY.o \
	obj/READ_FORCING.o \
	obj/READ_IN_VAR.o \
	obj/COMPUTE_NUM_FSTEPS.o \
	obj/CALC_WEIGHTS.o \
	obj/INTERP_FORCING.o \
	obj/LOCALTIME.o \
	obj/COSZENITH.o \
	obj/FINTERP.o \
	obj/ZTERP.o \
	obj/QDATAP.o \
	obj/E.o \
	obj/DQSDT.o \
	obj/SFLXALL_SRC.o \
	obj/WRITE_OUTPUT.o \
	obj/WRITE_PE.o \
	obj/WRITE_OUT_VAR.o \
	obj/WRITE_RESTART.o \
	obj/CLOSE_OUTPUT.o \
	obj/CLOSE_FORCING.o \
	obj/fland2.o \
	obj/fst2sac1.o \
	obj/fst2sac2.o \
	obj/sac2frz1.o \
	obj/frz2sac1.o \
	obj/smc2sac1.o \
	obj/h2o2sac1.o \
	obj/ex_sac1.o \
	obj/sac1.o \
	obj/begin.o \
	obj/chek.o \
	obj/datain.o \
	obj/guess.o \
	obj/meta.o \
	obj/nwsnow.o \
	obj/obtain.o \
	obj/pemb.o \
	obj/snowot.o \
	obj/snowtw.o \
	obj/subs.o \
	obj/surfac.o \
	obj/vapor.o \
	obj/water.o \
	obj/windf.o 

# Tunable parameters 
#
# COMPILER      Compiler name
# FFLAGS	Flags to the fortran compiler
# LIBS		A list of libraries to use 
# PROGRAM	Name of the executable

#COMPILER = ifort -r8 -i4 -fp_port -align
COMPILER = pgf90 -Mdalign -Mextend
FFLAGS = -c -g
## LIBS on tsunami
#LIBS = -L$(LIB_NETCDF) -lnetcdf
# LIBS on flood -- extra "f"
LIBS = -L$(LIB_NETCDF) -lnetcdff

HDRS = -I$(INC_NETCDF)
#PROGRAM = ulm
#PROGRAM = ulm.flood.full
PROGRAM = ulm.swm
#PROGRAM = ulm.single
#PROGRAM = ulmflood
#PROGRAM = ulmfloodmaj
#PROGRAM = ulmfloodmajcolo
#PROGRAM = ulmfloodminorfull
#PROGRAM = ulm.swm

##tsunami: no additional flags needed
LIBS2 = 
#flood: two more flags needed after program name
#LIBS2 = -lnetcdf -lcurl
#PROGRAM = ulm_test
#PROGRAM = ulm_wcrit
#PROGRAM = ulm_psnow

$(PROGRAM): $(OBJS)
	$(COMPILER) $(OBJS) $(LIBS) $(HDRS) -o $(PROGRAM) $(LIBS2)

obj/driverMod.o: driver/driverMod.f90
	$(COMPILER) $(FFLAGS) driver/driverMod.f90 -o obj/driverMod.o

obj/MAIN_DRIVER.o: driver/MAIN_DRIVER.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/MAIN_DRIVER.f90 -o obj/MAIN_DRIVER.o

obj/READCNTL.o: driver/READCNTL.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READCNTL.f90 -o obj/READCNTL.o

obj/DATE_TIME.o: driver/DATE_TIME.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/DATE_TIME.f90 -o obj/DATE_TIME.o

obj/GET_GRID.o: driver/GET_GRID.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/GET_GRID.f90 -o obj/GET_GRID.o

obj/GET_GRID_ASC.o: driver/GET_GRID_ASC.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/GET_GRID_ASC.f90 -o obj/GET_GRID_ASC.o

obj/GET_GRID_NETCDF.o: driver/GET_GRID_NETCDF.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/GET_GRID_NETCDF.f90 -o obj/GET_GRID_NETCDF.o

obj/ALLOCATE_ARRAYS.o: driver/ALLOCATE_ARRAYS.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/ALLOCATE_ARRAYS.f90 -o obj/ALLOCATE_ARRAYS.o

obj/READ_LSC.o: driver/READ_LSC.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_LSC.f90 -o obj/READ_LSC.o

obj/READ_INT_ASC.o: driver/READ_INT_ASC.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_INT_ASC.f90 -o obj/READ_INT_ASC.o

obj/READ_REAL_ASC.o: driver/READ_REAL_ASC.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_REAL_ASC.f90 -o obj/READ_REAL_ASC.o

obj/READ_INT_SPATIAL.o: driver/READ_INT_SPATIAL.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_INT_SPATIAL.f90 -o obj/READ_INT_SPATIAL.o

obj/READ_REAL_SPATIAL.o: driver/READ_REAL_SPATIAL.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_REAL_SPATIAL.f90 -o obj/READ_REAL_SPATIAL.o

obj/READ_LSC_NETCDF.o: driver/READ_LSC_NETCDF.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_LSC_NETCDF.f90 -o obj/READ_LSC_NETCDF.o

obj/READ_SNOWBANDS.o: driver/READ_SNOWBANDS.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_SNOWBANDS.f90 -o obj/READ_SNOWBANDS.o

obj/READ_INITIAL.o: driver/READ_INITIAL.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_INITIAL.f90 -o obj/READ_INITIAL.o

obj/CHECK_INITIAL.o: driver/CHECK_INITIAL.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/CHECK_INITIAL.f90 -o obj/CHECK_INITIAL.o

obj/OPEN_FORCING.o: driver/OPEN_FORCING.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/OPEN_FORCING.f90 -o obj/OPEN_FORCING.o

obj/OPEN_RESTART.o: driver/OPEN_RESTART.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/OPEN_RESTART.f90 -o obj/OPEN_RESTART.o

obj/OPEN_OUTPUT.o: driver/OPEN_OUTPUT.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/OPEN_OUTPUT.f90 -o obj/OPEN_OUTPUT.o

obj/INTERP_MONTHLY.o: driver/INTERP_MONTHLY.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/INTERP_MONTHLY.f90 -o obj/INTERP_MONTHLY.o

obj/READ_FORCING.o: driver/READ_FORCING.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_FORCING.f90 -o obj/READ_FORCING.o

obj/READ_IN_VAR.o: driver/READ_IN_VAR.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/READ_IN_VAR.f90 -o obj/READ_IN_VAR.o

obj/COMPUTE_NUM_FSTEPS.o: driver/COMPUTE_NUM_FSTEPS.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/COMPUTE_NUM_FSTEPS.f90 -o obj/COMPUTE_NUM_FSTEPS.o

obj/CALC_WEIGHTS.o: driver/CALC_WEIGHTS.f
	$(COMPILER) $(FFLAGS) driver/CALC_WEIGHTS.f -o obj/CALC_WEIGHTS.o

obj/INTERP_FORCING.o: driver/INTERP_FORCING.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/INTERP_FORCING.f90 -o obj/INTERP_FORCING.o

obj/LOCALTIME.o: driver/LOCALTIME.f
	$(COMPILER) $(FFLAGS) driver/LOCALTIME.f -o obj/LOCALTIME.o

obj/COSZENITH.o: driver/COSZENITH.f
	$(COMPILER) $(FFLAGS) driver/COSZENITH.f -o obj/COSZENITH.o

obj/FINTERP.o: driver/FINTERP.f90
	$(COMPILER) $(FFLAGS) driver/FINTERP.f90 -o obj/FINTERP.o

obj/ZTERP.o: driver/ZTERP.f
	$(COMPILER) $(FFLAGS) driver/ZTERP.f -o obj/ZTERP.o

obj/QDATAP.o: driver/QDATAP.f
	$(COMPILER) $(FFLAGS) driver/QDATAP.f -o obj/QDATAP.o

obj/E.o: driver/E.f
	$(COMPILER) $(FFLAGS) driver/E.f -o obj/E.o

obj/DQSDT.o: driver/DQSDT.f
	$(COMPILER) $(FFLAGS) driver/DQSDT.f -o obj/DQSDT.o

obj/SFLXALL_SRC.o: physics/SFLXALL_SRC.f
	$(COMPILER) $(FFLAGS) physics/SFLXALL_SRC.f -o obj/SFLXALL_SRC.o

obj/WRITE_OUTPUT.o: driver/WRITE_OUTPUT.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/WRITE_OUTPUT.f90 -o obj/WRITE_OUTPUT.o

obj/WRITE_PE.o: driver/WRITE_PE.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/WRITE_PE.f90 -o obj/WRITE_PE.o

obj/WRITE_OUT_VAR.o: driver/WRITE_OUT_VAR.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/WRITE_OUT_VAR.f90 -o obj/WRITE_OUT_VAR.o

obj/WRITE_RESTART.o: driver/WRITE_RESTART.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/WRITE_RESTART.f90 -o obj/WRITE_RESTART.o

obj/CLOSE_OUTPUT.o: driver/CLOSE_OUTPUT.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/CLOSE_OUTPUT.f90 -o obj/CLOSE_OUTPUT.o

obj/CLOSE_FORCING.o: driver/CLOSE_FORCING.f90 obj/driverMod.o
	$(COMPILER) $(FFLAGS) driver/CLOSE_FORCING.f90 -o obj/CLOSE_FORCING.o

obj/ex_sac1.o: physics/ex_sac1.f
	$(COMPILER) $(FFLAGS) physics/ex_sac1.f -o obj/ex_sac1.o

obj/sac1.o: physics/sac1.f
	$(COMPILER) $(FFLAGS) physics/sac1.f -o obj/sac1.o

obj/fland2.o: physics/fland2.f
	$(COMPILER) $(FFLAGS) physics/fland2.f -o obj/fland2.o

obj/fst2sac1.o: physics/fst2sac1.f
	$(COMPILER) $(FFLAGS) physics/fst2sac1.f -o obj/fst2sac1.o

obj/fst2sac2.o: physics/fst2sac2.f
	$(COMPILER) $(FFLAGS) physics/fst2sac2.f -o obj/fst2sac2.o

obj/frz2sac1.o: physics/frz2sac1.f
	$(COMPILER) $(FFLAGS) physics/frz2sac1.f -o obj/frz2sac1.o

obj/smc2sac1.o: physics/smc2sac1.f
	$(COMPILER) $(FFLAGS) physics/smc2sac1.f -o obj/smc2sac1.o

obj/h2o2sac1.o: physics/h2o2sac1.f
	$(COMPILER) $(FFLAGS) physics/h2o2sac1.f -o obj/h2o2sac1.o

obj/sac2frz1.o: physics/sac2frz1.f
	$(COMPILER) $(FFLAGS) physics/sac2frz1.f -o obj/sac2frz1.o 

obj/begin.o: physics/begin.f
	$(COMPILER) $(FFLAGS) physics/begin.f -o obj/begin.o

obj/chek.o: physics/chek.f
	$(COMPILER) $(FFLAGS) physics/chek.f -o obj/chek.o

obj/datain.o: physics/datain.f
	$(COMPILER) $(FFLAGS) physics/datain.f -o obj/datain.o

obj/guess.o: physics/guess.f
	$(COMPILER) $(FFLAGS) physics/guess.f -o obj/guess.o

obj/meta.o: physics/meta.f
	$(COMPILER) $(FFLAGS) physics/meta.f -o obj/meta.o

obj/nwsnow.o: physics/nwsnow.f
	$(COMPILER) $(FFLAGS) physics/nwsnow.f -o obj/nwsnow.o

obj/obtain.o: physics/obtain.f
	$(COMPILER) $(FFLAGS) physics/obtain.f -o obj/obtain.o

obj/pemb.o: physics/pemb.f
	$(COMPILER) $(FFLAGS) physics/pemb.f -o obj/pemb.o

obj/snowot.o: physics/snowot.f
	$(COMPILER) $(FFLAGS) physics/snowot.f -o obj/snowot.o

obj/snowtw.o: physics/snowtw.f
	$(COMPILER) $(FFLAGS) physics/snowtw.f -o obj/snowtw.o

obj/subs.o: physics/subs.f
	$(COMPILER) $(FFLAGS) physics/subs.f -o obj/subs.o

obj/surfac.o: physics/surfac.f
	$(COMPILER) $(FFLAGS) physics/surfac.f -o obj/surfac.o

obj/vapor.o: physics/vapor.f
	$(COMPILER) $(FFLAGS) physics/vapor.f -o obj/vapor.o

obj/water.o: physics/water.f
	$(COMPILER) $(FFLAGS) physics/water.f -o obj/water.o

obj/windf.o: physics/windf.f
	$(COMPILER) $(FFLAGS) physics/windf.f -o obj/windf.o

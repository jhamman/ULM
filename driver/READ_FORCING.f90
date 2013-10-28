SUBROUTINE READ_FORCING(TSTEP,mySWdown,myLWdown,myTair,myQair,myRainf,mySnowf,myPSurf,myWind)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! This subroutine reads one time step's worth of the atm forcing data,
  ! for the entire region, for the given timestep (TSTEP).
  ! It assumes the file is already open, and takes the FORCING_NCID of the open file as an input parameter.

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_FORCING.f90,v 1.3 2007/10/24 23:50:10 vicadmin Exp $"/

  ! Define local variables
  INTEGER TSTEP
  REAL    mySWdown(landlen)
  REAL    myLWdown(landlen)
  REAL    myTair(landlen)
  REAL    myQair(landlen)
  REAL    myRainf(landlen)
  REAL    myLSRainf(landlen)
  REAL    myCRainf(landlen)
  REAL    mySnowf(landlen)
  REAL    myPSurf(landlen)
  REAL    myWind(landlen)
  REAL    myWind_E(landlen)
  REAL    myWind_N(landlen)
  INTEGER varid
  INTEGER land_idx
  LOGICAL LBAND

  LBAND = .FALSE.

!write(*,*)'LWIND2D',LWIND2D,'LLSRAINF',LLSRAINF
!write(*,*)'COMP_FORCING',COMP_FORCING

  ! Get data
  status = NF_INQ_VARID(FORCING_NCID,'SWdown',varid)
  CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,mySWdown,2,COMP_FORCING,LBAND)

  status = NF_INQ_VARID(FORCING_NCID,'LWdown',varid)
  CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myLWdown,2,COMP_FORCING,LBAND)

  status = NF_INQ_VARID(FORCING_NCID,'Tair',varid)
  CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myTair,2,COMP_FORCING,LBAND)

  status = NF_INQ_VARID(FORCING_NCID,'Qair',varid)
  CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myQair,2,COMP_FORCING,LBAND)

  IF (LLSRAINF) THEN

    status = NF_INQ_VARID(FORCING_NCID,'LSRainf',varid)
    CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myLSRainf,2,COMP_FORCING,LBAND)

    status = NF_INQ_VARID(FORCING_NCID,'CRainf',varid)
    CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myCRainf,2,COMP_FORCING,LBAND)
    myRainf = myLSRainf + myCRainf

  ELSE

    status = NF_INQ_VARID(FORCING_NCID,'Rainf',varid)
    CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myRainf,2,COMP_FORCING,LBAND)

  ENDIF

  status = NF_INQ_VARID(FORCING_NCID,'Snowf',varid)
  CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,mySnowf,2,COMP_FORCING,LBAND)

  status = NF_INQ_VARID(FORCING_NCID,'PSurf',varid)
  CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myPSurf,2,COMP_FORCING,LBAND)

  IF (LWIND2D) THEN

    status = NF_INQ_VARID(FORCING_NCID,'Wind_E',varid)
    CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myWind_E,2,COMP_FORCING,LBAND)

    status = NF_INQ_VARID(FORCING_NCID,'Wind_N',varid)
    CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myWind_N,2,COMP_FORCING,LBAND)
    DO land_idx = 1,landlen
      myWind(land_idx) = SQRT( myWind_E(land_idx)*myWind_E(land_idx) + myWind_N(land_idx)*myWind_N(land_idx) )
    END DO

  ELSE

    status = NF_INQ_VARID(FORCING_NCID,'Wind',varid)
    CALL READ_IN_VAR(FORCING_NCID,varid,TSTEP,myWind,2,COMP_FORCING,LBAND)

  ENDIF

END

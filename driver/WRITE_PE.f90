SUBROUTINE WRITE_PE(FORCING_STEP)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! WRITES VARS FOR 1 TIMESTEP TO THE PE FILE

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: WRITE_PE.f90,v 1.2 2007/10/27 22:01:28 vicadmin Exp $"/

  ! Define local variables
  INTEGER FORCING_STEP
  INTEGER varid
  INTEGER K,L
  INTEGER start,count
  REAL TIME_tmp
  LOGICAL LBAND

  LBAND = .TRUE.
  start = FORCING_STEP
  count = 1

  ! We only care about the PE file here

  K = 7
!! Hack to avoid writing the file

  if (K ==2 ) then

  ! Write Time,Timestep
  TIME_tmp = (FORCING_STEP-1)*FORCING_DT
  status = NF_INQ_VARID(OUT_NCIDS(K),'time',varid)
  status = NF_PUT_VARA_REAL(OUT_NCIDS(K),varid,start,count,TIME_tmp)
  status = NF_INQ_VARID(OUT_NCIDS(K),'timestp',varid)
  status = NF_PUT_VARA_INT(OUT_NCIDS(K),varid,start,count,FORCING_STEP-1)

  ! Initialize var index L
  L = 1

  CALL WRITE_OUT_VAR(OUT_NCIDS(K),VARIDS(K,L),FORCING_STEP,PotEvapPE,3,COMP_OUTPUT,LBAND)

  endif 
!! End hack
END

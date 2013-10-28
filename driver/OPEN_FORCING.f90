SUBROUTINE OPEN_FORCING()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! OPENS FORCING FILE

  ! Modifications:
  ! 2008-May-15 Prints error message and quits if can't open input file.	TJB

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: OPEN_FORCING.f90,v 1.3 2008/05/16 00:33:25 vicadmin Exp $"/

  ! Define local variables
  INTEGER ndims,nvars,ngatts,natts,unlimited
  INTEGER xid,yid,landid,tstepid
  INTEGER xlen_forcing,ylen_forcing,landlen_forcing
  INTEGER timestepvarid
  CHARACTER*20 varname
  INTEGER I

  FORCING_DT = 0
  Z = Z1

  ! Open forcing file and get data dimensions
  status = NF_OPEN(FORCFILE, 0, FORCING_NCID)
  IF (status .ne. NF_NOERR) THEN
    WRITE(*,*)'ERROR: cannot open forcing file',FORCFILE
    STOP
  ENDIF
  status = NF_INQ(FORCING_NCID, ndims, nvars, ngatts, unlimited)
  ! Check for compression, signified by presence of "land" dimension
  status = NF_INQ_DIMID(FORCING_NCID,'land',landid)
  IF (status .eq. NF_NOERR) THEN
    COMP_FORCING=.TRUE.
    status = NF_INQ_DIMLEN(FORCING_NCID,landid,landlen_forcing)
  ELSE
    status = NF_INQ_DIMID(FORCING_NCID,'x',xid)
    status = NF_INQ_DIMID(FORCING_NCID,'y',yid)
    status = NF_INQ_DIMLEN(FORCING_NCID,xid,xlen_forcing)
    status = NF_INQ_DIMLEN(FORCING_NCID,yid,ylen_forcing)
    landlen_forcing = 0
  END IF
  status = NF_INQ_DIMID(FORCING_NCID,'tstep',tstepid)
  status = NF_INQ_DIMLEN(FORCING_NCID,tstepid,tsteplen)

  status = NF_INQ_VARID(FORCING_NCID,'timestep',timestepvarid)
  IF (status .ne. NF_NOERR) THEN
   status = NF_INQ_VARID(FORCING_NCID,'timestp',timestepvarid)
  END IF
  status = NF_GET_ATT_INT(FORCING_NCID,timestepvarid,'tstep_sec',FORCING_DT)
  status = NF_GET_ATT_TEXT(FORCING_NCID,timestepvarid,'time_origin',FORCING_START_TIME)
  FORCING_DT_REAL = FORCING_DT

  IF (COMP_FORCING) THEN
    IF (landlen_forcing /= landlen) THEN
      WRITE(*,*)'ERROR: landlen',landlen_forcing, &
        ' (from forcing file ',TRIM(FORCING),') /= landlen ',landlen, &
        ' (from lsc file ',TRIM(LSC),')'
      STOP
    END IF
  ELSE
    IF (ylen_forcing /= ylen) THEN
      WRITE(*,*)'ERROR: ylen',ylen_forcing, &
        ' (from forcing file ',TRIM(FORCING),') /= ylen ',ylen, &
        ' (from lsc file ',TRIM(LSC),')'
      STOP
    END IF
    IF (xlen_forcing /= xlen) THEN
      WRITE(*,*)'ERROR: xlen',xlen_forcing, &
        ' (from forcing file ',TRIM(FORCING),') /= xlen ',xlen, &
        ' (from lsc file ',TRIM(LSC),')'
      STOP
    END IF
  END IF

  IF (FORCING_DT <= 0) THEN
    WRITE(*,*)'ERROR: FORCING_DT',FORCING_DT
    STOP
  END IF

  ! check presence of 2D wind field, large-scale & convective rainfall
  LWIND2D=.FALSE.
  LLSRAINF=.FALSE.
  status = NF_INQ_NVARS(FORCING_NCID, nvars)
  DO I=1,nvars
    status = NF_INQ_VARNAME(FORCING_NCID,I,varname)
    IF(varname(1:7).EQ.'LSRainf')THEN
      LLSRAINF=.TRUE.
    ENDIF
    IF(varname(1:6).EQ.'Wind_E')THEN
      LWIND2D=.TRUE.
    ENDIF
  ENDDO

END

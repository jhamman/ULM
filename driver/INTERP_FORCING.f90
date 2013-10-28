SUBROUTINE INTERP_FORCING(day_of_year,tstep_of_day,MODEL_STEP,I)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! INTERPOLATES FORCING DATA TO MODEL TIME STEP

  ! Modifications:
  ! 2007-Nov-12 Removed SOLNET computation (it was commented out anyway)	TJB
  ! 2008-May-05 Uncommented call to ZTERP and added handling of CZMODEL,
  !             for compatibility with NOAH 2.8.				Ben Livneh
  ! 2008-Jul-15 Commented out call to ZTERP and handling of CZMODEL again.	Ben Livneh

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: INTERP_FORCING.f90,v 1.6 2008/07/30 22:33:03 vicadmin Exp $"/

  ! Define local variables
  INTEGER, INTENT(IN) :: day_of_year,tstep_of_day,MODEL_STEP,I
  REAL    :: GMT, LHOUR, BTIME, ETIME, MTIME, WEIGHT1, WEIGHT2
  INTEGER :: ZONE

  GMT=MODEL_STEP*MODEL_DT/3600.+ (tstep_of_day-1)*FORCING_DT/3600
  IF(GMT.GE.24)GMT=GMT-24
  CALL LOCALTIME (GMT,LON(X(I),Y(I)),LHOUR,ZONE)
  CALL COSZENITH(LON(X(I),Y(I)),LAT(X(I),Y(I)),LHOUR,ZONE,day_of_year,CZMODEL(I))

  CALL FINTERP(0,'I',0,SWdownnow(I),SWdownplus(I),0,MADTT,MODEL_STEP,DT_SOLDNB(I))
  CALL FINTERP(0,'I',0,LWdownnow(I),LWdownplus(I),0,MADTT,MODEL_STEP,DT_LWDN(I))
  CALL FINTERP(0,'I',0,Tairnow(I),Tairplus(I),0,MADTT,MODEL_STEP,DT_TAIR(I))
  CALL FINTERP(0,'I',0,Qairnow(I),Qairplus(I),0,MADTT,MODEL_STEP,DT_SPFH(I))
  CALL FINTERP(0,'I',0,Snowfnow(I),Snowfplus(I),0,MADTT,MODEL_STEP,DT_SNOW(I))
  CALL FINTERP(0,'I',0,Rainfnow(I),Rainfplus(I),0,MADTT,MODEL_STEP,DT_RAIN(I))
  CALL FINTERP(0,'I',0,PSurfnow(I),PSurfplus(I),0,MADTT,MODEL_STEP,DT_PSFC(I))
  CALL FINTERP(0,'I',0,Windnow(I),Windplus(I),0,MADTT,MODEL_STEP,DT_WIND(I))

  BTIME=(tstep_of_day-1)*FORCING_DT/3600
  ETIME=tstep_of_day*FORCING_DT/3600
  IF(ETIME.GE.24)ETIME=ETIME-24
  MTIME=GMT

!  CALL ZTERP(1.0,LAT(X(I),Y(I)),LON(X(I),Y(I)),BTIME,ETIME,MTIME,day_of_year,WEIGHT1,WEIGHT2)

!  DT_SOLDN(I) = SWdownnow(I)*WEIGHT1 + SWdownplus(I)*WEIGHT2

!  ! Arbitrary cutoff scheme from Brian Cosgrove
!  IF ((DT_SOLDN(I).GT.SWdownnow(I)).AND.(DT_SOLDN(I).GT.SWdownplus(I))) THEN
!    IF (CZMODEL(I).LE.0.1) THEN
!      DT_SOLDN(I) = DT_SOLDNB(I)
!    ENDIF
!  ENDIF

  DT_SOLDN(I) = DT_SOLDNB(I)

  DT_PRCP(I)  = DT_RAIN(I)+DT_SNOW(I)

END

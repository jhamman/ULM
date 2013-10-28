        SUBROUTINE ZTERP (IFLAG,LAT,LON,BTIME,ETIME,
     &          MBTIME,JULIANB,weight1,weight2)

C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: ZTERP.f,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

**********************************************************************************
*
*       PROGRAMMER: Brian Cosgrove    ORGANIZATION: NASA/GSFC   DATE: 10/2/98
*
*
*       ABSTRACT:
*
*       This subroutine is based, in part, on modified subroutines from
*       Jean C. Morrill of the GSWP project.  The program temporally interpolates
*       time averge or instantaneous data to that needed by the model
*       at the current timestep.  It does this by combining a linear interpolation
*       approach with a solar zenith angle approach in a fashion suitable for use
*       with data such as short wave radiation values.
*       It cannot be used with input data points which are more
*       than 24 hours apart.  The program outputs two weights which can then be
*       applied to the original data to arrive at the interpolated data.  If
*       IFLAG=0, then WEIGHT1 is the weight which should be applied to the
*       time averaged data (from the time period which the model is currently in)
*       to arrive at the interpolated value and weight 2 is not used at all.  If
*       IFLAG=1, then WEIGHT1 should be applied to the original instantaneous
*       data located just prior to the model time step, and WEIGHT2 should be
*       applied to the original instantaneous data located just after the model
*       time step.
*       i.e.  (IF IFLAG=0)   interpolated data = (WEIGHT1 * time averaged data from
*                                                time period that model is currently in)
*       i.e.  (IF IFLAG=1) interp. data = (WEIGHT1*past data)+(WEIGHT2*future data)
**********************************************************************************
*
*
*       PROGRAM HISTORY LOG:
*       10/2/98 Brian Cosgrove
*
**********************************************************************************
*
*       INPUT ARGUMENT LIST:
*
*       IFLAG   -Int., flag specifing whether input data is
*               time averaged (IFLAG=0) or
*               instantaneous (IFLAG=1)
*       LAT     -Real, latitude (deg) of current data point
*       LON     -Real, longitude (deg) of current data point
*       BTIME   -Real, beginning time of orig. avg. data (IFLAG=0) or
*               time of original instantaneous data point which
*               is located at or just prior to the current model
*                time step (IFLAG=1).  Expects GMT time in hour
*               fractions.  i.e., 6:30 Z would be 6.5.
*       ETIME   -Real, endding time of orig. avg. data (IFLAG=0) or
*               time of original instantaneous data point which
*               is located at or just after the current model
*               time step (IFLAG=1).    Expects GMT time in hour
*               fractions.  i.e., 6:30 Z would be 6.5.
*       MBTIME  -Real, time of current model time step.
*               Expects GMT time in hour fractions
*               i.e., 6:30 Z would be 6.5.
*       JULIANB -Int., Julian day upon which BTIME falls
*
*
**********************************************************************************
*
*       OUTPUT ARGUMENT LIST:
*
*       WEIGHT1 -Real, weight applied to original time averaged
*               data (IFLAG=0) or
*               weight applied to orig instantaneous data point
*               located just prior to the current model time step
*       WEIGHT2 -Real, weight applied to orig instantaneous
*               data point located just after the current model
*               time step (IFLAG=1)
*               If IFLAG=0, then this weight is meaningless and
*               should not be used
*
**********************************************************************************

        REAL LAT,LON,MBTIME,INTERVAL
        REAL BTIME,ETIME,IFLAG,WEIGHT1,WEIGHT2
        REAL CZMODEL,CZAVGDATA,TOTANGLE,AVGANGLE
        REAL CZENDDATA,CZBEGDATA,LHOUR,GMT
        REAL LMBTIME,LETIME,LBTIME,WEIGHTE,WEIGHTB
        REAL DIFFEB,DIFFBM,DIFFEM
        INTEGER ZONEETIME,ZONEBTIME,ZONEMBTIME
        INTEGER ZONE,JULIANB,JULIANE,JULIANMB,JULIANTEMP,JULIANB1,
     &          JULIANE1,JULIANMB1





C       This section contains hardwired data that will be supplied by main program.
C       These values were chosen arbitrarily and exist simply to check the
C       functioning of the program.
C
C       Initialize variables
C

        I=1
        TOTANGLE=0
        WEIGHTE=-999.999
        WEIGHTB=-999.999
        WEIGHT1=-999.999
        WEIGHT2=-999.999
        CZBEGDATA=-999.999
        CZENDDATA=-999.999
        CZMODEL=-999.999
        GMT=BTIME
        JULIANE=JULIANB
        JULIANMB=JULIANB
        JULIANTEMP=JULIANB

        IF (MBTIME.LT.BTIME) JULIANMB=JULIANMB+1
        IF (ETIME.LE.BTIME) JULIANE=JULIANE+1

C       First case, IFLAG=0 (Time average input, instantaneous output)
C
C       Compute time interval, here arbitrarily divided into 36 parts

        IF (IFLAG.EQ.0) THEN

                CALL LOCALTIME (MBTIME,LON,LHOUR,ZONE)
                CALL COSZENITH(LON,LAT,LHOUR,ZONE,JULIANMB,CZMODEL)
                IF (CZMODEL.EQ.0) THEN
                        WEIGHT1=0
                        GOTO 818
                ENDIF

                IF (ETIME.GT.BTIME) THEN
                        INTERVAL = ((ETIME-BTIME)/36.0)
                ELSEIF (ETIME.LT.BTIME) THEN
                        INTERVAL = (((24-BTIME)+(ETIME))/36.0)
                ELSE
                        INTERVAL=24.0/36.0
                ENDIF

C       Compute cosine of zenith angle for each time interval
C
                DO WHILE (I.LE.37)
                        IF ((GMT+INTERVAL).LT.24) THEN
                                GMT=GMT+INTERVAL
                        ELSE
                                GMT=(INTERVAL-(24-GMT))
                                JULIANTEMP=JULIANTEMP+1
                        ENDIF
                        CALL LOCALTIME (GMT,LON,LHOUR,ZONE)
                        CALL COSZENITH(LON,LAT,LHOUR,ZONE,
     &                     JULIANTEMP,CZAVGDATA)
                        TOTANGLE=TOTANGLE+CZAVGDATA
                        I=I+1

                ENDDO

C       Compute average cosine of zenith angle and also
C       weight which will be applied to original data (WEIGHT1)
C
                AVGANGLE=(TOTANGLE/37.0)
                IF (AVGANGLE.EQ.0) THEN
                        WEIGHT1=0
                ELSE
                WEIGHT1=(CZMODEL/AVGANGLE)
                ENDIF

        ENDIF


C       Second case:  IFLAG=1 (instantaneous input and output)
        IF (IFLAG.eq.1) THEN
C
C       Compute local times and cosine (zenith angle)
C
                CALL LOCALTIME (BTIME,LON,LBTIME,ZONEBTIME)
C               IF (LBTIME.GT.BTIME) JULIANB=JULIANB-1
                IF (LBTIME.GT.BTIME) THEN
                 JULIANB1=JULIANB-1
                ELSE
                 JULIANB1=JULIANB
                ENDIF
                CALL LOCALTIME (ETIME,LON,LETIME,ZONEETIME)
C               IF (LETIME.GT.ETIME) JULIANE=JULIANE-1
                IF (LETIME.GT.ETIME) THEN
                 JULIANE1=JULIANE-1
                ELSE
                 JULIANE1=JULIANE
                ENDIF
                CALL LOCALTIME (MBTIME,LON,LMBTIME,ZONEMBTIME)
C               IF (LMBTIME.GT.MBTIME) JULIANMB=JULIANMB-1
                IF (LMBTIME.GT.MBTIME) THEN
                 JULIANMB1=JULIANMB-1
                ELSE
                 JULIANMB1=JULIANMB
                ENDIF

                CALL COSZENITH (LON,LAT,LBTIME,ZONEBTIME,
C    &              JULIANB,CZBEGDATA)
     &              JULIANB1,CZBEGDATA)
                CALL COSZENITH (LON,LAT,LETIME,ZONEETIME,
C    &              JULIANE,CZENDDATA)
     &              JULIANE1,CZENDDATA)
                CALL COSZENITH (LON,LAT,LMBTIME,ZONEMBTIME,
C    &              JULIANMB,CZMODEL)
     &              JULIANMB1,CZMODEL)
c       print*,'CZBEGDATA=',CZBEGDATA,'LBTIME=',LBTIME,LON
c       print*,'CZENDDATA=',CZENDDATA,'LETIME=',LETIME,LON
c       print*,'CZMODEL=',CZMODEL,'LMBTIME=',LMBTIME,LON

C       Decision tree to deal with contingencies
C       If COS(zenith angle at current model time =0, weight =0
C       If COS(zenith angle =0 at beg. and end times, PROBLEM, STOP
C       Otherwise use beginning and ending data to calculate weight
C
                IF (CZMODEL.EQ.0) THEN
                        WEIGHT1=0
                        WEIGHT2=0
                ELSE
                IF ((CZBEGDATA.NE.0).or.(CZENDDATA.NE.0)) THEN
                        IF (CZBEGDATA.EQ.0) THEN
                                WEIGHT1=0
                                WEIGHT2=(CZMODEL/CZENDDATA)
                        ENDIF
                        IF (CZENDDATA.EQ.0) THEN
                                WEIGHT1=(CZMODEL/CZBEGDATA)
                                WEIGHT2=0
                        ENDIF
                        IF ((CZENDDATA.NE.0).and.(CZBEGDATA.NE.0)) THEN

                        IF (BTIME.LE.MBTIME) THEN
                                DIFFBM=MBTIME-BTIME
                        ELSE
                                DIFFBM=24-BTIME+MBTIME
                        ENDIF

                        IF (ETIME.GE.MBTIME) THEN
                                DIFFEM=ETIME-MBTIME
                        ELSE
                                DIFFEM=24-MBTIME+ETIME
                        ENDIF

                        IF (ETIME.GT.BTIME) THEN
                                DIFFEB=ETIME-BTIME
                        ELSEIF (ETIME.EQ.BTIME) THEN

                                DIFFEB=24
                        ELSE
                                DIFFEB=24-BTIME+ETIME
                        ENDIF
                        WEIGHTE=(DIFFBM/DIFFEB)
                        WEIGHTB=(DIFFEM/DIFFEB)

                        WEIGHT1=((CZMODEL/CZBEGDATA)*WEIGHTB)
                        WEIGHT2=((CZMODEL/CZENDDATA)*WEIGHTE)

                        ENDIF
                ELSE
                        print *,'NO DATA TO INTERPOLATE TO/FROM'
                        print *,'BEGINNING AND ENDING DATA BOTH = 0'
                        print *,'STOPPING!!!'
        print*,'CZBEGDATA=',CZBEGDATA,'LBTIME=',LBTIME,LON
        print*,'CZENDDATA=',CZENDDATA,'LETIME=',LETIME,LON
        print*,'CZMODEL=',CZMODEL,'LMBTIME=',LMBTIME,LON
c                       STOP
                ENDIF
                ENDIF




        ENDIF
 818    CONTINUE
        RETURN
        END

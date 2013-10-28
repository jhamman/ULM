      SUBROUTINE COSZENITH (LON,LATD,LHOUR,ZONE,JULIAN,CZENITH)
C
C
C     The purpose is to calculate the following:
C
C        1)  Day angle (GAMMA)
C
C        2)  Solar DEClination
C
C        3)  Equation of time
C
C        4)  Local apparent time
C
C        5)  Hour angle
C
C        6)  Cosine of zenith angle
C
C     All equations come from "An Introduction to
C     Solar Radition" By Muhammad Iqbal, 1983.
C
C
C-----------------------------------------------------------------------
      implicit none

C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: COSZENITH.f,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

C---------------------------Local variable------------------------------
C
      INTEGER
     $     ZONE,                ! time zone (1-24) GMT=12
     $     JULIAN               ! julian day

C
      REAL
     $     CZENITH,             ! cosine of zenith angle (radians)
     $     DECd,                ! solar declination (degrees)
     $     DEC,                 ! solar declination (radians)
     $     et,                  ! equation of time (minutes)
     $     GAMMA,               ! day angle (radians)
     $     LATime,              ! local apparent time
     $     LCORR,               ! longitudical correction
     $     LHOUR,               ! local standard time
     $     LON,                 ! local longitude (deg)
     $     LLAT,                ! local latitude in radians
     $     LATD ,               ! local latitude in degrees
     $     LS,                  ! standard longitude (deg)
     $     OMEGAD,              ! omega in degrees
     $     OMEGA ,              ! omega in radians
     $     PI,                  ! universal constant PI [-]
     $     ZENITH,              ! zenith angle(radians)
     $     ZEND                 ! zenith angle(degress)
C
C     Neither ZENITH nor ZEND are necessary for this program.
C     I originally used them as checks, and left them here in
C     case anyone else had a use for them.
************************************************************************
C
C     1)  Day angle GAMMA (radians) page 3

      PI= 3.141592              ! universal constant PI
      GAMMA=2*PI*(JULIAN-1)/365.
************************************************************************
C     2) Solar declination (assumed constant for a 24 hour period)  page 7
C     in radians
C
      DEC=(0.006918-0.399912*COS(GAMMA)+0.070257*SIN(GAMMA)
     $     -0.006758*COS(2*GAMMA)+0.000907*SIN(2*GAMMA)
     $     -0.002697*COS(3*GAMMA)+0.00148*SIN(3*GAMMA))
      DECd=DEC*(180./PI)
C
C     maximum error 0.0006 rad (<3'), leads to error of less than 1/2 degree
C     in ZENITH angle
************************************************************************^M
C     3)  Equation of time  page 11

      et=(0.000075+0.001868*COS(GAMMA)-0.032077*SIN(GAMMA)
     $     -0.014615*COS(2*GAMMA)-0.04089*SIN(2*GAMMA))*229.18
C
************************************************************************^M
C     4) Local apparent time  page 13
C
C     LS     standard longitude (nearest 15 degree meridian)
C     LON     local longitude
C     LHOUR  local standard time
C     LATIME local apparent time
C     LCORR  longitudunal correction (minutes)
C
      LS=((ZONE-1)*15)-180.
      LCORR=4.*(LS-LON)*(-1)
      LATIME=LHOUR+LCORR/60.+et/60.
      IF (LATIME.LT.0.) LATIME=LATIME+24
      IF (LATIME.GT.24.) LATIME=LATIME-24
************************************************************************
C     5) Hour angle OMEGA  page 15
C
C     hour angle is zero at noon, postive in the morning
C     It ranges from 180 to -180
C
C     OMEGAD is OMEGA in degrees, OMEGA is in radians
C
      OMEGAD=(LATime-12.)*(-15.)
      OMEGA=OMEGAD*PI/180.
C
************************************************************************
C     6)  Zenith angle page 15
C
C     CZENITH cosine of zenith angle (radians)
C     LATD=local latitude in degrees
C     LLAT=local latitude in radians
C
      LLAT=LATD*PI/180.
      CZENITH=SIN(DEC)*SIN(LLAT)+COS(DEC)*COS(LLAT)*COS(OMEGA)
      CZENITH=AMAX1(0.,CZENITH)
      ZENITH=ASIN(CZENITH)
      ZEND=180.*ZENITH/PI

************************************************************************
C
      RETURN
      END

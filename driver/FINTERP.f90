      SUBROUTINE finterp(ip,trp_flag,val0,val1,val2,val3,madtt,istep,vout)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

!======================================================================
! Description: INTERPOLATION OF ATMOSPHERIC FORCING DATA
!
! Revision 1.2  2002/09/20  20:55:32  dirmeyer
! Introduced a new interpolation type "P" which specifies a temporal
! disaggregation of interpolated precipitation to improve partitioning
! of infiltration/runoff when forcing data interval is large compared to
! the typical length of a convective (or total) rain event.
!
! Revision 1.1  2002/08/16  15:55:18  guo
! Initial revision
!
!======================================================================
 
      !** This routine interpolates atmospheric adtt seconds period forcing
      !   data to dtt seconds valtrp time-step for one adtt 
      !   interval. Interpolation is performed based on the value of 
      !   'trp_flag':
      !
      !   "L" or "l" = value represents average over interval ending at current time
      !   "N" or "n" = value represents average over interval beginning at current time
      !   "C" or "c" = value represents average over interval centered on current time
      !   "I" or "i" = instantaneous value at current time (linear interpolation)
      !   "P" or "p" = PDF disaggregation applied in time (for precip)
      !   "0" (zero) = no interpolation, centered on current time
      !   Otherwise  = no interpolation, applied beginning at current time
      !
      !
      !
      !************************************************************************
      !* NOTE!!!!!!!
      !*   For "L", "N", and "C" to conserve the mean over the interval after 
      !*   interpolation, 'madtt' MUST BE A MULTIPLE OF 2!
      !************************************************************************
      !
      !   Algorithm for conserving interpolation scheme - courtesy of 
      !   Dag Lohmann, NOAA/NCEP/EMC
      !
      !   - P. Dirmeyer, November 2001
      !===================================================================
 
      IMPLICIT NONE
  
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: FINTERP.f90,v 1.2 2007/10/04 20:44:20 vicadmin Exp $"/

! Input 
      INTEGER,INTENT(IN)           :: ip    ! Grid point in question
      INTEGER,INTENT(IN)           :: istep    ! step number
      CHARACTER(len=1),INTENT(IN)  :: trp_flag ! Type of interpolation

      REAL,INTENT(IN)              :: val0   ! last forcing value (for cases N & C)
      REAL,INTENT(IN)              :: val1   ! current forcing value (all cases)
      REAL,INTENT(IN)              :: val2   ! next forcing value (all cases but default)
      REAL,INTENT(IN)              :: val3   ! ueber-next forcing value (for cases C & L)

      INTEGER,INTENT(IN)           :: madtt  ! Number of valtrp timesteps in
                                             ! a forcing timestep.
! Output
      REAL, DIMENSION(madtt)       :: valtrp ! Interpolated forcing data vector
      REAL                         :: vout

! Local
      REAL               :: fac0 ! Weight for last value 
      REAL               :: fac1 ! Weight for current value
      REAL               :: fac2 ! Weight for next value
      REAL               :: rmadtt ! real madtt
      REAL               :: denom  ! denominator of scaling factor for 
                                   ! conserving interpolation
      REAL               :: numer  ! numerator of scaling factor for 
                                   ! conserving interpolation
      INTEGER            :: i
      INTEGER            :: j
      INTEGER            :: ip1  ! i + 1
      INTEGER            :: im1  ! i - 1

      REAL               :: rtdist(120) ! Weights for precip disag PDF
      REAL               :: factor      ! Scaling factor for brevity/intensity
      REAL               :: expo        ! Exponent to fit slope of log-log relationship
      REAL               :: rtsum       ! Ensure PDF weights add to 1.0
      REAL               :: rtsum0      ! Sum from last time interval
      REAL               :: p0, p1      ! Define edges of each time bin
      REAL               :: rintr1      ! Log-log tail target value

      LOGICAL            :: iniflag
      SAVE iniflag,rtdist
      DATA iniflag/.true./


      !===================================================================
      rmadtt = float(madtt)

!
!>>> Generate disaggregation PDF - Just do this once
!
      IF (iniflag) THEN
         write(6,4030)
 4030    format('Initializing precipitation disaggregation PDF:')
         factor = 66.6666  ! For 33% dry time steps in a rainy forcing interval
!         factor = 29.63  ! For 0% dry time steps in a rainy forcing interval
         expo = -3.0      
         rtsum = 0.0
         rintr1 = factor * float(madtt) ** expo

         DO j = 1, madtt
           rtsum0 = rtsum
           p0 = float(j-1)/float(madtt)
           p1 = float(j)/float(madtt)
           rtdist(j) = rintr1*expo/((expo+1)*(p1-p0)*100)* ((100*p1/rintr1)**((expo+1)/expo)- (100*p0/rintr1)**((expo+1)/expo))
           rtsum = rtsum + rtdist(j)
           IF (rtsum > 1.0) THEN
             IF (rtsum0 < 1.0) THEN
               rtdist(j) = rtdist(j) + 1.0 - rtsum   ! Stick any remaining weight in last interval
             ELSE
               rtdist(j) = 0.0
             ENDIF
           ENDIF
           write(6,4050) j,rtdist(j)
 4050      format('Time disag: ',i3,': ',f11.6)
         ENDDO
         IF (rtsum < 1.0) THEN
             rtdist(1) = rtdist(1) + 1.0 - rtsum  ! Ensure integral of PDF = 1.0
         ENDIF
      ENDIF
      iniflag = .false.

      IF (trp_flag == 'I' .or. trp_flag == 'i') THEN
      !** Current value valid at midpoint of interval
        DO j = 1, madtt
          fac1     = float(madtt+1-j)/madtt
          fac2     = 1.0-fac1
          valtrp(j)   = val1*fac1+val2*fac2
        ENDDO

      ELSEIF (trp_flag == 'N' .or. trp_flag == 'n') THEN
      !** Current value is average over next interval
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.0*rmadtt-abs(float(2*j-madtt-1)))/(rmadtt*2.0)
          fac0 = max(1.0-float(j*2+madtt-1)/(rmadtt*2.0),0.0)
          fac2 = max(1.0-float((madtt+1-j)*2+madtt-1)/(rmadtt*2.0),0.0)
          denom = 0.5*(val0+val2)+3.0*val1
          numer = 4.0*val1
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val0*fac0+val1*fac1+val2*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'C' .or. trp_flag == 'c') THEN
      !** Current value is average centered on current time
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.0*rmadtt-(float(2*j-1)))/(rmadtt*2.0)
          fac0 = 0.0
          fac2 = 1.0-fac1
          IF (j > madtt/2) THEN
            denom = 0.5*(val1+val3)+3.0*val2
            numer = 4.0*val2
          ELSE
            denom = 0.5*(val0+val2)+3.0*val1
            numer = 4.0*val1
          ENDIF
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val0*fac0+val1*fac1+val2*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'L' .or. trp_flag == 'l') THEN
      !** Current value is average for period ending at current time
        IF (mod(madtt,2) /= 0) STOP 'drv_interp mod(madtt,2) /= 0'
        DO j = 1, madtt
          fac1 = (2.0*rmadtt-abs(float(2*j-madtt-1)))/(rmadtt*2.0)
          fac0 = max(1.0-float(j*2+madtt-1)/(rmadtt*2.0),0.0)
          fac2 = max(1.0-float((madtt+1-j)*2+madtt-1)/(rmadtt*2.0),0.0)
          denom = 0.5*(val1+val3)+3.0*val2
          numer = 4.0*val2
          IF (denom > EPSILON(denom)) THEN
             valtrp(j) = (val1*fac0+val2*fac1+val3*fac2) * numer / denom
          ELSE
             valtrp(j) = 0.0
          ENDIF
        ENDDO

      ELSEIF (trp_flag == 'P' .or. trp_flag == 'p') THEN
      !** Current value is from a PDF that disagregates in time over the forcing interval
        DO j = 1, madtt
          valtrp(j) = MAX(val2*rtdist(j)*float(madtt),0.0)
        ENDDO

      ELSEIF (trp_flag == '0') THEN
      !** Current value is applied centered on current time without interpolation
        DO j = 1, madtt
          IF (j > madtt/2) THEN
            valtrp(j) = val2
          ELSE
            valtrp(j) = val1
          ENDIF
        ENDDO

      ELSE
      !** Current value is applied over this interval without interpolation
        DO j = 1, madtt
          valtrp(j) = val1
        ENDDO

      ENDIF

      vout=valtrp(istep)

      END SUBROUTINE finterp

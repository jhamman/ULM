SUBROUTINE READ_LSC()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS MODEL PARAMETERS

  ! Modifications:
  ! 2007-Nov-15 Added quality control on SHDFAC.				TJB
  ! 2008-May-05 Added SHDMAX for compatibility with NOAH 2.8.			Ben Livneh
  ! 2008-Oct-07 Included SAC parameters for use as unified model         Ben Livneh
  ! 2009-Feb-03 Made NSOIL = 2 for any SAC model case                   B.L.
  ! 2009-Feb-03 Made a temporary definition of SOILDEPTH_ACCUM for      B.L.
  ! SAC model cases, as a function of porosity and SAC params.
  ! 2009-Feb-03 Moved SOILTYP and SAC param call to top for depth computation
  ! 2009-Feb-04 Changed the dimension of TEMP array for 2-d SAC vars to 2-d  B.L.
  ! 

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_LSC.f90,v 1.3 2008/05/12 22:58:55 vicadmin Exp $"/

  ! Define local variables
  INTEGER J,I,K
  INTEGER land_idx
  INTEGER TEMP_INT(xlen,ylen),TEMP_DEPTH_INT(MAXNSOIL,xlen,ylen)
  REAL    TEMP(xlen,ylen),TEMP_DEPTH(MAXNSOIL,xlen,ylen),TEMP_TIME(NMONTHS,xlen,ylen)
  REAL    cumulative_depth(landlen)
  REAL    tmp_PCTIM,tmp_ADIMP,tmp_RIVA,tmp_SIDE,tmp_RSERV
  REAL    tmp_WCRIT,tmp_RICHARDS
  REAL    tmp_PSNOW1,tmp_PSNOW2
  REAL    MAXSMC(19)
  REAL    WLTSMC(19)
  REAL    REFSMC(19)
  CHARACTER*200 FILENAME_TMP
  CHARACTER*2   K_STR

  DATA MAXSMC/0.395, 0.421, 0.434, 0.476, 0.476, 0.439,&
       0.404, 0.464, 0.465, 0.406, 0.468, 0.457,&
       0.464, 0.000, 0.200, 0.421, 0.457, 0.200,&
       0.395/
  DATA WLTSMC/0.023, 0.028, 0.047, 0.084, 0.084, 0.066,&
       0.069, 0.120, 0.103, 0.100, 0.126, 0.135,&
       0.069, 0.000, 0.012, 0.028, 0.135, 0.012,&
       0.023/
  DATA REFSMC/0.196, 0.248, 0.282, 0.332, 0.332, 0.301,&
       0.293, 0.368, 0.361, 0.320, 0.388, 0.389,&
       0.319, 0.000, 0.116, 0.248, 0.389, 0.116,&
       0.196/

  IF (trim(PARAM_TYPE) == 'netcdf') THEN

    ! NetCDF file - contains ROW,COL,LAT,LON,CELLID and NOAH parameters that vary spatially
    CALL READ_LSC_NETCDF

  ELSE

    ! Compute ROW, COL, LAT, LON, CELLID, NSOIL
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          ROW(J,I) = I
          COL(J,I) = J
          LAT(J,I) = YLLCORNER + (I-0.5)*cellsize
          LON(J,I) = XLLCORNER + (J-0.5)*cellsize
          CELLID(J,I) = land_idx
          land_idx = land_idx + 1
        ELSE
          ROW(J,I) = NODATA_INT
          COL(J,I) = NODATA_INT
          LAT(J,I) = NODATA
          LON(J,I) = NODATA
          CELLID(J,I) = NODATA_INT
        END IF
      END DO
    END DO

    ! Get NOAH parameters that vary spatially
    ! For all SAC cases make NSOIL = 2


    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_INT_ASC(SOILTYP_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP_INT)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_INT_SPATIAL(SOILTYP_FILE,1,xlen,ylen,TEMP_INT)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          SOILTYP(land_idx) = TEMP_INT(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO
  
    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_INT_ASC(NSOIL_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP_INT)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_INT_SPATIAL(NSOIL_FILE,1,xlen,ylen,TEMP_INT)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
           NSOIL(land_idx) = TEMP_INT(J,I)
           land_idx = land_idx + 1
        END IF
      END DO
    END DO

    ! Get SAC parameters that vary spatially; needed for depth computation
  
    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(UZTWM_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(UZTWM_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          UZTWM_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(UZFWM_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(UZFWM_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          UZFWM_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(UZK_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(UZK_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          UZK_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(ZPERC_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(ZPERC_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          ZPERC_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(REXP_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(REXP_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          REXP_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(LZTWM_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(LZTWM_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          LZTWM_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(LZFPM_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(LZFPM_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          LZFPM_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(LZFSM_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(LZFSM_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          LZFSM_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(LZPK_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(LZPK_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          LZPK_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(LZSK_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(LZSK_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          LZSK_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(PFREE_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(PFREE_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          PFREE_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

  ! Get PE scale factor

!  IF (trim(PARAM_TYPE) == 'ascii') THEN
!    CALL READ_REAL_ASC(PE_SCALE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
!  ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
!    CALL READ_REAL_SPATIAL(PE_SCALE,1,xlen,ylen,TEMP)
!  END IF
!  land_idx = 1
!  DO I = 1,ylen
!    DO J = 1,xlen
!      IF (LANDMASK(J,I) == 1) THEN
!        PESCALE_2d(land_idx) = TEMP(J,I)
!        land_idx = land_idx + 1
!      END IF
!    END DO
!  END DO

  ! Get PE adjustment factors

!  IF (trim(PARAM_TYPE) == 'ascii') THEN
!    DO K = 1,NMONTHS
!      WRITE(K_STR,FMT='(I2.2)') K
!      WRITE(FILENAME_TMP,*) TRIM(PE_ADJ)//K_STR
!      CALL READ_REAL_ASC(FILENAME_TMP,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
!      land_idx = 1
!      DO I = 1,ylen
!        DO J = 1,xlen
!          IF (LANDMASK(J,I) == 1) THEN
!            PEADJ_2d(K,land_idx) = TEMP(J,I)
!            land_idx = land_idx + 1
!          END IF
!        END DO
!      END DO
!    END DO
!  ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
!    CALL READ_REAL_SPATIAL(PE_ADJ,NMONTHS,xlen,ylen,TEMP_TIME)
!    land_idx = 1
!    DO I = 1,ylen
!      DO J = 1,xlen
!        IF (LANDMASK(J,I) == 1) THEN
!          DO K = 1,NMONTHS
!            PEADJ_2d(K,land_idx) = TEMP_TIME(K,J,I)
!          END DO
!          land_idx = land_idx + 1
!        END IF
!      END DO
!    END DO
!  END IF


  ! Now compute soil depths based on SAC zone capacities

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      DO land_idx = 1,landlen
        cumulative_depth(land_idx) = 0
      END DO
      DO K = 1,MAXNSOIL
        WRITE(K_STR,FMT='(I2.2)') K
        WRITE(FILENAME_TMP,*) TRIM(SOILDEPTH_FILE)//TRIM(K_STR)
        CALL READ_REAL_ASC(FILENAME_TMP,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
        land_idx = 1
        DO I = 1,ylen
          DO J = 1,xlen
            IF (LANDMASK(J,I) == 1) THEN
               ! SAC cases -- make same definition twice (redundant) to keep iteration framework clean
               ! depths in metres to be consistent with native Noah
               !write(*,*)'filename_tmp',FILENAME_TMP,'TEMP(15,1)',TEMP(15,1)
               SOILDEPTH(land_idx,K) = TEMP(J,I)
               cumulative_depth(land_idx) = cumulative_depth(land_idx) + SOILDEPTH(land_idx,K)
               SOILDEPTH_ACCUM(land_idx,K) = cumulative_depth(land_idx)
               land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(SOILDEPTH_FILE,MAXNSOIL,xlen,ylen,TEMP_DEPTH,LREAL8)
      land_idx = 1
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            cumulative_depth(land_idx) = 0
            DO K = 1,MAXNSOIL
              SOILDEPTH(land_idx,K) = TEMP_DEPTH(K,J,I)
              cumulative_depth(land_idx) = cumulative_depth(land_idx) + SOILDEPTH(land_idx,K)
              SOILDEPTH_ACCUM(land_idx,K) = cumulative_depth(land_idx)
            END DO
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END IF
   ! write(*,*)'DEFINED SOILDEPTHS',SOILDEPTH(1,:),TEMP(15,1)
    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_INT_ASC(SLOPETYP_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP_INT)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_INT_SPATIAL(SLOPETYP_FILE,1,xlen,ylen,TEMP_INT)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          SLOPETYP(land_idx) = TEMP_INT(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(TBOT_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(TBOT_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          TBOT(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_INT_ASC(VEGTYP_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP_INT)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_INT_SPATIAL(VEGTYP_FILE,1,xlen,ylen,TEMP_INT)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          VEGTYP(land_idx) = TEMP_INT(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      DO K = 1,NMONTHS
        WRITE(K_STR,FMT='(I2.2)') K
        WRITE(FILENAME_TMP,*) TRIM(SHDFAC_FILE)//K_STR
        CALL READ_REAL_ASC(FILENAME_TMP,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
        land_idx = 1
        DO I = 1,ylen
          DO J = 1,xlen
            IF (LANDMASK(J,I) == 1) THEN
              SHDFAC(land_idx,K) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(SHDFAC_FILE,NMONTHS,xlen,ylen,TEMP_TIME,LREAL8)
      land_idx = 1
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            DO K = 1,NMONTHS
              SHDFAC(land_idx,K) = TEMP_TIME(K,J,I)
            END DO
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END IF

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      DO K = 1,NMONTHS
        WRITE(K_STR,FMT='(I2.2)') K
        WRITE(FILENAME_TMP,*) TRIM(ALBEDO_FILE)//K_STR
        CALL READ_REAL_ASC(FILENAME_TMP,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
        land_idx = 1
        DO I = 1,ylen
          DO J = 1,xlen
            IF (LANDMASK(J,I) == 1) THEN
              ALBEDO(land_idx,K) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(ALBEDO_FILE,NMONTHS,xlen,ylen,TEMP_TIME,LREAL8)
      land_idx = 1
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            DO K = 1,NMONTHS
              ALBEDO(land_idx,K) = TEMP_TIME(K,J,I)
            END DO
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END IF

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(SNOALB_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(SNOALB_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          SNOALB(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    IF (trim(PARAM_TYPE) == 'ascii') THEN
      CALL READ_REAL_ASC(ELEV_FILE,1,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)
    ELSE IF (trim(PARAM_TYPE) == 'bin') THEN
      CALL READ_REAL_SPATIAL(ELEV_FILE,1,xlen,ylen,TEMP,LREAL8)
    END IF
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          ELEV_2d(land_idx) = TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

  END IF

  ! The remaining SAC parameters are constant
  OPEN(99,FILE=SAC_CONST,STATUS='OLD')
  READ (99,*)
  READ (99,*) tmp_PCTIM,tmp_ADIMP,tmp_RIVA,tmp_SIDE,tmp_RSERV,tmp_WCRIT, &
       tmp_PSNOW1,tmp_PSNOW2,tmp_RICHARDS
  CLOSE(99)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        PCTIM_2d(land_idx) = tmp_PCTIM
        ADIMP_2d(land_idx) = tmp_ADIMP
        RIVA_2d(land_idx) = tmp_RIVA
        SIDE_2d(land_idx) = tmp_SIDE
        RSERV_2d(land_idx) = tmp_RSERV
        WCRIT_2d(land_idx) = tmp_WCRIT
        PSNOW1_2d(land_idx) = tmp_PSNOW1
        PSNOW2_2d(land_idx) = tmp_PSNOW2
        RICHARDS_2d(land_idx) = tmp_RICHARDS
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Check SHDFAC
  DO I = 1,landlen
    DO K = 1,NMONTHS
      IF (SHDFAC(I,K) < 0) THEN
        SHDFAC(I,K) = 0
      END IF
    END DO
  END DO

  ! Compute SHDMIN and SHDMAX
  SHDMAX = -999.0
  SHDMIN = 999.0
  DO I = 1,landlen
    DO K = 1,NMONTHS
      SHDMAX(I) = max(SHDMAX(I),SHDFAC(I,K))
      SHDMIN(I) = min(SHDMIN(I),SHDFAC(I,K))
    END DO
  END DO

  ! Compute ICE
  ICE = ICE1

  ! Compute PTU
  PTU = 0.0

END

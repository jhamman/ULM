SUBROUTINE WRITE_RESTART()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! WRITES RESTART FILE

  ! Modifications:
  ! 2007-Nov-06 Added SnowSoilT (T12).						Ben Livneh
  ! 2007-Nov-12 Added LSTSNW.							TJB
  ! 2008-May-05 Removed SnowSoilT (T12), for compatibility with NOAH 2.8.	Ben Livneh
  ! 2008-Jul-15 Added SnowTProf (TPACK).					Ben Livneh
  ! 2008-Jul-24 Added PACH20.							Ben Livneh
                                                                                
  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE
                                                                                
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: WRITE_RESTART.f90,v 1.6 2008/07/30 22:33:03 vicadmin Exp $"/

  ! Define local variables
  INTEGER varid
  INTEGER J,I,K,L,start2d(2),count2d(2),start3d(3),count3d(3),start4d(4),count4d(4)
  REAL    TEMP(xlen,ylen)
  INTEGER TEMP_INT(xlen,ylen)
  INTEGER land_idx

  ! Write Soil Moisture Content
  status = NF_INQ_VARID(RESTART_NCID,'SMC',varid)
  DO L = 1, MAXNSOIL
    DO K = 1, NBANDS
      start4d(1) = 1
      start4d(2) = 1
      start4d(3) = K
      start4d(4) = L
      count4d(1) = xlen
      count4d(2) = ylen
      count4d(3) = 1
      count4d(4) = 1
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            TEMP(J,I) = SMC(land_idx,K,L)
            land_idx = land_idx + 1
          ELSE
            TEMP(J,I) = NODATA
          END IF
        END DO
      END DO
      status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start4d,count4d,TEMP)
    END DO
  END DO

  ! Write Liquid Soil Moisture Content
  status = NF_INQ_VARID(RESTART_NCID,'SH2O',varid)
  DO L = 1, MAXNSOIL
    DO K = 1, NBANDS
      start4d(1) = 1
      start4d(2) = 1
      start4d(3) = K
      start4d(4) = L
      count4d(1) = xlen
      count4d(2) = ylen
      count4d(3) = 1
      count4d(4) = 1
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            TEMP(J,I) = SH2O(land_idx,K,L)
            land_idx = land_idx + 1
          ELSE
            TEMP(J,I) = NODATA
          END IF
        END DO
      END DO
      status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start4d,count4d,TEMP)
    END DO
  END DO

  ! Write Soil Temperature
  status = NF_INQ_VARID(RESTART_NCID,'STC',varid)
  DO L = 1, MAXNSOIL
    DO K = 1, NBANDS
      start4d(1) = 1
      start4d(2) = 1
      start4d(3) = K
      start4d(4) = L
      count4d(1) = xlen
      count4d(2) = ylen
      count4d(3) = 1
      count4d(4) = 1
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            TEMP(J,I) = STC(land_idx,K,L)
            land_idx = land_idx + 1
          ELSE
            TEMP(J,I) = NODATA
          END IF
        END DO
      END DO
      status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start4d,count4d,TEMP)
    END DO
  END DO

  ! Write SAC states
  status = NF_INQ_VARID(RESTART_NCID,'SACST',varid)
  DO L = 1, 6
    DO K = 1, NBANDS
      start4d(1) = 1
      start4d(2) = 1
      start4d(3) = K
      start4d(4) = L
      count4d(1) = xlen
      count4d(2) = ylen
      count4d(3) = 1
      count4d(4) = 1
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            TEMP(J,I) = SACST(land_idx,K,L)
            land_idx = land_idx + 1
          ELSE
            TEMP(J,I) = NODATA
          END IF
        END DO
      END DO
      status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start4d,count4d,TEMP)
    END DO
  END DO

  ! Write Frozen SAC states
  status = NF_INQ_VARID(RESTART_NCID,'FRZST',varid)
  DO L = 1, 10
    DO K = 1, NBANDS
      start4d(1) = 1
      start4d(2) = 1
      start4d(3) = K
      start4d(4) = L
      count4d(1) = xlen
      count4d(2) = ylen
      count4d(3) = 1
      count4d(4) = 1
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            TEMP(J,I) = FRZST(land_idx,K,L)
            land_idx = land_idx + 1
          ELSE
            TEMP(J,I) = NODATA
          END IF
        END DO
      END DO
      status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start4d,count4d,TEMP)
    END DO
  END DO

  ! Write Skin Temperature
  status = NF_INQ_VARID(RESTART_NCID,'T1',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = T1(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Snow Pack Temperature
  status = NF_INQ_VARID(RESTART_NCID,'TPACK',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = TPACK(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Snow Pack Liquid Water Content
  status = NF_INQ_VARID(RESTART_NCID,'PACH20',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = PACH20(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Canopy Moisture
  status = NF_INQ_VARID(RESTART_NCID,'CMC',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = CMC(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Soil Moisture State
  status = NF_INQ_VARID(RESTART_NCID,'UZTWC',varid)
  DO L = 1, NBANDS
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = UZTWC_2d(land_idx,L)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = 1.e20
        END IF
      END DO
    END DO
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = L
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  status = NF_INQ_VARID(RESTART_NCID,'UZFWC',varid)
  DO L = 1, NBANDS
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = UZFWC_2d(land_idx,L)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = 1.e20
        END IF
      END DO
    END DO
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = L
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  status = NF_INQ_VARID(RESTART_NCID,'LZTWC',varid)
  DO L = 1, NBANDS
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = LZTWC_2d(land_idx,L)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = 1.e20
        END IF
      END DO
    END DO
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = L
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  status = NF_INQ_VARID(RESTART_NCID,'LZFSC',varid)
  DO L = 1, NBANDS
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = LZFSC_2d(land_idx,L)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = 1.e20
        END IF
      END DO
    END DO
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = L
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  status = NF_INQ_VARID(RESTART_NCID,'LZFPC',varid)
  DO L = 1, NBANDS
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = LZFPC_2d(land_idx,L)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = 1.e20
        END IF
      END DO
    END DO
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = L
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  status = NF_INQ_VARID(RESTART_NCID,'ADIMC',varid)
  DO L = 1, NBANDS
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = ADIMC_2d(land_idx,L)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = 1.e20
        END IF
      END DO
    END DO
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = L
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Snow Pack Depth
  status = NF_INQ_VARID(RESTART_NCID,'SNOWH',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = SNOWH(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Snow Pack Water Equivalent
  status = NF_INQ_VARID(RESTART_NCID,'SNEQV',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = SNEQV(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Snow Cover Extent
  status = NF_INQ_VARID(RESTART_NCID,'SNCOVR',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = SNCOVR(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Last Snow Counter
  status = NF_INQ_VARID(RESTART_NCID,'LSTSNW',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP_INT(J,I) = LSTSNW(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP_INT(J,I) = NODATA_INT
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_INT(RESTART_NCID,varid,start3d,count3d,TEMP_INT)
  END DO

  ! Write Exchange Coeff for heat/moisture
  status = NF_INQ_VARID(RESTART_NCID,'CH',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = CH(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Write Exchange Coeff for momentum
  status = NF_INQ_VARID(RESTART_NCID,'CM',varid)
  DO K = 1, NBANDS
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    land_idx = 1
    DO I = 1, ylen
      DO J = 1, xlen
        IF (LANDMASK(J,I) == 1) THEN
          TEMP(J,I) = CM(land_idx,K)
          land_idx = land_idx + 1
        ELSE
          TEMP(J,I) = NODATA
        END IF
      END DO
    END DO
    status = NF_PUT_VARA_REAL(RESTART_NCID,varid,start3d,count3d,TEMP)
  END DO

  ! Close file
  status = NF_CLOSE(RESTART_NCID)

END

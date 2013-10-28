SUBROUTINE READ_LSC_NETCDF()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS NETCDF-FORMAT LAND SURFACE CHARACTERISTICS

  ! Modifications:
  ! 2008-Oct-07 Included SAC parameters for use as unified model         Ben Livneh

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_LSC_NETCDF.f90,v 1.2 2007/10/04 20:44:20 vicadmin Exp $"/

  ! Define local variables
  INTEGER ndims,natts
  INTEGER varid
  CHARACTER*20 name
  INTEGER xtype, dimids(2)
  INTEGER J,I,K,NT,start2d(2),count2d(2),start3d(3),count3d(3),start,count
  REAL    TEMP(xlen,ylen)
  INTEGER TEMP_INT(xlen,ylen)
  INTEGER land_idx
  REAL    cumulative_depth(landlen)

  ! Get Col
  status = NF_INQ_VARID(LSC_NCID,'col',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,COL)

  ! Get Row
  status = NF_INQ_VARID(LSC_NCID,'row',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,ROW)

  ! Get Lon
  status = NF_INQ_VARID(LSC_NCID,'lon',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,LON)

  ! Get Lat
  status = NF_INQ_VARID(LSC_NCID,'lat',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,LAT)

  ! Get Cellid
  status = NF_INQ_VARID(LSC_NCID,'CellID',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,CELLID)

  ! Get NSOIL
  status = NF_INQ_VARID(LSC_NCID,'NSOIL',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,TEMP_INT)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        NSOIL(land_idx) = TEMP_INT(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Get SOILDEPTH
  status = NF_INQ_VARID(LSC_NCID,'SOILDEPTH',varid)
  cumulative_depth = 0
  DO K = 1, MAXNSOIL
    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = K
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1
    status = NF_GET_VARA_REAL(LSC_NCID,varid,start3d,count3d,TEMP)
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          SOILDEPTH(land_idx,K) = TEMP(J,I)
          cumulative_depth(land_idx) = cumulative_depth(land_idx) + SOILDEPTH(land_idx,K)
          SOILDEPTH_ACCUM(land_idx,K) = cumulative_depth(land_idx)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO
  END DO

  ! Get SOILTYP
  status = NF_INQ_VARID(LSC_NCID,'SOILTYP',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,TEMP_INT)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        SOILTYP(land_idx) = TEMP_INT(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Get SLOPETYP
  status = NF_INQ_VARID(LSC_NCID,'SLOPETYP',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,TEMP_INT)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        SLOPETYP(land_idx) = TEMP_INT(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Get TBOT
  status = NF_INQ_VARID(LSC_NCID,'TBOT',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        TBOT(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Get VEGTYP
  status = NF_INQ_VARID(LSC_NCID,'VEGTYP',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,TEMP_INT)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        VEGTYP(land_idx) = TEMP_INT(J,I)
        IF (VEGTYP(land_idx) .GT. 13) THEN
          VEGTYP(land_idx) = 13
        END IF
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Get ALBEDO, SHDFAC
  DO NT = 1,NMONTHS

    start3d(1) = 1
    start3d(2) = 1
    start3d(3) = NT
    count3d(1) = xlen
    count3d(2) = ylen
    count3d(3) = 1

    status = NF_INQ_VARID(LSC_NCID,'SHDFAC',varid)
    status = NF_GET_VARA_REAL(LSC_NCID,varid,start3d,count3d,TEMP)
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          SHDFAC(land_idx,NT)=TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

    status = NF_INQ_VARID(LSC_NCID,'ALBEDO',varid)
    status = NF_GET_VARA_REAL(LSC_NCID,varid,start3d,count3d,TEMP)
    land_idx = 1
    DO I = 1,ylen
      DO J = 1,xlen
        IF (LANDMASK(J,I) == 1) THEN
          ALBEDO(land_idx,NT)=TEMP(J,I)
          land_idx = land_idx + 1
        END IF
      END DO
    END DO

  ENDDO

  ! Get SNOALB
  status = NF_INQ_VARID(LSC_NCID,'SNOALB',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        SNOALB(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  ! Get SAC parameters

  status = NF_INQ_VARID(LSC_NCID,'UZTWM',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        UZTWM_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'UZFWM',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        UZFWM_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'UZK',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        UZK_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'PCTIM',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        PCTIM_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'ADIMP',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        ADIMP_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'RIVA',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        RIVA_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'ZPERC',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        ZPERC_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'REXP',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        REXP_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'LZTWM',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        LZTWM_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'LZFSM',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        LZFSM_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'LZFPM',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        LZFPM_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'LZSK',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        LZSK_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'LZPK',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        LZPK_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'PFREE',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        PFREE_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'SIDE',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        SIDE_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

  status = NF_INQ_VARID(LSC_NCID,'RSERV',varid)
  status = NF_GET_VAR_REAL(LSC_NCID,varid,TEMP)
  land_idx = 1
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        RSERV_2d(land_idx) = TEMP(J,I)
        land_idx = land_idx + 1
      END IF
    END DO
  END DO

END

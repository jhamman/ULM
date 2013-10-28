SUBROUTINE READ_IN_VAR(ncid, varid, step, data, ndims, LCOMP, LBAND)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS A SINGLE VARIABLE FOR 1 TIMESTEP FROM AN INPUT FILE

  ! INPUTS:
  ! ncid:  filehandle of the input netcdf file (obtained via ncopen function)
  ! varid: id number of the variable we're reading, within the netcdf file (obtained via nc_inq_var function)
  ! step:  index of the current timestep that's being read
  ! data:  array containing data to be read;
  !        inside this subroutine, the data array has the maximum dimensions, just to be able
  !        to store the data from any of the possible output variables.  The actual data being
  !        passed into this subroutine from outside might only use a subset of these dimensions,
  !        i.e. when calling this subroutine, you can pass it a variable with fewer dimensions.
  !        You can tell this subroutine which dimensions to use, via the ndims, LCOMP, and LBAND
  !        arguments (see below).
  ! ndims, LCOMP, LBAND:  determine the number of dimensions of the variable for input and output.
  !
  !        ndims: number of dimensions, with the horizontal dimension always counted as 2 dimensions
  !               (i.e. x and y) even though internally the variable is always stored with just 1
  !               horizontal dimension.  This may seem confusing, but it works for now.
  !
  !        LCOMP: logical flag, determining whether to write the output variable using 2 horizontal
  !               dimensions (x,y) or 1 horizontal dimension (land).  Using (x,y) is "uncompressed",
  !               while using (land) is "compressed".  Using (land) to index the grid cells allows
  !               us to skip cells that are outside the catchment boundaries (determined by the land
  !               mask), and can potentially save a substantial amount of space.
  !
  !        LBAND: logical flag, used in the case of 3d variables to determine whether the 3rd dimension
  !               is the soil layer or elevation band dimension.
  !
  !        The valid combinations of these arguments are shown below:
  !
  !        ndims  LCOMP  LBAND  assumed_input_dimensions*      assumed_output_dimensions*
  !        2      true   n/a    data(landlen)                  data(landlen)
  !        2      false  n/a    data(xlen,ylen)                data(landlen)
  !        3      true   false  data(landlen,MAXNSOIL)         data(landlen,MAXNSOIL)
  !        3      false  false  data(xlen,ylen,MAXNSOIL)       data(landlen,MAXNSOIL)
  !        3      true   true   data(landlen,NBANDS)           data(landlen,NBANDS)
  !        3      false  true   data(xlen,ylen,NBANDS)         data(landlen,NBANDS)
  !        4      true   n/a    data(landlen,NBANDS,MAXNSOIL)  data(landlen,NBANDS,MAXNSOIL)
  !        4      false  n/a    data(xlen,ylen,NBANDS,MAXNSOIL)data(landlen,NBANDS,MAXNSOIL)
  !
  !       *NOTE: the dimensions here don't include time; the time dimension should be present in the
  !        output file's definition of the variable we're writing.
  !

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_IN_VAR.f90,v 1.2 2007/10/04 20:44:20 vicadmin Exp $"/

  ! Define local variables
  INTEGER ncid
  INTEGER varid
  INTEGER step
  INTEGER ndims
  INTEGER I,J,K,L,land_idx
  INTEGER start2d(3),count2d(3),start2d_cmp(2),count2d_cmp(2)
  INTEGER start3d(4),count3d(4),start3d_cmp(3),count3d_cmp(3)
  INTEGER start4d(5),count4d(5),start4d_cmp(4),count4d_cmp(4)
  INTEGER TEMPDIM
  REAL temp(xlen,ylen)
  REAL temp_cmp(landlen)
  REAL data(landlen,(NBANDS+MAXNSOIL),MAXNSOIL)
  LOGICAL LCOMP
  LOGICAL LBAND

  IF (ndims == 2) THEN

    IF (LCOMP) THEN
      start2d_cmp(1) = 1
      start2d_cmp(2) = step
      count2d_cmp(1) = landlen
      count2d_cmp(2) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start2d_cmp,count2d_cmp,temp_cmp)
      data(:,1,1) = temp_cmp
    ELSE
      start2d(1) = 1
      start2d(2) = 1
      start2d(3) = step
      count2d(1) = xlen
      count2d(2) = ylen
      count2d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start2d,count2d,temp)
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            data(land_idx,1,1) = temp(J,I)
          END IF
        END DO
      END DO
    END IF

  ELSE IF (ndims == 3) THEN

    IF (LBAND == .FALSE.) THEN
      TEMPDIM = MAXNSOIL
    ELSE
      TEMPDIM = NBANDS
    END IF

    IF (LCOMP) THEN
      DO K = 1,TEMPDIM
        start3d_cmp(1) = 1
        start3d_cmp(2) = K
        start3d_cmp(3) = step
        count3d_cmp(1) = landlen
        count3d_cmp(2) = 1
        count3d_cmp(3) = 1
        status = NF_GET_VARA_REAL(ncid,varid,start3d_cmp,count3d_cmp,temp_cmp)
        data(:,K,1) = temp_cmp(:)
      END DO
    ELSE
      DO K = 1,TEMPDIM
        start3d(1) = 1
        start3d(2) = 1
        start3d(3) = K
        start3d(4) = step
        count3d(1) = xlen
        count3d(2) = ylen
        count3d(3) = 1
        count3d(4) = 1
        status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,temp)
        land_idx = 0
        DO I = 1,ylen
          DO J = 1,xlen
            IF (LANDMASK(J,I) == 1) THEN
              land_idx = land_idx + 1
              data(land_idx,K,1) = temp(J,I)
            END IF
          END DO
        END DO
      END DO
    END IF

  ELSE

    IF (LCOMP) THEN
      DO L = 1,MAXNSOIL
        DO K = 1,NBANDS
          start4d_cmp(1) = 1
          start4d_cmp(2) = K
          start4d_cmp(3) = L
          start4d_cmp(4) = step
          count4d_cmp(1) = landlen
          count4d_cmp(2) = 1
          count4d_cmp(3) = 1
          count4d_cmp(4) = 1
          status = NF_GET_VARA_REAL(ncid,varid,start4d_cmp,count4d_cmp,temp_cmp)
          data(:,K,L) = temp_cmp(:)
        END DO
      END DO
    ELSE
      DO L = 1,MAXNSOIL
        DO K = 1,NBANDS
          start4d(1) = 1
          start4d(2) = 1
          start4d(3) = K
          start4d(4) = L
          start4d(5) = step
          count4d(1) = xlen
          count4d(2) = ylen
          count4d(3) = 1
          count4d(4) = 1
          count4d(5) = 1
          status = NF_GET_VARA_REAL(ncid,varid,start4d,count4d,temp)
          land_idx = 0
          DO I = 1,ylen
            DO J = 1,xlen
              IF (LANDMASK(J,I) == 1) THEN
                land_idx = land_idx + 1
                data(land_idx,K,L) = temp(J,I)
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF

  END IF

END

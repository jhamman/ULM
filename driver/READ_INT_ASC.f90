SUBROUTINE READ_INT_ASC(filename,zlen,xlen,ylen,XLLCORNER,YLLCORNER,cellsize,TEMP)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS MODEL PARAMETER FILE (ascii arcmap format)

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_INT_ASC.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

  ! Define local variables
  CHARACTER*200 filename
  INTEGER zlen,xlen,ylen
  REAL    XLLCORNER,YLLCORNER,cellsize
  INTEGER TEMP(zlen,xlen,ylen),TEMP2(zlen,xlen,ylen)
  INTEGER K,J,I
  CHARACTER*15 tmp_str
  INTEGER tmp_xlen,tmp_ylen
  REAL    tmp_XLLCORNER,tmp_YLLCORNER,tmp_cellsize,tmp_NODATA_value

  OPEN(99,FILE=filename,STATUS='OLD')
  READ (99,*) tmp_str,tmp_xlen
  READ (99,*) tmp_str,tmp_ylen
  READ (99,*) tmp_str,tmp_XLLCORNER
  READ (99,*) tmp_str,tmp_YLLCORNER
  READ (99,*) tmp_str,tmp_cellsize
  READ (99,*) tmp_str,tmp_NODATA_value
  IF (tmp_xlen /= xlen) THEN
    WRITE(*,*) 'ERROR: NCOLS',tmp_xlen,'from',trim(filename),'does not match NCOLS',xlen,'from mask'
    STOP
  END IF
  IF (tmp_ylen /= ylen) THEN
    WRITE(*,*) 'ERROR: NROWS',tmp_ylen,'from',trim(filename),'does not match NROWS',ylen,'from mask'
    STOP
  END IF
  IF (tmp_XLLCORNER /= XLLCORNER) THEN
    WRITE(*,*) 'ERROR: XLLCORNER',tmp_XLLCORNER,'from',trim(filename),'does not match XLLCORNER',XLLCORNER,'from mask'
    STOP
  END IF
  IF (tmp_YLLCORNER /= YLLCORNER) THEN
    WRITE(*,*) 'ERROR: YLLCORNER',tmp_YLLCORNER,'from',trim(filename),'does not match YLLCORNER',YLLCORNER,'from mask'
    STOP
  END IF
  IF (tmp_cellsize /= cellsize) THEN
    WRITE(*,*) 'ERROR: cellsize',cellsize,'from',trim(filename),'does not match cellsize',cellsize,'from mask'
    STOP
  END IF
  READ (99,*) (((TEMP2(K,J,I),K=1,zlen),J=1,xlen),I=1,ylen)
  CLOSE(99)

  ! This assumes that ascii file presents the top row (i.e.
  ! the row with the maximum index) first, and proceeds in
  ! descending order.
  DO I=1,ylen
    DO J=1,xlen
      DO K=1,zlen
        TEMP(K,J,ylen+1-I) = TEMP2(K,J,I)
      END DO
    END DO
  END DO

END

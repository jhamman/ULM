      SUBROUTINE READ_REAL_SPATIAL(filename,nz,nx,ny,temp,LREAL8)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

      IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_REAL_SPATIAL.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

      INTEGER     nx, ny, nz, i, j, k
      REAL        temp(nz,nx,ny)
      REAL*8      temp8(nz,nx,ny)
      CHARACTER*200 filename
      INTEGER     ios
      LOGICAL     LREAL8
      
      OPEN(50,FILE = filename,IOSTAT=ios,STATUS = 'OLD', FORM='unformatted')
      IF (ios .NE. 0) THEN
         WRITE(*,*) 'File '//filename//' does not exist, Subroutine READ_REAL_SPATIAL'
         STOP
      END IF
      IF (LREAL8 .EQ. .TRUE.) THEN
        READ(50) (((temp8(i,j,k),i=1,nz),j=1,nx),k=1,ny)
        temp = temp8
      ELSE
        READ(50) (((temp(i,j,k),i=1,nz),j=1,nx),k=1,ny)
      END IF
      CLOSE(50)

!      WRITE(*,*) (((temp(i,j,k),i=1,1),j=1,15),k=1,15)

      END

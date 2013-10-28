      SUBROUTINE READ_INT_SPATIAL(filename,nz,nx,ny,temp)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

      IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_INT_SPATIAL.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

      INTEGER nx, ny, nz, i, j, k
      INTEGER temp(nz,nx,ny)
      CHARACTER*200 filename
      INTEGER ios
      
      OPEN(50,FILE = filename,IOSTAT=ios,STATUS = 'OLD', FORM='unformatted')
      IF (ios .NE. 0) THEN
         WRITE(*,*) 'File '//filename//' does not exist, Subroutine READ_INT_SPATIAL'
         STOP
      END IF
      READ(50) (((temp(i,j,k),i=1,nz),j=1,nx),k=1,ny)
      CLOSE(50)

!      WRITE(*,*) (((temp(i,j,k),i=1,1),j=1,15),k=1,15)

      END

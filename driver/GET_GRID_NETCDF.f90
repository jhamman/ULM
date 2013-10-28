SUBROUTINE GET_GRID_NETCDF()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! OPENS NETCDF-FORMAT LSC FILE AND GETS DATA DIMENSIONS AND LANDMASK

  ! Modifications
  ! 2008-May-15 Prints error message and quits if can't open input file.	TJB

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: GET_GRID_NETCDF.f90,v 1.3 2008/05/16 00:33:25 vicadmin Exp $"/

  ! Define local variables
  INTEGER ndims,nvars,ngatts,natts,unlimited
  INTEGER xdimid,ydimid,landid,zdimid,landcoverdimid,tstepdimid
  CHARACTER*20 name
  INTEGER xtype,dimids(2)
  INTEGER varid

  ! Open LSC file and get data dimensions
  ! Assume that LSC file is NOT compressed
  status = NF_OPEN(LSC, 0, LSC_NCID)
  IF (status .ne. NF_NOERR) THEN
    WRITE(*,*)'ERROR: cannot open param file', LSC
    STOP
  END IF
  status = NF_INQ(LSC_NCID, ndims, nvars, ngatts, unlimited)
  status = NF_INQ_DIMID(LSC_NCID,'z',zdimid)
  status = NF_INQ_DIMID(LSC_NCID,'y',ydimid)
  status = NF_INQ_DIMID(LSC_NCID,'x',xdimid)
  status = NF_INQ_DIMID(LSC_NCID,'t',tstepdimid)
  status = NF_INQ_DIMLEN(LSC_NCID,zdimid,MAXNSOIL)
  status = NF_INQ_DIMLEN(LSC_NCID,ydimid,ylen)
  status = NF_INQ_DIMLEN(LSC_NCID,xdimid,xlen)
  status = NF_INQ_DIMLEN(LSC_NCID,tstepdimid,NMONTHS)

  ! Get LANDMASK
  ALLOCATE (LANDMASK(xlen,ylen))
  status = NF_INQ_VARID(LSC_NCID,'land',varid)
  status = NF_GET_VAR_INT(LSC_NCID,varid,LANDMASK)

END

SUBROUTINE OPEN_RESTART()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! CREATES RESTART FILE

  ! Modifications:
  ! 2007-Nov-06 Added SnowSoilT (T12).						Ben Livneh
  ! 2007-Nov-12 Added LSTSNW.							TJB
  ! 2008-May-05 Removed T12, for compatibility with NOAH 2.8.			Ben Livneh
  ! 2008-May-15 Prints error message and quits if can't open output file.	TJB
  ! 2008-Jul-15 Added SnowSoilT (TPACK).					Ben Livneh
  ! 2008-Jul-24 Added PACH20.							Ben Livneh
  ! 2008-Oct-07 Included SAC parameters for use as unified model         Ben Livneh
  ! 2011-May-2 Added the _FillValue attribute to compliment missing_value for nco tools B.L.

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: OPEN_RESTART.f90,v 1.10 2008/07/30 22:33:03 vicadmin Exp $"/

  ! Define local variables
  INTEGER ndims,nvars,ngatts,natts,unlimited
  INTEGER xdimid,ydimid,leveldimid,banddimid,sacdimid,frzdimid
  INTEGER rowvarid,colvarid,latvarid,lonvarid,cellidvarid,landmaskvarid,varid
  CHARACTER*20 name
  INTEGER xtype, dims2d(2),dims3d(3),dims4d(4),dims4dsac(4),dims4dfrz(4)
  INTEGER I,J,NT

  ! Create RESTART file
  ! Data will not be compressed
  status = NF_CREATE(RESTFILE,0,RESTART_NCID)
  IF (status .ne. NF_NOERR) THEN
    WRITE(*,*)'ERROR: cannot open restart file',RESTFILE
    STOP
  END IF
  write(*,*)'maxnsoil restart',maxnsoil
  status = NF_DEF_DIM(RESTART_NCID,'x',xlen,xdimid)
  status = NF_DEF_DIM(RESTART_NCID,'y',ylen,ydimid)
  status = NF_DEF_DIM(RESTART_NCID,'z',MAXNSOIL,leveldimid)
  status = NF_DEF_DIM(RESTART_NCID,'band',NBANDS,banddimid)
  status = NF_DEF_DIM(RESTART_NCID,'zone',6,sacdimid)
  status = NF_DEF_DIM(RESTART_NCID,'--',10,frzdimid)
  dims2d(1) = xdimid
  dims2d(2) = ydimid
  dims3d(1) = xdimid
  dims3d(2) = ydimid
  dims3d(3) = banddimid
  dims4d(1) = xdimid
  dims4d(2) = ydimid
  dims4d(3) = banddimid
  dims4d(4) = leveldimid
  dims4dsac(1) = xdimid
  dims4dsac(2) = ydimid
  dims4dsac(3) = banddimid
  dims4dsac(4) = sacdimid
  dims4dfrz(1) = xdimid
  dims4dfrz(2) = ydimid
  dims4dfrz(3) = banddimid
  dims4dfrz(4) = frzdimid

  status = NF_PUT_ATT_TEXT(RESTART_NCID,NF_GLOBAL,'Conventions',7,'GDT 1.3')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,NF_GLOBAL,'file_name',31,RESTART)

  status = NF_DEF_VAR(RESTART_NCID,'col',NF_INT,2,dims2d,colvarid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,colvarid,'units',1,'-')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,colvarid,'long_name',11,'Grid Column')
  status = NF_PUT_ATT_INT(RESTART_NCID,colvarid,'missing_value',NF_INT,1, NODATA_INT)
  status = NF_PUT_ATT_INT(RESTART_NCID,colvarid,'_FillValue',NF_INT,1, NODATA_INT)

  status = NF_DEF_VAR(RESTART_NCID,'row',NF_INT,2,dims2d,rowvarid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,rowvarid,'units',1,'-')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,rowvarid,'long_name',8,'Grid Row')
  status = NF_PUT_ATT_INT(RESTART_NCID,rowvarid,'missing_value',NF_INT,1, NODATA_INT)
  status = NF_PUT_ATT_INT(RESTART_NCID,rowvarid,'_FillValue',NF_INT,1, NODATA_INT)

  status = NF_DEF_VAR(RESTART_NCID,'lon',NF_REAL,2,dims2d,lonvarid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,lonvarid,'units',12,'Degrees East')
  status = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'valid_min',NF_REAL,1,-180.)
  status = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'valid_max',NF_REAL,1,180.)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,lonvarid,'long_name',9,'Longitude')
  status = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'_FillValue',NF_REAL,1, NODATA)

  status = NF_DEF_VAR(RESTART_NCID,'lat',NF_REAL,2,dims2d,latvarid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,latvarid,'units',13,'Degrees North')
  status = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'valid_min',NF_REAL,1,-90.)
  status = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'valid_max',NF_REAL,1,90.)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,latvarid,'long_name',8,'Latitude')
  status = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'_FillValue',NF_REAL,1, NODATA)

  status = NF_DEF_VAR(RESTART_NCID,'land',NF_INT,2,dims2d,landmaskvarid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,landmaskvarid,'units',1,'-')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,landmaskvarid,'long_name',8,'Landmask')
  status = NF_PUT_ATT_INT(RESTART_NCID,landmaskvarid,'missing_value',NF_INT,1, NODATA_INT)
  status = NF_PUT_ATT_INT(RESTART_NCID,landmaskvarid,'_FillValue',NF_INT,1, NODATA_INT)

  status = NF_DEF_VAR(RESTART_NCID,'CellID',NF_INT,2,dims2d,cellidvarid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,cellidvarid,'units',1,'-')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,cellidvarid,'long_name',7,'Cell ID')
  status = NF_PUT_ATT_INT(RESTART_NCID,cellidvarid,'missing_value',NF_INT,1, NODATA_INT)
  status = NF_PUT_ATT_INT(RESTART_NCID,cellidvarid,'_FillValue',NF_INT,1, NODATA_INT)

  status = NF_DEF_VAR(RESTART_NCID,'SMC',NF_REAL,4,dims4d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',5,'m3/m3')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',10,'z band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',28,'Soil Layer Moisture Fraction')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',10,'z band y x')

  status = NF_DEF_VAR(RESTART_NCID,'SH2O',NF_REAL,4,dims4d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',5,'m3/m3')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',10,'z band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',29,'Liquid Soil Moisture Fraction')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',10,'z band y x')

  status = NF_DEF_VAR(RESTART_NCID,'STC',NF_REAL,4,dims4d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'K')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',10,'z band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',22,'Soil Layer Temperature')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',10,'z band y x')

  status = NF_DEF_VAR(RESTART_NCID,'SACST',NF_REAL,4,dims4dsac,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',13,'zone band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',18,'SAC Moisture state')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',13,'zone band y x')

  status = NF_DEF_VAR(RESTART_NCID,'FRZST',NF_REAL,4,dims4dfrz,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',8,'K and mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',11,'-- band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',27,'Soil temp and frozen states')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',11,'-- band y x')

  status = NF_DEF_VAR(RESTART_NCID,'T1',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'K')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Surface Temperature')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'TPACK',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'K')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',20,'Snowpack Temperature')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'PACH20',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',21,'Snowpack Liquid Water')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'CMC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',5,'kg/m2')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Canopy Interception')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')


  status = NF_DEF_VAR(RESTART_NCID,'UZTWC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',33,'Upper zone tension water contents')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'UZFWC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',30,'Upper zone free water contents')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'LZTWC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',33,'Lower zone tension water contents')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'LZFSC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',37,'Lower zone free supplemental contents')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'LZFPC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',32,'Lower zone free primary contents')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'ADIMC',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',40,'Tension water contents of the ADIMP area')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'SNOWH',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'m')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',10,'Snow Depth')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'SNEQV',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'m')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Snow Water Equivalent')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'SNCOVR',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'-')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',17,'Snow Cover Extent')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'LSTSNW',NF_INT,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'-')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',17,'Last snow counter')
  status = NF_PUT_ATT_INT(RESTART_NCID,varid,'missing_value',NF_INT,1, NODATA_INT)
  status = NF_PUT_ATT_INT(RESTART_NCID,varid,'_FillValue',NF_INT,1, NODATA_INT)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'CH',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',3,'m/s')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Heat Exchange Coeff')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_DEF_VAR(RESTART_NCID,'CM',NF_REAL,3,dims3d,varid)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',3,'m/s')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Momentum Exchange Coeff')
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
  status = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')

  status = NF_ENDDEF(RESTART_NCID)

  ! Write row, col, etc.

  status = NF_PUT_VAR_INT(RESTART_NCID,rowvarid,ROW)
  status = NF_PUT_VAR_INT(RESTART_NCID,colvarid,COL)
  status = NF_PUT_VAR_REAL(RESTART_NCID,lonvarid,LON)
  status = NF_PUT_VAR_REAL(RESTART_NCID,latvarid,LAT)
  status = NF_PUT_VAR_INT(RESTART_NCID,landmaskvarid,LANDMASK)
  status = NF_PUT_VAR_INT(RESTART_NCID,cellidvarid,CELLID)
!hack
!  status = NF_CLOSE(RESTART_NCID)
!  stop
END

SUBROUTINE READ_INITIAL()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS INITIAL CONDITIONS

  ! Modifications:
  ! 2007-Nov-06 Added T12.							Ben Livneh
  ! 2007-Nov-12 Added LSTSNW.							TJB
  ! 2007-Nov-15 Fixed initialization for zero-area bands.			TJB
  ! 2008-May-05 Removed T12, for compatibility with NOAH 2.8.                   Ben Livneh
  ! 2008-May-15 Prints error message and quits if can't open input file.        TJB
  ! 2008-Jul-15 Added TPACK.							Ben Livneh
  ! 2008-Jul-24 Added PACH20.							Ben Livneh
  ! 2008-Aug-12 Extended MAXSMC, etc to match REDPRM of SFLX 2.8.		Ben Livneh
  ! 2008-Aug-14 Added FLIMIT as constraint to initialize soil moisture 
  ! 2008-Oct-07 Included SAC parameters for use as unified model               Ben Livneh
  ! 2009-Feb-03 Created separate initialization for SAC extension              Ben Livneh
  ! 2009-Feb-15 Added Anderson PEMB snow model parameters for unified model    Ben Livneh

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_INITIAL.f90,v 1.10 2008/08/12 18:45:51 vicadmin Exp $"/

  ! Define local variables
  INTEGER ncid,ndims,nvars,ngatts,unlimited
  INTEGER xdimid,ydimid,zdimid,banddimid
  INTEGER xlen_initial,ylen_initial,MAXNSOIL_initial,NBANDS_initial
  INTEGER varid
  INTEGER J,I,K,L,start3d(3),count3d(3),start4d(4),count4d(4),NL,q
  INTEGER landmask_initial(xlen,ylen)
  INTEGER land_idx
  REAL    TEMP(xlen,ylen),SMCWLT1,SMLOW
  INTEGER TEMP_INT(xlen,ylen)
  REAL    MAXSMC(30)
  REAL    WLTSMC(20)
  REAL    REFSMC(20)
  REAL    MU(20)
  REAL    KSAT(20)
  REAL    SATPSI(20)
  REAL    BB(30)
  REAL    FLIMIT(30)
  REAL    SMCWLT(landlen)
  REAL    DZUPPER,DZLOWER
  INTEGER NUP,NLOW
  LOGICAL LAYERDEF
  REAL    PI,N1,DS,BETA1,DELT,SWLTSAND,SWLT,SFLD,SMAX
  PARAMETER(PI = 3.141592653)
  PARAMETER(N1 = 1.6)
  PARAMETER(DS = 0.0000025)
  PARAMETER(BETA1 = 16)
  PARAMETER(DELT = 86400)
  PARAMETER(SWLTSAND = 0.04)
  !  DATA MAXSMC/0.395, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.457, 0.464, 0.000, 0.200, 0.421, 0.457, 0.200, 0.395, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
  !DATA SATPSI/0.0350, 0.0363, 0.1413, 0.7586, 0.7586, 0.3548, 0.1349, 0.6166, 0.2630, 0.0977, 0.3236, 0.4677, 0.3548, 0.0000, 0.0350, 0.0363, 0.4677, 0.0350, 0.0350, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
  !DATA BB/4.05, Q4.26, 4.74, 5.33, 5.33, 5.25, 6.77, 8.72, 8.17, 10.73, 10.39, 11.55, 5.25, 0.00, 4.05, 4.26, 11.55, 4.05, 4.05,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
  ! Updated terms from Koren 2003 derivation -- only covers up to class 12 -- use sand for water
  DATA MAXSMC/0.37308, 0.38568, 0.41592, 0.46758, 0.47766, 0.43482, 0.41592, 0.4764, 0.44868, 0.42348, 0.48144, 0.46128, 0.464, 0.37308, 0.200, 0.421, 0.457, 0.200, 0.395, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
  DATA WLTSMC/0.03469064, 0.05199094, 0.08743051, 0.14637683, 0.10712489, 0.13941739, 0.15698002, 0.24386303, 0.21203782, 0.20755672, 0.28488226, 0.28290603, 0.069, 0.03469064, 0.012, 0.028, 0.135, 0.012, 0.023, 0.000/
  DATA REFSMC/0.15229854, 0.19015085, 0.26621888, 0.34851191, 0.34354197, 0.29455855, 0.28586507, 0.40984752, 0.35636059, 0.32561179, 0.43177247, 0.40382873, 0.319, 0.15229854, 0.116, 0.248, 0.389, 0.116, 0.196, 0.000/
  DATA FLIMIT/0.590, 0.900, 0.850, 0.740, 0.740, 0.740, 0.850, 0.800, 0.860, 0.860, 0.900, 0.900, 0.800, 0.900, 0.800, 0.900, 0.900, 0.800, 0.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
  DATA SATPSI/0.48094695,0.65051019,1.34286003,4.63205997,5.89793145,2.11235130,1.34286003,5.72247662,2.94468452,1.60962573,6.45723810,3.98286611,0.3548,0.48094695,0.0350,0.0363,0.4677,0.0350,0.0350,0.0000/
  DATA BB/3.387, 3.864, 4.500, 4.977, 3.705, 5.772, 7.203, 8.316, 8.316, 9.588, 10.383, 12.132, 5.25, 3.387, 4.05, 4.26, 11.55, 4.05, 4.05,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
  DATA MU/0.28513201, 0.23306613, 0.14960345, 0.10230300, 0.12465128, 0.13427207, 0.11844550, 0.03895119, 0.06705860, 0.07388107, 0.02396340, 0.03051445, 0.28513201, 0.23306613, 0.14960345, 0.10230300, 0.12465128, 0.13427207, 0.13427207, 0.000/
  DATA KSAT/1.7600E-4,1.4078E-5,5.2304E-6,2.8089E-6,2.8089E-6,3.3770E-6,4.4518E-6,2.0348E-6,2.4464E-6, 7.2199E-6,1.3444E-6,9.7394E-7,1.7600E-4,1.4078E-5,5.2304E-6,2.8089E-6,2.8089E-6,3.3770E-6,4.4518E-6, 0.000/
  ! Required for ULM -- LAYERDEF flag determines whether to (true) use input SAC storages
  ! and define a consistent set of physical layer depths; OR (false) to use input physical
  ! layer depths and to define a suitable set of SAC storages. 
  ! They key here is that max and min values are consistent (e.g. (UZTWM + UZFWM)/DZ*1000 + SMCWLT = SMCMAX)
  SMLOW = 0.5
  LAYERDEF = .TRUE.
  ! For either case, first define SAC parameters
  DO I = 1,landlen
!Fixed parameters
     SACPAR(I,1) = UZTWM_2d(I)
     SACPAR(I,2) = UZFWM_2d(I)
     SACPAR(I,3) = UZK_2d(I)
     SACPAR(I,4) = PCTIM_2d(I)
     SACPAR(I,5) = ADIMP_2d(I)
!     SACPAR(I,6) = RIVA_2d(I) Not considering riparian vegetation ET
     SACPAR(I,6) = 0.0
     SACPAR(I,7) = ZPERC_2d(I)
     SACPAR(I,8) = REXP_2d(I)
     SACPAR(I,9) = LZTWM_2d(I)  
     SACPAR(I,10) = LZFSM_2d(I)
     SACPAR(I,11) = LZFPM_2d(I)
     SACPAR(I,12) = LZSK_2d(I)
     SACPAR(I,13) = LZPK_2d(I)
     SACPAR(I,14) = PFREE_2d(I)
     SACPAR(I,15) = SIDE_2d(I)
     SACPAR(I,16) = RSERV_2d(I)
     ! Compute wilting point to be consistent with SFLX
     !     SMCWLT1 = MAXSMC(SOILTYP(I)) * (200.0/SATPSI(SOILTYP(I)))**(-1.0/BB(SOILTYP(I)))
     !     SMCWLT(I) = SMCWLT1 - SMLOW * SMCWLT1
     SMCWLT(I) = WLTSMC(SOILTYP(I))
     if (i==1) then
        write(*,*)'smcwlt(i) smlow smcmax'
        write(*,*)smcwlt(i),smlow,maxsmc(soiltyp(i))
     endif
! Frozen soil parameters 
! 1)stxt 2)tbot 3)rsmax 4)cksl 5)zbot 6)rtup 7)rtlow 8)psisat 9)swlt 10-13)zsoil
     FRZPAR(I,1) = SOILTYP(I)
     FRZPAR(I,2) = TBOT(I)
     FRZPAR(I,3) = 0.58
     FRZPAR(I,4) = 8.0
     FRZPAR(I,5) = SOILDEPTH_ACCUM(I,MAXNSOIL) * -1.0
     FRZPAR(I,6) = 1.0
     FRZPAR(I,7) = 1.0
     FRZPAR(I,8) = SATPSI(SOILTYP(I))
     FRZPAR(I,9) = SMCWLT(I)
     FRZPAR(I,10) = SOILDEPTH(I,1) * -1.0
     FRZPAR(I,11) = SOILDEPTH(I,2) * -1.0
     FRZPAR(I,12) =  SOILDEPTH(I,3) * -1.0
     FRZPAR(I,13) =  SOILDEPTH(I,4) * -1.0
  ENDDO

     ! No state case
  IF (INITIAL == '') THEN

     ! DEFAULT VALUES
     
     WRITE (*,*) 'No initial conditions file supplied; using defaults'
! Assign upper and lower zones to physical layers for soil moisture
     NUP=2
     NLOW=2
     WRITE (*,*) 'Setting # upper layers to',NUP,'Lower layers to',NLOW
     DO I = 1,landlen
        if(i==1)write(*,*)'soiltype',soiltyp(i),wltsmc(soiltyp(i))
        SWLT = WLTSMC(SOILTYP(I))
        SFLD = REFSMC(SOILTYP(I))
        SMAX = MAXSMC(SOILTYP(I))
        if(i==1)write(*,*)'swlt sfld smax',swlt,sfld,smax
        IF (LAYERDEF.EQ..TRUE..AND.MODEL_TYPE.EQ.1) THEN
! Define 4 soil layer depths
           DZUPPER = (UZTWM_2d(I)+UZFWM_2d(I))/((SMAX-SWLT)*1000)
           DZLOWER = (LZTWM_2d(I)+LZFPM_2d(I)+LZFSM_2d(I))/((SMAX-SWLT)*1000)
           if(i==1)write(*,*)dzupper,dzlower,swlt,smax,sfld
           SOILDEPTH(I,1) = DZUPPER * 0.25
           SOILDEPTH(I,2) = DZUPPER * 0.75
           SOILDEPTH(I,3) = DZLOWER * 0.375
           SOILDEPTH(I,4) = DZLOWER * 0.625
           if (SOILDEPTH(I,1) < 0) then
              write(*,*)'neg',UZTWM_2d(I),UZFWM_2d(I)
              write(*,*)'dzupper',DZUPPER,SMAX,SWLT
           endif
           FRZPAR(I,10) = SOILDEPTH(I,1) * -1.0
           FRZPAR(I,11) = SOILDEPTH(I,2) * -1.0
           FRZPAR(I,12) = SOILDEPTH(I,3) * -1.0
           FRZPAR(I,13) = SOILDEPTH(I,4) * -1.0
           if (i==1) then
              write(*,*)'Z1 Z2 Z3 Z4',SOILDEPTH(I,1),SOILDEPTH(I,2)
              write(*,*)SOILDEPTH(I,3),SOILDEPTH(I,4)
              write(*,*)uztwm_2d(i),uzfwm_2d(i)
              write(*,*)lztwm_2d(i),lzfpm_2d(i),lzfsm_2d(i)
              write(*,*)'tbot',tbot(i)
           endif
        ENDIF
        IF (LAYERDEF.EQ..FALSE..AND.MODEL_TYPE.EQ.1) THEN
! Define SAC storages based on layer depths
          UZTWM_2d(I) = (SFLD-SWLT)*(SOILDEPTH(I,1)+SOILDEPTH(I,2)) 
          UZFWM_2d(I) = (SMAX-SFLD)*(SOILDEPTH(I,1)+SOILDEPTH(I,2)) 
          UZK_2d(I) = 1 - (SFLD/SMAX)**N1
          LZTWM_2d(I) = (SFLD-SWLT)*(SOILDEPTH(I,3)+SOILDEPTH(I,4)) 
          LZFSM_2d(I) = (SMAX-SFLD)*(SOILDEPTH(I,3)+SOILDEPTH(I,4))*((SWLT/SMAX)**N1) 
          LZFSM_2d(I) = (SMAX-SFLD)*(SOILDEPTH(I,3)+SOILDEPTH(I,4))*(1 - (SWLT/SMAX)**N1) 
          LZSK_2d(I) = (1 - (SFLD/SMAX)**N1)/(1 + 2*(1-SWLT))
          LZPK_2d(I) = 1 - EXP(-1.0*((PI**2)*KSAT(SOILTYP(I))*1000*5.5*(DS**2)*(SOILDEPTH(I,3)+SOILDEPTH(I,4))*DELT)/MU(SOILTYP(I)))
          PFREE_2d(I) = (SWLT/SMAX)**N1
          ZPERC_2d(I) = ((LZTWM+LZFSM*(1-LZSK)) + (LZFPM*(1-LZPK))) / (LZFSM*LZSK+LZFPM*LZPK)
          REXP_2d(I) = (SWLT/(SWLTSAND - 0.001))**0.5
! Redefine storages
          SACPAR(I,1) = UZTWM_2d(I)
          SACPAR(I,2) = UZFWM_2d(I)
          SACPAR(I,3) = UZK_2d(I)
          SACPAR(I,4) = PCTIM_2d(I)
          SACPAR(I,5) = ADIMP_2d(I)
          !     SACPAR(I,6) = RIVA_2d(I) Not considering riparian vegetation ET
          SACPAR(I,6) = 0.0
          SACPAR(I,7) = ZPERC_2d(I)
          SACPAR(I,8) = REXP_2d(I)
          SACPAR(I,9) = LZTWM_2d(I)  
          SACPAR(I,10) = LZFSM_2d(I)
          SACPAR(I,11) = LZFPM_2d(I)
          SACPAR(I,12) = LZSK_2d(I)
          SACPAR(I,13) = LZPK_2d(I)
          SACPAR(I,14) = PFREE_2d(I)
          SACPAR(I,15) = SIDE_2d(I)
          SACPAR(I,16) = RSERV_2d(I)
       ENDIF
        DO J = 1, NBANDS
           IF (band_area(I,J) > 0) THEN
! Contents -- can add multipliers, these are arbitrary (0 <= m <= 1)   
              UZTWC_2d(I,J) = UZTWM_2d(I) * 0.75
              UZFWC_2d(I,J) = UZFWM_2d(I) * 0.50
              LZTWC_2d(I,J) = LZTWM_2d(I) * 0.75
              LZFSC_2d(I,J) = LZFSM_2d(I) * 0.25
              LZFPC_2d(I,J) = LZFPM_2d(I) * 0.50
              ADIMC_2d(I,J) = (UZTWM_2d(I) + LZTWM_2d(I))
              if(adimc_2d(i,j)<uztwc_2d(i,j)) adimc_2d(i,j)=uztwc_2d(i,j)
! Store in state array
              SACST(I,J,1) = UZTWC_2d(I,J)
              SACST(I,J,2) = UZFWC_2d(I,J)
              SACST(I,J,3) = LZTWC_2d(I,J)
              SACST(I,J,4) = LZFSC_2d(I,J)
              SACST(I,J,5) = LZFPC_2d(I,J)
              SACST(I,J,6) = ADIMC_2d(I,J)
! Need to rewrite frozen state allocation, once the indices are sorted out
! (i.e. the SAC-HT code is hard-coded to 5 physical layers
              FRZST(I,J,1) = TBOT(I)
              FRZST(I,J,2) = TBOT(I)
              FRZST(I,J,3) = TBOT(I)
              FRZST(I,J,4) = TBOT(I)
              FRZST(I,J,5) = 0
              FRZST(I,J,6) = SACST(I,J,1)
              FRZST(I,J,7) = SACST(I,J,2)
              FRZST(I,J,8) = SACST(I,J,3)
              FRZST(I,J,9) = SACST(I,J,4)
              FRZST(I,J,10)= SACST(I,J,5)
              IF(I==1) THEN
                 write(*,*)'swlt',smcwlt(I)
              endif
              DO K = 1,MAXNSOIL
! Divide SAC SMC's by 1000 (to get meters) and by respective layer 
! depths to keep as a voumetric frac, add SMCWLT as bottom boundary
                 IF (MODEL_TYPE == 1 .OR. MODEL_TYPE == 3) THEN
                    IF (K.LE.NUP) THEN
                       SMC(I,J,K) = ((SACST(I,J,1) + SACST(I,J,2))/1000)/(SOILDEPTH(I,1)+SOILDEPTH(I,2)) + SMCWLT(I)
                    ELSE
                       SMC(I,J,K) = ((SACST(I,J,3) + SACST(I,J,4) + SACST(I,J,5))/1000)/(SOILDEPTH(I,3)+SOILDEPTH(I,4)) + SMCWLT(I)
                    ENDIF
! **Important** Assume warm-season start-up, no frozen soil yet
                    SH2O(I,J,K) = SMC(I,J,K)
                    if(i==1) then
                       write(*,*)'smc',K,smc(i,j,k)
                    endif
                 ELSE
                    SMC(I,J,K)  = MAXSMC(SOILTYP(I)) * FLIMIT(SOILTYP(I))
                    SH2O(I,J,K) = SMC(I,J,K)
                 END IF
                 STC(I,J,K)  = TBOT(I)
              END DO
              if (i==1) then
                 do q = 1,5
                    write(*,*)'sacst',q,sacst(i,j,q)
                 enddo
              endif
              CMC(I,J)    = 0.0
              SNOWH(I,J)  = 0.0
              SNEQV(I,J)  = 0.0
              SNCOVR(I,J) = 0.0
              LSTSNW(I,J) = 0
              CH(I,J)     = 1.0e-4
              CM(I,J)     = 1.0e-4
              T1(I,J)     = TBOT(I)
              TPACK(I,J)  = 273.15
              PACH20(I,J) = 0.0
              DO NL = 1,MAXNSNOW             
                 DSNOW(I,J,NL)    = 0.0
                 PSNOW(I,J,NL)    = 0.0
                 RTTSNOW(I,J,NL)  = 0.0
                 RTTDTSNOW(I,J,NL)= 0.0
                 WTSNOW(I,J,NL)   = 0.0
                 WTDTSNOW(I,J,NL) = 0.0
                 TTSNOW(I,J,NL)   = 0.0
                 TTDTSNOW(I,J,NL) = 0.0          
              END DO
           ELSE
              DO K = 1,MAXNSOIL
                 SMC(I,J,K)  = NODATA
                 SH2O(I,J,K) = NODATA
                 STC(I,J,K)  = NODATA
              END DO
              UZTWC_2d(I,J) = NODATA
              UZFWC_2d(I,J) = NODATA
              LZTWC_2d(I,J) = NODATA
              LZFPC_2d(I,J) = NODATA
              LZFSC_2d(I,J) = NODATA
              ADIMC_2d(I,J) = NODATA
              CMC(I,J)    = NODATA
              SNOWH(I,J)  = NODATA
              SNEQV(I,J)  = NODATA
              SNCOVR(I,J) = NODATA
              LSTSNW(I,J) = NODATA_INT
              CH(I,J)     = NODATA
              CM(I,J)     = NODATA
              T1(I,J)     = NODATA
              TPACK(I,J)  = NODATA
              PACH20(I,J) = NODATA
              DO NL = 1,MAXNSNOW             
                 DSNOW(I,J,NL)    = NODATA
                 PSNOW(I,J,NL)    = NODATA
                 RTTSNOW(I,J,NL)  = NODATA
                 RTTDTSNOW(I,J,NL)= NODATA
                 WTSNOW(I,J,NL)   = NODATA
                 WTDTSNOW(I,J,NL) = NODATA
                 TTSNOW(I,J,NL)   = NODATA
                 TTDTSNOW(I,J,NL) = NODATA          
              END DO
           END IF
        END DO
     END DO
  ELSE

    ! SPECIFIED VALUES

    WRITE (*,*) 'Reading initial conditions from ',trim(INITIAL)

    status = NF_OPEN(INITIAL, 0, ncid)
    IF (status .ne. NF_NOERR) THEN
      WRITE(*,*)'ERROR: cannot open initial condition file',INITIAL
      STOP
    END IF
    status = NF_INQ(ncid, ndims, nvars, ngatts, unlimited)
    status = NF_INQ_DIMID(ncid,'z',zdimid)
    status = NF_INQ_DIMID(ncid,'band',banddimid)
    status = NF_INQ_DIMID(ncid,'y',ydimid)
    status = NF_INQ_DIMID(ncid,'x',xdimid)
    status = NF_INQ_DIMLEN(ncid,zdimid,MAXNSOIL_initial)
    status = NF_INQ_DIMLEN(ncid,banddimid,NBANDS_initial)
    status = NF_INQ_DIMLEN(ncid,ydimid,ylen_initial)
    status = NF_INQ_DIMLEN(ncid,xdimid,xlen_initial)

write(*,*)'MAXNSOIL_initial',MAXNSOIL_initial
write(*,*)'NBANDS_initial',NBANDS_initial
write(*,*)'ylen_initial',ylen_initial
write(*,*)'xlen_initial',xlen_initial

    ! Validate data dimensions
    IF (MAXNSOIL_initial /= MAXNSOIL) THEN
      WRITE(*,*)'ERROR: MAXNSOIL',MAXNSOIL, &
        ' (from lsc file ',TRIM(LSC),') /= MAXNSOIL ',MAXNSOIL_initial, &
        ' (from initial condition file ',TRIM(INITIAL),')'
      STOP
    END IF
    IF (NBANDS_initial /= NBANDS) THEN
      WRITE(*,*)'ERROR: NBANDS',NBANDS, &
        ' (from lsc file ',TRIM(LSC),') /= NBANDS ',NBANDS_initial, &
        ' (from initial condition file ',TRIM(INITIAL),')'
      STOP
    END IF
    IF (ylen_initial /= ylen) THEN
      WRITE(*,*)'ERROR: ylen',ylen, &
        ' (from lsc file ',TRIM(LSC),') /= ylen ',ylen_initial, &
        ' (from initial condition file ',TRIM(INITIAL),')'
      STOP
    END IF
    IF (xlen_initial /= xlen) THEN
      WRITE(*,*)'ERROR: xlen',xlen, &
        ' (from lsc file ',TRIM(LSC),') /= xlen ',xlen_initial, &
        ' (from initial condition file ',TRIM(INITIAL),')'
      STOP
    END IF
    
    ! Get landmask and compare to landmask from LSC file
    status = NF_INQ_VARID(ncid,'land',varid)
    status = NF_GET_VAR_INT(ncid,varid,landmask_initial)
    DO I=1,ylen
      DO J=1,xlen
        IF (landmask_initial(J,I) /= LANDMASK(J,I)) THEN
          WRITE(*,*)'ERROR: landmask from lsc file ',TRIM(LSC), &
          ' /= landmask from initial condition file ',TRIM(INITIAL)
          STOP
        END IF
      END DO
    END DO

    ! Get Soil Moisture Content
    status = NF_INQ_VARID(ncid,'SMC',varid)
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
        status = NF_GET_VARA_REAL(ncid,varid,start4d,count4d,TEMP)
        land_idx = 1
        DO I = 1, ylen
          DO J = 1, xlen
            IF (LANDMASK(J,I) == 1) THEN
              SMC(land_idx,K,L) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    END DO

    ! Get Soil Liquid Moisture Content
    status = NF_INQ_VARID(ncid,'SH2O',varid)
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
        status = NF_GET_VARA_REAL(ncid,varid,start4d,count4d,TEMP)
        land_idx = 1
        DO I = 1, ylen
          DO J = 1, xlen
            IF (LANDMASK(J,I) == 1) THEN
              SH2O(land_idx,K,L) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    END DO

    ! Get Soil Temperature
    status = NF_INQ_VARID(ncid,'STC',varid)
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
        status = NF_GET_VARA_REAL(ncid,varid,start4d,count4d,TEMP)
        land_idx = 1
        DO I = 1, ylen
          DO J = 1, xlen
            IF (LANDMASK(J,I) == 1) THEN
              STC(land_idx,K,L) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    END DO

    ! Get SAC states
    status = NF_INQ_VARID(ncid,'SACST',varid)
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
        status = NF_GET_VARA_REAL(ncid,varid,start4d,count4d,TEMP)
        land_idx = 1
        DO I = 1, ylen
          DO J = 1, xlen
            IF (LANDMASK(J,I) == 1) THEN
              SACST(land_idx,K,L) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    END DO

    ! Get SAC Frozen states
    status = NF_INQ_VARID(ncid,'FRZST',varid)
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
        status = NF_GET_VARA_REAL(ncid,varid,start4d,count4d,TEMP)
        land_idx = 1
        DO I = 1, ylen
          DO J = 1, xlen
            IF (LANDMASK(J,I) == 1) THEN
              FRZST(land_idx,K,L) = TEMP(J,I)
              land_idx = land_idx + 1
            END IF
          END DO
        END DO
      END DO
    END DO

    ! Get Skin Temperature
    status = NF_INQ_VARID(ncid,'T1',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            T1(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Snow Pack Temperature
    status = NF_INQ_VARID(ncid,'TPACK',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            TPACK(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Snow Pack Liquid Water Content
    status = NF_INQ_VARID(ncid,'PACH20',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            PACH20(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Canopy Moisture
    status = NF_INQ_VARID(ncid,'CMC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            CMC(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get SAC Parameters
    status = NF_INQ_VARID(ncid,'UZTWC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            UZTWC_2d(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    status = NF_INQ_VARID(ncid,'UZFWC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            UZFWC_2d(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    status = NF_INQ_VARID(ncid,'LZTWC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            LZTWC_2d(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    status = NF_INQ_VARID(ncid,'LZFSC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            LZFSC_2d(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    status = NF_INQ_VARID(ncid,'LZFPC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            LZFPC_2d(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    status = NF_INQ_VARID(ncid,'ADIMC',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            ADIMC_2d(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Snow Pack Depth
    status = NF_INQ_VARID(ncid,'SNOWH',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            SNOWH(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Snow Pack Water Equivalent
    status = NF_INQ_VARID(ncid,'SNEQV',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            SNEQV(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Snow Cover Extent
    status = NF_INQ_VARID(ncid,'SNCOVR',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            SNCOVR(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Get Last Snow Counter
    status = NF_INQ_VARID(ncid,'LSTSNW',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_INT(ncid,varid,start3d,count3d,TEMP_INT)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            LSTSNW(land_idx,K) = TEMP_INT(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

    ! Exchange coefficients
    status = NF_INQ_VARID(ncid,'CH',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            CH(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO
    status = NF_INQ_VARID(ncid,'CM',varid)
    DO K = 1, NBANDS
      start3d(1) = 1
      start3d(2) = 1
      start3d(3) = K
      count3d(1) = xlen
      count3d(2) = ylen
      count3d(3) = 1
      status = NF_GET_VARA_REAL(ncid,varid,start3d,count3d,TEMP)
      land_idx = 1
      DO I = 1, ylen
        DO J = 1, xlen
          IF (LANDMASK(J,I) == 1) THEN
            CM(land_idx,K) = TEMP(J,I)
            land_idx = land_idx + 1
          END IF
        END DO
      END DO
    END DO

! Align soil depths with SAC-storages or vice versa -- needed for soil water balance
     NUP=2
     NLOW=2
     WRITE (*,*) 'Setting # upper layers to',NUP,'Lower layers to',NLOW
     DO I = 1,landlen
        if(i==1)write(*,*)'soiltype',soiltyp(i),wltsmc(soiltyp(i))
        SWLT = WLTSMC(SOILTYP(I))
        SFLD = REFSMC(SOILTYP(I))
        SMAX = MAXSMC(SOILTYP(I))
        if(i==1)write(*,*)'swlt sfld smax',swlt,sfld,smax
        IF (LAYERDEF.EQ..TRUE..AND.MODEL_TYPE.EQ.1) THEN
! Define 4 soil layer depths
           DZUPPER = (UZTWM_2d(I)+UZFWM_2d(I))/((SMAX-SWLT)*1000)
           DZLOWER = (LZTWM_2d(I)+LZFPM_2d(I)+LZFSM_2d(I))/((SMAX-SWLT)*1000)
           if(i==1)write(*,*)dzupper,dzlower,swlt,smax,sfld
           SOILDEPTH(I,1) = DZUPPER * 0.25
           SOILDEPTH(I,2) = DZUPPER * 0.75
           SOILDEPTH(I,3) = DZLOWER * 0.375
           SOILDEPTH(I,4) = DZLOWER * 0.625
           if (SOILDEPTH(I,1) < 0) then
              write(*,*)'neg',UZTWM_2d(I),UZFWM_2d(I)
              write(*,*)'dzupper',DZUPPER,SMAX,SWLT
           endif
           FRZPAR(I,10) = SOILDEPTH(I,1) * -1.0
           FRZPAR(I,11) = SOILDEPTH(I,2) * -1.0
           FRZPAR(I,12) = SOILDEPTH(I,3) * -1.0
           FRZPAR(I,13) = SOILDEPTH(I,4) * -1.0
           if (i==1) then
              write(*,*)'Z1 Z2 Z3 Z4',SOILDEPTH(I,1),SOILDEPTH(I,2)
              write(*,*)SOILDEPTH(I,3),SOILDEPTH(I,4)
              write(*,*)uztwm_2d(i),uzfwm_2d(i)
              write(*,*)lztwm_2d(i),lzfpm_2d(i),lzfsm_2d(i)
              write(*,*)'tbot',tbot(i)
           endif
        ENDIF
        ! define SAC-storages based on soil depths -- rare case
        IF (LAYERDEF.EQ..FALSE..AND.MODEL_TYPE.EQ.1) THEN
! Define SAC storages based on layer depths
          UZTWM_2d(I) = (SFLD-SWLT)*(SOILDEPTH(I,1)+SOILDEPTH(I,2)) 
          UZFWM_2d(I) = (SMAX-SFLD)*(SOILDEPTH(I,1)+SOILDEPTH(I,2)) 
          UZK_2d(I) = 1 - (SFLD/SMAX)**N1
          LZTWM_2d(I) = (SFLD-SWLT)*(SOILDEPTH(I,3)+SOILDEPTH(I,4)) 
          LZFSM_2d(I) = (SMAX-SFLD)*(SOILDEPTH(I,3)+SOILDEPTH(I,4))*((SWLT/SMAX)**N1) 
          LZFSM_2d(I) = (SMAX-SFLD)*(SOILDEPTH(I,3)+SOILDEPTH(I,4))*(1 - (SWLT/SMAX)**N1) 
          LZSK_2d(I) = (1 - (SFLD/SMAX)**N1)/(1 + 2*(1-SWLT))
          LZPK_2d(I) = 1 - EXP(-1.0*((PI**2)*KSAT(SOILTYP(I))*1000*5.5*(DS**2)*(SOILDEPTH(I,3)+SOILDEPTH(I,4))*DELT)/MU(SOILTYP(I)))
          PFREE_2d(I) = (SWLT/SMAX)**N1
          ZPERC_2d(I) = ((LZTWM+LZFSM*(1-LZSK)) + (LZFPM*(1-LZPK))) / (LZFSM*LZSK+LZFPM*LZPK)
          REXP_2d(I) = (SWLT/(SWLTSAND - 0.001))**0.5
! Redefine storages
          SACPAR(I,1) = UZTWM_2d(I)
          SACPAR(I,2) = UZFWM_2d(I)
          SACPAR(I,3) = UZK_2d(I)
          SACPAR(I,4) = PCTIM_2d(I)
          SACPAR(I,5) = ADIMP_2d(I)
          !     SACPAR(I,6) = RIVA_2d(I) Not considering riparian vegetation ET
          SACPAR(I,6) = 0.0
          SACPAR(I,7) = ZPERC_2d(I)
          SACPAR(I,8) = REXP_2d(I)
          SACPAR(I,9) = LZTWM_2d(I)  
          SACPAR(I,10) = LZFSM_2d(I)
          SACPAR(I,11) = LZFPM_2d(I)
          SACPAR(I,12) = LZSK_2d(I)
          SACPAR(I,13) = LZPK_2d(I)
          SACPAR(I,14) = PFREE_2d(I)
          SACPAR(I,15) = SIDE_2d(I)
          SACPAR(I,16) = RSERV_2d(I)
       ENDIF
        DO J = 1, NBANDS
           IF (band_area(I,J) > 0) THEN
              if(adimc_2d(i,j)<uztwc_2d(i,j)) adimc_2d(i,j)=uztwc_2d(i,j)
! Store in state array
              SACST(I,J,1) = UZTWC_2d(I,J)
              SACST(I,J,2) = UZFWC_2d(I,J)
              SACST(I,J,3) = LZTWC_2d(I,J)
              SACST(I,J,4) = LZFSC_2d(I,J)
              SACST(I,J,5) = LZFPC_2d(I,J)
              SACST(I,J,6) = ADIMC_2d(I,J)
! Need to rewrite frozen state allocation, once the indices are sorted out
! (i.e. the SAC-HT code is hard-coded to 5 physical layers
              FRZST(I,J,1) = TBOT(I)
              FRZST(I,J,2) = TBOT(I)
              FRZST(I,J,3) = TBOT(I)
              FRZST(I,J,4) = TBOT(I)
              FRZST(I,J,5) = 0
              FRZST(I,J,6) = SACST(I,J,1)
              FRZST(I,J,7) = SACST(I,J,2)
              FRZST(I,J,8) = SACST(I,J,3)
              FRZST(I,J,9) = SACST(I,J,4)
              FRZST(I,J,10)= SACST(I,J,5)
              IF(I==1) THEN
                 write(*,*)'swlt',smcwlt(I)
              endif
              DO K = 1,MAXNSOIL
! Divide SAC SMC's by 1000 (to get meters) and by respective layer 
! depths to keep as a voumetric frac, add SMCWLT as bottom boundary
                 IF (MODEL_TYPE == 1 .OR. MODEL_TYPE == 3) THEN
                    IF (K.LE.NUP) THEN
                       SMC(I,J,K) = ((SACST(I,J,1) + SACST(I,J,2))/1000)/(SOILDEPTH(I,1)+SOILDEPTH(I,2)) + SMCWLT(I)
                    ELSE
                       SMC(I,J,K) = ((SACST(I,J,3) + SACST(I,J,4) + SACST(I,J,5))/1000)/(SOILDEPTH(I,3)+SOILDEPTH(I,4)) + SMCWLT(I)
                    ENDIF
! **Important** Assume warm-season start-up, no frozen soil yet
                    SH2O(I,J,K) = SMC(I,J,K)
                    if(i==1) then
                       write(*,*)'smc',K,smc(i,j,k)
                    endif
                 ELSE
                    SMC(I,J,K)  = MAXSMC(SOILTYP(I)) * FLIMIT(SOILTYP(I))
                    SH2O(I,J,K) = SMC(I,J,K)
                 END IF
                 STC(I,J,K)  = TBOT(I)
              END DO
              if (i==1) then
                 do q = 1,5
                    write(*,*)'sacst',q,sacst(i,j,q)
                 enddo
              endif
              DO NL = 1,MAXNSNOW             
                 DSNOW(I,J,NL)    = 0.0
                 PSNOW(I,J,NL)    = 0.0
                 RTTSNOW(I,J,NL)  = 0.0
                 RTTDTSNOW(I,J,NL)= 0.0
                 WTSNOW(I,J,NL)   = 0.0
                 WTDTSNOW(I,J,NL) = 0.0
                 TTSNOW(I,J,NL)   = 0.0
                 TTDTSNOW(I,J,NL) = 0.0          
              END DO
           ELSE
              DO K = 1,MAXNSOIL
                 SMC(I,J,K)  = NODATA
                 SH2O(I,J,K) = NODATA
                 STC(I,J,K)  = NODATA
              END DO
              UZTWC_2d(I,J) = NODATA
              UZFWC_2d(I,J) = NODATA
              LZTWC_2d(I,J) = NODATA
              LZFPC_2d(I,J) = NODATA
              LZFSC_2d(I,J) = NODATA
              ADIMC_2d(I,J) = NODATA
              CMC(I,J)    = NODATA
              SNOWH(I,J)  = NODATA
              SNEQV(I,J)  = NODATA
              SNCOVR(I,J) = NODATA
              LSTSNW(I,J) = NODATA_INT
              CH(I,J)     = NODATA
              CM(I,J)     = NODATA
              T1(I,J)     = NODATA
              TPACK(I,J)  = NODATA
              PACH20(I,J) = NODATA
              DO NL = 1,MAXNSNOW             
                 DSNOW(I,J,NL)    = NODATA
                 PSNOW(I,J,NL)    = NODATA
                 RTTSNOW(I,J,NL)  = NODATA
                 RTTDTSNOW(I,J,NL)= NODATA
                 WTSNOW(I,J,NL)   = NODATA
                 WTDTSNOW(I,J,NL) = NODATA
                 TTSNOW(I,J,NL)   = NODATA
                 TTDTSNOW(I,J,NL) = NODATA          
              END DO
           END IF
        END DO
     END DO

 END IF

END

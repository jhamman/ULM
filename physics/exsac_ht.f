      SUBROUTINE EXSAC_HT(NSOLD,DTM,PCP,TMP,ETP,LAND_ID,iBAND,IPR,
C     SAC PARAMETERS
     &     UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
     &     REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
     &     SIDE,RSERV,EFC,
C     SAC State variables  
     &     UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,
C     FROZEN GROUND VARIABLES
     &     SNEQV,COVER,SNOWH,FRZST,FRZPAR,
     &     NSOIL,NUPL,NSAC,FRZDUP,FRZDBT,
     &     FROST,TSINT,SWINT,SWHINT,
     &     SACST_PRV,STEP1,
C     SAC OUTPUTS
     &     SMC,SH2O,QS,QG,Q,TET)
      
C      IMPLICIT NONE

C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: ex_sac1.f,v 1.1 2006/09/01 21:59:43 vicadmin Exp $"/


C     ...THIS SUBROUTINE IS THE EXECUTION ROUTINE FOR SMFLX MODEL...
C

      INTEGER NSOLD
      REAL    DTM
      REAL    PCP
      REAL    TMP
      REAL    ETP
      REAL    QS
      REAL    QG
      REAL    Q
      REAL    ETA
      REAL    TA
      REAL    SNEQV,COVER,AESC,SNOWH,SH,STXT,TBOT,RSMAX,ZBOT
      REAL    RTUP,RTLW,PSISAT,SWLT,ZO,Z1,Z2,Z3,Z4,SUPM,SLWM
      REAL    SMAX,BRT,QUARTZ,STYPE,ZSOIL_PRE(5),ZSOIL_POST(5)
      REAL    LZTWM,LZFSM,LZFPM,UZTWM,UZFWM,UZK,LZPK,LZSK
      REAL    LZTWH,LZFSH,LZFPH,UZTWH,UZFWH,ADIMP
      REAL    LZTWC,LZFSC,LZFPC,UZTWC,UZFWC,ADIMC
      REAL    TOTAL_S1, TOTAL_S2
      REAL    DT
      REAL    DS
      REAL    SACPAR(17),SACST(6),SACST_PRV(6),FRZPAR(14),FRZST(10)
      REAL    SMC(5),SH2O(5),DSINT(5),DSINTW(5),DEPTH_ARRAY(5)
      REAL    FROST,TSINT(5),SWINT(5),SWHINT(5)
      LOGICAL STEP1

      COMMON/FSMCO1/FGCO(6),RSUM(7),PPE,PSC,PTA,PWE
      COMMON/FSUMS1/SROT,SIMPVT,SRODT,SROST,SINTFT,SGWFP,SGWFS,SRECHT,
     &              SETT,SE1,SE3,SE4,SE5

C    TURN OFF FROZEN GROUND PROCESS

      IFRZE = 0

C     COMPUTE TOTAL INITIAL STORAGE

      TOTAL_S1 = UZTWC + UZFWC + LZTWC + LZFSC + LZFPC + ADIMC

C     COMPUTE SURFACE MOISTURE FLUXES

      EDMND = ETP
      PXV = PCP

C     COMPUTE OTHER NECESSARY QUANTITIES FOR SAC_HT
      TA = TMP
      WE = SNEQV
      AESC = COVER
      SH = SNOWH*100 ! Convert SNOWH into cm for hrt1.f
      DT = DTM/86400.0
      IVERS = 1
      DTFRZ = 1800     ! frozen ground time step
      IDTFRZ = 2      ! number of frz-grnd step per timestep
C     OK Mesonet depths
      DEPTH_ARRAY = (/0.05,0.25,0.60,0.75,1.00/)
      DO I = 1, 5
         DSINT(I) = DEPTH_ARRAY(I)
         DSINTW(I) = DEPTH_ARRAY(I)
      ENDDO
      NDSINT = SIZE(DEPTH_ARRAY)
      NDINTW = SIZE(DEPTH_ARRAY)
      NORMALIZE = 1.0   ! Flag to normalize temp and moisture profile

C Frost depth parameters FRZDUP,FRZDBT,FROST computed in FLAND 
C Option to normalize soil profile to desired profile provided
C by NORMALIZE variable, embodied in SWINT, SWHINT, DSINT, NDSINT,
C DSINTW, and NDINTW variables.   



C     ORGANIZE PARAMETERS AND DATA INTO ARRAYS FOR PASSAGE TO FLAND1

C     SAC PARAMETERS
      SACPAR(1) = UZTWM
      SACPAR(2) = UZFWM
      SACPAR(3) = UZK
      SACPAR(4) = PCTIM
      SACPAR(5) = ADIMP    
      SACPAR(6) = RIVA
      SACPAR(7) = ZPERC
      SACPAR(8) = REXP
      SACPAR(9) = LZTWM    
      SACPAR(10) = LZFSM
      SACPAR(11) = LZFPM
      SACPAR(12) = LZSK
      SACPAR(13) = LZPK
      SACPAR(14) = PFREE
      SACPAR(15) = SIDE
      SACPAR(16) = RSERV
      SACPAR(17) = EFC

C     SAC states, current and from previous time-step   
      SACST(1) = UZTWC
      SACST(2) = UZFWC
      SACST(3) = LZTWC    
      SACST(4) = LZFSC
      SACST(5) = LZFPC
      SACST(6) = ADIMC 

C     SAC frozen parameters and states passed from driver code

C      DO I = 1,10
C         if (ipr == 1) then
C            if (i.gt.5) then
C               write(*,*)'prefland',I,FRZST(I),SACST(I-5)
C            else
C               write(*,*)'prefland',I,FRZST(I)
C            endif
C         endif
C      ENDDO

      if (ipr ==1) write(*,*)'ABOVE CALL TO FLAND1'
      if (ipr ==1) write(*,*)'PXV',PXV
      if (ipr ==1) write(*,*)'EDMND',EDMND
      if (ipr ==1) write(*,*)'TA',TA
      if (ipr ==1) write(*,*)'WE',WE
      if (ipr ==1) write(*,*)'AESC',AESC
      if (ipr ==1) write(*,*)'SH',SH
      if (ipr ==1) write(*,*)'DT',DT
      if (ipr ==1) then
         DO I = 1, 6
            write(*,*)'SACST:',i,SACST(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, 10
            write(*,*)'FRZST:',i,FRZST(i)
         ENDDO
      endif      
      if (ipr ==1) then
         DO I = 1, 17
            write(*,*)'SACPAR:',i,SACPAR(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, 14
            write(*,*)'FRZPAR:',i,FRZPAR(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'NSOIL',NSOIL
      if (ipr ==1) write(*,*)'NUPL',NUPL
      if (ipr ==1) write(*,*)'NSAC',NSAC
      if (ipr ==1) write(*,*)'IVERS',IVERS
      if (ipr ==1) write(*,*)'SURF',SURF
      if (ipr ==1) write(*,*)'GRND',GRND
      if (ipr ==1) write(*,*)'TET',TET
      if (ipr ==1) then
         DO I = 1, NSOIL
            write(*,*)'SMC:',i,SMC(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, NSOIL
            write(*,*)'SH2O:',i,SH2O(i)
         ENDDO
      endif      
      if (ipr ==1) then
         DO I = 1, 6
            write(*,*)'SACST_PRV:',i,SACST_PRV(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'DTFRZ',DTFRZ
      if (ipr ==1) write(*,*)'IDTFRZ',IDTFRZ
      if (ipr ==1) write(*,*)'FRZDUP',FRZDUP
      if (ipr ==1) write(*,*)'FRZDBT',FRZDBT
      if (ipr ==1) write(*,*)'FROST',FROST
      if (ipr ==1) then
         DO I = 1, NDSINT
            write(*,*)'TSINT:',i,TSINT(i)
         ENDDO
      endif
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'SWINT:',i,SWINT(i)
         ENDDO
      endif
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'SWHINT:',i,SWHINT(i)
         ENDDO
      endif
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'SWHINT:',i,SWHINT(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, NDSINT
            write(*,*)'DSINT:',i,DSINT(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'NDSINT',NDSINT
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'DSINTW:',i,DSINTW(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'NDINTW',NDINTW
      if (ipr ==1) write(*,*)'NORMALIZE',NORMALIZE


C SAC frozen states FRZST(*) are passed directly to FLAND1, initialized
C in driver code. States from last time-step are SACST_PRV(*) are stored
C using driver code.
c 6, 10, 10, 7, 1
      if (ipr==1.and.step1.eq..true.) then
         write(99,'2(f10.2,f10.3),2f10.2,10f10.3,f10.4,
     &        5f10.2,2f10.1,f10.3,f10.4,8f10.2') uztwc,uzfwc,lztwc,
     &        lzfsc,lzfpc,adimc,sperc,rimp,sdro,ssurf,sif,bfs,bfp,tci,
     &        edmnd,tet,pxv,uztwh,uzfwh,lztwh,lzfsh,lzfph,we,sh,aesc,ta,
     &        frost,frzst(1),frzst(2),frzst(3),frzst(4),frzst(5),frzdup,
     &        frzdbt
      endif
      step1 = .false.

      CALL FLAND1(PXV,EDMND,TA,WE,AESC,SH,DT,SACST,FRZST,SACPAR,
     +     FRZPAR,NSOIL,NUPL,NSAC,IVERS,SURF,GRND,TET,SMC,SH2O,
     +     SACST_PRV,DTFRZ,IDTFRZ,FRZDUP,FRZDBT,FROST,
     +     TSINT, SWINT, SWHINT, DSINT, NDSINT, DSINTW, NDINTW,
c swap out line to run unmodified source code
c     +     NORMALIZE,LAND_ID,iband,ipr)
     +     NORMALIZE)

      if (ipr ==1) write(*,*)'BELOW CALL TO FLAND1'
      if (ipr ==1) write(*,*)'PXV',PXV
      if (ipr ==1) write(*,*)'EDMND',EDMND
      if (ipr ==1) write(*,*)'TA',TA
      if (ipr ==1) write(*,*)'WE',WE
      if (ipr ==1) write(*,*)'AESC',AESC
      if (ipr ==1) write(*,*)'SH',SH
      if (ipr ==1) write(*,*)'DT',DT
      if (ipr ==1) then
         DO I = 1, 6
            write(*,*)'SACST:',i,SACST(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, 10
            write(*,*)'FRZST:',i,FRZST(i)
         ENDDO
      endif      
      if (ipr ==1) then
         DO I = 1, 17
            write(*,*)'SACPAR:',i,SACPAR(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, 14
            write(*,*)'FRZPAR:',i,FRZPAR(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'NSOIL',NSOIL
      if (ipr ==1) write(*,*)'NUPL',NUPL
      if (ipr ==1) write(*,*)'NSAC',NSAC
      if (ipr ==1) write(*,*)'IVERS',IVERS
      if (ipr ==1) write(*,*)'SURF',SURF
      if (ipr ==1) write(*,*)'GRND',GRND
      if (ipr ==1) write(*,*)'TET',TET
      if (ipr ==1) then
         DO I = 1, NSOIL
            write(*,*)'SMC:',i,SMC(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, NSOIL
            write(*,*)'SH2O:',i,SH2O(i)
         ENDDO
      endif      
      if (ipr ==1) then
         DO I = 1, 6
            write(*,*)'SACST_PRV:',i,SACST_PRV(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'DTFRZ',DTFRZ
      if (ipr ==1) write(*,*)'IDTFRZ',IDTFRZ
      if (ipr ==1) write(*,*)'FRZDUP',FRZDUP
      if (ipr ==1) write(*,*)'FRZDBT',FRZDBT
      if (ipr ==1) write(*,*)'FROST',FROST
      if (ipr ==1) then
         DO I = 1, NDSINT
            write(*,*)'TSINT:',i,TSINT(i)
         ENDDO
      endif
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'SWINT:',i,SWINT(i)
         ENDDO
      endif
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'SWHINT:',i,SWHINT(i)
         ENDDO
      endif
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'SWHINT:',i,SWHINT(i)
         ENDDO
      endif
      if (ipr ==1) then
         DO I = 1, NDSINT
            write(*,*)'DSINT:',i,DSINT(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'NDSINT',NDSINT
      if (ipr==1) then
         DO I = 1, NDINTW
            write(*,*)'DSINTW:',i,DSINTW(i)
         ENDDO
      endif
      if (ipr ==1) write(*,*)'NDINTW',NDINTW
      if (ipr ==1) write(*,*)'NORMALIZE',NORMALIZE

C     Re-assign states
      UZTWC = SACST(1)
      UZFWC = SACST(2)
      LZTWC = SACST(3)    
      LZFSC = SACST(4)
      LZFPC = SACST(5)
      ADIMC = SACST(6)
      SACST_PRV(1) = UZTWC
      SACST_PRV(2) = UZFWC
      SACST_PRV(3) = LZTWC    
      SACST_PRV(4) = LZFSC
      SACST_PRV(5) = LZFPC
      SACST_PRV(6) = ADIMC  

      DO I = 1,10
         if (ipr == 1) then
            if (i.gt.5) then
               write(*,*)'postfland',I,FRZST(I),SACST(I-5)
            else
               write(*,*)'postfland',I,FRZST(I)
            endif
         endif
      ENDDO
  
C     COMPUTE FINAL TOTAL STORAGE AND WATER BALANCE

      QS = SURF
      QG = GRND
      Q = SURF + GRND

      TOTAL_S2 = UZTWC + UZFWC + LZTWC + LZFSC + LZFPC + ADIMC
      DS = (TOTAL_S2 - TOTAL_S1)

      BAL = P1-ETA-QS-QG-DS
C	  WRITE(*,*)''
C	  WRITE(*,*)'FRZST1-3',FRZST(1),FRZST(2),FRZST(3)
C          WRITE(*,*)'FRZST4-6',FRZST(4),FRZST(5),FRZST(6)
C          WRITE(*,*)'FRZST7-9',FRZST(7),FRZST(8),FRZST(9)
C          WRITE(*,*)'FRZST10',FRZST(10)
C	  WRITE(*,*)'SMC1-2',SMC(1),SMC(2)
C	  WRITE(*,*)'SMC3-4',SMC(3),SMC(4)
C	  WRITE(*,*)'SH2O1-2',SH2O(1),SH2O(2)
C	  WRITE(*,*)'SH2O3-4',SH2O(3),SH2O(4)	 
C      PRINT*,'exsac1 -',BAL,P1,ETA,Q,QS,QG,DS,TOTAL_S1,TOTAL_S2

CC----------------------------------------------------------------------
CC Old call to SAC
C      CALL SAC1(DT,P1,EP1,TCI,ROIMP,SDRO,SSUR,SIF,BFS,BFP,ETA,
CC     SAC FROZEN GROUND VARIABLES
C     &            IFRZE,TA,LWE,WE,ISC,AESC,
CC     SAC PARAMETERS
C     &            UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
C     &            REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
C     &            SIDE,RSERV,
CC     SAC State variables  ',
C     &            UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC)
C
C     QS = ROIMP + SDRO + SSUR + SIF
C     QG = BFS + BFP
C     Q  = TCI
C
CC----------------------------------------------------------------------

      RETURN
      END

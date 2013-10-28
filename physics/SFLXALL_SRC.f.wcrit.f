      SUBROUTINE SFLX (
     C  ICE,DT,ZLVL,NSOIL,SLDPTH,
     F  LWDN,SOLDN,SOLNET,SFCPRS,PRCP,SFCTMP,Q2,SFCSPD,
     I  TH2,Q2SAT,DQSDT2,
     S  VEGTYP,SOILTYP,SLOPETYP,SHDFAC,SHDMAX,SHDMIN,PTU,ALB,SNOALB,TBOT,
     H  CMC,T1,TPACK,PACH20,STC,SMC,SH2O,SNOWH,SNEQV,ALBEDO,CH,CM,
     S  FRZST,FRZPAR,SACST,SACPAR,
     O  ETA,SHEAT,
C State and output variables for PEMB
     &  NSNOW,DSNOW,PSNOW,RTTSNOW,RTTDTSNOW,WTSNOW,WTDTSNOW,TTSNOW,
     &  TTDTSNOW,
C ----------------------------------------------------------------------
C OUTPUTS, DIAGNOSTICS, PARAMETERS BELOW GENERALLY NOT NECESSARY WHEN
C COUPLED WITH E.G. A NWP MODEL (SUCH AS THE NOAA/NWS/NCEP MESOSCALE ETA
C MODEL).  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES. 
C ----------------------------------------------------------------------
     O  EC,EDIR,ET,ETT,ESNOW,DRIP,DEW,
     O  BETA,ETP,SSOIL,
     O  FLX1,FLX2,FLX3,
     O  SNOMLT,SNCOVR,
     O  RUNOFF1,RUNOFF2,RUNOFF3,
     O  RC,PC,RSMIN,XLAI1,RCS,RCT,RCQ,RCSOIL,FX,
     D  SOILW,SOILM,RGL,HS,ERFLAG,PRFLAG,CELLID,MSTEP,
     P  SMCWLT,SMCDRY,SMCREF,SMCMAX,NROOT,CZMODEL,LSTSNW1,MODEL_TYPE,
     P  LVRAIN,WCRIT)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SFLX - VERSION 2.7 - July 6th 2004
C ----------------------------------------------------------------------
C SUB-DRIVER FOR "NOAH/OSU LSM" FAMILY OF PHYSICS SUBROUTINES FOR A
C SOIL/VEG/SNOWPACK LAND-SURFACE MODEL TO UPDATE SOIL MOISTURE, SOIL
C ICE, SOIL TEMPERATURE, SKIN TEMPERATURE, SNOWPACK WATER CONTENT,
C SNOWDEPTH, AND ALL TERMS OF THE SURFACE ENERGY BALANCE AND SURFACE
C WATER BALANCE (EXCLUDING INPUT ATMOSPHERIC FORCINGS OF DOWNWARD
C RADIATION AND PRECIP)
C ----------------------------------------------------------------------
C SFLX ARGUMENT LIST KEY:
C ----------------------------------------------------------------------
C  C  CONFIGURATION INFORMATION
C  F  FORCING DATA
C  I  OTHER (INPUT) FORCING DATA
C  S  SURFACE CHARACTERISTICS
C  H  HISTORY (STATE) VARIABLES
C  O  OUTPUT VARIABLES
C  D  DIAGNOSTIC OUTPUT
C ----------------------------------------------------------------------
C 1. CONFIGURATION INFORMATION (C):
C ----------------------------------------------------------------------
C   ICE	       SEA-ICE FLAG  (=1: SEA-ICE, =0: LAND)
C   DT	       TIMESTEP (SEC) (DT SHOULD NOT EXCEED 3600 SECS, RECOMMEND
C                1800 SECS OR LESS)
C   ZLVL       HEIGHT (M) ABOVE GROUND OF ATMOSPHERIC FORCING VARIABLES
C   NSOIL      NUMBER OF SOIL LAYERS (AT LEAST 2, AND NOT GREATER THAN
C                PARAMETER NSOLD SET BELOW)
C   SLDPTH     THE THICKNESS OF EACH SOIL LAYER (M)
C ----------------------------------------------------------------------
C 2. FORCING DATA (F):
C ----------------------------------------------------------------------
C   LWDN       LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET LONGWAVE)
C   SOLDN      SOLAR DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET SOLAR)
C   SFCPRS     PRESSURE AT HEIGHT ZLVL ABOVE GROUND (PASCALS)
C   PRCP       PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE)
C   SFCTMP     AIR TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
C   TH2        AIR POTENTIAL TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
C   Q2         MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
C ----------------------------------------------------------------------
C 3. OTHER FORCING (INPUT) DATA (I):
C ----------------------------------------------------------------------
C   SFCSPD     WIND SPEED (M S-1) AT HEIGHT ZLVL ABOVE GROUND
C   Q2SAT      SAT MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
C   DQSDT2     SLOPE OF SAT SPECIFIC HUMIDITY CURVE AT T=SFCTMP
C                (KG KG-1 K-1)
C ----------------------------------------------------------------------
C 4. CANOPY/SOIL CHARACTERISTICS (S):
C ----------------------------------------------------------------------
C   VEGTYP     VEGETATION TYPE (INTEGER INDEX)
C   SOILTYP    SOIL TYPE (INTEGER INDEX)
C   SLOPETYP   CLASS OF SFC SLOPE (INTEGER INDEX)
C   SHDFAC     AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C                (FRACTION= 0.0-1.0)
C   SHDMIN     MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C                (FRACTION= 0.0-1.0) <= SHDFAC
C   PTU        PHOTO THERMAL UNIT (PLANT PHENOLOGY FOR ANNUALS/CROPS)
C                (NOT YET USED, BUT PASSED TO REDPRM FOR FUTURE USE IN
C                VEG PARMS)
C   ALB        BACKROUND SNOW-FREE SURFACE ALBEDO (FRACTION), FOR JULIAN
C                DAY OF YEAR (USUALLY FROM TEMPORAL INTERPOLATION OF
C                MONTHLY MEAN VALUES' CALLING PROG MAY OR MAY NOT
C                INCLUDE DIURNAL SUN ANGLE EFFECT)
C   SNOALB     UPPER BOUND ON MAXIMUM ALBEDO OVER DEEP SNOW (E.G. FROM
C                ROBINSON AND KUKLA, 1985, J. CLIM. & APPL. METEOR.)
C   TBOT       BOTTOM SOIL TEMPERATURE (LOCAL YEARLY-MEAN SFC AIR
C                TEMPERATURE)
C ----------------------------------------------------------------------
C 5. HISTORY (STATE) VARIABLES (H):
C ----------------------------------------------------------------------
C  CMC         CANOPY MOISTURE CONTENT (M)
C  T1          GROUND/CANOPY/SNOWPACK) EFFECTIVE SKIN TEMPERATURE (K)
C  TPACK       EFFECTIVE SNOWPACK TEMPERATURE (K)
C  PACH20      SNOWPACK LIQUID WATER CONTENT (M)
C  STC(NSOIL)  SOIL TEMP (K)
C  SMC(NSOIL)  TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
C  SH2O(NSOIL) UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
C                NOTE: FROZEN SOIL MOISTURE = SMC - SH2O
C  SNOWH       ACTUAL SNOW DEPTH (M)
C  SNEQV       LIQUID WATER-EQUIVALENT SNOW DEPTH (M)
C                NOTE: SNOW DENSITY = SNEQV/SNOWH
C  ALBEDO      SURFACE ALBEDO INCLUDING SNOW EFFECT (UNITLESS FRACTION)
C                =SNOW-FREE ALBEDO (ALB) WHEN SNEQV=0, OR
C                =FCT(MSNOALB,ALB,VEGTYP,SHDFAC,SHDMIN) WHEN SNEQV>0
C  CH          SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
C                (M S-1); NOTE: CH IS TECHNICALLY A CONDUCTANCE SINCE
C                IT HAS BEEN MULTIPLIED BY WIND SPEED.
C  CM          SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM (M S-1); NOTE:
C                CM IS TECHNICALLY A CONDUCTANCE SINCE IT HAS BEEN
C                MULTIPLIED BY WIND SPEED.  CM IS NOT NEEDED IN SFLX
C ----------------------------------------------------------------------
C 6. SAC-SMA PARAMETERS (S):
C ----------------------------------------------------------------------
C  UZTW   Upper zone tension water (M:max, C:content)
C  UZFW   Upper zone free water (M:max, C:content)
C  LZTW   Lower zone tension water (M:max, C:content)
C  LZFP   Lower zone free primary water (M:max, C:content)
C  LZFS   Lower zone free supplementary water (M:max, C:content)
C  UZK    Interflow depletion rate from uppler layer
C  LZSK   Depletion rate of lower zone supplemental storage
C  LZPK   Depletion rate of lower zone primary storage
C  ZPERC  Ratio of max and min percolation rates
C  REXP   Shape parameter of percolation curve
C  PFREE  Percolation fraction that goes directly to lower free water
C  PCTIM  Permanent impervious area fraction
C  ADIMP  Maximum frac. of additional impervious area due to saturation
C  RIVA   Riparian vegetation area fraction
C  SIDE   Ratio of deep percolation from lower layer free storages
C  RSERV  Fraction of lower free water not transferable to tension water
C ----------------------------------------------------------------------
C 7. OUTPUT (O):
C ----------------------------------------------------------------------
C OUTPUT VARIABLES NECESSARY FOR A COUPLED NUMERICAL WEATHER PREDICTION
C MODEL, E.G. NOAA/NWS/NCEP MESOSCALE ETA MODEL.  FOR THIS APPLICATION,
C THE REMAINING OUTPUT/DIAGNOSTIC/PARAMETER BLOCKS BELOW ARE NOT
C NECESSARY.  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES.
C   ETA        ACTUAL LATENT HEAT FLUX (W M-2: NEGATIVE, IF UP FROM
C	         SURFACE)
C   SHEAT      SENSIBLE HEAT FLUX (W M-2: NEGATIVE, IF UPWARD FROM
C	         SURFACE)
C ----------------------------------------------------------------------
C   EC         CANOPY WATER EVAPORATION (W M-2)
C   EDIR       DIRECT SOIL EVAPORATION (W M-2)
C   ET(NSOIL)  PLANT TRANSPIRATION FROM A PARTICULAR ROOT (SOIL) LAYER
C                 (W M-2)
C   ETT        TOTAL PLANT TRANSPIRATION (W M-2)
C   ESNOW      SUBLIMATION FROM SNOWPACK (W M-2)
C   DRIP       THROUGH-FALL OF PRECIP AND/OR DEW IN EXCESS OF CANOPY
C                WATER-HOLDING CAPACITY (M)
C   DEW        DEWFALL (OR FROSTFALL FOR T<273.15) (M)
C ----------------------------------------------------------------------
C   BETA       RATIO OF ACTUAL/POTENTIAL EVAP (DIMENSIONLESS)
C   ETP        POTENTIAL EVAPORATION (W M-2)
C   SSOIL      SOIL HEAT FLUX (W M-2: NEGATIVE IF DOWNWARD FROM SURFACE)
C ----------------------------------------------------------------------
C   FLX1       PRECIP-SNOW SFC (W M-2)
C   FLX2       FREEZING RAIN LATENT HEAT FLUX (W M-2)
C   FLX3       PHASE-CHANGE HEAT FLUX FROM SNOWMELT (W M-2)
C ----------------------------------------------------------------------
C   SNOMLT     SNOW MELT (M) (WATER EQUIVALENT)
C   SNCOVR     FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1)
C ----------------------------------------------------------------------
C   RUNOFF1    SURFACE RUNOFF (M S-1), NOT INFILTRATING THE SURFACE
C   RUNOFF2    SUBSURFACE RUNOFF (M S-1), DRAINAGE OUT BOTTOM OF LAST
C                SOIL LAYER
C   RUNOFF3    NUMERICAL TRUNCTATION IN EXCESS OF POROSITY (SMCMAX)
C                FOR A GIVEN SOIL LAYER AT THE END OF A TIME STEP
C ----------------------------------------------------------------------
C   RC         CANOPY RESISTANCE (S M-1)
C   PC         PLANT COEFFICIENT (UNITLESS FRACTION, 0-1) WHERE PC*ETP
C                = ACTUAL TRANSPIRATION
C   XLAI       LEAF AREA INDEX (DIMENSIONLESS)
C   RSMIN      MINIMUM CANOPY RESISTANCE (S M-1)
C   RCS        INCOMING SOLAR RC FACTOR (DIMENSIONLESS)
C   RCT        AIR TEMPERATURE RC FACTOR (DIMENSIONLESS)
C   RCQ        ATMOS VAPOR PRESSURE DEFICIT RC FACTOR (DIMENSIONLESS)
C   RCSOIL     SOIL MOISTURE RC FACTOR (DIMENSIONLESS)
C ----------------------------------------------------------------------
C 8. DIAGNOSTIC OUTPUT (D):
C ----------------------------------------------------------------------
C   SOILW      AVAILABLE SOIL MOISTURE IN ROOT ZONE (UNITLESS FRACTION
C	         BETWEEN SMCWLT AND SMCMAX)
C   SOILM      TOTAL SOIL COLUMN MOISTURE CONTENT (FROZEN+UNFROZEN) (M) 
C ----------------------------------------------------------------------
C 9. PARAMETERS (P):
C ----------------------------------------------------------------------
C   SMCWLT     WILTING POINT (VOLUMETRIC)
C   SMCDRY     DRY SOIL MOISTURE THRESHOLD WHERE DIRECT EVAP FRM TOP
C                LAYER ENDS (VOLUMETRIC)
C   SMCREF     SOIL MOISTURE THRESHOLD WHERE TRANSPIRATION BEGINS TO
C                STRESS (VOLUMETRIC)
C   SMCMAX     POROSITY, I.E. SATURATED VALUE OF SOIL MOISTURE
C                (VOLUMETRIC)
C   NROOT      NUMBER OF ROOT LAYERS, A FUNCTION OF VEG TYPE, DETERMINED
C              IN SUBROUTINE REDPRM.
C ----------------------------------------------------------------------
C 10. Ben's parameters for Albedo decay calculation
C    FFROZP    FROZEN PRECIPITATION FLAG
C    T1       EFFECTIVE SNOW/GROUND SURFACE TERM (K)
     
C ----------------------------------------------------------------------- 
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

C ----------------------------------------------------------------------
C DECLARATIONS - LOGICAL
C ----------------------------------------------------------------------
      LOGICAL FRZGRA
      LOGICAL SATURATED
      LOGICAL SNOWNG
      LOGICAL LVRAIN

C ----------------------------------------------------------------------
C DECLARATIONS - INTEGER
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER K,I
      INTEGER KZ
      INTEGER NSOIL
      INTEGER NROOT
      INTEGER SLOPETYP
      INTEGER SOILTYP
      INTEGER VEGTYP

C ----------------------------------------------------------------------
C DECLARATIONS - REAL
C ----------------------------------------------------------------------
      REAL ALBEDO
      REAL ALB
      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CH
      REAL CM
      REAL CMC
      REAL CMCMAX
      REAL CP
      REAL CSNOW
      REAL CSOIL
      REAL CZIL
      REAL DEW
      REAL DF1
      REAL DF1H
      REAL DF1A
      REAL DKSAT
      REAL DT
      REAL DWSAT
      REAL DQSDT2
      REAL DSOIL
      REAL DTOT
      REAL DRIP
      REAL EC
      REAL EDIR
      REAL ESNOW
      REAL ET(NSOIL)
      REAL ETT
      REAL FRCSNO
      REAL FRCSOI
      REAL EPSCA
      REAL ETA
      REAL ETP
      REAL FDOWN
      REAL F1
      REAL FLX1
      REAL FLX2
      REAL FLX3
      REAL FXEXP
      REAL FX
      REAL FRZX
      REAL SHEAT
      REAL HS
      REAL KDT
      REAL LWDN
      REAL LVH2O
      REAL PC
      REAL PRCP
      REAL PTU
      REAL PRCP1
      REAL PSISAT
      REAL Q2
      REAL Q2SAT
      REAL QUARTZ
      REAL R
      REAL RCH
      REAL REFKDT
      REAL RR
      REAL RTDIS(NSOLD)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RGL
      REAL RUNOFF3
      REAL RSMAX
      REAL RC
      REAL RSMIN
      REAL RCQ
      REAL RCS
      REAL RCSOIL
      REAL RCT
      REAL RSNOW
      REAL SNDENS
      REAL SNCOND 
      REAL SSOIL
      REAL SBETA
      REAL SFCPRS
      REAL SFCSPD
      REAL SFCTMP
      REAL SHDFAC
      REAL SHDMAX
      REAL SHDMIN
      REAL SH2O(NSOIL)
      REAL SLDPTH(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL SMC(NSOIL)
      REAL SNEQV
      REAL SNCOVR
      REAL SNOWH
      REAL SN_NEW
      REAL SLOPE
      REAL SNUP
      REAL SALP
      REAL SNOALB
      REAL STC(NSOIL)
      REAL SNOMLT
      REAL SOLDN
      REAL SOILM
      REAL SOILW
      REAL SOILWM
      REAL SOILWW
      REAL T1
      REAL T1V
      REAL T24
      REAL T2V
      REAL TBOT
      REAL TH2
      REAL TH2V
      REAL TOPT
      REAL TFREEZ
      REAL TSNOW
      REAL XLAI
      REAL XLAI1
      REAL ZLVL
      REAL ZBOT
      REAL Z0
      REAL ZSOIL(NSOLD)
      REAL EMISS
      REAL EMISSNOW

      REAL FFROZP
      REAL SOLNET
      REAL LSUBS
      REAL CZMODEL,PI,ALBEDO0,ALBEDO1,ALBEDO2,ALB1

      INTEGER LSTSNW1
      INTEGER LSTSNW	
      INTEGER ERFLAG,PRFLAG,CELLID
      INTEGER MODEL_TYPE,MSTEP
      REAL RIB
      REAL TPACK
      REAL PACH20
      REAL UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,LZFSM
      REAL LZFPM,LZSK,LZPK,PFREE,SIDE,RSERV,UZTWC,UZFWC,LZTWC,LZFSC
      REAL LZFPC,ADIMC
      REAL NSNOW,DSNOW,PSNOW,RTTSNOW,RTTDTSNOW,WTSNOW,WTDTSNOW,TTSNOW
      REAL TTDTSNOW,TW
      REAL FRZST(10)
      REAL FRZPAR(13)
      REAL SACST(6)
      REAL SACPAR(16)
      REAL FROST
      real smctot0,sh2otot0,runtot,smctot,sh2otot,runtot0,smcdif,sh2odif
      real wcrit

C ----------------------------------------------------------------------
C DECLARATIONS - PARAMETERS
C ----------------------------------------------------------------------
      PARAMETER(TFREEZ = 273.15)
      PARAMETER(LVH2O = 2.501E+6)
      PARAMETER(LSUBS = 2.83E+6)
      PARAMETER(R = 287.04)
      PARAMETER(CP = 1004.5)

C ----------------------------------------------------------------------
C   INITIALIZATION
C ----------------------------------------------------------------------
      RUNOFF1 = 0.0
      RUNOFF2 = 0.0
      RUNOFF3 = 0.0
      SNOMLT = 0.0
      ERFLAG = 0

C      SFCTMP = 283.15
C ----------------------------------------------------------------------
C  THE VARIABLE "ICE" IS A FLAG DENOTING SEA-ICE CASE 
C ----------------------------------------------------------------------
      IF (ICE .EQ. 1) THEN

C ----------------------------------------------------------------------
C SEA-ICE LAYERS ARE EQUAL THICKNESS AND SUM TO 3 METERS
C ----------------------------------------------------------------------
        DO KZ = 1,NSOIL
          ZSOIL(KZ) = -3.*FLOAT(KZ)/FLOAT(NSOIL)
        END DO

      ELSE

C ----------------------------------------------------------------------
C CALCULATE DEPTH (NEGATIVE) BELOW GROUND FROM TOP SKIN SFC TO BOTTOM OF
C   EACH SOIL LAYER.  NOTE:  SIGN OF ZSOIL IS NEGATIVE (DENOTING BELOW
C   GROUND)
C ----------------------------------------------------------------------
        ZSOIL(1) = -SLDPTH(1)
        DO KZ = 2,NSOIL
           ZSOIL(KZ) = -SLDPTH(KZ)+ZSOIL(KZ-1)
        END DO

      ENDIF
      smctot0=0.
      sh2otot0=0.
      runtot0=0.
      runtot=(runoff1+runoff2)*-1
      do i = 1,nsoil
         smctot0 = smctot0 + smc(i)*sldpth(i)
         sh2otot0 = sh2otot0 + sh2o(i)*sldpth(i)
      enddo    
      if (prflag==2) then
      write(*,*)'-----------sflx top------------' 
      write(*,*)'prcp           edir            ec'
      write(*,*)prcp,edir,ec
      write(*,*)'ett            eta             etp'
      write(*,*)ett,eta,etp
      write(*,*)'runoff          swe'
      write(*,*)runtot0,sneqv
      write(*,*)'smctot0        sh2otot0        cmc'
      write(*,*)smctot0,sh2otot0,cmc
      write(*,*)'--------------------------------' 
      endif
C ----------------------------------------------------------------------
C NEXT IS CRUCIAL CALL TO SET THE LAND-SURFACE PARAMETERS, INCLUDING
C SOIL-TYPE AND VEG-TYPE DEPENDENT PARAMETERS.
C ----------------------------------------------------------------------

      CALL REDPRM (VEGTYP,SOILTYP,SLOPETYP,
     +     CFACTR,CMCMAX,RSMAX,TOPT,REFKDT,KDT,SBETA,
     O     SHDFAC,SHDMAX,SHDMIN,RSMIN,RGL,HS,ZBOT,FRZX,PSISAT,SLOPE,
     +     SNUP,SALP,BEXP,DKSAT,DWSAT,SMCMAX,SMCWLT,SMCREF,
     O     SMCDRY,F1,QUARTZ,FXEXP,RTDIS,SLDPTH,ZSOIL,
     +     NROOT,NSOIL,Z0,CZIL,XLAI,CSOIL,PTU,EMISS,EMISSNOW,MODEL_TYPE)

      XLAI1=XLAI

C ----------------------------------------------------------------------
C  INITIALIZE PRECIPITATION LOGICALS.
C ----------------------------------------------------------------------
      SNOWNG = .FALSE.
      FRZGRA = .FALSE.

C ----------------------------------------------------------------------
C IF SEA-ICE CASE, ASSIGN DEFAULT WATER-EQUIV SNOW ON TOP
C ----------------------------------------------------------------------
      IF (ICE .EQ. 1) THEN
        SNEQV = 0.01
        SNOWH = 0.05
        SNDENS = SNEQV/SNOWH
      ENDIF

C ----------------------------------------------------------------------
C IF INPUT SNOWPACK IS NONZERO, THEN COMPUTE SNOW DENSITY "SNDENS" AND
C   SNOW THERMAL CONDUCTIVITY "SNCOND" (NOTE THAT CSNOW IS A FUNCTION
C   SUBROUTINE)
C ----------------------------------------------------------------------
      IF (SNEQV .EQ. 0.0) THEN
        SNDENS = 0.0
        SNOWH = 0.0
        SNCOND = 1.0
      ELSE
        SNDENS = SNEQV/SNOWH
        SNCOND = CSNOW(SNDENS) 
      ENDIF

C ----------------------------------------------------------------------
C DETERMINE IF IT'S PRECIPITATING AND WHAT KIND OF PRECIP IT IS.
C IF IT'S PRCPING AND THE AIR TEMP IS COLDER THAN 0 C, IT'S SNOWING!
C IF IT'S PRCPING AND THE AIR TEMP IS WARMER THAN 0 C, BUT THE GRND
C TEMP IS COLDER THAN 0 C, FREEZING RAIN IS PRESUMED TO BE FALLING.
C ----------------------------------------------------------------------
      IF (PRCP .GT. 0.0) THEN
C         CALL WETBULB(TW,SFCTMP,Q2,SFCPRS)
        IF (SFCTMP .LE. TFREEZ) THEN
C        IF (FFROZP .GT. 0.5) THEN
C         IF (TW.LE.(TFREEZ + 1)) THEN
            SNOWNG = .TRUE.
C     WRITE(*,*)'SFCTMP PRCP',SFCTMP,PRCP
C            WRITE(*,*)'SNOWNG FLAG'
         ELSE
            IF (T1 .LE. TFREEZ) FRZGRA = .TRUE.
C            WRITE(*,*)'FRZGRA FLAG'
         ENDIF
      ENDIF

C ----------------------------------------------------------------------
C IF EITHER PRCP FLAG IS SET, DETERMINE NEW SNOWFALL (CONVERTING PRCP
C RATE FROM KG M-2 S-1 TO A LIQUID EQUIV SNOW DEPTH IN METERS) AND ADD
C IT TO THE EXISTING SNOWPACK.
C NOTE THAT SINCE ALL PRECIP IS ADDED TO SNOWPACK, NO PRECIP INFILTRATES
C INTO THE SOIL SO THAT PRCP1 IS SET TO ZERO.
C ----------------------------------------------------------------------
C Note for the Anderson PEMB model, the below additions are only utlized
C on the first time step with snow, to initialize the pack.  Otherwise
C PEMB internal arrays hold pack ice/water/temperature values and precip 
C type is determined within the TR_19 code based on wet-bulb temperature. 
C ----------------------------------------------------------------------
C For the case of zero snow cover (i.e. initial conditions) if there is
C no SWE present, do not allow SWE to form if only the FRZGRA flag is
C set, since this is unrealistic and will not persist (e.g. wet-bulb 
C temp > 1 deg C, but ground is frozen)

      IF (SNEQV == 0.0.AND.SNOWNG.EQ..FALSE.) THEN
         FRZGRA = .FALSE.
      END IF
      IF ( (SNOWNG) .OR. (FRZGRA) ) THEN
        SN_NEW = PRCP * DT * 0.001
        SNEQV = SNEQV + SN_NEW
C		WRITE(*,*)'SNEQV SN_NEW',SNEQV,SN_NEW
        PRCP1 = 0.0

C ----------------------------------------------------------------------
C UPDATE SNOW DENSITY BASED ON NEW SNOWFALL, USING OLD AND NEW SNOW.
C UPDATE SNOW THERMAL CONDUCTIVITY
C ----------------------------------------------------------------------
        CALL SNOW_NEW (SFCTMP,SN_NEW,SNOWH,SNDENS)

        SNCOND = CSNOW (SNDENS) 
      ELSE

C ----------------------------------------------------------------------
C PRECIP IS LIQUID (RAIN), HENCE SAVE IN THE PRECIP VARIABLE THAT
C LATER CAN WHOLELY OR PARTIALLY INFILTRATE THE SOIL (ALONG WITH 
C ANY CANOPY "DRIP" ADDED TO THIS LATER)
C ----------------------------------------------------------------------
        PRCP1 = PRCP

      ENDIF

C ----------------------------------------------------------------------
C DETERMINE SNOWCOVER AND ALBEDO OVER LAND.
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C DETERMINE SNOWCOVER AND ALBEDO OVER LAND.
C ----------------------------------------------------------------------
CHELIN    diurnal variation of albedo
      ALB1=ALB
c     IF(ALB.LE.0.5)THEN
c        PI= 3.141592
c        ALBEDO0=-18.0*(0.5-ACOS(CZMODEL)/PI)
c        ALBEDO0=EXP(ALBEDO0)
c        ALBEDO1=(ALB1-0.054313)/0.945687
c        
c        ALBEDO2=ALBEDO1+(1-ALBEDO1)*ALBEDO0
c        ALBEDO2=MIN(ALBEDO2,3*ALB)
         ALBEDO2=ALB1*(2.287-3.374*CZMODEL+3.619*CZMODEL*CZMODEL
     1           -1.603**CZMODEL*CZMODEL*CZMODEL)
c weighting so diurnal adjustment only applied for direct beam
         ALB1=0.45*ALB+0.55*ALBEDO2
         IF(ALB1.LT.0)ALB1=0
c     ENDIF
 
      IF (ICE .EQ. 0) THEN

C ----------------------------------------------------------------------
C IF SNOW DEPTH=0, SET SNOW FRACTION=0, ALBEDO=SNOW FREE ALBEDO.
C ----------------------------------------------------------------------
         IF (SNEQV .EQ. 0.0) THEN
            SNCOVR = 0.0
            ALBEDO = ALB1
         ELSE
C ----------------------------------------------------------------------
C DETERMINE SNOW FRACTIONAL COVERAGE.
C DETERMINE SURFACE ALBEDO MODIFICATION DUE TO SNOWDEPTH STATE.
C ----------------------------------------------------------------------
            CALL SNFRAC (SNEQV,SNUP,SALP,SNOWH,SNCOVR)

C Noah default albedo algorithm

C            CALL ALCALC (ALB1,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO)
C ----------------------------------------------------------------------
C Ben Livneh 2007; Alternate albedo scheme, considering a higher initial
C albedo, with an exponential decay function

       CALL ALBCPU (ALB1,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO,
     &                   SNOWNG,TFREEZ,DT,SNEQV,T1,LSTSNW1,VEGTYP)

C ----------------------------------------------------------------------
         ENDIF

      ELSE
C ----------------------------------------------------------------------
C SNOW COVER, ALBEDO OVER SEA-ICE
C ----------------------------------------------------------------------
         SNCOVR = 1.0
C   changed in version 2.6 on June 2nd 2003
C     ALBEDO = 0.60
         ALBEDO = 0.65
      ENDIF

C ----------------------------------------------------------------------
C THERMAL CONDUCTIVITY FOR SEA-ICE CASE
C ----------------------------------------------------------------------
      IF (ICE .EQ. 1) THEN
         DF1 = 2.2
      ELSE

C ----------------------------------------------------------------------
C NEXT CALCULATE THE SUBSURFACE HEAT FLUX, WHICH FIRST REQUIRES
C CALCULATION OF THE THERMAL DIFFUSIVITY.  TREATMENT OF THE
C LATTER FOLLOWS THAT ON PAGES 148-149 FROM "HEAT TRANSFER IN 
C COLD CLIMATES", BY V. J. LUNARDINI (PUBLISHED IN 1981 
C BY VAN NOSTRAND REINHOLD CO.) I.E. TREATMENT OF TWO CONTIGUOUS 
C "PLANE PARALLEL" MEDIUMS (NAMELY HERE THE FIRST SOIL LAYER 
C AND THE SNOWPACK LAYER, IF ANY). THIS DIFFUSIVITY TREATMENT 
C BEHAVES WELL FOR BOTH ZERO AND NONZERO SNOWPACK, INCLUDING THE 
C LIMIT OF VERY THIN SNOWPACK.  THIS TREATMENT ALSO ELIMINATES
C THE NEED TO IMPOSE AN ARBITRARY UPPER BOUND ON SUBSURFACE 
C HEAT FLUX WHEN THE SNOWPACK BECOMES EXTREMELY THIN.
C ----------------------------------------------------------------------
C FIRST CALCULATE THERMAL DIFFUSIVITY OF TOP SOIL LAYER, USING
C BOTH THE FROZEN AND LIQUID SOIL MOISTURE, FOLLOWING THE 
C SOIL THERMAL DIFFUSIVITY FUNCTION OF PETERS-LIDARD ET AL.
C (1998,JAS, VOL 55, 1209-1224), WHICH REQUIRES THE SPECIFYING
C THE QUARTZ CONTENT OF THE GIVEN SOIL CLASS (SEE ROUTINE REDPRM)
C ----------------------------------------------------------------------
        CALL TDFCND (DF1,SMC(1),QUARTZ,SMCMAX,SH2O(1))

C ----------------------------------------------------------------------
C NEXT ADD SUBSURFACE HEAT FLUX REDUCTION EFFECT FROM THE 
C OVERLYING GREEN CANOPY, ADAPTED FROM SECTION 2.1.2 OF 
C PETERS-LIDARD ET AL. (1997, JGR, VOL 102(D4))
C ----------------------------------------------------------------------
        DF1 = DF1 * EXP(SBETA*SHDFAC)
      ENDIF

C ----------------------------------------------------------------------
C FINALLY "PLANE PARALLEL" SNOWPACK EFFECT FOLLOWING 
C V.J. LINARDINI REFERENCE CITED ABOVE. NOTE THAT DTOT IS
C COMBINED DEPTH OF SNOWDEPTH AND THICKNESS OF FIRST SOIL LAYER
C ----------------------------------------------------------------------
      DSOIL = -(0.5 * ZSOIL(1))

      IF (SNEQV .EQ. 0.) THEN
         SSOIL = DF1 * (T1 - STC(1) ) / DSOIL
c         write(*,*)'sflx swe  ssoil df1 t1 stc1 dsoil'
c         write(*,*)ssoil,df1,t1,stc(1),dsoil
      ELSE
        DTOT = SNOWH + DSOIL
        FRCSNO = SNOWH/DTOT
        FRCSOI = DSOIL/DTOT
C
C 1. HARMONIC MEAN (SERIES FLOW)
        DF1H = (SNCOND*DF1)/(FRCSOI*SNCOND+FRCSNO*DF1)
C 2. ARITHMETIC MEAN (PARALLEL FLOW)
        DF1A = FRCSNO*SNCOND + FRCSOI*DF1
C 3. GEOMETRIC MEAN (INTERMEDIATE BETWEEN HARMONIC AND ARITHMETIC MEAN)
        DF1 = DF1A*SNCOVR + DF1*(1.0-SNCOVR)

C ----------------------------------------------------------------------
C CALCULATE SUBSURFACE HEAT FLUX, SSOIL, FROM FINAL THERMAL DIFFUSIVITY
C OF SURFACE MEDIUMS, DF1 ABOVE, AND SKIN TEMPERATURE AND TOP 
C MID-LAYER SOIL TEMPERATURE
C ----------------------------------------------------------------------
        SSOIL = DF1 * (T1 - STC(1) ) / DTOT
c        write(*,*)'sflx swe0  ssoil df1 t1 stc1 dtot'
c        write(*,*)ssoil,df1,t1,stc(1),dtot
      ENDIF

C ----------------------------------------------------------------------
C DETERMINE SURFACE ROUGHNESS OVER SNOWPACK USING SNOW CONDITION FROM
C THE PREVIOUS TIMESTEP.
C ----------------------------------------------------------------------
      IF (SNCOVR .GT. 0.) THEN
        CALL SNOWZ0 (SNCOVR,Z0)
      ENDIF

C ----------------------------------------------------------------------
C NEXT CALL ROUTINE SFCDIF TO CALCULATE THE SFC EXCHANGE COEF (CH) FOR
C HEAT AND MOISTURE.
C
C NOTE !!!
C COMMENT OUT CALL SFCDIF, IF SFCDIF ALREADY CALLED IN CALLING PROGRAM
C (SUCH AS IN COUPLED ATMOSPHERIC MODEL).
C
C NOTE !!!
C DO NOT CALL SFCDIF UNTIL AFTER ABOVE CALL TO REDPRM, IN CASE
C ALTERNATIVE VALUES OF ROUGHNESS LENGTH (Z0) AND ZILINTINKEVICH COEF
C (CZIL) ARE SET THERE VIA NAMELIST I/O.
C
C NOTE !!!
C ROUTINE SFCDIF RETURNS A CH THAT REPRESENTS THE WIND SPD TIMES THE
C "ORIGINAL" NONDIMENSIONAL "Ch" TYPICAL IN LITERATURE.  HENCE THE CH
C RETURNED FROM SFCDIF HAS UNITS OF M/S.  THE IMPORTANT COMPANION
C COEFFICIENT OF CH, CARRIED HERE AS "RCH", IS THE CH FROM SFCDIF TIMES
C AIR DENSITY AND PARAMETER "CP".  "RCH" IS COMPUTED IN "CALL PENMAN".
C RCH RATHER THAN CH IS THE COEFF USUALLY INVOKED LATER IN EQNS.
C
C NOTE !!!
C SFCDIF ALSO RETURNS THE SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM, CM,
C ALSO KNOWN AS THE SURFACE DRAGE COEFFICIENT, BUT CM IS NOT USED HERE.
C ----------------------------------------------------------------------
C CALC VIRTUAL TEMPS AND VIRTUAL POTENTIAL TEMPS NEEDED BY SUBROUTINES
C SFCDIF AND PENMAN.
C ----------------------------------------------------------------------
      T2V = SFCTMP * (1.0 + 0.61 * Q2 )
C ----------------------------------------------------------------------
C COMMENT OUT BELOW 2 LINES IF CALL SFCDIF IS COMMENTED OUT, I.E. IN THE
C COUPLED MODEL.
C ----------------------------------------------------------------------
      T1V = T1 * (1.0 + 0.61 * Q2)
      TH2V = TH2 * (1.0 + 0.61 * Q2)

      CALL SFCDIF (ZLVL,Z0,T1V,TH2V,SFCSPD,CZIL,CM,CH)


C -------- Slater Fix for NOAH  (changed by Youlong, March 2008)-------
      RIB = 9.8*ZLVL*(TH2V-T1V)/(SFCSPD*SFCSPD)
      IF (RIB .gt. 0.) THEN
      CH = CH * MAX(1. - RIB/0.5, 0.05)
      CM = CM * MAX(1. - RIB/0.5, 0.05)
      ENDIF	
C ----------------------------------------------------------------------
C CALCULATE TOTAL DOWNWARD RADIATION (SOLAR PLUS LONGWAVE) NEEDED IN
C PENMAN EP SUBROUTINE THAT FOLLOWS
C ----------------------------------------------------------------------
      FDOWN = SOLDN*(1.0-ALBEDO) + LWDN
C      FDOWN = SOLNET + LWDN

C ----------------------------------------------------------------------
C CALL PENMAN SUBROUTINE TO CALCULATE POTENTIAL EVAPORATION (ETP), AND
C OTHER PARTIAL PRODUCTS AND SUMS SAVE IN COMMON/RITE FOR LATER
C CALCULATIONS.
C ----------------------------------------------------------------------
       CALL PENMAN (SFCTMP,SFCPRS,CH,T2V,TH2,PRCP,FDOWN,T24,SSOIL,
     &              Q2,Q2SAT,ETP,RCH,EPSCA,RR,SNOWNG,FRZGRA,
     &              DQSDT2,FLX2,EMISS,EMISSNOW,SNCOVR)

      smctot=0.
      sh2otot=0.
      runtot=0.
      runtot=(runoff1+runoff2)*-1
      do i = 1,nsoil
         smctot = smctot + smc(i)*sldpth(i)
         sh2otot = sh2otot + sh2o(i)*sldpth(i)
      enddo  
      smcdif=smctot0-smctot
      sh2odif=sh2otot0-sh2otot
      if (prflag==2) then
      write(*,*)'-----------sflx BELOW PENMAN-----' 
      write(*,*)'prcp           edir            ec'
      write(*,*)prcp,edir,ec
      write(*,*)'ett            eta             etp'
      write(*,*)ett,eta,etp
      write(*,*)'runoff          swe'
      write(*,*)runtot,sneqv
      write(*,*)'smcdif         sh2odif         cmc'
      write(*,*)smcdif,sh2odif,cmc
      write(*,*)'--------------------------------' 
      endif
C ----------------------------------------------------------------------
C CALL CANRES TO CALCULATE THE CANOPY RESISTANCE AND CONVERT IT INTO PC
C IF NONZERO GREENNESS FRACTION
C ----------------------------------------------------------------------
chelin      IF (SHDFAC .GT. 0.) THEN
      
C ----------------------------------------------------------------------
C  FROZEN GROUND EXTENSION: TOTAL SOIL WATER "SMC" WAS REPLACED 
C  BY UNFROZEN SOIL WATER "SH2O" IN CALL TO CANRES BELOW
C ----------------------------------------------------------------------
        CALL CANRES (SOLDN,CH,SFCTMP,Q2,SFCPRS,SH2O,ZSOIL,NSOIL,
     &               SMCWLT,SMCREF,RSMIN,RC,PC,NROOT,Q2SAT,DQSDT2,
     &               TOPT,RSMAX,RGL,HS,XLAI,STC,
     &               RCS,RCT,RCQ,RCSOIL,EMISS,EMISSNOW,SNCOVR,RTDIS)

chelin      ENDIF

C ----------------------------------------------------------------------
C NOW DECIDE MAJOR PATHWAY BRANCH TO TAKE DEPENDING ON WHETHER SNOWPACK
C EXISTS OR NOT:
C ----------------------------------------------------------------------
        ESNOW = 0.0
        LVRAIN = .FALSE.
      IF (SNEQV .EQ. 0.0) THEN
        CALL NOPAC (ETP,ETA,PRCP,SMC,SMCMAX,SMCWLT,
     &        SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,SHDFAC,
     &        SBETA,Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,
     &        STC,EPSCA,BEXP,PC,RCH,RR,CFACTR,
     &        SH2O,SLOPE,KDT,FRZX,PSISAT,ZSOIL,
     &        DKSAT,DWSAT,TBOT,ZBOT,RUNOFF1,RUNOFF2,RUNOFF3,
     &        EDIR,EC,ET,ETT,NROOT,ICE,RTDIS,QUARTZ,FXEXP,FX,
     &        CSOIL,FRZST,FRZPAR,SACST,SACPAR,BETA,DRIP,DEW,
     &        FLX1,FLX2,FLX3,EMISS,MODEL_TYPE,PRFLAG,CELLID,MSTEP,wcrit)
      ELSE
         IF (MODEL_TYPE == 0.OR. MODEL_TYPE == 1) THEN
C            CALL SNOPAC (ETP,ETA,PRCP,PRCP1,SNOWNG,SMC,SMCMAX,SMCWLT,
C     &        SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
C     &        SBETA,DF1,
C     &        Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
C     &        SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,SNEQV,SNDENS,
C     &        SNOWH,SH2O,SLOPE,KDT,FRZX,PSISAT,SNUP,
C     &        ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
C     &        RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
C     &        ICE,RTDIS,QUARTZ,FXEXP,FX,CSOIL,
C     &        BETA,DRIP,DEW,FLX1,FLX2,FLX3,ESNOW,EMISS,EMISSNOW,
C     &        MODEL_TYPE,FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP)
            CALL SNWPAC (ETP,ETA,PRCP,PRCP1,SNOWNG,FRZGRA,SMC,SMCMAX,
     &           SMCWLT,SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &           SBETA,DF1,PACH20,CH,SFCSPD,
     &           Q2,T1,TPACK,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &           SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,SNEQV,SNDENS,
     &           SNOWH,SH2O,SLOPE,KDT,FRZX,PSISAT,SNUP,
     &           ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &           RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &           ICE,RTDIS,QUARTZ,FXEXP,FX,CSOIL,BETA,DRIP,
     &           DEW,FLX1,FLX2,FLX3,ESNOW,EMISS,EMISSNOW,MODEL_TYPE,
     &           FRZST,FRZPAR,SACST,SACPAR,ERFLAG,PRFLAG,CELLID,MSTEP,
     &           LVRAIN,wcrit)
         ELSE
            WRITE(*,*)'ABOVE TR_19'
            CALL TR_19 (ETP,ETA,PRCP,PRCP1,SNOWNG,FRZGRA,SMC,SMCMAX,
     &           SMCWLT,SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &           SBETA,DF1,PACH20,CH,SFCSPD,
     &           Q2,T1,TPACK,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &           SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,SNEQV,SNDENS,
     &           SNOWH,SH2O,SLOPE,KDT,FRZX,PSISAT,SNUP,
     &           ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &           RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &           ICE,RTDIS,QUARTZ,FXEXP,FX,CSOIL,UZTWM,UZFWM,UZK,
     &           PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,LZFSM,LZFPM,LZSK,
     &           LZPK,PFREE,SIDE,RSERV,UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,
     &           ADIMC,ERFLAG,BETA,DRIP,DEW,FLX1,FLX2,FLX3,ESNOW,
     &           EMISS,EMISSNOW,ZLVL,SOLDN,LWDN,SNOALB,SN_NEW,
C     State and output variables for PEMB model
     &           NSNOW,DSNOW,PSNOW,RTTSNOW,RTTDTSNOW,WTSNOW,WTDTSNOW,
     &           TTSNOW,TTDTSNOW)
            WRITE(*,*)'BELOW TR_19'
         ENDIF
C     ESNOW = ETA
      ENDIF

C ----------------------------------------------------------------------
C UPDATE WATER BALANCE, ADD ANY DEW/FROST FALL TO PRECIP. CHANGE IT 
C FROM M/S TO KG M-2 S-1 AS IT WAS PASSED FROM DRIVER
C ----------------------------------------------------------------------

c      PRCP = PRCP + (DEW * 1000)
      smctot=0.
      sh2otot=0.
      runtot=0.
      runtot=(runoff1+runoff2)*-1
      do i = 1,nsoil
         smctot = smctot + smc(i)*sldpth(i)
         sh2otot = sh2otot + sh2o(i)*sldpth(i)
      enddo  
      smcdif=smctot0-smctot
      sh2odif=sh2otot0-sh2otot
      if (prflag==2) then
      write(*,*)'-----------sflx BELOW NOPAC------' 
      write(*,*)'prcp           edir            ec'
      write(*,*)prcp,edir,ec
      write(*,*)'ett            eta             etp'
      write(*,*)ett,eta,etp
      write(*,*)'runoff          swe'
      write(*,*)runtot,sneqv
      write(*,*)'smcdif         sh2odif         cmc'
      write(*,*)smcdif,sh2odif,cmc
      write(*,*)'--------------------------------' 
      endif
C ----------------------------------------------------------------------
C   PREPARE SENSIBLE HEAT (H) FOR RETURN TO PARENT MODEL
C ----------------------------------------------------------------------
      SHEAT = -(CH * CP * SFCPRS)/(R * T2V) * ( TH2 - T1 )
          
C ----------------------------------------------------------------------
C  CONVERT UNITS AND/OR SIGN OF TOTAL EVAP (ETA), POTENTIAL EVAP (ETP),
C  SUBSurfACE HEAT FLUX (S), AND RUNOFFS FOR WHAT PARENT MODEL EXPECTS
C  CONVERT ETA FROM KG M-2 S-1 TO W M-2
C ----------------------------------------------------------------------

CBL Moved conversion from M/S to KG M-2 S-1 to this point to account for
CBL both snow covered and snow-free cases at once (e.g. * 1000)
      EDIR = (EDIR * 1000) * LVH2O
      EC = (EC * 1000) * LVH2O
      DO K=1,4
        ET(K) = (ET(K) * 1000) * LVH2O
      ENDDO
      ETT = (ETT * 1000) * LVH2O
CBL multiply by the 'L' that was used to obtain ESNOW
      IF (LVRAIN .EQ. .TRUE.) THEN
         ESNOW = (ESNOW * 1000) * LVH2O
      ELSE
         ESNOW = (ESNOW * 1000) * LSUBS
      ENDIF
CBL      ETP = ETP*((1.-SNCOVR)*LVH2O + SNCOVR*LSUBS)
      ETP = ETP*LVH2O
      IF (ETP .GT. 0.) THEN
        ETA = EDIR + EC + ETT + ESNOW
      ELSE
        ETA = ETP * (1 - SNCOVR) + ESNOW
      ENDIF
      BETA = ETA/ETP
      if (prflag==1) then
      write(*,*)'sflxlv edir ec ett esnow eta etp'
      write(*,*)edir,ec,ett,esnow,eta,etp
      endif
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C CONVERT THE SIGN OF SOIL HEAT FLUX SO THAT:
C   SSOIL>0: WARM THE SURFACE  (NIGHT TIME)
C   SSOIL<0: COOL THE SURFACE  (DAY TIME)
C ----------------------------------------------------------------------
      SSOIL = -1.0*SSOIL      

C ----------------------------------------------------------------------
C  CONVERT RUNOFF3 (INTERNAL LAYER RUNOFF FROM SUPERSAT) FROM M TO M S-1
C  AND ADD TO SUBSURFACE RUNOFF/DRAINAGE/BASEFLOW
C ----------------------------------------------------------------------
      RUNOFF3 = RUNOFF3/DT
      RUNOFF2 = RUNOFF2+RUNOFF3

C ----------------------------------------------------------------------
C TOTAL COLUMN SOIL MOISTURE IN METERS (SOILM) AND ROOT-ZONE 
C SOIL MOISTURE AVAILABILITY (FRACTION) RELATIVE TO POROSITY/SATURATION
C ----------------------------------------------------------------------
      SOILM = -1.0*SMC(1)*ZSOIL(1)
      DO K = 2,NSOIL
        SOILM = SOILM+SMC(K)*(ZSOIL(K-1)-ZSOIL(K))
      END DO
      SOILWM = -1.0*(SMCMAX-SMCWLT)*ZSOIL(1)
      SOILWW = -1.0*(SMC(1)-SMCWLT)*ZSOIL(1)
      DO K = 2,NROOT
        SOILWM = SOILWM+(SMCMAX-SMCWLT)*(ZSOIL(K-1)-ZSOIL(K))
        SOILWW = SOILWW+(SMC(K)-SMCWLT)*(ZSOIL(K-1)-ZSOIL(K))
      END DO
      SOILW = SOILWW/SOILWM


C ----------------------------------------------------------------------
C END SUBROUTINE SFLX
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE ALCALC (ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1)
C   ALB     SNOWFREE ALBEDO
C   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO
C   SHDFAC    AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SHDMIN    MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SNCOVR  FRACTIONAL SNOW COVER
C   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT
C   TSNOW   SNOW SURFACE TEMPERATURE (K)
C ----------------------------------------------------------------------
      REAL ALB, SNOALB, SHDFAC, SHDMIN, SNCOVR, ALBEDO, TSNOW
      REAL TM
      
C ----------------------------------------------------------------------
C SNOALB IS ARGUMENT REPRESENTING MAXIMUM ALBEDO OVER DEEP SNOW,
C AS PASSED INTO SFLX, AND ADAPTED FROM THE SATELLITE-BASED MAXIMUM 
C SNOW ALBEDO FIELDS PROVIDED BY D. ROBINSON AND G. KUKLA 
C (1985, JCAM, VOL 24, 402-411)
C ----------------------------------------------------------------------
C         changed in version 2.6 on June 2nd 2003
C          ALBEDO = ALB + (1.0-(SHDFAC-SHDMIN))*SNCOVR*(SNOALB-ALB) 
          ALBEDO = ALB + SNCOVR*(SNOALB-ALB)
          IF (ALBEDO .GT. SNOALB) ALBEDO=SNOALB

C     BASE FORMULATION (DICKINSON ET AL., 1986, COGLEY ET AL., 1990)
C          IF (TSNOW.LE.263.16) THEN
C            ALBEDO=SNOALB
C          ELSE
C            IF (TSNOW.LT.273.16) THEN
C              TM=0.1*(TSNOW-263.16)
C              ALBEDO=0.5*((0.9-0.2*(TM**3))+(0.8-0.16*(TM**3)))
C            ELSE
C              ALBEDO=0.67
C            ENDIF
C          ENDIF

C     ISBA FORMULATION (VERSEGHY, 1991; BAKER ET AL., 1990)
C          IF (TSNOW.LT.273.16) THEN
C            ALBEDO=SNOALB-0.008*DT/86400
C          ELSE
C            ALBEDO=(SNOALB-0.5)*EXP(-0.24*DT/86400)+0.5
C          ENDIF

C ----------------------------------------------------------------------
C END SUBROUTINE ALCALC
C ----------------------------------------------------------------------
      RETURN
      END

C ---------------------------------------------------------------------
C BETA CODE TO PROVIDE AN ALBEDO DECAY FUNTION SIMILAR TO VIC MODEL
C Ben Livneh JUNE 2007
C
C Modifications:
C 2008-Sep-20 Added minimum snow albedo MINALB; this is necessary
C             when using the NCEP "compromise" scheme of initial
C             snow albedo = avg of 0.85 and satellite-derived scene
C             albedo.							BL via TJB
C 2008-Sep-20 Changed COEF to 1.0 so that initial snow albedo is 0.85.	TJB
C ---------------------------------------------------------------------

      SUBROUTINE ALBCPU (ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO,
     &                   SNOWNG,TFREEZ,DT,SNEQV,T12,LSTSNW1,VEGTYP)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1)
C   ALB     SNOWFREE ALBEDO
C   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO
C   SHDFAC  AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SHDMIN  MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SNCOVR  FRACTIONAL SNOW COVER
C   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT
C   TSNOW   SNOW SURFACE TEMPERATURE (K)
C   SNOWNG  FLAG INDICATING WHETHER SNOW IS FALLING
C   TFREEZ  FREEZING TEMPERATURE
C   DT      TIME STEP IN SECONDS
C   SNEQV   LIQUID WATER-EQUIVALENT SNOW DEPTH (M)
C   LSTSNW  NUMBER OF TIMESTEPS SINCE LAST SNOWFALL
C   LSTSNW1 INITIAL NUMBER OF TIMESTEPS SINCE LAST SNOWFALL
C   T12     EFFECTIVE SNOW/GROUND TEMP
C ----------------------------------------------------------------------
      REAL ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,ALBEDO,TSNOW
      REAL TFREEZ,DT,SNACCA,SNACCB,SNTHWA,SNTHWB,SNEQV,T12
C  ----------------- SNOALB2 adjusted maximum albedo
C  ----------------- COEF is adjustable parameter --------------------
      REAL SNOALB1,SNOALB2,COEF 
      REAL MINALB

      INTEGER VEGTYP	
      INTEGER LSTSNW1,LSTSNW
      LOGICAL SNOWNG

C      PARAMETER (COEF=0.5)
      PARAMETER (COEF=1.0)
      PARAMETER
     &     (SNACCA=0.94,SNACCB=0.58,SNTHWA=0.82,SNTHWB=0.46)
      PARAMETER (MINALB=0.4)

C ----------------------------------------------------------------------
C SNOALB IS CONSIDERED AS THE MAXIMUM SNOW ALBEDO FOR NEW SNOW, AT 
C A VALUE OF 85%. SNOW ALBEDO CURVE DEFAULTS ARE FROM BRAS P.263. SHOULD 
C NOT BE CHANGED EXCEPT FOR SERIOUS PROBLEMS WITH SNOW MELT.
C TO IMPLEMENT ACCUMULATIN PARAMETERS, SNACCA AND SNACCB, ASSERT THAT IT 
C IS INDEED ACCUMULATION SEASON. I.E. THAT SNOW SURFACE TEMP IS BELOW
C ZERO AND THE DATE FALLS BETWEEN OCTOBER AND FEBRUARY 
C ----------------------------------------------------------------------
	   SNOALB1 = SNOALB+COEF*(0.85-SNOALB)
C           IF(VEGTYP.EQ.1.OR.VEGTYP.EQ.5) SNOALB1=0.85 
	   SNOALB2=SNOALB1	   
C ---------------- Initial LSTSNW --------------------------------------
		LSTSNW=LSTSNW1
          IF (SNOWNG) THEN
             SNOALB2=SNOALB1
             LSTSNW=0
          ELSE
            LSTSNW=LSTSNW+1
            IF (SNEQV.GT.0.0) THEN
              IF (T12.LT.TFREEZ) THEN
               SNOALB2=SNOALB1*(SNACCA**((LSTSNW*DT/86400.0)**SNACCB))
              ELSE
              SNOALB2 =SNOALB1*(SNTHWA**((LSTSNW*DT/86400.0)**SNTHWB))
              ENDIF
            ENDIF
          ENDIF

	  IF (SNOALB2 .LT. MINALB) THEN
	    SNOALB2 = MINALB
	  ENDIF
  
	   ALBEDO = ALB + SNCOVR*(SNOALB2-ALB)
           IF (ALBEDO .GT. SNOALB2) ALBEDO=SNOALB2

	  LSTSNW1=LSTSNW
C ----------------------------------------------------------------------
C END SUBROUTINE ALBCPU
C ----------------------------------------------------------------------
      RETURN
      END



      SUBROUTINE ALBCPU27 (ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO,
     &                   LSTSNW,SNOWNG,TFREEZ,DT,SNEQV,T12)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1)
C   ALB     SNOWFREE ALBEDO
C   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO
C   SHDFAC  AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SHDMIN  MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SNCOVR  FRACTIONAL SNOW COVER
C   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT
C   TSNOW   SNOW SURFACE TEMPERATURE (K)
C   SNOWNG  FLAG INDICATING WHETHER SNOW IS FALLING
C   TFREEZ  FREEZING TEMPERATURE
C   DT      TIME STEP IN SECONDS
C   SNEQV   LIQUID WATER-EQUIVALENT SNOW DEPTH (M)
C   LSTSNW  NUMBER OF TIMESTEPS SINCE LAST SNOWFALL
C   T12     EFFECTIVE SNOW/GROUND TEMP
C ----------------------------------------------------------------------
      REAL ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,ALBEDO,TSNOW
      REAL TFREEZ,DT,SNACCA,SNACCB,SNTHWA,SNTHWB,SNEQV,T12

      INTEGER LSTSNW

      LOGICAL SNOWNG

      PARAMETER (SNACCA=0.94,SNACCB=0.58,SNTHWA=0.82,SNTHWB=0.46)

C ----------------------------------------------------------------------
C SNOALB IS CONSIDERED AS THE MAXIMUM SNOW ALBEDO FOR NEW SNOW, AT
C A VALUE OF 85%. SNOW ALBEDO CURVE DEFAULTS ARE FROM BRAS P.263. SHOULD
C NOT BE CHANGED EXCEPT FOR SERIOUS PROBLEMS WITH SNOW MELT.
C Parameters appropriate for accumulation or ablation are chosen based on
C Whether T12 is above freezing.
C ----------------------------------------------------------------------
      SNOALB = 0.85
      IF (SNOWNG) THEN
        ALBEDO=SNOALB
        LSTSNW=0
      ELSE
        LSTSNW=LSTSNW+1
        IF (SNEQV.GT.0.0) THEN
          IF (T12.LE.TFREEZ) THEN
            ALBEDO=ALB*(1-SNCOVR)+SNCOVR*SNOALB*(SNACCA**((LSTSNW*DT/86400)**SNACCB))
          ELSE
            ALBEDO=ALB*(1-SNCOVR)+SNCOVR*SNOALB*(SNTHWA**((LSTSNW*DT/86400)**SNTHWB))
          ENDIF
          IF (ALBEDO .GT. SNOALB) ALBEDO=SNOALB
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C END SUBROUTINE ALBCPU27
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE CANRES (SOLAR,CH,SFCTMP,Q2,SFCPRS,SMC,ZSOIL,NSOIL,
     &                   SMCWLT,SMCREF,RSMIN,RC,PC,NROOT,Q2SAT,DQSDT2, 
     &                   TOPT,RSMAX,RGL,HS,XLAI,STC,
     &                   RCS,RCT,RCQ,RCSOIL,EMISS,EMISSNOW,SNCOVR,RTDIS)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE CANRES                    
C ----------------------------------------------------------------------
C CALCULATE CANOPY RESISTANCE WHICH DEPENDS ON INCOMING SOLAR RADIATION,
C AIR TEMPERATURE, ATMOSPHERIC WATER VAPOR PRESSURE DEFICIT AT THE
C LOWEST MODEL LEVEL, AND SOIL MOISTURE (PREFERABLY UNFROZEN SOIL
C MOISTURE RATHER THAN TOTAL)
C ----------------------------------------------------------------------
C SOURCE:  JARVIS (1976), NOILHAN AND PLANTON (1989, MWR), JACQUEMIN AND
C NOILHAN (1990, BLM)
C SEE ALSO:  CHEN ET AL (1996, JGR, VOL 101(D3), 7251-7268), EQNS 12-14
C AND TABLE 2 OF SEC. 3.1.2         
C ----------------------------------------------------------------------
C INPUT:
C   SOLAR   INCOMING SOLAR RADIATION
C   CH      SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
C   SFCTMP  AIR TEMPERATURE AT 1ST LEVEL ABOVE GROUND
C   Q2      AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
C   Q2SAT   SATURATION AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
C   DQSDT2  SLOPE OF SATURATION HUMIDITY FUNCTION WRT TEMP
C   SFCPRS  SURFACE PRESSURE
C   SMC     VOLUMETRIC SOIL MOISTURE 
C   ZSOIL   SOIL DEPTH (NEGATIVE SIGN, AS IT IS BELOW GROUND)
C   NSOIL   NO. OF SOIL LAYERS
C   NROOT   NO. OF SOIL LAYERS IN ROOT ZONE (1.LE.NROOT.LE.NSOIL)
C   XLAI    LEAF AREA INDEX
C   SMCWLT  WILTING POINT
C   SMCREF  REFERENCE SOIL MOISTURE (WHERE SOIL WATER DEFICIT STRESS
C             SETS IN)
C RSMIN, RSMAX, TOPT, RGL, HS ARE CANOPY STRESS PARAMETERS SET IN
C   SURBOUTINE REDPRM
C OUTPUT:
C   PC  PLANT COEFFICIENT
C   RC  CANOPY RESISTANCE
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER K
      INTEGER NROOT
      INTEGER NSOIL

      REAL CH
      REAL CP
      REAL DELTA
      REAL DQSDT2
      REAL FF
      REAL GX
      REAL GTX
      REAL HS
      REAL P
      REAL PART(NSOLD) 
      REAL RTDIS(NSOLD) 
      REAL PC
      REAL Q2
      REAL Q2SAT
      REAL RC
      REAL RSMIN
      REAL RCQ
      REAL RCS
      REAL RCSOIL
      REAL RCT
      REAL RD
      REAL RGL
      REAL RR
      REAL RSMAX
      REAL SFCPRS
      REAL SFCTMP
      REAL SIGMA
      REAL EMISS
      REAL EMISSNOW
      REAL SNCOVR
      REAL SLV
      REAL SLVS
      REAL SMC(NSOIL)
      REAL STC(NSOIL)
      REAL SMCREF
      REAL SMCWLT
      REAL SOLAR
      REAL TOPT
      REAL SLVCP
      REAL ST1
      REAL TAIR4
      REAL XLAI
      REAL ZSOIL(NSOIL)
      REAL RCQMIN
      REAL STCRT
      REAL DSTCRT
      REAL RCQEXP
      REAL RCQFAC

      PARAMETER(CP = 1004.5)
      PARAMETER(RD = 287.04)
      PARAMETER(SIGMA = 5.67E-8)
      PARAMETER(SLV = 2.501000E6)
      PARAMETER(SLVS = 2.83E6)

C ----------------------------------------------------------------------
C INITIALIZE CANOPY RESISTANCE MULTIPLIER TERMS.
C ----------------------------------------------------------------------
      RCS = 0.0
      RCT = 0.0
      RCQ = 0.0
      RCSOIL = 0.0
      RC = 0.0

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO INCOMING SOLAR RADIATION
C ----------------------------------------------------------------------
      FF = 0.55*2.0*SOLAR/(RGL*XLAI)
      RCS = (FF + RSMIN/RSMAX) / (1.0 + FF)
      RCS = MAX(RCS,0.0001)

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO AIR TEMPERATURE AT FIRST MODEL LEVEL ABOVE GROUND
C RCT EXPRESSION FROM NOILHAN AND PLANTON (1989, MWR).
C ----------------------------------------------------------------------
      RCT = 1.0 - 0.0016*((TOPT-SFCTMP)**2.0)
      RCT = MAX(RCT,0.0001)

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO VAPOR PRESSURE DEFICIT AT FIRST MODEL LEVEL.
C RCQ EXPRESSION FROM SSIB 
C ----------------------------------------------------------------------
c     RCQ = 1.0/(1.0+HS*(Q2SAT-Q2))
c     RCQ = MAX(RCQ,0.01)
C RCQ expression from original OSU model (via Noilhan and Planton),
C extended:
C when RCQEXP=1.0, and RCQMIN=0.25, this is the original NP expression
      RCQMIN = 0.25
c     RCQMIN = 0.1
c     RCQEXP = 1.0
      RCQEXP = 0.5
      RCQFAC = RCQMIN**(1./RCQEXP)
      RCQ = (1.0-HS*(Q2SAT-Q2))
      IF (RCQ .GT. RCQFAC) THEN
        RCQ = RCQ**RCQEXP
      ELSE
        RCQ = RCQMIN
      ENDIF
c M.Ek, Sep 2004:
C RCQ expression from original OSU model (via Noilhan and Planton),
C but with RCQMIN=0.0001 (RCQMIN=0.25, original NP expression)
c      RCQMIN = 0.25
c     RCQMIN = 0.0001
c     RCQ = (1.0-HS*(Q2SAT-Q2))
c     IF (RCQ .GT. 1.0) RCQ = 1.0
c     IF (RCQ .LT. RCQMIN) RCQ = RCQMIN

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO SOIL MOISTURE AVAILABILITY.
C DETERMINE CONTRIBUTION FROM EACH SOIL LAYER, THEN ADD THEM UP.
C ----------------------------------------------------------------------
      GX = (SMC(1) - SMCWLT) / (SMCREF - SMCWLT)
C calculate root zone soil temp
       STCRT=0
       DSTCRT=0
      DO K=1,NROOT
       STCRT=STCRT+ZSOIL(K)*STC(K)
       DSTCRT=DSTCRT+ZSOIL(K)
      ENDDO

       STCRT=STCRT/DSTCRT
        
      GTX =  MAX(0.0,1.-0.0016* MAX(298.- STCRT,0.0)**2)
c     print*,"STCRT=",STCRT,"GTX=",GTX
      IF (GX .GT. 1.) GX = 1.
      IF (GX .LT. 0.) GX = 0.

C ----------------------------------------------------------------------
C USE SOIL DEPTH AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
C     PART(1) = (ZSOIL(1)/ZSOIL(NROOT)) * GX
C ----------------------------------------------------------------------
C USE ROOT DISTRIBUTION AS WEIGHTING FACTOR
       PART(1) = RTDIS(1) * GX*GTX
C ----------------------------------------------------------------------
      DO K = 2,NROOT
        GX = (SMC(K) - SMCWLT) / (SMCREF - SMCWLT)
        GTX= MAX(0.0,1.-0.0016* MAX(298.- STCRT,0.0)**2)
        IF (GX .GT. 1.) GX = 1.
        IF (GX .LT. 0.) GX = 0.
C ----------------------------------------------------------------------
C USE SOIL DEPTH AS WEIGHTING FACTOR        
C ----------------------------------------------------------------------
C       PART(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT)) * GX
C ----------------------------------------------------------------------
C USE ROOT DISTRIBUTION AS WEIGHTING FACTOR
        PART(K) = RTDIS(K) * GX * GTX 
C ----------------------------------------------------------------------
      END DO

      DO K = 1,NROOT
        RCSOIL = RCSOIL+PART(K)
      END DO
      RCSOIL = MAX(RCSOIL,0.0001)

C ----------------------------------------------------------------------
C DETERMINE CANOPY RESISTANCE DUE TO ALL FACTORS.  CONVERT CANOPY
C RESISTANCE (RC) TO PLANT COEFFICIENT (PC) TO BE USED WITH POTENTIAL
C EVAP IN DETERMINING ACTUAL EVAP.  PC IS DETERMINED BY:
C   PC * LINERIZED PENMAN POTENTIAL EVAP =
C   PENMAN-MONTEITH ACTUAL EVAPORATION (CONTAINING RC TERM).
C ----------------------------------------------------------------------
      RC = RSMIN/(XLAI*RCS*RCT*RCQ*RCSOIL)

c      TAIR4 = SFCTMP**4.
c      ST1 = (4.*SIGMA*RD)/CP
c      SLVCP = SLV/CP
c      RR = ST1*TAIR4/(SFCPRS*CH) + 1.0
      RR = (4.*(EMISSNOW*SNCOVR+EMISS*(1.0-SNCOVR))*SIGMA*RD/CP)*
     &     (SFCTMP**4.)/(SFCPRS*CH) + 1.0
C      RR = (4.*SIGMA*RD/CP)*(SFCTMP**4.)/(SFCPRS*CH) + 1.0
      DELTA = (((1.0-SNCOVR)*SLV+SNCOVR*SLVS)/CP)*DQSDT2

      PC = (RR+DELTA)/(RR*(1.+RC*CH)+DELTA)

C ----------------------------------------------------------------------
C END SUBROUTINE CANRES
C ----------------------------------------------------------------------
      RETURN
      END

      FUNCTION CSNOW (DSNOW)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C FUNCTION CSNOW
C ----------------------------------------------------------------------
C CALCULATE SNOW TERMAL CONDUCTIVITY
C ----------------------------------------------------------------------
      REAL C
      REAL DSNOW
      REAL CSNOW
      REAL UNIT

      PARAMETER(UNIT = 0.11631) 
                                         
C ----------------------------------------------------------------------
C CSNOW IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C)
C BASIC VERSION IS DYACHKOVA EQUATION (1960), FOR RANGE 0.1-0.4
C ----------------------------------------------------------------------
      C=0.328*10**(2.25*DSNOW)
      CSNOW=UNIT*C

C ----------------------------------------------------------------------
C DE VAUX EQUATION (1933), IN RANGE 0.1-0.6
C ----------------------------------------------------------------------
C      CSNOW=0.0293*(1.+100.*DSNOW**2)
      
C ----------------------------------------------------------------------
C E. ANDERSEN FROM FLERCHINGER
C ----------------------------------------------------------------------
C      CSNOW=0.021+2.51*DSNOW**2        
      
C ----------------------------------------------------------------------
C END FUNCTION CSNOW
C ----------------------------------------------------------------------
      RETURN                                                      
      END
      SUBROUTINE DEVAP (EDIR1,ETP1,SMC,ZSOIL,SHDFAC,SMCMAX,BEXP,
c      FUNCTION DEVAP (ETP1,SMC,ZSOIL,SHDFAC,SMCMAX,BEXP,
     &                DKSAT,DWSAT,SMCDRY,SMCREF,SMCWLT,FXEXP,FX)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE DEVAP
C FUNCTION DEVAP
C ----------------------------------------------------------------------
C CALCULATE DIRECT SOIL EVAPORATION
C ----------------------------------------------------------------------
      REAL BEXP
c      REAL DEVAP
      REAL EDIR1
      REAL DKSAT
      REAL DWSAT
      REAL ETP1
      REAL FX
      REAL FXEXP
      REAL SHDFAC
      REAL SMC
      REAL SMCDRY
      REAL SMCMAX
      REAL ZSOIL
      REAL SMCREF
      REAL SMCWLT
      REAL SRATIO

C ----------------------------------------------------------------------
C DIRECT EVAP A FUNCTION OF RELATIVE SOIL MOISTURE AVAILABILITY, LINEAR
C WHEN FXEXP=1.
C FX > 1 REPRESENTS DEMAND CONTROL
C FX < 1 REPRESENTS FLUX CONTROL
C ----------------------------------------------------------------------
      SRATIO = (SMC - SMCDRY) / (SMCMAX - SMCDRY)
      IF (SRATIO .GT. 0.) THEN
        FX = SRATIO**FXEXP
        FX = MAX ( MIN ( FX, 1. ) ,0. )
      ELSE
        FX = 0.
      ENDIF

C ----------------------------------------------------------------------
C ALLOW FOR THE DIRECT-EVAP-REDUCING EFFECT OF SHADE
C ----------------------------------------------------------------------
c      DEVAP = FX * ( 1.0 - SHDFAC ) * ETP1
      EDIR1 = FX * ( 1.0 - SHDFAC ) * ETP1

C ----------------------------------------------------------------------
C END SUBROUTINE DEVAP
C END FUNCTION DEVAP
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE EVAPO (ETA1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
     &                  SH2O,STC,
     &                  SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &                  SMCREF,SHDFAC,CMCMAX,
     &                  SMCDRY,CFACTR,
     &                  EDIR1,EC1,ET1,ETT1,SFCTMP,Q2,NROOT,RTDIS,
     &                  FXEXP,FX,MODEL_TYPE,SMCFLAG)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE EVAPO
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER
C UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH
C PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
C FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND
C CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K
      INTEGER NSOIL
      INTEGER NROOT

      REAL BEXP
      REAL CFACTR
      REAL CMC
      REAL CMC2MS
      REAL CMCMAX
c      REAL DEVAP
      REAL DKSAT
      REAL DT
      REAL DWSAT
      REAL EC1
      REAL EDIR1
      REAL ET1(NSOIL)
      REAL ETA1
      REAL ETP1
      REAL ETT1
      REAL FXEXP
      REAL FX
      REAL PC
      REAL Q2
      REAL RTDIS(NSOIL)
      REAL SFCTMP
      REAL SHDFAC
      REAL SMC(NSOIL)
      REAL STC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL ZSOIL(NSOIL)
      INTEGER MODEL_TYPE,SMCFLAG

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE IF THE POTENTIAL EVAPOTRANSPIRATION IS
C GREATER THAN ZERO.
C ----------------------------------------------------------------------

      EC1 = 0.
      DO K = 1,NSOIL
        ET1(K) = 0.
      END DO
      ETT1 = 0.

      IF (ETP1 .GT. 0.0) THEN

C ----------------------------------------------------------------------
C RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE.  CALL THIS FUNCTION
C ONLY IF VEG COVER NOT COMPLETE.
C FROZEN GROUND VERSION:  SH2O STATES REPLACE SMC STATES.
C ----------------------------------------------------------------------
         IF (SHDFAC .LT. 1.) THEN
C Only compute this EDIR if using Noah SM scheme
            EDIR1 = 0.
           IF (MODEL_TYPE.EQ.0) THEN
               CALL DEVAP (EDIR1,ETP1,SH2O(1),ZSOIL(1),SHDFAC,SMCMAX,
c     EDIR = DEVAP(ETP1,SH2O(1),ZSOIL(1),SHDFAC,SMCMAX,
     &              BEXP,DKSAT,DWSAT,SMCDRY,SMCREF,SMCWLT,FXEXP,FX)
            ENDIF
         ENDIF

C ----------------------------------------------------------------------
C INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION,
C AND ACCUMULATE IT FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
        IF (SHDFAC.GT.0.0) THEN

         CALL TRANSP(ET1,NSOIL,ETP1,SH2O,STC,CMC,ZSOIL,SHDFAC,SMCWLT,
     &                 CMCMAX,PC,CFACTR,SMCREF,SFCTMP,Q2,NROOT,RTDIS)

          DO K = 1,NSOIL
            ETT1 = ETT1 + ET1(K)
          END DO

C ----------------------------------------------------------------------
C CALCULATE CANOPY EVAPORATION.
C IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR CMC=0.0.
C ----------------------------------------------------------------------
          IF (CMC .GT. 0.0) THEN
            EC1 = SHDFAC * ( ( CMC / CMCMAX ) ** CFACTR ) * ETP1
          ELSE
            EC1 = 0.0
          ENDIF

C ----------------------------------------------------------------------
C EC SHOULD BE LIMITED BY THE TOTAL AMOUNT OF AVAILABLE WATER ON THE
C CANOPY.  -F.CHEN, 18-OCT-1994
C ----------------------------------------------------------------------
          CMC2MS = CMC / DT
          EC1 = MIN ( CMC2MS, EC1 )
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C TOTAL UP EVAP AND TRANSP TYPES TO OBTAIN ACTUAL EVAPOTRANSP
C ----------------------------------------------------------------------
      ETA1 = EDIR1 + ETT1 + EC1

C ----------------------------------------------------------------------
C END SUBROUTINE EVAPO
C ----------------------------------------------------------------------
      RETURN
      END

      FUNCTION FRH2O (TKELV,SMC,SH2O,SMCMAX,BEXP,PSIS)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C FUNCTION FRH2O
C ----------------------------------------------------------------------
C CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT IF
C TEMPERATURE IS BELOW 273.15K (T0).  REQUIRES NEWTON-TYPE ITERATION TO
C SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF KOREN ET AL
C (1999, JGR, VOL 104(D16), 19569-19585).
C ----------------------------------------------------------------------
C NEW VERSION (JUNE 2001): MUCH FASTER AND MORE ACCURATE NEWTON
C ITERATION ACHIEVED BY FIRST TAKING LOG OF EQN CITED ABOVE -- LESS THAN
C 4 (TYPICALLY 1 OR 2) ITERATIONS ACHIEVES CONVERGENCE.  ALSO, EXPLICIT
C 1-STEP SOLUTION OPTION FOR SPECIAL CASE OF PARAMETER CK=0, WHICH
C REDUCES THE ORIGINAL IMPLICIT EQUATION TO A SIMPLER EXPLICIT FORM,
C KNOWN AS THE "FLERCHINGER EQN". IMPROVED HANDLING OF SOLUTION IN THE
C LIMIT OF FREEZING POINT TEMPERATURE T0.
C ----------------------------------------------------------------------
C INPUT:
C
C   TKELV.........TEMPERATURE (Kelvin)
C   SMC...........TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC)
C   SH2O..........LIQUID SOIL MOISTURE CONTENT (VOLUMETRIC)
C   SMCMAX........SATURATION SOIL MOISTURE CONTENT (FROM REDPRM)
C   B.............SOIL TYPE "B" PARAMETER (FROM REDPRM)
C   PSIS..........SATURATED SOIL MATRIC POTENTIAL (FROM REDPRM)
C
C OUTPUT:
C   FRH2O.........SUPERCOOLED LIQUID WATER CONTENT
C ----------------------------------------------------------------------
      REAL BEXP
      REAL BLIM
      REAL BX
      REAL CK
      REAL DENOM
      REAL DF
      REAL DH2O
      REAL DICE
      REAL DSWL
      REAL ERROR
      REAL FK
      REAL FRH2O
      REAL GS
      REAL HLICE
      REAL PSIS
      REAL SH2O
      REAL SMC
      REAL SMCMAX
      REAL SWL
      REAL SWLK
      REAL TKELV
      REAL T0

      INTEGER NLOG
      INTEGER KCOUNT

      PARAMETER(CK = 8.0)
C      PARAMETER(CK = 0.0)
      PARAMETER(BLIM = 5.5)
      PARAMETER(ERROR = 0.005)

      PARAMETER(HLICE = 3.335E5)
      PARAMETER(GS = 9.81)
      PARAMETER(DICE = 920.0)
      PARAMETER(DH2O = 1000.0)
      PARAMETER(T0 = 273.15)

C ----------------------------------------------------------------------
C LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)
C SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT IS
C NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES.
C ----------------------------------------------------------------------
      BX = BEXP
      IF (BEXP .GT. BLIM) BX = BLIM

C ----------------------------------------------------------------------
C INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
C ----------------------------------------------------------------------
      NLOG=0
      KCOUNT=0

C ----------------------------------------------------------------------
C  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC
C ----------------------------------------------------------------------
      IF (TKELV .GT. (T0 - 1.E-3)) THEN
 	FRH2O = SMC
      ELSE
        IF (CK .NE. 0.0) THEN

C ----------------------------------------------------------------------
C OPTION 1: ITERATED SOLUTION FOR NONZERO CK
C IN KOREN ET AL, JGR, 1999, EQN 17
C ----------------------------------------------------------------------
C INITIAL GUESS FOR SWL (frozen content)
C ----------------------------------------------------------------------
          SWL = SMC-SH2O

C ----------------------------------------------------------------------
C KEEP WITHIN BOUNDS.
C ----------------------------------------------------------------------
          IF (SWL .GT. (SMC-0.02)) SWL = SMC-0.02
          IF (SWL .LT. 0.) SWL = 0.

C ----------------------------------------------------------------------
C  START OF ITERATIONS
C ----------------------------------------------------------------------
          DO WHILE ( (NLOG .LT. 10) .AND. (KCOUNT .EQ. 0) )
            NLOG = NLOG+1
            DF = ALOG(( PSIS*GS/HLICE ) * ( ( 1.+CK*SWL )**2. ) *
     &        ( SMCMAX/(SMC-SWL) )**BX) - ALOG(-(TKELV-T0)/TKELV)
            DENOM = 2. * CK / ( 1.+CK*SWL ) + BX / ( SMC - SWL )
            SWLK = SWL - DF/DENOM
C ----------------------------------------------------------------------
C BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
C ----------------------------------------------------------------------
            IF (SWLK .GT. (SMC-0.02)) SWLK = SMC - 0.02
            IF (SWLK .LT. 0.) SWLK = 0.

C ----------------------------------------------------------------------
C MATHEMATICAL SOLUTION BOUNDS APPLIED.
C ----------------------------------------------------------------------
            DSWL = ABS(SWLK-SWL)
            SWL = SWLK

C ----------------------------------------------------------------------
C IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
C WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
C ----------------------------------------------------------------------
            IF ( DSWL .LE. ERROR )  THEN
 	      KCOUNT = KCOUNT+1
            ENDIF
          END DO

C ----------------------------------------------------------------------
C  END OF ITERATIONS
C ----------------------------------------------------------------------
C BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
C ----------------------------------------------------------------------
          FRH2O = SMC - SWL

C ----------------------------------------------------------------------
C END OPTION 1
C ----------------------------------------------------------------------
        ENDIF

C ----------------------------------------------------------------------
C OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
C IN KOREN ET AL., JGR, 1999, EQN 17
C APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
C ----------------------------------------------------------------------
        IF (KCOUNT .EQ. 0) THEN
          Print*,'Flerchinger used in NEW version. Iterations=',NLOG
 	  FK = (((HLICE/(GS*(-PSIS)))*
     &      ((TKELV-T0)/TKELV))**(-1/BX))*SMCMAX
 	  IF (FK .LT. 0.02) FK = 0.02
 	  FRH2O = MIN (FK, SMC)
C ----------------------------------------------------------------------
C END OPTION 2
C ----------------------------------------------------------------------
        ENDIF

      ENDIF

C ----------------------------------------------------------------------
C END FUNCTION FRH2O
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE HRT (RHSTS,STC,SMC,SMCMAX,NSOIL,ZSOIL,YY,ZZ1,
     &                TBOT,ZBOT,PSISAT,SH2O,DT,BEXP,
     &                F1,DF1,QUARTZ,CSOIL,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE HRT
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C THERMAL DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
C COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      LOGICAL ITAVG

      INTEGER I
      INTEGER K
      INTEGER NSOIL

C ----------------------------------------------------------------------
C DECLARE WORK ARRAYS NEEDED IN TRI-DIAGONAL IMPLICIT SOLVER
C ----------------------------------------------------------------------
      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

C ----------------------------------------------------------------------
C DECLARATIONS
C ----------------------------------------------------------------------
      REAL BEXP
      REAL CAIR
      REAL CH2O
      REAL CICE
      REAL CSOIL
      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DF1
      REAL DF1N
      REAL DF1K
      REAL DT
      REAL DTSDZ
      REAL DTSDZ2
      REAL F1
      REAL HCPCT
      REAL PSISAT
      REAL QUARTZ
      REAL QTOT
      REAL RHSTS(NSOIL)
      REAL SSOIL
      REAL SICE
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCMAX
      REAL SNKSRC
      REAL STC(NSOIL)
      REAL T0
      REAL TAVG
      REAL TBK
      REAL TBK1
      REAL TBOT
      REAL ZBOT
      REAL TSNSR
      REAL TSURF
      REAL YY
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      PARAMETER(T0 = 273.15)

C ----------------------------------------------------------------------
C SET SPECIFIC HEAT CAPACITIES OF AIR, WATER, ICE, SOIL MINERAL       
C ----------------------------------------------------------------------
      PARAMETER(CAIR = 1004.0)
      PARAMETER(CH2O = 4.2E6)
      PARAMETER(CICE = 2.106E6)
C NOTE: CSOIL NOW SET IN ROUTINE REDPRM AND PASSED IN
C      PARAMETER(CSOIL = 1.26E6)

C ----------------------------------------------------------------------
C INITIALIZE LOGICAL FOR SOIL LAYER TEMPERATURE AVERAGING.
C ----------------------------------------------------------------------
      ITAVG = .TRUE.
C      ITAVG = .FALSE.

C ----------------------------------------------------------------------
C BEGIN SECTION FOR TOP SOIL LAYER
C ----------------------------------------------------------------------
C CALC THE HEAT CAPACITY OF THE TOP SOIL LAYER
C ----------------------------------------------------------------------
      HCPCT = SH2O(1)*CH2O + (1.0-SMCMAX)*CSOIL + (SMCMAX-SMC(1))*CAIR
     &        + ( SMC(1) - SH2O(1) )*CICE

C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1.0 / ( -0.5 * ZSOIL(2) )
      AI(1) = 0.0
      CI(1) = (DF1 * DDZ) / (ZSOIL(1) * HCPCT)
      BI(1) = -CI(1) + DF1 / (0.5 * ZSOIL(1) * ZSOIL(1)*HCPCT*ZZ1)

C ----------------------------------------------------------------------
C CALCULATE THE VERTICAL SOIL TEMP GRADIENT BTWN THE 1ST AND 2ND SOIL
C LAYERS.  THEN CALCULATE THE SUBSURFACE HEAT FLUX. USE THE TEMP
C GRADIENT AND SUBSFC HEAT FLUX TO CALC "RIGHT-HAND SIDE TENDENCY
C TERMS", OR "RHSTS", FOR TOP SOIL LAYER.
C ----------------------------------------------------------------------
      DTSDZ = (STC(1) - STC(2)) / (-0.5 * ZSOIL(2))
      SSOIL = DF1 * (STC(1) - YY) / (0.5 * ZSOIL(1) * ZZ1)
      RHSTS(1) = (DF1 * DTSDZ - SSOIL) / (ZSOIL(1) * HCPCT)

C ----------------------------------------------------------------------
C NEXT CAPTURE THE VERTICAL DIFFERENCE OF THE HEAT FLUX AT TOP AND
C BOTTOM OF FIRST SOIL LAYER FOR USE IN HEAT FLUX CONSTRAINT APPLIED TO
C POTENTIAL SOIL FREEZING/THAWING IN ROUTINE SNKSRC.
C ----------------------------------------------------------------------
      QTOT = SSOIL - DF1*DTSDZ

C ----------------------------------------------------------------------
C IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):
C SET TEMP "TSURF" AT TOP OF SOIL COLUMN (FOR USE IN FREEZING SOIL
C PHYSICS LATER IN FUNCTION SUBROUTINE SNKSRC).  IF SNOWPACK CONTENT IS
C ZERO, THEN TSURF EXPRESSION BELOW GIVES TSURF = SKIN TEMP.  IF
C SNOWPACK IS NONZERO (HENCE ARGUMENT ZZ1=1), THEN TSURF EXPRESSION
C BELOW YIELDS SOIL COLUMN TOP TEMPERATURE UNDER SNOWPACK.  THEN
C CALCULATE TEMPERATURE AT BOTTOM INTERFACE OF 1ST SOIL LAYER FOR USE
C LATER IN FUNCTION SUBROUTINE SNKSRC
C ----------------------------------------------------------------------
      IF (ITAVG) THEN 
        TSURF = (YY + (ZZ1-1) * STC(1)) / ZZ1
        CALL TBND (STC(1),STC(2),ZSOIL,ZBOT,1,NSOIL,TBK)
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE FROZEN WATER CONTENT IN 1ST SOIL LAYER. 
C ----------------------------------------------------------------------
      SICE = SMC(1) - SH2O(1)

C ----------------------------------------------------------------------
C IF FROZEN WATER PRESENT OR ANY OF LAYER-1 MID-POINT OR BOUNDING
C INTERFACE TEMPERATURES BELOW FREEZING, THEN CALL SNKSRC TO
C COMPUTE HEAT SOURCE/SINK (AND CHANGE IN FROZEN WATER CONTENT)
C DUE TO POSSIBLE SOIL WATER PHASE CHANGE
C ----------------------------------------------------------------------
      IF ( (SICE   .GT. 0.) .OR. (TSURF .LT. T0) .OR.
     &     (STC(1) .LT. T0) .OR. (TBK   .LT. T0) ) THEN

        IF (ITAVG) THEN 
          CALL TMPAVG(TAVG,TSURF,STC(1),TBK,ZSOIL,NSOIL,1)
        ELSE
          TAVG = STC(1)
        ENDIF
        TSNSR = SNKSRC (TAVG,SMC(1),SH2O(1), 
     &    ZSOIL,NSOIL,SMCMAX,PSISAT,BEXP,DT,1,QTOT)

        RHSTS(1) = RHSTS(1) - TSNSR / ( ZSOIL(1) * HCPCT )
      ENDIF
 
C ----------------------------------------------------------------------
C THIS ENDS SECTION FOR TOP SOIL LAYER.
C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0

C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS
C (EXCEPT SUBSFC OR "GROUND" HEAT FLUX NOT REPEATED IN LOWER LAYERS)
C ----------------------------------------------------------------------
      DF1K = DF1
      DO K = 2,NSOIL

C ----------------------------------------------------------------------
C CALCULATE HEAT CAPACITY FOR THIS SOIL LAYER.
C ----------------------------------------------------------------------
        HCPCT = SH2O(K)*CH2O +(1.0-SMCMAX)*CSOIL +(SMCMAX-SMC(K))*CAIR
     &        + ( SMC(K) - SH2O(K) )*CICE

        IF (K .NE. NSOIL) THEN
C ----------------------------------------------------------------------
C THIS SECTION FOR LAYER 2 OR GREATER, BUT NOT LAST LAYER.
C ----------------------------------------------------------------------
C CALCULATE THERMAL DIFFUSIVITY FOR THIS LAYER.
C ----------------------------------------------------------------------
          CALL TDFCND (DF1N,SMC(K),QUARTZ,SMCMAX,SH2O(K))

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER
C ----------------------------------------------------------------------
          DENOM = 0.5 * ( ZSOIL(K-1) - ZSOIL(K+1) )
          DTSDZ2 = ( STC(K) - STC(K+1) ) / DENOM

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
C ----------------------------------------------------------------------
          DDZ2 = 2. / (ZSOIL(K-1) - ZSOIL(K+1))
          CI(K) = -DF1N * DDZ2 / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)

C ----------------------------------------------------------------------
C IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE
C TEMP AT BOTTOM OF LAYER.
C ----------------------------------------------------------------------
          IF (ITAVG) THEN 
            CALL TBND (STC(K),STC(K+1),ZSOIL,ZBOT,K,NSOIL,TBK1)
          ENDIF
        ELSE

C ----------------------------------------------------------------------
C SPECIAL CASE OF BOTTOM SOIL LAYER:  CALCULATE THERMAL DIFFUSIVITY FOR
C BOTTOM LAYER.
C ----------------------------------------------------------------------
          CALL TDFCND (DF1N,SMC(K),QUARTZ,SMCMAX,SH2O(K))

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU BOTTOM LAYER.
C ----------------------------------------------------------------------
          DENOM = .5 * (ZSOIL(K-1) + ZSOIL(K)) - ZBOT
          DTSDZ2 = (STC(K)-TBOT) / DENOM

C ----------------------------------------------------------------------
C SET MATRIX COEF, CI TO ZERO IF BOTTOM LAYER.
C ----------------------------------------------------------------------
          CI(K) = 0.

C ----------------------------------------------------------------------
C IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE
C TEMP AT BOTTOM OF LAST LAYER.
C ----------------------------------------------------------------------
          IF (ITAVG) THEN 
            CALL TBND (STC(K),TBOT,ZSOIL,ZBOT,K,NSOIL,TBK1)
          ENDIF 

        ENDIF
C ----------------------------------------------------------------------
C THIS ENDS SPECIAL LOOP FOR BOTTOM LAYER.
C ----------------------------------------------------------------------
C CALCULATE RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT.
C ----------------------------------------------------------------------
        DENOM = ( ZSOIL(K) - ZSOIL(K-1) ) * HCPCT
        RHSTS(K) = ( DF1N * DTSDZ2 - DF1K * DTSDZ ) / DENOM
        QTOT = -1.0*DENOM*RHSTS(K)
        SICE = SMC(K) - SH2O(K)

        IF ( (SICE .GT. 0.) .OR. (TBK .LT. T0) .OR.
     &     (STC(K) .LT. T0) .OR. (TBK1 .LT. T0) ) THEN

          IF (ITAVG) THEN 
            CALL TMPAVG(TAVG,TBK,STC(K),TBK1,ZSOIL,NSOIL,K)
          ELSE
            TAVG = STC(K)
          ENDIF
          TSNSR = SNKSRC(TAVG,SMC(K),SH2O(K),ZSOIL,NSOIL,
     &                   SMCMAX,PSISAT,BEXP,DT,K,QTOT)
          RHSTS(K) = RHSTS(K) - TSNSR / DENOM
        ENDIF 

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
C ----------------------------------------------------------------------
        AI(K) = - DF1K * DDZ / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        BI(K) = -(AI(K) + CI(K))

C ----------------------------------------------------------------------
C RESET VALUES OF DF1, DTSDZ, DDZ, AND TBK FOR LOOP TO NEXT SOIL LAYER.
C ----------------------------------------------------------------------
        TBK   = TBK1
        DF1K  = DF1N
        DTSDZ = DTSDZ2
        DDZ   = DDZ2
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE HRT
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE HRTICE (RHSTS,STC,NSOIL,ZSOIL,YY,ZZ1,DF1,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE HRTICE
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C THERMAL DIFFUSION EQUATION IN THE CASE OF SEA-ICE PACK.  ALSO TO
C COMPUTE (PREPARE) THE MATRIX COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX
C OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER K
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DF1
      REAL DTSDZ
      REAL DTSDZ2
      REAL HCPCT
      REAL RHSTS(NSOIL)
      REAL SSOIL
      REAL STC(NSOIL)
      REAL TBOT
      REAL YY
      REAL ZBOT
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      DATA TBOT /271.16/

C ----------------------------------------------------------------------
C SET A NOMINAL UNIVERSAL VALUE OF THE SEA-ICE SPECIFIC HEAT CAPACITY,
C HCPCT = 1880.0*917.0.
C ----------------------------------------------------------------------
      PARAMETER(HCPCT = 1.72396E+6)

C ----------------------------------------------------------------------
C THE INPUT ARGUMENT DF1 IS A UNIVERSALLY CONSTANT VALUE OF SEA-ICE
C THERMAL DIFFUSIVITY, SET IN ROUTINE SNOPAC AS DF1 = 2.2.
C ----------------------------------------------------------------------
C SET ICE PACK DEPTH.  USE TBOT AS ICE PACK LOWER BOUNDARY TEMPERATURE
C (THAT OF UNFROZEN SEA WATER AT BOTTOM OF SEA ICE PACK).  ASSUME ICE
C PACK IS OF N=NSOIL LAYERS SPANNING A UNIFORM CONSTANT ICE PACK
C THICKNESS AS DEFINED BY ZSOIL(NSOIL) IN ROUTINE SFLX.
C ----------------------------------------------------------------------
      ZBOT = ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1.0 / ( -0.5 * ZSOIL(2) )
      AI(1) = 0.0
      CI(1) = (DF1 * DDZ) / (ZSOIL(1) * HCPCT)
      BI(1) = -CI(1) + DF1/(0.5 * ZSOIL(1) * ZSOIL(1) * HCPCT * ZZ1)

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT BTWN THE TOP AND 2ND SOIL LAYERS.
C RECALC/ADJUST THE SOIL HEAT FLUX.  USE THE GRADIENT AND FLUX TO CALC
C RHSTS FOR THE TOP SOIL LAYER.
C ----------------------------------------------------------------------
      DTSDZ = ( STC(1) - STC(2) ) / ( -0.5 * ZSOIL(2) )
      SSOIL = DF1 * ( STC(1) - YY ) / ( 0.5 * ZSOIL(1) * ZZ1 )
      RHSTS(1) = ( DF1 * DTSDZ - SSOIL ) / ( ZSOIL(1) * HCPCT )

C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0

C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        IF (K .NE. NSOIL) THEN

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER.
C ----------------------------------------------------------------------
          DENOM = 0.5 * ( ZSOIL(K-1) - ZSOIL(K+1) )
          DTSDZ2 = ( STC(K) - STC(K+1) ) / DENOM

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT.
C ----------------------------------------------------------------------
          DDZ2 = 2. / (ZSOIL(K-1) - ZSOIL(K+1))
          CI(K) = -DF1 * DDZ2 / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        ELSE

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU THE LOWEST LAYER.
C ----------------------------------------------------------------------
          DTSDZ2 = (STC(K)-TBOT)/(.5 * (ZSOIL(K-1) + ZSOIL(K))-ZBOT)

C ----------------------------------------------------------------------
C SET MATRIX COEF, CI TO ZERO.
C ----------------------------------------------------------------------
          CI(K) = 0.
        ENDIF

C ----------------------------------------------------------------------
C CALC RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT.
C ----------------------------------------------------------------------
        DENOM = ( ZSOIL(K) - ZSOIL(K-1) ) * HCPCT
        RHSTS(K) = ( DF1 * DTSDZ2 - DF1 * DTSDZ ) / DENOM

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
C ----------------------------------------------------------------------
        AI(K) = - DF1 * DDZ / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        BI(K) = -(AI(K) + CI(K))

C ----------------------------------------------------------------------
C RESET VALUES OF DTSDZ AND DDZ FOR LOOP TO NEXT SOIL LYR.
C ----------------------------------------------------------------------
        DTSDZ = DTSDZ2
        DDZ   = DDZ2

      END DO
C ----------------------------------------------------------------------
C END SUBROUTINE HRTICE
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE HSTEP (STCOUT,STCIN,RHSTS,DT,NSOIL,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE HSTEP
C ----------------------------------------------------------------------
C CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER K
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL CIin(NSOLD)
      REAL DT
      REAL RHSTS(NSOIL)
      REAL RHSTSin(NSOIL)
      REAL STCIN(NSOIL)
      REAL STCOUT(NSOIL)

C ----------------------------------------------------------------------
C CREATE FINITE DIFFERENCE VALUES FOR USE IN ROSR12 ROUTINE
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTS(K) = RHSTS(K) * DT
        AI(K) = AI(K) * DT
        BI(K) = 1. + BI(K) * DT
        CI(K) = CI(K) * DT
      END DO

C ----------------------------------------------------------------------
C COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
         RHSTSin(K) = RHSTS(K)
      END DO
      DO K = 1,NSOLD
        CIin(K) = CI(K)
      END DO

C ----------------------------------------------------------------------
C SOLVE THE TRI-DIAGONAL MATRIX EQUATION
C ----------------------------------------------------------------------
      CALL ROSR12(CI,AI,BI,CIin,RHSTSin,RHSTS,NSOIL)

C ----------------------------------------------------------------------
C CALC/UPDATE THE SOIL TEMPS USING MATRIX SOLUTION
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        STCOUT(K) = STCIN(K) + CI(K)
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE HSTEP
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE NOPAC(ETP,ETA,PRCP,SMC,SMCMAX,SMCWLT,
     &     SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,SHDFAC,
     &     SBETA,Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,
     &     STC,EPSCA,BEXP,PC,RCH,RR,CFACTR, 
     &     SH2O,SLOPE,KDT,FRZFACT,PSISAT,ZSOIL,
     &     DKSAT,DWSAT,TBOT,ZBOT,RUNOFF1,RUNOFF2,RUNOFF3,
     &     EDIR,EC,ET,ETT,NROOT,ICE,RTDIS,QUARTZ,FXEXP,FX,
     &     CSOIL,FRZST,FRZPAR,SACST,SACPAR,BETA,DRIP,DEW,
     &     FLX1,FLX2,FLX3,EMISS,MODEL_TYPE,PRFLAG,CELLID,MSTEP,wcrit)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE NOPAC
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES AND UPDATE SOIL MOISTURE
C CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN NO SNOW PACK IS
C PRESENT.
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER NROOT
      INTEGER NSOIL

      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL CP
      REAL CSOIL
      REAL DEW
      REAL DF1
      REAL DKSAT
      REAL DRIP
      REAL DT
      REAL DWSAT
      REAL EC
      REAL EDIR
      REAL EPSCA
      REAL ETA
      REAL ETA1
      REAL ETA0
      REAL ETP
      REAL ETP1
      REAL ET(NSOIL)
      REAL ETT
      REAL FDOWN
      REAL F1
      REAL FXEXP
      REAL FX
      REAL FLX1
      REAL FLX2
      REAL FLX3
      REAL FRZFACT
      REAL KDT
      REAL PC
      REAL PRCP
      REAL PRCP1
      REAL PSISAT
      REAL Q2
      REAL QUARTZ
      REAL RCH
      REAL RR
      REAL RTDIS(NSOIL)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL SSOIL
      REAL SBETA
      REAL SFCTMP
      REAL SHDFAC
      REAL SH2O(NSOIL)
      REAL SIGMA
      REAL EMISS
      REAL SLOPE
      REAL SMC(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL STC(NSOIL)
      REAL T1
      REAL T24
      REAL TBOT
      REAL TH2
      REAL YY
      REAL YYNUM
      REAL ZBOT
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      REAL EC1
      REAL EDIR1
      REAL ET1(NSOIL)
      REAL ETT1
      REAL UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,LZFSM
      REAL LZFPM,LZSK,LZPK,PFREE,SIDE,RSERV,UZTWC,UZFWC,LZTWC,LZFSC
      REAL LZFPC,ADIMC
      REAL FRZST(10)
      REAL FRZPAR(13)
      REAL SACST(6)
      REAL SACPAR(16)
      REAL FROST,wcrit
      INTEGER SMCFLAG

      INTEGER K,PRFLAG,CELLID,MSTEP
      INTEGER MODEL_TYPE

      PARAMETER(CP = 1004.5)
      PARAMETER(SIGMA = 5.67E-8)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C CONVERT ETP FROM KG M-2 S-1 TO MS-1 AND INITIALIZE DEW.
C ----------------------------------------------------------------------
      PRCP1 = PRCP * 0.001
      ETP1 = ETP * 0.001
      DEW = 0.0
      if (prflag==1) then
      write(*,*)'nopac prcp1 prcp',prcp1,prcp
      endif
      EDIR = 0.
      EDIR1 = 0.
      EC = 0.
      EC1 = 0.
      DO K = 1,NSOIL
        ET(K) = 0.
        ET1(K) = 0.
      END DO
      ETT = 0.
      ETT1 = 0.
      IF (ETP .GT. 0.0) THEN

C ----------------------------------------------------------------------
C CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1'.
C ----------------------------------------------------------------------
           CALL EVAPO (ETA1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
     &                 SH2O,STC,SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &                 SMCREF,SHDFAC,CMCMAX,
     &                 SMCDRY,CFACTR, 
     &                 EDIR1,EC1,ET1,ETT1,SFCTMP,Q2,NROOT,RTDIS,
     &                 FXEXP,FX,MODEL_TYPE,SMCFLAG)
           IF (MODEL_TYPE.EQ.0) THEN
           CALL SMFLX (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &                 SH2O,SLOPE,KDT,FRZFACT,
     &                 SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &                 SHDFAC,CMCMAX,
     &                 RUNOFF1,RUNOFF2,RUNOFF3, 
     &                 EDIR1,EC1,ET1,DRIP)
           END IF
C ETP1 NEEDS TO BE CONVERTED TO MM/TIMESTEP FOR SAC, ETA NEEDS TO BE IN
C SAME UNITS AS QS,QG,ETC...
           IF (MODEL_TYPE.EQ.1) THEN
           CALL SMFLX_SAC (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &             SH2O,SLOPE,KDT,FRZFACT,0.0,SMCREF,
     &             SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &             SHDFAC,CMCMAX,SFCTMP,ETP1,ETA1,
     &             RUNOFF1,RUNOFF2,RUNOFF3,EDIR1,EC1,ET1,DRIP, 
     &             FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP,wcrit)
           END IF

CBL Enumerate the intermediate term ETA1 , so as to allow for proper 
CBL units forET efficiency computation (ETA/ETP)and as well, correct
CBL unit conversion in SFLX KG M-2 S-1

           ETA1 = EDIR1 + ETT1 + EC1

C ----------------------------------------------------------------------
C       CONVERT MODELED EVAPOTRANSPIRATION FM  M S-1  TO  KG M-2 S-1
C ----------------------------------------------------------------------
c        ETA = ETA1 * 1000.0

C ----------------------------------------------------------------------
c        EDIR = EDIR1 * 1000.0
c        EC = EC1 * 1000.0
c        ETT = ETT1 * 1000.0
c        ET(1) = ET1(1) * 1000.0
c        ET(2) = ET1(2) * 1000.0
c        ET(3) = ET1(3) * 1000.0
c        ET(4) = ET1(4) * 1000.0
C ----------------------------------------------------------------------

      ELSE

C ----------------------------------------------------------------------
C IF ETP < 0, ASSUME DEW FORMS (TRANSFORM ETP1 INTO DEW AND REINITIALIZE
C ETP1 TO ZERO).
C ----------------------------------------------------------------------
        DEW = -ETP1
        if (prflag==1) then
        write(*,*)'nopac dew',dew
        endif
c        ETP1 = 0.0

C ----------------------------------------------------------------------
C CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1' AND ADD DEW AMOUNT.
C ----------------------------------------------------------------------
        PRCP1 = PRCP1 + DEW
        if (prflag==1) then
        write(*,*)'nopac prcp1 dew',prcp1,dew
        endif
C
c      CALL EVAPO (ETA1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
c     &            SH2O,SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
c     &            SMCREF,SHDFAC,CMCMAX,
c     &            SMCDRY,CFACTR, 
c     &            EDIR1,EC1,ET1,ETT,SFCTMP,Q2,NROOT,RTDIS,FXEXP)
        IF (MODEL_TYPE.EQ.0) THEN
      CALL SMFLX (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &            SH2O,SLOPE,KDT,FRZFACT,
     &            SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &            SHDFAC,CMCMAX,
     &            RUNOFF1,RUNOFF2,RUNOFF3, 
     &            EDIR1,EC1,ET1,DRIP)
      END IF
      IF (MODEL_TYPE.EQ.1) THEN
      CALL SMFLX_SAC (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &        SH2O,SLOPE,KDT,FRZFACT,0.0,SMCREF,
     &        SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &        SHDFAC,CMCMAX,SFCTMP,ETP1,ETA1,
     &        RUNOFF1,RUNOFF2,RUNOFF3,EDIR1,EC1,ET1,DRIP, 
     &        FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP,wcrit)
      END IF
C ----------------------------------------------------------------------
C CONVERT MODELED EVAPOTRANSPIRATION FROM 'M S-1' TO 'KG M-2 S-1'.
C ----------------------------------------------------------------------
c        ETA = ETA1 * 1000.0

C ----------------------------------------------------------------------
c        EDIR = EDIR1 * 1000.0
c        EC = EC1 * 1000.0
c        ETT = ETT1 * 1000.0
c        ET(1) = ET1(1) * 1000.0
c        ET(2) = ET1(2) * 1000.0
c        ET(3) = ET1(3) * 1000.0
c        ET(4) = ET1(4) * 1000.0
C ----------------------------------------------------------------------

      ENDIF

C ----------------------------------------------------------------------
C       CONVERT MODELED EVAPOTRANSPIRATION FM  M S-1  TO  KG M-2 S-1
CBL Move this conversion to SFLX so that snow and snow-free cases handled
CBL equally and together to avoid confusion or mistakes, ET components 
CBL will all be in M/S throughout the code until the bottom of SFLX.
C ----------------------------------------------------------------------
CBL        ETA = ETA1 * 1000.0
CBL
C ----------------------------------------------------------------------
CBL Re-enumerate ett1
      ETT = 0.0
      EDIR = EDIR1
      EC = EC1
      DO K = 1,NSOIL
        ET(K) = ET1(K)
c        ET(1) = ET1(1).0
c        ET(2) = ET1(2).0
c        ET(3) = ET1(3).0
c        ET(4) = ET1(4).0
        ETT = ETT + ET1(K)
      ENDDO
C ----------------------------------------------------------------------
      if (prflag==1) then
      write(*,*)'nopac prcp edir ec ett eta'
      write(*,*)edir,prcp1,ec,ett,eta
      endif
C ----------------------------------------------------------------------
C BASED ON ETP AND E VALUES, DETERMINE BETA
C ----------------------------------------------------------------------
CBL use ETA0 insteat of ETA1 to keep units straight

      ETA0 = ETA1 * 1000
      IF ( ETP .LE. 0.0 ) THEN
        BETA = 0.0
        IF ( ETP .LT. 0.0 ) THEN
          BETA = 1.0
c          ETA = ETP
        ENDIF
      ELSE
CBL        BETA = ETA / ETP
CBL use ETA0 instead
         BETA = ETA0 / ETP
      ENDIF
      if (prflag==1) then
      write(*,*)'nopac beta eta etp',beta,eta,etp
      endif

C ----------------------------------------------------------------------
C GET SOIL THERMAL DIFFUXIVITY/CONDUCTIVITY FOR TOP SOIL LYR,
C CALC. ADJUSTED TOP LYR SOIL TEMP AND ADJUSTED SOIL FLUX, THEN
C CALL SHFLX TO COMPUTE/UPDATE SOIL HEAT FLUX AND SOIL TEMPS.
C ----------------------------------------------------------------------
      CALL TDFCND (DF1,SMC(1),QUARTZ,SMCMAX,SH2O(1))

C ----------------------------------------------------------------------
C VEGETATION GREENNESS FRACTION REDUCTION IN SUBSURFACE HEAT FLUX 
C VIA REDUCTION FACTOR, WHICH IS CONVENIENT TO APPLY HERE TO THERMAL 
C DIFFUSIVITY THAT IS LATER USED IN HRT TO COMPUTE SUB SFC HEAT FLUX
C (SEE ADDITIONAL COMMENTS ON VEG EFFECT SUB-SFC HEAT FLX IN 
C ROUTINE SFLX)
C ----------------------------------------------------------------------
      DF1 = DF1 * EXP(SBETA*SHDFAC)

C ----------------------------------------------------------------------
C COMPUTE INTERMEDIATE TERMS PASSED TO ROUTINE HRT (VIA ROUTINE 
C SHFLX BELOW) FOR USE IN COMPUTING SUBSURFACE HEAT FLUX IN HRT
C ----------------------------------------------------------------------
      YYNUM = FDOWN - EMISS*SIGMA * T24
      YY = SFCTMP + (YYNUM/RCH+TH2-SFCTMP-BETA*EPSCA) / RR
c      write(*,*)'nopac above shflx'
c      write(*,*)'yy sfctmp yynum rch th2 beta epsca rr'
c      write(*,*)yy,sfctmp,yynum,rch,th2,beta,epsca,rr

      ZZ1 = DF1 / ( -0.5 * ZSOIL(1) * RCH * RR ) + 1.0

      CALL SHFLX (SSOIL,STC,SMC,SMCMAX,NSOIL,T1,DT,YY,ZZ1,ZSOIL,
     &            TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &            QUARTZ,CSOIL)

C ----------------------------------------------------------------------
C SET FLX1 AND FLX3 (SNOPACK PHASE CHANGE HEAT FLUXES) TO ZERO SINCE
C THEY ARE NOT USED HERE IN SNOPAC.  FLX2 (FREEZING RAIN HEAT FLUX) WAS
C SIMILARLY INITIALIZED IN THE PENMAN ROUTINE.
C ----------------------------------------------------------------------
      FLX1 = 0.0
      FLX3 = 0.0

C ----------------------------------------------------------------------
C END SUBROUTINE NOPAC
C ----------------------------------------------------------------------
      RETURN
      END

C ------------ This is partial snow cover PENMAN treatment
C              made by Mike Ek
C  IF SNEQV.GT.0.0,AND THEN SNCOVR=1,THIS TREATMENT IS THE SAME 
C  AS NCAR TREATMENT 
C  YOLONG ----- Oct 2007
C --------------------------------------------------------------------
      SUBROUTINE PENMAN (SFCTMP,SFCPRS,CH,T2V,TH2,PRCP,FDOWN,T24,SSOIL,
     &                   Q2,Q2SAT,ETP,RCH,EPSCA,RR,SNOWNG,FRZGRA,
     &                   DQSDT2,FLX2,EMISS,EMISSNOW,SNCOVR)
      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE PENMAN
C ----------------------------------------------------------------------
C CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT.  VARIOUS
C PARTIAL SUMS/PRODUCTS ARE ALSO CALCULATED AND PASSED BACK TO THE
C CALLING ROUTINE FOR LATER USE.
C ----------------------------------------------------------------------
      LOGICAL SNOWNG
      LOGICAL FRZGRA

      REAL A
      REAL BETA
      REAL CH
      REAL CP
      REAL CPH2O
      REAL CPICE
      REAL DELTA
      REAL DQSDT2
      REAL ELCP
      REAL EPSCA
      REAL ETP
      REAL FDOWN
      REAL FLX2
      REAL FNET
      REAL LSUBC
      REAL LSUBF
      REAL LSUBS
      REAL PRCP
      REAL Q2
      REAL Q2SAT
      REAL R
      REAL RAD
      REAL RCH
      REAL RHO
      REAL RR
      REAL SSOIL
      REAL SFCPRS
      REAL SFCTMP
      REAL SIGMA
      REAL EMISS
      REAL EMISSNOW
      REAL SNCOVR
      REAL T24
      REAL T2V
      REAL TH2
      REAL SNOWC 

      REAL EMISSI
      REAL ELCP1
      REAL LVS

      PARAMETER(CP = 1004.6)
      PARAMETER(CPH2O = 4.218E+3)
      PARAMETER(CPICE = 2.106E+3)
      PARAMETER(R = 287.04)
      PARAMETER(ELCP = 2.4888E+3)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(LSUBC = 2.501000E+6)
      PARAMETER(LSUBS = 2.83E+6)
      PARAMETER(SIGMA = 5.67E-8)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C ----------------------------------------------------------------------
C account for partial snow cover fraction
C comment out for this sentence for EMC version
        SNOWC=SNCOVR
C        IF (SNOWC.GT.0.0) SNOWC=1.0       ! ONLY FOR NCAR VERSION
        EMISSI=(1.0-SNOWC)*EMISS + SNOWC*EMISSNOW
        ELCP1=(1.0-SNOWC)*ELCP + SNOWC*ELCP*LSUBS/LSUBC
        LVS=(1.0-SNOWC)*LSUBC + SNOWC*LSUBS

        FLX2 = 0.0

C ----------------------------------------------------------------------
C PREPARE PARTIAL QUANTITIES FOR PENMAN EQUATION.
C ----------------------------------------------------------------------
      DELTA = ELCP1 * DQSDT2
      T24 = SFCTMP * SFCTMP * SFCTMP * SFCTMP
      RR = EMISSI*T24 * 6.48E-8 / (SFCPRS * CH) + 1.0
      RHO = SFCPRS / (R * T2V)
      RCH = RHO * CP * CH

C ----------------------------------------------------------------------
C ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT
C EFFECTS CAUSED BY FALLING PRECIPITATION.
C ----------------------------------------------------------------------
      IF (.NOT. SNOWNG) THEN
        IF (PRCP .GT. 0.0) RR = RR + CPH2O*PRCP/RCH
          ELSE
         RR = RR + CPICE*PRCP/RCH
      ENDIF

      FNET = FDOWN -  EMISSI*SIGMA * T24- SSOIL
c      write(*,*)'penman fnet fdown emissi sigma t24 ssoil'
c      write(*,*)fnet,fdown,emissi,sigma,t24,ssoil
C ----------------------------------------------------------------------
C INCLUDE THE LATENT HEAT EFFECTS OF FRZNG RAIN CONVERTING TO ICE ON
C IMPACT IN THE CALCULATION OF FLX2 AND FNET.
C ----------------------------------------------------------------------
      IF (FRZGRA) THEN
        FLX2 = -LSUBF * PRCP
        FNET = FNET - FLX2
      ENDIF

C ----------------------------------------------------------------------
C FINISH PENMAN EQUATION CALCULATIONS.
C ----------------------------------------------------------------------
      RAD = FNET/RCH + TH2 - SFCTMP
c      write(*,*)'penman rad fnet rch th2 sfctmp'
c      write(*,*)rad,fnet,rch,th2,sfctmp
      A = ELCP1 * (Q2SAT - Q2)
      EPSCA = (A*RR + RAD*DELTA) / (DELTA + RR)
c      write(*,*)'penman epsca a rr rad delta'
c      write(*,*)epsca,a,rr,rad,delta
      ETP = EPSCA * RCH / LVS
c      write(*,*)'penman etp epsca rch lvs'
c      write(*,*)etp,epsca,rch,lvs

C ----------------------------------------------------------------------
C END SUBROUTINE PENMAN
C --------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE ROSR12 (P,A,B,C,D,DELTA,NSOIL)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE ROSR12
C ----------------------------------------------------------------------
C INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
C ###                                            ### ###  ###   ###  ###
C #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
C #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
C # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
C # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
C # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
C # .                                          .   # #  .   # = #   .  #
C # .                                          .   # #  .   #   #   .  #
C # .                                          .   # #  .   #   #   .  #
C # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
C # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
C # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
C ###                                            ### ###  ###   ###  ###
C ----------------------------------------------------------------------
      INTEGER K
      INTEGER KK
      INTEGER NSOIL
      
      REAL A(NSOIL)
      REAL B(NSOIL)
      REAL C(NSOIL)
      REAL D(NSOIL)
      REAL DELTA(NSOIL)
      REAL P(NSOIL)
      
C ----------------------------------------------------------------------
C INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
C ----------------------------------------------------------------------
      C(NSOIL) = 0.0

C ----------------------------------------------------------------------
C SOLVE THE COEFS FOR THE 1ST SOIL LAYER
C ----------------------------------------------------------------------
      P(1) = -C(1) / B(1)
      DELTA(1) = D(1) / B(1)

C ----------------------------------------------------------------------
C SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        P(K) = -C(K) * ( 1.0 / (B(K) + A (K) * P(K-1)) )
        DELTA(K) = (D(K)-A(K)*DELTA(K-1))*(1.0/(B(K)+A(K)*P(K-1)))
      END DO

C ----------------------------------------------------------------------
C SET P TO DELTA FOR LOWEST SOIL LAYER
C ----------------------------------------------------------------------
      P(NSOIL) = DELTA(NSOIL)

C ----------------------------------------------------------------------
C ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
         KK = NSOIL - K + 1
         P(KK) = P(KK) * P(KK+1) + DELTA(KK)
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE ROSR12
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SFCDIF_ALT (ZLM,Z0,THZ0,THLM,SFCSPD,CZIL,AKMS,AKHS)
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SFCDIF_ALT
C ----------------------------------------------------------------------
C CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
C SEE CHEN ET AL (1997, BLM)
C ----------------------------------------------------------------------
      
      REAL WWST, WWST2, G, VKRM, EXCM, BETA, BTG, ELFC, WOLD, WNEW
      REAL PIHF, EPSU2, EPSUST, EPSIT, EPSA, ZTMIN, ZTMAX, HPBL, SQVISC
      REAL RIC, RRIC, FHNEU, RFC, RFAC, ZZ, PSLMU, PSLMS, PSLHU, PSLHS
      REAL XX, PSPMU, YY, PSPMS, PSPHU, PSPHS, ZLM, Z0, THZ0, THLM
      REAL SFCSPD, CZIL, AKMS, AKHS, ZILFC, ZU, ZT, RDZ, CXCH
      REAL DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT
      REAL RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4
      REAL XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN, RLMA
CCC   ......REAL ZTFC
      
      INTEGER ITRMX, ILECH, ITR
      
      PARAMETER
     &     (WWST=1.2,WWST2=WWST*WWST,G=9.8,VKRM=0.40,EXCM=0.001
     &     ,BETA=1./270.,BTG=BETA*G,ELFC=VKRM*BTG
     &     ,WOLD=.15,WNEW=1.-WOLD,ITRMX=05,PIHF=3.14159265/2.)
C ----------------------------------------------------------------------
      PARAMETER
     &     (EPSU2=1.E-4,EPSUST=0.07,EPSIT=1.E-4,EPSA=1.E-8
     &     ,ZTMIN=-5.,ZTMAX=1.,HPBL=1000.0
     &     ,SQVISC=258.2)
C ----------------------------------------------------------------------
      PARAMETER
     &     (RIC=0.183,RRIC=1.0/RIC,FHNEU=0.8,RFC=0.191
     &     ,RFAC=RIC/(FHNEU*RFC*RFC))
      
C ----------------------------------------------------------------------
C NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
C ----------------------------------------------------------------------
C LECH'S SURFACE FUNCTIONS
C ----------------------------------------------------------------------
      PSLMU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLMS(ZZ)=ZZ*RRIC-2.076*(1.-1./(ZZ+1.))
      PSLHU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLHS(ZZ)=ZZ*RFAC-2.076*(1.-1./(ZZ+1.))
      
C ----------------------------------------------------------------------
C PAULSON'S SURFACE FUNCTIONS
C ----------------------------------------------------------------------
      PSPMU(XX)=-2.*log((XX+1.)*0.5)-log((XX*XX+1.)*0.5)+2.*ATAN(XX)
     &     -PIHF
      PSPMS(YY)=5.*YY
      PSPHU(XX)=-2.*log((XX*XX+1.)*0.5)
      PSPHS(YY)=5.*YY

C ----------------------------------------------------------------------
C THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
C OVER SOLID SURFACE (LAND, SEA-ICE).  
C ----------------------------------------------------------------------
      ILECH=0
      
C ----------------------------------------------------------------------
C     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
C     C......ZTFC=0.1
C     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
C ----------------------------------------------------------------------
      ZILFC=-CZIL*VKRM*SQVISC
      
C ----------------------------------------------------------------------
      ZU=Z0
C     C.......ZT=Z0*ZTFC
      RDZ=1./ZLM
      CXCH=EXCM*RDZ
      DTHV=THLM-THZ0     
      DU2=MAX(SFCSPD*SFCSPD,EPSU2)

C ----------------------------------------------------------------------
C BELJARS CORRECTION OF USTAR
C ----------------------------------------------------------------------
      BTGH=BTG*HPBL
ccc   If statements to avoid TANGENT LINEAR problems near zero
      IF (BTGH*AKHS*DTHV .NE. 0.0) THEN
         WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
      ELSE
         WSTAR2=0.0
      ENDIF
      USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

C ----------------------------------------------------------------------
C ZILITINKEVITCH APPROACH FOR ZT
C ----------------------------------------------------------------------
      ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0

C ----------------------------------------------------------------------
      ZSLU=ZLM+ZU
      ZSLT=ZLM+ZT
C     PRINT*,'ZSLT=',ZSLT
C     PRINT*,'ZLM=',ZLM
C     PRINT*,'ZT=',ZT
C     
      RLOGU=log(ZSLU/ZU)
      RLOGT=log(ZSLT/ZT)
C     
      RLMO=ELFC*AKHS*DTHV/USTAR**3
C     PRINT*,'RLMO=',RLMO
C     PRINT*,'ELFC=',ELFC
C     PRINT*,'AKHS=',AKHS
C     PRINT*,'DTHV=',DTHV
C     PRINT*,'USTAR=',USTAR
      
      DO ITR=1,ITRMX
C ----------------------------------------------------------------------
C 1./MONIN-OBUKKHOV LENGTH-SCALE
C ----------------------------------------------------------------------

         ZETALT=MAX(ZSLT*RLMO,ZTMIN)
         RLMO=ZETALT/ZSLT
         ZETALU=ZSLU*RLMO
         ZETAU=ZU*RLMO
         ZETAT=ZT*RLMO
 
CC Bound zeta (z/L) at 1.0 so as to reduce estimate errors via the findings
CC of Yague 2006, which show ~ constant stability function when z/L > 1
C         ZETALT=MAX(ZSLT*RLMO,ZTMIN)
C         RLMO=ZETALT/ZSLT
C         ZETALT = MIN(ZETALT,1.0)
C         ZETALU= MIN(ZSLU*RLMO,1.0)
C         ZETAU= MIN(ZU*RLMO,1.0)
C         ZETAT= MIN(ZT*RLMO,1.0)
        
         IF(ILECH.EQ.0) THEN
            IF(RLMO.LT.0.)THEN
               XLU4=1.-16.*ZETALU
               XLT4=1.-16.*ZETALT
               XU4 =1.-16.*ZETAU
               XT4 =1.-16.*ZETAT
               
               XLU=SQRT(SQRT(XLU4))
               XLT=SQRT(SQRT(XLT4))
               XU =SQRT(SQRT(XU4))
               XT =SQRT(SQRT(XT4))
               
               PSMZ=PSPMU(XU)
C     PRINT*,'-----------1------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSPMU(ZETAU)=',PSPMU(ZETAU)
C     PRINT*,'XU=',XU
C     PRINT*,'------------------------'
               SIMM=PSPMU(XLU)-PSMZ+RLOGU
               PSHZ=PSPHU(XT)
               SIMH=PSPHU(XLT)-PSHZ+RLOGT
            ELSE
               ZETALU=MIN(ZETALU,ZTMAX)
               ZETALT=MIN(ZETALT,ZTMAX)
               PSMZ=PSPMS(ZETAU)
C     PRINT*,'-----------2------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSPMS(ZETAU)=',PSPMS(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSPMS(ZETALU)-PSMZ+RLOGU
               PSHZ=PSPHS(ZETAT)
               SIMH=PSPHS(ZETALT)-PSHZ+RLOGT
C ------- avoid numerical problem ------------------------------------       
               IF(SIMM.EQ.0.0) SIMM=-1.0E-6
               IF(SIMH.EQ.0.0) SIMH=1.0E-6
            ENDIF
         ELSE
C ----------------------------------------------------------------------
C LECH'S FUNCTIONS
C ----------------------------------------------------------------------
            IF(RLMO.LT.0.)THEN
               PSMZ=PSLMU(ZETAU)
C     PRINT*,'-----------3------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSLMU(ZETAU)=',PSLMU(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSLMU(ZETALU)-PSMZ+RLOGU
               PSHZ=PSLHU(ZETAT)
               SIMH=PSLHU(ZETALT)-PSHZ+RLOGT
            ELSE
               ZETALU=MIN(ZETALU,ZTMAX)
               ZETALT=MIN(ZETALT,ZTMAX)
C     
               PSMZ=PSLMS(ZETAU)
C     PRINT*,'-----------4------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSLMS(ZETAU)=',PSLMS(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSLMS(ZETALU)-PSMZ+RLOGU
               PSHZ=PSLHS(ZETAT)
               SIMH=PSLHS(ZETALT)-PSHZ+RLOGT
            ENDIF
         ENDIF
C ----------------------------------------------------------------------
C BELJAARS CORRECTION FOR USTAR
C ----------------------------------------------------------------------
         USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

C ----------------------------------------------------------------------
C ZILITINKEVITCH FIX FOR ZT
C ----------------------------------------------------------------------
         ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0
         
         ZSLT=ZLM+ZT
         RLOGT=log(ZSLT/ZT)
C-----------------------------------------------------------------------
         USTARK=USTAR*VKRM
         AKMS=MAX(USTARK/SIMM,CXCH)
         AKHS=MAX(USTARK/SIMH,CXCH)
C       PRINT*,'AKMS=',AKMS
C       PRINT*,'SIMM=',SIMM
C       PRINT*,'AKHS=',AKHS
C       PRINT*,'SIMH=',SIMH
C-----------------------------------------------------------------------
C IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
C-----------------------------------------------------------------------
         IF (BTGH*AKHS*DTHV .NE. 0.0) THEN
            WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
         ELSE
            WSTAR2=0.0
         ENDIF
         RLMN=ELFC*AKHS*DTHV/USTAR**3
C-----------------------------------------------------------------------
         RLMA=RLMO*WOLD+RLMN*WNEW
C-----------------------------------------------------------------------
C     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
C-----------------------------------------------------------------------
         RLMO=RLMA
C-----------------------------------------------------------------------
      END DO

C     PRINT*,'----------------------------'
C     PRINT*,'SFCDIF OUTPUT !  ! ! ! ! ! ! ! !  !   !    !'
C     
C     PRINT*,'ZLM=',ZLM
C     PRINT*,'Z0=',Z0
C     PRINT*,'THZ0=',THZ0
C     PRINT*,'THLM=',THLM
C     PRINT*,'AKMS=',AKMS
C     PRINT*,'AKHS=',AKHS
C     PRINT*,'----------------------------'
C     
C ----------------------------------------------------------------------
C END SUBROUTINE SFCDIF_ALT
C ----------------------------------------------------------------------
      RETURN
      END


      SUBROUTINE SFCDIF (ZLM,Z0,THZ0,THLM,SFCSPD,CZIL,AKMS,AKHS)
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SFCDIF
C ----------------------------------------------------------------------
C CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
C SEE CHEN ET AL (1997, BLM)
C ----------------------------------------------------------------------
      
      REAL WWST, WWST2, G, VKRM, EXCM, BETA, BTG, ELFC, WOLD, WNEW
      REAL PIHF, EPSU2, EPSUST, EPSIT, EPSA, ZTMIN, ZTMAX, HPBL, SQVISC
      REAL RIC, RRIC, FHNEU, RFC, RFAC, ZZ, PSLMU, PSLMS, PSLHU, PSLHS
      REAL XX, PSPMU, YY, PSPMS, PSPHU, PSPHS, ZLM, Z0, THZ0, THLM
      REAL SFCSPD, CZIL, AKMS, AKHS, ZILFC, ZU, ZT, RDZ, CXCH
      REAL DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT
      REAL RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4
      REAL XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN, RLMA
CCC   ......REAL ZTFC
      
      INTEGER ITRMX, ILECH, ITR
      
      PARAMETER
     &     (WWST=1.2,WWST2=WWST*WWST,G=9.8,VKRM=0.40,EXCM=0.001
     &     ,BETA=1./270.,BTG=BETA*G,ELFC=VKRM*BTG
     &     ,WOLD=.15,WNEW=1.-WOLD,ITRMX=05,PIHF=3.14159265/2.)
C ----------------------------------------------------------------------
      PARAMETER
     &     (EPSU2=1.E-4,EPSUST=0.07,EPSIT=1.E-4,EPSA=1.E-8
     &     ,ZTMIN=-5.,ZTMAX=1.,HPBL=1000.0
     &     ,SQVISC=258.2)
C ----------------------------------------------------------------------
      PARAMETER
     &     (RIC=0.183,RRIC=1.0/RIC,FHNEU=0.8,RFC=0.191
     &     ,RFAC=RIC/(FHNEU*RFC*RFC))
      
C ----------------------------------------------------------------------
C NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
C ----------------------------------------------------------------------
C LECH'S SURFACE FUNCTIONS
C ----------------------------------------------------------------------
      PSLMU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLMS(ZZ)=ZZ*RRIC-2.076*(1.-1./(ZZ+1.))
      PSLHU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLHS(ZZ)=ZZ*RFAC-2.076*(1.-1./(ZZ+1.))
      
C ----------------------------------------------------------------------
C PAULSON'S SURFACE FUNCTIONS
C ----------------------------------------------------------------------
      PSPMU(XX)=-2.*log((XX+1.)*0.5)-log((XX*XX+1.)*0.5)+2.*ATAN(XX)
     &     -PIHF
      PSPMS(YY)=5.*YY
      PSPHU(XX)=-2.*log((XX*XX+1.)*0.5)
      PSPHS(YY)=5.*YY

C ----------------------------------------------------------------------
C THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
C OVER SOLID SURFACE (LAND, SEA-ICE).  
C ----------------------------------------------------------------------
      ILECH=0
      
C ----------------------------------------------------------------------
C     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
C     C......ZTFC=0.1
C     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
C ----------------------------------------------------------------------
      ZILFC=-CZIL*VKRM*SQVISC
      
C ----------------------------------------------------------------------
      ZU=Z0
C     C.......ZT=Z0*ZTFC
      RDZ=1./ZLM
      CXCH=EXCM*RDZ
      DTHV=THLM-THZ0     
      DU2=MAX(SFCSPD*SFCSPD,EPSU2)

C ----------------------------------------------------------------------
C BELJARS CORRECTION OF USTAR
C ----------------------------------------------------------------------
      BTGH=BTG*HPBL
ccc   If statements to avoid TANGENT LINEAR problems near zero
      IF (BTGH*AKHS*DTHV .NE. 0.0) THEN
         WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
      ELSE
         WSTAR2=0.0
      ENDIF
      USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

C ----------------------------------------------------------------------
C ZILITINKEVITCH APPROACH FOR ZT
C ----------------------------------------------------------------------
      ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0

C ----------------------------------------------------------------------
      ZSLU=ZLM+ZU
      ZSLT=ZLM+ZT
C     PRINT*,'ZSLT=',ZSLT
C     PRINT*,'ZLM=',ZLM
C     PRINT*,'ZT=',ZT
C     
      RLOGU=log(ZSLU/ZU)
      RLOGT=log(ZSLT/ZT)
C     
      RLMO=ELFC*AKHS*DTHV/USTAR**3
C     PRINT*,'RLMO=',RLMO
C     PRINT*,'ELFC=',ELFC
C     PRINT*,'AKHS=',AKHS
C     PRINT*,'DTHV=',DTHV
C     PRINT*,'USTAR=',USTAR
      
      DO ITR=1,ITRMX
C ----------------------------------------------------------------------
C 1./MONIN-OBUKKHOV LENGTH-SCALE
C ----------------------------------------------------------------------
         ZETALT=MAX(ZSLT*RLMO,ZTMIN)
         RLMO=ZETALT/ZSLT
         ZETALU=ZSLU*RLMO
         ZETAU=ZU*RLMO
         ZETAT=ZT*RLMO
         
         IF(ILECH.EQ.0) THEN
            IF(RLMO.LT.0.)THEN
               XLU4=1.-16.*ZETALU
               XLT4=1.-16.*ZETALT
               XU4 =1.-16.*ZETAU
               XT4 =1.-16.*ZETAT
               
               XLU=SQRT(SQRT(XLU4))
               XLT=SQRT(SQRT(XLT4))
               XU =SQRT(SQRT(XU4))
               XT =SQRT(SQRT(XT4))
               
               PSMZ=PSPMU(XU)
C     PRINT*,'-----------1------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSPMU(ZETAU)=',PSPMU(ZETAU)
C     PRINT*,'XU=',XU
C     PRINT*,'------------------------'
               SIMM=PSPMU(XLU)-PSMZ+RLOGU
               PSHZ=PSPHU(XT)
               SIMH=PSPHU(XLT)-PSHZ+RLOGT
            ELSE
               ZETALU=MIN(ZETALU,ZTMAX)
               ZETALT=MIN(ZETALT,ZTMAX)
               PSMZ=PSPMS(ZETAU)
C     PRINT*,'-----------2------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSPMS(ZETAU)=',PSPMS(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
	       IF (ZETALU .EQ. ZTMAX) PSMZ=PSPMS(ZTMAX*ZU/ZSLU)
               SIMM=PSPMS(ZETALU)-PSMZ+RLOGU
               PSHZ=PSPHS(ZETAT)
	       IF (ZETALT .EQ. ZTMAX) PSHZ=PSPHS(ZTMAX*ZT/(ZLM+ZT))
               SIMH=PSPHS(ZETALT)-PSHZ+RLOGT
C ------- avoid numerical problem ------------------------------------
                  IF(SIMM.LE.0.0.OR.SIMH.LE.0.0) THEN
                  WRITE(*,*) 'SIMM/SIMH<=0.0 and STOP to RUN Noah '
                  STOP
                  ENDIF     
CY               IF(SIMM.EQ.0.0) SIMM=-1.0E-6
CY               IF(SIMH.EQ.0.0) SIMH=1.0E-6
            ENDIF
         ELSE
C ----------------------------------------------------------------------
C LECH'S FUNCTIONS
C ----------------------------------------------------------------------
            IF(RLMO.LT.0.)THEN
               PSMZ=PSLMU(ZETAU)
C     PRINT*,'-----------3------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSLMU(ZETAU)=',PSLMU(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSLMU(ZETALU)-PSMZ+RLOGU
               PSHZ=PSLHU(ZETAT)
               SIMH=PSLHU(ZETALT)-PSHZ+RLOGT
            ELSE
               ZETALU=MIN(ZETALU,ZTMAX)
               ZETALT=MIN(ZETALT,ZTMAX)
C     
               PSMZ=PSLMS(ZETAU)
C     PRINT*,'-----------4------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSLMS(ZETAU)=',PSLMS(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSLMS(ZETALU)-PSMZ+RLOGU
               PSHZ=PSLHS(ZETAT)
               SIMH=PSLHS(ZETALT)-PSHZ+RLOGT
            ENDIF
         ENDIF
C ----------------------------------------------------------------------
C BELJAARS CORRECTION FOR USTAR
C ----------------------------------------------------------------------
         USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

C ----------------------------------------------------------------------
C ZILITINKEVITCH FIX FOR ZT
C ----------------------------------------------------------------------
         ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0
         
         ZSLT=ZLM+ZT
         RLOGT=log(ZSLT/ZT)
C-----------------------------------------------------------------------
         USTARK=USTAR*VKRM
         AKMS=MAX(USTARK/SIMM,CXCH)
         AKHS=MAX(USTARK/SIMH,CXCH)
C       PRINT*,'AKMS=',AKMS
C       PRINT*,'SIMM=',SIMM
C       PRINT*,'AKHS=',AKHS
C       PRINT*,'SIMH=',SIMH
C-----------------------------------------------------------------------
C IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
C-----------------------------------------------------------------------
         IF (BTGH*AKHS*DTHV .NE. 0.0) THEN
            WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
         ELSE
            WSTAR2=0.0
         ENDIF
         RLMN=ELFC*AKHS*DTHV/USTAR**3
C-----------------------------------------------------------------------
         RLMA=RLMO*WOLD+RLMN*WNEW
C-----------------------------------------------------------------------
C     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
C-----------------------------------------------------------------------
         RLMO=RLMA
C-----------------------------------------------------------------------
      END DO

C     PRINT*,'----------------------------'
C     PRINT*,'SFCDIF OUTPUT !  ! ! ! ! ! ! ! !  !   !    !'
C     
C     PRINT*,'ZLM=',ZLM
C     PRINT*,'Z0=',Z0
C     PRINT*,'THZ0=',THZ0
C     PRINT*,'THLM=',THLM
C     PRINT*,'AKMS=',AKMS
C     PRINT*,'AKHS=',AKHS
C     PRINT*,'----------------------------'
C     
C ----------------------------------------------------------------------
C END SUBROUTINE SFCDIF
C ----------------------------------------------------------------------
      RETURN
      END


      SUBROUTINE SHFLX (SSOIL,STC,SMC,SMCMAX,NSOIL,T1,DT,YY,ZZ1,ZSOIL,
     &                  TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &                  QUARTZ,CSOIL)
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SHFLX
C ----------------------------------------------------------------------
C UPDATE THE TEMPERATURE STATE OF THE SOIL COLUMN BASED ON THE THERMAL
C DIFFUSION EQUATION AND UPDATE THE FROZEN SOIL MOISTURE CONTENT BASED
C ON THE TEMPERATURE.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER ICE
      INTEGER IFRZ
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL BEXP
      REAL CSOIL
      REAL DF1
      REAL DT
      REAL F1
      REAL PSISAT
      REAL QUARTZ
      REAL RHSTS(NSOLD)
      REAL SSOIL
      REAL SH2O(NSOIL)
      REAL SMC(NSOIL)
      REAL SMCMAX
      REAL SMCWLT
      REAL STC(NSOIL)
      REAL STCF(NSOLD)
      REAL T0
      REAL T1
      REAL TBOT
      REAL YY
      REAL ZBOT
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      PARAMETER(T0 = 273.15)



C ----------------------------------------------------------------------
C HRT ROUTINE CALCS THE RIGHT HAND SIDE OF THE SOIL TEMP DIF EQN
C ----------------------------------------------------------------------
      IF (ICE.EQ.1) THEN

C ----------------------------------------------------------------------
C SEA-ICE CASE
C ----------------------------------------------------------------------
         CALL HRTICE (RHSTS,STC,NSOIL,ZSOIL,YY,ZZ1,DF1,AI,BI,CI)

         CALL HSTEP (STCF,STC,RHSTS,DT,NSOIL,AI,BI,CI)
         
      ELSE

C ----------------------------------------------------------------------
C LAND-MASS CASE
C ----------------------------------------------------------------------
c         write(*,*)'shflx0 rhsts stc1 smc1 yy'
c         write(*,*)rhsts(1),stc(1),smc(1),yy
c         write(*,*)'zz1 sh2o1 f1 df1 csoil'
c         write(*,*)zz1,sh2o(1),f1,df1,csoil
c         write(*,*)'tbot ai bi ci'
c         write(*,*)tbot,ai(1),bi(1),ci(1)

         CALL HRT (RHSTS,STC,SMC,SMCMAX,NSOIL,ZSOIL,YY,ZZ1,TBOT,
     &             ZBOT,PSISAT,SH2O,DT,
     &             BEXP,F1,DF1,QUARTZ,CSOIL,AI,BI,CI)

c         write(*,*)'shflx1 rhsts stc1 smc1 yy'
c         write(*,*)rhsts(1),stc(1),smc(1),yy
c         write(*,*)'zz1 sh2o1 f1 df1 csoil'
c         write(*,*)zz1,sh2o(1),f1,df1,csoil
c         write(*,*)'tbot ai bi ci'
c         write(*,*)tbot,ai(1),bi(1),ci(1)

         CALL HSTEP (STCF,STC,RHSTS,DT,NSOIL,AI,BI,CI)

c         write(*,*)'shflx2 rhsts stc1 smc1 yy'
c         write(*,*)rhsts(1),stc(1),smc(1),yy
c         write(*,*)'zz1 sh2o1 f1 df1 csoil'
c         write(*,*)zz1,sh2o(1),f1,df1,csoil
c         write(*,*)'tbot ai bi ci'
c         write(*,*)tbot,ai(1),bi(1),ci(1)

      ENDIF

      DO I = 1,NSOIL
         STC(I) = STCF(I)
      END DO
      
C ----------------------------------------------------------------------
C IN THE NO SNOWPACK CASE (VIA ROUTINE NOPAC BRANCH,) UPDATE THE GRND
C (SKIN) TEMPERATURE HERE IN RESPONSE TO THE UPDATED SOIL TEMPERATURE 
C PROFILE ABOVE.  (NOTE: INSPECTION OF ROUTINE SNOPAC SHOWS THAT T1
C BELOW IS A DUMMY VARIABLE ONLY, AS SKIN TEMPERATURE IS UPDATED
C DIFFERENTLY IN ROUTINE SNOPAC) 
C ----------------------------------------------------------------------
      T1 = (YY + (ZZ1 - 1.0) * STC(1)) / ZZ1

c      write(*,*)'shflx3 t1 yy zz1 stc1'
c      write(*,*)t1,yy,zz1,stc(1)
            
C ----------------------------------------------------------------------
C CALCULATE SURFACE SOIL HEAT FLUX
C ----------------------------------------------------------------------
      SSOIL = DF1 * (STC(1) - T1) / (0.5 * ZSOIL(1))

C ----------------------------------------------------------------------
C END SUBROUTINE SHFLX
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SMFLX_SAC (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &     SH2O,SLOPE,KDT,FRZFACT,SNCOVR,SMCREF,
     &     SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &     SHDFAC,CMCMAX,SFCTMP,ETP1,ETA1,
     &     RUNOFF1,RUNOFF2,RUNOFF3,EDIR1,EC1,ET1,DRIP,
     &     FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP,wcrit)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SMFLX_SAC
C ----------------------------------------------------------------------
C Solves the soil water balance using the Sacramento SMA (Burnash, 1973)
C The SAC-HT (heat transfer) code here was modfied from the version 
C received from Victor Koren, such that the frozen soil physics are done
C within the Noah model (SHFLX), but effects of frozen soil such as rate
C of infiltration and differences in liquid and total water are handled
C correctly within SAC-HT via the SMC,SH2O physical layers (volume frac)
C and xZxxC,xZxxH SAC storages (mm). Conversion between these quantities
C is done before and after the call to SAC.
C
C During snow-free conditions, PET and precipitation are scaled by canopy
C coverage (SHDFAC) so as not to double-count either quantity here.
C Effective precipitation is a combination of precipitation, drip water 
C from canopy, and dewfall.
C
C During snow cover, PET is scaled based on snow cover, while other ET 
C components are passed to this subroutine already rescaled by SNCOVR in
C the SNXPAC subroutine(s).  In that case, effective precip combines 
C melt water + (dewfall + precip) over (1-SNCOVR) area.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL BEXP
      REAL CMC
      REAL CMCMAX
      REAL DKSAT
      REAL DRIP
      REAL DT
      REAL DUMMY
      REAL DWSAT
      REAL EC1
      REAL EDIR1
      REAL ET1(NSOIL)
      REAL EXCESS
      REAL FRZFACT
      REAL KDT
      REAL PCPDRP
      REAL PRCP1
      REAL RHSCT
      REAL RHSTT(NSOLD)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL SHDFAC
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SICE(NSOLD)
      REAL SH2OA(NSOLD)
      REAL SH2OFG(NSOLD)
      REAL SLOPE
      REAL SMCMAX
      REAL SMCWLT
      REAL TRHSCT
      REAL ZSOIL(NSOIL)
      REAL EDMND
      REAL RAIM
      REAL TMP
	
      REAL FAC2	
      REAL FLIMIT

      REAL SFCTMP,ETA1,ETP1,TCI
      REAL FRZST(10)
      REAL FRZPAR(13)
      REAL SACST(6)
      REAL SACST_PRV(6)
      REAL SACPAR(16)
      real sacstdif(6)
      REAL FROST,PXV,TA,DWT,DWF,SURF,GRND,TET,DTDAY,SNCOVR,temp1,EDMND0
      INTEGER IVERS,NUPL,NSAC,NUP,PRFLAG,CELLID,MSTEP
      real sactot,sactot0,sacdif,frztot,frztot0,frzdif,runtot,runtot0
      real smctot,smctot0,smcdif,sh2otot,sh2otot0,sh2odif
      real tempet1,difsum,diftot,difliq,stot0,sliq0,stot,sliq
      REAL ROIMP,SDRO,SSUR,SIF,BFS,BFP,WE,QS,QG,Q,DS,BAL,DSH
      REAL critsmc,wcrit,scale,SMCREF,fcrit
      INTEGER IFRZE,LWE,ISC

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE.
C ----------------------------------------------------------------------
      if (prflag==1) then
      write(*,*)'smflx_sac prcp1 sncovr',prcp1,sncovr
      endif
C ----------------------------------------------------------------------
C ADD EXCESS MOISTURE FROM CANOPY TO THE SOIL, VIA THE 
C COMPUTATION OF THE RIGHT HAND SIDE OF THE CANOPY EQN TERM ( RHSCT )
C ----------------------------------------------------------------------
      RHSCT = SHDFAC * PRCP1 - EC1

C ----------------------------------------------------------------------
C CONVERT RHSCT (A RATE) TO TRHSCT (AN AMOUNT) AND ADD IT TO EXISTING
C CMC.  IF RESULTING AMT EXCEEDS MAX CAPACITY, IT BECOMES DRIP AND WILL
C FALL TO THE GRND.
C ----------------------------------------------------------------------
      DRIP = 0.
      TRHSCT = DT * RHSCT
      CMC = CMC + TRHSCT
      IF (CMC .GT. CMCMAX) DRIP = CMC - CMCMAX
      CMC = MIN(CMC,CMCMAX)      
      IF (CMC .LT. 1.E-20) CMC=0.0

C ----------------------------------------------------------------------
C PCPDRP IS THE COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE
C SOIL
C ----------------------------------------------------------------------

      PCPDRP = (1. - SHDFAC) * PRCP1 + DRIP / DT
      if (prflag==1) then
      write(*,*)'smflx shdfac pcpdrp drip',shdfac,pcpdrp,drip
      endif

C ----------------------------------------------------------------------
C CONVERT MOISTURE FLUXES FROM M/S TO MM OVER TIMESTEP FOR SAC,
C REDUCE PET TO ACCOUNT FOR AVAILABLE SURFACE AREA, SO AS NOT 
C TO DOUBLE-COUNT UNDER CANOPY OR SNOW, CONSIDER FROZEN SOIL (IVERS=1)
C REDUCE ETP BY A FACTOR TO ET MAPS
C ----------------------------------------------------------------------

      EDMND0 = ETP1 * DT * (1 - SNCOVR) * (1 - SHDFAC) * 1000
cbl      EDMND = MAX(EDMND0,0.0) * 0.5
      EDMND = MAX(EDMND0,0.0)
cbl Employ the critical moisture content concept, such that soil
cbl evaporation and transpiration rates decrease linearly from
cbl unrestricted, at the critical point to 0 the wilting point.
cbl SMCREF (field capacity) computed in subroutine REDPRM with 
cbl other SMC constants to ensure consistency.
cbl      wcrit = 0.333 ---read from input
      critsmc = smcwlt + (smcref-smcwlt) * wcrit
c      wcrit = (smc(1) - smcwlt) / (smcref - smcwlt)
      fcrit = (smc(1) - smcwlt) / (critsmc - smcwlt)
      if(fcrit.lt.1.0.and.prflag==1)then
         write(*,*)'WCRIT',wcrit
         write(*,*)'smcwlt critsmc smcref'
         write(*,*)smcwlt,critsmc,smcref
         write(*,*)'smc1 fcrit',smc(1),fcrit
      endif

      EDMND=MIN(EDMND,(EDMND*fcrit))
      PXV = PCPDRP * DT * 1000
      TA = SFCTMP
      DTDAY = DT / 86400
      IVERS = 1

      NUPL = 2
      NSAC = NSOIL
c      write(*,*)'cellid',CELLID
      
C      IF (PRFLAG == 1) THEN
C         write(*,*) '1.SMC:',(SMC(I),I=1,4)
C         write(*,*) '1.SH2O:',(SH2O(I),I=1,4)
C         write(*,*) '1.SACST:',(SACST(I),I=1,5)
C         write(*,*) '1.FRZST:',(FRZST(I+5),I=1,5)
C      ENDIF

      smctot0=0.
      sh2otot0=0.
      runtot0=0.
      runtot=(runoff1+runoff2)
      do i = 1,nsoil
         smctot0 = smctot0 + smc(i)
         sh2otot0 = sh2otot0 + sh2o(i)
      enddo  
c      sacdif = sactot - sactot0
c      frzdif = frztot - frztot0
      if (prflag==1) then
      write(*,*)'-----------smflx TOP SMFLX------' 
      write(*,*)'prcp           edir            ec'
      write(*,*)pcpdrp*dt,edir1,ec1
      write(*,*)'ett            eta             etp'
      write(*,*)tempet1,eta1,edmnd
      write(*,*)'smctot0        sh2otot0'
      write(*,*)smctot0,sh2otot0
c      write(*,*)'sacdif        frzdif'
c      write(*,*)sacdif,frzdif
      write(*,*)'---------------------------------' 
      endif
      sactot0=0.
      frztot0=0.
      runtot0=0.
      stot0 = 0.
      sliq0 = 0.
      do i = 1,5
         sactot0 = sactot0 + sacst(i)
         frztot0 = frztot0 + frzst(i+5)
         if(prflag==1) write(*,*)'sac frz',sacst(i),frzst(i+5),sactot0
      enddo
      stot0 = stot0 + 1000 * (smc(1) - smcwlt) * (0.0 - zsoil(1))
      sliq0 = sliq0 + 1000 * (sh2o(1) - smcwlt) * (0.0 - zsoil(1))
      if(prflag==1) write(*,*)'smctot 1',stot0,smc(1),smcwlt,zsoil(1)
      do i = 2,nsoil
         stot0 = stot0+1000*(smc(i) - smcwlt) * (zsoil(i-1) - zsoil(i))
         sliq0 = sliq0+1000*(sh2o(i) - smcwlt) * (zsoil(i-1) - zsoil(i))
         if (prflag==1) write(*,*)'smctot',i,stot0,smc(i)         
      enddo
      diftot = sactot0 - stot0
      difliq = frztot0 - sliq0
      if (prflag==1) then
      write(*,*)'stot0 sliq0',stot0,sliq0
      write(*,*)'sactot0 frztot0',sactot0,frztot0
      write(*,*)'diftot difliq',diftot,difliq
      endif

C     Convert physical layers to SAC storages
      CALL FRZ2SAC1(PXV,ZSOIL,SMC,SH2O,NUPL,NSAC,SMCWLT,1.0,
     +     1.0,FRZST,SACST,FROST,PRFLAG,CELLID)

      DO I = 1,6
         SACST_PRV(I) = SACST(I)
      ENDDO

      if(prflag==1)write(*,*)'below frz2sac'
      sactot0=0.
      frztot0=0.
      runtot0=0.
      stot0 = 0.
      sliq0 = 0.
      do i = 1,5
         sactot0 = sactot0 + sacst(i)
         frztot0 = frztot0 + frzst(i+5)
         if(prflag==1) write(*,*)'sac frz',sacst(i),frzst(i+5),sactot0
      enddo
      stot0 = stot0 + 1000 * (smc(1) - smcwlt) * (0.0 - zsoil(1))
      sliq0 = sliq0 + 1000 * (sh2o(1) - smcwlt) * (0.0 - zsoil(1))
      if(prflag==1) write(*,*)'smctot 1',stot0,smc(1),smcwlt,zsoil(1)
      do i = 2,nsoil
         stot0 = stot0+1000*(smc(i) - smcwlt) * (zsoil(i-1) - zsoil(i))
         sliq0 = sliq0+1000*(sh2o(i) - smcwlt) * (zsoil(i-1) - zsoil(i))
         if (prflag==1) write(*,*)'smctot',i,stot0,smc(i)         
      enddo
      diftot = sactot0 - stot0
      difliq = frztot0 - sliq0
      if (prflag==1) then
      write(*,*)'stot0 sliq0',stot0,sliq0
      write(*,*)'sactot0 frztot0',sactot0,frztot0
      write(*,*)'diftot difliq',diftot,difliq
      endif

C Remove root zone uptake (transpiration) and update SAC storages
C Make sure water is available.  Assume this moiture comes from tension
C water, as this is where other evap is extracted from SAC storages

      CALL EXTRANSP (ET1,NSOIL,SMC,SH2O,SMCWLT,SACST,FRZST,DT,
     &     NUPL,NSAC,PRFLAG,CELLID,MSTEP,smcref,critsmc)

c      DO I = 1,6
c         SACST_PRV(I) = SACST(I)
c      ENDDO

      tempet1=0
      do i = 1,nsoil
         tempet1 = tempet1 + et1(i)
      enddo

c      IF (PRFLAG==1) THEN
c         write(*,*)'pxv,edmnd,ta,dt',pxv,edmnd,ta,dt
c         write(*,*)'sacst',(sacst(i),i=1,6)
c         write(*,*)'frzst',(frzst(i),i=6,10)
c         write(*,*)'sacpar',(sacpar(i),i=1,16)
c         write(*,*)'frzpar',(frzpar(i),i=1,13)
c         write(*,*)'nsoil nupl nsac ivers',nsoil,nupl,nsac,ivers
c         write(*,*)'surf grnd tet',surf,grnd,tet
c         write(*,*) 'SMC:',(SMC(I),I=1,4)
c         write(*,*) 'SH2O:',(SH2O(I),I=1,4)
c         write(*,*) 'SACST_PRV:',(SACST_PRV(I),I=1,5)
c      ENDIF

      sactot0=0.
      frztot0=0.
      runtot0=0.
      runtot0=(surf+grnd)
      do i = 1,5
         sactot0 = sactot0 + sacst(i)
      enddo  
      do i = 1,5
         frztot0 = frztot0 + frzst(i+5)
      enddo  
      if (prflag==1) then
      write(*,*)'-----------smflx ABOVE FLAND2------' 
      write(*,*)'prcp           edir            ec'
      write(*,*)pxv,edir1,ec1
      write(*,*)'ett            eta             etp'
      write(*,*)tempet1,tet,edmnd
      write(*,*)'surf           grnd            runoff'
      write(*,*)surf,grnd,runtot0
      write(*,*)'sactot        frztot          tci'
      write(*,*)sactot0,frztot0,tci
      write(*,*)'---------------------------------'
      endif
      
C      write(*,*)'cell id',cellid,'model step',mstep

      CALL FLAND2(PXV,EDMND,TA,DTDAY,SACST,FRZST,SACPAR,
     +     FRZPAR,NSOIL,NUPL,NSAC,IVERS,SURF,GRND,TCI,TET,
     +     SMC,SH2O,SACST_PRV,cellid,prflag)

C      CALL SAC1(DTDAY,PXV,EDMND,TCI,ROIMP,SDRO,SSUR,SIF,BFS,BFP,TET,
C     &     IFRZE,TA,LWE,WE,ISC,SNCOVR,SACPAR,SACST)


      difsum = 0
      do i = 1,6
         sacstdif(i) = sacst(i) - sacst_prv(i)
         if (prflag==1) write(*,*)'sacstdif',i,sacstdif(i)
         difsum = difsum + sacstdif(i)
      enddo
      if (prflag==1) write(*,*)'difsum',difsum

      sactot=0.
      frztot=0.
      runtot=0.
      runtot=(surf+grnd)
      do i = 1,5
         sactot = sactot + sacst(i)
      enddo
      do i = 1,5
         frztot = frztot + frzst(i+5)
      enddo  
      sacdif = sactot - sactot0
      frzdif = frztot - frztot0
      if (prflag==1) then
      write(*,*)'-----------smflx BELOW FLAND2------' 
      write(*,*)'prcp           edir            ec'
      write(*,*)pxv,edir1,ec1
      write(*,*)'ett            eta             etp'
      write(*,*)tempet1,tet,edmnd
      write(*,*)'surf           grnd            runoff'
      write(*,*)surf,grnd,runtot
      write(*,*)'sactot        frztot          tci'
      write(*,*)sactot,frztot,tci
      write(*,*)'sacdif        frzdif'
      write(*,*)sacdif,frzdif
      write(*,*)'---------------------------------' 
      endif
      DS = sactot - sactot0
      BAL = PXV-TET-SURF-GRND-DS
      if(prflag==1)then
      write(*,"(A,9(xf14.6))") 'bal sac-ht ',BAL,PXV,TET,RUNTOT,SURF,
     &     GRND,DS,sactot0,sactot     
      endif
      
      QS = ROIMP + SDRO + SSUR + SIF
      QG = BFS + BFP
      Q  = TCI
      DS = sactot - sactot0
      BAL = PXV-TET-QS-QG-DS
c      write(*,"(A,9(xf14.6))")'bal exsac1 ',BAL,PXV,TET,Q,QS,
c     &     QG,DS,sactot0,sactot
      
C     Convert SAC storages back to physical layers
      DWT=SACST(1)-SACST_PRV(1)
      DWF=SACST(2)-SACST_PRV(2)
      NUP=1
      NUPL=2
      if(prflag==1)write(*,*)'upper dwt dwf',dwt,dwf
      CALL SAC2FRZ1(DWT,DWF,SMC,SH2O,NUP,NUPL,ZSOIL,SMCMAX,SMCWLT,
     &     PRFLAG)


C  LOWER ZONE STORAGES
      DWT=SACST(3)-SACST_PRV(3)
      DWF=SACST(4)+SACST(5)-SACST_PRV(4)-SACST_PRV(5)
      NUP=NUPL+1
      if(prflag==1)write(*,*)'lower dwt dwf',dwt,dwf
      CALL SAC2FRZ1(DWT,DWF,SMC,SH2O,NUP,NSAC,ZSOIL,SMCMAX,SMCWLT,
     &     PRFLAG)

      temp1=0
      do i = 1,nsoil
         temp1 = temp1 + et1(i)
      enddo

C ----------------------------------------------------------------------
C CONVERT MOISTURE FLUXES BACK TO MM, UPDATE ET
C ----------------------------------------------------------------------

      EDIR1 = TET / (DT * 1000)
      RUNOFF1 = SURF / (DT * 1000)
      RUNOFF2 = GRND / (DT * 1000)
      ETA1 = ETA1 + EDIR1
      if (prflag==1) then
      write(*,*)'smflx edir1 ett eta1'
      write(*,*)edir1,temp1,eta1
      endif

      smctot=0.
      sh2otot=0.
      runtot=0.
      sactot=0.
      frztot=0.
      stot=0.
      sliq=0.
      runtot=(runoff1+runoff2)
      do i = 1,nsoil
         smctot = smctot + smc(i)
         sh2otot = sh2otot + sh2o(i)
      enddo  
      smcdif = smctot - smctot0
      sh2odif = sh2otot - sh2otot0
      if (prflag==1) then
      write(*,*)'-----------smflx BOTTOM SMFLX------' 
      write(*,*)'prcp           edir            ec'
      write(*,*)'runoff1        runoff2         runoff'
      write(*,*)runoff1,runoff2,runtot
      write(*,*)'smctot        sh2otot'
      write(*,*)smctot,sh2otot
      write(*,*)'smcdif        sh2odif'
      write(*,*)smcdif,sh2odif
      write(*,*)'---------------------------------' 
      endif
      do i = 1,5
         sactot = sactot + sacst(i)
         frztot = frztot + frzst(i+5)
      enddo
      stot = stot + 1000 * (smc(1) - smcwlt) * (0.0 - zsoil(1))
      sliq = sliq + 1000 * (sh2o(1) - smcwlt) * (0.0 - zsoil(1))
      if (prflag==1) write(*,*)'smctot',1,stot,smc(i)  
      do i = 2,nsoil
         stot = stot+1000*(smc(i) - smcwlt) * (zsoil(i-1) - zsoil(i))
         sliq = sliq+1000*(sh2o(i) - smcwlt) * (zsoil(i-1) - zsoil(i))
         if (prflag==1) write(*,*)'smctot',i,stot,smc(i)  
      enddo
      diftot = sactot - stot
      difliq = frztot - sliq
      if (prflag==1) then
      write(*,*)'stot sliq',stot,sliq
      write(*,*)'sactot frztot',sactot,frztot
      write(*,*)'diftot difliq',diftot,difliq
      endif

      DS = (stot - stot0) / 1000 
      DSH = (sliq - sliq0) / 1000
      BAL = (pcpdrp-edir1-ec1-tempet1-runtot)*dt-ds
      if(prflag==1)then
      write(*,"(A,9(xf16.9))") 'bal-smc ',BAL,pcpdrp*dt,edir1*dt,
     &     ec1*dt,tempet1*dt,runtot*dt,ds,dsh
      endif

C ----------------------------------------------------------------------
C STORE ICE CONTENT AT EACH SOIL LAYER BEFORE CALLING SRT & SSTEP
C ----------------------------------------------------------------------

      DO I = 1,NSOIL
        SICE(I) = SMC(I) - SH2O(I)
       END DO
   
C ----------------------------------------------------------------------
C END SUBROUTINE SMFLX_SAC
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE EXTRANSP (ET1,NSOIL,SMC,SH2O,SMCWLT,SACST,FRZST,DT,
     &     NUPL,NSAC,PRFLAG,CELLID,MSTEP,smcref,critsmc)
      
      IMPLICIT NONE
C ----------------------------------------------------------------------
C SUBROUTINE EXTRANSP
C ----------------------------------------------------------------------
C Remove root zone uptake (transpiration) and update SAC storages
C Make sure water is available.  Assume this moiture comes from tension
C water, as this is where other evap is extracted from SAC storages
C ----------------------------------------------------------------------
      INTEGER NSOIL,I,PRFLAG,NUPL,NSAC,CELLID,MSTEP
      REAL DWTUP,DWTLW
      REAL ET1(NSOIL)
      REAL ET1MM(NSOIL)
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SACST(6)
      REAL FRZST(10)
      REAL SMCWLT,DT
      REAL UZTWH,LZTWH,UZFWH,LZFPH,LZFSH
      REAL UZTWC,LZTWC,UZFWC,LZFPC,LZFSC,smcref,critsmc,fcrit

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE.
C ----------------------------------------------------------------------

C Enumerate transpired water to be removed now from SAC storages
C Since SAC-storages in mm, need to convert ET1 from M/s to MM
C over entire time step.

      DO I = 1,NSOIL
         fcrit = (smc(i) - smcwlt) / (critsmc - smcwlt)
         if(fcrit.lt.1.0.and.prflag==1)then
            write(*,*)'smc(i) fcrit',i,smc(i),fcrit
         endif
         ET1(I) = MIN(ET1(I),(ET1(I)*fcrit))
         ET1MM(I) = ET1(I) * 1000 * DT
      ENDDO

      UZTWH = FRZST(6)
      UZFWH = FRZST(7)
      UZTWC = SACST(1)
      UZFWC = SACST(2)
      
C Only remove free water for transpiration
      DO I = 1,NUPL
         IF (ET1MM(I).GT.0.0) THEN
            IF (PRFLAG==1) write(*,*)'ET1MM(I) UZTWH',ET1MM(I),UZTWH
c            IF (ET1MM(I).GT.UZTWH) THEN
            IF (ET1MM(I).GT.UZFWH) THEN
c               write(*,*)'error upper soil moisture threatened'
c               write(*,*)'cell',CELLID,'time',MSTEP,ET1MM(I)
c               write(*,*)'zone',I,' transpiration',ET1MM(I)
c               write(*,*)'uztwh uzfwh',uztwh,uzfwh
c               write(*,*)'lztwh lzfph lzfsh',lztwh,lzfph,lzfsh
               ET1MM(I) = 0.0
c     stop
            ELSE
C               UZTWH = UZTWH - ET1MM(I)
               UZFWH = UZFWH - ET1MM(I)
               UZFWC = UZFWC - ET1MM(I)
            ENDIF
         ENDIF
      ENDDO

      LZTWH = FRZST(8)
      LZFSH = FRZST(9)
      LZFPH = FRZST(10)
      LZTWC = SACST(3)
      LZFSC = SACST(4)
      LZFPC = SACST(5)

C Only remove free water for transpiration
      DO I = (NUPL+1),NSAC
         IF (ET1MM(I).GT.0.0) THEN
            IF (PRFLAG==1) write(*,*)'ET1MM(I) LZTWH',ET1MM(I),LZTWH
C            IF (ET1MM(I).GT.LZTWH) THEN
C            IF (ET1MM(I).GT.LZFPH) THEN
            IF ((ET1MM(I).GT.LZFPH).AND.(ET1MM(I).GT.LZFSH)) THEN
               ET1MM(I) = 0.0
c               write(*,*)'error lower soil moisture threatened'
c               write(*,*)'uztwh uzfwh',uztwh,uzfwh
c               write(*,*)'lztwh lzfph lzfsh',lztwh,lzfph,lzfsh
c     stop
            ELSE
C               LZTWH = LZTWH - ET1MM(I)
C               LZFPH = LZFPH - ET1MM(I)
C               LZFPC = LZFPC - ET1MM(I)
               IF (ET1MM(I).GT.LZFPH) THEN
                  LZFSH = LZFSH - ET1MM(I)
                  LZFSC = LZFSC - ET1MM(I)
               ELSE
                  LZFPH = LZFPH - ET1MM(I)
                  LZFPC = LZFPC - ET1MM(I)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

C Save changes to states in appropriate arrays

      FRZST(6) = UZTWH
      FRZST(7) = UZFWH
      SACST(1) = UZTWC
      SACST(2) = UZFWC
      FRZST(8) = LZTWH
      FRZST(10)= LZFPH
      SACST(3) = LZTWC
      SACST(5) = LZFPC

C Update transpiration if threatened wilting point

      DO I = 1,NSOIL
         ET1(I) = ET1MM(I) / (1000 * DT)
      ENDDO

cC Remove root zone water from transpiration so long as soil moisture does
cC not fall below wilting point
c
c      DO I = 1,NSOIL
c         IF (ET1(I).GT.0.0) THEN
c            IF (ET1(I).GT.(SH2O(I)-SMCWLT)) THEN
c               ET1(I) = MAX(0.0,(SH2O(I)-SMCWLT))
c               SH2O(I) = SMCWLT
c               SMC(I) = SMC(I) - ET1(I)
c            ELSE
c               SH2O(I) = SH2O(I) - ET1(I)
c               SMC(I) = SMC(I) - ET1(I)
c            ENDIF
c         ENDIF
cc         write(*,*)'i smc sh2o et1',i,smc(i),sh2o(i),et1(i)
c      ENDDO
c
cC Enumerate transpired water to be removed now from SAC storages
cC Since SAC-storages in mm, need to convert ET1 from M/s to MM
c      DWTUP = (ET1(1) + ET1(2)) * 1000 * DT
c      DWTLW = (ET1(3) + ET1(4)) * 1000 * DT
c
c      UZTWH = FRZST(6)
c      UZFWH = FRZST(7)
c      UZTWC = SACST(1)
c      UZFWC = SACST(2)
c      
c      IF (DWTUP.GT.0.0) THEN
c         IF (PRFLAG==1) write(*,*)'DWTUP UZTWH',DWTUP,UZTWH
c         IF (DWTUP.GT.UZTWH) THEN
c           DWTUP = DWTUP - UZTWH
c            UZTWH = 0
c            UZFWH = UZFWH - DWTUP
c            IF (UZFWH.LT.0) THEN
c               write(*,*)'error upper soil moisture exhausted'
c               write(*,*)'uztwh uzfwh',uztwh,uzfwh
c               write(*,*)'lztwh lzfph lzfsh',lztwh,lzfph,lzfsh
c               stop
c            ENDIF
c         ELSE
c            UZTWH = UZTWH - DWTUP
c         ENDIF
c      ENDIF
c
c      LZTWH = FRZST(8)
c      LZFPH = FRZST(10)
c      LZTWC = SACST(3)
c      LZFPC = SACST(5)
c      LZFSC = SACST(4)
c      LZFSH = FRZST(9)
c
c      IF (DWTLW.GT.0.0) THEN
c         IF(PRFLAG==1) write(*,*)'DWTLW LZTWH',DWTUP,UZTWH
c         IF (DWTLW.GT.LZTWH) THEN
c            DWTLW = DWTLW - LZTWH
c            LZTWH = 0
c            LZFPH = LZFPH - DWTLW
c            IF (LZFPH.LT.0) THEN
c               write(*,*)'error lower soil moisture exhausted'
c               write(*,*)'uztwh uzfwh',uztwh,uzfwh
c               write(*,*)'lztwh lzfph lzfsh',lztwh,lzfph,lzfsh
c               stop
c            ENDIF
c         ELSE
c            LZTWH = LZTWH - DWTLW
c         ENDIF
c      ENDIF
c
cC Save changes to states in appropriate arrays
c
c      FRZST(6) = UZTWH
c      FRZST(7) = UZFWH
c      SACST(1) = UZTWC
c      SACST(2) = UZFWC
c      FRZST(8) = LZTWH
c      FRZST(10)= LZFPH
c      SACST(3) = LZTWC
c      SACST(5) = LZFPC


C ----------------------------------------------------------------------
C SUBROUTINE EXTRANSP
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SMFLX (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &                 SH2O,SLOPE,KDT,FRZFACT,
     &                 SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &                 SHDFAC,CMCMAX,
     &                 RUNOFF1,RUNOFF2,RUNOFF3,
     &                 EDIR1,EC1,ET1,DRIP)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SMFLX
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER
C UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH
C PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
C FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND
C CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL BEXP
      REAL CMC
      REAL CMCMAX
      REAL DKSAT
      REAL DRIP
      REAL DT
      REAL DUMMY
      REAL DWSAT
      REAL EC1
      REAL EDIR1
      REAL ET1(NSOIL)
      REAL EXCESS
      REAL FRZFACT
      REAL KDT
      REAL PCPDRP
      REAL PRCP1
      REAL RHSCT
      REAL RHSTT(NSOLD)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL SHDFAC
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SICE(NSOLD)
      REAL SH2OA(NSOLD)
      REAL SH2OFG(NSOLD)
      REAL SLOPE
      REAL SMCMAX
      REAL SMCWLT
      REAL TRHSCT
      REAL ZSOIL(NSOIL)
	
      REAL FAC2	
      REAL FLIMIT

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE.
C ----------------------------------------------------------------------
      DUMMY = 0.

C ----------------------------------------------------------------------
C COMPUTE THE RIGHT HAND SIDE OF THE CANOPY EQN TERM ( RHSCT )
C ----------------------------------------------------------------------
      RHSCT = SHDFAC * PRCP1 - EC1

C ----------------------------------------------------------------------
C CONVERT RHSCT (A RATE) TO TRHSCT (AN AMOUNT) AND ADD IT TO EXISTING
C CMC.  IF RESULTING AMT EXCEEDS MAX CAPACITY, IT BECOMES DRIP AND WILL
C FALL TO THE GRND.
C ----------------------------------------------------------------------
      DRIP = 0.
      TRHSCT = DT * RHSCT
      EXCESS = CMC + TRHSCT
      IF (EXCESS .GT. CMCMAX) DRIP = EXCESS - CMCMAX

C ----------------------------------------------------------------------
C PCPDRP IS THE COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE
C SOIL
C ----------------------------------------------------------------------
      PCPDRP = (1. - SHDFAC) * PRCP1 + DRIP / DT

C ----------------------------------------------------------------------
C STORE ICE CONTENT AT EACH SOIL LAYER BEFORE CALLING SRT & SSTEP
C ----------------------------------------------------------------------
      DO I = 1,NSOIL
        SICE(I) = SMC(I) - SH2O(I)
      END DO
            
C ----------------------------------------------------------------------
C CALL SUBROUTINES SRT AND SSTEP TO SOLVE THE SOIL MOISTURE
C TENDENCY EQUATIONS. 
C
C IF THE INFILTRATING PRECIP RATE IS NONTRIVIAL,
C   (WE CONSIDER NONTRIVIAL TO BE A PRECIP TOTAL OVER THE TIME STEP 
C    EXCEEDING ONE ONE-THOUSANDTH OF THE WATER HOLDING CAPACITY OF 
C    THE FIRST SOIL LAYER)
C THEN CALL THE SRT/SSTEP SUBROUTINE PAIR TWICE IN THE MANNER OF 
C   TIME SCHEME "F" (IMPLICIT STATE, AVERAGED COEFFICIENT)
C   OF SECTION 2 OF KALNAY AND KANAMITSU (1988, MWR, VOL 116, 
C   PAGES 1945-1958)TO MINIMIZE 2-DELTA-T OSCILLATIONS IN THE 
C   SOIL MOISTURE VALUE OF THE TOP SOIL LAYER THAT CAN ARISE BECAUSE
C   OF THE EXTREME NONLINEAR DEPENDENCE OF THE SOIL HYDRAULIC 
C   DIFFUSIVITY COEFFICIENT AND THE HYDRAULIC CONDUCTIVITY ON THE
C   SOIL MOISTURE STATE
C OTHERWISE CALL THE SRT/SSTEP SUBROUTINE PAIR ONCE IN THE MANNER OF
C   TIME SCHEME "D" (IMPLICIT STATE, EXPLICIT COEFFICIENT) 
C   OF SECTION 2 OF KALNAY AND KANAMITSU
C PCPDRP IS UNITS OF KG/M**2/S OR MM/S, ZSOIL IS NEGATIVE DEPTH IN M 
C ----------------------------------------------------------------------
CY  According to Dr. Ken Mitchell's suggestion, add the second contraint  
CY  to remove numerical instability of runoff and soil moisture
CY  FLIMIT is a limit value for FAC2
	FAC2=0.0
	DO I=1,NSOIL
	FAC2=MAX(FAC2,SH2O(I)/SMCMAX)
	ENDDO

	CALL FAC2MIT(SMCMAX,FLIMIT)

	 IF (((PCPDRP*DT) .GT. (0.0001*1000.0*(-ZSOIL(1)))*SMCMAX).OR.
     &    (FAC2.GT.FLIMIT)) THEN
	
C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C SMC STATES REPLACED BY SH2O STATES IN SRT SUBR.  SH2O & SICE STATES
C INCLUDED IN SSTEP SUBR.  FROZEN GROUND CORRECTION FACTOR, FRZFACT
C ADDED.  ALL WATER BALANCE CALCULATIONS USING UNFROZEN WATER
C ----------------------------------------------------------------------
        CALL SRT (RHSTT,EDIR1,ET1,SH2O,SH2O,NSOIL,PCPDRP,ZSOIL,
     &            DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1, 
     &            RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZFACT,SICE,AI,BI,CI)
         
        CALL SSTEP (SH2OFG,SH2O,DUMMY,RHSTT,RHSCT,DT,NSOIL,SMCMAX,
     &              CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,AI,BI,CI)
         
        DO K = 1,NSOIL
          SH2OA(K) = (SH2O(K) + SH2OFG(K)) * 0.5
        END DO
        
        CALL SRT (RHSTT,EDIR1,ET1,SH2O,SH2OA,NSOIL,PCPDRP,ZSOIL,
     &            DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1,
     &            RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZFACT,SICE,AI,BI,CI)
         
        CALL SSTEP (SH2O,SH2O,CMC,RHSTT,RHSCT,DT,NSOIL,SMCMAX,
     &              CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,AI,BI,CI)
         
      ELSE
         
        CALL SRT (RHSTT,EDIR1,ET1,SH2O,SH2O,NSOIL,PCPDRP,ZSOIL,
     &            DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1,
     &            RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZFACT,SICE,AI,BI,CI)

        CALL SSTEP (SH2O,SH2O,CMC,RHSTT,RHSCT,DT,NSOIL,SMCMAX,
     &              CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,AI,BI,CI)
         
      ENDIF
      
c      RUNOF = RUNOFF

C ----------------------------------------------------------------------
C END SUBROUTINE SMFLX
C ----------------------------------------------------------------------
      RETURN
      END

C ----------------------------------------------------------------------
      SUBROUTINE FAC2MIT(SMCMAX,FLIMIT)
      IMPLICIT NONE		
      REAL SMCMAX,FLIMIT

	 FLIMIT=0.90

      IF(SMCMAX.EQ.0.395) THEN
	FLIMIT=0.59
           ELSE IF(SMCMAX.EQ.0.434.OR.SMCMAX.EQ.0.404) THEN
    	 FLIMIT=0.85
	        ELSE IF(SMCMAX.EQ.0.465.OR.SMCMAX.EQ.0.406) THEN
    	 FLIMIT=0.86
 	    ELSE IF(SMCMAX.EQ.0.476.OR.SMCMAX.EQ.0.439) THEN
    	 FLIMIT=0.74
               ELSE IF(SMCMAX.EQ.0.200.OR.SMCMAX.EQ.0.464) THEN
	FLIMIT=0.80
          ENDIF

	RETURN
	END
     
C ----------------------------------------------------------------------
C END SUBROUTINE FAC2MIT
C ----------------------------------------------------------------------


      SUBROUTINE SNFRAC (SNEQV,SNUP,SALP,SNOWH,SNCOVR)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SNFRAC
C ----------------------------------------------------------------------
C CALCULATE SNOW FRACTION (0 -> 1)
C SNEQV   SNOW WATER EQUIVALENT (M)
C SNUP    THRESHOLD SNEQV DEPTH ABOVE WHICH SNCOVR=1
C SALP    TUNING PARAMETER
C SNCOVR  FRACTIONAL SNOW COVER
C ----------------------------------------------------------------------
      REAL SNEQV, SNUP, SALP, SNCOVR, RSNOW, Z0N, SNOWH
      
C ----------------------------------------------------------------------
C SNUP IS VEG-CLASS DEPENDENT SNOWDEPTH THRESHHOLD (SET IN ROUTINE
C REDPRM) ABOVE WHICH SNCOVR=1.
C ----------------------------------------------------------------------
          IF (SNEQV .LT. SNUP) THEN
            RSNOW = SNEQV/SNUP
            SNCOVR = 1. - ( EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))
          ELSE
            SNCOVR = 1.0
          ENDIF

          Z0N=0.035 
C     FORMULATION OF DICKINSON ET AL. 1986

C        SNCOVR=SNOWH/(SNOWH + 5*Z0N)

C     FORMULATION OF MARSHALL ET AL. 1994
C        SNCOVR=SNEQV/(SNEQV + 2*Z0N)

C ----------------------------------------------------------------------
C END SUBROUTINE SNFRAC
C ----------------------------------------------------------------------
      RETURN
      END

      FUNCTION SNKSRC (TAVG,SMC,SH2O,ZSOIL,NSOIL,
     &                 SMCMAX,PSISAT,BEXP,DT,K,QTOT) 
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C FUNCTION SNKSRC
C ----------------------------------------------------------------------
C CALCULATE SINK/SOURCE TERM OF THE TERMAL DIFFUSION EQUATION. (SH2O) IS
C AVAILABLE LIQUED WATER.
C ----------------------------------------------------------------------
      INTEGER K
      INTEGER NSOIL
      
      REAL BEXP
      REAL DF
      REAL DH2O
      REAL DT
      REAL DZ
      REAL DZH
      REAL FREE
      REAL FRH2O
      REAL HLICE
      REAL PSISAT
      REAL QTOT
      REAL SH2O
      REAL SMC
      REAL SMCMAX
      REAL SNKSRC
      REAL T0
      REAL TAVG
      REAL TDN
      REAL TM
      REAL TUP
      REAL TZ
      REAL X0
      REAL XDN
      REAL XH2O
      REAL XUP
      REAL ZSOIL (NSOIL)

      PARAMETER(DH2O = 1.0000E3)
      PARAMETER(HLICE = 3.3350E5)
      PARAMETER(T0 = 2.7315E2)
      
      IF (K .EQ. 1) THEN
        DZ = -ZSOIL(1)
      ELSE
        DZ = ZSOIL(K-1)-ZSOIL(K)
      ENDIF

C ----------------------------------------------------------------------
C VIA FUNCTION FRH2O, COMPUTE POTENTIAL OR 'EQUILIBRIUM' UNFROZEN
C SUPERCOOLED FREE WATER FOR GIVEN SOIL TYPE AND SOIL LAYER TEMPERATURE.
C FUNCTION FRH20 INVOKES EQN (17) FROM V. KOREN ET AL (1999, JGR, VOL.
C 104, PG 19573).  (ASIDE:  LATTER EQN IN JOURNAL IN CENTIGRADE UNITS.
C ROUTINE FRH2O USE FORM OF EQN IN KELVIN UNITS.)
C ----------------------------------------------------------------------
      FREE = FRH2O(TAVG,SMC,SH2O,SMCMAX,BEXP,PSISAT)

C ----------------------------------------------------------------------
C IN NEXT BLOCK OF CODE, INVOKE EQN 18 OF V. KOREN ET AL (1999, JGR,
C VOL. 104, PG 19573.)  THAT IS, FIRST ESTIMATE THE NEW AMOUNTOF LIQUID
C WATER, 'XH2O', IMPLIED BY THE SUM OF (1) THE LIQUID WATER AT THE BEGIN
C OF CURRENT TIME STEP, AND (2) THE FREEZE OF THAW CHANGE IN LIQUID
C WATER IMPLIED BY THE HEAT FLUX 'QTOT' PASSED IN FROM ROUTINE HRT.
C SECOND, DETERMINE IF XH2O NEEDS TO BE BOUNDED BY 'FREE' (EQUIL AMT) OR
C IF 'FREE' NEEDS TO BE BOUNDED BY XH2O.
C ----------------------------------------------------------------------
      XH2O = SH2O + QTOT*DT/(DH2O*HLICE*DZ)

C ----------------------------------------------------------------------
C FIRST, IF FREEZING AND REMAINING LIQUID LESS THAN LOWER BOUND, THEN
C REDUCE EXTENT OF FREEZING, THEREBY LETTING SOME OR ALL OF HEAT FLUX
C QTOT COOL THE SOIL TEMP LATER IN ROUTINE HRT.
C ----------------------------------------------------------------------
      IF ( XH2O .LT. SH2O .AND. XH2O .LT. FREE) THEN 
        IF ( FREE .GT. SH2O ) THEN
          XH2O = SH2O
        ELSE
          XH2O = FREE
        ENDIF
      ENDIF
              
C ----------------------------------------------------------------------
C SECOND, IF THAWING AND THE INCREASE IN LIQUID WATER GREATER THAN UPPER
C BOUND, THEN REDUCE EXTENT OF THAW, THEREBY LETTING SOME OR ALL OF HEAT
C FLUX QTOT WARM THE SOIL TEMP LATER IN ROUTINE HRT.
C ----------------------------------------------------------------------
      IF ( XH2O .GT. SH2O .AND. XH2O .GT. FREE )  THEN
        IF ( FREE .LT. SH2O ) THEN
          XH2O = SH2O
        ELSE
          XH2O = FREE
        ENDIF
      ENDIF 

      IF (XH2O .LT. 0.) XH2O = 0.
      IF (XH2O .GT. SMC) XH2O = SMC

C ----------------------------------------------------------------------
C CALCULATE PHASE-CHANGE HEAT SOURCE/SINK TERM FOR USE IN ROUTINE HRT
C AND UPDATE LIQUID WATER TO REFLCET FINAL FREEZE/THAW INCREMENT.
C ----------------------------------------------------------------------
      SNKSRC = -DH2O*HLICE*DZ*(XH2O-SH2O)/DT
      SH2O = XH2O
      
C ----------------------------------------------------------------------
C END FUNCTION SNKSRC
C ----------------------------------------------------------------------
77    RETURN
      END
      SUBROUTINE SNOPAC (ETP,ETA,PRCP,PRCP1,SNOWNG,SMC,SMCMAX,SMCWLT,
     &     SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &     SBETA,DF1,
     &     Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &     SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,ESD,SNDENS,
     &     SNOWH,SH2O,SLOPE,KDT,FRZFACT,PSISAT,SNUP,
     &     ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &     RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &     ICE,RTDIS,QUARTZ,FXEXP,FX,CSOIL,
     &     BETA,DRIP,DEW,FLX1,FLX2,FLX3,ESNOW,EMISS,EMISSNOW,MODEL_TYPE,
     &     FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP,wcrit)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SNOPAC
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES & UPDATE SOIL MOISTURE
C CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN A SNOW PACK IS
C PRESENT.
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER NROOT
      INTEGER NSOIL

      LOGICAL SNOWNG

      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL CP
      REAL CPH2O
      REAL CPICE
      REAL CSOIL
      REAL DENOM
      REAL DEW
      REAL DF1
      REAL DKSAT
      REAL DRIP
      REAL DSOIL
      REAL DTOT
      REAL DT
      REAL DWSAT
      REAL EC
      REAL EDIR
      REAL EPSCA
      REAL ESD
      REAL ESDMIN
      REAL EXPSNO
      REAL EXPSOI
      REAL ETA
      REAL ETA1
      REAL ETP
      REAL ETP1
      REAL ETP2
      REAL ET(NSOIL)
      REAL ETT
      REAL EX
      REAL EXPFAC
      REAL FDOWN
      REAL FXEXP
      REAL FX
      REAL FLX1
      REAL FLX2
      REAL FLX3
      REAL F1
      REAL KDT
      REAL LSUBF
      REAL LSUBC
      REAL LSUBS
      REAL PC
      REAL PRCP
      REAL PRCP1
      REAL Q2
      REAL RCH
      REAL RR
      REAL RTDIS(NSOIL)
      REAL SSOIL
      REAL SBETA
      REAL SSOIL1
      REAL SFCTMP
      REAL SHDFAC
      REAL SIGMA
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL SNOMLT
      REAL SNOWH
      REAL STC(NSOIL)
      REAL T1
      REAL T11
      REAL T12
      REAL T12A
      REAL T12B
      REAL T24
      REAL TBOT
      REAL ZBOT
      REAL TH2
      REAL YY
      REAL ZSOIL(NSOIL)
      REAL ZZ1
      REAL TFREEZ
      REAL SALP
      REAL SFCPRS
      REAL SLOPE
      REAL FRZFACT
      REAL PSISAT
      REAL SNUP
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL QUARTZ
      REAL SNDENS
      REAL SNCOND
      REAL RSNOW
      REAL SNCOVR
      REAL EMISS
      REAL EMISSNOW
      REAL QSAT
      REAL ETP3
      REAL SEH
      REAL T14
      REAL CSNOW

      REAL EC1
      REAL EDIR1
      REAL ET1(NSOIL)
      REAL ETT1

      REAL ETNS
      REAL ETNS1
      REAL ESNOW
      REAL ESNOW1
      REAL ESNOW2
      REAL ETANRG

      REAL SNOEXP

      INTEGER K
      INTEGER MODEL_TYPE
      REAL FRZST(10)
      REAL FRZPAR(13)
      REAL SACST(6)
      REAL SACPAR(16)
      REAL FROST,wcrit
      INTEGER PRFLAG,CELLID,MSTEP,SMCFLAG

      PARAMETER(CP = 1004.5)
      PARAMETER(CPH2O = 4.218E+3)
      PARAMETER(CPICE = 2.106E+3)
      PARAMETER(ESDMIN = 1.E-6)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(LSUBC = 2.501000E+6)
      PARAMETER(LSUBS = 2.83E+6)
      PARAMETER(SIGMA = 5.67E-8)
      PARAMETER(TFREEZ = 273.15)

C     DATA SNOEXP /1.0/
      DATA SNOEXP /2.0/

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C CONVERT POTENTIAL EVAP (ETP) FROM KG M-2 S-1 TO M S-1 AND THEN TO AN
C AMOUNT (M) GIVEN TIMESTEP (DT) AND CALL IT AN EFFECTIVE SNOWPACK
C REDUCTION AMOUNT, ETP2 (M).  THIS IS THE AMOUNT THE SNOWPACK WOULD BE
C REDUCED DUE TO EVAPORATION FROM THE SNOW SFC DURING THE TIMESTEP.
C EVAPORATION WILL PROCEED AT THE POTENTIAL RATE UNLESS THE SNOW DEPTH
C IS LESS THAN THE EXPECTED SNOWPACK REDUCTION.
C IF SEAICE (ICE=1), BETA REMAINS=1.
C ----------------------------------------------------------------------
      PRCP1 = PRCP1*0.001

c      ETP2 = ETP * 0.001Q * DT
      BETA = 1.0
      IF (ICE .NE. 1) THEN
        IF (ESD .LT. ETP2) THEN
c          BETA = ESD / ETP2
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
      EDIR = 0.0
      EDIR1 = 0.0
      EC = 0.0
      EC1 = 0.0
      DO K = 1,NSOIL
        ET(K) = 0.0
        ET1(K) = 0.0
      ENDDO
      ETT = 0.0
      ETT1 = 0.0
      ETNS = 0.0
      ETNS1 = 0.0
      ESNOW = 0.0
      ESNOW1 = 0.0
      ESNOW2 = 0.0
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C IF ETP<0 (DOWNWARD) THEN DEWFALL (=FROSTFALL IN THIS CASE).
C ----------------------------------------------------------------------
      DEW = 0.0
      ETP1 = ETP*0.001
      IF (ETP .LT. 0.0) THEN
c        DEW = -ETP * 0.001
        DEW = -ETP1
c        ESNOW2 = ETP * 0.001 * DT
        ESNOW2 = ETP1 * DT
        ETANRG = ETP*((1.-SNCOVR)*LSUBC + SNCOVR*LSUBS)
c      ENDIF
      ELSE
C ----------------------------------------------------------------------
c      ETP1 = 0.0
c        ETP1 = ETP*0.001
        IF (ICE .NE. 1) THEN
          IF (SNCOVR .LT. 1.) THEN
c          CALL EVAPO (ETA1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
            CALL EVAPO (ETNS1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
     &            SH2O,STC,SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &            SMCREF,SHDFAC,CMCMAX,
     &            SMCDRY,CFACTR,
     &            EDIR1,EC1,ET1,ETT1,SFCTMP,Q2,NROOT,RTDIS,
     &            FXEXP,FX,MODEL_TYPE,SMCFLAG)
c        ENDIF
C ----------------------------------------------------------------------
            EDIR1 = EDIR1*(1.-SNCOVR)
            EC1 = EC1*(1.-SNCOVR)
            DO K = 1,NSOIL
              ET1(K) = ET1(K)*(1.-SNCOVR)
            END DO
            ETT1 = ETT1*(1.-SNCOVR)
            ETNS1 = ETNS1*(1.-SNCOVR)
C ----------------------------------------------------------------------
            EDIR = EDIR1 * 1000.0
            EC = EC1 * 1000.0
            DO K = 1,NSOIL
              ET(K) = ET1(K) * 1000.0
            END DO
            ETT = ETT1 * 1000.0
            ETNS = ETNS1 * 1000.0
C ----------------------------------------------------------------------
          ENDIF
          ESNOW = ETP*SNCOVR
c          ESNOW1 = ETP*0.001
          ESNOW1 = ESNOW*0.001
          ESNOW2 = ESNOW1*DT
          ETANRG = ESNOW*LSUBS + ETNS*LSUBC
        ENDIF
      ENDIF
   
C ----------------------------------------------------------------------
C IF PRECIP IS FALLING, CALCULATE HEAT FLUX FROM SNOW SFC TO NEWLY
C ACCUMULATING PRECIP.  NOTE THAT THIS REFLECTS THE FLUX APPROPRIATE FOR
C THE NOT-YET-UPDATED SKIN TEMPERATURE (T1).  ASSUMES TEMPERATURE OF THE
C SNOWFALL STRIKING THE GOUND IS =SFCTMP (LOWEST MODEL LEVEL AIR TEMP).
C ----------------------------------------------------------------------
      FLX1 = 0.0
      IF (SNOWNG) THEN
        FLX1 = CPICE * PRCP * (T1 - SFCTMP)
      ELSE
        IF (PRCP .GT. 0.0) FLX1 = CPH2O * PRCP * (T1 - SFCTMP)
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE AN 'EFFECTIVE SNOW-GRND SFC TEMP' (T12) BASED ON HEAT FLUXES
C BETWEEN THE SNOW PACK AND THE SOIL AND ON NET RADIATION.
C INCLUDE FLX1 (PRECIP-SNOW SFC) AND FLX2 (FREEZING RAIN LATENT HEAT)
C FLUXES.
C FLX2 REFLECTS FREEZING RAIN LATENT HEAT FLUX USING T1 CALCULATED IN
C PENMAN.
C ----------------------------------------------------------------------
      DSOIL = -(0.5 * ZSOIL(1))
      DTOT = SNOWH + DSOIL
      DENOM = 1.0 + DF1 / (DTOT * RR * RCH)
c      T12A = ( (FDOWN-FLX1-FLX2-SIGMA*T24)/RCH
c     &       + TH2 - SFCTMP - BETA*EPSCA ) / RR
C      T12A = ( (FDOWN-FLX1-FLX2-SIGMA*T24)/RCH
C     &       + TH2 - SFCTMP - ETANRG/RCH ) / RR
      T12A = ( (FDOWN-FLX1-FLX2-
     &       (EMISSNOW*SNCOVR+EMISS*(1.0-SNCOVR))*SIGMA*T24)/RCH
     &       + TH2 - SFCTMP - ETANRG/RCH ) / RR

      T12B = DF1 * STC(1) / (DTOT * RR * RCH)
      T12 = (SFCTMP + T12A + T12B) / DENOM      

C ----------------------------------------------------------------------
C IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS AT OR BELOW FREEZING, NO SNOW
C MELT WILL OCCUR.  SET THE SKIN TEMP TO THIS EFFECTIVE TEMP.  REDUCE
C (BY SUBLIMINATION ) OR INCREASE (BY FROST) THE DEPTH OF THE SNOWPACK,
C DEPENDING ON SIGN OF ETP.
C UPDATE SOIL HEAT FLUX (SSOIL) USING NEW SKIN TEMPERATURE (T1)
C SINCE NO SNOWMELT, SET ACCUMULATED SNOWMELT TO ZERO, SET 'EFFECTIVE'
C PRECIP FROM SNOWMELT TO ZERO, SET PHASE-CHANGE HEAT FLUX FROM SNOWMELT
C TO ZERO.
C ----------------------------------------------------------------------
      IF (T12 .LE. TFREEZ) THEN
        T1 = T12
        SSOIL = DF1 * (T1 - STC(1)) / DTOT
c        ESD = MAX(0.0, ESD-ETP2)
        ESD = MAX(0.0, ESD-ESNOW2)
        FLX3 = 0.0
        EX = 0.0
        SNOMLT = 0.0

      ELSE
C ----------------------------------------------------------------------
C IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS ABOVE FREEZING, SNOW MELT
C WILL OCCUR.  CALL THE SNOW MELT RATE,EX AND AMT, SNOMLT.  REVISE THE
C EFFECTIVE SNOW DEPTH.  REVISE THE SKIN TEMP BECAUSE IT WOULD HAVE CHGD
C DUE TO THE LATENT HEAT RELEASED BY THE MELTING. CALC THE LATENT HEAT
C RELEASED, FLX3. SET THE EFFECTIVE PRECIP, PRCP1 TO THE SNOW MELT RATE,
C EX FOR USE IN SMFLX.  ADJUSTMENT TO T1 TO ACCOUNT FOR SNOW PATCHES.
C CALCULATE QSAT VALID AT FREEZING POINT.  NOTE THAT ESAT (SATURATION
C VAPOR PRESSURE) VALUE OF 6.11E+2 USED HERE IS THAT VALID AT FRZZING
C POINT.  NOTE THAT ETP FROM CALL PENMAN IN SFLX IS IGNORED HERE IN
C FAVOR OF BULK ETP OVER 'OPEN WATER' AT FREEZING TEMP.
C UPDATE SOIL HEAT FLUX (S) USING NEW SKIN TEMPERATURE (T1)
C ----------------------------------------------------------------------
C        T1 = TFREEZ * SNCOVR + T12 * (1.0 - SNCOVR)
C mek Feb2004
C non-linear weighting of snow vs non-snow covered portions of gridbox
C so with SNOEXP = 2.0 (>1), surface skin temperature is higher than for
C the linear case (SNOEXP = 1).
 
        T1 = TFREEZ * SNCOVR**SNOEXP + T12 * (1.0 - SNCOVR**SNOEXP) 

C        QSAT = (0.622*6.11E2)/(SFCPRS-0.378*6.11E2)
C        ETP = RCH*(QSAT-Q2)/CP
C        ETP2 = ETP*0.001*DT
        BETA = 1.0
        SSOIL = DF1 * (T1 - STC(1)) / DTOT
	
C ----------------------------------------------------------------------
C IF POTENTIAL EVAP (SUBLIMATION) GREATER THAN DEPTH OF SNOWPACK.
C BETA<1
C SNOWPACK HAS SUBLIMATED AWAY, SET DEPTH TO ZERO.
C ----------------------------------------------------------------------
c        IF (ESD .LE. ETP2) THEN
c        IF (ESD .LE. ESNOW2) THEN

        IF (ESD-ESNOW2 .LE. ESDMIN) THEN
           FLX3=0.0
c          BETA = ESD / ETP2
           ESD = 0.0
           EX = 0.0
           SNOMLT = 0.0

        ELSE
C ----------------------------------------------------------------------
C POTENTIAL EVAP (SUBLIMATION) LESS THAN DEPTH OF SNOWPACK, RETAIN
C   BETA=1.
C SNOWPACK (ESD) REDUCED BY POTENTIAL EVAP RATE
C ETP3 (CONVERT TO FLUX)
C ----------------------------------------------------------------------
c     ESD = ESD-ETP2
           ESD = ESD-ESNOW2
c     ETP3 = ETP*LSUBC
           SEH = RCH*(T1-TH2)
           T14 = T1*T1
           T14 = T14*T14
c     FLX3 = FDOWN - FLX1 - FLX2 - SIGMA*T14 - SSOIL - SEH - ETP3
C     FLX3 = FDOWN - FLX1 - FLX2 - SIGMA*T14 - SSOIL - SEH - ETANRG
           FLX3 = FDOWN - FLX1 - FLX2 -
     &          (EMISSNOW*SNCOVR+EMISS*(1.0-SNCOVR))*SIGMA*T14
     &          - SSOIL - SEH - ETANRG
           
           IF (FLX3 .LE .0.0) FLX3 = 0.0
           EX = FLX3*0.001/LSUBF
           
C     ----------------------------------------------------------------------
C     SNOWMELT REDUCTION DEPENDING ON SNOW COVER
C     IF SNOW COVER LESS THAN 5% NO SNOWMELT REDUCTION
C     ***NOTE:  DOES 'IF' BELOW FAIL TO MATCH THE MELT WATER WITH THE MELT
C     ENERGY?
C     ----------------------------------------------------------------------
c     IF (SNCOVR .GT. 0.05) EX = EX * SNCOVR
           SNOMLT = EX * DT
           
C     ----------------------------------------------------------------------
C     ESDMIN REPRESENTS A SNOWPACK DEPTH THRESHOLD VALUE BELOW WHICH WE
C     CHOOSE NOT TO RETAIN ANY SNOWPACK, AND INSTEAD INCLUDE IT IN SNOWMELT.
C     ----------------------------------------------------------------------
           IF (ESD-SNOMLT .GE. ESDMIN) THEN
              ESD = ESD - SNOMLT
              
           ELSE
C     ----------------------------------------------------------------------
C     SNOWMELT EXCEEDS SNOW DEPTH
C     ----------------------------------------------------------------------
              EX = ESD/DT
              FLX3 = EX*1000.0*LSUBF
              SNOMLT = ESD
              ESD = 0.0
              
           ENDIF
C     ----------------------------------------------------------------------
C     END OF 'ESD .LE. ETP2' IF-BLOCK
C     ----------------------------------------------------------------------
        ENDIF
        
        PRCP1 = PRCP1 + EX
        
C     ----------------------------------------------------------------------
C     END OF 'T12 .LE. TFREEZ' IF-BLOCK
C     ----------------------------------------------------------------------
      ENDIF
      
C ----------------------------------------------------------------------
C FINAL BETA NOW IN HAND, SO COMPUTE EVAPORATION.  EVAP EQUALS ETP
C UNLESS BETA<1.
C ----------------------------------------------------------------------
c      ETA = BETA*ETP

C ----------------------------------------------------------------------
C SET THE EFFECTIVE POTNL EVAPOTRANSP (ETP1) TO ZERO SINCE THIS IS SNOW
C CASE, SO SURFACE EVAP NOT CALCULATED FROM EDIR, EC, OR ETT IN SMFLX
C (BELOW).
C IF SEAICE (ICE=1) SKIP CALL TO SMFLX.
C SMFLX RETURNS UPDATED SOIL MOISTURE VALUES.  IN THIS, THE SNOW PACK
C CASE, ETA1 IS NOT USED IN CALCULATION OF EVAP.
C ----------------------------------------------------------------------
c      ETP1 = 0.0
      IF (ICE .NE. 1) THEN
c     CALL EVAPO (ETA1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
c     &              SH2O,SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
c     &              SMCREF,SHDFAC,CMCMAX,
c     &              SMCDRY,CFACTR,
c     &              EDIR1,EC1,ET1,ETT1,SFCTMP,Q2,NROOT,RTDIS,FXEXP)
         IF (MODEL_TYPE.EQ.0) THEN
            CALL SMFLX (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &           SH2O,SLOPE,KDT,FRZFACT,
     &           SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &           SHDFAC,CMCMAX,
     &           RUNOFF1,RUNOFF2,RUNOFF3,
     &           EDIR1,EC1,ET1,
     &           DRIP)
         ENDIF
         IF (MODEL_TYPE.EQ.1) THEN
            CALL SMFLX_SAC (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &           SH2O,SLOPE,KDT,FRZFACT,SNCOVR,SMCREF,
     &           SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &           SHDFAC,CMCMAX,SFCTMP,ETP1,ETA1,
     &           RUNOFF1,RUNOFF2,RUNOFF3,EDIR1,EC1,ET1,DRIP, 
     &           FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP,wcrit)
         ENDIF
      ENDIF

C ----------------------------------------------------------------------
c        EDIR = EDIR1 * 1000.0
c        EC = EC1 * 1000.0
c        ETT = ETT1 * 1000.0
c        ET(1) = ET1(1) * 1000.0
c        ET(2) = ET1(2) * 1000.0
c        ET(3) = ET1(3) * 1000.0
c        ET(4) = ET1(4) * 1000.0
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C BEFORE CALL SHFLX IN THIS SNOWPACK CASE, SET ZZ1 AND YY ARGUMENTS TO
C SPECIAL VALUES THAT ENSURE THAT GROUND HEAT FLUX CALCULATED IN SHFLX
C MATCHES THAT ALREADY COMPUTER FOR BELOW THE SNOWPACK, THUS THE SFC
C HEAT FLUX TO BE COMPUTED IN SHFLX WILL EFFECTIVELY BE THE FLUX AT THE
C SNOW TOP SURFACE.  T11 IS A DUMMY ARGUEMENT SO WE WILL NOT USE THE
C SKIN TEMP VALUE AS REVISED BY SHFLX.
C ----------------------------------------------------------------------
      ZZ1 = 1.0
      YY = STC(1)-0.5*SSOIL*ZSOIL(1)*ZZ1/DF1
      T11 = T1

C ----------------------------------------------------------------------
C SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX 
C (SSOIL1) AND THE SKIN TEMP (T11) OUTPUT FROM THIS SHFLX CALL ARE NOT
C USED  IN ANY SUBSEQUENT CALCULATIONS. RATHER, THEY ARE DUMMY VARIABLES
C HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT FLUX ARE
C UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
C ----------------------------------------------------------------------
      CALL SHFLX (SSOIL1,STC,SMC,SMCMAX,NSOIL,T11,DT,YY,ZZ1,ZSOIL,
     &            TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &            QUARTZ,CSOIL)
      
C ----------------------------------------------------------------------
C SNOW DEPTH AND DENSITY ADJUSTMENT BASED ON SNOW COMPACTION.  YY IS
C ASSUMED TO BE THE SOIL TEMPERTURE AT THE TOP OF THE SOIL COLUMN.
C ----------------------------------------------------------------------
      IF (ESD .GT. 0.) THEN
        CALL SNOWPACK (ESD,DT,SNOWH,SNDENS,T1,YY)
      ELSE
        ESD = 0.
        SNOWH = 0.
        SNDENS = 0.
        SNCOND = 1.
        SNCOVR = 0.
      ENDIF

C ----------------------------------------------------------------------
C END SUBROUTINE SNOPAC
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SNWPAC (ETP,ETA,PRCP,PRCP1,SNOWNG,FRZGRA,SMC,SMCMAX,
     &     SMCWLT,SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &     SBETA,DF1,PACH20,CH,SFCSPD,
     &     Q2,T1,TPACK,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &     SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,ESD,SNDENS,
     &     SNOWH,SH2O,SLOPE,KDT,FRZFACT,PSISAT,SNUP,
     &     ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &     RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &     ICE,RTDIS,QUARTZ,FXEXP,FX,CSOIL,BETA,DRIP,
     &     DEW,FLX1,FLX2,FLX3,ESNOW,EMISS,EMISSNOW,MODEL_TYPE,
     &     FRZST,FRZPAR,SACST,SACPAR,ERFLAG,PRFLAG,CELLID,MSTEP,
     &     LVRAIN,wcrit)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SNWPAC
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES & UPDATE SOIL MOISTURE
C CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN A SNOW PACK IS
C PRESENT.
C
C Modifications:
C 2008-Aug-11 Replaced TFREEZ with T1 in initial guess for T brackets.	Ben Livneh
C 2008-Aug-12 Added ERFLAG to report errors to parent function.         Ben Livneh
C 2008-Sep-08 Changed initial lower bound on T1 iteration to be
C             relative to T1_0, the previous time step's final T1.      TJB
C 2008-Sep-20 Added condition on TLOW to prevent it from being higher
C             than TFREEZ-BOUND.					TJB
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER NROOT
      INTEGER NSOIL
      INTEGER LSTSNW
      INTEGER ITRMAX
      INTEGER ITER
      INTEGER K
      INTEGER J,I
      INTEGER ERFLAG

      LOGICAL SNOWNG
      LOGICAL FRZGRA
      LOGICAL LVRAIN

      REAL BEXP,BETA,CC,CFACTR,CMC,CMCMAX,CHICE,CP,CPH2O,CPICE,CSOIL,ESD
      REAL DENOM,DEW,DF1,DKSAT,DRIP,DSOIL,DTOT,DT,DWSAT,EC,EDIR,EPSCA
      REAL ESDMIN,EXPSNO,EXPSNOI,ETA,ETA1,ETP,ETP1,ETP2,ET(NSOIL),SNOMLT
      REAL ETT,EX,EXPFAC,FX,EMISS,EMISSNOW,FDOWN,FXEXP,FLX1,FLX2,FLX3
      REAL F1,KDT,LSUBF,LSUBC,LSUBS,PC,PRCP,PRCP1,Q2,RCH,RR,RTDIS(NSOIL)
      REAL SSOIL,SBETA,SSOIL1,SFCTMP,SHDFAC,SIGMA,SMC(NSOIL),SH2O(NSOIL)
      REAL SMCDRY,SMCMAX,SMCREF,SMCWLT,SNOWLT,SNOWH,STC(NSOIL),T12,T14
      REAL T1,T1A,T1B,TH2,T2,TBOT,ZBOT,YY,ZSOIL(NSOIL),ZZ1,TFREEZ,SALP
      REAL SFCPRS,SLOPE,FRZFACT,PSISAT,SNUP,RUNOFF1,RUNOFF2,RUNOFF3,T11
      REAL QUARTZ,SNDENS,SNCOND,RSNOW,SNCOVR,QSAT,ETP3,SEH,CSNOW,EC1
      REAL EDIR1,ET1(NSOIL),ETT1,ETNS,ETNS1,ESNOW,ESNOW1,ESNOW2,ETANRG
      REAL ESSNOW
      REAL LHF
      REAL QNET
      REAL FRCH20
      REAL IRRSAT
      REAL EACT
      REAL VMF
      REAL AIRDEN
      REAL CH
      REAL SFCSPD
      REAL CORECT
      REAL TOPSWE
      REAL PACSWE
      REAL MAXSWE
      REAL TOPCC
      REAL PACCC
      REAL TOPH20
      REAL PACH20
      REAL MAXH20
      REAL RAIN
      REAL SNOW
      REAL SNOWCC
      REAL DELSWE
      REAL DELCC
      REAL RFZNRG
      REAL RFZH20
      REAL MLTNRG
      REAL A
      REAL B
      REAL C
      REAL D
      REAL E
      REAL P
      REAL Q
      REAL R
      REAL S
      REAL FA
      REAL FB
      REAL FC
      REAL T1HIGH
      REAL T1LOW
      REAL TOL
      REAL TOL1
      REAL XM
      REAL EPS
      REAL BOUNDS
      REAL T12A
      REAL T12B
      REAL T24
      REAL TPACK
      REAL T1_0
      REAL FRZST(10)
      REAL FRZPAR(13)
      REAL SACST(6)
      REAL SACPAR(16)
      REAL FROST
      REAL SNOEXP
      REAL DEWNS
      REAL DEWS
      REAL ettemp
      REAL PRCP0
      REAL PRCP2
      INTEGER MODEL_TYPE,PRFLAG,CELLID,MSTEP,SMCFLAG
      real swe0,swe,swe1,swedif,bal,swe2,wcrit

      PARAMETER(CP = 1004.5)
      PARAMETER(CPH2O = 4.218E+3)
      PARAMETER(CPICE = 2.106E+3)
      PARAMETER(CHICE = 2.100E+3)
      PARAMETER(ESDMIN = 1.E-6)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(LSUBC = 2.501E+6)
      PARAMETER(LSUBS = 2.83E+6)
      PARAMETER(SIGMA = 5.67E-8)
      PARAMETER(TFREEZ = 273.15)
      PARAMETER(IRRSAT = 0.04)
      PARAMETER(MAXSWE = 0.10)
      PARAMETER(EPS = 3E-8)

c      DATA SNOEXP /1.0/
      DATA SNOEXP /2.0/

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C CONVERT POTENTIAL EVAP (ETP) FROM KG M-2 S-1 TO M S-1 AND THEN TO AN
C AMOUNT (M) GIVEN TIMESTEP (DT) AND CALL IT AN EFFECTIVE SNOWPACK
C REDUCTION AMOUNT, ETP2 (M).  THIS IS THE AMOUNT THE SNOWPACK WOULD BE
C REDUCED DUE TO EVAPORATION FROM THE SNOW SFC DURING THE TIMESTEP.
C EVAPORATION WILL PROCEED AT THE POTENTIAL RATE UNLESS THE SNOW DEPTH
C IS LESS THAN THE EXPECTED SNOWPACK REDUCTION.
C IF SEAICE (ICE=1), BETA REMAINS=1.
C ----------------------------------------------------------------------
C Compute snow melt and pack temperature via an energy balance approach 
C that allows for liquid water refreeze within the pack.  Use cold 
C content approach to transfer new snow to the pack and consider a 
C pseudo-two layer pack when SWE depth exceeds the prescribed threshold
C MAXSWE.  Parts of layering approach and iterative solution for pack 
C temperature were adapted from the VIC model through the Brent method 
C (1973: Algorithms for Minimization without Derivatives).
C Author: Ben Livneh -- Univ. of Wash.
C Last modified: 12/9/2009
C ----------------------------------------------------------------------
      swe0 = 0
      swe0 = esd
      if (prflag==2) write(*,*)'snwpac swe0 prcp prcp1',swe0,prcp,prcp1
C Any precip over non-snow area will be passed to SMFLX as an input (M/S)
cbl PRCP0 = total precip volume over time step
cbl PRCP1 = precip rate over non-snow covered area.
cbl This sum is passed to soil, so it can also include 
cbl dew/frost over non-snow cover as well as total melt
cbl PRCP2 = total precip volume over snow
      PRCP1 = (PRCP / 1000) * (1 - SNCOVR)
      EDIR = 0.0
      EDIR1 = 0.0
      EC = 0.0
      EC1 = 0.0
      DO K = 1,NSOIL
        ET(K) = 0.0
        ET1(K) = 0.0
      ENDDO
      ETT = 0.0
      ETT1 = 0.0
      ETNS = 0.0
      ETNS1 = 0.0
      ESNOW = 0.0
      ESNOW1 = 0.0
      ESNOW2 = 0.0
C ----------------------------------------------------------------------
C IF ETP<0 (DOWNWARD) THEN DEWFALL (=FROSTFALL IN THIS CASE).
C When ETP<0 all ET components are zero, and EDMND will be constrained
C to zero in SMFLX if it is negative. Dew (M/S) added below to prcp.
C ----------------------------------------------------------------------
      DEW = 0.0
      DEWNS = 0.0
      ETP1 = ETP*0.001
      IF (ETP .LT. 0.0) THEN
        DEWNS = -ETP1 * (1 - SNCOVR)
        ESNOW2 = ETP1 * DT
        ETANRG = ETP*((1.-SNCOVR)*LSUBC + SNCOVR*LSUBS)
      ELSE
C ----------------------------------------------------------------------
C When ETP>0, ET components are computed/scaled by SNCOVR (exception: 
C EDIR will be 0 for SAC case, to be recomputed by the ETP passed to 
C SMFLX, which will be scaled by SNCOVR therein)
C ----------------------------------------------------------------------
        IF (ICE .NE. 1) THEN
          IF (SNCOVR .LT. 1.) THEN
            CALL EVAPO (ETNS1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
     &            SH2O,STC,SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &            SMCREF,SHDFAC,CMCMAX,
     &            SMCDRY,CFACTR,
     &            EDIR1,EC1,ET1,ETT1,SFCTMP,Q2,NROOT,RTDIS,
     &            FXEXP,FX,MODEL_TYPE,SMCFLAG)
C ----------------------------------------------------------------------
            EDIR1 = EDIR1*(1.-SNCOVR)
            EC1 = EC1*(1.-SNCOVR)
            DO K = 1,NSOIL
              ET1(K) = ET1(K)*(1.-SNCOVR)
            END DO
            ETT1 = ETT1*(1.-SNCOVR)
            ETNS1 = ETNS1*(1.-SNCOVR)
C ----------------------------------------------------------------------
CBL Move this conversion from M/s to KG M-2 S-1 back to SFLX to avoid
CBL confusion and mistakes, as these values were passed to SMFLX 
CBL Therefore, had to supplement ETANRG calculation below with *1000's     
CBL            EDIR = EDIR1 * 1000.0
CBL            EC = EC1 * 1000.0
CBL            DO K = 1,NSOIL
CBL              ET(K) = ET1(K) * 1000.0
CBL            END DO
CBL            ETT = ETT1 * 1000.0
CBL            ETNS = ETNS1 * 1000.0
C ----------------------------------------------------------------------
          ENDIF
          ESNOW = ETP*SNCOVR
          ESNOW1 = ESNOW*0.001
          ESNOW2 = ESNOW1*DT
CBL          ETANRG = ESNOW*LSUBS + ETNS*LSUBC
          ETANRG = ESNOW*LSUBS + ETNS*1000*LSUBC
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C IF PRECIP IS FALLING, CALCULATE HEAT FLUX FROM SNOW SFC TO NEWLY
C ACCUMULATING PRECIP.  NOTE THAT THIS REFLECTS THE FLUX APPROPRIATE FOR
C THE NOT-YET-UPDATED SKIN TEMPERATURE (T1).  ASSUMES TEMPERATURE OF THE
C SNOWFALL STRIKING THE GOUND IS =SFCTMP (LOWEST MODEL LEVEL AIR TEMP).
C ----------------------------------------------------------------------
      FLX1 = 0.0
      IF (SNOWNG) THEN
         FLX1 = CPICE * PRCP * (T1 - SFCTMP)
      ELSE
         IF (PRCP .GT. 0.0) FLX1 = CPH2O * PRCP * (T1 - SFCTMP)
      ENDIF


C ----------------------------------------------------------------------
C Compute effective snow-ground temperature, T12, for the case when only
C partial snowcover and use this temperature to find a bulk surface temp
C ----------------------------------------------------------------------

      DSOIL = -(0.5 * ZSOIL(1))
      DTOT = SNOWH + DSOIL
      DENOM = 1.0 + DF1 / (DTOT * RR * RCH)
      T12A = ( (FDOWN-FLX1-FLX2-
     &       (EMISSNOW*SNCOVR+EMISS*(1.0-SNCOVR))*SIGMA*T24)/RCH
     &       + TH2 - SFCTMP - ETANRG/RCH ) / RR

      T12B = DF1 * STC(1) / (DTOT * RR * RCH)
      T12 = (SFCTMP + T12A + T12B) / DENOM      
      BETA = 1.0
	
C ----------------------------------------------------------------------
C IF POTENTIAL EVAP (SUBLIMATION) GREATER THAN DEPTH OF SNOWPACK.
C BETA<1
C SNOWPACK HAS SUBLIMATED AWAY, SET DEPTH TO ZERO.
C ----------------------------------------------------------------------

C      IF ((ESD-ESNOW2) .LE. ESDMIN) THEN
CCBL change to reflect an ablated pack
C         SNOMLT = MAX(0.0,ESD)
C         ESD = 0.0
C         EX = 0.0
C         WRITE(*,*)'MIN SNOW FLAG',esd,esnow2,snomlt
C
C      ELSE
C -----------------------------------------------------------------------
C Rebuild the snowpack by first separating new fallen precip and liquid 
C water, then add their cold contents sequentially, also account for 
C dew/frost fall in the process
C -----------------------------------------------------------------------
         PACSWE = 0.0
         TOPSWE = 0.0
         PACCC = 0.0
         TOPCC = 0.0
         SNOWCC = 0.0
         SNOW = 0.0
         RAIN = 0.0
         RFZNRG = 0.0
         RFZH20 = 0.0
         SNOMLT = 0.0
         FLX3 = 0.0
         PRCP0 = 0.0
         PRCP2 = 0.0

C Must convert PRCP0 to same units as snow (ESD =  M) from a rate to an amount
         PRCP0 = ((PRCP * DT) / 1000)
         PRCP2 = ((PRCP * DT) / 1000) * SNCOVR
         if (prflag==2) then
         write(*,*)'snwpac prcp0 sncovr',prcp0,sncovr
         write(*,*)'topswe esd pach20',topswe,esd,pach20
         endif
CBL       PRCP = PRCP + DEW
         TOPSWE = ESD - PACH20
         if (prflag==2) write(*,*)'topswe esd pach20',topswe,esd,pach20
C Remove newly added snow and redistribute it *correctly*
         IF (SNOWNG.OR.FRZGRA) THEN
            SNOW = PRCP0
            TOPSWE = TOPSWE - SNOW
            SNOWCC = CHICE*SNOW*(MIN(SFCTMP,TFREEZ)-TFREEZ)
            if (prflag==2) write(*,*)'topswe-snow',topswe,snow
            swe1 = esd - prcp0
            if (prflag==2) write(*,*)'swe1 = esd - prcp0',swe1,esd,prcp0
         ELSE
C This would be a rain on snow case where (SWE>0;SFCTMP>TFREEZ;T1>TFREEZ)
C Therefore, rain has not yet been added to the pack in SFLX
            if (prflag==2) then
            write(*,*)'adding rain...'
            write(*,*)'sfctmp t1',sfctmp,t1
            endif
            SNOWCC = 0.0
            SNOW = 0.0
            RAIN = PRCP0
            swe1 = esd
            if (prflag==2) write(*,*)'swe1 = esd',swe1,esd,prcp0
            PACH20 = PACH20 + PRCP0
         END IF
         if (prflag==2) write(*,*)'topswe esd pach20',topswe,esd,pach20
         swe2 = topswe + pach20 - prcp0
         if (prflag==2) then
            write(*,*)'swe2 topswe pach20',swe2,topswe,pach20,prcp0
         endif
C     Test to see if 2 layers needed
         IF (TOPSWE.GT.MAXSWE) THEN
            PACSWE = TOPSWE - MAXSWE
            TOPSWE = MAXSWE
         END IF
         
C     Compute cold contents based on current temperaturse
         TOPCC = CHICE*TOPSWE*(T1-TFREEZ)
         PACCC = CHICE*PACSWE*(TPACK-TFREEZ)
       
C     Now add snowfall and test if surface layer depth is exceeded by new snow
         IF (SNOW.GT.(MAXSWE-TOPSWE)) THEN
            DELSWE = (SNOW + TOPSWE) - MAXSWE
            IF (DELSWE.GT.TOPSWE) THEN
               DELCC = TOPCC + ((SNOW - MAXSWE)/SNOW)*SNOWCC
            ELSE
               DELCC = (DELSWE/TOPSWE)*TOPCC
            END IF
            TOPSWE = MAXSWE
            TOPCC = TOPCC + SNOWCC - DELCC
            PACSWE = PACSWE + DELSWE
            PACCC = PACCC + DELCC
         ELSE
            TOPSWE = TOPSWE + SNOW
            TOPCC = TOPCC + SNOWCC
         END IF
         ESD = TOPSWE + PACSWE
         swe2 = esd + pach20
         if (prflag==2) then
         write(*,*)'added snowfall swe2',swe2,esd,topswe,pach20,snow      
         endif
C     Compute initial layer temperatures, T1_0, based on cold contents
         IF (TOPSWE.GT.0.0) THEN
            T1_0 = TOPCC/(CHICE*TOPSWE) + TFREEZ
         ELSE
            T1_0 = TFREEZ
         END IF
         IF (PACSWE.GT.0.0) THEN
            TPACK = PACCC/(CHICE*PACSWE) + TFREEZ
         ELSE
         END IF

C Determine liquid water within the pack layer(s). PACH20 now changes 
C meaning to be the liquid water content of the pack layer
         TOPH20 = ((PACH20 - RAIN)/ESD)*TOPSWE + RAIN
         PACH20 = ((PACH20 - RAIN)/ESD)*PACSWE 
         swe2 = esd + pach20 + toph20
         if (prflag==2) then
         write(*,*)'enumerated h20 swe2',swe2,esd,pach20,toph20    
         endif
C mek Feb2004
C non-linear weighting of snow vs non-snow covered portions of gridbox
C so with SNOEXP = 2.0 (>1); take into account partial snow cover.
C After energy balance, proceed with 'snow only' calculations for T1
C after which reinvoke m. ek nonlinear weighting of T1 to be returned
C as the 'bulk' surface temperature

         IF (T12.GT.TFREEZ) THEN
            T1 = TFREEZ * SNCOVR**SNOEXP + T12 * (1.0 - SNCOVR**SNOEXP)
         ELSE
            T1 = TFREEZ
         END IF

C Now solve the snowpackenergy balance (SEB) with an initial guess of a
C ripe snowpack (snow surface = TFREEZ) to obtain the net energy available 
C to melt/refreeze the pack, QNET.
  
         CALL SEB (T1,QNET,TOPSWE,RCH,CH,DT,Q2,SFCPRS,SFCTMP,SFCSPD,
     &        STC,DTOT,DF1,T1_0,TH2,FDOWN,FLX1,FLX2,SSOIL,SEH,DELCC,
     &        ETANRG,ESNOW2,TOPH20,PACH20,PACCC,RFZNRG,EMISSNOW,EMISS,
     &        SNCOVR,NSOIL)
         

C Test the net energy available to the pack, QNET (which accounts for 
C refreeze energy, RFZNRG of pack water) is balanced, or in surplus or 
C deficit. First case is an energy surplus, within which refreeze will 
C occur if the residual of QNET and RFZNRG is negative, otherwise, the
C second case is that melt will occur.  The third case is an energy
C deficit, in which a colder temperature is solved for iteratively and
C all water from the surface layer is refrozen.  At the end, test if 
C the cold content of pack layer is suffice to refreez pack water.

         IF (QNET.GE.0.0) THEN
            T1 = TFREEZ
C First case is that the energy needed to refreeze available water 
C exceeds a melt energy defecit
            IF ((QNET - RFZNRG).LT.0.0) THEN
               RFZH20 = (ABS((QNET-RFZNRG))*DT)/(LSUBF*1000)
               TOPSWE = TOPSWE + RFZH20
               TOPH20 = TOPH20 - RFZH20
               FLX3 = 0.0
               SNOMLT = 0.0
               swe2=topswe+toph20+pacswe+pach20
               if (prflag==2) then
               write(*,*)'case 1a swe2 snomlt',swe2,snomlt
               write(*,*)'topswe toph20 pacswe pach20'
               write(*,*)topswe,toph20,pacswe,pach20
               endif
C Second case is that there is a melt energy surplus, such that no 
C refreeze takes place and some melt will occur, use melt energy, 
C MLTNRG, to track its dissipation
            ELSE
               RFZH20 = 0.0
               FLX3 = QNET - RFZNRG
               MLTNRG = FLX3
               SNOMLT = ((MLTNRG)*DT)/(LSUBF*1000)
C Test if snow melt is contained within surface layer SWE, if not, then
C ripen a pack layer (if one exists) through depletion of cold content 
C and then melt accordingly
               IF (SNOMLT.LT.TOPSWE) THEN
                  TOPSWE = TOPSWE - SNOMLT
                  TOPH20 = TOPH20 + SNOMLT
                  MLTNRG = 0.0
CBL Since snomelt is transfered here as pack water, set snomlt to zero
CBL pack water is handled down below when checking maximum pack water limits
                  swe2=topswe+toph20+pacswe+pach20
                  if (prflag==2) then
                  write(*,*)'case 1bi swe2 snomlt',swe2,snomlt
                  write(*,*)'topswe toph20 pacswe pach20'
                  write(*,*)topswe,toph20,pacswe,pach20
                  endif
                  SNOMLT = 0.0
C If a pack layer exists, use melt energy to heat, then melt the pack 
C layer very rare case...
               ELSE
CBL All top-layer swe becomes liquid, which accounts for only part of the 
CBL melt energy conversion needed.
                  TOPH20 = TOPH20 + TOPSWE
                  IF   (MLTNRG.LT.((TOPSWE*LSUBF*1000/DT)+PACCC+
     &                 (PACSWE*LSUBF*1000/DT))) THEN
                     MLTNRG = MLTNRG - (TOPSWE*LSUBF*1000)/DT
                     IF (MLTNRG.LT.PACCC) THEN     
                        TPACK = TPACK+(MLTNRG/ABS(PACCC))*(TFREEZ-TPACK)
                        PACCC = CHICE*PACSWE*(TPACK-TFREEZ)
                     ELSE
                        MLTNRG = MLTNRG - PACCC
                        PACCC = 0.0
                        TPACK = TFREEZ
                        PACSWE = PACSWE - (MLTNRG*DT)/(LSUBF*1000)
                        PACH20 = PACH20 + (MLTNRG*DT)/(LSUBF*1000)
                        MLTNRG = 0.0
C If melt energy doesn't exceed CC of pack layer, warm the pack 
C accordingly
                     END IF
                     TOPSWE = 0.0
                     TOPCC = 0.0
                     swe2=topswe+toph20+pacswe+pach20
                     if (prflag==2) then
                     write(*,*)'case 1bii swe2 snomlt',swe2,snomlt
                     write(*,*)'topswe toph20 pacswe pach20'
                     write(*,*)topswe,toph20,pacswe,pach20
                     endif
CBL Set snomlt to zero, since the above conversion of frozen SWE into liquid
CBL water has accounted for the melt energy conversion.  The additional 
CBL liquid water will be handled down below
                     SNOMLT = 0.0
                  ELSE
C Complete melt situation; runoff pack liquid water
                     PACH20 = PACH20 + PACSWE
                     PACSWE = 0.0
                     PACCC = 0.0
CBL                     SNOMLT = TOPH20 + PACH20
CBL Since the above PACSWE was added to the PACH20, then after the energetics
CBL all the excess water (TOPH20+PACH20>MAXH20) will become snomlt, so set
CBL SNOMLT to zero here, so as not to double count it.
                     SNOMLT = 0.0
                     FLX3 = (SNOMLT*LSUBF*1000)/DT
                     T1 = TFREEZ
                     TPACK = T1
                     swe2=topswe+toph20+pacswe+pach20
                     if (prflag==2) then
                     write(*,*)'case 1biii swe2 snomlt',swe2,snomlt
                     write(*,*)'topswe toph20 pacswe pach20'
                     write(*,*)topswe,toph20,pacswe,pach20
                     endif
                  END IF
                  TOPSWE = 0.0
                  TOPCC = 0.0
c                  TOPH20 = 0.0
                  swe2=topswe+toph20+pacswe+pach20
                  if (prflag==2) then
                  write(*,*)'case 1b swe2 snomlt',swe2,snomlt
                  write(*,*)'topswe toph20 pacswe pach20'
                  write(*,*)topswe,toph20,pacswe,pach20
                  endif
               END IF
            END IF
            swe2 = topswe + pacswe +toph20 + pach20
            if (prflag==2) then
            write(*,*)'warming case swe2 snomlt',swe2,snomlt     
            write(*,*)'topswe toph20 pacswe pach20'
            write(*,*)topswe,toph20,pacswe,pach20  
            endif
         ELSE
C Third case, residual of QNET and RFZNRG is negative, meaning T1<TFREEZ 
C and will be solved iteratively, using the Brent (1973) method. Adapted
C from "Numerical Recipes in Fortran 77" - sec. 9.3. Initial guesses are
C bracketed by a variable BOUNDS and convergence tolerance and maximum 
C number of iterations (TOL AND ITRMAX) can be specified

            ITRMAX = 40
            ITER = 0
            TOL = 0.01
C     BOUNDS = 100.0
            BOUNDS = 50.0
C     T1LOW = T1 - BOUNDS
            T1LOW = T1_0 - BOUNDS
            T1HIGH = T1
            IF (T1LOW > T1HIGH - BOUNDS) THEN
               T1LOW = T1HIGH - BOUNDS
            ENDIF
            A = T1LOW
            B = T1HIGH
            CALL SEB (A,FA,TOPSWE,RCH,CH,DT,Q2,SFCPRS,SFCTMP,SFCSPD,
     &           STC,DTOT,DF1,T1_0,TH2,FDOWN,FLX1,FLX2,SSOIL,SEH,DELCC,
     &           ETANRG,ESNOW2,TOPH20,PACH20,PACCC,RFZNRG,EMISSNOW,
     &           EMISS,SNCOVR,NSOIL)  
            
            CALL SEB (B,FB,TOPSWE,RCH,CH,DT,Q2,SFCPRS,SFCTMP,SFCSPD,
     &           STC,DTOT,DF1,T1_0,TH2,FDOWN,FLX1,FLX2,SSOIL,SEH,DELCC,
     &           ETANRG,ESNOW2,TOPH20,PACH20,PACCC,RFZNRG,EMISSNOW,
     &           EMISS,SNCOVR,NSOIL)

            IF(FA*FB.GT.0.0) THEN
               WRITE(*,*)'SNWPAC: SNOW PACK T',A,B, 
     &              'FALLS OUTSIDE LOW AND HIGH BOUNDS',FA,FB
               ERFLAG = 1
               GOTO 20
            END IF
            
            C = B
            FC = FB
 10         IF (ITER.LT.ITRMAX) THEN
               ITER = ITER + 1
               IF (FB*FC.GT.0.0) THEN
C     Rename a,b,c and bounding interval, d
                  C = A
                  FC = FA
                  D = B - A
                  E = D
               END IF
               IF (ABS(FC).LT.ABS(FB)) THEN
                  A = B
                  B = C
                  C = A
                  FA = FB
                  FB = FC
                  FC = FA
               END IF
               TOL1 = 2*EPS*ABS(B) + 0.5*TOL
C     Convergence check
               XM = 0.5*(C-B)
               IF ((ABS(XM).LE.TOL1).OR.(FB.EQ.0.)) THEN
C     WRITE(*,*)'Converged T1 B',T1,B
                  T1 = B
                  GOTO 20
               END IF
               IF ((ABS(E).GE.TOL1).AND.(ABS(FA).GT.ABS(FB))) THEN
C     Attempt inverse quadratic interpolation
                  S = FB/FA
                  IF (A.EQ.C) THEN
                     P = 2.*XM*S
                     Q = 1. - S
                  ELSE
                     Q = FA/FC
                     R = FB/FC
                     P = S*(2.*XM*Q*(Q-R) - (B-A)*(R-1.))
                     Q = (Q-1.)*(R-1.)*(S-1.)
                  END IF
C     Check whether in bounds
                  IF (P.GT.0.) THEN
                     Q = Q * -1.
                  END IF
                  P = ABS(P)
                  IF ((2.*P).LT.MIN((3.*XM*Q-ABS(TOL1*Q)),ABS(E*Q))) THEN
C     Accept intermpolation
C     WRITE(*,*)(2*P),(3.*XM*Q - ABS(TOL1*Q)),ABS(E*Q)
                     E = D
                     D = P/Q
                  ELSE
C     Interpolation failed, use bisection
                     D = XM
                     E = D
                  END IF
               ELSE
C     Bounds decreasing too slowly, use bisectoin
                  D = XM
                  E = D
               END IF
C     Move last best guess to A
               A = B
               FA = FB
C     Evaluate new trial root
               IF (ABS(D).GT.TOL1) THEN
                  B = B + D
               ELSE
                  B = B + SIGN(TOL1,XM)
               END IF
               CALL SEB (B,FB,TOPSWE,RCH,CH,DT,Q2,SFCPRS,SFCTMP,SFCSPD,
     &              STC,DTOT,DF1,T1_0,TH2,FDOWN,FLX1,FLX2,SSOIL,SEH,
     &              DELCC,ETANRG,ESNOW2,TOPH20,PACH20,PACCC,RFZNRG,
     &              EMISSNOW,EMISS,SNCOVR,NSOIL)     
               
C     WRITE(*,*)'**BOTTOM OF BRENT** B FB',B,FB
               GOTO 10
            END IF
 20         RFZH20 = TOPH20
            TOPH20 = 0.0
            T1 = B
            TOPSWE = TOPSWE + RFZH20
            SNOMLT = 0.
            FLX3 = 0.
            swe2 = topswe+pacswe
            if (prflag==2) write(*,*)'bottom of brent swe2',swe2,snomlt              
         END IF
         
C     Refreeze any pack layer water with pack cold content (J) if available
         RFZNRG = PACH20*LSUBF*1000
         
C     First case that refreezable pack water is exceeded by cold content
         IF (PACCC.LT.(RFZNRG*-1)) THEN
            PACSWE = PACSWE + PACH20
            PACH20 = 0.
            IF (PACSWE.GT.0.) THEN
               PACCC = PACSWE*CHICE*(TPACK-TFREEZ) + RFZNRG
               TPACK = PACCC/(CHICE*PACSWE) + TFREEZ
            ELSE
               TPACK = TFREEZ
            END IF
         ELSE
C     Second case: pack cold content doesn't exceed refreezable pack water
            TPACK = TFREEZ
            PACH20 = PACH20 - (ABS(PACCC))/(LSUBF*1000)
            PACSWE = PACSWE + (ABS(PACCC))/(LSUBF*1000)     
            PACCC = 0.0
         END IF
         
C Use the irreducible saturation water, MAXH20, to update the pack water.
C An amount of water equal to 4% of the pore volume of the pack, IRRSAT
C was chosen and converted to a fraction of total swe, FRCH20, based on 
C experimental data, such as Denoth, 2003 (Cld Rgn Sci&Tech 37,227-32, 
C showing MAXH20 asymptotes to 4% of snow pore volume,compute pore space 
C using density from call to SNOWPACK. Send the residual water to snomlt
         ESD = TOPSWE + PACSWE   
         if (prflag==2) write(*,*)'esd=topswe+pacswe',esd,topswe,pacswe
C     MAXH20 = ESD * FRCH20
C     Test Denoth estimate of irreducible saturation = 4% of pore volume
cbl Special thanks to Rick Steed
         ZZ1 = 1.0
         YY = STC(1)-0.5*SSOIL*ZSOIL(1)*ZZ1/DF1
         CALL SNOWPACK (ESD,DT,SNOWH,SNDENS,T1,YY)
         FRCH20 = IRRSAT*((1 - SNDENS/0.9)/SNDENS)
C     FRCH20 = 0.04
         MAXH20 = ESD * FRCH20
         IF ((PACH20+TOPH20).GE.MAXH20) THEN
            SNOMLT = SNOMLT + (PACH20+TOPH20-MAXH20)
            if (prflag==2) then
            write(*,*)'maxh20 exceeded snomlt maxh20',snomlt,maxh20   
            write(*,*)'maxh20 pach20 toph20',maxh20,pach20,toph20
            endif
            PACH20 = MAXH20
            TOPH20 = 0.0
            ESD = ESD + PACH20
            if (prflag==2) write(*,*)'esd+=pach20',esd,pach20
         ELSE
            PACH20 = PACH20 + TOPH20
            ESD = ESD + PACH20
         END IF
         
         swe2 = esd
         if (prflag==2) then
            write(*,*)'enumerate water and density swe2',swe2,snomlt   
         endif
CBL ESNOW will be the evap over the snow potion, which may have a
CBL  different sign from ETP (rarely), but with MUST be updated in SFLX.
         DEWS = 0.0
C     BL convert ESNOW2(VMF) from M to M/S (ESNOW) for transfer       
         ESNOW2 = (ESNOW2 * SNCOVR) 
         ESNOW = ESNOW2 / DT
         IF (ESNOW .LT. 0.) THEN
C     BL Frost-fall over the snowpack
            DEWS = -ESNOW
            ESD = ESD + (ESNOW2 * -1.)
            if (prflag==2) then
               write(*,*)'enumerated frostfall esd snomlt',esd,snomlt
            endif
         ELSE
C     BL Vapor flux away from snowpack, compute ESNOW
            IF ((ESD-ESNOW2) .LE. ESDMIN) THEN
C     BL change to reflect an ablated pack
               ESNOW2 = MAX(0.0,ESD)
               ESNOW = ESNOW2 / DT
               SNOMLT = SNOMLT + MAX(0.0,ESD)
               ESD = 0.0
               EX = 0.0
               if (prflag==2) then
               write(*,*)'enumerated submimation esd snomlt',esd,snomlt
               endif
            ELSE
               ESD = ESD - ESNOW2
               if (prflag==2) then
               write(*,*)'enumerated submimation esd snomlt',esd,snomlt
            endif
            ENDIF
         ENDIF
         
C     BL Check for T1 > TFREEZ to ensure correct latent heat flag applied
C     BL in SFLX
         IF (T1 .GT. TFREEZ) LVRAIN = .TRUE.
         
C     Redefine T1 to be its intended, bulk skin temperature for the partial
C     snow cover case, if it was invoked before the initial call to SEB.
         
         IF (T12.GT.TFREEZ) THEN
            T1 = T1 * SNCOVR**SNOEXP + T12 * (1.0 - SNCOVR**SNOEXP)       
         END IF
         
C----------------------------------------------------------------------
C     END OF 'ESD .LE. ETP2' IF-BLOCK
C----------------------------------------------------------------------
C      ENDIF
      
      EX = SNOMLT/DT
      PRCP1 = PRCP1 + (EX * SNCOVR)
      
      swe2 = esd 
      if (prflag==2) write(*,*)'enumerated vapor swe2',swe2       
c     swedif = swe - swe
      swedif = swe2 - swe1
      bal = prcp0 - swedif - esnow2 - snomlt
      if (prflag==2) then
      write(*,*)'prcp0 swedif',prcp0,swedif
      write(*,*)'swe1 swe2',swe1,swe2
      write(*,*)'esnow2 snomlt',esnow2,snomlt
      endif
      bal=prcp0-esnow2-snomlt
      if (prflag==2) write(*,*)'change in i/o storage',bal,swedif
      bal = prcp0 - swedif - esnow2 - snomlt
      if(prflag==2) then
      write(*,*)'bal-swe',bal,prcp0,swedif,esnow2,snomlt,swe1,swe2
      endif
      if (prflag==2) write(*,*)'snwpac prcp1 snomlt',prcp1,ex

C Add any dewfall (M/S) and melt (M/S) over non snow-covered area to PRCP1
      PRCP1 = PRCP1 + DEWNS 

C DEW will is the combination of snowfree area (etp) and snow area (esnow)
      DEW = DEWNS + DEWS
C ----------------------------------------------------------------------
C FINAL BETA NOW IN HAND, SO COMPUTE EVAPORATION.  EVAP EQUALS ETP
C UNLESS BETA<1.
C ----------------------------------------------------------------------
c      ETA = BETA*ETP
      if (prflag==2)write(*,*)'snwpac prcp1 dewns snomlt',prcp1,dewns,ex
C ----------------------------------------------------------------------
C SET THE EFFECTIVE POTNL EVAPOTRANSP (ETP1) TO ZERO SINCE THIS IS SNOW
C CASE, SO SURFACE EVAP NOT CALCULATED FROM EDIR, EC, OR ETT IN SMFLX
C (BELOW).
C IF SEAICE (ICE=1) SKIP CALL TO SMFLX.
C SMFLX RETURNS UPDATED SOIL MOISTURE VALUES.  IN THIS, THE SNOW PACK
C CASE, ETA1 IS NOT USED IN CALCULATION OF EVAP.
C ----------------------------------------------------------------------
c      ETP1 = 0.0
      IF (ICE .NE. 1) THEN
         IF (MODEL_TYPE.EQ.0) THEN
            CALL SMFLX (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &           SH2O,SLOPE,KDT,FRZFACT,
     &           SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &           SHDFAC,CMCMAX,
     &           RUNOFF1,RUNOFF2,RUNOFF3,
     &           EDIR1,EC1,ET1,
     &           DRIP)
         ENDIF
         IF (MODEL_TYPE.EQ.1) THEN
            CALL SMFLX_SAC (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &           SH2O,SLOPE,KDT,FRZFACT,SNCOVR,SMCREF,
     &           SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &           SHDFAC,CMCMAX,SFCTMP,ETP1,ETA1,
     &           RUNOFF1,RUNOFF2,RUNOFF3,EDIR1,EC1,ET1,DRIP, 
     &           FRZST,FRZPAR,SACST,SACPAR,PRFLAG,CELLID,MSTEP,wcrit)
         ENDIF
      ENDIF
      
    
      if (prflag==2) write(*,*)'esnow2 etp1 dew',esnow2,etp1,dew
      ettemp = 0
      do i = 1,nsoil
         ettemp = ettemp + et1(i)
      enddo
      if (prflag==2) then
      write(*,*)'snwpac edir ec ett esnow esnow2 eta'
      write(*,*)edir1,ec,ettemp,esnow,esnow2,eta1
      endif

C ----------------------------------------------------------------------
C BEFORE CALL SHFLX IN THIS SNOWPACK CASE, SET ZZ1 AND YY ARGUMENTS TO
C SPECIAL VALUES THAT ENSURE THAT GROUND HEAT FLUX CALCULATED IN SHFLX
C MATCHES THAT ALREADY COMPUTER FOR BELOW THE SNOWPACK, THUS THE SFC
C HEAT FLUX TO BE COMPUTED IN SHFLX WILL EFFECTIVELY BE THE FLUX AT THE
C SNOW TOP SURFACE.  T11 IS A DUMMY ARGUEMENT SO WE WILL NOT USE THE
C SKIN TEMP VALUE AS REVISED BY SHFLX.
C ----------------------------------------------------------------------
      ZZ1 = 1.0
      YY = STC(1)-0.5*SSOIL*ZSOIL(1)*ZZ1/DF1
      T11 = T12

C ----------------------------------------------------------------------
C SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX 
C (SSOIL1) AND THE SKIN TEMP (T11) OUTPUT FROM THIS SHFLX CALL ARE NOT
C USED  IN ANY SUBSEQUENT CALCULATIONS. RATHER, THEY ARE DUMMY VARIABLES
C HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT FLUX ARE
C UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
C ----------------------------------------------------------------------
      CALL SHFLX (SSOIL1,STC,SMC,SMCMAX,NSOIL,T11,DT,YY,ZZ1,ZSOIL,
     &            TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &            QUARTZ,CSOIL)

      
C ----------------------------------------------------------------------
C SNOW DEPTH AND DENSITY ADJUSTMENT BASED ON SNOW COMPACTION.  YY IS
C ASSUMED TO BE THE SOIL TEMPERTURE AT THE TOP OF THE SOIL COLUMN.
C ----------------------------------------------------------------------
      IF (ESD .GT. 0.) THEN
        CALL SNOWPACK (ESD,DT,SNOWH,SNDENS,T1,YY)
C        WRITE(*,*)'SNDENS',SNDENS
      ELSE
        ESD = 0.
        SNOWH = 0.
        SNDENS = 0.
        SNCOND = 1.
        SNCOVR = 0.
      ENDIF

C ----------------------------------------------------------------------
C END SUBROUTINE SNWPAC
C ----------------------------------------------------------------------
      RETURN
      END

       SUBROUTINE SEB (TSURF,QNET,TOPSWE,RCH,CH,DT,Q2,SFCPRS,SFCTMP,
     &	SFCSPD,STC,DTOT,DF1,T1_0,TH2,FDOWN,FLX1,FLX2,SSOIL,SEH,DELCC,
     &  ETANRG,ESNOW2,TOPH20,PACH20,PACCC,RFZNRG,EMISSNOW,EMISS,SNCOVR,
     &  NSOIL)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SEB
C ----------------------------------------------------------------------
C CALCULATE THE SNOWPACK ENERGY BALANCE AND DEFINE THE AMOUNT OF MELT ENERGY
C QNET AVAILABLE BASED PRIMARILY ON THE SURFACE TEMPERATURE OF THE PACK
C ----------------------------------------------------------------------
      INTEGER NSOIL

      REAL TSURF
      REAL QNET
      REAL TOPSWE
      REAL RCH
      REAL CH
      REAL DT
      REAL Q2
      REAL SFCPRS
      REAL SFCTMP
      REAL SFCSPD
      REAL STC(NSOIL)
      REAL DTOT
      REAL DF1
      REAL T1
      REAL TH2
      REAL FDOWN
      REAL FLX1
      REAL FLX2
      REAL SSOIL
      REAL SEH
      REAL DELCC
      REAL ETANRG
      REAL ESNOW2
      REAL TFREEZ
      REAL SIGMA
      REAL LSUBC
      REAL LSUBS
      REAL LSUBF
      REAL CHICE
      REAL TCELS
      REAL COEFF1
      REAL ESSNOW
      REAL EACT
      REAL AIRDEN
      REAL VMF
      REAL LHF
      REAL SIGT4
      REAL TOPH20
      REAL PACH20
      REAL PACCC 
      REAL RFZNRG
      REAL EMISSNOW
      REAL EMISS
      REAL SNCOVR
      REAL T1_0

      PARAMETER(TFREEZ = 273.15)
      PARAMETER(SIGMA = 5.67E-8)
      PARAMETER(LSUBC = 2.501E+6)
      PARAMETER(LSUBS = 2.83E+6)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(CHICE = 2.100E+3)
	  
      TCELS = TSURF - TFREEZ
	  
C Calculate the saturated vapor pressure of the pack, ESSNOW
      COEFF1 = 1 + 0.00972*TCELS + 0.000042*TCELS**2
      ESSNOW = 0.6107*(2.7182818**((17.269*TCELS)/(237.3 + TCELS)))*1000
      IF (TSURF.LT.TFREEZ) THEN
        ESSNOW = ESSNOW*COEFF1
       END IF

C Variable for actual vapor pressure
       EACT = (Q2/(Q2+0.622))*SFCPRS
       AIRDEN = SFCPRS/(287*(SFCTMP+0.622*Q2))
         
C Compute vapor mass flux, VMF (M), and divide through by the density of water
       VMF = (AIRDEN*(0.622/SFCPRS)*(EACT-ESSNOW)*CH)/(1000*SFCSPD)
C         WRITE(*,*)'vmf adns prs wnd',VMF,AIRDEN,SFCPRS,SFCSPD
C         WRITE(*,*)'eact essnow ch',EACT,ESSNOW,CH
       IF (TSURF.GE.TFREEZ) THEN
         LHF = LSUBC*VMF*1000
       ELSE
         LHF = LSUBS*VMF*1000
       END IF

       ESNOW2 = VMF
       ETANRG = LHF
	
C Sensible, ground, heat fluxes and long wave emitted
       SEH = RCH * (TSURF - TH2)
       SSOIL = DF1 * (TSURF - STC(1))/ DTOT
       SIGT4 =(EMISSNOW*SNCOVR+EMISS*(1.0-SNCOVR))*SIGMA*TSURF**4

C Change in cold content from the initial condition
       DELCC = (CHICE * TOPSWE * (TSURF - T1_0))/DT
	
C Compute refreeze energy, RFZNRG, needed for pack layer
       RFZNRG = (LSUBF*TOPH20*1000)/DT
C   RFZNRG = (LSUBF*TOPH20*1000)/DT + PACCC + (LSUBF*PACH20*1000)/DT	

C Compute net melt energy
       QNET = FDOWN-SIGT4-SSOIL-SEH-ETANRG-FLX1-FLX2-DELCC+RFZNRG

C ----------------------------------------------------------------------
C END SUBROUTINE SEB
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE TR_19 (ETP,ETA,PRCP,PRCP1,SNOWNG,FRZGRA,SMC,SMCMAX,
     &               SMCWLT,SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &               SBETA,DF1,PACH20,CH,SFCSPD,
     &               Q2,T1,TPACK,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &               SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,SNEQV,SNDENS,
     &               SNOWH,SH2O,SLOPE,KDT,FRZX,PSISAT,SNUP,
     &               ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &               RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &               ICE,RTDIS,QUARTZ,FXEXP,FX,CSOIL,UZTWM,UZFWM,UZK,
     &               PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,LZFSM,LZFPM,LZSK,
     &               LZPK,PFREE,SIDE,RSERV,UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,
     &               ADIMC,ERFLAG,BETA,DRIP,DEW,FLX1,FLX2,FLX3,ESNOW,
     &               EMISS,EMISSNOW,ZLVL,SOLDN,LWDN,SNOALB,SN_NEW,
C State and output variables
     &     NSNOW,DSNOW,PSNOW,RTTSNOW,RTTDTSNOW,WTSNOW,WTDTSNOW,TTSNOW,
     &     TTDTSNOW)


      IMPLICIT NONE

      INTEGER ICE
      INTEGER NROOT
      INTEGER NSOIL
      INTEGER LSTSNW
      INTEGER ITRMAX
      INTEGER ITER
      INTEGER K
      INTEGER J
      INTEGER ERFLAG

      LOGICAL SNOWNG
      LOGICAL FRZGRA

      REAL BEXP,BETA,CC,CFACTR,CMC,CMCMAX,CHICE,CP,CPH2O,CPICE,CSOIL,ESD
      REAL DENOM,DEW,DF1,DKSAT,DRIP,DSOIL,DTOT,DT,DWSAT,EC,EDIR,EPSCA
      REAL ESDMIN,EXPSNO,EXPSNOI,ETA,ETA1,ETP,ETP1,ETP2,ET(NSOIL),SNOMLT
      REAL ETT,EX,EXPFAC,FX,EMISS,EMISSNOW,FDOWN,FXEXP,FLX1,FLX2,FLX3
      REAL F1,KDT,LSUBF,LSUBC,LSUBS,PC,PRCP,PRCP1,Q2,RCH,RR,RTDIS(NSOIL)
      REAL SSOIL,SBETA,SSOIL1,SFCTMP,SHDFAC,SIGMA,SMC(NSOIL),SH2O(NSOIL)
      REAL SMCDRY,SMCMAX,SMCREF,SMCWLT,SNOWLT,SNOWH,STC(NSOIL),T12,T14
      REAL T1,T1A,T1B,TH2,T2,TBOT,ZBOT,YY,ZSOIL(NSOIL),ZZ1,TFREEZ,SALP
      REAL SFCPRS,SLOPE,FRZFACT,PSISAT,SNUP,RUNOFF1,RUNOFF2,RUNOFF3,T11
      REAL QUARTZ,SNDENS,SNCOND,RSNOW,SNCOVR,QSAT,ETP3,SEH,CSNOW,EC1
      REAL EDIR1,ET1(NSOIL),ETT1,ETNS,ETNS1,ESNOW,ESNOW1,ESNOW2,ETANRG
      REAL ESSNOW,SNEQV,FRZX,SOLDN,LWDN,SNOALB,SN_NEW
      REAL LHF
      REAL QNET
      REAL FRCH20
      REAL IRRSAT
      REAL EACT
      REAL VMF
      REAL AIRDEN
      REAL CH
      REAL SFCSPD
      REAL CORECT
      REAL TOPSWE
      REAL PACSWE
      REAL MAXSWE
      REAL TOPCC
      REAL PACCC
      REAL TOPH20
      REAL PACH20
      REAL MAXH20
      REAL RAIN
      REAL SNOW
      REAL SNOWCC
      REAL DELSWE
      REAL DELCC
      REAL RFZNRG
      REAL RFZH20
      REAL MLTNRG
      REAL A
      REAL B
      REAL C
      REAL D
      REAL E
      REAL P
      REAL Q
      REAL R
      REAL S
      REAL FA
      REAL FB
      REAL FC
      REAL T1HIGH
      REAL T1LOW
      REAL TOL
      REAL TOL1
      REAL XM
      REAL EPS
      REAL BOUNDS
      REAL T12A
      REAL T12B
      REAL T24
      REAL TPACK
      REAL T1_0
      REAL UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,LZFSM
      REAL LZFPM,LZSK,LZPK,PFREE,SIDE,RSERV,UZTWC,UZFWC,LZTWC,LZFSC
      REAL LZFPC,ADIMC

      REAL SNOEXP,ZLVL
      INTEGER MODEL_TYPE,ITOPT,ITS,IPUNCH,IGRAD,IQAE,IFU,NN,NEWPACK
      INTEGER NSNOW,I,N
      REAL DELTAT,DTOUT,THEDA,TOLER,GRMAX,DTONE,TEXP,THICK,CTHICK
      REAL ADJQA,HEIGHT,D_O,PA,X,DTG,TCG,DCG,SFC,RCF,PLWHC,PLWMAX
      REAL PLWDEN,FUCOEF,COEFKE,CV,RICRIT,ZO,C1,C2,C3,C4,DMETA,C5
      REAL CW1,CW2,CW3,CW4,G1,G2,G3,MM,DEN,TEMP,THICK0
      REAL DSNOW(NSNOW),PSNOW(NSNOW),RTTSNOW(NSNOW),RTTDTSNOW(NSNOW)
      REAL WTSNOW(NSNOW),WTDTSNOW(NSNOW),TTSNOW(NSNOW)
      REAL TTDTSNOW(NSNOW),SCOUTSNOW,VAPOURSNOW,QGSNOW,WATER0

      DATA SNOEXP /2.0/

      PARAMETER(TFREEZ = 273.15)
C Metamorphism and water transmission parameters
C Compaction parameters
      PARAMETER (C1 = 0.21,C2 = 21)
C Destructive metamorphism parameters
      PARAMETER (C3 = 0.01,C4 = 0.04, DMETA = 0.15)
C Liquid water metamorphism & transmission
      PARAMETER (C5 = 2.0,CW1 = 10.,CW2 = 1.0,CW3 = 5.0, CW4 = 450.)
C Grain size params: grain size = G1+G2*(DEN**2)+G3*(DEN**4)
      PARAMETER (G1 = 0.16,G2 = 0.0,G3 = 110.)
                       
C ----------------------------------------------------------------------
C SUBROUTINE TR_19
C ----------------------------------------------------------------------
C CALCULATE THE SNOWPACK ENERGY BALANCE, PACK TEMPERATURE, AND LIQUID
C WATER CONTENT ACCORDING TO ANDERSONS POINT ENEGY BALANCE (Ref:
C Anderson, E. A. 1976. A Point Energy and Mass Balance Model of a Snow 
C Cover. A Technical Report NWS 19, NOAA).
C A multi-layer snow model (<100 layers) that uses an implicit finite-
C difference scheme via the Newton-Raphson method to solve the snowpack 
C energy balance for each layer. I/O of original code was modified to 
C run at sub-hourly timesteps using the given forcing variables, rather
C than reading inputs from cards and tape.
C Misc. parameter values were a mix from Eric Andersons suggestions and 
C estimations from Danville, VT
C Ben Livneh 01/25/2009 Univ. of Washington
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Define input parameters for PEMB model.  Parameters relating to model
C start time and model ending time have been removed, since these are
C handled by driver code framework.  Other parameters related to punch
C cards and input snow-surface-temperature were removed from model code
C ----------------------------------------------------------------------

Change value for accumulation and melt season based on SNOW SURF TEMP
      ITOPT = 0
      ITS = 1
      IPUNCH = 0
      DELTAT = DT
      THEDA = 0.5
C Tolerance in deg C of iterative solution
      TOLER = 0.05 
      IGRAD = 0
      GRMAX = 2.0
      DTONE = 3.0
      TEXP = 0.5
C May wish to change IQAE to 1, to get model to estimate longwave
      IQAE = 0
C Layer thickness in centimeters
      THICK = 2.5
      CTHICK = 0.05
      ADJQA = 1.0
C Theoretical wind function with stability adjustment (IFU=1); empirical 
C wind (IFU=0) computed each time step; may wish to change to IFU = 1.
      IFU = 0
C Height of wind, temp, and vapor pressure measurements (meters)
      HEIGHT = ZLVL
C Diffusion coefficient for water vapor
      D_O = 0.9
      PA = SFCPRS / 100
      X = 14
C Depth of soil temperature measurement
      DTG = (-ZSOIL(NSOIL))*100
C Thermal conductivity of the soil in cal/cm/sec/deg C
      TCG = 0.001
C Diffusion coeff for water vapor in the soil (cm2/sec)
      DCG = 0.2
C Precip multiplication factor during snow
      SFC = 0.76
C Precip multiplication factor during rain
      RCF = 1.0
C Liquid water holding capacity (basic value)
      PLWHC = 0.03
C Liquid water holding capacity (zero density)
      PLWMAX = 0.10 
C Density above which basic liquid water capacity is used
      PLWDEN = 0.20
C Wind function coeff (mm/mb/km)
      FUCOEF = 0.002
      COEFKE = 0.006
C Extinction coefficient (=CV*DENSITY*SQRT(1/GRAIN SIZE)
      CV = 1.2
      RICRIT = 0.4
C Roughness height in cm for determining bulk transfer coefficient
      ZO = ZLVL

C For a new snowpack, define number of layers, liquid wc, density, temp
      NEWPACK = 0

      IF ( (SNOWNG) .OR. (FRZGRA) ) THEN
C New snowpack case         
         IF ((SNEQV - SN_NEW) == 0.0) THEN 
            NEWPACK = 1
            NSNOW = 1
            THICK0 = SN_NEW * 100
            DEN = SNDENS
            TEMP = MIN((SFCTMP-273.15),0.0)
            WATER0 = 0
         ENDIF
      ENDIF

C Arrays for (up to) 100 snow layers, temperature, density, liquid w/c, 
C SWE as well as global variables: number of layers (N). Move 
C definitions from inside PEMB to MAIN_DRIVER
      WRITE(*,*)'ABOVE PEMB'
      CALL PEMB(ZO,MODEL_TYPE,ITOPT,ITS,IPUNCH,IGRAD,IQAE,IFU,NN,
     &     DELTAT,DTOUT,THEDA,TOLER,GRMAX,DTONE,TEXP,THICK,CTHICK,
     &     ADJQA,HEIGHT,D_O,PA,X,DTG,TCG,DCG,SFC,RCF,PLWHC,PLWMAX,
     &     PLWDEN,FUCOEF,COEFKE,CV,RICRIT,C1,C2,C3,C4,DMETA,C5,
     &     CW1,CW2,CW3,CW4,G1,G2,G3,MM,DEN,TEMP,STC(1),NEWPACK,THICK0,
     &     WATER0,SFCTMP,Q2,SFCPRS,SFCSPD,SOLDN,LWDN,
     &     SNOALB,PRCP,SNOWNG,
C State and output variables
     &     NSNOW,DSNOW,PSNOW,RTTSNOW,RTTDTSNOW,WTSNOW,WTDTSNOW,TTSNOW,
     &     TTDTSNOW,SCOUTSNOW,VAPOURSNOW,QGSNOW)
      WRITE(*,*)'BELOW PEMB'
C ----------------------------------------------------------------------
C Obtain sums over N layers and convert necessary units
C SNOMLT (m), ESD (m), ESNOW(mm/s), SSOIL (W/m2)
C ---------------------------------------------------------------------- 
      SSOIL = QGSNOW * (41840/DT)
      SNOMLT = SCOUTSNOW * 0.001

C Accumulate layer's ice depth and liquid water depth SWE (ESD)
      ESD = 0
      DO I=1,NSNOW
         ESD = ESD + DSNOW(I)*PSNOW(I) + WTDTSNOW(I)
      END DO
      ESD = ESD *0.01

      IF (NSNOW.EQ.1) THEN
         TPACK = TTDTSNOW(1)
      ELSE
         IF (NSNOW.GT.1) THEN
            DO N = 2, NSNOW
               TPACK = TPACK + TTDTSNOW(N)
            END DO
         END IF
          TPACK = TPACK / (NSNOW - 1)
       END IF

C ----------------------------------------------------------------------
C ESNOW represents the total vapour transfer to (positive) and from
C (negative) the snowpack; e.g. evap and sublimation are negative
C ----------------------------------------------------------------------
      ESNOW = (VAPOURSNOW * 10) / DT
C ----------------------------------------------------------------------
C Partial snow cover, composite skin temperature with T12
C ----------------------------------------------------------------------
      DSOIL = -(0.5 * ZSOIL(1))
      DTOT = SNOWH + DSOIL
      DENOM = 1.0 + DF1 / (DTOT * RR * RCH)
      T12A = ( (FDOWN-FLX1-FLX2-
     &       (EMISSNOW*SNCOVR+EMISS*(1.0-SNCOVR))*SIGMA*T24)/RCH
     &       + TH2 - SFCTMP - ETANRG/RCH ) / RR

      T12B = DF1 * STC(1) / (DTOT * RR * RCH)
      T12 = (SFCTMP + T12A + T12B) / DENOM    

      IF (T12.GT.TFREEZ) THEN
         T1 = TTDTSNOW(1) * SNCOVR**SNOEXP + T12 * (1.0 - SNCOVR**SNOEXP)
      ELSE
         T1 = TTDTSNOW(1)
      ENDIF

C ----------------------------------------------------------------------
C SMFLX RETURNS UPDATED SOIL MOISTURE VALUES.  IN THIS, THE SNOW PACK
C CASE, ETA1 IS NOT USED IN CALCULATION OF EVAP.
C ----------------------------------------------------------------------
c      ETP1 = 0.0
      IF (ICE .NE. 1) THEN
         IF (MODEL_TYPE.EQ.0.OR.MODEL_TYPE == 2) THEN
            CALL SMFLX (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &           SH2O,SLOPE,KDT,FRZFACT,
     &           SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &           SHDFAC,CMCMAX,
     &           RUNOFF1,RUNOFF2,RUNOFF3,
     &           EDIR1,EC1,ET1,DRIP)
         ENDIF
         IF (MODEL_TYPE.EQ.1.OR.MODEL_TYPE == 3) THEN
            CALL SMFLX_SAC (SMC,NSOIL,CMC,DT,PRCP1,ZSOIL,
     &           SH2O,SLOPE,KDT,FRZFACT,SNCOVR,SMCREF,
     &           SMCMAX,BEXP,SMCWLT,DKSAT,DWSAT,
     &           SHDFAC,CMCMAX,SFCTMP,ETP1,ETA1,
     &           RUNOFF1,RUNOFF2,RUNOFF3, 
     &           EDIR1,EC1,ET1,DRIP,UZTWM,UZFWM,UZK,PCTIM,ADIMP,
     &           RIVA,ZPERC,REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
     &           SIDE,RSERV,UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC)
         ENDIF
      ENDIF


C ----------------------------------------------------------------------
C BEFORE CALL SHFLX IN THIS SNOWPACK CASE, SET ZZ1 AND YY ARGUMENTS TO
C SPECIAL VALUES THAT ENSURE THAT GROUND HEAT FLUX CALCULATED IN SHFLX
C MATCHES THAT ALREADY COMPUTED FOR BELOW THE SNOWPACK, THUS THE SFC
C HEAT FLUX TO BE COMPUTED IN SHFLX WILL EFFECTIVELY BE THE FLUX AT THE
C SNOW TOP SURFACE.  T11 IS A DUMMY ARGUEMENT SO WE WILL NOT USE THE
C SKIN TEMP VALUE AS REVISED BY SHFLX.
C ----------------------------------------------------------------------
      ZZ1 = 1.0
      YY = STC(1)-0.5*SSOIL*ZSOIL(1)*ZZ1/DF1
      T11 = TPACK

C ----------------------------------------------------------------------
C SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX 
C (SSOIL1) AND THE SKIN TEMP (T11) OUTPUT FROM THIS SHFLX CALL ARE NOT
C USED  IN ANY SUBSEQUENT CALCULATIONS. RATHER, THEY ARE DUMMY VARIABLES
C HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT FLUX ARE
C UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
C ----------------------------------------------------------------------
      CALL SHFLX (SSOIL1,STC,SMC,SMCMAX,NSOIL,T11,DT,YY,ZZ1,ZSOIL,
     &            TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &            QUARTZ,CSOIL)


C ----------------------------------------------------------------------
C END SUBROUTINE TR_19
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SNOWPACK (ESD,DTSEC,SNOWH,SNDENS,TSNOW,TSOIL)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SNOWPACK
C ----------------------------------------------------------------------
C CALCULATE COMPACTION OF SNOWPACK UNDER CONDITIONS OF INCREASING SNOW
C DENSITY, AS OBTAINED FROM AN APPROXIMATE SOLUTION OF E. ANDERSON'S
C DIFFERENTIAL EQUATION (3.29), NOAA TECHNICAL REPORT NWS 19, BY VICTOR
C KOREN, 03/25/95.
C ----------------------------------------------------------------------
C ESD     WATER EQUIVALENT OF SNOW (M)
C DTSEC   TIME STEP (SEC)
C SNOWH   SNOW DEPTH (M)
C SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
C TSNOW   SNOW SURFACE TEMPERATURE (K)
C TSOIL   SOIL SURFACE TEMPERATURE (K)
C
C SUBROUTINE WILL RETURN NEW VALUES OF SNOWH AND SNDENS
C ----------------------------------------------------------------------
      INTEGER IPOL, J

      REAL BFAC,C1,C2,SNDENS,DSX,DTHR,DTSEC,DW,SNOWHC,SNOWH,PEXP,TAVGC,
     &     TSNOW,TSNOWC,TSOIL,TSOILC,ESD,ESDC,ESDCX,G,KN

      PARAMETER(C1 = 0.01, C2=21.0, G=9.81, KN=4000.0)

C ----------------------------------------------------------------------
C CONVERSION INTO SIMULATION UNITS
C ----------------------------------------------------------------------
      SNOWHC = SNOWH*100.
      ESDC = ESD*100.
      DTHR = DTSEC/3600.
      TSNOWC = TSNOW-273.15
      TSOILC = TSOIL-273.15

C ----------------------------------------------------------------------
C CALCULATING OF AVERAGE TEMPERATURE OF SNOW PACK
C ----------------------------------------------------------------------
      TAVGC = 0.5*(TSNOWC+TSOILC)                                    

C ----------------------------------------------------------------------
C CALCULATING OF SNOW DEPTH AND DENSITY AS A RESULT OF COMPACTION
C  SNDENS=DS0*(EXP(BFAC*ESD)-1.)/(BFAC*ESD)
C  BFAC=DTHR*C1*EXP(0.08*TAVGC-C2*DS0)
C NOTE: BFAC*ESD IN SNDENS EQN ABOVE HAS TO BE CAREFULLY TREATED
C NUMERICALLY BELOW:
C   C1 IS THE FRACTIONAL INCREASE IN DENSITY (1/(CM*HR)) 
C   C2 IS A CONSTANT (CM3/G) KOJIMA ESTIMATED AS 21 CMS/G
C ----------------------------------------------------------------------
      IF (ESDC .GT. 1.E-2) THEN
        ESDCX = ESDC
      ELSE
        ESDCX = 1.E-2
      ENDIF
      BFAC = DTHR*C1*EXP(0.08*TAVGC-C2*SNDENS)

C      DSX = SNDENS*((DEXP(BFAC*ESDC)-1.)/(BFAC*ESDC))
C ----------------------------------------------------------------------
C THE FUNCTION OF THE FORM (e**x-1)/x IMBEDDED IN ABOVE EXPRESSION
C FOR DSX WAS CAUSING NUMERICAL DIFFICULTIES WHEN THE DENOMINATOR "x"
C (I.E. BFAC*ESDC) BECAME ZERO OR APPROACHED ZERO (DESPITE THE FACT THAT
C THE ANALYTICAL FUNCTION (e**x-1)/x HAS A WELL DEFINED LIMIT AS 
C "x" APPROACHES ZERO), HENCE BELOW WE REPLACE THE (e**x-1)/x 
C EXPRESSION WITH AN EQUIVALENT, NUMERICALLY WELL-BEHAVED 
C POLYNOMIAL EXPANSION.
C
C NUMBER OF TERMS OF POLYNOMIAL EXPANSION, AND HENCE ITS ACCURACY, 
C IS GOVERNED BY ITERATION LIMIT "IPOL".
C      IPOL GREATER THAN 9 ONLY MAKES A DIFFERENCE ON DOUBLE
C            PRECISION (RELATIVE ERRORS GIVEN IN PERCENT %).
C       IPOL=9, FOR REL.ERROR <~ 1.6 E-6 % (8 SIGNIFICANT DIGITS)
C       IPOL=8, FOR REL.ERROR <~ 1.8 E-5 % (7 SIGNIFICANT DIGITS)
C       IPOL=7, FOR REL.ERROR <~ 1.8 E-4 % ...
C ----------------------------------------------------------------------
      IPOL = 4
      PEXP = 0.
      DO J = IPOL,1,-1
C        PEXP = (1. + PEXP)*BFAC*ESDC/REAL(J+1) 
        PEXP = (1. + PEXP)*BFAC*ESDCX/REAL(J+1) 
      END DO
      PEXP = PEXP + 1.

      DSX = SNDENS*(PEXP)
C ----------------------------------------------------------------------
C ABOVE LINE ENDS POLYNOMIAL SUBSTITUTION
C ----------------------------------------------------------------------
C     END OF KOREAN FORMULATION

C     BASE FORMULATION (COGLEY ET AL., 1990)
C     CONVERT DENSITY FROM G/CM3 TO KG/M3
C       DSM=SNDENS*1000.0
 
C       DSX=DSM+DTSEC*0.5*DSM*G*ESD/
C    &      (1E7*EXP(-0.02*DSM+KN/(TAVGC+273.16)-14.643))
 
C     CONVERT DENSITY FROM KG/M3 TO G/CM3
C       DSX=DSX/1000.0

C     END OF COGLEY ET AL. FORMULATION 

C ----------------------------------------------------------------------
C SET UPPER/LOWER LIMIT ON SNOW DENSITY
C ----------------------------------------------------------------------
      IF (DSX .GT. 0.40) DSX = 0.40
      IF (DSX .LT. 0.05) DSX = 0.05
      SNDENS = DSX
C ----------------------------------------------------------------------
C UPDATE OF SNOW DEPTH AND DENSITY DEPENDING ON LIQUID WATER DURING
C SNOWMELT.  ASSUMED THAT 13% OF LIQUID WATER CAN BE STORED IN SNOW PER
C DAY DURING SNOWMELT TILL SNOW DENSITY 0.40.
C ----------------------------------------------------------------------
      IF (TSNOWC .GE. 0.) THEN
        DW = 0.13*DTHR/24.
        SNDENS = SNDENS*(1.-DW)+DW
        IF (SNDENS .GT. 0.40) SNDENS = 0.40
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE SNOW DEPTH (CM) FROM SNOW WATER EQUIVALENT AND SNOW DENSITY.
C CHANGE SNOW DEPTH UNITS TO METERS
C ----------------------------------------------------------------------
      SNOWHC = ESDC/SNDENS
      SNOWH = SNOWHC*0.01

C ----------------------------------------------------------------------
C END SUBROUTINE SNOWPACK
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE SNOWZ0 (SNCOVR,Z0)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SNOWZ0
C ----------------------------------------------------------------------
C CALCULATE TOTAL ROUGHNESS LENGTH OVER SNOW
C SNCOVR  FRACTIONAL SNOW COVER
C Z0      ROUGHNESS LENGTH (m)
C Z0S     SNOW ROUGHNESS LENGTH:=0.001 (m)
C ----------------------------------------------------------------------
      REAL SNCOVR, Z0, Z0S
c      PARAMETER (Z0S=0.001)
      
C CURRENT NOAH LSM CONDITION - MBEK, 09-OCT-2001
      Z0S = Z0
C
      Z0 = (1-SNCOVR)*Z0 + SNCOVR*Z0S
C ----------------------------------------------------------------------
C END SUBROUTINE SNOWZ0
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SNOW_NEW (TEMP,NEWSN,SNOWH,SNDENS)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SNOW_NEW
C ----------------------------------------------------------------------
C CALCULATE SNOW DEPTH AND DENSITITY TO ACCOUNT FOR THE NEW SNOWFALL.
C NEW VALUES OF SNOW DEPTH & DENSITY RETURNED.
C
C TEMP    AIR TEMPERATURE (K)
C NEWSN   NEW SNOWFALL (M)
C SNOWH   SNOW DEPTH (M)
C SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
C ----------------------------------------------------------------------
      REAL SNDENS
      REAL DSNEW
      REAL SNOWHC
      REAL HNEWC
      REAL SNOWH
      REAL NEWSN
      REAL NEWSNC
      REAL TEMP 
      REAL TEMPC
      
C ----------------------------------------------------------------------
C CONVERSION INTO SIMULATION UNITS      
C ----------------------------------------------------------------------
      SNOWHC = SNOWH*100.
      NEWSNC = NEWSN*100.
      TEMPC = TEMP-273.15
      
C ----------------------------------------------------------------------
C CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE
C EQUATION FROM GOTTLIB L. 'A GENERAL RUNOFF MODEL FOR SNOWCOVERED
C AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE,
C VEMADOLEN, SWEDEN, 1980, 172-177PP.
C-----------------------------------------------------------------------
      IF (TEMPC .LE. -15.) THEN
        DSNEW = 0.05
      ELSE                                                      
        DSNEW = 0.05+0.0017*(TEMPC+15.)**1.5
      ENDIF
      
C ----------------------------------------------------------------------
C ADJUSTMENT OF SNOW DENSITY DEPENDING ON NEW SNOWFALL      
C ----------------------------------------------------------------------
      HNEWC = NEWSNC/DSNEW
      SNDENS = (SNOWHC*SNDENS+HNEWC*DSNEW)/(SNOWHC+HNEWC)
      SNOWHC = SNOWHC+HNEWC
      SNOWH = SNOWHC*0.01
      
C ----------------------------------------------------------------------
C END SUBROUTINE SNOW_NEW
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SRT (RHSTT,EDIR,ET,SH2O,SH2OA,NSOIL,PCPDRP,
     &                ZSOIL,DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1, 
     &                RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZX,SICE,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SRT
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C WATER DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
C COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER CVFRZ      
      INTEGER IALP1
      INTEGER IOHINF
      INTEGER J
      INTEGER JJ      
      INTEGER K
      INTEGER KS
      INTEGER NSOIL

      REAL ACRT
      REAL AI(NSOLD)
      REAL BEXP
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL DD
      REAL DDT
      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DENOM2
      REAL DICE
      REAL DKSAT
      REAL DMAX(NSOLD)
      REAL DSMDZ
      REAL DSMDZ2
      REAL DT
      REAL DT1
      REAL DWSAT
      REAL EDIR
      REAL ET(NSOIL)
      REAL FCR
      REAL FRZX
      REAL INFMAX
      REAL KDT
      REAL MXSMC
      REAL MXSMC2
      REAL NUMER
      REAL PCPDRP
      REAL PDDUM
      REAL PX
      REAL RHSTT(NSOIL)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL SH2O(NSOIL)
      REAL SH2OA(NSOIL)
      REAL SICE(NSOIL)
      REAL SICEMAX
      REAL SLOPE
      REAL SLOPX
      REAL SMCAV
      REAL SMCMAX
      REAL SMCWLT
      REAL SSTT
      REAL SUM
      REAL VAL
      REAL WCND
      REAL WCND2
      REAL WDF
      REAL WDF2
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C REFERENCE FROZEN GROUND PARAMETER, CVFRZ, IS A SHAPE PARAMETER OF
C AREAL DISTRIBUTION FUNCTION OF SOIL ICE CONTENT WHICH EQUALS 1/CV.
C CV IS A COEFFICIENT OF SPATIAL VARIATION OF SOIL ICE CONTENT.  BASED
C ON FIELD DATA CV DEPENDS ON AREAL MEAN OF FROZEN DEPTH, AND IT CLOSE
C TO CONSTANT = 0.6 IF AREAL MEAN FROZEN DEPTH IS ABOVE 20 CM.  THAT IS
C WHY PARAMETER CVFRZ = 3 (INT{1/0.6*0.6}).
C CURRENT LOGIC DOESN'T ALLOW CVFRZ BE BIGGER THAN 3
C ----------------------------------------------------------------------
        PARAMETER(CVFRZ = 3)
        
C ----------------------------------------------------------------------
C DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF.  INCLUDE THE
C INFILTRATION FORMULE FROM SCHAAKE AND KOREN MODEL.
C MODIFIED BY Q DUAN
C ----------------------------------------------------------------------
      IOHINF=1

C ----------------------------------------------------------------------
C LET SICEMAX BE THE GREATEST, IF ANY, FROZEN WATER CONTENT WITHIN SOIL
C LAYERS.
C ----------------------------------------------------------------------
      SICEMAX = 0.0
      DO KS=1,NSOIL
       IF (SICE(KS) .GT. SICEMAX) SICEMAX = SICE(KS)
      END DO

C ----------------------------------------------------------------------
C DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF
C ----------------------------------------------------------------------
      PDDUM = PCPDRP
      RUNOFF1 = 0.0
      IF (PCPDRP .NE. 0.0) THEN

C ----------------------------------------------------------------------
C MODIFIED BY Q. DUAN, 5/16/94
C ----------------------------------------------------------------------
C        IF (IOHINF .EQ. 1) THEN

        DT1 = DT/86400.
        SMCAV = SMCMAX - SMCWLT
        DMAX(1)=-ZSOIL(1)*SMCAV

C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C ----------------------------------------------------------------------
        DICE = -ZSOIL(1) * SICE(1)
          
        DMAX(1)=DMAX(1)*(1.0 - (SH2OA(1)+SICE(1)-SMCWLT)/SMCAV)
        DD=DMAX(1)

        DO KS=2,NSOIL
          
C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C ----------------------------------------------------------------------
          DICE = DICE + ( ZSOIL(KS-1) - ZSOIL(KS) ) * SICE(KS)
         
          DMAX(KS) = (ZSOIL(KS-1)-ZSOIL(KS))*SMCAV
          DMAX(KS) = DMAX(KS)*(1.0 - (SH2OA(KS)+SICE(KS)-SMCWLT)/SMCAV)
          DD = DD+DMAX(KS)
        END DO

C ----------------------------------------------------------------------
C VAL = (1.-EXP(-KDT*SQRT(DT1)))
C IN BELOW, REMOVE THE SQRT IN ABOVE
C ----------------------------------------------------------------------
        VAL = (1.-EXP(-KDT*DT1))
        DDT = DD*VAL
        PX = PCPDRP*DT  
        IF (PX .LT. 0.0) PX = 0.0
        INFMAX = (PX*(DDT/(PX+DDT)))/DT
          
C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C REDUCTION OF INFILTRATION BASED ON FROZEN GROUND PARAMETERS
C ----------------------------------------------------------------------
        FCR = 1. 
        IF (DICE .GT. 1.E-2) THEN 
          ACRT = CVFRZ * FRZX / DICE 
          SUM = 1.
          IALP1 = CVFRZ - 1 
          DO J = 1,IALP1
            K = 1
            DO JJ = J+1,IALP1
              K = K * JJ
            END DO
            SUM = SUM + (ACRT ** ( CVFRZ-J)) / FLOAT (K) 
          END DO
          FCR = 1. - EXP(-ACRT) * SUM 
        ENDIF 
        INFMAX = INFMAX * FCR

C ----------------------------------------------------------------------
C CORRECTION OF INFILTRATION LIMITATION:
C IF INFMAX .LE. HYDROLIC CONDUCTIVITY ASSIGN INFMAX THE VALUE OF
C HYDROLIC CONDUCTIVITY
C ----------------------------------------------------------------------
C         MXSMC = MAX ( SH2OA(1), SH2OA(2) ) 
        MXSMC = SH2OA(1)

        CALL WDFCND (WDF,WCND,MXSMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &               SICEMAX)

        INFMAX = MAX(INFMAX,WCND)
        INFMAX = MIN(INFMAX,PX)

        IF (PCPDRP .GT. INFMAX) THEN
          RUNOFF1 = PCPDRP - INFMAX
          PDDUM = INFMAX
        ENDIF

      ENDIF

C ----------------------------------------------------------------------
C TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN LINE
C BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
C 'MXSMC = MAX(SH2OA(1), SH2OA(2))'
C ----------------------------------------------------------------------
      MXSMC = SH2OA(1)

      CALL WDFCND (WDF,WCND,MXSMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &             SICEMAX)
 
C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1. / ( -.5 * ZSOIL(2) )
      AI(1) = 0.0
      BI(1) = WDF * DDZ / ( -ZSOIL(1) )
      CI(1) = -BI(1)

C ----------------------------------------------------------------------
C CALC RHSTT FOR THE TOP LAYER AFTER CALC'NG THE VERTICAL SOIL MOISTURE
C GRADIENT BTWN THE TOP AND NEXT TO TOP LAYERS.
C ----------------------------------------------------------------------
      DSMDZ = ( SH2O(1) - SH2O(2) ) / ( -.5 * ZSOIL(2) )
      RHSTT(1) = (WDF * DSMDZ + WCND - PDDUM + EDIR + ET(1))/ZSOIL(1)
      SSTT = WDF * DSMDZ + WCND + EDIR + ET(1)

C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0

C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        DENOM2 = (ZSOIL(K-1) - ZSOIL(K))
        IF (K .NE. NSOIL) THEN
          SLOPX = 1.

C ----------------------------------------------------------------------
C AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN
C LINE BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
C 'MXSMC2 = MAX (SH2OA(K), SH2OA(K+1))'
C ----------------------------------------------------------------------
          MXSMC2 = SH2OA(K)

          CALL WDFCND (WDF2,WCND2,MXSMC2,SMCMAX,BEXP,DKSAT,DWSAT,
     &                 SICEMAX)

C ----------------------------------------------------------------------
C CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALC'NG RHSTT
C ----------------------------------------------------------------------
          DENOM = (ZSOIL(K-1) - ZSOIL(K+1))
          DSMDZ2 = (SH2O(K) - SH2O(K+1)) / (DENOM * 0.5)

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
C ----------------------------------------------------------------------
          DDZ2 = 2.0 / DENOM
          CI(K) = -WDF2 * DDZ2 / DENOM2
        ELSE

C ----------------------------------------------------------------------
C SLOPE OF BOTTOM LAYER IS INTRODUCED
C ----------------------------------------------------------------------
          SLOPX = SLOPE

C ----------------------------------------------------------------------
C RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR
C THIS LAYER
C ----------------------------------------------------------------------
          CALL WDFCND (WDF2,WCND2,SH2OA(NSOIL),SMCMAX,BEXP,DKSAT,DWSAT,
     &                 SICEMAX)

C ----------------------------------------------------------------------
C CALC A PARTIAL PRODUCT FOR LATER USE IN CALC'NG RHSTT
C ----------------------------------------------------------------------
          DSMDZ2 = 0.0

C ----------------------------------------------------------------------
C SET MATRIX COEF CI TO ZERO
C ----------------------------------------------------------------------
          CI(K) = 0.0
        ENDIF

C ----------------------------------------------------------------------
C CALC RHSTT FOR THIS LAYER AFTER CALC'NG ITS NUMERATOR
C ----------------------------------------------------------------------
        NUMER = (WDF2 * DSMDZ2) + SLOPX * WCND2 - (WDF * DSMDZ) 
     &    - WCND + ET(K)
        RHSTT(K) = NUMER / (-DENOM2)

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
C ----------------------------------------------------------------------
        AI(K) = -WDF * DDZ / DENOM2
        BI(K) = -( AI(K) + CI(K) )

C ----------------------------------------------------------------------
C RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
C RUNOFF2:  SUB-SURFACE OR BASEFLOW RUNOFF
C ----------------------------------------------------------------------
        IF (K .EQ. NSOIL) THEN
          RUNOFF2 = SLOPX * WCND2
        ENDIF

        IF (K .NE. NSOIL) THEN
          WDF = WDF2
          WCND = WCND2
          DSMDZ = DSMDZ2
          DDZ = DDZ2
        ENDIF
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE SRT
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE SSTEP (SH2OOUT,SH2OIN,CMC,RHSTT,RHSCT,DT,
     &                  NSOIL,SMCMAX,CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,
     &                  AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SSTEP
C ----------------------------------------------------------------------
C CALCULATE/UPDATE SOIL MOISTURE CONTENT VALUES AND CANOPY MOISTURE
C CONTENT VALUES.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K 
      INTEGER KK11
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL CIin(NSOLD)
      REAL CMC
      REAL CMCMAX
      REAL DDZ
      REAL DT
      REAL RHSCT
      REAL RHSTT(NSOIL)
      REAL RHSTTin(NSOIL)
      REAL RUNOFF3
      REAL SH2OIN(NSOIL)
      REAL SH2OOUT(NSOIL)
      REAL SICE(NSOIL)
      REAL SMC(NSOIL)
      REAL SMCMAX
      REAL STOT
      REAL WPLUS
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C CREATE 'AMOUNT' VALUES OF VARIABLES TO BE INPUT TO THE
C TRI-DIAGONAL MATRIX ROUTINE.
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTT(K) = RHSTT(K) * DT
        AI(K) = AI(K) * DT
        BI(K) = 1. + BI(K) * DT
        CI(K) = CI(K) * DT
      END DO

C ----------------------------------------------------------------------
C COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTTin(K) = RHSTT(K)
      END DO
      DO K = 1,NSOLD
        CIin(K) = CI(K)
      END DO

C ----------------------------------------------------------------------
C CALL ROSR12 TO SOLVE THE TRI-DIAGONAL MATRIX
C ----------------------------------------------------------------------
      CALL ROSR12 (CI,AI,BI,CIin,RHSTTin,RHSTT,NSOIL)

C ----------------------------------------------------------------------
C SUM THE PREVIOUS SMC VALUE AND THE MATRIX SOLUTION TO GET A
C NEW VALUE.  MIN ALLOWABLE VALUE OF SMC WILL BE 0.02.
C RUNOFF3: RUNOFF WITHIN SOIL LAYERS
C ----------------------------------------------------------------------
      WPLUS = 0.0
      RUNOFF3 = 0.
      DDZ = -ZSOIL(1)
      
      DO K = 1,NSOIL
        IF (K .NE. 1) DDZ = ZSOIL(K - 1) - ZSOIL(K)
        SH2OOUT(K) = SH2OIN(K) + CI(K) + WPLUS / DDZ

        STOT = SH2OOUT(K) + SICE(K)
        IF (STOT .GT. SMCMAX) THEN
          IF (K .EQ. 1) THEN
            DDZ = -ZSOIL(1)
          ELSE
            KK11 = K - 1
            DDZ = -ZSOIL(K) + ZSOIL(KK11)
          ENDIF
          WPLUS = (STOT-SMCMAX) * DDZ
        ELSE
          WPLUS = 0.
        ENDIF
        SMC(K) = MAX ( MIN(STOT,SMCMAX),0.02 )
        SH2OOUT(K) = MAX((SMC(K)-SICE(K)),0.0)
      END DO

      RUNOFF3 = WPLUS

C ----------------------------------------------------------------------
C UPDATE CANOPY WATER CONTENT/INTERCEPTION (CMC).  CONVERT RHSCT TO 
C AN 'AMOUNT' VALUE AND ADD TO PREVIOUS CMC VALUE TO GET NEW CMC.
C ----------------------------------------------------------------------
      CMC = CMC + DT * RHSCT
      IF (CMC .LT. 1.E-20) CMC=0.0
      CMC = MIN(CMC,CMCMAX)

C ----------------------------------------------------------------------
C END SUBROUTINE SSTEP
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE TBND (TU,TB,ZSOIL,ZBOT,K,NSOIL,TBND1)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TBND
C ----------------------------------------------------------------------
C CALCULATE TEMPERATURE ON THE BOUNDARY OF THE LAYER BY INTERPOLATION OF
C THE MIDDLE LAYER TEMPERATURES
C ----------------------------------------------------------------------
      INTEGER NSOIL
      INTEGER K

      REAL TBND1
      REAL T0
      REAL TU
      REAL TB
      REAL ZB
      REAL ZBOT
      REAL ZUP
      REAL ZSOIL (NSOIL)

      PARAMETER(T0 = 273.15)

C ----------------------------------------------------------------------
C USE SURFACE TEMPERATURE ON THE TOP OF THE FIRST LAYER
C ----------------------------------------------------------------------
      IF (K .EQ. 1) THEN
        ZUP = 0.
      ELSE
        ZUP = ZSOIL(K-1)
      ENDIF

C ----------------------------------------------------------------------
C USE DEPTH OF THE CONSTANT BOTTOM TEMPERATURE WHEN INTERPOLATE
C TEMPERATURE INTO THE LAST LAYER BOUNDARY
C ----------------------------------------------------------------------
      IF (K .EQ. NSOIL) THEN
        ZB = 2.*ZBOT-ZSOIL(K)
      ELSE
        ZB = ZSOIL(K+1)
      ENDIF

C ----------------------------------------------------------------------
C LINEAR INTERPOLATION BETWEEN THE AVERAGE LAYER TEMPERATURES
C ----------------------------------------------------------------------
      TBND1 = TU+(TB-TU)*(ZUP-ZSOIL(K))/(ZUP-ZB)
      
C ----------------------------------------------------------------------
C END SUBROUTINE TBND
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE TDFCND ( DF, SMC, QZ,  SMCMAX, SH2O)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TDFCND
C ----------------------------------------------------------------------
C CALCULATE THERMAL DIFFUSIVITY AND CONDUCTIVITY OF THE SOIL FOR A GIVEN
C POINT AND TIME.
C ----------------------------------------------------------------------
C PETERS-LIDARD APPROACH (PETERS-LIDARD et al., 1998)
C June 2001 CHANGES: FROZEN SOIL CONDITION.
C ----------------------------------------------------------------------
       REAL DF
       REAL GAMMD
       REAL THKDRY
       REAL AKE
       REAL THKICE
       REAL THKO
       REAL THKQTZ
       REAL THKSAT
       REAL THKS
       REAL THKW
       REAL QZ
       REAL SATRATIO
       REAL SH2O
       REAL SMC
       REAL SMCMAX
       REAL XU
       REAL XUNFROZ

C ----------------------------------------------------------------------
C WE NOW GET QUARTZ AS AN INPUT ARGUMENT (SET IN ROUTINE REDPRM):
C      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52, 
C     &             0.35, 0.60, 0.40, 0.82/
C ----------------------------------------------------------------------
C IF THE SOIL HAS ANY MOISTURE CONTENT COMPUTE A PARTIAL SUM/PRODUCT
C OTHERWISE USE A CONSTANT VALUE WHICH WORKS WELL WITH MOST SOILS
C ----------------------------------------------------------------------
C  THKW ......WATER THERMAL CONDUCTIVITY
C  THKQTZ ....THERMAL CONDUCTIVITY FOR QUARTZ
C  THKO ......THERMAL CONDUCTIVITY FOR OTHER SOIL COMPONENTS
C  THKS ......THERMAL CONDUCTIVITY FOR THE SOLIDS COMBINED(QUARTZ+OTHER)
C  THKICE ....ICE THERMAL CONDUCTIVITY
C  SMCMAX ....POROSITY (= SMCMAX)
C  QZ .........QUARTZ CONTENT (SOIL TYPE DEPENDENT)
C ----------------------------------------------------------------------
C USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).
C
C                                  PABLO GRUNMANN, 08/17/98
C REFS.:
C      FAROUKI, O.T.,1986: THERMAL PROPERTIES OF SOILS. SERIES ON ROCK 
C              AND SOIL MECHANICS, VOL. 11, TRANS TECH, 136 PP.
C      JOHANSEN, O., 1975: THERMAL CONDUCTIVITY OF SOILS. PH.D. THESIS,
C              UNIVERSITY OF TRONDHEIM,
C      PETERS-LIDARD, C. D., ET AL., 1998: THE EFFECT OF SOIL THERMAL 
C              CONDUCTIVITY PARAMETERIZATION ON SURFACE ENERGY FLUXES
C              AND TEMPERATURES. JOURNAL OF THE ATMOSPHERIC SCIENCES,
C              VOL. 55, PP. 1209-1224.
C ----------------------------------------------------------------------
C NEEDS PARAMETERS
C POROSITY(SOIL TYPE):
C      POROS = SMCMAX
C SATURATION RATIO:
      SATRATIO = SMC/SMCMAX

C PARAMETERS  W/(M.K)
      THKICE = 2.2
      THKW = 0.57
      THKO = 2.0
C      IF (QZ .LE. 0.2) THKO = 3.0
      THKQTZ = 7.7
C SOLIDS' CONDUCTIVITY      
      THKS = (THKQTZ**QZ)*(THKO**(1.- QZ))

C UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN))
      XUNFROZ = (SH2O + 1.E-9) / (SMC + 1.E-9)

C UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
      XU=XUNFROZ*SMCMAX 
C SATURATED THERMAL CONDUCTIVITY
      THKSAT = THKS**(1.-SMCMAX)*THKICE**(SMCMAX-XU)*THKW**(XU)

C DRY DENSITY IN KG/M3
      GAMMD = (1. - SMCMAX)*2700.

C DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
      THKDRY = (0.135*GAMMD + 64.7)/(2700. - 0.947*GAMMD)

      IF ( (SH2O + 0.0005) .LT. SMC ) THEN
C FROZEN
              AKE = SATRATIO
      ELSE
C UNFROZEN
C RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
          IF ( SATRATIO .GT. 0.1 ) THEN

C KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT 
C LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
C (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).

              AKE = LOG10(SATRATIO) + 1.0

          ELSE

C USE K = KDRY
              AKE = 0.0

          ENDIF
      ENDIF

C  THERMAL CONDUCTIVITY

       DF = AKE*(THKSAT - THKDRY) + THKDRY

C ----------------------------------------------------------------------
C END SUBROUTINE TDFCND
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE TMPAVG (TAVG,TUP,TM,TDN,ZSOIL,NSOIL,K) 
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE TMPAVG
C ----------------------------------------------------------------------
C CALCULATE SOIL LAYER AVERAGE TEMPERATURE (TAVG) IN FREEZING/THAWING
C LAYER USING UP, DOWN, AND MIDDLE LAYER TEMPERATURES (TUP, TDN, TM),
C WHERE TUP IS AT TOP BOUNDARY OF LAYER, TDN IS AT BOTTOM BOUNDARY OF
C LAYER.  TM IS LAYER PROGNOSTIC STATE TEMPERATURE.
C ----------------------------------------------------------------------
      INTEGER K
      INTEGER NSOIL

      REAL DZ
      REAL DZH
      REAL T0
      REAL TAVG
      REAL TDN
      REAL TM
      REAL TUP
      REAL X0
      REAL XDN
      REAL XUP
      REAL ZSOIL (NSOIL)

      PARAMETER(T0 = 2.7315E2)

C ----------------------------------------------------------------------
      IF (K .EQ. 1) THEN
        DZ = -ZSOIL(1)
      ELSE
        DZ = ZSOIL(K-1)-ZSOIL(K)
      ENDIF

      DZH=DZ*0.5

      IF (TUP .LT. T0) THEN
        IF (TM .LT. T0) THEN
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP, TM, TDN < T0
C ----------------------------------------------------------------------
            TAVG = (TUP + 2.0*TM + TDN)/ 4.0            
          ELSE
C ----------------------------------------------------------------------
C TUP & TM < T0,  TDN >= T0
C ----------------------------------------------------------------------
            X0 = (T0 - TM) * DZH / (TDN - TM)
            TAVG = 0.5 * (TUP*DZH+TM*(DZH+X0)+T0*(2.*DZH-X0)) / DZ
          ENDIF      
        ELSE
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP < T0, TM >= T0, TDN < T0
C ----------------------------------------------------------------------
            XUP  = (T0-TUP) * DZH / (TM-TUP)
            XDN  = DZH - (T0-TM) * DZH / (TDN-TM)
            TAVG = 0.5 * (TUP*XUP+T0*(2.*DZ-XUP-XDN)+TDN*XDN) / DZ
          ELSE
C ----------------------------------------------------------------------
C TUP < T0, TM >= T0, TDN >= T0
C ----------------------------------------------------------------------
            XUP  = (T0-TUP) * DZH / (TM-TUP)
            TAVG = 0.5 * (TUP*XUP+T0*(2.*DZ-XUP)) / DZ
          ENDIF   
        ENDIF
      ELSE
        IF (TM .LT. T0) THEN
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP >= T0, TM < T0, TDN < T0
C ----------------------------------------------------------------------
            XUP  = DZH - (T0-TUP) * DZH / (TM-TUP)
            TAVG = 0.5 * (T0*(DZ-XUP)+TM*(DZH+XUP)+TDN*DZH) / DZ
          ELSE
C ----------------------------------------------------------------------
C TUP >= T0, TM < T0, TDN >= T0
C ----------------------------------------------------------------------
            XUP  = DZH - (T0-TUP) * DZH / (TM-TUP)
            XDN  = (T0-TM) * DZH / (TDN-TM)
            TAVG = 0.5 * (T0*(2.*DZ-XUP-XDN)+TM*(XUP+XDN)) / DZ
          ENDIF   
        ELSE
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP >= T0, TM >= T0, TDN < T0
C ----------------------------------------------------------------------
            XDN  = DZH - (T0-TM) * DZH / (TDN-TM)
            TAVG = (T0*(DZ-XDN)+0.5*(T0+TDN)*XDN) / DZ                 
          ELSE
C ----------------------------------------------------------------------
C TUP >= T0, TM >= T0, TDN >= T0
C ----------------------------------------------------------------------
            TAVG = (TUP + 2.0*TM + TDN) / 4.0
          ENDIF
        ENDIF
      ENDIF
C ----------------------------------------------------------------------
C END SUBROUTINE TMPAVG
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE TRANSP (ET1,NSOIL,ETP1,SMC,STC,CMC,ZSOIL,SHDFAC,SMCWLT,
     &                   CMCMAX,PC,CFACTR,SMCREF,SFCTMP,Q2,NROOT,RTDIS)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TRANSP
C ----------------------------------------------------------------------
C CALCULATE TRANSPIRATION FOR THE VEG CLASS.
C ----------------------------------------------------------------------
      INTEGER I
      INTEGER K
      INTEGER NSOIL
      INTEGER NROOT

      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL DENOM
      REAL ET1(NSOIL)
      REAL ETP1
      REAL ETP1A
      REAL GX (7)
C.....REAL PART(NSOIL)
      REAL PC
      REAL Q2
      REAL RTDIS(NSOIL)
      REAL GTX(NSOIL)
      REAL RTX
      REAL SFCTMP
      REAL SGX
      REAL SHDFAC
      REAL SMC(NSOIL)
      REAL STC(NSOIL)
      REAL SMCREF
      REAL SMCWLT
      REAL ZSOIL(NSOIL)
      REAL STCRT
      REAL DSTCRT

C ----------------------------------------------------------------------
C INITIALIZE PLANT TRANSP TO ZERO FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        ET1(K) = 0.
      END DO

C ----------------------------------------------------------------------
C CALCULATE AN 'ADJUSTED' POTENTIAL TRANSPIRATION
C IF STATEMENT BELOW TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
C NOTE: GX AND OTHER TERMS BELOW REDISTRIBUTE TRANSPIRATION BY LAYER,
C ET(K), AS A FUNCTION OF SOIL MOISTURE AVAILABILITY, WHILE PRESERVING
C TOTAL ETP1A.
C ----------------------------------------------------------------------
      IF (CMC .NE. 0.0) THEN
        ETP1A = SHDFAC * PC * ETP1 * (1.0 - (CMC /CMCMAX) ** CFACTR)
      ELSE
        ETP1A = SHDFAC * PC * ETP1
      ENDIF
      
chelin      SGX = 0.0
chelin        DO I = 1,NROOT
chelin          GX(I) = ( SMC(I) - SMCWLT ) / ( SMCREF - SMCWLT )
chelin          GX(I) = MAX ( MIN ( GX(I), 1. ), 0. )
chelin        SGX = SGX + GX (I)
chelin        END DO
chelin      SGX = SGX / NROOT
      
chelin      DENOM = 0.
chelin      DO I = 1,NROOT
chelin        RTX = RTDIS(I) + GX(I) - SGX
chelin        GX(I) = GX(I) * MAX ( RTX, 0. )
chelin        DENOM = DENOM + GX(I)
chelin      END DO
chelin      IF (DENOM .LE. 0.0) DENOM = 1.

chelin      DO I = 1,NROOT
chelin        ET1(I) = ETP1A * GX(I) / DENOM
chelin      END DO
C seasonal factor for reduced root water uptake
C calculate root zone soil temp
       STCRT=0
       DSTCRT=0
      DO K=1,NROOT
       STCRT=STCRT+ZSOIL(K)*STC(K)
       DSTCRT=DSTCRT+ZSOIL(K)
      ENDDO
       STCRT=STCRT/DSTCRT

        DO K = 1,NROOT
         GTX(K)= MAX(0.0,1.-0.0016* MAX(298.- STCRT,0.0)**2)
        ENDDO

C ----------------------------------------------------------------------
C ABOVE CODE ASSUMES A VERTICALLY UNIFORM ROOT DISTRIBUTION
C CODE BELOW TESTS A VARIABLE ROOT DISTRIBUTION
C ----------------------------------------------------------------------
C      ET1(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * GX * ETP1A
C      ET1(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * ETP1A
C ----------------------------------------------------------------------
C USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
       ET1(1) = RTDIS(1) * ETP1A * GTX(1)
C      ET1(1) = ETP1A * PART(1)
C ----------------------------------------------------------------------
C LOOP DOWN THRU THE SOIL LAYERS REPEATING THE OPERATION ABOVE,
C BUT USING THE THICKNESS OF THE SOIL LAYER (RATHER THAN THE
C ABSOLUTE DEPTH OF EACH LAYER) IN THE FINAL CALCULATION.
C ----------------------------------------------------------------------
       DO K = 2,NROOT
c        GX = ( SMC(K) - SMCWLT ) / ( SMCREF - SMCWLT )
c        GX = MAX ( MIN ( GX, 1. ), 0. )
C TEST CANOPY RESISTANCE
C        GX = 1.0
C        ET1(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*GX*ETP1A
C        ET1(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*ETP1A
C ----------------------------------------------------------------------
C USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
         ET1(K) = RTDIS(K)*ETP1A * GTX(K)
C        ET1(K) = ETP1A*PART(K)
       END DO      
C ----------------------------------------------------------------------
C END SUBROUTINE TRANSP
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE WDFCND (WDF,WCND,SMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &                   SICEMAX)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE WDFCND
C ----------------------------------------------------------------------
C CALCULATE SOIL WATER DIFFUSIVITY AND SOIL HYDRAULIC CONDUCTIVITY.
C ----------------------------------------------------------------------
      REAL BEXP
      REAL DKSAT
      REAL DWSAT
      REAL EXPON
      REAL FACTR1
      REAL FACTR2
      REAL SICEMAX
      REAL SMC
      REAL SMCMAX
      REAL VKwgt
      REAL WCND
      REAL WDF

C ----------------------------------------------------------------------
C     CALC THE RATIO OF THE ACTUAL TO THE MAX PSBL SOIL H2O CONTENT
C ----------------------------------------------------------------------
      SMC = SMC
      SMCMAX = SMCMAX
      FACTR1 = 0.2 / SMCMAX
      FACTR2 = SMC / SMCMAX

C ----------------------------------------------------------------------
C PREP AN EXPNTL COEF AND CALC THE SOIL WATER DIFFUSIVITY
C ----------------------------------------------------------------------
      EXPON = BEXP + 2.0
      WDF = DWSAT * FACTR2 ** EXPON

C ----------------------------------------------------------------------
C FROZEN SOIL HYDRAULIC DIFFUSIVITY.  VERY SENSITIVE TO THE VERTICAL
C GRADIENT OF UNFROZEN WATER. THE LATTER GRADIENT CAN BECOME VERY
C EXTREME IN FREEZING/THAWING SITUATIONS, AND GIVEN THE RELATIVELY 
C FEW AND THICK SOIL LAYERS, THIS GRADIENT SUFFERES SERIOUS 
C TRUNCTION ERRORS YIELDING ERRONEOUSLY HIGH VERTICAL TRANSPORTS OF
C UNFROZEN WATER IN BOTH DIRECTIONS FROM HUGE HYDRAULIC DIFFUSIVITY.  
C THEREFORE, WE FOUND WE HAD TO ARBITRARILY CONSTRAIN WDF 
C --
C VERSION D_10CM: ........  FACTR1 = 0.2/SMCMAX
C WEIGHTED APPROACH...................... PABLO GRUNMANN, 28_SEP_1999.
C ----------------------------------------------------------------------
      IF (SICEMAX .GT. 0.0)  THEN
        VKWGT = 1./(1.+(500.*SICEMAX)**3.)
        WDF = VKWGT*WDF + (1.- VKWGT)*DWSAT*FACTR1**EXPON
      ENDIF

C ----------------------------------------------------------------------
C RESET THE EXPNTL COEF AND CALC THE HYDRAULIC CONDUCTIVITY
C ----------------------------------------------------------------------
      EXPON = (2.0 * BEXP) + 3.0
      WCND = DKSAT * FACTR2 ** EXPON

C ----------------------------------------------------------------------
C END SUBROUTINE WDFCND
C ----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE REDPRM (
     &     VEGTYP,SOILTYP,SLOPETYP,
     &     CFACTR,CMCMAX,RSMAX,TOPT,REFKDT,KDT,SBETA,
     &     SHDFAC,SHDMAX,SHDMIN,RSMIN,RGL,HS,ZBOT,FRZX,
     &     PSISAT,SLOPE,SNUP,SALP,BEXP,DKSAT,DWSAT,
     &     SMCMAX,SMCWLT,SMCREF,SMCDRY,F1,QUARTZ,FXEXP,
     &     RTDIS,SLDPTH,ZSOIL,NROOT,NSOIL,Z0,CZIL,LAI,CSOIL,
     &     PTU,EMISS,EMISSNOW,MODEL_TYPE)

      IMPLICIT NONE
C ----------------------------------------------------------------------
C SUBROUTINE REDPRM
C ----------------------------------------------------------------------
C INTERNALLY SET (DEFAULT VALUESS), OR OPTIONALLY READ-IN VIA NAMELIST
C I/O, ALL SOIL AND VEGETATION PARAMETERS REQUIRED FOR THE EXECUSION OF
C THE NOAH LSM.
C
C OPTIONAL NON-DEFAULT PARAMETERS CAN BE READ IN, ACCOMMODATING UP TO 30
C SOIL, VEG, OR SLOPE CLASSES, IF THE DEFAULT MAX NUMBER OF SOIL, VEG,
C AND/OR SLOPE TYPES IS RESET.
C
C FUTURE UPGRADES OF ROUTINE REDPRM MUST EXPAND TO INCORPORATE SOME OF
C THE EMPIRICAL PARAMETERS OF THE FROZEN SOIL AND SNOWPACK PHYSICS (SUCH
C AS IN ROUTINES FRH2O, SNOWPACK, AND SNOW_NEW) NOT YET SET IN THIS
C REDPRM ROUTINE, BUT RATHER SET IN LOWER LEVEL SUBROUTINES.
C
C SET MAXIMUM NUMBER OF SOIL-, VEG-, AND SLOPETYP IN DATA STATEMENT.
C ----------------------------------------------------------------------
      INTEGER MAX_SLOPETYP
      INTEGER MAX_SOILTYP
      INTEGER MAX_VEGTYP

      PARAMETER(MAX_SLOPETYP = 30)
      PARAMETER(MAX_SOILTYP = 30)
      PARAMETER(MAX_VEGTYP = 30)

C ----------------------------------------------------------------------
C NUMBER OF DEFINED SOIL-, VEG-, AND SLOPETYPS USED.
C ----------------------------------------------------------------------
      INTEGER DEFINED_VEG
      INTEGER DEFINED_SOIL
      INTEGER DEFINED_SLOPE

      DATA DEFINED_VEG/27/
      DATA DEFINED_SOIL/19/
      DATA DEFINED_SLOPE/9/

C ----------------------------------------------------------------------
C  SET-UP SOIL PARAMETERS FOR GIVEN SOIL TYPE
C  INPUT: SOLTYP: SOIL TYPE (INTEGER INDEX)
C  OUTPUT: SOIL PARAMETERS:
C    MAXSMC: MAX SOIL MOISTURE CONTENT (POROSITY)
C    REFSMC: REFERENCE SOIL MOISTURE (ONSET OF SOIL MOISTURE
C            STRESS IN TRANSPIRATION)
C    WLTSMC: WILTING PT SOIL MOISTURE CONTENTS
C    DRYSMC: AIR DRY SOIL MOIST CONTENT LIMITS
C    SATPSI: SATURATED SOIL POTENTIAL
C    SATDK:  SATURATED SOIL HYDRAULIC CONDUCTIVITY
C    BB:     THE 'B' PARAMETER
C    SATDW:  SATURATED SOIL DIFFUSIVITY
C    F11:    USED TO COMPUTE SOIL DIFFUSIVITY/CONDUCTIVITY
C    QUARTZ:  SOIL QUARTZ CONTENT
C ----------------------------------------------------------------------
C SOIL  STATSGO
C TYPE  CLASS
C ----  -------
C   1   SAND
C   2   LOAMY SAND
C   3   SANDY LOAM
C   4   SILT LOAM
C   5   SILT
C   6   LOAM
C   7   SANDY CLAY LOAM
C   8   SILTY CLAY LOAM
C   9   CLAY LOAM
C  10   SANDY CLAY
C  11   SILTY CLAY
C  12   CLAY
C  13   ORGANIC MATERIAL
C  14   WATER
C  15   BEDROCK
C  16   OTHER(land-ice)
C  17   PLAYA
C  18   LAVA
C  19   WHITE SAND
C ----------------------------------------------------------------------

      REAL BB(MAX_SOILTYP)
      REAL BB2(MAX_SOILTYP)
      REAL DRYSMC(MAX_SOILTYP)
      REAL F11(MAX_SOILTYP)
      REAL MAXSMC(MAX_SOILTYP)
      REAL MAXSMC2(MAX_SOILTYP)
      REAL REFSMC(MAX_SOILTYP)
      REAL REFSMC2(MAX_SOILTYP)
      REAL SATPSI(MAX_SOILTYP)
      REAL SATPSI2(MAX_SOILTYP)
      REAL SATDK(MAX_SOILTYP)
      REAL SATDW(MAX_SOILTYP)
      REAL WLTSMC(MAX_SOILTYP)
      REAL WLTSMC2(MAX_SOILTYP)
      REAL QTZ(MAX_SOILTYP)

      REAL BEXP
      REAL DKSAT
      REAL DWSAT
      REAL F1
      REAL PTU
      REAL QUARTZ
      REAL REFSMC1
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL WLTSMC1
      INTEGER MODEL_TYPE
      

C ----------------------------------------------------------------------
C SOIL TEXTURE-RELATED ARRAYS.
CBL Note: all arrays with suffix "2" were derived from Koren, 2003
CBL relationships, which use Cosby et al., 1984, Campbell et al., 1974
CBL and Armstrong, 1978 to derive soil properties from sand and clay
CBL contents.  These will be used when SAC is envoked for consistency.
C ----------------------------------------------------------------------
      DATA MAXSMC2/0.37308, 0.38568, 0.41592, 0.46758, 0.47766, 0.43482,
     &            0.41592, 0.47640, 0.44868, 0.42348, 0.48144, 0.46128,
     &            0.464, 0.000, 0.200, 0.421, 0.457, 0.200,
     &            0.395, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      DATA MAXSMC/0.395, 0.421, 0.434, 0.476, 0.476, 0.439,
     &            0.404, 0.464, 0.465, 0.406, 0.468, 0.457,
     &            0.464, 0.000, 0.200, 0.421, 0.457, 0.200,
     &            0.395, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
C ----------------------------------------------------------------------
      DATA SATPSI2/0.48094695, 0.65051019, 1.34286003, 4.63205997, 
     &             5.89793145, 2.11235130, 1.34286003, 5.72247662,
     &             2.94468452, 1.60962573, 6.45723810, 3.98286611,
     &            0.3548, 0.0000, 0.0350, 0.0363, 0.4677, 0.0350,
     &            0.0350, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     &            0.000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
      DATA SATPSI/0.0350, 0.0363, 0.1413, 0.7586, 0.7586, 0.3548,
     &            0.1349, 0.6166, 0.2630, 0.0977, 0.3236, 0.4677,
     &            0.3548, 0.0000, 0.0350, 0.0363, 0.4677, 0.0350,
     &            0.0350, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     &            0.000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000/
C ----------------------------------------------------------------------
      DATA SATDK /1.7600E-4, 1.4078E-5, 5.2304E-6, 2.8089E-6, 2.8089E-6,
     &            3.3770E-6, 4.4518E-6, 2.0348E-6, 2.4464E-6, 7.2199E-6,
     &            1.3444E-6, 9.7394E-7, 3.3770E-6,       0.0, 1.4078E-5,
     &            1.4078E-5, 9.7394E-7, 1.4078E-5, 1.7600E-4,       0.0,
     &                  0.0,       0.0,       0.0,       0.0,       0.0,
     &                  0.0,       0.0,       0.0,       0.0,       0.0/
C ----------------------------------------------------------------------
      DATA BB2    /3.387, 3.864, 4.500, 4.977, 3.705, 5.772,
     &            7.203, 8.316, 8.316, 9.588, 10.383, 12.132,
     &            5.25,  0.00,  4.05,  4.26, 11.55,  4.05,
     &            4.05,  0.00,  0.00,  0.00,  0.00,  0.00,
     &            0.00,  0.00,  0.00,  0.00,  0.00,  0.00/
      DATA BB    /4.05,  4.26,  4.74,  5.33,  5.33,  5.25,
     &            6.77,  8.72,  8.17, 10.73, 10.39, 11.55,
     &            5.25,  0.00,  4.05,  4.26, 11.55,  4.05,
     &            4.05,  0.00,  0.00,  0.00,  0.00,  0.00,
     &            0.00,  0.00,  0.00,  0.00,  0.00,  0.00/
C ----------------------------------------------------------------------
      DATA QTZ   /0.92, 0.82, 0.60, 0.25, 0.10, 0.40,
     &            0.60, 0.10, 0.35, 0.52, 0.10, 0.25,
     &            0.05, 0.00, 0.07, 0.25, 0.60, 0.52,
     &            0.92, 0.00, 0.00, 0.00, 0.00, 0.00,
     &            0.00, 0.00, 0.00, 0.00, 0.00, 0.00/

C ----------------------------------------------------------------------
C THE FOLLOWING 5 PARAMETERS ARE DERIVED LATER IN REDPRM.F FROM THE SOIL
C DATA, AND ARE JUST GIVEN HERE FOR REFERENCE AND TO FORCE STATIC
C STORAGE ALLOCATION. -DAG LOHMANN, FEB. 2001
C ----------------------------------------------------------------------
C !!!!!!!!!!!!!! The following values in the table are NOT used
C !!!!!!!!!!!!!! and are just given for reference
      DATA REFSMC/0.196, 0.248, 0.282, 0.332, 0.332, 0.301,
     &            0.293, 0.368, 0.361, 0.320, 0.388, 0.389,
     &            0.319, 0.000, 0.116, 0.248, 0.389, 0.116,
     &            0.196, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      DATA REFSMC2/0.15229854, 0.19015085, 0.26621888, 0.34851191,
     &     0.34354197, 0.29455855, 0.28586507, 0.40984752, 0.35636059, 
     &     0.32561179, 0.43177247, 0.40382873, 0.319, 0.000, 0.116, 
     &     0.248, 0.389, 0.116, 0.196, 0.000, 0.196, 0.248, 0.282, 
     &     0.332, 0.332, 0.301, 0.000, 0.000, 0.000, 0.000/
C !!!!!!!!!!!!!! The following values in the table are NOT used
C !!!!!!!!!!!!!! and are just given for reference
      DATA WLTSMC/0.023, 0.028, 0.047, 0.084, 0.084, 0.066,
     &            0.069, 0.120, 0.103, 0.100, 0.126, 0.135,
     &            0.069, 0.000, 0.012, 0.028, 0.135, 0.012,
     &            0.023, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      DATA WLTSMC2/0.03469064, 0.05199094, 0.08743051, 0.14637683, 
     &     0.10712489, 0.13941739, 0.15698002, 0.24386303, 0.21203782,
     &     0.20755672, 0.28488226, 0.28290603, 0.069, 0.000, 0.012, 
     &     0.028, 0.135, 0.012, 0.023, 0.000, 0.000, 0.000, 0.000,
     &     0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
C !!!!!!!!!!!!!! The following values in the table are NOT used
C !!!!!!!!!!!!!! and are just given for reference
      DATA DRYSMC/0.023, 0.028, 0.047, 0.084, 0.084, 0.066,
     &            0.069, 0.120, 0.103, 0.100, 0.126, 0.135,
     &            0.069, 0.000, 0.012, 0.028, 0.135, 0.012,
     &            0.023, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/

C !!!!!!!!!!!!!! The following values in the table are NOT used
C !!!!!!!!!!!!!! and are just given for reference
      DATA SATDW /0.632E-4, 0.517E-5, 0.807E-5, 0.239E-4, 0.239E-4,
     &            0.143E-4, 0.101E-4, 0.236E-4, 0.113E-4, 0.186E-4,
     &            0.966E-5, 0.115E-4, 0.136E-4,      0.0, 0.998E-5,
     &            0.517E-5, 0.115E-4, 0.998E-5, 0.632E-4,      0.0,
     &                 0.0,      0.0,      0.0,      0.0,      0.0,
     &                 0.0,      0.0,      0.0,      0.0,      0.0/
C !!!!!!!!!!!!!! The following values in the table are NOT used
C !!!!!!!!!!!!!! and are just given for reference
      DATA F11  /-1.090, -1.041, -0.568,  0.162,  0.162, -0.327,
     &           -1.535, -1.118, -1.297, -3.211, -1.916, -2.258,
     &           -0.201,  0.000, -2.287, -1.041, -2.258, -2.287,
     &           -1.090,  0.000,  0.000,  0.000,  0.000,  0.000,
     &            0.000,  0.000,  0.000,  0.000,  0.000,  0.000/

C ----------------------------------------------------------------------
C SET-UP VEGETATION PARAMETERS FOR A GIVEN VEGETAION TYPE:
C INPUT: VEGTYP = VEGETATION TYPE (INTEGER INDEX)
C OUPUT: VEGETATION PARAMETERS
C   SHDFAC: VEGETATION GREENNESS FRACTION
C   RSMIN:  MIMIMUM STOMATAL RESISTANCE
C   RGL:    PARAMETER USED IN SOLAR RAD TERM OF
C           CANOPY RESISTANCE FUNCTION
C   HS:     PARAMETER USED IN VAPOR PRESSURE DEFICIT TERM OF
C           CANOPY RESISTANCE FUNCTION
C   SNUP:   THRESHOLD SNOW DEPTH (IN WATER EQUIVALENT M) THAT
C           IMPLIES 100% SNOW COVER
C ----------------------------------------------------------------------
C UMD VEGETATION TYPE
C   1   Evergreen Needleleaf Forest 
C   2   Evergreen Broadleaf Forest
C   3   Deciduous Needleleaf Forest
C   4   Deciduous Broadleaf Forest
C   5   Mixed Cover
C   6   Woodland 
C   7   Wooded Grassland
C   8   Closed Shrubland
C   9   Open Shrubland
C  10   Grassland
C  11   Cropland
C  12   Bare Ground
C  13   Urban and Built-Up
C ----------------------------------------------------------------------

      INTEGER NROOT
      INTEGER NROOT_DATA(MAX_VEGTYP)
      INTEGER NROOT_DATA_SAC(MAX_VEGTYP)
      INTEGER NROOT_DATA_NOAH(MAX_VEGTYP)
      REAL FRZFACT
      REAL HS
      REAL HSTBL(MAX_VEGTYP)
      REAL LAI
      REAL LAI_DATA(MAX_VEGTYP)
      REAL LAI_MAX(MAX_VEGTYP)
      REAL LAI_MIN(MAX_VEGTYP)
      REAL ROOT_A(MAX_VEGTYP)
      REAL ROOT_B(MAX_VEGTYP)
      REAL RTA
      REAL RTB
      REAL PSISAT
      REAL RSMIN
      REAL RGL
      REAL RGLTBL(MAX_VEGTYP)
      REAL RSMTBL(MAX_VEGTYP)
      REAL SHDFAC
      REAL SHDMAX
      REAL SHDMIN
      REAL SNUP
      REAL SNUPX(MAX_VEGTYP)
      REAL Z0
      REAL Z0_DATA(MAX_VEGTYP)
      REAL EMISS_DATA(MAX_VEGTYP)
      REAL EMISSNOW_DATA
      REAL EMISS
      REAL EMISSNOW
      REAL RMAX_DATA(MAX_VEGTYP)
      REAL CROOT_DATA(MAX_VEGTYP)
      REAL D50_DATA(MAX_VEGTYP)


C ----------------------------------------------------------------------
C VEGETATION CLASS-RELATED ARRAYS
C ----------------------------------------------------------------------
cbl      DATA NROOT_DATA_SAC /2,2,2,2,2,2,2,2,2,2,
cbl     &     2,2,1,0,0,0,0,0,0,0,
cbl     &     0,0,0,0,0,0,0,0,0,0/
      DATA NROOT_DATA_SAC /4,4,4,4,4,4,4,3,3,3,
     &     3,3,1,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0/
      DATA NROOT_DATA_NOAH /4,4,4,4,4,4,4,3,3,3,
     &     3,3,1,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0/
      DATA NROOT_DATA /4,4,4,4,4,4,4,3,3,3,
     &     3,3,1,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0/
      DATA RMAX_DATA /2.5,2.5,2.5,2.5,2.0,2.0,1.8,1.8,1.6,1.5,
     &     1.0,0.0,0.0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0/
      DATA CROOT_DATA /-1.8,-1.8,-1.8,-1.8,-1.7,-1.7,-1.4,-1.4,-1.2,
     &     -1.1,-1.1,-1.1,-1.1,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
     &     -1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0/
      DATA D50_DATA /0.22,0.22,0.22,0.22,0.18,0.18,0.06,0.06,0.08,0.08,
     &     0.08,0.04,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
     &     0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01/

c     DATA RSMTBL /150.0, 150.0, 100.0, 100.0, 125.0,  70.0,
c    &              70.0, 300.0, 400.0,  40.0,  40.0, 400.0,
c    &             200.0,   0.0,   0.0,   0.0,   0.0,   0.0,
c    &               0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
c    &               0.0,   0.0,   0.0,   0.0,   0.0,   0.0/
c test new table
      DATA RSMTBL /300.0, 300.0, 300.0, 175.0, 175.0,  70.0,
     &              70.0, 225.0, 225.0,  35.0,  35.0, 400.0,
     &             200.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     &               0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     &               0.0,   0.0,   0.0,   0.0,   0.0,   0.0/
      DATA RGLTBL / 30.0,  30.0,  30.0,  30.0,  30.0,  65.0,
     &              65.0, 100.0, 100.0, 100.0, 100.0, 100.0,
     &             100.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     &               0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     &               0.0,   0.0,   0.0,   0.0,   0.0,   0.0/
      DATA HSTBL /47.35, 41.69, 47.35, 54.53, 51.93, 54.53,
     &            54.53, 42.00, 42.00, 36.35, 36.35, 42.00,
     &            42.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     &             0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     &             0.00,  0.00,  0.00,  0.00,  0.00,  0.00/
      DATA SNUPX /0.040, 0.040, 0.040, 0.040, 0.040, 0.040,
     &            0.040, 0.020, 0.020, 0.020, 0.020, 0.013,
     &            0.020, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
c      DATA SNUPX /0.080, 0.080, 0.080, 0.080, 0.080, 0.080,
c     &            0.080, 0.040, 0.040, 0.040, 0.040, 0.025,
c     &            0.040, 0.000, 0.000, 0.000, 0.000, 0.000,
c     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
c     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      DATA Z0_DATA /1.089, 2.653, 0.854, 0.826, 0.563, 0.856,
     &              0.856, 0.238, 0.065, 0.035, 0.035, 0.011,
     &              1.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &              0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &              0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      DATA LAI_DATA /3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     &               3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     &               3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &               0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

c     DATA LAI_MAX /6.00, 6.00, 6.00, 5.99, 5.98, 5.70,
c    &              4.63, 5.07, 6.00, 4.00, 5.98, 0.74,
c    &              4.57, 0.00, 0.00, 0.00, 0.00, 0.00,
c    &              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
c    &              0.00, 0.00, 0.00, 0.00, 0.00, 0.00/ 

c helin new table 2/15/2006
      DATA LAI_MAX /6.00, 6.00, 6.00, 5.99, 5.98, 5.70,
     &              3.50, 5.07, 6.00, 2.64, 3.00, 0.74,
     &              4.57, 0.00, 0.00, 0.00, 0.00, 0.00,
     &              0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     &              0.00, 0.00, 0.00, 0.00, 0.00, 0.00/

      DATA LAI_MIN /5.00, 5.00, 1.00, 1.00, 2.88, 3.36,
     &              1.98, 1.39, 0.64, 0.65, 0.78, 0.06,
     &              1.20, 0.00, 0.00, 0.00, 0.00, 0.00,
     &              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
     &              0.00, 0.00, 0.00, 0.00, 0.00, 0.00/ 
C two parameters for root distribution 

      DATA ROOT_A /6.706, 7.344, 7.066, 5.990, 4.453, 8.235,
     &             7.604, 6.326, 7.718,10.740, 5.558, 5.558,
     &             5.558, 0.000, 0.000, 0.000, 0.000, 0.000,
     &             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &             0.000, 0.000, 0.000, 0.000, 0.000, 0.000/ 

      DATA ROOT_B /2.175, 1.303, 1.953, 1.955, 1.631, 1.627,
     &             2.300, 1.567, 1.262, 2.608, 2.614, 2.614,
     &             2.614, 0.000, 0.000, 0.000, 0.000, 0.000,
     &             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &             0.000, 0.000, 0.000, 0.000, 0.000, 0.000/

      DATA EMISS_DATA /1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     &               1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     &               1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     &               1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     &               1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
C  ----------------- Emissivity for SNOW -------------------------------
C      DATA EMISSNOW_DATA /0.90/
       DATA EMISSNOW_DATA /1.00/ 
C ----------------------------------------------------------------------
C CLASS PARAMETER 'SLOPETYP' WAS INCLUDED TO ESTIMATE LINEAR RESERVOIR
C COEFFICIENT 'SLOPE' TO THE BASEFLOW RUNOFF OUT OF THE BOTTOM LAYER.
C LOWEST CLASS (SLOPETYP=0) MEANS HIGHEST SLOPE PARAMETER = 1.
C DEFINITION OF SLOPETYP FROM 'ZOBLER' SLOPE TYPE:
C SLOPE CLASS  PERCENT SLOPE
C 1            0-8
C 2            8-30
C 3            > 30
C 4            0-30
C 5            0-8 & > 30
C 6            8-30 & > 30
C 7            0-8, 8-30, > 30
C 9            GLACIAL ICE
C BLANK        OCEAN/SEA
C ----------------------------------------------------------------------
C NOTE:
C CLASS 9 FROM 'ZOBLER' FILE SHOULD BE REPLACED BY 8 AND 'BLANK' 9
C ----------------------------------------------------------------------
      REAL SLOPE
      REAL SLOPE_DATA(MAX_SLOPETYP)

      DATA SLOPE_DATA /0.1,  0.6, 1.0, 0.35, 0.55, 0.8,
     &                 0.63, 0.0, 0.0, 0.0,  0.0,  0.0,
     &                 0.0 , 0.0, 0.0, 0.0,  0.0,  0.0,
     &                 0.0 , 0.0, 0.0, 0.0,  0.0,  0.0,
     &                 0.0 , 0.0, 0.0, 0.0,  0.0,  0.0/

C ----------------------------------------------------------------------
C SET NAMELIST FILE NAME
C ----------------------------------------------------------------------
      CHARACTER*50 NAMELIST_NAME

C ----------------------------------------------------------------------
C SET UNIVERSAL PARAMETERS (NOT DEPENDENT ON SOIL, VEG, SLOPE TYPE)
C ----------------------------------------------------------------------
      INTEGER I
      INTEGER NSOIL
      INTEGER SLOPETYP
      INTEGER SOILTYP
      INTEGER VEGTYP

      INTEGER BARE
      DATA BARE /12/

      LOGICAL LPARAM
      DATA LPARAM /.TRUE./

      LOGICAL LFIRST
      DATA LFIRST /.TRUE./

C ----------------------------------------------------------------------
C PARAMETER USED TO CALCULATE ROUGHNESS LENGTH OF HEAT.
C ----------------------------------------------------------------------
      REAL CZIL
      REAL CZIL_DATA
C   changed in version 2.6 June 2nd 2003
C      DATA CZIL_DATA /0.2/
c     DATA CZIL_DATA /0.1/
c test 0.05
      DATA CZIL_DATA /0.05/

C ----------------------------------------------------------------------
C PARAMETER USED TO CALUCULATE VEGETATION EFFECT ON SOIL HEAT FLUX.
C ----------------------------------------------------------------------
      REAL SBETA
      REAL SBETA_DATA
      DATA SBETA_DATA /-2.0/

C ----------------------------------------------------------------------
C BARE SOIL EVAPORATION EXPONENT USED IN DEVAP.
C ----------------------------------------------------------------------
      REAL FXEXP
      REAL FXEXP_DATA
      DATA FXEXP_DATA /2.0/

C ----------------------------------------------------------------------
C SOIL HEAT CAPACITY [J M-3 K-1]
C ----------------------------------------------------------------------
      REAL CSOIL
      REAL CSOIL_DATA
C      DATA CSOIL_DATA /1.26E+6/
      DATA CSOIL_DATA /2.00E+6/

C ----------------------------------------------------------------------
C SPECIFY SNOW DISTRIBUTION SHAPE PARAMETER SALP - SHAPE PARAMETER OF
C DISTRIBUTION FUNCTION OF SNOW COVER. FROM ANDERSON'S DATA (HYDRO-17)
C BEST FIT IS WHEN SALP = 2.6
C ----------------------------------------------------------------------
      REAL SALP
      REAL SALP_DATA
C     changed for version 2.6 June 2nd 2003
C      DATA SALP_DATA /2.6/
      DATA SALP_DATA /4.0/

C ----------------------------------------------------------------------
C KDT IS DEFINED BY REFERENCE REFKDT AND DKSAT; REFDK=2.E-6 IS THE SAT.
C DK. VALUE FOR THE SOIL TYPE 2
C ----------------------------------------------------------------------
      REAL REFDK
      REAL REFDK_DATA
      DATA REFDK_DATA /2.0E-6/

      REAL REFKDT
      REAL REFKDT_DATA
      DATA REFKDT_DATA /3.0/

      REAL FRZX
      REAL KDT

C ----------------------------------------------------------------------
C FROZEN GROUND PARAMETER, FRZK, DEFINITION: ICE CONTENT THRESHOLD ABOVE
C WHICH FROZEN SOIL IS IMPERMEABLE REFERENCE VALUE OF THIS PARAMETER FOR
C THE LIGHT CLAY SOIL (TYPE=3) FRZK = 0.15 M.
C ----------------------------------------------------------------------
      REAL FRZK
      REAL FRZK_DATA
      DATA FRZK_DATA /0.15/

      REAL RTDIS1(NSOIL)
      REAL RTDIS(NSOIL)
      REAL SLDPTH(NSOIL)
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C SET TWO CANOPY WATER PARAMETERS.
C ----------------------------------------------------------------------
      REAL CFACTR
      REAL CFACTR_DATA
      DATA CFACTR_DATA /0.5/

      REAL CMCMAX
      REAL CMCMAX_DATA
      DATA CMCMAX_DATA /0.5E-3/

C ----------------------------------------------------------------------
C SET MAX. STOMATAL RESISTANCE.
C ----------------------------------------------------------------------
      REAL RSMAX
      REAL RSMAX_DATA
      DATA RSMAX_DATA /5000.0/

C ----------------------------------------------------------------------
C SET OPTIMUM TRANSPIRATION AIR TEMPERATURE.
C ----------------------------------------------------------------------
      REAL TOPT
      REAL TOPT_DATA
      DATA TOPT_DATA /298.0/

C ----------------------------------------------------------------------
C SPECIFY DEPTH[M] OF LOWER BOUNDARY SOIL TEMPERATURE.
C ----------------------------------------------------------------------
      REAL ZBOT
      REAL ZBOT_DATA
C     changed for version 2.5.2
C      DATA ZBOT_DATA /-3.0/
      DATA ZBOT_DATA /-8.0/

C ----------------------------------------------------------------------
C SET TWO SOIL MOISTURE WILT, SOIL MOISTURE REFERENCE PARAMETERS
C ----------------------------------------------------------------------
      REAL SMLOW
      REAL SMLOW_DATA
      DATA SMLOW_DATA /0.5/

      REAL SMHIGH
      REAL SMHIGH_DATA
C     changed in 2.6 from 3 to 6 on June 2nd 2003
      DATA SMHIGH_DATA /2.0/
      REAL RTACC
      REAL RMAX,CROOT,D50

C ----------------------------------------------------------------------
C NAMELIST DEFINITION:
C ----------------------------------------------------------------------
      NAMELIST /SOIL_VEG/ SLOPE_DATA, RSMTBL, RGLTBL, HSTBL, SNUPX,
     &  BB, DRYSMC, F11, MAXSMC, REFSMC, SATPSI, SATDK, SATDW,
     &  WLTSMC, QTZ, LPARAM, ZBOT_DATA, SALP_DATA, CFACTR_DATA,
     &  CMCMAX_DATA, SBETA_DATA, RSMAX_DATA, TOPT_DATA,
     &  REFDK_DATA, FRZK_DATA, BARE, DEFINED_VEG, DEFINED_SOIL,
     &  DEFINED_SLOPE, FXEXP_DATA, NROOT_DATA, REFKDT_DATA, Z0_DATA,
     &  CZIL_DATA, LAI_DATA, CSOIL_DATA, RMAX_DATA, CROOT_DATA, D50_DATA

C ----------------------------------------------------------------------
C READ NAMELIST FILE TO OVERRIDE DEFAULT PARAMETERS ONLY ONCE.
C NAMELIST_NAME must be 50 characters or less.
C ----------------------------------------------------------------------
      IF (LFIRST) THEN
c        WRITE(*,*) 'READ NAMELIST'
C        OPEN(58, FILE = 'namelist_filename.txt')
C         READ(58,'(A)') NAMELIST_NAME
C         CLOSE(58)
C         WRITE(*,*) 'Namelist Filename is ', NAMELIST_NAME
C         OPEN(59, FILE = NAMELIST_NAME)
C 50      CONTINUE
C         READ(59, SOIL_VEG, END=100)
C         IF (LPARAM) GOTO 50
C 100     CONTINUE
C         CLOSE(59)
c         WRITE(*,NML=SOIL_VEG)
         LFIRST = .FALSE.
         IF (DEFINED_SOIL .GT. MAX_SOILTYP) THEN
            WRITE(*,*) 'Warning: DEFINED_SOIL too large in namelist'
            STOP 222
         ENDIF
         IF (DEFINED_VEG .GT. MAX_VEGTYP) THEN
            WRITE(*,*) 'Warning: DEFINED_VEG too large in namelist'
            STOP 222
         ENDIF
         IF (DEFINED_SLOPE .GT. MAX_SLOPETYP) THEN
            WRITE(*,*) 'Warning: DEFINED_SLOPE too large in namelist'
            STOP 222
         ENDIF
         
         SMLOW = SMLOW_DATA
         SMHIGH = SMHIGH_DATA
         
         DO I = 1,DEFINED_SOIL
            SATDW(I)  = BB(I)*SATDK(I)*(SATPSI(I)/MAXSMC(I))
            F11(I) = ALOG10(SATPSI(I)) + BB(I)*ALOG10(MAXSMC(I)) + 2.0
            if (model_type == 0) then
               REFSMC1 = MAXSMC(I)*(5.79E-9/SATDK(I))
     &              **(1.0/(2.0*BB(I)+3.0))
               REFSMC(I) = REFSMC1 + (MAXSMC(I)-REFSMC1) / SMHIGH
               WLTSMC1 = MAXSMC(I) * (200.0/SATPSI(I))**(-1.0/BB(I))
               WLTSMC(I) = WLTSMC1 - SMLOW * WLTSMC1
               DRYSMC(I) = WLTSMC(I)
            else
               WLTSMC(I) = WLTSMC2(I)
               DRYSMC(I) = WLTSMC2(I)
            endif
            
C     ----------------------------------------------------------------------
C     CURRENT VERSION DRYSMC VALUES THAT EQUATE TO WLTSMC.
C     FUTURE VERSION COULD LET DRYSMC BE INDEPENDENTLY SET VIA NAMELIST.
C     ----------------------------------------------------------------------
         END DO
         
C     ----------------------------------------------------------------------
C     END LFIRST BLOCK
C     ----------------------------------------------------------------------
      ENDIF
      
      IF (SOILTYP .GT. DEFINED_SOIL) THEN
        WRITE(*,*) 'Warning: too many soil types'
        STOP 333
      ENDIF
      IF (VEGTYP .GT. DEFINED_VEG) THEN
        WRITE(*,*) 'Warning: too many veg types'
        STOP 333
      ENDIF
      IF (SLOPETYP .GT. DEFINED_SLOPE) THEN
        WRITE(*,*) 'Warning: too many slope types'
        STOP 333
      ENDIF

C ----------------------------------------------------------------------
C SET-UP UNIVERSAL PARAMETERS (NOT DEPENDENT ON SOILTYP, VEGTYP OR
C SLOPETYP)
C ----------------------------------------------------------------------
      ZBOT = ZBOT_DATA
      SALP = SALP_DATA
      CFACTR = CFACTR_DATA
      CMCMAX = CMCMAX_DATA
      SBETA = SBETA_DATA
      RSMAX = RSMAX_DATA
      TOPT = TOPT_DATA
      REFDK = REFDK_DATA
      FRZK = FRZK_DATA
      FXEXP = FXEXP_DATA
      REFKDT = REFKDT_DATA
      CZIL = CZIL_DATA
      CSOIL = CSOIL_DATA

C ----------------------------------------------------------------------
C  SET-UP SOIL PARAMETERS
C ----------------------------------------------------------------------
      BEXP = BB(SOILTYP)
CYoulong Change ------ multiply by 5.5 according to MOPEX calibration ---
      DKSAT = 5.5*SATDK(SOILTYP)
      DWSAT = SATDW(SOILTYP)
      F1 = F11(SOILTYP)
CY      KDT = REFKDT * DKSAT/REFDK
CYoulong Change ------ from MOPEX calibration ----------------------------
      KDT = 2.59-0.044*DKSAT
      PSISAT = SATPSI(SOILTYP)
      QUARTZ = QTZ(SOILTYP)
      if (model_type == 0) then
         SMCDRY = DRYSMC(SOILTYP)
         SMCMAX = MAXSMC(SOILTYP)
         SMCREF = REFSMC(SOILTYP)
         SMCWLT = WLTSMC(SOILTYP)
      else
         SMCDRY = WLTSMC2(SOILTYP)
         SMCMAX = MAXSMC2(SOILTYP)
         SMCREF = REFSMC2(SOILTYP)
         SMCWLT = WLTSMC2(SOILTYP)
      endif
      FRZFACT = (SMCMAX / SMCREF) * (0.412 / 0.468)
      !write(*,*)'wlt max',smcwlt,smcmax

C ----------------------------------------------------------------------
C TO ADJUST FRZK PARAMETER TO ACTUAL SOIL TYPE: FRZK * FRZFACT
C ----------------------------------------------------------------------
      FRZX = FRZK * FRZFACT

C ----------------------------------------------------------------------
C SET-UP VEGETATION PARAMETERS
C ----------------------------------------------------------------------
      NROOT = NROOT_DATA(VEGTYP)
c      IF (MODEL_TYPE == 1.OR.MODEL_TYPE == 3) THEN
c         NROOT = NROOT_DATA_SAC(VEGTYP)
c      ELSE
c         NROOT = NROOT_DATA_NOAH(VEGTYP)
c      ENDIF
      SNUP = SNUPX(VEGTYP)
      RSMIN = RSMTBL(VEGTYP)
      RGL = RGLTBL(VEGTYP)
      HS = HSTBL(VEGTYP)
      Z0 = Z0_DATA(VEGTYP)
C     LAI = LAI_DATA(VEGTYP)
      IF(SHDMAX-SHDMIN.NE.0)THEN
      LAI = LAI_MIN(VEGTYP)+((SHDFAC-SHDMIN)/(SHDMAX-SHDMIN))*
     &      (LAI_MAX(VEGTYP)-LAI_MIN(VEGTYP))
      ELSE
      LAI=LAI_MAX(VEGTYP)
      ENDIF
      LAI=MAX(LAI_MIN(VEGTYP),LAI)
      RTA = ROOT_A(VEGTYP)
      RTB = ROOT_B(VEGTYP)
c     write(99,*)LAI_MIN(VEGTYP),LAI_MAX(VEGTYP),SHDFAC,SHDMIN,SHDMAX,LAI
      EMISS = EMISS_DATA(VEGTYP)
      EMISSNOW = EMISSNOW_DATA
      IF (VEGTYP .EQ. BARE) SHDFAC = 0.0

      IF (NROOT .GT. NSOIL) THEN
        WRITE(*,*) 'Warning: too many root layers'
        STOP 333
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE ROOT DISTRIBUTION.  PRESENT VERSION ASSUMES UNIFORM
C DISTRIBUTION BASED ON SOIL LAYER DEPTHS.
C ----------------------------------------------------------------------
c     DO I = 1,NROOT
c       RTDIS(I) = -SLDPTH(I)/ZSOIL(NROOT)
c     END DO
        RTDIS1(1)  = 1.0 - 0.5*(EXP(-RTA*(0.5*SLDPTH(1)))
     &               +EXP(-RTB*(0.5*SLDPTH(1))))
      DO I = 2,NROOT
        RTDIS1(I) = 1.0 - 0.5*(EXP(-RTA*(0.5*SLDPTH(I)-
     &              ZSOIL(NROOT-1)))
     &              +EXP(-RTB*(0.5*SLDPTH(I)-ZSOIL(NROOT-1))))
      END DO

        RTDIS(1) = RTDIS1(1)

      DO I = 2,NROOT-1
        RTDIS(I) = RTDIS1(I)-RTDIS1(I-1)
      END DO

        RTACC = 0
      DO I = 1,NROOT-1
        RTACC = RTDIS(I)+RTACC
      ENDDO

        RTDIS(NROOT)=1-RTACC

      DO I = 1, NROOT
       RTDIS(I) = MIN (RTDIS(I), 1.)
       RTDIS(I) = MAX (RTDIS(I), 0.)
      ENDDO
c no root distribution function used
      DO I = 1,NROOT
        RTDIS(I) = -SLDPTH(I)/ZSOIL(NROOT)
cbl        write(*,*)'rtdis',i,rtdis(i)
      END DO

      if (model_type==1) then
         RMAX = RMAX_DATA(VEGTYP)
         CROOT = CROOT_DATA(VEGTYP)
         D50 = D50_DATA(VEGTYP)
         NROOT = 0
         do i =1,nsoil
            rtdis1(i) = RMAX/(1+(-ZSOIL(I)/D50)**CROOT)
cbl         write(*,*)'ZSOIL',I,-ZSOIL(I)
cbl         write(*,*)rmax,croot,d50
cbl         write(*,*)'rtdis1',i,rtdis1(i)
            if (I.NE.1.AND.RMAX.GT.-zsoil(i-1)) then
               NROOT = I
            endif
         enddo
cbl      write(*,*)'nroot',nroot

         rtdis(1) = rtdis1(1)/RMAX 
cbl      write(*,*)'RTDIS',1.0,RTDIS(1)
         do i = 2,nsoil
            rtdis(i) = (RTDIS1(I) - RTDIS1(I-1)) / RMAX
cbl         write(*,*)'RTDIS',i,rtdis(i)
         enddo
      endif
C ----------------------------------------------------------------------
C  SET-UP SLOPE PARAMETER
C ----------------------------------------------------------------------
      SLOPE = SLOPE_DATA(SLOPETYP)

C ----------------------------------------------------------------------
C END SUBROUTINE REDPRM
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE WETBULB (TW,TA,Q2,PA)                                           
C***********************************************************************
C     COMPUTES WET-BULB(TW)FROM DRY-BULB(TA) AND VAPOR PRESSURE(EA)         
C***********************************************************************
      EA = ((Q2/(Q2+0.622))*SFCPRS)/100
      TAC= TA-273.16                                      
      EAS= 2.7489E8*EXP(-4278.63/(TAC+242.792))
      DELTA= (4278.63/((TAC+242.792)*(TAC+242.792)))*EAS
      DO 100 I=1,3
      TWC= DELTA*TAC+6.6E-4*PA*TAC+7.59E-7*PA*TAC*TAC+EA-EAS
      TWC= TWC/(DELTA+6.6E-4*PA+7.59E-7*PA*TAC)
      TAV= (TAC+TWC)*0.5
      EAV= 2.7489E8*EXP(-4278.63/(TAV+242.792))
      DELTA= (4278.63/((TAV+242.792)*(TAV+242.792)))*EAV
  100 CONTINUE
C     CONVERT WET-BULB TO DEGREES KELVIN
      TW= TWC+273.16
      RETURN
      END


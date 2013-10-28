C     POINT MODEL OF THE ENERGY AND MASS BALANCE OF A SNOW COVER.               
C     ....NOTE....THIS MODEL IS NOT DESIGNED FOR SIMULATION OF VERY SHALLOW     
C         INTERMITTENT SNOW COVERS.  A SNOW COVER MUST EXIST AT THE TIME        
C         COMPUTATIONS START.                                                   
C*******************************************************************************
C     INPUT SUMMARY                                                             
C     CARD NO.     FORMAT    REMARKS                                            
C*******************************************************************************
C     1           I5    =0,SOLVE FOR SNOW SURFACE TEMPERATURE.               
C     ITOPT             =1,INPUT SNOW SURFACE TEMPERATURE,                   
C                             COMPUTE OTHER LAYERS. MELT ONLY OCCURS            
C                             WHEN SURFACE TEMPERATURE EQUALS ZERO              
C                             DEGREES CELSIUS.                                  
C     ITS            I5    METHOD OF INPUTING SNOW SURFACE TEMPERATURE          
C                          DATA WHEN IT IS NEEDED AS AN INPUT VARIABLE.         
C                          INPUT VARIABLES MUST NOT CONTAIN MISSING DATA        
C                          =0,INPUT FROM TAPE.                                  
C                          =1,INPUT FROM CARDS                                  
C                          .GT.1 -- USE A SELECTED PATTERN. IN THIS             
C                              CASE THE PROGRAM STARTS ON DAY ONE               
C                              AND RUNS FOR A MAXIMUM OF ONE MONTH. THE         
C                              AVAILABLE PATTERNS ARE.                          
C                              =2,SUDDEN INSTANTANEOUS CHANGE FOLLOWED          
C                                 BY A CONSTANT TEMPERATURE.                    
C                              =3,SIN-FUNCTION VARIATION                        
C      IPUNCH        I5    =1, PUNCH EXCESS LIQUID-WATER (NET MELT PLUS RAIN).  
C                          =0, DO NOT PUNCH                                     
C      TITLE    5X,15A4    GENERAL INFORMATION FOR THE RUN.                     
C*******************************************************************************
C     CARD TWO ONLY NEEDED IF SNOW SURFACE TEMPERATURE IS TO                    
C        BE AN INPUT VARIABLE.                                                  
C        2A          I5    =1,INPUT METEOROLOGICAL VARIABLES                    
C        NOBS                 AND VERIFICATION DATA ARE NOT NEEDED.             
C                             THUS ONLY THE EFFECT OF SURFACE                   
C                             TEMPERATURE (AND SOIL TEMPERATURE) ON             
C                             INTERNAL SNOW COVER TEMPERATURE IS                
C                             TO BE EXAMINED.                                   
C                          =0,METEOROLOGICAL AND VERIFICATION DATA ARE          
C                             NEEDED. WHEN A SELECTED PATTERN IS USED           
C                             OR WHEN SNOW SURFACE TEMPERATURE IS INPUT FROM    
C                             CARDS, THIS CONTROL PARAMETER IS ALWAYS=1.        
C        TSMAX     F5.0    ONLY NEEDED IF A SELECTED PATTERN IS USED.           
C                          VALUE DEPENDS ON THE TYPE OF PATTERN                 
C                          SELECTED.                                            
C                              1.MAGNITUDE OF INSTANTEOUS CHANGE IN             
C                                  DEGREES CELSIUS.                             
C                              2.AMPLITUDE OF SIN-VARIATION IN DEGREES          
C                                  CELSIUS.                                     
C        CYCLE     F5.0    ONLY NEEDED IF SIN-VARIATION IS USED.                
C                          TIME IN DAYS FOR ONE COMPLETE CYCLE.                 
C        LENGTH      I5    ONLY NEEDED FOR SELECTED PATTERN RUNS. LENGTH        
C                             OF THE RUN IN DAYS. 31 IS THE MAXIMUM.            
C                             THE LONGER OF 2 DAYS OR 2 CYCLES IS ADDED AT THE  
C                             BEGINNING OF THE SIN VARIATION PATTERN WITH NO    
C                             OUTPUT.                                           
C        TI        F5.0    ONLY NEEDED FOR SELECTED PATTERN RUNS.  INITIAL      
C                             UNIFORM SNOW COVER TEMPERATURE IN DEGREES CELSIUS.
C*******************************************************************************
C     THIS CARD NEEDED WHENEVER A SELECTED SURFACE TEMPERATURE                  
C        PATTERN IS NOT BEING USED.                                             
C        2 IMO       I5    INITIAL MONTH TO BE USED.                            
C          IDA       I5    INITIAL DAY                                          
C          IYR       I5    INITIAL YEAR (4 DIGITS)                              
C          LMO       I5    LAST MONTH TO BE USED                                
C          LDA       I5    LAST DAY                                             
C          LYR       I5    LAST YEAR (4 DIGITS)                                 
C*******************************************************************************
C     NEEDED PROGRAM INFORMATION.                                               
C        3         F5.0    COMPUTATIONAL TIME INTERVAL IN HOURS.                
C        DELTAT               1,2,3,4,6,8,12,AND 24 HOURS CAN BE USED.          
C        DTOUT     F5.0    OUTPUT TIME INTERVAL IN HOURS. THIS                  
C                              VALUE MUST BE A MULTIPLE OF THE                  
C                              COMPUTATIONAL TIME INTERVAL. MAXIMUM IS          
C                              24 HOURS.                                        
C        THEDA     F5.1    WEIGHTING FACTOR THAT WEIGHTS THE SPACTIAL           
C                              DERIVATIVES WITH RESPECT TO TIME.                
C                              MUST BE IN THE RANGE 0.0 TO 1.0.                 
C                              IF=0.0 THEN AN EXPLICIT FORMULATION.             
C                              IF=1.0 THEN FULLY IMPLICIT.                      
C        TOLER     F5.2    TOLERANCE FOR ITERATIVE SOLUTION IN                  
C                              DEGREES CELSIUS. ALL CORRECTIONS                 
C                              MUST BE LESS THAN THIS VALUE BEFORE              
C                              COMPUTATIONS PROCEED. LIQUID-WATER               
C                              TOLERANCE IS COMPUTED FROM THIS VALUE            
C                              FOR EACH LAYER.                                  
C        IGRAD       I5    =1,USE TEMPERATURE GRADIENT RATIO AND PREVIOUS       
C                             CHANGE IN HEAT STORAGE TO TRY TO IMPROVE THE      
C                             FIRST GUESS FOR THE UPPER THREE LAYERS OVER THE   
C                             REGULAR FIRST GUESS METHOD.                       
C                          =0, USE REGULAR FIRST GUESS METHOD. (BASED ON A      
C                             LINEAR PROJECTION OF THE PREVIOUS CHANGE FOR      
C                             INTERNAL LAYERS.  THE SURFACE TEMPERATURE,        
C                             COMPUTED WITH HEAT STORAGE=0.0, IS USED FOR THE   
C                             SURFACE LAYER.)                                   
C        GRMAX     F5.0    MAXIMUM GRADIENT RATIO TO BE ALLOWED.                
C                              2.0 IS RECOMMENDED.                              
C        DTONE     F5.0    PROJECTED TEMPERATURE CHANGE(DEGREES CELSIUS) IN     
C                             SURFACE LAYERS BELOW WHICH THE BASIC COMPUTATIONAL
C                             TIME INTERVAL WILL BE USED.  IF PROJECTED CHANGE  
C                             EXCEEDS THIS VALUE THE PERIOD WILL BE SUB-DIVIDED.
C                             FORTRAN IDENTIFIER IS(DTONE).                     
C                  *** NOTE *** SUB-DIVIDING OF THE TIME INTERVAL CANNOT BE     
C                            USED WHEN A THEORETICALLY BASED WIND FUNCTION IS   
C                            USED.  THIS IS BECAUSE F(U) VARIES THUS ENERGY     
C                            BALANCE COMPONENTS CANNOT BE COMPUTED IN SUBROUTINE
C                            CHECK.  IN THIS CASE DTONE IS SET TO 99.0.         
C         TEXP     F5.1    EXPONENT IN SUB-DIVIDING EQUATION. FORTRAN IDENTIFIER
C                             IS (TEXP).  THE EQUATION IS.                      
C                                  NINC=((CHANGE-DTONE)**TEXP)+2                
C                             WHERE NINC IS THE NUMBER OF TIME INCREMENTS TO    
C                                USE FOR THAT PERIOD.                           
C         IQAE       I5    =1, USE ESTIMATED ATMOSPHERIC LONGWAVE RADIATION     
C                          =0, USE OBSERVED ATMOSPHERIC LONGWAVE.               
C                                  (ONLY NEEDED IF METEOROLOGICAL VARIABLES     
C                                       ARE BEING USED.)                        
C         THICK    F5.1    DESIRED THICKNESS OF THE SURFACE LAYERS OF THE       
C                             SNOW COVER IN CENTIMETERS.(USED DOWN TO 30 CM.)   
C                             THIS IS AN INDEX. ACTUAL THICKNESS WILL VARY      
C                             FROM 55 TO 155 PERCENT OF THE DESIRED VALUE.      
C                             (NOT NEEDED WITH SELECTED PATTERN RUNS OR         
C                               WHEN SURFACE TEMPERATURE IS INPUT FROM CARDS.)  
C         CTHICK   F5.2    COEFFICIENT (CTHICK) THAT INCREASES THE DESIRED      
C                             THICKNESS WITH DEPTH AFTER THE TOP 30 CM.         
C                             EQUATION IS.                                      
C                                  DESIRED THICKNESS=SURFACE LAYER THICKNESS    
C                                          +CTHICK*(Z-30.0)                     
C                             WHERE. Z=DEPTH BELOW SURFACE.                     
C                          RECOMMENDED VALUE FOR CTHICK IS 0.05                 
C                             (NOT NEEDED WITH SELECTED PATTERN RUNS OR         
C                               WHEN SURFACE TEMPERATURE IS INPUT FROM CARDS.)  
C         ADJQA    F5.2    ATMOSPHERIC LONGWAVE RADIATION ADJUSTMENT FACTOR.    
C         IFU      I5      =1, USE THEORETICAL WIND FUNCTION WITH               
C                                                STABILITY ADJUSTMENT.          
C                          =0, USE EMPIRICAL WIND FUNCTION.                     
C         HEIGHT   F5.1    HEIGHT OF WIND,TEMP., AND VAPOR PRESSURE             
C                                    MEASUREMENTS (METERS)                      
C*******************************************************************************
C     PARAMETER VALUES.                                                         
C        4         F5.2    DIFFUSION COEFFICIENT FOR WATER VAPOR                
C        D_O                IN SNOW AT ZERO DEGREES CELSIUS                      
C                          AND 1000 MB. PRESSURE IN CM2/SEC.                    
C        PA        F5.0    AVERAGE STATION PRESSURE IN MILLIBARS.               
C        X         F5.1    EXPONENT IN DIFFUSION COEFFIENT EQUATION.            
C        DTG       F5.0    DEPTH OF SOIL TEMPERATURE MEASUREMENT IN CM.         
C        TCG       F5.3    THERMAL CONDUCTIVITY OF THE SOIL IN                  
C                              CAL/CM/SEC/DEG.C.                                
C        DCG       F5.2    DIFFUSION COEFFICIENT FOR WATER VAPOR                
C                              IN THE SOIL IN CM2/SEC.                          
C        SFC       F5.2    PRECIPITATION MULTIPLICATION FACTOR DURING SNOW.     
C                             ADJUSTS FOR THE DIFFERENCE BETWEEN THE INCREASE   
C                             IN WATER-EQUIVALENT AT THE SITE AND THE           
C                             PRECIPITATION MEASURED IN THE GAGE.               
C        RCF       F5.2    PRECIPITATION MULTIPLICATION FACTOR DURING RAIN.     
C        PLWHC     F5.2    PERCENT(DECIMAL) LIQUID-WATER HOLDING CAPACITY OF    
C                             THE SNOW. (THIS IS THE BASIC VALUE.)              
C        PLWMAX    F5.2    PERCENT(DECIMAL) LIQUID-WATER-HOLDING CAPACITY       
C                             AT ZERO DENSITY.                                  
C        PLWDEN    F5.2    DENSITY ABOVE WHICH THE BASIC LIQUID-WATER-HOLDING   
C                             CAPACITY IS USED.  THE CAPACITY VARIES LINEARLY   
C                             BETWEEN ZERO AND THIS VALUE.                      
C        FUCOEF   F5.5    WIND FUNCTION COEFFICIENT (MM/MB/KM)                 
C                                  EMPIRICAL F(U) ONLY.                         
C        COEFKE    F5.4    EFFECTIVE THERMAL CONDUCTIVITY COEFFICIENT IN FORM-- 
C                                 KE=0.00005+COEFFICIENT*(DENSITY**2)           
C        CV        F5.2    EXTINCTION COEFFICIENT PARAMETER (CV)                
C                            EXTINCTION COEFF.=CV*DENSITY*SQRT(1/GRAIN SIZE)    
C        RICRIT    F5.2    CRITICAL RICHARDSON NUMBER.  (RI VALUE WHEN          
C                              STABILITY CORRECTION GOES TO ZERO.)              
C        ZO        F5.2    ROUGHNESS HEIGHT (Z0) IN CENTIMETERS.  (FOR USE IN   
C                                 DETERMINING THE BULK TRANSFER COEFFICIENT     
C                                 UNDER NEUTRAL CONDITIONS.)                    
C     ***NOTE***IF THE DIFFUSION COEFFICIENT FOR SNOW IS ZERO THEN VAPOR        
C            TRANSFER WITHIN THE SNOW COVER IS NEGLECTED.                       
C*******************************************************************************
C     THIS CARD CONTAINS METAMORPHISM AND WATER TRANSMISSION PARAMETERS.        
C        4A                COMPACTION PARAMETERS.                               
C                  F5.4       C1                                                
C                  F5.1       C2                                                
C                          DESTRUCTIVE METAMORPHISM PARAMETERS.                 
C                  F5.3       C3                                                
C                  F5.2       C4                                                
C      DMETA       F5.2       DENSITY LIMIT                                     
C                          LIQUID-WATER METAMORPHISM                            
C                  F5.1       C5                                                
C                          LIQUID-WATER TRANSMISSION PARAMETERS                 
C                  F5.1       CW1 (MAXIMUM VALUE IS 10.0)                       
C                 F10.5       CW2                                               
C                 2F5.1       CW3,CW4                                           
C                          GRAIN SIZE PARAMETERS FOR DETERMINING THE            
C                               EXTINCTION COEFFICIENT.                         
C                               GRAIN SIZE=G1+G2*(DEN**2)+G3*(DEN**4)           
C                 3F5.0    G1,G2,G3                                             
C*******************************************************************************
C     CARDS 5 AND 6 DEFINE INITIAL SNOW COVER CONDITIONS.--SUBROUTINE BEGIN     
C*******************************************************************************
C   SUBROUTINE BEGIN
C        5    N      I5    TOTAL NUMBER OF LAYERS IN THE INITIAL                
C             MM             SNOW COVER. 100 IS MAXIMUM.                      
C                    I5    UNITS OF LIQUID-WATER CONTENT ON CARD 6.             
C                            =1, MILLIMETERS.                                   
C                            =0, PERCENT (DECIMAL)                              
C*******************************************************************************
C        6  NN       I5    LAYER NUMBER                                         
C           THICK 5.0     THICKNESS OF THE LAYER IN CENTIMETERS.               
C           DEN   F5.2     DENSITY OF THE LAYER--DECIMAL. (ICE PLUS LIQUID)     
C           TEMP  F5.0     MEAN TEMPERATURE OF THE LAYER IN DEGREES             
C                              CELSIUS.                                         
C           WATER F5.2     LIQUID-WATER CONTENT OF THE LAYER.                   
C     REPEAT CARD 6 AS NEEDED. A CARD IS NEEDED WHENEVER INITIAL                
C        CONDITIONS CHANGE. THAT IS, EACH CARD DEFINES INITIAL                  
C        VALUES FOR LAYER NP+1 THROUGH THE INDICATED LAYER,WHERE                
C        NP IS THE LAYER NUMBER ON THE PREVIOUS CARD. THE LAYER                 
C        NUMBER ON THE LAST CARD 6 MUST MATCH WITH THE TOTAL                    
C        NUMBER OF LAYERS.                                                      
C*******************************************************************************
C...NOTE...REPEAT CARDS 7-9 FOR EACH MONTH THAT IS INCLUDED IN THE RUN.  DATAIN
C*******************************************************************************
C   SUBROUTINE DATAIN
C     SOIL TEMPERATURE DATA                                                     
C        7 MO        I2    MONTH NUMBER                                         
C          JYR    1X,I2    YEAR (2 DIGITS)                                      
C                    THE REMAINING THREE VALUES ARE REPEATED AS A               
C                    GROUP FIVE TIMES.                                          
C          ID(I)     I5    DAY NUMBER                                           
C          IH(I)     I5    HOUR                                                 
C          GT(I)    F5.0   SOIL TEMPERATURE IN DEGREES FAHRENHEIT               
C     SOIL TEMPERATURE VALUES ARE NEEDED FOR THE LAST HOUR                      
C        OF THE MONTH AND ALL OTHER HOURS WHEN CHANGES OCCUR.                   
C        THUS REPEAT CARD 7 AS NECESSARY.                                       
C*******************************************************************************
C   SUBROUTINE DATAIN
C     CARD 8 IS ONLY NEEDED IF SNOW SURFACE TEMPERATURES                        
C        ARE TO BE INPUT FROM CARDS.                                            
C        8 MO        I2    MONTH NUMBER                                         
C          IDAY   1X,I2    DAY NUMBER                                           
C          JYR    1X,I2    YEAR (2 DIGITS)                                      
C          ST    24F3.0    SNOW SURFACE TEMPERATURE FOR EACH HOUR               
C                              IN DEGREES FAHRENHEIT.                           
C     REPEAT CARD 8 FOR EACH DAY OF THE MONTH.                                  
C*******************************************************************************
C   SUBROUTINE DATAIN
C     DENSITY OF NEW SNOW -- CARD 9 ONLY NEEDED IF METEOROLOGICAL               
C        AND VERIFICATION DATA ARE NEEDED. ONLY INPUT NEW SNOW                  
C        DENSITY VALUES WHEN THE NORMAL COMPUTATIONS OF THE                     
C        DENSITY OF SNOWFALL ARE TO BE OVER-RIDDEN. THERE MUST                  
C        BE AT LEAST ONE CARD 9 PER MONTH WHENEVER METEOROLOGICAL               
C        AND VERIFICATION DATA ARE NEEDED.                                      
C        9 MO        I2    MONTH NUMBER                                         
C          JYR       I2    YEAR (2 DIGITS)                                      
C          IEND      I1    =1, THIS IS THE LAST NEW SNOW DENSITY CARD.          
C                          =0, THIS IS NOT THE LAST CARD FOR THE MONTH.         
C                   THE NEXT FOUR VALUES ARE REPEATED AS A GROUP FIVE TIMES.    
C          ID(I)    I4     DAY NUMBER OF FIRST HOUR TO BE CHANGED.              
C          IH(I)    I3     FIRST HOUR                                           
C          NH(I)  F5.2     NEW SNOW DENSITY -- DECIMAL                          
C          SD(I)    I3     NUMBER OF HOURS TO BE CHANGED.                       
C        THIS CAN BE USED TO CHANGE THE FORM OF THE PRECIPITATION ALSO.         
C                  =1.0, INDICATES RAIN                                         
C                  (.GT.0.0).AND.(.LT.0.90) INDICATES SNOW.                     
C*******************************************************************************
C*******************************************************************************
C
C Other variables
C     PNS     Precipitation flag (Not-Snow = 1; for snow 0<PNS<1)
C     SUMPX   Sum of precipitation (cm) over snow period -- reset after 
C                ablation is complete
C     WT(N)   Layer liquid water (cm?) at current time step
C     WTDT(N) Layer liquid water at next timestep
C     TT(N)   Layer temperature (Kelvin)
C     TTDT(N) Layer temperature at next timestep
C     DNS     Depth of new snow (NWSNOW.f)
C     TA      Air temperature (Kelvin)
C     VAPOUR  Vapour exchange b/w pack and air (down is positive)
C
C***********************************************************************
      SUBROUTINE PEMB (ZO,MODEL_TYPE,ITOPT,ITS,IPUNCH,IGRAD,IQAE,IFU,NN,
     &     DELTAT,DTOUT,THEDA,TOLER,GRMAX,DTONE,TEXP,THICK,CTHICK,
     &     ADJQA,HEIGHT,D_O,PA,X,DTG,TCG,DCG,SFC,RCF,PLWHC,PLWMAX,
     &     PLWDEN,FUCOEF,COEFKE,CV,RICRIT,C1,C2,C3,C4,DMETA,C5,
     &     CW1,CW2,CW3,CW4,G1,G2,G3,MM,DEN,TEMP,STC,NEWPACK,THICK0,
     &     WATER0,SFCTMP,Q2,SFCPRS,SFCSPD,SOLDN,LWDN,
     &     SNOALB,PRCP,SNOWNG,
C State and output variables
     &     N,D,P,RTT,RTTDT,WT,WTDT,TT,TTDT,SCOUT,VAPOUR,QG)

      DIMENSION D(100),P(100),TT(100),WT(100),NOKNOW(100),ND(12),               
     & TTDT(100),WTDT(100),RTT(100),RTTDT(100),TITLE(15)
C      COMMON/SPCASE/ITS                                                         
C      COMMON/WPM/CW1,CW2,CW3,CW4                                               
C      COMMON/WSPEED/IFU,UMSEC,FUCOEF,HEIGHT                                     
C      COMMON/EFFTC/COEFKE                                                       
C      COMMON/EXTINC/CV,G1,G2,G3                                                 
C      COMMON/CRITRI/RICRIT,ZO                                                   
C      DATA ND/31,28,31,30,31,30,31,31,30,31,30,31/                              
      NOBS=0                                                                    
C     DEFINE THE TYPE OF RUN  
C     ITOPT,ITS = SNOW SURF TEMP SOLVE FLAGS
C     IPUNCH = LIQUID WATER FLAG
C      READ 908,ITOPT,ITS,IPUNCH,TITLE                                           
C  908 FORMAT (3I5, 5X,15A4)                                                     
      IF (ITOPT.NE.0) IPUNCH=0                                                  
      IF (ITS.LT.1)ITS=0                                                        
      IF (ITOPT.NE.1) GO TO 106     
C     ONLY NEEDED IF SURFACE TEMP TO BE INPUT
C      READ 901,NOBS,TSMAX,CYCLE,LENGTH,TI                                       
C  901 FORMAT (I5,2F5.0,I5,F5.0)                                                 
      IF (TI.GT.0.0) STOP 10                                                    
      IF (ITS.GT.0)NOBS=1                                                       
  106 IF (ITOPT.NE.1) ITS=0                                                     
      IF (ITOPT.NE.1) NOBS=0                                                    
C     INITIAL MONTH,DAY,YEAR; LAST DAY,MONTH,YEAR
C      IF (ITS.LT.2)READ 900,IMO,IDA,IYR,LMO,LDA,LYR                             
C  900 FORMAT (6I5)           
C     DELTAT(HOURS),DTOUT(HOURS),THEDA-IMPLICIT/EXPLICIT FACTOR
C     TOLER(DEG.C),IGRAD
C      READ 902,DELTAT,DTOUT,THEDA,TOLER,IGRAD,GRMAX,DTONE,TEXP,IQAE,            
C     1THICK,CTHICK,ADJQA,IFU,HEIGHT                                             
C  902 FORMAT (2F5.0,F5.1,F5.2,I5,2F5.0,F5.1,I5,F5.1,2F5.2,I5,F5.1)              
      IF (IFU.EQ.1) DTONE=99.0                                                  
C      IDT=DELTAT+0.001                                                          
C      IF ((IDT*(24/IDT)).NE.24)STOP 11                                          
C      IRT=(DTOUT/DELTAT)+0.001                                                  
C      IOUT=IDT*IRT                                                              
C      IF (IOUT.GT.24) STOP 12                                                   
      IF ((THEDA.LT.0.0).OR.(THEDA.GT.1.0))STOP 13                              
C      READ 903,DO,PA,X,DTG,TCG,DCG,SCF,RCF,PLWHC,PLWMAX,PLWDEN,FUCOEF,          
C     1COEFKE,CV,RICRIT,ZO                                                       
C  903 FORMAT (F5.2,F5.0,F5.1,F5.0,F5.3,6F5.2,F5.5,F5.4,3F5.2)                   
      IF (IFU.LT.1) GO TO 107                                                   
      CWN=0.16/((ALOG(HEIGHT*100.0/ZO))**2)                                     
      AIR=0.00035*PA/273.16                                                     
      FUCOEF=(AIR*0.622E6*CWN)/PA                                               
  107 CONTINUE                                                                  
C      READ 917,C1,C2,C3,C4,DMETA,C5,CW1,CW2,CW3,CW4,G1,G2,G3                    
C  917 FORMAT (F5.4,F5.1,F5.3,F5.2,F5.2,F5.1,F5.1,F10.5,2F5.1,3F5.0)             
C Time stuff and input snow temp stuff
C      IF (C5.LT.1.0) C5=1.0                                                     
C      IF (ITS.EQ.0) GO TO 100                                                   
C      C1=0.0                                                                    
C      C3=0.0                                                                    
C      C5=0.0                                                                    
C      IF (ITS.LT.2) GO TO 100                                                   
C      IF (LENGTH.GT.31) LENGTH=31                                               
C      IMO=1                                                                     
C      IDA=1                                                                     
C      IYR=1                                                                     
C      LMO=1                                                                     
C      LDA=LENGTH                                                                
C      ICYCLE=0                                                                  
C      IF (ITS.NE.3) GO TO 104                                                   
C      ICYCLE=2.0*CYCLE+0.1                                                      
C      IF (ICYCLE.LT.2) ICYCLE=2                                                 
C      LDA=LDA+ICYCLE                                                            
C      IF (LDA.GT.31) LDA=31                                                     
C      ICYCLE=ICYCLE*24                                                          
C  104 LYR=1                                                                
  102 CONTINUE                                                                  
C*******************************************************************************
C     INITIALIZE THE SNOW COVER: for a new pack, compute water equivalent 
      WRITE(*,*)'ABOVE BEGIN, NEWPACK',NEWPACK
      IF (NEWPACK == 1) THEN
         CALL BEGIN(N,D,P,TT,WT,NOKNOW,ITS,TI,WEI,THICK0,DEN,TEMP,WATER0)      
         SURDEN=P(1)
         SUMPX=0.0                               
      ENDIF                                                          
C      MONTH=IMO                                                                 
C      IYEAR=IYR                                                                 
C      IDAY=IDA                                                                  
C  101 IHOUR=0                                                                   
C*******************************************************************************
C     BEGIN MONTHLY LOOP.                                                       
C      LEAPYR=0                                                                  
C      IF (((IYEAR/4)*4).EQ.IYEAR)LEAPYR=1                                       
C      LAST=ND(MONTH)                                                            
C      IF ((MONTH.EQ.2).AND.(LEAPYR.EQ.1))LAST=29                                
C     INPUT HOURLY DATA FOR THE MONTH.   
C  LAST:#DAYS IN MONTH      
C          
C Subroutine to input soil temperature and forcings (originally from tape)
C and perform unit conversions. Replace monthly (31) and hourly (744) 
C arrays for each quantity with a single value from parent model.
C
C      CALL DATAIN(LAST,ITOPT,ITS,NOBS,TSMAX,TT(1),CYCLE,IMO,IYR,MONTH,          
C     1     IYEAR,IQAE)
C
C Temperatures in Kelvin
      TG = STC
      TAD = SFCTMP
C Vapor pressure in millibars
      EAD = ((Q2/(Q2+0.622))*SFCPRS)/100
C Wind speed m/s
      UA = SFCSPD
C Solar radiation incoming, reflected, longwave incoming (langleys cal/cm2)
C Conversion: W/m**-2 = (DT/41840)*ly 
      QID = SOLDN * (DELTAT/41840)
      QRD = QID * SNOALB
      QAD = LWDN * (DELTAT/41840)
C Snow cover outflow (mm)
      POD = 9999.0
C Precipitation (mm) and type; if snow, compute density later (denns = 2)
      PXD = PRCP * DELTAT
      WRITE(*,*)'PXD TG TAD',PXD,TG,TAD
      WRITE(*,*)'EAD UA QID',EAD,UA,QID
      IF (SNOWNG) THEN
         DENNS = 2.0
      ELSE
         DENNS = 1.0
      ENDIF        
C*******************************************************************************
C     BEGIN DELTA TIME COMPUTATION LOOP.                                        
C  105 IHOUR=IHOUR+IDT                                                           
C      IF (IHOUR.LE.24) GO TO 110                                                
C      IDAY=IDAY+1                                                               
C      IHOUR=IHOUR-24                                                            
C     GET VALUES FOR INPUT VARIABLES AND VERIFICATION DATA FOR
C        THE TIME PERIOD.
      WRITE(*,*)'PEMB:ABOVE OBTAIN N layers',N
  110 CALL OBTAIN(DELTAT,IDAY,IHOUR,PNS,QI,QR,QA,FU,TA,EA,
     1PX,TGT,TGTDT,TT(1),TTDT(1),ITOPT,NOBS,TO,PO,WESC,WEPW,DEPTH,
     2STAKE,PA,TSNOW,SCF,RCF,NOKNOW(1),ADJQA,
     3TG,DENNS,TAD,EAD,UA,QID,QRD,QAD,POD,PXD,NEWPACK
     4IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO)
      IF(NOBS.NE.1) CALL NWSNOW(N,TT,WT,TTDT,WTDT,NOKNOW,D,P,RTT,RTTDT,
     1PX,PNS,TSNOW,SUMPX,THICK,ITOPT)
      IF (N.EQ.0) GO TO 140
C     GET A FIRST GUESS FOR TEMPERATURE AND LIQUID-WATER
C        FOR EACH LAYER FOR TIME T+DT.
      WRITE(*,*)'PEMB:ABOVE GUESS WT TT',WT(1),TT(1)
      CALL GUESS(N,TT,TTDT,WT,WTDT,NOKNOW,D,P,ITOPT,QI,QR,QA,
     1FU,TA,EA,PX,D_O,PA,X,THEDA,TOLER,DELTAT,MONTH,IDAY,IYEAR,
     2IHOUR,RTT,RTTDT,IGRAD,GRMAX,DTONE,TEXP,NINC,ITS,TSMAX,EXTINC,
     3CV,G1,G2,G3,QG,IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE)
      WRITE(*,*)'PEMB:BELOW GUESS WT TT',WT(1),TT(1)
C     COMPUTE TEMPERATURE AND LIQUID-WATER FOR EACH LAYER
C        FOR TIME T+DT.
      IF (N.GT.1) GO TO 111
C     SNOW COVER CONSISTS OF A SINGLE LAYER
      IGUESS=0
      WRITE(*,*)'PEMB:ABOVE SURFAC WT TT',WT(1),TT(1)
      CALL SURFAC(N,TT(1),TTDT(1),WT(1),WTDT(1),D,P,ITOPT,QI,QR,
     1QA,FU,TA,EA,PX,D_O,PA,X,DTG,TGT,TGTDT,TCG,DCG,THEDA,
     2NOKNOW(1),TOLER,DELTAT,ITER,MONTH,IDAY,IYEAR,IHOUR,IGUESS,QG,
     3EXTINC,CV,G1,G2,G3,IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE,ITS)
      WRITE(*,*)'PEMB:BELOW SURFAC WT TT',WT(1),TT(1)
      GO TO 112
C     SNOW COVER CONSISTS OF MORE THAN ONE LAYER.
      WRITE(*,*)'PEMB:ABOVE SNOWTW'
  111 CALL SNOWTW(N,TT,TTDT,WT,WTDT,D,P,ITOPT,QI,QR,QA,FU,
     1TA,EA,PX,D_O,PA,X,DTG,TGT,TGTDT,TCG,DCG,THEDA,NOKNOW,TOLER,
     2DELTAT,ITER,MONTH,IDAY,IYEAR,IHOUR,NINC,QG,ITS,EXTINC,CV,G1,G2,G3,
     3IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE)
 112  WRITE(*,*)'PEMB:ABOVE CHEK'
      IF (ITS.EQ.0)  CALL CHEK(N,D,P,TT,TTDT,WT,WTDT,NOKNOW,SCOUT,
     1THICK,CTHICK,WE,TDEPTH,QI,QR,QA,DELTAT,PA,FU,TA,EA,PX,TGT,TGTDT,
     2     TCG,DCG,DTG,D_O,X,ITOPT,ITS,COEFKE)
      WRITE(*,*)'PEMB:BELOW CHEK N',N
      WRITE(*,*)'PEMB: GO TO 130 BELOW CHEK'
      IF (N.LT.0) GO TO 130
C     RETAIN TEMPERATURE CHANGE DUE TO HEAT TRANSFER FOR EACH LAYER.
      WRITE(*,*)'PEMB:ABOVE RETAIN'
      CALL RETAIN(N,TT,TTDT,RTT,RTTDT)
      IF ((D_O.GT.0.0).OR.(ITS.EQ.0)) CALL VAPOR(N,P,D_O,PA,X,TT,TTDT,D,
     1VAPOUR,DTG,DCG,TGT,TGTDT,EA,FU,DELTAT,WTDT(1),NOBS,WTDT(N),
     2NOKNOW(N),SOILVT)
      IF (NOBS.EQ.1) GO TO 113
  130 WRITE(*,*)'PEMB: AT 130'
      WRITE(*,*)'DELTAT',DELTAT,'N',N,'P',P,'D',D,'WTDT',WTDT
      WRITE(*,*)'TTDT',TTDT,'NOKNOW',NOKNOW,'PLWHC',PLWHC
      WRITE(*,*)'PLWMAX',PLWMAX,'PLWDEN',LWDEN,'PX',PX
      WRITE(*,*)'SCOUT',SCOUT,'MONTH',MONTH,'IDAY',IDAY,'IYEAR',IYEAR
      WRITE(*,*)'IHOUR',IHOUR,'IPUNCH',IPUNCH,'CW1',CW1,'CW2',CW2
      WRITE(*,*)'CW3',CW3,'CW4',CW4
C May need to make WLAG a global array which describes the amount of 
C liquid water that has been lagged
      CALL WATER(DELTAT,N,P,D,WTDT,TTDT,NOKNOW,PLWHC,PLWMAX,PLWDEN,PX,
     &     SCOUT,MONTH,IDAY,IYEAR,IHOUR,IPUNCH,CW1,CW2,CW3,CW4)
      WRITE(*,*)'PEMB:BELOW WATER'
      IF (N.LT.0) GO TO 140
  113 CALL META(N,P,D,WTDT,TTDT,PLWHC,C1,C2,C3,C4,C5,DELTAT,DMETA)
C     CHECK IF PRINTER OUTPUT IS WANTED FOR THIS TIME PERIOD.                   
C      MHR=(IDAY-1)*24+IHOUR                                                     
C      IF (((MHR/IOUT)*IOUT).NE.MHR) GO TO 125                                   
C     A PRINT IS WANTED FOR THIS TIME PERIOD.                                   
      IF (ITS.LT.2) GO TO 120                                                   
C      IF (MHR.LE.ICYCLE) GO TO 150                                              
C     OUTPUT FOR SELECTED SNOW SURFACE TEMPERATURE PATTERNS.                    
C      CALL PATERN (N,TTDT,WTDT,ITS,TSMAX,CYCLE,D,P,IDAY,                        
C     1IHOUR,ITER,D_O,TI,SURDEN)                                                  
      GO TO 150                                                                 
C     OUTPUT FOR OTHER CASES                
C Remove standard output based on time
 120  CONTINUE
C  120 CALL SNOWOT(N,TTDT,WTDT,D,P,MONTH,IDAY,IYEAR,IHOUR,ITER,WE,TDEPTH)        
  125 IF (NOBS.EQ.1) GO TO 150                                                  
C Remove verification check done in STATDA
  140 CONTINUE  
C  140 CALL STATDA(TT(1),TTDT(1),TO,SCOUT,PO,WE,WESC,WEPW,VAPOUR,TDEPTH,         
C     1DEPTH,STAKE,MONTH,IDAY,IYEAR,IHOUR,DELTAT,N,SOILVT,TA)                    
C     END OF DELTA TIME LOOP, CHECK FOR END OF RUN AND MONTH.                   
C*******************************************************************************
  150 CONTINUE
C Remove time controls  
C  150 IF (IHOUR.LT.24) GO TO 105                                                
C      IF (MONTH.NE.LMO) GO TO 160                                               
C      IF ((IYEAR.EQ.LYR).AND.(IDAY.EQ.LDA)) GO TO 200                           
C      GO TO 105                                                                 
C  160 IF (IDAY.LT.LAST) GO TO 105                                               
C      MONTH=MONTH+1                                                             
C      IDAY=1                                                                    
C      IF (MONTH.LE.12) GO TO 101                                                
C      MONTH=1                                                                   
C      IYEAR=IYEAR+1                                                             
C      GO TO 101                                                                 
C     END OF THE RUN                                                            
C*******************************************************************************
C  200 IF (NOBS.NE.1) CALL FINAL(IMO,IDA,IYR,LMO,LDA,LYR,WEI,SUMPX,WE)           
      RETURN                                                                   
      END                                                                       

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
C        DO                IN SNOW AT ZERO DEGREES CELSIUS                      
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
      DIMENSION D(100),P(100),TT(100),WT(100),NOKNOW(100),ND(12),               
     1TTDT(100),WTDT(100),RTT(100),RTTDT(100),TITLE(15)                         
      COMMON/SPCASE/ITS                                                         
      COMMON/WPM/CW1,CW2,CW3,CW4                                                
      COMMON/WSPEED/IFU,UMSEC,FUCOEF,HEIGHT                                     
      COMMON/EFFTC/COEFKE                                                       
      COMMON/EXTINC/CV,G1,G2,G3                                                 
      COMMON/CRITRI/RICRIT,ZO                                                   
      DATA ND/31,28,31,30,31,30,31,31,30,31,30,31/                              
      NOBS=0                                                                    
C     DEFINE THE TYPE OF RUN  
C     ITOPT,ITS = SNOW SURF TEMP SOLVE FLAGS
C     IPUNCH = LIQUID WATER FLAG
      READ 908,ITOPT,ITS,IPUNCH,TITLE                                           
  908 FORMAT (3I5, 5X,15A4)                                                     
      IF (ITOPT.NE.0) IPUNCH=0                                                  
      IF (ITS.LT.1)ITS=0                                                        
      IF (ITOPT.NE.1) GO TO 106     
C     ONLY NEEDED IF SURFACE TEMP TO BE INPUT
      READ 901,NOBS,TSMAX,CYCLE,LENGTH,TI                                       
  901 FORMAT (I5,2F5.0,I5,F5.0)                                                 
      IF (TI.GT.0.0) STOP 10                                                    
      IF (ITS.GT.0)NOBS=1                                                       
  106 IF (ITOPT.NE.1) ITS=0                                                     
      IF (ITOPT.NE.1) NOBS=0                                                    
C     INITIAL MONTH,DAY,YEAR; LAST DAY,MONTH,YEAR
      IF (ITS.LT.2)READ 900,IMO,IDA,IYR,LMO,LDA,LYR                             
  900 FORMAT (6I5)           
C     DELTAT(HOURS),DTOUT(HOURS),THEDA-IMPLICIT/EXPLICIT FACTOR
C     TOLER(DEG.C),IGRAD
      READ 902,DELTAT,DTOUT,THEDA,TOLER,IGRAD,GRMAX,DTONE,TEXP,IQAE,            
     1THICK,CTHICK,ADJQA,IFU,HEIGHT                                             
  902 FORMAT (2F5.0,F5.1,F5.2,I5,2F5.0,F5.1,I5,F5.1,2F5.2,I5,F5.1)              
      IF (IFU.EQ.1) DTONE=99.0                                                  
      IDT=DELTAT+0.001                                                          
      IF ((IDT*(24/IDT)).NE.24)STOP 11                                          
      IRT=(DTOUT/DELTAT)+0.001                                                  
      IOUT=IDT*IRT                                                              
      IF (IOUT.GT.24) STOP 12                                                   
      IF ((THEDA.LT.0.0).OR.(THEDA.GT.1.0))STOP 13                              
      READ 903,DO,PA,X,DTG,TCG,DCG,SCF,RCF,PLWHC,PLWMAX,PLWDEN,FUCOEF,          
     1COEFKE,CV,RICRIT,ZO                                                       
  903 FORMAT (F5.2,F5.0,F5.1,F5.0,F5.3,6F5.2,F5.5,F5.4,3F5.2)                   
      IF (IFU.LT.1) GO TO 107                                                   
      CWN=0.16/((ALOG(HEIGHT*100.0/ZO))**2)                                     
      AIR=0.00035*PA/273.16                                                     
      FUCOEF=(AIR*0.622E6*CWN)/PA                                               
  107 CONTINUE                                                                  
      READ 917,C1,C2,C3,C4,DMETA,C5,CW1,CW2,CW3,CW4,G1,G2,G3                    
  917 FORMAT (F5.4,F5.1,F5.3,F5.2,F5.2,F5.1,F5.1,F10.5,2F5.1,3F5.0)             
      IF (C5.LT.1.0) C5=1.0                                                     
      IF (ITS.EQ.0) GO TO 100                                                   
      C1=0.0                                                                    
      C3=0.0                                                                    
      C5=0.0                                                                    
      IF (ITS.LT.2) GO TO 100                                                   
      IF (LENGTH.GT.31) LENGTH=31                                               
      IMO=1                                                                     
      IDA=1                                                                     
      IYR=1                                                                     
      LMO=1                                                                     
      LDA=LENGTH                                                                
      ICYCLE=0                                                                  
      IF (ITS.NE.3) GO TO 104                                                   
      ICYCLE=2.0*CYCLE+0.1                                                      
      IF (ICYCLE.LT.2) ICYCLE=2                                                 
      LDA=LDA+ICYCLE                                                            
      IF (LDA.GT.31) LDA=31                                                     
      ICYCLE=ICYCLE*24                                                          
  104 LYR=1                                                                     
C*******************************************************************************
C     PRINT TITLE PAGE.                                                         
  100 DO 102 I=1,2                                                              
      PRINT 904                                                                 
  904 FORMAT(1H1)                                                               
      DO 103 J=1,10                                                             
  103 PRINT 905                                                                 
  905 FORMAT(1H0)                                                               
      PRINT 906,IMO,IDA,IYR,LMO,LDA,LYR                                         
  906 FORMAT (1H ,25X,46HPOINT SNOW COVER ENERGY AND MASS BALANCE MODEL,        
     15X,I2,1H/,I2,1H/,I4,1X,2HTO,I3,1H/,I2,1H/,I4)                             
      PRINT 909,TITLE                                                           
  909 FORMAT (1H0,30X,15A4)                                                     
      PRINT 907,DELTAT,THEDA,TOLER                                              
  907 FORMAT (1H0,34HBASIC COMPUTATIONAL TIME INTERVAL=,F3.0,1X,5HHOURS,        
     15X,6HTHEDA=,F4.2,5X,32HITERATIVE TERMINATION TOLERANCE=,F5.3,1X,          
     216HDEGREES CELSIUS.)                                                      
      PRINT 910,DTONE,TEXP                                                      
  910 FORMAT (1H , 5X,92HBASIC TIME INTERVAL IS SUB-DIVIDED WHEN THE TEM        
     1PERATURE CHANGE IN THE SURFACE LAYERS EXCEEDS,F5.1,1X,6HDEG. C,5X,        
     29HEXPONENT=,F4.2)                                                         
      PRINT 911                                                                 
  911 FORMAT (1H0,50X,16HPARAMETER VALUES)                                      
      PRINT 912,TCG,DCG,DTG                                                     
  912 FORMAT (1H ,29HSOIL -- THERMAL CONDUCTIVITY=,F6.4,1X,16HCAL/CM/SEC        
     1/DEG.C,5X,22HDIFFUSION COEFFICIENT=,F5.3,1X,7HCM2/SEC,5X,                 
     221HDEPTH OF MEASUREMENT=,F3.0,1X,3HCM.)                                   
      PRINT 913,DO,X,PA                                                         
  913 FORMAT (1H ,67HSNOW -- VAPOR DIFFUSION COEFFICIENT (1000 MB,ZERO D        
     1EGREES CELSIUS)=,F5.3,1X,7HCM2/SEC,5X,9HEXPONENT=,F4.1,4X,                
     217HSTATION PRESSURE=,F5.0,1X,3HMB.)                                       
      PRINT 926,COEFKE                                                          
  926 FORMAT (1H ,8X,39HEFFECTIVE THERMAL CONDUCTIVITY=0.00005+,                
     1F6.4,13H*(DENSITY**2))                                                    
      PRINT 927,CV                                                              
  927 FORMAT (1H ,8X,17HEXTINCTION COEFF=,F5.2,27H*DENSITY*SQRT(1/GRAIN         
     1SIZE))                                                                    
      PRINT 928,G1,G2,G3                                                        
  928 FORMAT (1H ,20X,11HGRAIN SIZE=,F6.4,1H+,F5.2,14H*(DENSITY**2)+,           
     1F5.0,13H*(DENSITY**4))                                                    
      IF (NOBS.EQ.1) GO TO 102                                                  
      PRINT 914,THICK,CTHICK                                                    
  914 FORMAT (1H ,43HDESIRED THICKNESS OF THE UPPER SNOW LAYERS=,F4.1,          
     11X,3HCM.,5X,40HTHICKNESS COEFFICIENT FOR DEEPER LAYERS=,F5.3)             
      PRINT 915,PLWMAX,PLWDEN,PLWHC                                             
  915 FORMAT(1H ,46HPERCENT(DECIMAL) LIQUID-WATER-HOLDING CAPACITY,3X,          
     18HMAXIMUM=,F3.2,5X,41HBASIC CAPACITY(FOR DENSITIES IN EXCESS OF,          
     2F4.2,8H) EQUALS,F4.2)                                                     
      PRINT 918,C3,C4,DMETA                                                     
  918 FORMAT (1H ,50HMETAMORPHISM PARAMETERS--DESTRUCTIVE METAMORPHISM ,        
     13HC3=,F5.3,1X,3HC4=,F4.3,1X,14HDENSITY LIMIT=,F3.2)                       
      PRINT 919,C1,C2,C5                                                        
  919 FORMAT (1H ,25X,14HCOMPACTION C1=,F7.5,1X,3HC2=,F4.1,5X,                  
     116HLIQUID-WATER C5=,F4.1)                                                 
      PRINT 922,CW1,CW2,CW3,CW4                                                 
  922 FORMAT (1H ,42HLIQUID-WATER TRANSMISSION PARAMETERS--CW1=,F5.2,3X,        
     14HCW2=,F10.7,3X,4HCW3=,F4.1,3X,4HCW4=,F7.2)                               
      PRINT 923,HEIGHT,FUCOEF                                                   
  923 FORMAT (1H ,45HWIND FUNCTION INFORMATION--INSTRUMENT HEIGHT=,             
     1F4.1,1X,6HMETERS,2X,12HCOEFFICIENT=,F7.5,1X,9HMM/MB/KM.)                  
      IF (IFU.LT.1) PRINT 924                                                   
  924 FORMAT (1H+,90X,24HEMPIRICAL WIND FUNCTION.)                              
      IF (IFU.GT.0) PRINT 925                                                   
  925 FORMAT (1H+,90X,38HTHEORETICAL WITH STABILITY ADJUSTMENT.)                
      IF (IFU.GT.0) PRINT 929,RICRIT,ZO                                         
  929 FORMAT (1H ,40X,27HCRITICAL RICHARDSON NUMBER=,F4.2,10X,                  
     13HZ0=,F5.3,1X,11HCENTIMETERS)                                             
      PRINT 920,RCF,SCF                                                         
  920 FORMAT (1H ,41HPRECIPITATION ADJUSTMENT FACTORS -- RAIN=,F4.2,2X,         
     15HSNOW=,F4.2)                                                             
      PRINT 921,ADJQA                                                           
  921 FORMAT (1H0,42HATMOSPHERIC LONGWAVE RADIATION ADJUSTMENT=,F5.2)           
      IF (IQAE.EQ.1) PRINT 916                                                  
  916 FORMAT (1H+,55X,55HESTIMATED ATMOSPHERIC LONGWAVE RADIATION DATA A        
     1RE USED.)                                                                 
  102 CONTINUE                                                                  
C*******************************************************************************
C     INITIALIZE THE SNOW COVER.                                                
      CALL BEGIN(N,D,P,TT,WT,NOKNOW,ITS,TI,WEI)                                 
      SURDEN=P(1)                                                               
      SUMPX=0.0                                                                 
      MONTH=IMO                                                                 
      IYEAR=IYR                                                                 
      IDAY=IDA                                                                  
  101 IHOUR=0                                                                   
C*******************************************************************************
C     BEGIN MONTHLY LOOP.                                                       
      LEAPYR=0                                                                  
      IF (((IYEAR/4)*4).EQ.IYEAR)LEAPYR=1                                       
      LAST=ND(MONTH)                                                            
      IF ((MONTH.EQ.2).AND.(LEAPYR.EQ.1))LAST=29                                
C     INPUT HOURLY DATA FOR THE MONTH.   
C  LAST:#DAYS IN MONTH                
      CALL DATAIN(LAST,ITOPT,ITS,NOBS,TSMAX,TT(1),CYCLE,IMO,IYR,MONTH,          
     1IYEAR,IQAE)                                                               
C*******************************************************************************
C     BEGIN DELTA TIME COMPUTATION LOOP.                                        
  105 IHOUR=IHOUR+IDT                                                           
      IF (IHOUR.LE.24) GO TO 110                                                
      IDAY=IDAY+1                                                               
      IHOUR=IHOUR-24                                                            
C     GET VALUES FOR INPUT VARIABLES AND VERIFICATION DATA FOR                  
C        THE TIME PERIOD.                                                       
  110 CALL OBTAIN(DELTAT,IDAY,IHOUR,PNS,QI,QR,QA,FU,TA,EA,                      
     1PX,TGT,TGTDT,TT(1),TTDT(1),ITOPT,NOBS,TO,PO,WESC,WEPW,DEPTH,              
     2STAKE,PA,TSNOW,SCF,RCF,NOKNOW(1),ADJQA)                                   
      IF(NOBS.NE.1) CALL NWSNOW(N,TT,WT,TTDT,WTDT,NOKNOW,D,P,RTT,RTTDT,         
     1PX,PNS,TSNOW,SUMPX,THICK,ITOPT)                                           
      IF (N.EQ.0) GO TO 140                                                     
C     GET A FIRST GUESS FOR TEMPERATURE AND LIQUID-WATER                        
C        FOR EACH LAYER FOR TIME T+DT.                                          
      CALL GUESS(N,TT,TTDT,WT,WTDT,NOKNOW,D,P,ITOPT,QI,QR,QA,                   
     1FU,TA,EA,PX,DO,PA,X,THEDA,TOLER,DELTAT,MONTH,IDAY,IYEAR,                  
     2IHOUR,RTT,RTTDT,IGRAD,GRMAX,DTONE,TEXP,NINC,ITS,TSMAX)                    
C     COMPUTE TEMPERATURE AND LIQUID-WATER FOR EACH LAYER                       
C        FOR TIME T+DT.                                                         
      IF (N.GT.1) GO TO 111                                                     
C     SNOW COVER CONSISTS OF A SINGLE LAYER                                     
      IGUESS=0                                                                  
      CALL SURFAC(N,TT(1),TTDT(1),WT(1),WTDT(1),D,P,ITOPT,QI,QR,                
     1QA,FU,TA,EA,PX,DO,PA,X,DTG,TGT,TGTDT,TCG,DCG,THEDA,                       
     2NOKNOW(1),TOLER,DELTAT,ITER,MONTH,IDAY,IYEAR,IHOUR,IGUESS)                
      GO TO 112                                                                 
C     SNOW COVER CONSISTS OF MORE THAN ONE LAYER.                               
  111 CALL SNOWTW(N,TT,TTDT,WT,WTDT,D,P,ITOPT,QI,QR,QA,FU,                      
     1TA,EA,PX,DO,PA,X,DTG,TGT,TGTDT,TCG,DCG,THEDA,NOKNOW,TOLER,                
     2DELTAT,ITER,MONTH,IDAY,IYEAR,IHOUR,NINC)                                  
  112 IF (ITS.EQ.0)  CALL CHECK(N,D,P,TT,TTDT,WT,WTDT,NOKNOW,SCOUT,             
     1THICK,CTHICK,WE,TDEPTH,QI,QR,QA,DELTAT,PA,FU,TA,EA,PX,TGT,TGTDT,          
     2TCG,DCG,DTG,DO,X,ITOPT)                                                   
      IF (N.LT.0) GO TO 130                                                     
C     RETAIN TEMPERATURE CHANGE DUE TO HEAT TRANSFER FOR EACH LAYER.            
      CALL RETAIN(N,TT,TTDT,RTT,RTTDT)                                          
      IF ((DO.GT.0.0).OR.(ITS.EQ.0)) CALL VAPOR(N,P,DO,PA,X,TT,TTDT,D,          
     1VAPOUR,DTG,DCG,TGT,TGTDT,EA,FU,DELTAT,WTDT(1),NOBS,WTDT(N),               
     2NOKNOW(N),SOILVT)                                                         
      IF (NOBS.EQ.1) GO TO 113                                                  
  130 CALL WATER(DELTAT,N,P,D,WTDT,TTDT,NOKNOW,PLWHC,PLWMAX,PLWDEN,PX,          
     1SCOUT,MONTH,IDAY,IYEAR,IHOUR,IPUNCH)                                      
      IF (N.LT.0) GO TO 140                                                     
  113 CALL META(N,P,D,WTDT,TTDT,PLWHC,C1,C2,C3,C4,C5,DELTAT,DMETA)              
C     CHECK IF PRINTER OUTPUT IS WANTED FOR THIS TIME PERIOD.                   
      MHR=(IDAY-1)*24+IHOUR                                                     
      IF (((MHR/IOUT)*IOUT).NE.MHR) GO TO 125                                   
C     A PRINT IS WANTED FOR THIS TIME PERIOD.                                   
      IF (ITS.LT.2) GO TO 120                                                   
      IF (MHR.LE.ICYCLE) GO TO 150                                              
C     OUTPUT FOR SELECTED SNOW SURFACE TEMPERATURE PATTERNS.                    
      CALL PATERN (N,TTDT,WTDT,ITS,TSMAX,CYCLE,D,P,IDAY,                        
     1IHOUR,ITER,DO,TI,SURDEN)                                                  
      GO TO 150                                                                 
C     OUTPUT FOR OTHER CASES                                                    
  120 CALL SNOWOT(N,TTDT,WTDT,D,P,MONTH,IDAY,IYEAR,IHOUR,ITER,WE,TDEPTH)        
  125 IF (NOBS.EQ.1) GO TO 150                                                  
  140 CALL STATDA(TT(1),TTDT(1),TO,SCOUT,PO,WE,WESC,WEPW,VAPOUR,TDEPTH,         
     1DEPTH,STAKE,MONTH,IDAY,IYEAR,IHOUR,DELTAT,N,SOILVT,TA)                    
C     END OF DELTA TIME LOOP, CHECK FOR END OF RUN AND MONTH.                   
C*******************************************************************************
  150 IF (IHOUR.LT.24) GO TO 105                                                
      IF (MONTH.NE.LMO) GO TO 160                                               
      IF ((IYEAR.EQ.LYR).AND.(IDAY.EQ.LDA)) GO TO 200                           
      GO TO 105                                                                 
  160 IF (IDAY.LT.LAST) GO TO 105                                               
      MONTH=MONTH+1                                                             
      IDAY=1                                                                    
      IF (MONTH.LE.12) GO TO 101                                                
      MONTH=1                                                                   
      IYEAR=IYEAR+1                                                             
      GO TO 101                                                                 
C     END OF THE RUN                                                            
C*******************************************************************************
  200 IF (NOBS.NE.1) CALL FINAL(IMO,IDA,IYR,LMO,LDA,LYR,WEI,SUMPX,WE)           
      STOP                                                                      
      END                                                                       

      SUBROUTINE DATAIN (LAST,ITOPT,ITS,NOBS,TSMAX,TST,CYCLE,                   
     1IMO,IYR,MONTH,IYEAR,IQAE)                                                 
C*******************************************************************************
C     INPUTS BASIC HOURLY DATA FROM TAPE (OR DISK) AND CARDS,                   
C        PLUS CONVERTS UNITS.                                                   
C*******************************************************************************
      COMMON/DADATA/WEP(31),DSS(31),WET(31),DSC(31),SURDEN(31)                  
      COMMON/HRDATA/TG(744),TS(744),DENNS(744),TAD(744),EAD(744),               
     1UA(744),QID(744),QRD(744),QAD(744),POD(744),PXD(744)                      
C*******************************************************************************
      DIMENSION ID(5),IH(5),GT(5),ST(24),SD(5),NH(5)                            
      DIMENSION TEMPOR(744)                                                     
      DATA THERMO,RAIN,SNOW,FICE/1HT,1HR,1HS,1HI/                               
      DATA B1/1H /                                                              
      LMHR=LAST*24                                                              
C*******************************************************************************
C     CARD SECTION - SOIL TEMPERATURE ALWAYS,SNOW SURFACE                       
C        TEMPERATURE AND DENSITY OF NEW SNOW UNDER CERTAIN OPTIONS.             
C*******************************************************************************
C     INPUT SOIL TEMPERATURE IN DEGREES FAHRENHEIT.                             
C        NEED VALUES FOR THE LAST HOUR OF THE MONTH AND                         
C        ALL OTHER HOURS WHEN CHANGES OCCUR.                                    
      MHR1=1                                                                    
  102 READ 900,MO,JYR,((ID(I),IH(I),GT(I)),I=1,5)                               
  900 FORMAT (I2,1X,I2,5(I5,I5,F5.0))                                           
      DO 100 I=1,5                                                              
      IF (ID(I).LT.1) STOP 30                                                   
      MHR2=(ID(I)-1)*24+IH(I)                                                   
      DO 101 MHR=MHR1,MHR2                                                      
      TG(MHR)=((5.0/9.0)*(GT(I)-32.0))+273.16                                   
  101 CONTINUE                                                                  
      IF (MHR2.EQ.LMHR) GO TO 105                                               
      MHR1=MHR2+1                                                               
  100 CONTINUE                                                                  
      GO TO 102                                                                 
C*******************************************************************************
C     INPUT SNOW SURFACE TEMPERATURE ITOPT=1,ITS=1.                             
  105 IF (ITOPT.NE.1) GO TO 120                                                 
      IF (ITS.NE.1) GO TO 110                                                   
C     NEED VALUES FOR EACH HOUR OF THE MONTH - ONE CARD PER DAY                 
C        UNITS ARE DEGREES FAHRENHEIT.                                          
      DO 106 IDA=1,LAST                                                         
      READ 901,MO,IDAY,JYR,ST                                                   
  901 FORMAT (I2,1X,I2,1X,I2,24F3.0)                                            
      IF (IDA.NE.IDAY) STOP 31                                                  
      MHR1=(IDA-1)*24                                                           
      DO 107 IHR=1,24                                                           
      MHR=MHR1+IHR                                                              
      TS(MHR)=((5.0/9.0)*(ST(IHR)-32.0))+273.16                                 
  107 CONTINUE                                                                 
  106 CONTINUE                                                                  
      GO TO 120                                                                 
  110 IF(ITS.LT.2) GO TO 120                                                    
C*******************************************************************************
C     PRESET SNOW SURFACE TEMPERATURE PATTERN-NO OBSERVED DATA ALLOWED.         
      NOBS=1                                                                    
      IF (ITS.GT.2) GO TO 112                                                   
C     ITS=2,SUDDEN CHANGE FROM INITIAL VALUE BY AN AMOUNT TSMAX.                
C        SURFACE TEMPERATURE THEN REMAINS CONSTANT                              
      DO 111 MHR=1,744                                                          
      TS(MHR)=TST+TSMAX                                                         
      IF (TS(MHR).GT.273.16)TS(MHR)=273.16                                      
  111 CONTINUE                                                                  
      TST=TST+TSMAX                                                             
      IF (TST.GT.273.16)TST=273.16                                              
      GO TO 120                                                                 
  112 IF (ITS.GT.3) STOP 32                                                     
C     ITS=3,SIN-VARIATION CHANGE AT THE SURFACE ABOUT TST.                      
C        AMPLITUDE=TSMAX,CYCLE=THE TIME IN DAYS FOR                             
C        ONE COMPLETE SURFACE TEMPERATURE VARIATION.                            
      PI2=6.2832                                                                
      DO 113 MHR=1,744                                                          
C     COMPUTE THE FRACTION OF THE CYCLE                                         
      HOUR=MHR                                                                  
      FRAC=HOUR/(CYCLE*24.0)                                                    
C     ONE CYCLE IS 2PI RADIANS.                                                 
C     CONVERT TO RADIANS                                                        
      FRAC=FRAC*PI2                                                             
      TS(MHR)=TST+TSMAX*SIN(FRAC)                                               
      IF (TS(MHR).GT.273.16)TS(MHR)=273.16                                      
  113 CONTINUE                                                                  
  120 IF (NOBS.NE.1) GO TO 124                                                  
      IF (ITS.EQ.0) GO TO 125                                                   
C     NO OBSERVED DATA ARE BEING USED                                           
      RETURN                                                                    
C*******************************************************************************
C     INPUT NEW SNOW DENSITY FOR HOURS WHEN THE COMPUTATION                     
C        IS TO BE OVER-RIDDEN. FIVE PERIODS PER CARD.                           
C        THIS CAN BE USED TO CHANGE THE FORM OF THE PRECIPITATION ALSO.         
C                  =1.0, INDICATES RAIN                                         
C                  (.GT.0.0).AND.(.LT.0.90) INDICATES SNOW.                     
C     PRESET ARRAY TO NEGATIVE.                                                 
  124 DO 121 MHR=1,744                                                          
  121 DENNS(MHR)=-1.0                                                           
  123 READ 902,MO,JYR,IEND,((ID(I),IH(I),NH(I),SD(I)),I=1,5)                    
  902 FORMAT (2I2,I1,5(I4,2I3,F5.2))                                            
      DO 122 I=1,5                                                              
      IF (ID(I).LT.1) GO TO 125                                                 
      MHR1=(ID(I)-1)*24+IH(I)                                                   
      MHR2=MHR1+NH(I)-1                                                         
      DO 119 MHR=MHR1,MHR2                                                      
      DENNS(MHR)=SD(I)                                                          
  119 CONTINUE                                                                  
  122 CONTINUE                                                                  
      IF (IEND.EQ.1)  GO TO 125                                                 
      GO TO 123                                                                 
C*******************************************************************************
C     REMAINING DATA FROM TAPE OR DISK                                          
C*******************************************************************************
C     FORMAT OF THE SEQUENTIAL FILE CONTAINING METEOROLOGICAL DATA IS AS        
C         FOLLOWS.                                                              
C              BINARY FILE COMPOSED OF 13 LOGICAL RECORDS PER MONTH.            
C                 FIRST 12 RECORDS CONTAIN 748 WORDS EACH.                      
C                   WORD 1=ALPHANUMERIC DATA TYPE CODE                          
C                        2=1 IF IN METRIC UNITS,  =0 IF IN RECORDED UNITS       
C                        3=MONTH NUMBER                                         
C                        4=YEAR (4 DIGITS)                                      
C                    5-748=HOURLY DATA FOR EACH HOUR OF THE MONTH.              
C                                                                               
C                        THE DATA CONTAINED IN THESE 12 RECORDS ARE             
C                            1. AIR TEMPERATURE                                 
C                            2. VAPOR PRESSURE OF THE AIR                       
C                            3. WIND SPEED                                      
C                            4. INCOMING SOLAR RADIATION                        
C                            5. REFLECTED SOLAR RADIATION                       
C                            6. OBSERVED INCOMING LONGWAVE RADIATION            
C                            7. ESTIMATED INCOMING LONGWAVE RADIATION           
C                            8. OBSERVED SNOW SURFACE TEMPERATURE               
C                            9. ALPHANUMERIC CHARACTER IDENTIFYING SNOW         
C                                  SURFACE TEMPERATURE INSTRUMENT IN USE.       
C                                     BLANK=INFRARED THERMOMETER                
C                                     T=THERMOCOUPLE  M=MISSING                 
C                           10. OBSERVED SNOW COVER OUTFLOW.                    
C                           11. PRECIPITATION                                   
C                           12. FORM OF PRECIPITATION (ALPHANUMERIC)            
C                                 R=RAIN, S=SNOW,  I=ICE                        
C                                                                               
C                 THE 13TH RECORD CONTAINS 159 WORDS.                           
C                   WORD 1=DAILY DATA CODE =PDO                                 
C                        2=1 IF IN METRIC UNITS, =0 IF IN RECORDED UNITS        
C                        3=MONTH NUMBER                                         
C                        4=YEAR (4 DIGITS)                                      
C                     5-35=SNOW PILLOW WATER-EQUIVALENT FOR EACH DAY            
C                    36-66=DEPTH FROM SNOW STAKES FOR EACH DAY                  
C                    67-97=SNOW COURSE WATER-EQUIVALENT FOR EACH DAY            
C                   98-128=SNOW COURSE DEPTH FOR EACH DAY                       
C                  129-159=SNOW SURFACE DENSITY                                 
C                                                                               
C              MISSING DATA IS INDICATED BY 9999.0.                             
C                                                                               
C    EACH SUCCEEDING MONTH USED IN A SIMULATION RUN CAN NOT BE STORED           
C         IN THE FILE PRIOR TO THE PREVIOUS MONTH.                              
C*******************************************************************************
C     POSITION TAPE AT FIRST MONTH.                                             
  125 IF ((MONTH.NE.IMO).OR.(IYEAR.NE.IYR)) GO TO 130                           
      REWIND 1                                                                  
  127 READ (1) TYPE,METRIC,JMO,JYR,TAD                                          
      IF ((JMO.EQ.MONTH).AND.(JYR.EQ.IYEAR)) GO TO 131                          
      DO 126 I=1,12                                                             
      READ (1)                                                                  
  126 CONTINUE                                                                  
      GO TO 127                                                                 
C*******************************************************************************
C     READ DATA FOR THE MONTH -- CONVERTED TO PROPER UNITS AND STORE.           
C*******************************************************************************
C     AIR TEMPERATURE -- DEGREES KELVIN.                                        
  130 READ(1)TYPE,METRIC,JMO,JYR,TAD                                            
      IF ((JMO.NE.MONTH).OR.(JYR.NE.IYEAR)) STOP 33                             
  131 DO 132 MHR=1,LMHR                                                         
      IF (METRIC.EQ.0)TAD(MHR)=(5.0/9.0)*(TAD(MHR)-32.0)                        
      TAD(MHR)=TAD(MHR)+273.16                                                  
  132 CONTINUE                                                                  
C*******************************************************************************
C     VAPOR PRESSURE -- MILLIBARS.                                              
      READ(1)TYPE,METRIC,JMO,JYR,EAD                                            
      IF (METRIC.EQ.1) GO TO 135                                                
      DO 133 MHR=1,LMHR                                                         
  133 EAD(MHR)=33.864*EAD(MHR)                                                  
C*******************************************************************************
C     WIND SPEED -- METERS/SECOND.                                              
  135 READ(1)TYPE,METRIC,JMO,JYR,UA                                             
      IF (METRIC.EQ.1) GO TO 140                                                
      DO 136 MHR=1,LMHR                                                         
  136 UA(MHR)=0.447*UA(MHR)                                                     
C*******************************************************************************
C     RADIATION DATA                                                            
C     INCOMING SOLAR RADIATION -- LANGLEYS(CAL/CM2).                            
  140 READ(1)TYPE,METRIC,JMO,JYR,QID                                            
C     REFLECTED SOLAR RADIATION -- LANGLEYS.                                    
      READ(1)TYPE,METRIC,JMO,JYR,QRD                                            
C     ATMOSPHERIC LONGWAVE RADIATION -- LANGLEYS.                               
      IF (IQAE.EQ.1) GO TO 141                                                  
C     USE MEASURED ATMOSPHERIC LONGWAVE.                                        
      READ(1)TYPE,METRIC,JMO,JYR,QAD                                            
  141 READ(1)TYPE,METRIC,JMO,JYR,TEMPOR                                         
      IF (IQAE.NE.1) GO TO 145                                                  
C     USE ESTIMATED ATMOSPHERIC LONGWAVE.                                       
      READ(1)TYPE,METRIC,JMO,JYR,QAD                                            
C*******************************************************************************
C     SNOW SURFACE TEMPERATURE -- DEGREES KELVIN                                
  145 IF (ITS.GT.0) GO TO 149                                                   
      READ(1)TYPE,METRIC,JMO,JYR,TS                                             
      DO 146 MHR=1,LMHR                                                         
      IF (TS(MHR).GT.9000.0) GO TO 146                                          
      IF (METRIC.EQ.0)TS(MHR)=(5.0/9.0)*(TS(MHR)-32.0)                          
      TS(MHR)=TS(MHR)+273.16                                                    
  146 CONTINUE                                                                  
      READ(1)TYPE,METRIC,JMO,JYR,TEMPOR                                         
C     IF TO BE USED FOR VERIFICATION,OFFSET THERMOCOUPLE VALUES                 
C             BY 100 DEGREES.                                                   
      IF (ITOPT.GT.0) GO TO 150                                                 
      DO 147 MHR=1,LMHR                                                         
      IF (TS(MHR).GT.9000.0) GO TO 147                                          
      IF (TEMPOR(MHR).EQ.THERMO)TS(MHR)=TS(MHR)+100.0                           
  147 CONTINUE                                                                  
      GO TO 150                                                                 
C     SNOW SURFACE TEMPERATURE NOT NEEDED FROM TAPE.                            
  149 READ(1)TYPE,METRIC,JMO,JYR,TEMPOR                                         
      READ(1)TYPE,METRIC,JMO,JYR,TEMPOR                                         
C*******************************************************************************
C     SNOW COVER OUTFLOW -- MILLIMETERS                                         
  150 READ(1)TYPE,METRIC,JMO,JYR,TEMPOR                                         
      DO 151 MHR=1,LMHR                                                         
      IF (TEMPOR(MHR).GT.9000.0) GO TO 153                                      
      IF (TEMPOR(MHR).LT.0.0) GO TO 153                                         
      IF (MHR.EQ.1) GO TO 152                                                   
      IF (TEMPOR(MHR-1).GE.0.0) GO TO 152                                       
  153 POD(MHR)=9999.0                                                           
      GO TO 151                                                                 
  152 POD(MHR)=TEMPOR(MHR)                                                      
      IF (METRIC.EQ.0)POD(MHR)=25.4*POD(MHR)                                    
  151 CONTINUE                                                                  
C*******************************************************************************
C     PRECIPITATION -- MILLIMETERS                                              
      READ(1)TYPE,METRIC,JMO,JYR,PXD                                            
      IF (METRIC.EQ.1) GO TO 160                                                
      DO 156 MHR=1,LMHR                                                         
  156 PXD(MHR)=25.4*PXD(MHR)                                                    
C     READ FORM OF PRECIPITATION -- USE DENSITY OF NEW SNOW                     
C             ARRAY TO COMMUNICATE TO REST OF PROGRAM.                          
C             DENSITY OF NEW SNOW = NEGATIVE,EITHER NO PRECIPITATION            
C                                      OR FORM NOT AVAILABLE.                   
C                                 = POSITIVE,BUT LESS THAN 1.0 -- VALID         
C                                      VALUE,OBTAINED FROM ICE OR WAS           
C                                      INPUT PREVIOUSLY.                        
C                                 = 1.0, RAIN                                   
C                                 = GREATER THAN ONE,SNOW -- DENSITY            
C                                      TO BE COMPUTED FROM                      
C                                      TEMPERATURE LATER.                       
  160 READ(1)TYPE,METRIC,JMO,JYR,TEMPOR                                         
      DO 161 MHR=1,LMHR                                                         
      IF (PXD(MHR).LT.0.001) GO TO 161                                          
      IF (TEMPOR(MHR).EQ.B1) GO TO 161                                          
      IF (DENNS(MHR).GT.0.0) GO TO 161                                          
      IF (TEMPOR(MHR).EQ.RAIN)DENNS(MHR)=1.0                                    
      IF (TEMPOR(MHR).EQ.SNOW)DENNS(MHR)=2.0                                    
      IF (TEMPOR(MHR).EQ.FICE)DENNS(MHR)=0.8                                    
  161 CONTINUE                                                                  
C*******************************************************************************
C     DAILY OBSERVATIONS -- USE FOR VERIFICATION.                               
      READ(1)TYPE,METRIC,JMO,JYR,WEP,DSS,WET,DSC,SURDEN                         
C     CONVERT UNITS -- DEPTH(CM) -- WATER EQUIVALENT(MM)                        
      IF (METRIC.EQ.1) GO TO 169                                                
      DO 165 KD=1,LAST                                                          
      IF (DSS(KD).GT.9000.0) GO TO 166                                          
      DSS(KD)=2.54*DSS(KD)                                                      
  166 IF (DSC(KD).GT.9000.0) GO TO 167                                          
      DSC(KD)=2.54*DSC(KD)                                                      
  167 IF (WEP(KD).GT.9000.0) GO TO 168                                          
      WEP(KD)=25.4*WEP(KD)                                                      
  168 IF (WET(KD).GT.9000.0) GO TO 165                                          
      WET(KD)=25.4*WET(KD)                                                      
  165 CONTINUE                                                                  
  169 CONTINUE                                                                  
C*******************************************************************************
      RETURN                                                                    
      END                                                                       

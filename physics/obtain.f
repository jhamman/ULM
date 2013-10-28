      SUBROUTINE OBTAIN(DELTAT,IDAY,IHOUR,PNS,QI,QR,QA,FU,TA,EA,                
     1PX,TGT,TGTDT,TST,TSTDT,ITOPT,NOBS,TO,PO,WESC,WEPW,                        
     2DEPTH,STAKE,PA,TSNOW,SCF,RCF,KNOWN,ADJQA,
     3TG,DENNS,TAD,EAD,UA,QID,QRD,QAD,POD,PXD,NEWPACK,
     4IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO)            
C*******************************************************************************
C     ESTABLISHES DATA VALUES FOR THE TIME PERIOD.                              
C*******************************************************************************
C Replace the common blocks below (originally from DATAIN) to direct passes on
C line 3 above
      COMMON/DADATA/WEP(31),DSS(31),WET(31),DSC(31),SURDEN(31)                  
C      COMMON/HRDATA/TG,TS,DENNS,TAD,EAD,               
C     1UA,QID,QRD,QAD,POD,PXD                      
C      COMMON/IFU,UMSEC,FUCOEF,HEIGHT                                     
C*******************************************************************************
C     COMPUTE FIRST AND LAST HOUR IN THIS TIME PERIOD.                          
      MHR2=(IDAY-1)*24+IHOUR                                                    
      IDT=DELTAT+0.01                                                           
      MHR1=MHR2-IDT+1                                                           
      MHRT=MHR1-1                                                               
      IF (MHRT.EQ.0) MHRT=1                                                     
      IF (ITOPT.NE.1) GO TO 100                                                 
C*******************************************************************************
C     BEGINNING AND END OF PERIOD SNOW SURFACE TEMPERATURE                      
C        FOR THE CASE WHEN ITOPT=1.                                             
      IF (KNOWN.GE.0) TST=TSTDT                                                 
      TSTDT=TS                                                            
C*******************************************************************************
C     BEGINNING AND ENDING SOIL TEMPERATURE                                     
  100 TGT=TG                                                              
      TGTDT=TG                                                            
      IF (NOBS.NE.1) GO TO 110                                                  
C*******************************************************************************
C     NO OBSERVED INPUT OR VERIFICATION DATA                                    
      IF (ITOPT.NE.1)STOP 41                                                    
      PNS=-1.0                                                                  
      QI=0.0                                                                    
      QR=0.0                                                                    
      QA=0.0                                                                    
      FU=0.0                                                                    
      TA=0.0                                                                    
      EA=0.0                                                                    
      PX=0.0                                                                    
      TSNOW=0.0                                                                 
      TO=9999.0                                                                 
      PO=9999.0                                                                 
      WESC=9999.0                                                               
      WEPW=9999.0                                                               
      DEPTH=9999.0                                                              
      STAKE=9999.0                                                              
      RETURN                                                                    
C*******************************************************************************
C     COMPUTE INPUT AND VERIFICATION VALUES FROM HOURLY DATA                    
C             COMPUTE VARIOUS TOTALS FOR THE PERIOD. NO CHECK MADE FOR          
C             MISSING DATA. USER MUST BE SURE ALL NEEDED DATA ARE               
C             AVAILABLE FOR THE PERIOD BEING USED.                              
  110 TA=0.0                                                                    
      EA=0.0                                                                    
      U=0.0                                                                     
      QI=0.0                                                                    
      QR=0.0                                                                    
      QA=0.0                                                                    
      PX=0.0                                                                    
      PO=0.0                                                                    
C      DO 111 MHR=MHR1,MHR2                                                      
      TA=TA+TAD                                                            
      EA=EA+EAD                                                            
      U=U+UA                                                               
      QI=QI+QID                                                            
      QR=QR+QRD                                                            
      QA=QA+QAD*ADJQA                                                      
      PX=PX+PXD                                                            
      PO=PO+POD                                                            
C  111 CONTINUE                                                                  
C     COMPUTE MEANS WHERE NEEDED                                                
C      TA=TA/DELTAT                                                              
C      EA=EA/DELTAT                                                              
C      U=U/DELTAT                                                                
      UMSEC=U                                                                   
C*******************************************************************************
C     MAKE NEEDED CHECKS AND CONVERSIONS                                        
      IF (QI.LT.0.1) GO TO 117                                                  
      A=QR/QI                                                                   
      IF (A.GT.0.90) QR=0.90*QI                                                 
      IF (A.LT.0.10) QR=0.10*QI                                                 
      GO TO 118                                                                 
  117 QR=0.0     
C Convert precipitation to centimeters
  118 PX=PX*0.1                                                                 
      IF (PO.GT.9000.0) PO=9999.0                                               
      IF (PX.LT.0.0001) PX=0.0                                                  
C     ONLY USE INFRARED RADIOMETER DATA FOR SURFACE TEMPERATURE                 
C             VERIFICATION.                                                     
      IF (ITOPT.EQ.1) GO TO 116                                                 
      TO=0.0                                                                    
C      DO 113 MHR=MHR1,MHR2                                                      
C      IF (TS(MHR).GT.300.0) GO TO 114                                           
C      TO=TO+TS(MHR)                                                             
C      GO TO 113                                                                 
C  114 TO=TO+9999.0                                                              
  113 CONTINUE                                                                  
C      IF (TO.GT.9000.0) GO TO 116                                               
C      TO=TO/DELTAT                                                              
C      TO=TO-273.16                                                              
C      GO TO 115                                                                 
  116 TO=9999.0                                                                 
C*******************************************************************************
C     DAILY VALUES                                                              
C  115 WESC=WET(IDAY)                                                            
C      WEPW=WEP(IDAY)                                                            
C      DEPTH=DSC(IDAY)                                                           
C      STAKE=DSS(IDAY)                                                           
C*******************************************************************************
C     DETERMINE THE WIND FUNCTION FOR THE PERIOD.                               
C               (IF EMPIRICAL WIND FUNCTION IS TO BE USED.)  
C Provide for TSTDT as needed below for a new pack to compute the Richardson #          
      IF (IFU.LT.1) THEN
         IF (NEWPACK == 1) THEN
            TSTDT = TST
         ENDIF
         CALL WINDF(DELTAT,UMSEC,TA,TST,TSTDT,FUCOEF,PA,HEIGHT,IFU,FU,
     &    RICRIT,ZO)
      ENDIF              
C*******************************************************************************
C     DETERMINE DENSITY OF NEW SNOW -- FINAL VALUES.                            
C             NEGATIVE IF NO PRECIPITATION.                                     
C             =1.0 IF RAIN.                                                     
C             .GT.0.0.AND.LT.0.90 IF SNOW.                                      
      PNS=-1.0                                                                  
      TSNOW=0.0                                                                 
      IF (PX.EQ.0.0) RETURN                                                     
C     COMPUTE PERIOD WET-BULB TEMPERATURE.                                      
      CALL WTBULB(TW,TA,EA,PA)                                                  
C     CHECK FOR FROZEN PRECIPITATION.   
      WRITE(*,*)'TW TA EA PA',TW,TA,EA,PA                                        
      ICE=0                                                                     
      IF (TW.LE.274.16) ICE=1                                                   
C     CONVERT FORM IF INPUT VALUES FOR ALL HOURS ARE REVERSED.                  
      HOURS=0.1                                                                 
C      DO 120 MHR=MHR1,MHR2                                                      
      IF (ICE.EQ.1) GO TO 121                                                   
      IF ((DENNS.GT.0.0).AND.(DENNS.NE.1.0)) HOURS=HOURS+1.0          
      GO TO 120                                                                 
  121 IF (DENNS.EQ.1.0) HOURS=HOURS+1.0                                    
  120 CONTINUE                                                                  
      IF (HOURS.LT.DELTAT) GO TO 125                                            
C     PRECIPIATION FORM IS OPPOSITE TO THAT INDICATED BY THE WET-BULB.          
C      IF (ICE.EQ.0) INDEX=1                                                     
C      IF (ICE.EQ.1) INDEX=0                                                     
C      ICE=INDEX                                                                 
  125 IF (ICE.EQ.1) GO TO 130                                                   
C*******************************************************************************
C     PRECIPITATION IS RAIN                                                     
      PNS=1.0                                                                   
      TSNOW=0.0                                                                 
      PX=RCF*PX                                                                 
      RETURN                                                                    
C*******************************************************************************
C     PRECIPITATION IS FROZEN -- COMPUTE DENSITY                                
  130 PNS=0.0                                                                   
C      DO 131 MHR=MHR1,MHR2                                                      
      DEN=DENNS                                                            
      IF ((DEN.GT.0.0).AND.(DEN.LT.0.90)) GO TO 133                             
C     COMPUTE DENSITY OF NEW SNOW BASED ON WET BULB TEMPERATURE.                
C             ALTA RELATIONSHIP.                                                
      CALL WTBULB(TW,TAD,EAD,PA)                                      
      IF (TW.LE.258.16) GO TO 132                                               
      DEN=0.05+0.0017*((TW-258.16)**1.5)                                        
      GO TO 133                                                                 
  132 DEN=0.05                                                                  
  133 PNS=PNS+DEN*PXD*0.1                                                  
      IF (TW.GT.273.16)TW=273.16                                                
      TSNOW=TSNOW+TW*PXD*0.1                                               
  131 CONTINUE                                                                  
      PNS=PNS/PX                                                                
      TSNOW=TSNOW/PX                                                            
      PX=SCF*PX                                                                 
C*******************************************************************************
      RETURN                                                                    
      END                                                                       

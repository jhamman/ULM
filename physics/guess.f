      SUBROUTINE GUESS(N,TT,TTDT,WT,WTDT,NOKNOW,D,P,ITOPT,                      
     1QI,QR,QA,FU,TA,EA,PX,D_O,PA,X,THEDA,TOLER,DELTAT,MONTH,                    
     2IDAY,IYEAR,IHOUR,RTT,RTTDT,IGRAD,GRMAX,DTONE,TEXP,NINC,ITS,TSMAX,
     3EXTINC,CV,G1,G2,G3,QG,IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE)         
C*******************************************************************************
C     COMPUTES THE FIRST GUESS FOR TEMPERATURE AND LIQUID-WATER AT TIME         
C        T+DT FOR EACH LAYER OF THE SNOW COVER.                                 
C*******************************************************************************
      REAL LF                                                                   
      COMMON/DHEAT/DHS                                                          
      DIMENSION TT(100),TTDT(100),WT(100),WTDT(100),NOKNOW(100),                
     1D(100),P(100),RTT(100),RTTDT(100),Z(100),AQI(100)                         
      DATA D1,D2,D3,D4,D5/1.0,1.0,1.0,1.0,1.0/                                  
C*******************************************************************************
C     ESTABLISH VALUES OF TEMPERATURE AND LIQUID WATER AT TIME T FROM           
C     THE VALUE AT THE END OF THE LAST TIME PERIOD. VALUES AT TIME T            
C     FOR NEWLY ESTABLISHED LAYERS HAVE BEEN DEFINED.                           
      IGUESS=1                                                                  
      NINC=1                                                                    
      IF (NOKNOW(1).LT.0) DHS=0.0                                               
      DO 100 I=1,N                                                              
      IF (NOKNOW(I).LT.0) GO TO 104                                             
      IF ((I.EQ.1).AND.(ITOPT.EQ.1)) GO TO 102                                  
      TT(I)=TTDT(I)                                                             
  102 WT(I)=WTDT(I)                                                             
      GO TO 100                                                                 
  104 RTT(I)=-1.0                                                               
  100 CONTINUE                                                                  
C     SET T+DT VALUE EQUAL TO TIME T VALUE TO START WITH. THEN                  
C        REFINE GUESS FROM THAT VALUE.                                          
      DO 101 I=1,N                                                              
      IF ((I.EQ.1).AND.(ITOPT.EQ.1)) GO TO 103                                  
      TTDT(I)=TT(I)                                                             
  103 WTDT(I)=WT(I)                                                             
  101 CONTINUE                                                                  
C*******************************************************************************
C     MAKE FIRST GUESS FOR THE SURFACE LAYER.                                   
      IF (NOKNOW(1).LT.0) NOKNOW(1)=NOKNOW(1)+2                                 
      IF (ITOPT.EQ.1) GO TO 110                                                 
      IF (IGRAD.NE.1) DHS=0.0     
      WRITE(*,*)'GUESS: ABOVE SURFAC WT TT',WT(1),TT(1)                           
      CALL SURFAC(N,TT(1),TTDT(1),WT(1),WTDT(1),D,P,ITOPT,QI,QR,                
     1QA,FU,TA,EA,PX,D_O,PA,X,D1,D2,D3,D4,D5,THEDA,NOKNOW(1),TOLER,              
     2DELTAT,ITER,MONTH,IDAY,IYEAR,IHOUR,IGUESS,QG,
     3EXTINC,CV,G1,G2,G3,IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE,ITS)
      WRITE(*,*)'GUESS: BELOW SURFAC WT TT',WT(1),TT(1)
      IF (NOKNOW(1).GT.0) GO TO 105                                             
      IF (TTDT(1).GT.273.16) TTDT(1)=273.16                                     
      GO TO 110                                                                 
  105 IF (WTDT(1).LT.0.0) WTDT(1)=0.0                                           
C*******************************************************************************
C     REMAINING LAYERS                                                          
  110 IF (N.EQ.1) RETURN                                                        
      CALL ZDEPTH(N,D,Z,TDEPTH)                                                 
      LF=79.7                                                                   
      PW=1.0                                                                    
      CALL ABSORB(N,QI,QR,Z,D,P,AQI,EXTINC,CV,G1,G2,G3)                        
      DO 120 I=2,N                                                              
      IF (NOKNOW(I).GE.0) GO TO 115                                             
C     LAYERS WHICH ARE ALL OR PARTLY NEW SNOW,JUST USE TIME T VALUE.            
      NOKNOW(I)=NOKNOW(I)+2                                                     
      GO TO 120                                                                 
  115 IF (NOKNOW(I).GT.0) GO TO 118                                             
C     TEMPERATURE UNKNOWN                                                       
C     FORWARD LINEAR PROJECTION OF PREVIOUS HEAT TRANSFER CHANGE                
      TTDT(I)=TT(I)+(RTTDT(I)-RTT(I))                                           
      IF (TTDT(I).GT.273.16)TTDT(I)=273.16                                      
      IF (IGRAD.NE.1) GO TO 120                                                 
C     USE GRADIENT RATIO TO TRY TO IMPROVE FIRST GUESS. (UPPER 3 LAYERS ONLY.)  
      IF (I.GT.3) GO TO 120                                                     
C     CHECK IF  ABOVE   LAYER EXSISTED DURING THE LAST TIME PERIOD.             
      IF (RTT(I-1).LT.0.0) GO TO 120                                            
      G1=RTT(I-1)-RTT(I)                                                        
      G2=RTTDT(I-1)-RTTDT(I)                                                    
      G3=TT(I-1)-TT(I)                                                          
      G4=TTDT(I-1)-TTDT(I)                                                      
      GPAST=(G1+G2)*0.5                                                         
      GNOW=(G3+G4)*0.5                                                          
      IF (ABS(GPAST).GT.0.5) GO TO 116                                          
      RATIO=0.0                                                                 
      GO TO 117                                                                 
  116 RATIO=GNOW/GPAST                                                          
      ARATIO=ABS(RATIO)                                                         
      IF (ARATIO.GT.GRMAX) RATIO=GRMAX*RATIO/ARATIO                             
C     GRMAX IS THE MAXIMUM GRADIENT RATIO WHICH WILL                            
C        BE USED TO PROJECT NEW FIRST GUESS. NECESSARY                          
C        BECAUSE SMALL GPAST CREATE LARGE RATIO.                                
  117 TTDT(I)=TT(I)+RATIO*(RTTDT(I)-RTT(I))                                     
      IF (TTDT(I).GT.273.16) TTDT(I)=273.16                                     
      GO TO 120                                                                 
C     LIQUID-WATER UNKNOWN                                                      
  118 WTDT(I)=WT(I)+(AQI(I)/(LF*PW))                                            
  120 CONTINUE                                                                  
C*******************************************************************************
C     CHECK TO DETERMINE IF THE COMPUTATIONAL TIME INTERVAL                     
C             SHOULD BE SUB-DIVIDED DURING THIS PERIOD.                         
C     SUB-DIVIDE WHEN THERE IS A LARGE CHANGE PROJECTED                         
C             IN THE UPPER TWO LAYERS.                                          
C     NO SUB-DIVIDING WHEN ONLY ONE LAYER.                                      
      IF (THEDA.GT.0.0) GO TO 130                                               
C  EXPLICIT -- USE STABILITY CRITERION, DT=0.25*(DZ**2)/THERMAL DIFFUSIVITY.    
      PMAX=P(1)                                                                 
      DMIN=D(1)                                                                 
      DO 132 I=2,N                                                              
      IF (P(I).GT.PMAX) PMAX=P(I)                                               
      IF (D(I).LT.DMIN) DMIN=D(I)                                               
  132 CONTINUE                                                                  
C       ASSUME KE=0.00005+0.0085*(DENSITY**2)                                   
C       ASSUME CI=0.4      (BOTH ASSUMPTIONS ARE CONSERVATIVE)                  
      ALPHA=(0.00005+0.0085*PMAX*PMAX)/(PMAX*0.4)                               
      DT=0.25*DMIN*DMIN/ALPHA                                                   
C      NINC=3600.0*DELTAT/DT+0.5   
C NINC uses DELTAT in seconds, thus make the change in its definition     
      NINC=DELTAT/DT+0.5                                                 
      IF (NINC.LT.2) NINC=2                                                     
      GO TO 131                                                                 
  130 CHANGE=0.0                                                                
      DO 121 I=1,2                                                              
  121 CHANGE=CHANGE+ABS(TTDT(I)-TT(I))                                          
      IF (ITS.NE.2) GO TO 123                                                   
      ITIME=(IDAY-1)*24+IHOUR                                                   
      IF (ITIME.GT.1) GO TO 123                                                 
C     FIRST PERIOD OF INSTANTANEOUS CHANGE.                                     
      CHANGE=CHANGE+ABS(TSMAX)                                                  
  123 IF (CHANGE.LE.DTONE)RETURN                                                
C     CHANGE LARGE ENOUGH TO SUB-DIVIDE.                                        
C     COMPUTE NUMBER OF INCREMENTS.                                             
      NINC=((CHANGE-DTONE)**TEXP)+2.01                                          
  131 TINC=1.0/NINC                                                             
C     ADJUST FIRST GUESS FOR TIME INCREMENT.                                    
      DO 122 I=1,N                                                              
      TTDT(I)=TT(I)+(TTDT(I)-TT(I))*TINC                                        
      WTDT(I)=WT(I)+(WTDT(I)-WT(I))*TINC                                        
      IF ((TT(I).EQ.273.16).AND.(WT(I).EQ.0.0)) GO TO 122                       
      IF(WT(I).GT.0.0) GO TO 125                                                
C     TEMPERATURE WAS UNKNOWN AT START OF PERIOD.                               
      IF (NOKNOW(I).EQ.0) GO TO 122                                             
C     LIQUID-WATER NOW UNKNOWN--CHANGE BACK TO TEMPERATURE.                     
      WTDT(I)=0.0                                                               
      NOKNOW(I)=0                                                               
      GO TO 122                                                                 
C     LIQUID-WATER WAS UNKNOWN AT START OF PERIOD.                              
  125 IF (NOKNOW(I).EQ.1) GO TO 122                                             
C     TEMPERATURE NOW UNKNOWN--CHANGE BACK TO LIQUID-WATER.                     
      TTDT(I)=273.16                                                            
      NOKNOW(I)=1                                                               
  122 CONTINUE                                                                  
C*******************************************************************************
      RETURN                                                                    
      END                                                                       

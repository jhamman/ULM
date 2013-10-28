      SUBROUTINE SNOWTW(N,TT,TTDT,WT,WTDT,D,P,ITOPT,QI,QR,QA,FU,TA,EA,          
     1PX,D_O,PA,X,DTG,TGT,TGTDT,TCG,DCG,THEDA,NOKNOW,TOLER,DELTAT,ITER,          
     2MONTH,IDAY,IYEAR,IHOUR,NINC,QG,ITS,EXTINC,CV,G1,G2,G3,
     3IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE)             
C*******************************************************************************
C     THIS SUBROUTINE COMPUTES SNOW COVER TEMPERATURES AND AMOUNT               
C        OF LIQUID-WATER FOR TIME T + DT. THE SNOW COVER MUST CONSIST           
C        OF TWO OR MORE LAYERS.                                                 
C*******************************************************************************
C     ITEMS IN THE ARGUMENT LIST ARE.                                           
C             1. NUMBER OF LAYERS IN THE SNOW COVER (N) INCLUDING               
C                SURFACE LAYER                                                  
C             2. TEMPERATURE OF EACH LAYER AT TIME T IN DEG. K (TT)             
C             3. TEMPERATURE OF EACH LAYER AT TIME T+DT IN DEG. K (TTDT)        
C                  INPUT (INITIAL GUESS)                                        
C                  OUTPUT (COMPUTED VALUE FROM ENERGY BALANCE)                  
C             4. LIQUID-WATER CONTENT OF EACH LAYER AT TIME T IN CM.(WT)        
C             5. LIQUID-WATER CONTENT OF EACH LAYER AT TIME T+DT                
C                IN CM. (WTDT)                                                  
C                  INPUT (INITIAL GUESS)                                        
C                  OUTPUT (COMPUTED VALUE FROM ENERGY BALANCE)                  
C             6. DEPTH OF EACH LAYER IN CM. (D)                                 
C             7. DENSITY OF THE ICE PORTION OF EACH LAYER, GM/CM3(P)            
C             ***NOTE*** THE FIRST LAYER IS THE SURFACE LAYER IN THE            
C             COMPUTER PROGRAM. THIS IS LAYER ZERO IN THE TEXT.                 
C             8. ITOPT =0,COMPUTE SNOW SURFACE TEMPERATURE                      
C                         INPUT TTDT(1) IS INITIAL GUESS                        
C                         OUTPUT TTDT(1) IS COMPUTED VALUE                      
C                      =1, DO NOT COMPUTE SURFACE TEMPERATURE,USE               
C                          INPUT VALUES. COMPUTE LIQUID WATER WHEN              
C                          TTDT(1) EQUALS ZERO DEGREES CELSIUS.                 
C             9. INCOMING SOLAR RADIATION IN CAL/CM2(QI)                        
C            10. REFLECTED SOLAR RADIATION IN CAL/CM2(QR)                       
C            11. ATMOSPHERIC LONGWAVE RADIATION IN CAL*CM-2(QA)                 
C            12. WIND FUNCTION IN CM./MB.(FU). COMPUTE FROM WIND                
C                SPEED IN SUBROUTINE WINDF.                                     
C            13. AIR TEMPERATURE IN DEG. K (TA)                                 
C            14. VAPOR PRESSURE OF THE AIR IN MILLIBARS (EA)                    
C            15. AMOUNT OF RAIN IN CM. (PX)                                     
C            ***NOTE*** ITEMS 9-15 ARE PERIOD VALUES.                           
C               9,10,11,15 ARE TOTALS -- 12,13,14 ARE AVERAGE VALUES.           
C            16. VAPOR DIFFUSION COEFFICIENT AT 1000 MB PRESSURE                
C                AND TEMPERATURE=273.16 DEG. K.(DO)CM2/SEC                      
C                IF (DO=0.0) THEN HEAT TRANSFER DUE TO VAPOR DIFFUSION          
C                   IS NEGLECTED.                                               
C            17. MEAN STATION ATMOSPHERIC PRESSURE IN MILLIBARS (PA)            
C            18. EXPONENT IN DIFFUSION COEF. EQUATION (X)                       
C            19. DEPTH INTO GROUND OF SOIL TEMPERATURE,CM(DTG)                  
C            20. SOIL TEMPERATURE AT TIME T IN DEG. K (TGT)                     
C            21. SOIL TEMPERATURE AT TIME T+DT IN DEG. K (TGTDT)                
C                THIS IS A KNOWN VALUE, IT IS NOT COMPUTED.                     
C            22. THERMAL CONDUCTIVITY OF THE SOIL IN                            
C                CAL/CM/SEC/DEG.K (TCG)                                         
C            23. VAPOR DIFFUSION COEFFICIENT FOR SOIL CM2/SEC (DCG)             
C            24. (THEDA)-WEIGHTING FACTOR THAT WEIGHTS THE SPACIAL              
C                DERIVATIVES WITH RESPECT TO TIME.                              
C            25. ARRAY INDICATING THE CURRENT UNKNOWN FOR EACH SNOW             
C                LAYER (NOKNOW). IF ELEMENT=0,TEMPERATURE IS                    
C                UNKNOWN--IF=1,LIQUID WATER IS UNKNOWN.                         
C            26. ITERATIVE SOLUTION TOLERANCE FOR TEMPERATURE IN DEG.K          
C                (TOLER)LIQUID-WATER TOLERANCE COMPUTED FROM THIS VALUE.        
C            27. TIME INTERVAL IN HOURS (DELTAT)                                
C            28. OUTPUT,THE NUMBER OF ITERATIONS THAT WERE MADE (ITER)          
C            29. MONTH NUMBER (MONTH)                                           
C            30. DAY (IDAY)                                                     
C            31. YEAR (IYEAR) - 4 DIGITS                                        
C            32. LAST HOUR OF THE PERIOD (IHOUR)                                
C            33. NUMBER OF TIME INCREMENTS TO USE FOR THE                       
C                TIME PERIOD.(NINC)                                             
C*******************************************************************************
C     NOTE.....DURING EACH CALL TO THIS SUBROUTINE THE NUMBER                   
C          OF LAYERS REMAINS FIXED. THE INITIAL DEPTHS AND                      
C          DENSITIES ARE ASSUMED TO BE APPLICABLE THROUGHOUT                    
C          THE PERIOD. HOWEVER,DURING SOME PERIODS THE STATE                    
C          OF SOME LAYERS MAY CHANGE DRASTICALLY, THE LAYERS                    
C          MAY EVEN DISAPPEAR. A REDEFINITION OF THE STATE OF                   
C          EACH LAYER WILL TAKE PLACE AFTER THIS SUBROUTINE IS                  
C          EXECUTED. SOME LAYERS MAY CONTAIN MORE LIQUID-                       
C          WATER THAN THE INITIAL WATER-EQUIVALENT OF THE LAYER.                
C          THIS WILL BE TAKEN CARE OF LATER BY REDUCING THE                     
C          ICE CONTENT OF SUBSEQUENT LAYERS TO MATCH THE LIQUID                 
C          EXCESS. THIS PROCEDURE IS REASONABLE AND MUCH MORE                   
C          COMPUTATIONLY EFFICIENT THAN REDEFINITION DURING                     
C          ITERATIVE LOOPS.                                                     
C*******************************************************************************
C     SUBROUTINE IS CURRENTLY DIMENSIONED FOR 100 LAYERS.                       
      REAL LS,LF                                                                
      COMMON/DEBUG/IDEBUG,IDBUG,IHR,LHR,LISTF                                   
C      COMMON/IFU,UMSEC,FUCOEF,HEIGHT                                     
      COMMON/STEMP/TOSIM                                                        
      DIMENSION TT(100),TTDT(100),WT(100),WTDT(100),D(100),P(100),              
     1NOKNOW(100),Z(100),AQI(100),TC(100),DTCDP(100),DPDZ(100),F1T(100),        
     2F2T(100),F3T(100),F4T(100),F5T(100),F6T(100),F7T(100),                    
     3DL1(100),DL2(100),F1TDT(100),F2TDT(100),F3TDT(100),F4TDT(100),            
     4F5TDT(100),F6TDT(100),F7TDT(100),CIT(100),CITDT(100),                     
     5DEDU(100,100),R(100),DU(100),TTI(100),WTI(100)                            
      DL4=D(N)*0.5                                                              
      DL5=D(N-1)*0.5                                                            
      TOSIM=0.0                                                                 
      TINC=1.0/NINC                                                             
      IF (NINC.EQ.1) GO TO 170                                                  
C     STORE TIME T VALUES.                                                      
      DO 171 I=1,N                                                              
      TTI(I)=TT(I)                                                              
      WTI(I)=WT(I)                                                              
  171 CONTINUE                                                                  
      FUI=FU                                                                    
      PXI=PX                                                                    
      NITER=0                                                                   
      FU=FU*TINC                                                                
      PX=PX*TINC                                                                
      QAI=QA                                                                    
C     DT IS DELTA TIME IN SECONDS FOR EACH TIME INCREMENT.                      
C  170 DT=3600.0*DELTAT*TINC                
C Change deltaT into seconds
  170 DT=DELTAT*TINC 
      INC=1                                                                     
C*******************************************************************************
C     OBTAIN OTHER PERIOD CONSTANTS                                             
      CALL VSTART(QIR,QI,QR,QA,SIGMA,DT,CW,LS,PW,PA,GAMMA,                      
     1THEDA,THEDAM,C1,LF,C2,D_O,X,C3,C4,E,RW,ITS)                             
      QA=QA*TINC                                                                
      QIR=QIR*TINC                                                              
      CALL ZDEPTH(N,D,Z,TDEPTH)                                                 
      CALL ABSORB(N,QI,QR,Z,D,P,AQI,EXTINC,CV,G1,G2,G3)                          
      DO 172 I=1,N                                                              
  172 AQI(I)=AQI(I)*TINC                                                        
      CALL COEFF(N,P,TC,DTCDP,COEFKE)                                             
      CALL CDPDZ(N,P,Z,DPDZ)                                                    
      TW=0.0                                                                    
      IF (PX.GT.0.0) CALL WTBULB(TW,TA,EA,PA)                                   
      IF ((PX.GT.0.0).AND.(TW.LT.273.16)) TW=273.16                             
C     OBTAIN OTHER VALUES NEEDED AT TIME T.                                     
  175 TT4=TT(1)*TT(1)*TT(1)*TT(1)                                               
      EOT=C1*EXP(-6141.9/TT(1))                                                 
      CALL VALUES(N,D_O,PA,X,C1,TC,TT,Z,F1T,F2T,F3T,F4T,                         
     1F5T,F6T,F7T,DL1,DL2)                                                      
      CALL SPHEAT(N,TT,C3,C4,CIT)                                               
      CALL GROUND(TGT,C1,TCG,LS,DCG,FGT,DCDTG)                                  
C*******************************************************************************
C     BEGIN ITERATIVE LOOP TO DETERMINE TEMPERATURE AND                         
C        LIQUID-WATER AT T+DT.                                                  
      MAX=10                                                                    
      IF (N.EQ.2) MAX=20                                                        
      DO 100 LOOP=1,MAX                                                         
      ITER=LOOP-1                                                               
      I=1                                                                       
C     OBTAIN VALUES FOR TIME T+DT. DL1 AND DL2 DO NOT CHANGE.                   
      CALL VALUES(N,D_O,PA,X,C1,TC,TTDT,Z,F1TDT,F2TDT,F3TDT,                     
     1F4TDT,F5TDT,F6TDT,F7TDT,DL1,DL2)                                          
      CALL SPHEAT(N,TTDT,C3,C4,CITDT)                                           
      CALL GROUND(TGTDT,C1,TCG,LS,DCG,FGTDT,DCDTG)                              
      IF (IFU.GT.0) CALL WINDF(DELTAT,UMSEC,TA,TT(1),TTDT(1),FUCOEF,PA,         
     1HEIGHT,IFU,FU,RICRIT,ZO)                                                  
C*******************************************************************************
C     COMPUTE RESIDUALS FOR EACH EQUATION,STORE NEGATIVE RESIDUAL.              
C*******************************************************************************
C     SURFACE LAYER EQUATION                                                    
      TTDT2=TTDT(1)*TTDT(1)                                                     
      TTDT3=TTDT2*TTDT(1)                                                       
      TTDT4=TTDT3*TTDT(1)                                                       
      EOTDT=C1*EXP(-6141.9/TTDT(1))                                             
      IF (IDEBUG.EQ.0) GO TO 988                                                
C***************************************                                        
C     DEBUG STATEMENTS.                                                         
C Remove time controls
c      IF (IDAY.NE.IDBUG) GO TO 988                                              
c      IF ((IHOUR.GE.IHR).AND.(IHOUR.LE.LHR)) GO TO 987                          
      GO TO 988                                                                 
  987 TO4=E*SIGMA*TT4                                                           
      TO4DT=E*SIGMA*TTDT4                                                       
      PRINT 990,MONTH,IDAY,IYEAR,IHOUR,LOOP,NINC,INC,DT                         
  990 FORMAT (1H1,I2,1H/,I2,1H/,I4,2X,5HHOUR=,I2,5X,5HLOOP=,I1,5X,              
     15HNINC=,I2,5X,4HINC=,I2,5X,3HDT=,F6.0)                                    
      PRINT 991,QIR,TO4,TO4DT,FU,TA,EA,EOT,EOTDT,PX,TW,TGT,TGTDT,FGT,           
     1FGTDT                                                                     
  991 FORMAT (1H ,14F9.4)                                                       
      PRINT 991,(TT(L),L=1,N)                                                   
      PRINT 991,(TTDT(L),L=1,N)                                                 
      PRINT 991,(WT(L),L=1,N)                                                   
      PRINT 991,(WTDT(L),L=1,N)                                                 
      PRINT 992,(NOKNOW(L),L=1,N)                                               
  992 FORMAT (1H ,14I9)                                                         
      PRINT 991,(CIT(L),L=1,N)                                                  
      PRINT 991,(CITDT(L),L=1,N)                                                
      PRINT 991,(D(L),L=1,N)                                                    
      PRINT 991,(Z(L),L=1,N)                                                    
      PRINT 991,(P(L),L=1,N)                                                    
      PRINT 991,(AQI(L),L=1,N)                                                  
      PRINT 991,(TC(L),L=1,N)                                                   
      IF (LISTF.EQ.0) GO TO 986                                                 
      PRINT 991,(DTCDP(L),L=1,N)                                                
      PRINT 991,(DPDZ(L),L=1,N)                                                 
      PRINT 991,(DL1(L),L=1,N)                                                  
      PRINT 991,(DL2(L),L=1,N)                                                  
      PRINT 991,(F1T(L),L=1,N)                                                  
      PRINT 991,(F1TDT(L),L=1,N)                                                
      PRINT 994,(F2T(L),L=1,N)                                                  
  994 FORMAT (1H ,14E9.3)                                                       
      PRINT 994,(F2TDT(L),L=1,N)                                                
      PRINT 994,(F3T(L),L=1,N)                                                  
      PRINT 994,(F3TDT(L),L=1,N)                                                
      PRINT 994,(F4T(L),L=1,N)                                                  
      PRINT 994,(F4TDT(L),L=1,N)                                                
      PRINT 994,(F5T(L),L=1,N)                                                  
      PRINT 994,(F5TDT(L),L=1,N)                                                
      PRINT 991,(F6T(L),L=1,N)                                                  
      PRINT 991,(F6TDT(L),L=1,N)                                                
      PRINT 991,(F7T(L),L=1,N)                                                  
      PRINT 991,(F7TDT(L),L=1,N)                                                
  986 PRINT 993                                                                 
  993 FORMAT (1H )                                                              
C     END DEBUG.                                                                
C***************************************                                        
  988 CONTINUE                                                                  
      IF ((TT(N).GT.273.159).AND.(TTDT(N).GT.273.159)) GO TO 101                
C     QG WHEN BOTTON LAYER IS NOT AT ZERO DEGREES CELSIUS.                      
      T1=(THEDAM*DT*F1T(N)*FGT*(TGT-TT(N)))/(FGT*DL4+F1T(N)*DTG)                
      T2=(THEDA*DT*F1TDT(N)*FGTDT*(TGTDT-TTDT(N)))/                             
     1(FGTDT*DL4+F1TDT(N)*DTG)                                                  
      QG=T1+T2                                                                  
      GO TO 102                                                                 
C     QG WHEN BOTTON LAYER IS AT ZERO DEGREES CELSIUS                           
  101 T3=THEDAM*DT*FGT*(TGT-TT(N))/DTG                                          
      T4=THEDA*DT*FGTDT*(TGTDT-TTDT(N))/DTG                                     
      QG=T3+T4                                                                  
  102 IF(ITOPT.EQ.1) GO TO 108                                                  
C     COMPUTE CHANGE IN HEAT STORAGE FOR THE REST OF THE SNOW COVER.            
      DQ=0.0                                                                    
      DO 103 J=2,N                                                              
      DQ=DQ+D(J)*P(J)*CIT(J)*TT(J)-D(J)*P(J)*CITDT(J)*TTDT(J)                   
     1+LF*PW*WT(J)-LF*PW*WTDT(J)                                                
  103 CONTINUE                                                                  
      DH=QIR-THEDAM*E*SIGMA*TT4-THEDA*E*SIGMA*TTDT4                             
     1+LS*PW*FU*GAMMA*TA-THEDAM*LS*PW*FU*GAMMA*TT(1)                            
     2-THEDA*LS*PW*FU*GAMMA*TTDT(1)+LS*PW*FU*EA                                 
     3-THEDAM*LS*PW*FU*EOT-THEDA*LS*PW*FU*EOTDT+                                
     4CW*PW*PX*TW-CW*PW*PX*273.16+QG+DQ                                         
      RESID=D(1)*P(1)*CITDT(1)*TTDT(1)-D(1)*P(1)*CIT(1)*TT(1)                   
     1+LF*PW*WTDT(1)-LF*PW*WT(1)-DH                                             
      R(I)=-RESID                                                               
      I=I+1                                                                     
  108 IF (N.EQ.2) GO TO 105                                                     
C*******************************************************************************
C     INTERMEDIATE LAYER EQUATIONS.                                             
      NN=N-1                                                                    
      DO 104 J=2,NN                                                             
      P1=(THEDAM*AQI(J))/(P(J)*D(J)*CIT(J))+(THEDA*AQI(J))/                     
     1(P(J)*D(J)*CITDT(J))+(THEDAM*2.0*DT*F1T(J)*F7T(J))/                       
     2(P(J)*CIT(J)*DL1(J)*DL2(J)*(DL1(J)+DL2(J)))+(THEDA*2.0*DT*                
     3F1TDT(J)*F7TDT(J))/(P(J)*CITDT(J)*DL1(J)*DL2(J)*                          
     4(DL1(J)+DL2(J)))+(THEDAM*DT*DTCDP(J)*DPDZ(J)*F6T(J))/                     
     5(2.0*P(J)*CIT(J)*DL1(J)*DL2(J))+(THEDA*DT*DTCDP(J)*                       
     6DPDZ(J)*F6TDT(J))/(2.0*P(J)*CITDT(J)*DL1(J)*DL2(J))                       
      P2=(THEDAM*DT*C2*LS*(((X/TT(J))*F2T(J))+F4T(J))*F6T(J)                    
     1*F6T(J))/(4.0*P(J)*CIT(J)*DL1(J)*DL1(J)*DL2(J)*DL2(J))+                   
     2(THEDA*DT*C2*LS*(((X/TTDT(J))*F2TDT(J))+F4TDT(J))                         
     3*F6TDT(J)*F6TDT(J))/(4.0*P(J)*CITDT(J)*DL1(J)*DL1(J)*                     
     4DL2(J)*DL2(J))                                                            
      DTN=P1+P2                                                                 
      RESID=TTDT(J)-TT(J)+(LF*WTDT(J)*PW)/(D(J)*P(J)*CITDT(J))                  
     1-(LF*WT(J)*PW)/(D(J)*P(J)*CIT(J))-DTN                                     
      R(I)=-RESID                                                               
      I=I+1                                                                     
  104 CONTINUE                                                                  
C*******************************************************************************
C     BOTTOM LAYER EQUATION.                                                    
  105 IF ((TT(N).GT.273.159).AND.(TTDT(N).GT.273.159)) GO TO 106                
C     BOTTOM LAYER NOT AT ZERO DEGREES CELSIUS.                                 
      DGT=T1/(P(N)*D(N)*CIT(N))+T2/(P(N)*D(N)*CITDT(N))                         
      GO TO 107                                                                 
C     BOTTOM LAYER AT ZERO DEGREES CELSIUS.                                     
  106 DGT=T3/(P(N)*D(N)*CIT(N))+T4/(P(N)*D(N)*CITDT(N))                         
  107 DN=(THEDAM*AQI(N))/(P(N)*D(N)*CIT(N))+(THEDA*AQI(N))/                     
     1(P(N)*D(N)*CITDT(N))+DGT+(THEDAM*DT*F1T(N)                                
     2*F1T(N-1)*(TT(N-1)-TT(N)))/(P(N)*D(N)*CIT(N)*                             
     3(F1T(N)*DL5+F1T(N-1)*DL4))+(THEDA*DT*F1TDT(N)*F1TDT(N-1)*                 
     4(TTDT(N-1)-TTDT(N)))/(P(N)*D(N)*CITDT(N)*                                 
     5(F1TDT(N)*DL5+F1TDT(N-1)*DL4))                                            
      RESID=TTDT(N)-TT(N)+(LF*WTDT(N)*PW)/                                      
     1(D(N)*P(N)*CITDT(N))-(LF*WT(N)*PW)/(D(N)*P(N)*CIT(N))-DN                  
      R(I)=-RESID                                                               
C*******************************************************************************
C     INITIALIZE PARTIAL DERIATIVE MATRIX TO ZERO.                              
      NUM=I                                                                     
C     NUM IS THE NUMBER OF CORRECTIONS THAT NEED TO BE COMPUTED.                
      DO 109 I=1,NUM                                                            
      DO 109 J=1,NUM                                                            
  109 DEDU(I,J)=0.0                                                             
C*******************************************************************************
C     COMPUTE DEDU VALUES FOR EACH EQUATION.                                    
      DO 110 K=1,NUM                                                            
C     K IS THE ROW NUMBER.                                                      
      IF (ITOPT.EQ.1) GO TO 117                                                 
      IF (K.GT.1) GO TO 117                                                     
C*******************************************************************************
C     COMPUTE PARTIALS FOR THE SURFACE LAYER                                    
      DO 111 J=1,NUM                                                            
      IF (NOKNOW(J).EQ.1) GO TO 112                                             
C     TEMPERATURE IS THE UNKNOWN.                                               
      IF (J.GT.1) GO TO 113                                                     
C     SURFACE LAYER                                                             
      DEODTO=D(1)*P(1)*C3+2.0*D(1)*P(1)*C4*TTDT(1)+                             
     14.0*THEDA*E*SIGMA*TTDT3+THEDA*LS*PW*FU*GAMMA                              
     2+(THEDA*LS*LS*PW*FU*EOTDT)/(RW*TTDT2)                                     
      DEDU(K,J)=DEODTO                                                          
      GO TO 111                                                                 
  113 IF (J.EQ.N) GO TO 114                                                     
C     INTERMEDIATE LAYERS                                                       
      DEODT=D(J)*P(J)*C3+2.0*D(J)*P(J)*C4*TTDT(J)                               
      DEDU(K,J)=DEODT                                                           
      GO TO 111                                                                 
C     BOTTOM LAYER                                                              
  114 IF ((TT(N).GT.273.159).AND.(TTDT(N).GT.273.159)) GO TO 115                
C     BOTTOM LAYER NOT AT ZERO DEGREES CELSIUS.                                 
      DEODTN=((THEDA*DT*FGTDT)/(FGTDT*DL4+F1TDT(N)*DTG))*                       
     1((DTG*LS*C2*F1TDT(N)*F3TDT(N)*(TGTDT-TTDT(N)))/                           
     2(FGTDT*DL4+F1TDT(N)*DTG)-LS*C2*F3TDT(N)*                                  
     3(TGTDT-TTDT(N))+F1TDT(N))                                                 
      GO TO 116                                                                 
C     BOTTOM LAYER AT ZERO DEGREES CELSIUS.                                     
  115 DEODTN=(THEDA*DT*FGTDT)/DTG                                               
  116 DEDU(K,J)=DEODTN                                                          
      GO TO 111                                                                 
C     LIQUID WATER IS THE UNKNOWN.                                              
  112 DEDU(K,J)=LF*PW                                                           
  111 CONTINUE                                                                  
      GO TO 110                                                                 
  117 IF (N.EQ.2) GO TO 130                                                     
      IF (K.EQ.NUM) GO TO 130                                                   
C*******************************************************************************
C     COMPUTE PARTIALS FOR INTERMEDIATE LAYERS.                                 
      J=K-1                                                                     
      I=K                                                                       
      IF (ITOPT.EQ.1) I=K+1                                                     
      IF (K.EQ.1) GO TO 120                                                     
C     K WILL ONLY EQUAL ONE THE FIRST TIME THROUGH THIS LOOP                    
C        IF ITOPT EQUALS ONE.                                                   
      IF (NOKNOW(I-1).EQ.1) GO TO 120                                           
      DEDTM=(-THEDA*2.0*DT*F1TDT(I))/(P(I)*CITDT(I)*                            
     1DL1(I)*(DL1(I)+DL2(I)))+(THEDA*DT*DTCDP(I)*                               
     2DPDZ(I))/(2.0*P(I)*CITDT(I)*DL1(I))+                                      
     3(THEDA*2.0*DT*C2*LS*F6TDT(I)*((X/TTDT(I))                                 
     4*F2TDT(I)+F4TDT(I)))/(4.0*P(I)*CITDT(I)*                                  
     5DL1(I)*DL1(I)*DL2(I))                                                     
      DEDU(K,J)=DEDTM                                                           
C     COMPUTE PARTIAL WHICH LIES ON THE DIAGONAL                                
  120 J=K                                                                       
      IF (NOKNOW(I).EQ.1) GO TO 121                                             
C     TEMPERATURE UNKNOWN                                                       
      P1  =1.0+(THEDA*C4*AQI(I))/(D(I)*P(I)*CITDT(I)*CITDT(I))+                 
     1((THEDA*2.0*DT)/(P(I)*DL1(I)*DL2(I)*(DL1(I)+DL2(I))*                      
     2CITDT(I)))*((F1TDT(I)*C4*F7TDT(I))/CITDT(I)-F7TDT(I)*LS*                  
     3C2*F3TDT(I)+F1TDT(I)*(DL1(I)+DL2(I)))+((THEDA*DT*DTCDP(I)*                
     4DPDZ(I))/(2.0*P(I)*DL1(I)*DL2(I)*CITDT(I)))*                              
     5((C4*F6TDT(I)/CITDT(I))-(DL2(I)-DL1(I)))                                  
      P2=((THEDA*DT*C2*LS)/                                                     
     1(4.0*P(I)*DL1(I)*DL1(I)*DL2(I)*DL2(I)*CITDT(I)))*                         
     2(((C4*F6TDT(I)*F6TDT(I))*((X/TTDT(I))*F2TDT(I)+F4TDT(I)))/                
     3CITDT(I)-(F6TDT(I)*F6TDT(I))*(F5TDT(I)+((X/TTDT(I))*                      
     4F3TDT(I))-((X*F2TDT(I))/(TTDT(I)*TTDT(I))))-2.0*F6TDT(I)*                 
     5(DL2(I)-DL1(I))*(((X/TTDT(I))*F2TDT(I))+F4TDT(I)))                        
      DEDT=P1+P2                                                                
      DEDU(K,J)=DEDT                                                            
      GO TO 125                                                                 
C     LIQUID-WATER UNKNOWN.                                                     
  121 DEDU(K,J)=(LF*PW)/(D(I)*P(I)*CITDT(I))                                    
  125 J=K+1                                                                     
      IF (NOKNOW(I+1).EQ.1) GO TO 110                                           
      DEDTP=((-THEDA*2.0*DT*F1TDT(I))/(P(I)*CITDT(I)*DL2(I)*                    
     1(DL1(I)+DL2(I))))-((THEDA*DT*DTCDP(I)*DPDZ(I))/                           
     2(2.0*P(I)*CITDT(I)*DL2(I)))-(((THEDA*2.0*DT*C2*LS*                        
     3F6TDT(I))/(4.0*P(I)*CITDT(I)*DL1(I)*DL2(I)*DL2(I)))                       
     4*(((X/TTDT(I))*F2TDT(I))+F4TDT(I)))                                       
      DEDU(K,J)=DEDTP                                                           
      GO TO 110                                                                 
C     END OF PARTIALS FOR INTERMEDIATE LAYERS.                                  
C*******************************************************************************
C     COMPUTE PARTIALS FOR THE BOTTOM LAYER                                     
  130 J=K-1                                                                     
      IF ((N.EQ.2).AND.(ITOPT.EQ.1)) GO TO 131                                  
      IF (NOKNOW(N-1).EQ.1) GO TO 131                                           
      DEDTNM=((THEDA*DT*F1TDT(N))/(P(N)*D(N)*CITDT(N)*                          
     1(F1TDT(N)*DL5+F1TDT(N-1)*DL4)))*(((F1TDT(N-1)*DL4*LS*C2*                  
     2F3TDT(N-1)*(TTDT(N-1)-TTDT(N)))/(F1TDT(N)*DL5+F1TDT(N-1)*DL4))            
     3-(LS*C2*F3TDT(N-1)*(TTDT(N-1)-TTDT(N)))-F1TDT(N-1))                       
      DEDU(K,J)=DEDTNM                                                          
  131 J=K                                                                       
      IF (NOKNOW(N).EQ.1) GO TO 132                                             
C     TEMPERATURE UNKNOWN.                                                      
      IF ((TT(N).GT.273.159).AND.(TTDT(N).GT.273.159)) GO TO 133                
C     BOTTOM LAYER NOT AT ZERO DEGREES CELSIUS.                                 
      DTGDTN=((THEDA*DT*FGTDT)/(P(N)*D(N)*CITDT(N)*(F1TDT(N)*                   
     1DTG+FGTDT*DL4)))*(((F1TDT(N)*DTG*LS*C2*F3TDT(N)*                          
     2(TGTDT-TTDT(N)))/(F1TDT(N)*DTG+FGTDT*DL4))+((F1TDT(N)*C4*                 
     3(TGTDT-TTDT(N)))/CITDT(N))-LS*C2*F3TDT(N)*                                
     4(TGTDT-TTDT(N))+F1TDT(N))                                                 
      GO TO 134                                                                 
C     BOTTOM LAYER AT ZERO DEGREES CELSIUS.                                     
  133 DTGDTN=((THEDA*DT*FGTDT)/(P(N)*D(N)*DTG*CITDT(N)))*                       
     1(((C4*(TGTDT-TTDT(N)))/CITDT(N))+1.0)                                     
  134 DENDTN=1.0+(THEDA*C4*AQI(N))/(P(N)*D(N)*CITDT(N)*CITDT(N))                
     1+((THEDA*DT*F1TDT(N-1))/(P(N)*D(N)*CITDT(N)*(F1TDT(N)                     
     2*DL5+F1TDT(N-1)*DL4)))*(((F1TDT(N)*DL5*LS*C2*F3TDT(N)*                    
     3(TTDT(N-1)-TTDT(N)))/(F1TDT(N)*DL5+F1TDT(N-1)*DL4))+                      
     4(F1TDT(N)*C4*(TTDT(N-1)-TTDT(N))/CITDT(N))-                               
     5(LS*C2*F3TDT(N)*(TTDT(N-1)-TTDT(N)))+F1TDT(N))+DTGDTN                     
      DEDU(K,J)=DENDTN                                                          
      GO TO 110                                                                 
C     LIQUID-WATER UNKNOWN                                                      
  132 DEDU(K,J)=(LF*PW)/(D(N)*P(N)*CITDT(N))                                    
  110 CONTINUE                                                                  
C     PARTIALS HAVE BEEN COMPUTED FOR EACH EQUATION AND EACH UNKNOWN.           
C*******************************************************************************
C     DETERMINE THE CORRECTIONS FOR EACH LAYER.                                 
      IF (IDEBUG.EQ.0) GO TO 985                                                
C***************************************                                        
C     DEBUG STATEMENTS.                                                         
C      IF (IDAY.NE.IDBUG) GO TO 985                                              
C      IF ((IHOUR.LT.IHR).OR.(IHOUR.GT.LHR)) GO TO 985                           
C Remove time controls
      IF (LISTF.EQ.0) GO TO 983                                                 
      DO 989 LL=1,N                                                             
  989 PRINT 991,(DEDU(LL,L),L=1,N)                                              
  983 PRINT 991,(R(L),L=1,N)                                                    
C     END DEBUG.                                                                
C***************************************                                        
  985 CONTINUE                                                                  
      CALL MATRIX(DEDU,R,DU,NUM)                                                
      IF (IDEBUG.EQ.0) GO TO 984                                                
C***************************************                                        
C     DEBUG STATEMENTS.                                                         
C      IF (IDAY.NE.IDBUG) GO TO 984                                              
C      IF ((IHOUR.LT.IHR).OR.(IHOUR.GT.LHR)) GO TO 984                           
C Remove time controls
      PRINT 991,(DU(L),L=1,N)                                                   
C     END DEBUG.                                                                
C***************************************                                        
  984 CONTINUE                                                                  
C*******************************************************************************
C     ADD CORRECTIONS AND CHECK TOLERANCES.                                     
      IF (ITOPT.NE.1) GO TO 139                                                 
      IF (LOOP.GT.1) GO TO 139                                                  
C     COMPUTE CHANGE IN LIQUID WATER IF ANY FOR ITOPT.EQ.1                      
      IF (TTDT(1).LT.273.159) GO TO 139                                         
      TAVG=(TT(1)+TTDT(1))*0.5                                                  
      EO=C1*EXP(-6141.9/TAVG)                                                   
      HEAT=AQI(1)+E*QA-E*SIGMA*TAVG*TAVG*TAVG*TAVG+LS*FU*                       
     1(GAMMA*(TA-TAVG)+(EA-EO))+CW*PW*PX*(TW-273.16)                            
C     ASSUME IF THE SNOW SURFACE TEMPERATURE AT THE END                         
C     OF THE PERIOD IS SET TO ZERO DEGREES CELSIUS THAT                         
C     MELT MUST OCCUR. IF ENERGY EXCHANGE IS NEGATIVE,                          
C     ASSUME MELT IS ZERO.                                                      
      IF (HEAT.LT.0.0) HEAT=0.0                                                 
      DELTAW=HEAT/(LF*PW)                                                       
  139 NOVER=0                                                                   
C     NOVER IS THE NUMBER OF LAYERS WHERE THE CORRECTION                        
C        EXCEEDS THE TOLERANCE.                                                 
      DO 138 K=1,NUM                                                            
      ABSDU=ABS(DU(K))                                                          
      IF (ABSDU.LE.30.0) GO TO 138                                              
      DU(K)=(DU(K)/ABSDU)*30.0                                                  
  138 CONTINUE                                                                  
      DO 140 K=1,NUM                                                            
      IF (K.GT.1) GO TO 150                                                     
      IF (ITOPT.EQ.1) GO TO 141                                                 
C*******************************************************************************
C     SURFACE LAYER                                                             
      IF (NOKNOW(1).EQ.1) GO TO 142                                             
C     TEMPERATURE UNKNOWN                                                       
      ABSDU=ABS(DU(K))                                                          
      IF (ABSDU.GT.TOLER)NOVER=NOVER+1                                          
      TTDT(1)=TTDT(1)+DU(K)                                                     
      IF (TTDT(1).LE.273.16) GO TO 140                                          
C     TEMPERATURE EXCEEDS ZERO DEGREES CELSIUS--LIQUID WATER NOW UNKNOWN        
      EXCESS=TTDT(1)-273.16                                                     
      TTDT(1)=273.16                                                            
      NOKNOW(1)=1                                                               
      WTDT(1)=((D(1)*P(1)*CITDT(1))/(LF*PW))*EXCESS                             
      GO TO 140                                                                 
C     LIQUID-WATER UNKNOWN                                                      
  142 WTOL=((D(1)*P(1)*CITDT(1))/(LF*PW))*TOLER                                 
      IF (WTOL.LT.0.0001) WTOL=0.0001                                           
      ABSDU=ABS(DU(K))                                                          
      IF (ABSDU.GT.WTOL)NOVER=NOVER+1                                           
      WTDT(1)=WTDT(1)+DU(K)                                                     
      IF (WTDT(1).GE.0.0) GO TO 140                                             
C     LIQUID-WATER NEGATIVE--TEMPERATURE NOW UNKNOWN.                           
      EXCESS=WTDT(1)                                                            
      WTDT(1)=0.0                                                               
      NOKNOW(1)=0                                                               
      TTDT(1)=273.16                                                            
      GO TO 140                                                                 
C     SURFACE LAYER--ITOPT.EQ.1 (ONLY FOR LOOP.EQ.1)                            
  141 IF (LOOP.GT.1) GO TO 150                                                  
      IF (TTDT(1).LT.273.159) GO TO 145                                         
      WTDT(1)=WTDT(1)+DELTAW                                                    
      GO TO 150                                                                 
C     IF SURFACE TEMPERATURE GOES BELOW ZERO DEGREES CELSIUS,                   
C        SET LIQUID-WATER TO ZERO.                                              
  145 WTDT(1)=0.0                                                               
C*******************************************************************************
C     REMAINING SNOW COVER LAYERS                                               
  150 I=K                                                                       
      IF (ITOPT.EQ.1)I=K+1                                                      
C     I IS THE LAYER NUMBER.                                                    
      IF (NOKNOW(I).EQ.1) GO TO 151                                             
C     TEMPERATURE UNKNOWN                                                       
      ABSDU=ABS(DU(K))                                                          
      IF (ABSDU.GT.TOLER) NOVER=NOVER+1                                         
      TTDT(I)=TTDT(I)+DU(K)                                                     
      IF (TTDT(I).LE.273.16) GO TO 140                                          
C     TEMPERATURE EXCEEDS ZERO DEGREES CELSIUS--LIQUID-WATER                    
C        NOW UNKNOWN.                                                           
      EXCESS=TTDT(I)-273.16                                                     
      TTDT(I)=273.16                                                            
      NOKNOW(I)=1                                                               
      WTDT(I)=((D(I)*P(I)*CITDT(I))/(LF*PW))*EXCESS                             
      GO TO 140                                                                 
C     LIQUID-WATER UNKNOWN                                                      
  151 WTOL=((D(I)*P(I)*CITDT(I))/(LF*PW))*TOLER                                 
      IF (WTOL.LT.0.0001) WTOL=0.0001                                           
      ABSDU=ABS(DU(K))                                                          
      IF (ABSDU.GT.WTOL) NOVER=NOVER+1                                          
      WTDT(I)=WTDT(I)+DU(K)                                                     
      IF (WTDT(I).GE.0.0) GO TO 140                                             
C     LIQUID-WATER NEGATIVE--TEMPERATURE NOW UNKNOWN.                           
      EXCESS=WTDT(I)                                                            
      WTDT(I)=0.0                                                               
      NOKNOW(I)=0                                                               
      TTDT(I)=273.16                                                            
  140 CONTINUE                                                                  
C     CORRECTIONS HAVE BEEN MADE AND TOLERANCES CHECKED                         
C*******************************************************************************
C     CHECK FOR TERMINATION                                                     
      IF (NOVER.EQ.0) GO TO 180                                                 
  100 CONTINUE                                                                  
C     END OF ITERATION LOOP                                                     
C*******************************************************************************
C     PRINT MESSAGE THAT SOLUTION DID NOT CONVERGE.                             
      PRINT 900,MONTH,IDAY,IYEAR,IHOUR,TOLER                                    
  900 FORMAT(1H1,28HSOLUTION DID NOT CONVERGE ON,I3,1H/,I2,1H/,I4,              
     12X,5HHOUR=,I2,1H.,5X,12HTOLERANCE IS,F5.2,1X,14HDEGREES KELVIN)           
C     PRINT HEADING                                                             
      PRINT 902                                                                 
  902 FORMAT (1H0,5HLAYER,4X,5HDEPTH,2X,9HTHICKNESS,3X,7HDENSITY,               
     17X,3HAQI,2X,8HCONDUCT.,4X,6HNOKNOW,8X,2HTT,6X,4HTTDT,8X,                  
     22HWT,6X,4HWTDT,2X,8HRESIDUAL,2X,10HCORRECTION)                            
      I=1                                                                       
      K=0                                                                       
C     PRINT SURFACE LAYER.                                                      
      IF (ITOPT.EQ.1) GO TO 160                                                 
      PRINT 901,I,Z(1),D(1),P(1),AQI(1),TC(1),NOKNOW(1),                        
     1TT(1),TTDT(1),WT(1),WTDT(1),R(1),DU(1)                                    
  901 FORMAT (1H ,I5,4F10.3,F10.5,I10,4F10.3,F10.4,F10.5)                       
      GO TO 161                                                                 
  160 PRINT 901,I,Z(1),D(1),P(1),AQI(1),TC(1),NOKNOW(1),TT(1),TTDT(1),          
     1WT(1),WTDT(1)                                                             
      K=1                                                                       
  161 DO 162 I=2,N                                                              
      J=I-K                                                                     
      PRINT 901,I,Z(I),D(I),P(I),AQI(I),TC(I),NOKNOW(I),TT(I),                  
     1TTDT(I),WT(I),WTDT(I),R(J),DU(J)                                          
  162 CONTINUE                                                                  
      TO4=E*SIGMA*TT4                                                           
      TO4DT=E*SIGMA*TTDT4                                                       
      PRINT 903,QIR,FU,TA,EA,EOT,EOTDT,PX,TW,QG,DQ,TO4,TO4DT,DT                 
  903 FORMAT (1H0,4HQIR=,F3.0,2X,3HFU=,F5.4,2X,3HTA=,F5.1,2X,3HEA=,             
     1F4.2,2X,3HEO=,2F5.2,2X,3HPX=,F4.2,2X,3HTW=,F5.1,2X,3HQG=,F5.2,            
     22X,3HDQ=,F5.2,2X,4HTO4=,2F5.1,2X,3HDT=,F5.0)                              
C*******************************************************************************
  180 TOSIM=TOSIM+TINC*0.5*(TTDT(1)+TT(1))                                      
      IF (NINC.EQ.1) RETURN                                                     
      NITER=NITER+ITER                                                          
C     CHECK FOR LAST TIME INCREMENT                                             
      IF (INC.EQ.NINC) GO TO 185                                                
      INC=INC+1                                                                 
C*******************************************************************************
C     COMPUTE FIRST GUESS VALUES FOR TIME T+DT FOR THE NEXT                     
C             TIME INCREMENT. USE LINEAR PROJECTION.                            
      DO 181 I=1,N                                                              
C     WHEN ITOPT=1,A LINEAR VARIATION IN SURFACE TEMPERATURE FROM               
C         TIME T TO T+DT IS USED.                                               
      IF (NOKNOW(I).EQ.1) GO TO 182                                             
C     TEMPERATURE UNKNOWN                                                       
      CHANGE=TTDT(I)-TT(I)                                                      
      TT(I)=TTDT(I)                                                             
      WT(I)=WTDT(I)                                                             
      TTDT(I)=TT(I)+CHANGE                                                      
      IF (TTDT(I).GT.273.16)TTDT(I)=273.16                                      
      GO TO 181                                                                 
C     LIQUID-WATER UNKNOWN.                                                     
  182 CHANGE=WTDT(I)-WT(I)                                                      
      WT(I)=WTDT(I)                                                             
      TT(I)=TTDT(I)                                                             
      WTDT(I)=WT(I)+CHANGE                                                      
      IF (WTDT(I).LT.0.0)WTDT(I)=0.0                                            
  181 CONTINUE                                                                  
      GO TO 175                                                                 
C*******************************************************************************
C     END OF THE TIME PERIOD                                                    
  185 ITER=NITER                                                                
      DO 186 I=1,N                                                              
      TT(I)=TTI(I)                                                              
      WT(I)=WTI(I)                                                              
  186 CONTINUE                                                                  
      FU=FUI                                                                    
      PX=PXI                                                                    
      QA=QAI                                                                    
      RETURN                                                                    
      END                                                                       

      SUBROUTINE ZDEPTH (N,D,Z,TDEPTH)                                          
C*******************************************************************************
C     COMPUTE DEPTH TO THE MID-POINT OF EACH LAYER FROM THE SURFACE.            
C        D= DEPTH OF EACH LAYER                                                 
C        Z= DEPTH FROM SURFACE TO MID-POINT                                     
C        N= NUMBER OF LAYERS, TDEPTH=TOTAL DEPTH IN CM.                         
C*******************************************************************************
      DIMENSION D(100),Z(100)                                                   
      TDEPTH=0.0                                                                
      DO 100 I=1,N                                                              
      Z(I)=TDEPTH+D(I)*0.5                                                      
      TDEPTH=TDEPTH+D(I)                                                        
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE VALUES (N,D_O,PA,X,C1,TC,T,Z,F1,F2,F3,F4,F5,F6,F7,DL1,          
     1DL2)                                                                      
C*******************************************************************************
C     COMPUTE VARIOUS VALUES NEEDED FOR THE FINITE DIFFERENCE SOLUTION          
C     DCDT= PARTIAL OF VAPOR CONCENTRATION WITH RESPECT TO TEMPERATURE.         
C     DCDT2= SECOND PARTIAL, DCDT3=THIRD PARTIAL, T=TEMPERATURE.                
C*******************************************************************************
      DIMENSION T(100),TC(100),Z(100),F1(100),F2(100),F3(100),F4(100),          
     1F5(100),F6(100),F7(100),DL1(100),DL2(100)                                 
      XM= X-1.0                                                                 
      DO 100 I=1,N                                                              
      EI= C1*EXP (-6141.9/T(I))                                                 
C     LS= LATENT HEAT OF SUBLIMATION= 677 CAL/GM.                               
C     RW= GAS CONSTANT FOR WATER VAPOR= 4615.0 MB*CM3/GM*DEG.K.                 
C     LS/RW= 6141.9 DEG.K--EI IS VAPOR PRESSURE OVER ICE IN MB.                 
C     ERG= DYNE*CM= 1.0E-3 MB*CM3= 2.38844E-8 CAL.                              
C        THUS CAL= 4.1868E4 MB*CM3                                              
C        (4.1686E4*677.0)/4615.0= 6141.9 DEG.K.                                 
      AA=6141.9/T(I)                                                            
      TP= EI/(4615.0 * T(I) * T(I))                                             
      DCDT= TP * (AA-1.0)                                                       
      TP= TP/T(I)                                                               
      DCDT2= TP* (AA**2-4.0*AA+2.0)                                             
      TP= TP/T(I)                                                               
      DCDT3= TP* (AA**3-9.0* (AA**2)+18.0*AA-6.0)                               
C     THE F VALUES ARE VARIOUS VALUES NEEDED IN THE SOLUTION                    
      DC= D_O*(1000.0/PA)*((T(I)/273.16)**X)                                     
      F1(I)= TC(I) + 677.0 * DC*DCDT                                            
C     TC= THERMAL CONDUCTIVITY, DC= DIFFUSION COEFFICIENT                       
      TX= T(I)**X                                                               
      TXM= T(I)**XM                                                             
      F2(I)= TX*DCDT                                                            
      F3(I)= TX*DCDT2 + X*DCDT*TXM                                              
C     ONLY F1 AND F3 ARE NEEDED FOR THE NTH LAYER AND SURFACE LAYER.            
      IF (I.EQ.1) GO TO 100                                                     
      IF (I.EQ.N)  GO TO 100                                                    
      F4(I)= TX*DCDT2                                                           
      F5(I)= TX*DCDT3 + X*DCDT2*TXM                                             
C     DL VALUES ARE DEPTH INCREMENTS BETWEEN LAYER MID-POINTS                   
      DL1(I)= Z(I)-Z(I-1)                                                       
      DL2(I)= Z(I+1)-Z(I)                                                       
      F6(I)= DL1(I)* T(I+1)+(DL2(I)-DL1(I)) * T(I)-DL2(I) * T(I-1)              
      F7(I)= DL1(I) * T(I+1)-(DL1(I)+DL2(I)) * T(I) + DL2(I) * T(I-1)           
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WTBULB (TW,TA,EA,PA)                                           
C*******************************************************************************
C     COMPUTES WET-BULB(TW)FROM DRY-BULB(TA) AND VAPOR PRESSURE(EA)             
C*******************************************************************************
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
C     CONVERT WET-BULB TO DEGREES KELVIN.                                       
      TW= TWC+273.16                                                            
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CDPDZ (N,P,Z,DPDZ)                                             
C*******************************************************************************
C     COMPUTES THE FINITE DIFFERENCE APPROXIMATION TO THE PARTIAL               
C        OF DENSITY WITH RESPECT TO DEPTH (DPDZ)                                
C        P=DENSITY,Z=DEPTH--DPDZ NOT NEEDED FOR SURFACE AND                     
C        BOTTOM SNOW LAYERS.                                                    
C*******************************************************************************
      DIMENSION DPDZ(100),P(100),Z(100)                                         
      DO 100 I=1,N                                                              
      IF ((I.EQ.1).OR.(I.EQ.N)) GO TO 100                                       
      DPDZ(I)= 0.5*(((P(I)-P(I-1))/(Z(I)-Z(I-1)))+                              
     1   ((P(I+1)-P(I))/(Z(I+1)-Z(I))))                                         
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COEFF (N,P,TC,DTCDP,COEFKE)                                        
C*******************************************************************************
C     COMPUTES THERMAL CONDUCTIVITY AND THE PARITAL OF CONDUCTIVITY WITH        
C        RESPECT TO DENSITY FOR EACH SNOW LAYER.  DENSITY(P) IS KNOWN.          
C*******************************************************************************
C      COMMON/EFFTC/COEFKE                                                       
      DIMENSION P(100), TC(100), DTCDP(100)                                     
      DO 100 I= 1,N                                                             
C             EXPRESSION FOR EFFECTIVE THERMAL CONDUCTIVITY-CAL/CM/SEC/DEG.K    
      TC(I)=0.00005+COEFKE*P(I)*P(I)                                            
      DTCDP(I)=2.0*COEFKE*P(I)                                                  
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE ABSORB (N,QI,QR,Z,D,P,AQI,EXTINC,CV,G1,G2,G3)                       
C*******************************************************************************
C     COMPUTES THE EXTINCTION COEFFICIENT FOR EACH LAYER PLUS THE               
C        AMOUNT OF SOLAR RADIATION ABSORDED IN EACH LAYER (AQI)                 
C*******************************************************************************
      DIMENSION Z(100), D(100), P(100), AQI(100)                                
C      COMMON/EXTINC/CV,G1,G2,G3                                                 
      TOTAL= QI-QR                                                              
C     TOTAL= ABSORDED RADIATION CAL/CM2                                         
      IF (TOTAL.GT.0.0) GO TO 101                                               
C     NO ABSORBED SOLAR RADIATION                                               
      DO 100 I= 1,N                                                             
  100 AQI(I)= 0.0                                                               
      RETURN                                                                    
  101 SUM= 0.0                                                                  
      DO 102 I= 1,N                                                             
C     EXC IS THE EXTINCTION COEFFICIENT -- TEMPORARY EXPRESSION                 
      DS=G1+G2*P(I)*P(I)+G3*P(I)*P(I)*P(I)*P(I)                                 
      EXC=CV*P(I)*SQRT(1.0/DS)                                                  
      PERCAB=1.0-EXP(-EXC*D(I))                                                 
      AQI(I)=(TOTAL-SUM)*PERCAB                                                 
      SUM= SUM + AQI(I)                                                         
  102 CONTINUE                                                                  
C     RESIDUAL IS ABSORBED IN THE BOTTOM SNOW LAYER                             
      RES= TOTAL-SUM                                                            
      IF (RES.LE.0.0) RETURN                                                    
      AQI(N)= AQI(N) + RES                                                      
      RETURN                                                                    
      END                                                                       
      SUBROUTINE VSTART(QIR,QI,QR,QA,SIGMA,DT,CW,LS,PW,PA,                     
     &   GAMMA,THEDA,THEDAM,C1,LF,C2,D_O,X,C3,C4,E,RW,ITS)                           
C*******************************************************************************
C     ESTABLISHES OR COMPUTES VARIOUS NUMBERS NEEDED TO SOLVE                   
C        THE FINITE DIFFENCE EQUATIONS.  THESE VALUES ARE                       
C        FOR THE ENTIRE PERIOD.                                                 
C*******************************************************************************
C      COMMON/SPCASE/ITS                                                         
      REAL LS,LF         
      E= 0.99                                                                   
      QIR= QI-QR+E*QA                                                           
C     EMISSIVITY OF SNOW ASSUMED EQUAL TO 0.99.                                 
C     SIGMA IS STEFAN-BOLTZMAN CONSTANT FOR THE TIME INTERVAL.                  
C        STEFAN-BOLTZMAN CONSTANT IS 1.355E-12 CAL/CM2/K4/SEC.                  
      SIGMA= 1.355E-12*DT                                                       
      CW=1.0                                                                    
      LS= 677.0                                                                 
      PW= 1.0                                                                   
C     GAMMA IS THE PSYCHROMETRIC CONSTANT IN MB/DEG.K                           
      GAMMA=(0.240*PA)/(0.622*LS)                                               
      THEDAM= 1.0-THEDA      
      C1= 3.5558E10                                                             
C     RW NEEDED IN CAL/GM/DEG K.                                                
      RW=0.110226                                                               
C     LF IS LATENT HEAT OF FUSION                                               
      LF=79.7                                                                   
      C2= (D_O*1000.0)/(PA*(273.16)**X)       
C      SPECFIC HEAT=0.5 FOR SELECTED PATTERN CASES.                             
      IF (ITS.GT.1) GO TO 200          
      C3= 0.0222384                                                             
      C4=0.00176  
      RETURN                                                                    
  200 C3=0.5    
      C4=0.0    
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SPHEAT (N,T,C3,C4,CI)                                          
C*******************************************************************************
C     COMPUTES THE SPECIFIC HEAT OF ICE FOR EACH LAYER                          
C*******************************************************************************
      DIMENSION T(100),CI(100)                                                  
C      SPECFIC HEAT=0.5 FOR SELECTED PATTERN CASES.                             
      IF (C4.LT.0.00001) GO TO 101                                              
      DO 100 I=1,N                                                              
      CI(I)=C3+C4*T(I)                                                          
  100 CONTINUE                                                                  
      RETURN                                                                    
  101 DO 102 I=1,N                                                              
  102 CI(I)=C3                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE GROUND (TG,C1,TCG,LS,DCG,FG,DCDT)                              
C*******************************************************************************
C     COMPUTES EFFECTIVE THERMAL CONDUCTIVITY OF THE SOIL INCLUDING             
C     VAPOR DIFFUSION.                                                          
C*******************************************************************************
      REAL LS,LV                                                                
      LV=597.3                                                                  
      AA= 6141.9/TG                                                             
      IF (TG.GT.273.16) GO TO 101                                               
      EG= C1*EXP(-6141.9/TG)                                                    
      GO TO 102                                                                 
  101 TGC= TG-273.16                                                            
      EG= 2.7489E8*EXP(-4278.63/(TGC+242.792))                                  
  102 DCDT= (EG/(4615.0*TG*TG))*(AA-1.0)                                        
C     IF THE SOIL IS BELOW ZERO DEGREES CELSIUS THEN THE LATENT HEAT OF         
C         SUBLIMATION IS USED.  OTHERWISE THE LATENT HEAT OF VAPORIZATION       
C         IS USED.                                                              
      IF (TG.GE.273.16) GO TO 103                                               
      FG= TCG+LS*DCG*DCDT                                                       
      RETURN                                                                    
  103 FG=TCG+LV*DCG*DCDT                                                        
      RETURN                                                                    
      END                                                                       
      SUBROUTINE MATRIX (A,R,C,N)                                               
C*******************************************************************************
C     SUBROUTINE FOR MATRIX SOLUTION BY GAUSS ELIMINATION.                      
C        ALSO SOLVES SIMPLE EQUATION WITH ONE UNKNOWN WHEN                      
C        A CORRECTION IS ONLY NEEDED FOR ONE LAYER.                             
C     SOLUTION IS FOR A THREE WIDE BANDED MATRIX WITH THE                       
C        FIRST ROW COMPLETELY FULL.  THE FIRST DOES NOT HAVE                    
C        TO BE COMPLETELY FULL.                                                 
C*******************************************************************************
      DIMENSION A(100,100),R(100),C(100)                                        
      IF (N.EQ.1) GO TO 15                                                      
      DO 10 I= 2,N                                                              
      J=I-1                                                                     
      IF (A(I,J).EQ.0.0) GO TO 10                                               
      FM= -A(I,J)/A(J,J)                                                        
      DO 11 K=I,N                                                               
   11 A(I,K)= A(I,K)+FM*A(J,K)                                                  
      R(I)=R(I)+FM*R(J)                                                         
   10 CONTINUE                                                                  
   15 C(N)= R(N)/A(N,N)                                                         
      IF (N.EQ.1) RETURN                                                        
      M=N-1                                                                     
      DO 12 J=1,M                                                               
      I=N-J                                                                     
      L=I+1                                                                     
      SUM=0.0                                                                   
      DO 13 K=L,N                                                               
   13 SUM=SUM+A(I,K)*C(K)                                                       
   12 C(I)= (R(I)-SUM)/A(I,I)                                                   
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RETAIN (N,TT,TTDT,RTT,RTTDT)                                   
C*******************************************************************************
C     RETAINS THE TEMPERATURE AT TIME T AND T+DT SO THAT THE                     
C        CHANGE IN TEMPERATURE DUE TO HEAT TRANSFER CAN BE USED                 
C        TO PROJECT A FIRST GUESS FOR THE NEXT TIME PERIOD.                     
C*******************************************************************************
      DIMENSION TT(100),TTDT(100),RTT(100),RTTDT(100)                           
      DO 100 I=1,N                                                              
      RTT(I)=TT(I)                                                              
      RTTDT(I)=TTDT(I)                                                          
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

      SUBROUTINE VAPOR(N,P,D_O,PA,X,TT,TTDT,D,VAPOUR,DTG,DCG,TGT,                
     1TGTDT,EA,FU,DELTAT,WSTDT,NOBS,WNTDT,NOKNOW,SOILVT)                        
C*******************************************************************************
C     COMPUTES THE CHANGE IN DENSITY DUE TO VAPOR TRANSFER. ALSO                
C        COMPUTES THE GAIN OR LOSS OF VAPOR AT THE SNOW-AIR INTERFACE           
C        (VAPOUR) AND ACROSS THE SOIL-SNOW INTERFACE(SOILVT).                   
C*******************************************************************************
C Vapor transfer (VAPOUR) is positive downward, thus evaporation or sublimation
C will be negative in sign.
      REAL LS,LV                                                                
      COMMON/SSVT/SUMSVT                                                        
      DIMENSION P(100),PT(100),TT(100),TTDT(100),D(100),Z(100),F2T(100),        
     1F3T(100),DTT(100),DT2T(100),F2TDT(100),F3TDT(100),                        
     2DTTDT(100),DT2TDT(100)                                                    
      DATA IFIRST/0/                                                            
      IF (IFIRST.GT.0) GO TO 107                                                
      SUMSVT=0.0                                                                
      IFIRST=1                                                                  
  107 C1=3.5558E10                                                              
      LS=677.0                                                                  
      LV=597.3                                                                  
C     SAVE DENSITY AT TIME T FOR USE LATER IN THE SUBROUTINE                    
      DO 100 I=1,N                                                              
  100 PT(I)=P(I)                                                                
C*******************************************************************************
C     ESTABLISH VALUES NEEDED IN THIS SUBROUTINE                                
      CALL ZDEPTH (N,D,Z,TDEPTH)                                                
      PW=1.0                                                                    
      C2=(D_O*(1000.0/PA))/(273.16**X)                                           
C Remove time controls in hours
C      DT=3600.0*DELTAT        
      DT=DELTAT                                                   
      CALL CDTDZ(N,TT,F2T,F3T,X,Z,DTT,DT2T)                                     
      CALL CDTDZ(N,TTDT,F2TDT,F3TDT,X,Z,DTTDT,DT2TDT)                           
C     COMPUTE AVERAGE SURFACE VAPOR PRESSURE OF THE SNOW COVER.                 
      EOT=C1*EXP(-6141.9/TT(1))                                                 
      EOTDT=C1*EXP(-6141.9/TTDT(1))                                             
      EO=(EOT+EOTDT)*0.5                                                        
      VG=0.0                                                                    
      SUMWE=0.0                                                                 
      SOILVT=0.0                                                                
      IF (N.EQ.1) GO TO 110                                                     
C*******************************************************************************
C     MORE THAN ONE LAYER,BOTTOM LAYER EXISTS.                                  
      DL4=D(N)*0.5                                                              
      DL5=D(N-1)*0.5                                                            
C     COMPUTE DCDT FOR SOIL.                                                    
      TCG=0.0                                                                   
      CALL GROUND(TGT,C1,TCG,LS,DCG,FG,DCGT)                                    
      CALL GROUND(TGTDT,C1,TCG,LS,DCG,FG,DCGTDT)                                
      IF ((TT(N).GT.273.159).AND.(TTDT(N).GT.273.159)) GO TO 101                
C     BOTTOM LAYER LESS THAN ZERO DEGREES CELSIUS.                              
      VG=((0.5*C2*F2T(N)*DCG*DCGT*(TGT-TT(N)))/(DCG*DCGT*                       
     1DL4+C2*F2T(N)*DTG))+((0.5*C2*F2TDT(N)*DCG*DCGTDT*                         
     2(TGTDT-TTDT(N)))/(DCG*DCGTDT*DL4+C2*F2TDT(N)*DTG))                        
      GO TO 102                                                                 
C     BOTTOM LAYER AT ZERO DEGREES CELSIUS.                                     
  101 VG=0.5*DCG*DCGT*(TGT-TT(N))/DTG+0.5*DCG*DCGTDT*                           
     1(TGTDT-TTDT(N))/DTG                                                       
  102 T1=((0.5*C2*F2T(N)*F2T(N-1)*(TT(N-1)-TT(N)))/                             
     1(F2T(N-1)*DL4+F2T(N)*DL5))+((0.5*C2*F2TDT(N)*F2TDT(N-1)*                  
     2(TTDT(N-1)-TTDT(N)))/(F2TDT(N-1)*DL4+F2TDT(N)*DL5))                       
C     CHANGE IN BOTTOM LAYER DENSITY.                                           
      P(N)=PT(N)+(DT/D(N))*T1                                                   
C     CHANGE IN BOTTOM DUE TO SOIL-SNOW VAPOR TRANSFER.                         
      SUMSVT=SUMSVT+(DT/PW)*VG                                                  
      IF ((TT(N).GT.273.159).AND.(TTDT(N).GT.273.159)) GO TO 103                
      P(N)=P(N)+(DT/D(N))*VG                                                    
      GO TO 104                                                                 
  103 SOILVT=(DT/PW)*VG                                                         
      SUMWE=SOILVT*PW                                                           
      WNTDT=WNTDT+SOILVT                                                        
      IF (WNTDT.LT.0.0) GO TO 106                                               
      NOKNOW=1                                                                  
      GO TO 104                                                                 
  106 P(N)=P(N)+(WNTDT*PW)/D(N)                                                 
      SUMWE=SUMWE+WNTDT*PW                                                      
      WNTDT=0.0                                                                 
  104 IF (N.LE.2) GO TO 110                                                     
C*******************************************************************************
C     INTERMEDIATE LAYERS EXIST.                                                
      NN=N-1                                                                    
      DO 105 I=2,NN                                                             
      P(I)=PT(I)+0.5*C2*DT*(F2T(I)*DT2T(I)+F2TDT(I)*                            
     1DT2TDT(I)+F3T(I)*DTT(I)*DTT(I)+F3TDT(I)*DTTDT(I)*DTTDT(I))                
  105 CONTINUE                                                                  
C*******************************************************************************
C     SURFACE LAYER DENSITY CHANGES.                                            
C     COMPUTE NET VAPOR CHANGE FOR ALL OTHER LAYERS.                            
  110 IF (N.EQ.1) GO TO 112                                                     
      DO 111 I=2,N                                                              
  111 SUMWE=D(I)*(P(I)-PT(I))+SUMWE                                             
C     COMPUTE SURFACE DENSITY CHANGE DUE TO VAPOR TRANSFER.                     
C     BETWEEN THE SURFACE LAYER AND THE REST OF THE SNOW COVER.                 
  112 P(1)=PT(1)-(SUMWE/D(1))+(DT*VG)/D(1)                                      
      VAPOUR=0.0                                                                
      IF (NOBS.EQ.1) RETURN                                                     
C*******************************************************************************
C     COMPUTE THE EFFECT OF VAPOR EXCHANGE BETWEEN THE                          
C        AIR AND SNOW ON THE SURFACE LAYER.                                     
      IF (EA.GT.EO) GO TO 120                                                   
      IF (EA.EQ.EO) GO TO 130                                                   
C     VAPOR LOSS - CHECK FOR EVAPORATION.                                       
      IF (WSTDT.LE.0.0) GO TO 115                                               
      PVT=FU*(EA-EO)*(LS/LV)                                                    
      IF (WSTDT.GT.(-PVT)) GO TO 113                                            
C     POTENTIAL VAPOR TRANSFER EXCEEDS LIQUID-WATER:SUBLIMATE THE               
C        REMAINDER.                                                             
      VAPOUR=-WSTDT                                                             
      WSTDT=0.0                                                                 
      PVT=(PVT-VAPOUR)*(LV/LS)                                                  
      GO TO 116                                                                 
C     LIQUID-WATER IN SURAFCE LAYER SATISFIES EVAPORATION DEMAND.               
  113 VAPOUR=PVT                                                                
      WSTDT=WSTDT+VAPOUR                                                        
C     EVAPORATION DOES NOT EFFECT DEPTH OR DENSITY.                             
      GO TO 130                                                                 
C     SUBLIMATION OCCURS.                                                       
  115 PVT=FU*(EA-EO)                                                            
C     SUBLIMATION DECREASES DEPTH                                               
  116 VAPOUR=VAPOUR+PVT                                                         
      D(1)=D(1)+(PVT*PW)/P(1)                                                   
      IF (D(1).GT.0.0) GO TO 130                                                
      D(1)=D(1)-(PVT*PW)/P(1)                                                   
      VAPOUR=VAPOUR-PVT                                                         
      GO TO 130                                                                 
C     VAPOR GAIN - CHECK IF SURFACE AT ZERO DEGREES CELSIUS.                    
  120 IF (TTDT(1).GT.273.159) GO TO 125                                         
C     FROST FORMATION -- ASSUME DENSITY THE SAME,DEPTH INCREASES.               
      VAPOUR=FU*(EA-EO)                                                         
      D(1)=D(1)+(VAPOUR*PW)/P(1)                                                
      GO TO 130                                                                 
C     CONDENSATION ON TO THE SURFACE                                            
  125 VAPOUR=FU*(EA-EO)*(LS/LV)                                                 
      WSTDT=WSTDT+VAPOUR                                                        
  130 CONTINUE                                                                  
C*******************************************************************************
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CDTDZ(N,T,F2,F3,X,Z,DTDZ,DTDZ2)                                
C*******************************************************************************
C     COMPUTES VALUES NEEDED BY SUBROUTINE VAPOR INCLUDING                      
C        THE FINITE DIFFERENCE APPROXIMATIONS TO THE                            
C        FIRST AND SECOND PARTIALS OF TEMPERATURE WITH RESPECT TO DEPTH.        
C*******************************************************************************
      DIMENSION T(100),F2(100),F3(100),Z(100),DTDZ(100),DTDZ2(100)              
      XM=X-1.0                                                                  
      C1=3.5558E10                                                              
      DO 100 I=1,N                                                              
      EI=C1*EXP(-6141.9/T(I))                                                   
      AA=6141.9/T(I)                                                            
      TP=EI/(4615.0*T(I)*T(I))                                                  
      DCDT=TP*(AA-1.0)                                                          
      TP=TP/T(I)                                                                
      DCDT2=TP*(AA**2-4.0*AA+2.0)                                               
      TX=T(I)**X                                                                
      TXM=T(I)**XM                                                              
      F2(I)=TX*DCDT                                                             
      F3(I)=TX*DCDT2+X*DCDT*TXM                                                 
C     REMAINING VALUES NOT NEEDED FOR LAYER N AND THE SURFACE LAYER.            
      IF (I.EQ.1) GO TO 100                                                     
      IF (I.EQ.N) GO TO 100                                                     
      DL1=Z(I)-Z(I-1)                                                           
      DL2=Z(I+1)-Z(I)                                                           
      F6=DL1*T(I+1)+(DL2-DL1)*T(I)-DL2*T(I-1)                                   
      F7=DL1*T(I+1)-(DL1+DL2)*T(I)+DL2*T(I-1)                                   
      DTDZ(I)=F6/(2.0*DL1*DL2)                                                  
      DTDZ2(I)=(2.0*F7)/(DL1*DL2*(DL1+DL2))                                     
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

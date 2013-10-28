      SUBROUTINE META(N,P,D,WTDT,TTDT,PLWHC,C1,C2,C3,C4,C5,DELTAT,DMETA)        
C*******************************************************************************
C     COMPUTES THE CHANGE IN DENSITY OF THE SNOW COVER CAUSED BY                
C        DESTRUCTIVE(EQUI-TEMPERATURE)METAMORPHISM,COMPACTION,AND THE           
C        PRESENCE OF LIQUID-WATER.                                              
C*******************************************************************************
      DIMENSION P(100),D(100),WTDT(100),TTDT(100)                               
C     WEIGHT IS THE WATER-EQUIVALENT(CM.)ABOVE THE LAYER.                       
      WEIGHT=0.0                                                                
      PICE=0.917                                                                
      IF((C1.GT.0.0).OR.(C3.GT.0.0))GO TO 109                                   
      IF(C5.GT.0.0)GO TO 109                                                    
      RETURN                                                                    
  109 DO 100 I=1,N                                                              
      WEL=P(I)*D(I)                                                             
C*******************************************************************************
C     DESTRUCTIVE METAMORPHISM TERM.                                            
      IF(P(I).GE.PICE) GO TO 101                                                
      TERM=EXP(-C4*(273.16-TTDT(I)))                                            
      IF (P(I).LE.DMETA) GO TO 107                                              
      T1=TERM*C3*EXP(-46.0*(P(I)-DMETA))                                        
      GO TO 102                                                                 
  107 T1=TERM*C3                                                                
      GO TO 102                                                                 
  101 T1=0.0                                                                    
C*******************************************************************************
C     COMPACTION TERM.                                                          
  102 IF(P(I).GE.PICE)GO TO 103                                                 
      TERM=EXP(-0.08*(273.16-TTDT(I)))                                          
      T2=WEIGHT*C1*EXP(-C2*P(I))*TERM                                           
      GO TO 104                                                                 
  103 T2=0.0                                                                    
C*******************************************************************************
C     LIQUID-WATER TERM.                                                        
  104 IF (WTDT(I).LE.0.0) GO TO 106                                             
      T1=C5*T1                                                                  
C*******************************************************************************
C     DENSIFICATION OF THE LAYER.  (WATER-EQUIVALENT STAYS THE SAME.)           
  106 P(I)=(1.0+DELTAT*(T1+T2))*P(I)                                            
      D(I)=WEL/P(I)                                                             
      WEIGHT=WEIGHT+WEL+WTDT(I)                                                 
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

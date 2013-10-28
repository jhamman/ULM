      SUBROUTINE NWSNOW(N,TT,WT,TTDT,WTDT,NOKNOW,D,P,RTT,RTTDT,PX,PNS,          
     &   TSNOW,SUMPX,THICK,ITOPT)                                               
C*******************************************************************************
C     CHECKS FOR NEW SNOWFALL AND ADDS NEW SNOW TO THE SNOW COVER.              
C*******************************************************************************
      DIMENSION TTDT(100),WTDT(100),NOKNOW(100),D(100),P(100),RTT(100),         
     &   RTTDT(100),TT(100),WT(100)                                             
      LF=79.7                                                                   
      PW=1.0                                                                    
      NL=0                                                                      
C*******************************************************************************
C     CHECK IF PRECIPITATION OCCURRED.                                          
      IF(PX.EQ.0.0)RETURN                                                       
C     CHECK FOR RAIN OR SNOW.                                                   
      IF(PNS.LT.1.0)GO TO 100                                                   
C*******************************************************************************
C     PRECIPITATION IS RAIN--CHECK IF SNOW IS PRESENT.                          
      IF(N.EQ.0)RETURN                                                          
      SUMPX=SUMPX+PX                                                            
      RETURN                                                                    
C*******************************************************************************
C     PRECIPITATION IS SNOW                                                     
  100 IF(ITOPT.NE.1)GO TO 101                                                   
      TST=TT(1)                                                                 
      TSTDT=TTDT(1)                                                             
      TTDT(1)=TST     
C Depth of new snow                                                          
  101 DNS=PX/PNS                                                                
      SUMPX=SUMPX+PX                                                            
      IF (N.GT.0) GO TO 107      
C No existing layers                                               
      DC=0.90*THICK                                                             
      TDNS=DNS                                                                  
      GO TO 112                                                                 
  107 IF (NOKNOW(1).GE.0) GO TO 103                                             
      DO 106 I=1,N                                                              
      TTDT(I)=TT(I)                                                             
      WTDT(I)=WT(I)                                                             
  106 CONTINUE                                                                  
C     CHECK IF NEW SNOW CAUSES THE SURFACE LAYER TO EXCEED THE                  
C        SPECIFIED DEPTH LIMIT.                                                 
  103 IF((D(1)+DNS).GT.(1.55*THICK))GO TO 110                                   
C*******************************************************************************
C     INCORPORATE NEW SNOW INTO OLD SURFACE LAYER.                              
      I=1                                                                       
  102 WEI=P(I)*D(I)                                                             
      WENS=PX                                                                   
      WEL=WEI+WENS                                                              
      D(I)=D(I)+DNS                                                             
      P(I)=WEL/D(I)                                                             
      TTDT(I)=(WENS*TSNOW+WEI*TTDT(I))/WEL                                      
      IF(NOKNOW(I).EQ.0)GO TO 105                                               
C     LIQUID-WATER AND SUBFREEZING TEMPERATURES MAY EXIST IN NEW LAYER.         
      FREEZE=((273.16-TTDT(I))*0.503*WEL)/(LF*PW)                               
      IF(WTDT(I).GT.FREEZE)GO TO 104                                            
      TTDT(I)=TTDT(I)+(WTDT(I)*LF*PW)/(0.503*WEL)                               
      WEL=WEL+WTDT(I)                                                           
      P(I)=WEL/D(I)                                                             
      WTDT(I)=0.0                                                               
      NOKNOW(I)=0                                                               
      GO TO 105                                                                 
  104 WTDT(I)=WTDT(I)-FREEZE                                                    
      WEL=WEL+FREEZE                                                            
      P(I)=WEL/D(I)                                                             
      TTDT(I)=273.16                                                            
C     SET UNKNOWN NEGATIVE TO INDICATE ALL OR PART OF THE LAYER IS NEW          
C        SNOW.                                                                  
  105 NOKNOW(I)=NOKNOW(I)-2                                                     
      TT(I)=TTDT(I)                                                             
      WT(I)=WTDT(I)                                                             
  114 PX=0.0                                                                    
      N=N+NL                                                                    
      IF (ITOPT.NE.1) RETURN                                                    
      TT(1)=TST                                                                 
      TTDT(1)=TSTDT                                                             
      RETURN                                                                    
C*******************************************************************************
C     SURFACE LAYER EXCEEDS LIMIT--SUB-DIVIDE LAYER.                            
  110 IF(D(1).LE.(0.90*THICK))GO TO 115                                         
C     MOVE 90 PERCENT OF OLD LAYER DOWN.                                        
      DO 111 J=1,N                                                              
      L=N-J+1                                                                   
      LL=L+1                                                                    
      P(LL)=P(L)                                                                
      TTDT(LL)=TTDT(L)                                                          
      WTDT(LL)=WTDT(L)                                                          
      D(LL)=D(L)                                                                
      RTT(LL)=RTT(L)                                                            
      RTTDT(LL)=RTTDT(L)                                                        
      NOKNOW(LL)=NOKNOW(L)                                                      
      IF (NOKNOW(LL).GE.0) GO TO 111                                            
      TT(LL)=TTDT(LL)                                                           
      WT(LL)=WTDT(LL)                                                           
  111 CONTINUE                                                                  
      N=N+1                                                                     
      D(2)=0.90*THICK                                                           
      WTDT(2)=(D(2)/D(1))*WTDT(1)                                               
      D(1)=D(1)-D(2)                                                            
      WTDT(1)=WTDT(1)-WTDT(2)                                                   
C     NOW SURFACE LAYER IS.LE.0.90*THICK.  NOW ADD NEW SNOW.                    
      IF((D(1)+DNS).GT.(1.55*THICK))GO TO 115                                   
      I=1                                                                       
      GO TO 102                                                                 
C*******************************************************************************
C     NEW LAYERS WHICH ARE ENTIRELY NEW SNOW ARE CREATED IN ADDITION TO         
C        A COMBINATION LAYER.                                                   
  115 TDNS=DNS                                                                  
      TPX=PX                                                                    
      DC=0.90*THICK                                                             
      DNS=DC-D(1)                                                               
      PX=(DNS/TDNS)*TPX                                                         
      TDNS=TDNS-DNS                                                             
      TPX=TPX-PX                                                                
C*******************************************************************************
C     COMPUTE THE NUMBER OF LAYERS OF COMPLETELY NEW SNOW.                      
  112 NL=TDNS/DC                                                                
      EXTRA=TDNS-NL*DC                                                          
      IF((DC+EXTRA).GT.(1.55*THICK))NL=NL+1                                     
      IF(NL.EQ.0)NL=1                                                           
      IF (N.EQ.0) GO TO 113                                                     
C     MOVE LAYERS DOWN TO MAKE ROOM FOR THE NEW SNOW LAYERS                     
C        NL IS THE NUMBER OF LAYERS COMPRISED OF ONLY NEW SNOW.                 
      DO 116 J=1,N                                                              
      L=N-J+1                                                                   
      LL=L+NL                                                                   
      D(LL)=D(L)                                                                
      P(LL)=P(L)                                                                
      TTDT(LL)=TTDT(L)                                                          
      WTDT(LL)=WTDT(L)                                                          
      RTT(LL)=RTT(L)                                                            
      RTTDT(LL)=RTTDT(L)                                                        
      NOKNOW(LL)=NOKNOW(L)                                                      
  116 CONTINUE                                                                  
  113 IF(NL.EQ.1)GO TO 120                                                      
C     FULL SIZE NEW SNOW LAYERS.                                                
      DO 117 J=2,NL                                                             
      D(J)=DC                                                                   
      P(J)=PNS                                                                  
      TT(J)=TSNOW                                                               
      WT(J)=0.0                                                                 
      NOKNOW(J)=-2                                                              
  117 CONTINUE                                                                  
  120 EXTRA=TDNS-(NL-1)*DC                                                      
      D(1)=EXTRA                                                                
      P(1)=PNS                                                                  
      TT(1)=TSNOW                                                               
      WT(1)=0.0                                                                 
      NOKNOW(1)=-2                                                              
      IF (N.EQ.0) GO TO 114                                                     
      I=NL+1                                                                    
      GO TO 102                                                                 
      END                                                                       

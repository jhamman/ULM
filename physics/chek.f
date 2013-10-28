      SUBROUTINE CHEK(N,D,P,TT,TTDT,WT,WTDT,NOKNOW,SCOUT,THICK,CTHICK,         
     1WE,TDEPTH,QI,QR,QA,DELTAT,PA,FU,TA,EA,PX,TGT,TGTDT,TCG,DCG,DTG,           
     2D_O,X,ITOPT,ITS,COEFKE)                                                          
C*******************************************************************************
C     CHECKS EACH LAYER TO MAKE SURE THAT TIME T+DT VALUES ARE                  
C          CONSISTENT WITH THE TIME T STATE OF THE LAYER. ALSO                  
C          CHECKS THE SIZE OF EACH LAYER TO MAKE SURE THE                       
C          THICKNESS IS WITHIN THE SPECIFIED LIMITS. CHANGES IN                 
C          DEPTH DURING THE REST OF THE TIME PERIOD IN                          
C          SUBROUTINES VAPOR AND META WILL BE SMALL,THUS                        
C          LIMITS CAN BE CHECKED AT THIS POINT.                                 
C*******************************************************************************
      REAL LF,LS                                                                
      COMMON/EBAL/QS,QLW,QH,QE,QG,QPX,DQ                                        
      COMMON/DHEAT/DHS                                                          
C      COMMON/EFFTC/COEFKE                                                       
      COMMON/STEMP/TOSIM                                                        
      DIMENSION D(100),P(100),TT(100),TTDT(100),WT(100),WTDT(100),              
     1NOKNOW(100),CIT(100),CITDT(100)                                           
      LF=79.7                                                                   
      PW=1.0                                                                    
C     COMPUTE ENERGY BALANCE COMPONENTS FOR THE PERIOD.                         
      IF(ITOPT.EQ.1) GO TO 111                                                  
      QS=QI-QR                                                                  
C      DT=3600*DELTAT   
C Replace time variable in seconds                                                         
      DT=DELTAT    
      D1=0.5                                                                    
      TG=(TGT+TGTDT)*0.5                                                        
      TN=(TT(N)+TTDT(N))*0.5                                                    
      CALL VSTART(QIR,QI,QR,QA,SIGMA,DT,CW,LS,PW,PA,GAMMA,D1,D2,C1,LF,          
     1C2,D_O,X,C3,C4,E,RW,ITS)                                                  
      TS=(TT(1)+TTDT(1))*0.5                                                    
      IF (N.GT.1) TS=TOSIM                                                      
      TS4=TS*TS*TS*TS                                                           
      QLW=E*(QA-SIGMA*TS4)                                                      
      QH=GAMMA*PW*LS*FU*(TA-TS)                                                 
      ES=C1*EXP(-6141.9/TS)                                                     
      QE=LS*PW*FU*(EA-ES)                                                       
      CALL GROUND(TG,C1,TCG,LS,DCG,FG,DCDT)                                     
      IF(TN.LT.273.16)GO TO 161                                                 
      QG=DT*FG*(TG-TN)/DTG                                                      
      GO TO 162                                                                 
  161 DN=C2*(TN**X)                                                             
      TCN=0.00005+COEFKE*P(N)*P(N)                                              
      CALL GROUND(TN,C1,TCN,LS,DN,FN,DCDT)                                      
      DN=D(N)*0.5                                                               
      QG=(DT*FN*FG*(TG-TN))/(FG*DN+FN*DTG)                                      
  162 CALL WTBULB(TW,TA,EA,PA)                                                  
      IF (TW.LT.273.16) TW=273.16                                               
      QPX=CW*PW*PX*(TW-273.16)                                                  
      CALL SPHEAT(N,TT,C3,C4,CIT)                                               
      CALL SPHEAT(N,TTDT,C3,C4,CITDT)                                           
      DQ=0.0                                                                    
      DO 163 I=1,N                                                              
      DQ=DQ+P(I)*D(I)*(CITDT(I)*TTDT(I)-CIT(I)*TT(I))                           
     1+LF*PW*(WTDT(I)-WT(I))                                                    
  163 CONTINUE                                                                  
      DHS=DQ-P(1)*D(1)*(CITDT(1)*TTDT(1)-CIT(1)*TT(1))                          
     1-LF*PW*(WTDT(1)-WT(1))                                                    
C*******************************************************************************
C     COMPUTE THE CHANGE IN LIQUID-WATER IN EACH LAYER. THIS AMOUNT             
C          NEEDS TO BE ADDED TO ICE CONTENT(NEGATIVE CHANGE) OR                 
C          SUBTRACTED FROM THE ICE CONTENT(POSITIVE CHANGE-MELT).               
C          VAPOR TRANSFER AND RAIN ARE NOT INCLUDED YET AND WILL                
C          BE ACCOUNTED FOR IN SUBSEQUENT SUBROUTINES.                          
  111 NL=N                                                                      
      ABOVE=0.0                                                                 
      EXCESS=0.0                                                                
      DO 100 I=1,N                                                              
      IF (NOKNOW(I).EQ.1) GO TO 103                                             
C*******************************************************************************
C     TEMPERATURE IS UNKNOWN--SOME EXCESS LIQUID-WATER FROM                     
C          THE ABOVE LAYER MUST BE FROZEN.                                      
      IF (EXCESS.EQ.0.0) GO TO 103                                              
C     COMPUTE AMOUNT TO BE FROZE IN ORDER TO RAISE THE TEMPERATURE.             
C          TO ZERO DEGREES CELSIUS. (USE SPECIFIC HEAT OF ICE=.503)             
      FREEZE=((273.16-TTDT(I))*0.503*P(I)*D(I))/(LF*PW)                         
      WEL=P(I)*D(I)                                                             
      IF (EXCESS.GT.FREEZE) GO TO 109                                           
C     EXCESS IS ALL FROZEN IN THIS LAYER.                                       
      TTDT(I)=TTDT(I)+(EXCESS*LF*PW)/(0.503*P(I)*D(I))                          
      WEL=WEL+EXCESS                                                            
      P(I)=WEL/D(I)                                                             
      EXCESS=0.0                                                                
      ABOVE=0.0                                                                 
      GO TO 103                                                                 
C     EXCESS EXCEEDS REFREEZE.                                                  
  109 TTDT(I)=273.16                                                            
      WEL=WEL+FREEZE                                                            
      P(I)=WEL/D(I)                                                             
      EXCESS=EXCESS-FREEZE                                                      
      ABOVE=ABOVE-FREEZE                                                        
      NOKNOW(I)=1                                                               
C*******************************************************************************
  103 WTDT(I)=WTDT(I)+EXCESS                                                    
      CHANGE=WTDT(I)-WT(I)-ABOVE                                                
      IF (ABS(CHANGE).LT.0.00001) GO TO 104                                     
      IF (CHANGE.GT.0.0) GO TO 101                                              
C     SOME LIQUID-WATER IS FROZEN. ADD TO ICE CONTENT OF LAYER.                 
      WEL=P(I)*D(I)                                                             
      WEL=WEL-CHANGE                                                            
      P(I)=WEL/D(I)                                                             
      GO TO 104                                                                 
C     MELT HAS OCCURRED,SUBTRACT FROM ICE CONTENT                               
  101 WEL=P(I)*D(I)                                                             
C     IF MELT EXCEEDS ICE CONTENT OF LAYER--LAYER IS GONE.                      
      IF (CHANGE.GE.WEL) GO TO 102                                              
C     LAYER REMAINS                                                             
      WEL=WEL-CHANGE                                                            
      D(I)=WEL/P(I)                                                             
      GO TO 104                                                                 
C     LAYER IS GONE                                                             
C          LIQUID-WATER IS ADDED TO THE LAYER BELOW.(EXCESS).THE ICE            
C          CONTENT PLUS WT OF THE LAYER CANNOT BE TAKEN FROM THE NEXT           
C          LAYER,BUT IS STILL ADDED TO THE LIQUID-WATER CONTENT OF              
C          THE NEXT LAYER (ABOVE).                                              
  102 NL=NL-1                                                                   
      ABOVE=ABOVE+WEL+WT(I)                                                     
      EXCESS=WTDT(I)                                                            
      D(I)=0.0                                                                  
      GO TO 100                                                                 
  104 ABOVE=0.0                                                                 
      EXCESS=0.0                                                                
  100 CONTINUE                                                                  
C*******************************************************************************
C     CHECK TO SEE IF ENTIRE SNOW COVER IS GONE.                                
      IF (NL.GT.0) GO TO 105                                                    
C     NEGATIVE N INDICATES TO SUBROUTINE STATDA THAT THE SNOW COVER HAS         
C         JUST DISAPPEARED.                                                     
      N=-1                                                                      
      WE=0.0                                                                    
      TDEPTH=0.0                                                                
      SCOUT=ABOVE                                                               
      RETURN                                                                    
C     ELIMINATE LAYERS WHICH ARE GONE.                                          
  105 IF (NL.EQ.N) GO TO 110                                                    
      DO 106 I=1,NL                                                             
      IF (D(I).GT.0.0) GO TO 106                                                
C     LAYER GONE. MOVE OTHER LAYERS UP.                                         
      NEXT=I+1                                                                  
      DO 107 J=NEXT,N                                                           
      L=J-1                                                                     
      D(L)=D(J)                                                                 
      P(L)=P(J)                                                                 
      TT(L)=TT(J)                                                               
      TTDT(L)=TTDT(J)                                                           
      WT(L)=WT(J)                                                               
      WTDT(L)=WTDT(J)                                                           
      NOKNOW(L)=NOKNOW(J)                                                       
  107 CONTINUE                                                                  
      N=N-1                                                                     
      IF (N.EQ.NL) GO TO 108                                                    
  106 CONTINUE                                                                  
  108 WTDT(N)=WTDT(N)+ABOVE                                                     
  110 CONTINUE                                                                  
      IF (N.EQ.1) GO TO 160                                                     
C*******************************************************************************
C     CHECK THICKNESS OF EACH LAYER. IF NOT WITHIN SPECIFIED                    
C          LIMITS,DIVIDE (IF TOO LARGE) OR ADD TO AN ADJACENT                   
C          LAYER (IF TOO SMALL).                                                
      TDEPTH=0.0                                                                
      I=1                                                                       
  150 ZZ=TDEPTH+D(I)*0.5                                                        
C*******************************************************************************
C     COMPUTED DESIRED THICKNESS FOR THE LAYER.                                 
      IF (ZZ.GT.30.0) GO TO 121                                                 
      DZ=THICK                                                                  
      GO TO 122                                                                 
  121 DZ=CTHICK*(ZZ-30.0)+THICK                                                 
C     CHECK ACTUAL THICKNESS AGAINST DESIRED THICKNESS.                         
  122 IF (D(I).GT.1.55*DZ) GO TO 125                                            
      IF (D(I).LT.0.55*DZ) GO TO 130                                            
C     THICKNESS IS WITHIN SPECIFIED LIMITS.                                     
      TDEPTH=TDEPTH+D(I)                                                        
      GO TO 120                                                                 
C*******************************************************************************
C     THICKNESS IS GREATER THAN SPECIFIED LIMTS.                                
C          SUB-DIVIDE LAYER,THUS CREATING A NEW LAYER.                          
C          PROPERTIES ARE THE SAME FOR EACH.                                    
C          MOVE OTHER LAYERS DOWN.                                              
  125 NEXT=I+1                                                                  
      IF (NEXT.GT.N) GO TO 127                                                  
      DO 126 J=NEXT,N                                                           
      L=N-J+NEXT                                                                
      LL=L+1                                                                    
      D(LL)=D(L)                                                                
      P(LL)=P(L)                                                                
      TT(LL)=TT(L)                                                              
      TTDT(LL)=TTDT(L)                                                          
      WT(LL)=WT(L)                                                              
      WTDT(LL)=WTDT(L)                                                          
      NOKNOW(LL)=NOKNOW(L)                                                      
  126 CONTINUE                                                                  
  127 TDEPTH=TDEPTH+D(I)                                                        
      D(NEXT)=D(I)*0.5                                                          
      P(NEXT)=P(I)                                                              
      TT(NEXT)=TT(I)                                                            
      TTDT(NEXT)=TTDT(I)                                                        
      WT(NEXT)=WT(I)*0.5                                                        
      WTDT(NEXT)=WTDT(I)*0.5                                                    
      NOKNOW(NEXT)=NOKNOW(I)                                                    
      D(I)=D(I)*0.5                                                             
      WT(I)=WT(I)*0.5                                                           
      WTDT(I)=WTDT(I)*0.5                                                       
      I=I+1                                                                     
      N=N+1                                                                     
      GO TO 120                                                                 
C*******************************************************************************
C     THICKNESS IS SMALLER THAN SPECIFIED LIMITS.                               
C          ADD TO THE SMALLEST ADJACENT LAYER, THUS                             
C          LOSING A LAYER. PROPERIES OF THE NEW LAYER                           
C          ARE THE WEIGHTED AVERAGE OF THE TWO FORMER                           
C          LAYERS. MOVE OTHER LAYERS UP.                                        
  130 IF (I.EQ.1) GO TO 131                                                     
      IF (I.EQ.N) GO TO 132                                                     
      NL=I-1                                                                    
      IF (D(I+1).LT.D(I-1)) NL=I+1                                              
      GO TO 135                                                                 
  131 NL=2                                                                      
      GO TO 135                                                                 
  132 NL=N-1                                                                    
C     NL IS THE SMALLEST ADJACENT LAYER.                                        
C          ADD LAYER I TO LAYER NL                                              
  135 WEI=P(I)*D(I)                                                             
      WENL=P(NL)*D(NL)                                                          
      WEL=WEI+WENL                                                              
      D(NL)=D(NL)+D(I)                                                          
      P(NL)=WEL/D(NL)                                                           
      WT(NL)=WT(NL)+WT(I)                                                       
      WTDT(NL)=WTDT(NL)+WTDT(I)                                                 
      TT(NL)=(WENL*TT(NL)+WEI*TT(I))/WEL                                        
      TTDT(NL)=(WENL*TTDT(NL)+WEI*TTDT(I))/WEL                                  
      IF (NOKNOW(I).EQ.NOKNOW(NL)) GO TO 140                                    
C     UNKNOWNS ARE DIFFERENT. COMPUTE THE UNKNOWN FOR THE NEW LAYER.            
      FREEZE=((273.16-TTDT(NL))*0.503*P(NL)*D(NL))/(LF*PW)                      
      IF (WTDT(NL).GT.FREEZE) GO TO 136                                         
C     TEMPERATURE IS UNKNOWN                                                    
      TTDT(NL)=TTDT(NL)+(WTDT(NL)*LF*PW)/(0.503*P(NL)*D(NL))                    
      WEL=WEL+WTDT(NL)                                                          
      P(NL)=WEL/D(NL)                                                           
      NOKNOW(NL)=0                                                              
      GO TO 140                                                                 
C     LIQUID-WATER IS UNKNOWN.                                                  
  136 WTDT(NL)=WTDT(NL)-FREEZE                                                  
      WEL=WEL+FREEZE                                                            
      P(NL)=WEL/D(NL)                                                           
      NOKNOW(NL)=1                                                              
  140 IF (NOKNOW(NL).EQ.1) TTDT(NL)=273.16                                      
      IF (NOKNOW(NL).EQ.0) WTDT(NL)=0.0                                         
C     MOVE OTHER LAYERS UP.                                                     
      IF (NL.LT.I) TDEPTH=TDEPTH+D(I)                                           
      NEXT=I+1                                                                  
      IF (NEXT.GT.N) GO TO 145                                                  
      DO 141 J=NEXT,N                                                           
      L=J-1                                                                     
      D(L)=D(J)                                                                 
      P(L)=P(J)                                                                 
      TT(L)=TT(J)                                                               
      TTDT(L)=TTDT(J)                                                           
      WT(L)=WT(J)                                                               
      WTDT(L)=WTDT(J)                                                           
      NOKNOW(L)=NOKNOW(J)                                                       
  141 CONTINUE                                                                  
  145 I=I-1                                                                     
      N=N-1                                                                     
C*******************************************************************************
C     CHECK FOR LAST LAYER.                                                     
  120 IF (I.EQ.N) GO TO 160                                                     
      IF (N.EQ.1) GO TO 160                                                     
      I=I+1                                                                     
      GO TO 150                                                                 
  160 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

      SUBROUTINE WATER(DELTAT,N,P,D,WTDT,TTDT,NOKNOW,PLWHC,PLWMAX,              
     &     PLWDEN,PX,SCOUT,MONTH,IDAY,IYEAR,IHOUR,IPUNCH,CW1,CW2,CW3,CW4)           
C*******************************************************************************
C     THIS SUBROUTINE DETERMINES THE SNOW COVER OUTFLOW DURING EACH PERIOD      
C     BASED ON THE CHANGE IN LIQUID-WATER DURING THE PERIOD AND                 
C     THE PHYSICAL CHARTERISTICS OF THE SNOW COVER.                             
C*******************************************************************************
      REAL LF                                                                   
      COMMON/WLAGS/WLAG(11),NLAG,STORE                                          
C      COMMON/WPM/CW1,CW2,CW3,CW4                                                
      DIMENSION P(100),D(100),WTDT(100),TTDT(100),NOKNOW(100)
      DIMENSION WNET(24)                                                        
C      DATA IFIRST/0/  
      IFIRST = 0
      WRITE(*,*)'WATER:TOP OF SBRTN N',N
      LF=79.7                                                                   
      PW=1.0                                                                    
      IF(N.LT.0) GO TO 135                                                      
      IF(IFIRST.GT.0) GO TO 100                                                 
C*******************************************************************************
C     INITIALIZE VARIABLES.                                                     
C     1. STORE IS THE AMOUNT OF LIQUID-WATER(THAT HAS ALREADY BEEN              
C      LAGGED) THAT IS IN STORAGE IN THE SNOW-COVER.                            
C     2. WLAG() IS THE AMOUNT OF LIQUID-WATER IN THE PROCESS OF                 
C      BEING LAGGED.  NLAG IS THE NUMBER OF ARRAY ELEMENTS USED.                
      STORE=0.0                                                                 
      IDT=DELTAT+0.01                                                           
      JDT=CW1+0.01                                                              
      NLAG=JDT+2                                                                
      IF (NLAG.GT.11) NLAG=11                                                   
      DO 101 J=1,NLAG                                                           
  101 WLAG(J)=0.0                                                               
      IFIRST=1                                                                  
C*******************************************************************************
C     DETERMINE THE EXCESS LIQUID-WATER(IN EXCESS OF LIQUID-WATER               
C      HOLDING CAPACITY) GENERATED DURING THIS TIME PERIOD.                     
C     AT THE SURFACE THE EXCESS IS THE AMOUNT OF RAIN.                          
  100 EXCESS=PX                                                                 
      IF (N.EQ.1) GO TO 102                                                     
C     EXCESS WATER IN BOTTOM LAYER IS NOT TO BE ROUTED.                         
      BOTTOM=0.0                                                                
      IF (NOKNOW(N).NE.1) GO TO 102                                             
      WEL=P(N)*D(N)                                                             
      IF (P(N).GE.PLWDEN) GO TO 103                                             
      PLW=(PLWMAX-PLWHC)*(PLWDEN-P(N))/PLWDEN+PLWHC                             
      GO TO 104                                                                 
  103 PLW=PLWHC                                                                 
  104 WMAX=PLW*WEL                                                              
      IF (WTDT(N).LE.WMAX) GO TO 102                                            
      BOTTOM=WTDT(N)-WMAX                                                       
  102 WE=0.0                                                                    
      TDEPTH=0.0                                                                
      DO 110 I=1,N                                                              
      WEL=P(I)*D(I)                                                             
      IF(NOKNOW(I).EQ.1) GO TO 115                                              
C     LAYER IS BELOW ZERO DEGREES CELSIUS.                                      
      IF(EXCESS.EQ.0.0) GO TO 112                                               
C     FREEZE SOME OF THE LIQUID-WATER.                                          
      FREEZE=((273.16-TTDT(I))*0.503*WEL)/(LF*PW)                               
      IF(EXCESS.GT.FREEZE) GO TO 111                                            
C     EXCESS IS ALL FROZEN IN THIS LAYER.                                       
      TTDT(I)=TTDT(I)+(EXCESS*LF*PW)/(0.503*WEL)                                
      WEL=WEL+EXCESS                                                            
      EXCESS=0.0                                                                
      P(I)=WEL/D(I)                                                             
      GO TO 112                                                                 
C     EXCESS EXCEEDS REFREEZE.                                                  
  111 TTDT(I)=273.16                                                            
      WEL=WEL+FREEZE                                                            
      P(I)=WEL/D(I)                                                             
      EXCESS=EXCESS-FREEZE                                                      
      NOKNOW(I)=1                                                               
  115 IF (P(I).GE.PLWDEN) GO TO 117                                             
      PLW=(PLWMAX-PLWHC)*(PLWDEN-P(I))/PLWDEN+PLWHC                             
      GO TO 118                                                                 
  117 PLW=PLWHC                                                                 
  118 WMAX=PLW*WEL                                                              
      W=WTDT(I)+EXCESS                                                          
      IF(W.GT.WMAX) GO TO 116                                                   
C     LIQUID-WATER HOLDING CAPACITY IS NOT SATISFIED.                           
      WTDT(I)=W                                                                 
      EXCESS=0.0                                                                
      GO TO 112                                                                 
C     LIQUID-WATER HOLDING CAPACITY IS EXCEEDED.                                
  116 WTDT(I)=WMAX                                                              
      EXCESS=W-WMAX                                                             
  112 WE=WE+WEL                                                                 
      TDEPTH=TDEPTH+D(I)                                                        
  110 CONTINUE                                                                  
C*******************************************************************************
      IF (IPUNCH.NE.1) GO TO 200                                                
C     STORE EXCESS WATER AND PUNCH ONCE A DAY.                                  
      IAP=IHOUR/IDT                                                             
      IF (IAP.GT.1) GO TO 201                                                   
      DO 202 I=1,24                                                             
 202  WNET(I)=0.0                                                               
 201  WNET(IAP)=EXCESS*10.0                                                     
      IF (IHOUR.NE.24) GO TO 200                                                
      NP=24/IDT                                                                 
      IF (IDT.EQ.1) NP=12                                                       
C      PUNCH 900,MONTH,IDAY,IYEAR,IDT,IDT,(WNET(I),I=1,NP)                       
C 900  FORMAT (I2,1H/,I2,1H/,I4,1H-,I2,1X,3HDT=,I2,1X,12F5.1)                    
      IF (IDT.NE.1) GO TO 200                                                   
      IHR=13                                                                    
      NP=24                                                                     
C      PUNCH 900,MONTH,IDAY,IYEAR,IHR,IDT,(WNET(I),I=IHR,NP)                     
 200  CONTINUE                                                                  
C*******************************************************************************
      IF (N.EQ.1) GO TO 133                                                     
      EXCESS=EXCESS-BOTTOM                                                      
C*******************************************************************************
C     ROUTE EXCESS WATER THROUGH THE SNOW COVER.                                
C        EMPIRICAL LAG AND ATTENUATION EQUATIONS.                               
C        ONE HOUR TIME STEP USED.                                               
      EXCESS=EXCESS/DELTAT                                                      
      DENSE=WE/TDEPTH                                                           
      SCOUT=0.0                                                                 
C Removed time constraint
C      DO 121 IHR=1,IDT                                                          
      OUTHR=0.0                                                                 
C     LAG-FUNCTION OF DEPTH,DENSITY,AND EXCESS WATER.                           
      IF(EXCESS.LT.0.001) GO TO 125                                             
      NI=((EXCESS*100.0)**0.3)+0.5                                              
      IF(NI.LT.1) NI=1                                                          
      FN=NI                                                                     
      FLMAX=CW1*(1.0-EXP(-0.0025*TDEPTH/DENSE))                                 
      DO 120 J=1,NI                                                             
      FJ=J                                                                      
      FLAG=FLMAX/(CW2*EXCESS*(FJ-0.5)/FN+1.0)                                   
      K=FLAG+1.0                                                                
      POR=K-FLAG                                                                
      WINC=1.0/FN                                                               
      WLAG(K)=WLAG(K)+EXCESS*WINC*POR                                           
      WLAG(K+1)=WLAG(K+1)+EXCESS*WINC*(1.0-POR)                                 
  120 CONTINUE                                                                  
      GO TO 130                                                                 
  125 WLAG(1)=WLAG(1)+EXCESS                                                    
C     ATTENUATION-FUNCTION OF DENSITY AND                                       
C        PREVIOUS OUTFLOW.                                                      
  130 IF((STORE+WLAG(1)).EQ.0.0) GO TO 131                                      
      R=1.0/(CW3*EXP(-CW4*WLAG(1)*DENSE/TDEPTH)+1.0)                            
      OUTHR=(STORE+WLAG(1))*R                                                   
      STORE=STORE+WLAG(1)-OUTHR                                                 
      SCOUT=SCOUT+OUTHR                                                         
      IF(STORE.GT.0.001) GO TO 131                                              
      OUTHR=OUTHR+STORE                                                         
      SCOUT=SCOUT+STORE                                                         
      STORE=0.0                                                                 
  131 NI=NLAG-1                                                                 
      DO 132 J=1,NI                                                             
  132 WLAG(J)=WLAG(J+1)                                                         
      WLAG(NLAG)=0.0                                                            
  121 CONTINUE                                                                  
      SCOUT=SCOUT+BOTTOM                                                        
      GO TO 140                                                                 
C     ONLY ONE LAYER EXISTS--OUTFLOW IS NOT DELAYED.                            
  133 SCOUT=EXCESS                                                              
C*******************************************************************************
C     SNOW COVER HAS JUST DISAPPEARED                                           
  135 DO 136 J=1,NLAG                                                           
      SCOUT=SCOUT+WLAG(J)                                                       
  136 WLAG(J)=0.0                                                               
      SCOUT=SCOUT+STORE                                                         
      STORE=0.0                                                                 
      IF (N.EQ.1) GO TO 140                                                     
      SCOUT=SCOUT+PX                                                            
C*******************************************************************************
C     CONVERT SNOW COVER OUTFLOW TO MILLIMETERS FOR OUTPUT.                     
  140 SCOUT=SCOUT*10.0                                                          
      RETURN                                                                    
      END                                                                       

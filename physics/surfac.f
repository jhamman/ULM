      SUBROUTINE SURFAC(N,TST,TSTDT,WST,WSTDT,D,P,ITOPT,QI,QR,                  
     1QA,FU,TA,EA,PX,D_O,PA,X,DTG,TGT,TGTDT,TCG,DCG,THEDA,NOKNOW,                
     2TOLER,DELTAT,ITER,MONTH,IDAY,IYEAR,IHOUR,IGUESS,QG,
     3EXTINC,CV,G1,G2,G3,IFU,UMSEC,FUCOEF,HEIGHT,RICRIT,ZO,COEFKE,ITS)             
C*******************************************************************************
C     COMPUTES TEMPERATURE AND LIQUID-WATER FOR THE SURFACE LAYER               
C        AT TIME T+DT FOR TWO CASES.                                            
C             1. TO PROVIDE AN INTIAL GUESS FOR THE SURFACE LAYER --            
C                   ASSUMES QG AND DQ ARE ZERO.                                 
C             2. COMPUTES BALANCE WHEN THE SNOW COVER CONSISTS                  
C                   OF ONLY ONE LAYER. WHEN ITOPT=1 AND THE                     
C                   TEMPERATURE AT T+DT IS 273.16, GROUND-MELT                  
C                   IS IGNORED.                                                 
C*******************************************************************************
      REAL LS,LF                                                                
      COMMON/DHEAT/DHS                                                          
C      COMMON/IFU,UMSEC,FUCOEF,HEIGHT                                     
      DIMENSION D(100),P(100),Z(100),AQI(100),TC(100),DTCDP(100)                
C      DT=3600.0*DELTAT                                                          
C Keep DT in seconds
      DT = DELTAT
      E=0.99
      WRITE(*,*)'SURFAC: ABOVE VSTART'
      CALL VSTART(QIR,QI,QR,QA,SIGMA,DT,CW,LS,PW,PA,GAMMA,THEDA,                
     1 THEDAM,C1,LF,C2,D_O,X,C3,C4,E,RW,ITS) 
      CALL ZDEPTH(N,D,Z,TDEPTH)                                                 
      CALL ABSORB(N,QI,QR,Z,D,P,AQI,EXTINC,CV,G1,G2,G3)                         
      CALL COEFF(N,P,TC,DTCDP,COEFKE)   
      WRITE(*,*)'SURFAC: BELOW COEFF'                                          
      TW=0.0     
      IF (PX.GT.0.0) CALL WTBULB(TW,TA,EA,PA)                                   
      IF (ITOPT.EQ.1) GO TO 125                                                 
      TT4=TST*TST*TST*TST                                                       
      EOT=C1*EXP(-6141.9/TST)                                                   
      DL4=D(N)*0.5                                                              
      DCDT=(EOT/(4615.0*TST*TST))*((LS/(RW*TST))-1.0)                           
      DSUR=((D_O*1000.0)/PA)*((TST/273.16)**X)                                   
      FST=TC(1)+LS*DSUR*DCDT                                                    
      CIT=C3+C4*TST                                                             
      CALL GROUND(TGT,C1,TCG,LS,DCG,FGT,DCDTG)                                  
C*******************************************************************************
C     BEGIN ITERATIVE LOOP                                                      
      DO 100 LOOP=1,10                                                          
      ITER=LOOP-1                                                               
      EOTDT=C1*EXP(-6141.9/TSTDT)   
      WRITE(*,*)'TGT',TGT,'TSTDT',TSTDT
      DCDT=(EOTDT/(4615.0*TSTDT*TSTDT))*((LS/(RW*TSTDT))-1.0)                   
      DSUR=((D_O*1000.0)/PA)*((TSTDT/273.16)**X)                                 
      FSTDT=TC(1)+LS*DSUR*DCDT                                                  
      CITDT=C3+C4*TSTDT                                                         
      CALL GROUND(TGTDT,C1,TCG,LS,DCG,FGTDT,DCDTG)                              
      IF (IFU.GT.0) CALL WINDF(DELTAT,UMSEC,TA,TST,TSTDT,FUCOEF,PA,             
     1HEIGHT,IFU,FU,RICRIT,ZO)                                                 
      TTDT2=TSTDT*TSTDT                                                         
      TTDT3=TTDT2*TSTDT                                                         
      TTDT4=TTDT3*TSTDT                                                         
C*******************************************************************************
C     COMPUTE RESIDUAL.                                                         
      IF (IGUESS.EQ.1) GO TO 101                                                
      DQ=0.0                                                                    
      IF ((TST.GT.273.159).AND.(TSTDT.GT.273.159)) GO TO 102                    
C     QG WHEN BOTTOM LAYER IS NOT AT ZERO DEGREES CELSIUS.                      
      T1=(THEDAM*DT*FST*FGT*(TGT-TST))/(FGT*DL4+FST*DTG)                        
      T2=(THEDA*DT*FSTDT*FGTDT*(TGTDT-TSTDT))/(FGTDT*DL4+FSTDT*DTG)             
      QG=T1+T2                                                                  
      GO TO 105                                                                 
C     QG WHEN BOTTOM LAYER IS AT ZERO DEGREES CELSIUS.                          
  102 T3=THEDAM*DT*FGT*(TGT-TST)/DTG                                            
      T4=THEDA*DT*FGTDT*(TGTDT-TSTDT)/DTG                                       
      QG=T3+T4                                                                  
      GO TO 105                                                                 
C     QG WHEN MAKING A FIRST GUESS.                                             
  101 QG=0.0                                                                    
      DQ=DHS                                                                    
  105 DH=QIR -THEDAM*E*SIGMA*TT4-THEDA*E*SIGMA*TTDT4                            
     1+LS*PW*FU*GAMMA*TA-THEDAM*LS*PW*FU*GAMMA*TST-THEDA                        
     2*LS*PW*FU*GAMMA*TSTDT+LS*PW*FU*EA-THEDAM*LS*PW*FU*EOT-                    
     3THEDA*LS*PW*FU*EOTDT+CW*PW*PX*TW-CW*PW*PX*273.16+QG-DQ                    
      RESID=D(1)*P(1)*CITDT*TSTDT-D(1)*P(1)*CIT*TST+LF                          
     1*PW*WSTDT-LF*PW*WST-DH    
      WRITE(*,*)'DH RESID',DH,RESID
C*******************************************************************************
C     COMPUTE THE PARTIAL DERIVATIVE FOR THE SURFACE LAYER                      
      IF (NOKNOW.EQ.1) GO TO 106                                                
C     TEMPERATURE IS THE UNKNOWN. WITH ONLY ONE LAYER, THAT                     
C        LAYER IS BOTH THE SURFACE LAYER AND THE BOTTOM LAYER.                  
      IF (IGUESS.EQ.1) GO TO 107                                                
      IF ((TST.GT.273.159).AND.(TSTDT.GT.273.159)) GO TO 108                    
C     BOTTOM LAYER NOT AT ZERO DEGREES CELSIUS.                                 
      DCDT2=(EOTDT/(4615.0*TTDT3))*(((LS/(RW*TSTDT))**2)                        
     1-((4.0*LS)/(RW*TSTDT))+2.0)                                               
      DVN=((TSTDT**X)*DCDT2)+(X*DCDT*(TSTDT**(X-1.0)))                          
      DEODTN=((THEDA*DT*FGTDT)/(FGTDT*DL4+FSTDT*DTG))*                          
     1((DTG*LS*FSTDT*DVN*(TGTDT-TSTDT))/(FGTDT*DL4+                             
     2FSTDT*DTG)-LS*C2*DVN*(TGTDT-TSTDT)+FSTDT)                                 
      GO TO 110                                                                 
C     BOTTOM LAYER AT ZERO DEGREES CELSIUS.                                     
  108 DEODTN=(THEDA*DT*FGTDT)/DTG                                               
      GO TO 110                                                                 
  107 DEODTN=0.0                                                                
  110 DEDU=D(1)*P(1)*C3+2.0*D(1)*P(1)*C4*TSTDT+4.0*THEDA*E                      
     1*SIGMA*TTDT3+THEDA*LS*PW*FU*GAMMA+(THEDA*LS*LS*PW*FU*                     
     2EOTDT)/(RW*TTDT2)+DEODTN                                                  
      GO TO 115                                                                 
C     LIQUID-WATER UNKNOWN                                                      
  106 DEDU=LF*PW        
      WRITE(*,*)'DEDU LF PW',DEDU,LF,PW
C*******************************************************************************
C     COMPUTE THE CORRECTION                                                    
  115 DU=-RESID/DEDU     
      WRITE(*,*)'DU RESID DEDU',DU,RESID,DEDU
C*******************************************************************************
C     ADD CORRECTION AND CHECK TOLERANCE                                        
      NOVER=0                                                                   
      IF (NOKNOW.EQ.1) GO TO 116                                                
C     TEMPERATURE UNKNOWN                                                       
      ABSDU=ABS(DU)                                                             
      IF (ABSDU.GT.TOLER)NOVER=1                                                
      TSTDT=TSTDT+DU                                                            
      IF (TSTDT.LE.273.16) GO TO 120                                            
C     TEMPERATURE EXCEEDS ZERO CELSIUS-LIQUID-WATER NOW UNKNOWN                 
      EXCESS=TSTDT-273.16                                                       
      TSTDT=273.16                                                              
      NOKNOW=1                                                                  
      WSTDT=((D(1)*P(1)*CITDT)/(LF*PW))*EXCESS                                  
      GO TO 120                                                                 
C     LIQUID-WATER UNKNOWN                                                      
  116 WTOL=((D(1)*P(1)*CITDT)/(LF*PW))*TOLER   
      WRITE(*,*)'D P CITDT',D(1),P(1),CITDT
      IF (WTOL.LT.0.0001) WTOL=0.0001                                           
      ABSDU=ABS(DU)
      WRITE(*,*)'ABSDU WTOL',ABSDU,WTOL
      IF (ABSDU.GT.WTOL) NOVER=1                                                
      WSTDT=WSTDT+DU                                                            
      IF (WSTDT.GE.0.0) GO TO 120                                               
C     LIQUID-WATER NEGATIVE - TEMPERATURE NOW UNKNOWN.                          
      EXCESS=WSTDT                                                              
      WSTDT=0.0                                                                 
      NOKNOW=0                                                                  
      TSTDT=273.16                                                              
  120 CONTINUE                                                                  
C     CORRECTION HAS BEEN MADE AND TOLERANCE CHECKED.                           
C*******************************************************************************
C     CHECK FOR TERMINATION                                                     
      IF (NOVER.EQ.0) RETURN                                                    
  100 CONTINUE                                                                  
C     END OF ITERATION LOOP                                                     
C*******************************************************************************
      IF (NOKNOW.EQ.1) GO TO 124                                                
      IF ((IGUESS.EQ.1).AND.(ABSDU.LT.0.5)) RETURN                              
  124 I=1                                                                       
      PRINT 900,MONTH,IDAY,IYEAR,IHOUR,TOLER,IGUESS                             
  900 FORMAT(1H1,28HSOLUTION DID NOT CONVERGE ON,I3,1H/,                        
     1I2,1H/,I4,2X,5HHOUR=,I2,1H.,5X,12HTOLERANCE IS,F5.2,                      
     21X,14HDEGREES KELVIN,5X,7HIGUESS=,I1)                                     
      PRINT 902                                                                 
  902 FORMAT(1H0,5HLAYER,4X,5HDEPTH,2X,9HTHICKNESS,3X,                          
     17HDENSITY,5X,3HAQI,4X,8HCONDUCT.,5X,6HNOKNOW,8X,2HTT,                     
     26X,4HTTDT,10X,2HWT,8X,4HWTDT,2X,8HRESIDUAL,2X,10HCORRECTION)               
      PRINT 901,I,Z(1),D(1),P(1),AQI(1),TC(1),NOKNOW,                           
     1TST,TSTDT,WST,WSTDT,RESID,DU                                              
  901 FORMAT (1H ,I5,2X,4F9.4,3x,F10.5,I10,4(1XF10.3),2F10.4)                            
      STOP                                                                      
C*******************************************************************************
C     ITOPT=1,WILL ONLY OCCUR WHEN IGUESS=0.                                    
  125 IF (TSTDT.LT.273.159) GO TO 126                                           
      TAVG=(TST+TSTDT)*0.5                                                      
      EO=C1*EXP(-6141.9/TAVG)                                                   
      HEAT=AQI(1)+E*QA-E*SIGMA*TAVG*TAVG*TAVG*TAVG+LS*FU*                       
     1(GAMMA*(TA-TAVG)+(EA-EO))+CW*PW*PX*(TW-273.16)                            
      IF (HEAT.LT.0.0) HEAT=0.0                                                 
      DELTAW=HEAT/(LF*PW)                                                       
      WSTDT=WSTDT+DELTAW                                                        
      GO TO 130                                                                 
  126 WSTDT=0.0                                                                 
  130 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

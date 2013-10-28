      SUBROUTINE PATERN(N,TTDT,WTDT,ITS,TSMAX,CYCLE,D,P,IDAY,IHOUR,             
     1ITER,DO,TI,SURDEN)                                                        
C*******************************************************************************
C     PERPARES PRINTED OUTPUT WHEN A SELECTED SURFACE TEMPERATURE               
C        PATTERN IS BEING USED.                                                 
C*******************************************************************************
      REAL LF                                                                   
      DIMENSION TTDT(100),WTDT(100),D(100),P(100),GEI1(7),GEI2(31),             
     1Z(100),TC(100),DTCDP(100),CI(100)                                         
C     GAUSS ERROR INTEGRAL 0.0 TO 0.10                                          
      DATA GEI1/0.0,0.0226,0.0451,0.0676,0.0901,0.1125,0.1349/                  
C     GAUSS ERROR INTEGRAL 0.0 TO 3.0                                           
      DATA GEI2/0.0,0.1125,0.2227,0.3286,0.4284,0.5205,0.6039,0.6778,           
     10.7421,0.7969,0.8427,0.8802,0.9103,0.9340,0.9523,0.9661,0.9763,           
     20.9827,0.9891,0.9922,0.9953,0.9967,0.9981,0.9986,0.9991,0.9996,           
     30.9997,0.9998,0.9999,1.0,1.0/                                             
C*******************************************************************************
      LAPSE=(IDAY-1)*24+IHOUR                                                   
      ETIME=LAPSE                                                               
      IF (ITS.GT.2) GO TO 100                                                   
      PRINT 900,TSMAX                                                           
  900 FORMAT (1H1,10X,23HSELECTED PATTERN OUTPUT,10X,                           
     123HINSTANTANEOUS CHANGE OF, F6.1,1X,16HDEGREES CELSIUS.)                  
      GO TO 101                                                                 
  100 PRINT 901,CYCLE,TSMAX                                                     
  901 FORMAT (1H1,10X,23HSELECTED PATTERN OUTPUT,10X,                           
     120HSIN-VARIATION-CYCLE=,F4.1,1X,16HDAYS--AMPLITUDE=,F5.1,1X,              
     216HDEFRGES CELSIUS.)                                                      
  101 PRINT 902, LAPSE,ITER,DO                                                  
  902 FORMAT(1H0,13HELAPSED TIME=,I3,1X,5HHOURS,5X,                             
     120HREQUIRED ITERATIONS=,I2,5X,25HDIFFUSION COEFFICIENT-DO=,F5.3,          
     21X,8HCM2/SEC.)                                                            
      PRINT 909                                                                 
  909 FORMAT (1H ,5X,105HTHE THEORETICAL TEMPERATURE IS FOR A SNOW COVER        
     1 WITH UNIFORM DENSITY AND A DIFFUSION COEFFICIENT OF ZERO.)               
      PRINT 903                                                                 
  903 FORMAT (1H0,12X,8HDEPTH TO,40X,5HWATER,10X,5HTOTAL,19X,                   
     111HTHEORETICAL)                                                           
      PRINT 904                                                                 
  904 FORMAT (1H ,5HLAYER,6X,9HMID-POINT,6X,9HTHICKNESS,3X,                     
     112HLIQUID-WATER,5X,10HEQUIVALENT,8X,7HDENSITY,4X,11HTEMPERATURE,          
     24X,11HTEMPERATURE,  5X,10HDIFFERENCE)                                     
      PRINT 905                                                                 
  905 FORMAT(1H ,17X,3HCM.,12X,3HCM.,12X,3HMM.,12X,3HMM.,21X,                   
     19HDEGREES C,6X,9HDEGREES C,6X,9HDEGREES C)                                
C*******************************************************************************
C     OBTAIN VALUES NEEDED FOR COMPUTATIONS.                                    
      CALL ZDEPTH(N,D,Z,TDEPTH)                                                 
C     THE THEORETICAL TEMPERATURE IS COMPUTED ASSUMING UNIFORM DENSITY.         
C         (TOTAL SNOW COVER DENSITY IS ASSUMED EQUAL TO THE INITIAL SURFACE     
C         LAYER DENSITY.)                                                       
      NL=1                                                                      
      PSAVE=P(1)                                                                
      P(1)=SURDEN                                                               
      CALL COEFF(NL,P,TC,DTCDP)                                                 
      P(1)=PSAVE                                                                
C      SPECFIC HEAT=0.5 FOR SELECTED PATTERN CASES.                             
      C3=0.5                                                                    
      C4=0.0                                                                    
      LF=79.7                                                                   
      PW=1.0                                                                    
      CIO=0.5                                                                   
      CALL SPHEAT (N,TTDT,C3,C4,CI)                                             
C     INITIALIZE SNOW-COVER WATER-EQUIVALENT AND HEAT DEFICIT.                  
      WE=0.0                                                                    
      HEAT=0.0                                                                  
      SUMD=0.0                                                                  
C*******************************************************************************
C     PRINT SURFACE LAYER                                                       
      I=1                                                                       
      WATER=WTDT(1)*10.0                                                        
      WEL=(P(1)*D(1)+WTDT(1))*10.0                                              
      WE=WE+WEL                                                                 
      DEN=(WEL*0.1)/D(1)                                                        
      TEMP=TTDT(1)-273.16                                                       
      THEO=TEMP                                                                 
      DIFF=0.0                                                                  
      AVGCI=(CIO+CI(1))*0.5                                                     
      HEAT=HEAT-TEMP*AVGCI*P(1)*D(1)                                            
      PRINT 906,I,Z(1),D(1),WATER,WEL,DEN,TEMP,THEO,DIFF                        
  906 FORMAT(1H ,I5,2F15.3,6F15.2)                                              
      TEMP1=TEMP                                                                
      ET=ETIME*3600.0                                                           
      FREQ=1.0/(CYCLE*24.0*3600.0)                                              
C*******************************************************************************
C     COMPUTE AND PRINT OTHER LAYERS                                            
      DO 110 I=2,N                                                              
      ZZ=Z(I)-Z(1)                                                              
      WATER=WTDT(I)*10.0                                                        
      WEL=(P(I)*D(I)+WTDT(I))*10.0                                              
      WE=WE+WEL                                                                 
      DEN=(WEL*0.1)/D(I)                                                        
      TEMP=TTDT(I)-273.16                                                       
      AVGCI=(CIO+CI(I))*0.5                                                     
      HEAT=HEAT-TEMP*AVGCI*P(I)*D(I)-WTDT(I)*LF*PW                              
C     COMPUTE THEORETICAL TEMPERATURE                                           
      IF(ITS.GT.2) GO TO 115                                                    
C     INSTANTANEOUS SURFACE CHANGE--UNIFORM INITIAL TEMPERATURE.                
      X=ZZ/(2.0*SQRT((TC(1)*ET)/(CI(I)*SURDEN)))                                
      IF(X.GT.3.0)GO TO 111                                                     
      IF (X.GT.0.10) GO TO 112                                                  
      IX=(X*50.0)+1.0                                                           
      FRAC=(X-((IX-1)*0.02))/0.02                                               
      GEI=GEI1(IX)+FRAC*(GEI1(IX+1)-GEI1(IX))                                   
      GO TO 114                                                                 
  112 IX=(X*10.0)+1.0                                                           
      FRAC=(X-((IX-1)*0.1))/0.1                                                 
      GEI=GEI2(IX)+FRAC*(GEI2(IX+1)-GEI2(IX))                                   
      GO TO 114                                                                 
  111 GEI=1.0                                                                   
  114 THEO=TEMP1-TSMAX*GEI                                                      
      GO TO 116                                                                 
C     SIN-VARIATION--UNIFORM INITIAL TEMPERATURE.                               
  115 X=ZZ*SQRT(3.1416*FREQ*CI(I)*SURDEN/TC(1))                                 
      ANGLE=6.2832*FREQ*ET-X                                                    
      IF (ANGLE.GT.0.0) GO TO 117                                               
      THEO=TI                                                                   
      GO TO 116                                                                 
  117 THEO=TI+TSMAX*EXP(-X)*SIN(ANGLE)                                          
  116 DIFF=TEMP-THEO                                                            
      SUMD=SUMD+ABS(DIFF)                                                       
      PRINT 906,I,Z(I),D(I),WATER,WEL,DEN,TEMP,THEO,DIFF                        
  110 CONTINUE                                                                  
C*******************************************************************************
C     COMPUTE AND PRINT TOTAL SNOW COVER VALUES.                                
      DEN=(WE*0.1)/TDEPTH                                                       
      PRINT 908,SUMD                                                            
  908 FORMAT (1H ,92X,28HSUM OF ABSOLUTE DIFFERENCES=,F5.2)                     
      PRINT 907,TDEPTH,WE,DEN,HEAT                                              
  907 FORMAT(1H0,24HTOTAL SNOW COVER--DEPTH=,F5.1,1X,3HCM.,                     
     15X,17HWATER-EQUIVALENT=,F5.0,1X,3HMM.,5X,8HDENSITY=,                      
     2F3.2,5X,13HHEAT DEFICIT=,F5.1,1X,8HCAL/CM2.)                              
      RETURN                                                                    
      END                                                                       

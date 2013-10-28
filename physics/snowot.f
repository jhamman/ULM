      SUBROUTINE SNOWOT (N,TTDT,WTDT,D,P,MONTH,IDAY,IYEAR,IHOUR,                
     1ITER,WE,TDEPTH)                                                           
C*******************************************************************************
C     PRINTS THE STATE OF THE SNOW COVER AT THE END OF EACH OUTPUT INTERVAL.    
C*******************************************************************************
C 2009-Feb-15 Removed print statements                    Ben Livneh
C
      REAL LF                                                                   
      COMMON/WLAGS/WLAG(11),NLAG,STORE                                          
      DIMENSION TTDT(100),WTDT(100),D(100),P(100),Z(100),CI(100)                
C      PRINT 900,MONTH,IDAY,IYEAR,IHOUR,ITER                                     
C  900 FORMAT (1H1,10X,26HSTATE OF THE SNOW COVER ON,I3,1H/,I2,1H/,              
C     1I4,3X,5HHOUR=,I2,20X,I2,1X,34HITERATIONS DURING THE LAST PERIOD.)         
C      PRINT 901                                                                 
C  901 FORMAT (1H0,7X,18HDEPTH TO MID-POINT,35X,5HWATER,10X,5HTOTAL)             
C      PRINT 902                                                                 
C  902 FORMAT (1H ,5HLAYER,6X,9HFROM SOIL,6X,9HTHICKNESS,3X,                     
C     112HLIQUID-WATER,5X,10HEQUIVALENT,8X,7HDENSITY,4X,11HTEMPERATURE)          
C      PRINT 903                                                                 
C  903 FORMAT (1H ,17X,3HCM.,12X,3HCM.,12X,3HMM.,12X,3HMM.,                      
C     121X,9HDEGREES C)                                                          
C*******************************************************************************
C     OBTAIN VALUES NEEDED FOR COMPUTATIONS.                                    
      CALL ZDEPTH(N,D,Z,TDEPTH)                                                 
      DO 102 I=1,N                                                              
  102 Z(I)=TDEPTH-Z(I)                                                          
      C3=0.0222384                                                              
      C4=0.00176                                                                
      LF=79.7                                                                   
      PW=1.0                                                                    
      CIO=0.503                                                                 
      CALL SPHEAT(N,TTDT,C3,C4,CI)                                              
C     INITIALIZE SNOW COVER WATER-EQUIVALENT AND HEAT DEFICIT.                  
      WE=0.0                                                                    
      HEAT=0.0                                                                  
C*******************************************************************************
C     PRINT EACH LAYER                                                          
      DO 100 I=1,N                                                              
      WATER=WTDT(I)*10.0                                                        
      WEL=(P(I)*D(I)+WTDT(I))*10.0                                              
      WE=WE+WEL                                                                 
      DEN=(WEL*0.1)/D(I)                                                        
      TEMP=TTDT(I)-273.16                                                       
      AVGCI=(CIO+CI(I))*0.5                                                     
      HEAT=HEAT-TEMP*AVGCI*P(I)*D(I)                                            
C      PRINT 904,I,Z(I),D(I),WATER,WEL,DEN,TEMP                                  
C  904 FORMAT (1H ,I5,2F15.3,2F15.2,F15.3,F15.2)                                 
  100 CONTINUE                                                                  
C*******************************************************************************
C     COMPUTE AND PRINT TOTAL SNOW COVER VALUES.                                
      DO 101 I=1,NLAG                                                           
  101 WE=WE+WLAG(I)*10.0                                                        
      WE=WE+STORE*10.0                                                          
      DEN=(WE*0.1)/TDEPTH                                                       
C      PRINT 905,TDEPTH,WE,DEN,HEAT                                              
C  905 FORMAT (1H0,24HTOTAL SNOW COVER--DEPTH=,F5.1,1X,3HCM.,                    
C     15X,17HWATER-EQUIVALENT=,F5.0,1X,3HMM.,5X,8HDENSITY=,                      
C     2F4.3,5X,13HHEAT DEFICIT=,F5.1,1X,8HCAL/CM2.)                              
      RETURN                                                                    
      END                                                                       

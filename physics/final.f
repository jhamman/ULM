      SUBROUTINE FINAL(IMO,IDA,IYR,LMO,LDA,LYR,WEI,SUMPX,WE)                    
C*******************************************************************************
C     THIS SUBROUTINE COMPUTES AND PRINT A RUN SUMMARY, INCLUDING A             
C        WATER BALANCE.                                                         
C*******************************************************************************
      COMMON/STATS/NCPO,OPO,SPO,SPO2,OSPO,OSPO2,PODIFF,NCTO,OTO,                
     1STO,STO2,OSTO,OSTO2,TODIFF,VAPOR,WATER,VSOIL,OPO2,OTO2,                   
     2 TQS,TQLW,TQH,TQE,TQG,TQPX,IRF(10),ITF(10)                                
      COMMON/SSVT/SUMSVT                                                        
      DIMENSION RL(9),TL(9)                                                     
      DATA RL/-0.2,-0.02,0.0,0.01,0.02,0.05,0.1,0.2,0.4/                        
      DATA TL/-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0/                          
      PRINT 912                                                                 
  912 FORMAT (1H1)                                                              
      DO 100 I=1,6                                                              
  100 PRINT 913                                                                 
  913 FORMAT (1H0)                                                              
C      PRINT 900,IMO,IDA,IYR,LMO,LDA,LYR                                         
C  900 FORMAT(1H0,40X,28HFINAL SUMMARY FOR THE PERIOD,I3,1H/,I2,1H/,I4,          
     11X,2HTO,I3,1H/,I2,1H/,I4,1H.)                                             
      IF((NCPO.EQ.0).AND.(NCTO.EQ.0))GO TO 110                                  
      IF(NCTO.EQ.0)GO TO 115                                                    
C*******************************************************************************
C     SNOW SURFACE TEMPERATURE STATISTICS.                                      
C      PRINT 911                                                                 
C  911 FORMAT (1H0,53HSNOW SURFACE TEMPERATURE--(UNITS ARE DEGREES CELSIU        
     1S))                                                                       
      PRINT 901                                                                 
  901 FORMAT (1H0,5HCASES,3X,9HOBS. MEAN,3X,9HSIM. MEAN,                        
     1             8X,4HBIAS,3X,9HRMS ERROR,5X,15HAVG. ABS. ERROR,7X,           
     2   23HCORRELATION COEFFICIENT,7X,13HBEST FIT LINE)                        
      CC=1.0                                                                    
      A=0.0                                                                     
      B=1.0                                                                     
      AVG=OTO/NCTO                                                              
      AVGS=STO/NCTO                                                             
      BIAS=AVGS-AVG                                                             
      AVEABS=TODIFF/NCTO                                                        
      RMS=SQRT(OSTO2/NCTO)                                                      
      IF(AVGS.GT.-0.05)B=99.99                                                  
      IF(AVG.GT.-0.05)B=0.0                                                     
      IF(B.NE.1.0)GO TO 111                                                     
      A=(OTO*STO2-STO*OSTO)/(NCTO*STO2-STO*STO)                                 
      B=(NCTO*OSTO-STO*OTO)/(NCTO*STO2-STO*STO)                                 
      D=SQRT((NCTO*STO2-STO*STO)*(NCTO*OTO2-OTO*OTO))                           
      CC=(NCTO*OSTO-STO*OTO)/D                                                  
  111 PRINT 902,NCTO,AVG,AVGS,BIAS,RMS,AVEABS,CC,A,B                            
  902 FORMAT (1H ,I5,4F12.1,F17.1,F28.3,                                        
     16X,4HOBS=,F5.1,1H+,F5.2,4H*SIM)                                           
  115 IF(NCPO.EQ.0)GO TO 110                                                    
C*******************************************************************************
C     SNOW COVER OUTFLOW STATISTICS.                                            
      PRINT 903                                                                 
  903 FORMAT (1H0,43HSNOW COVER OUTFLOW--(UNITS ARE MILLIMETERS))               
      PRINT 901                                                                 
      CC=1.0                                                                    
      A=0.0                                                                     
      B=1.0                                                                     
      AVG=OPO/NCPO                                                              
      AVGS=SPO/NCPO                                                             
      BIAS=AVGS-AVG                                                             
      AVEABS=PODIFF/NCPO                                                        
      RMS=SQRT(OSPO2/NCPO)                                                      
      IF(AVGS.LT.0.05)B=99.99                                                   
      IF(AVG.LT.0.05)B=0.0                                                      
      IF(B.NE.1.0)GO TO 116                                                     
      A=(OPO*SPO2-SPO*OSPO)/(NCPO*SPO2-SPO*SPO)                                 
      B=(NCPO*OSPO-SPO*OPO)/(NCPO*SPO2-SPO*SPO)                                 
      D=SQRT((NCPO*SPO2-SPO*SPO)*(NCPO*OPO2-OPO*OPO))                           
      CC=(NCPO*OSPO-SPO*OPO)/D                                                  
  116 PRINT 914,NCPO,AVG,AVGS,BIAS,RMS,AVEABS,CC,A,B                            
  914 FORMAT (1H ,I5,4F12.2,F17.2,F28.3,                                        
     16X,4HOBS=,F5.2,1H+,F5.2,4H*SIM)                                           
C*******************************************************************************
C     WATER BALANCE.                                                            
  110 WE=WE*0.1                                                                 
      BAL=WEI-WE+VAPOR+VSOIL+SUMPX-WATER                                        
      PRINT 904                                                                 
  904 FORMAT(1H0,46HWATER BALANCE SUMMARY--(UNITS ARE CENTIMETERS))             
      PRINT 905                                                                 
  905 FORMAT(1H0,14X,16HWATER-EQUIVALENT,11X,14HVAPOR TRANSFER,25X,             
     110HSNOW COVER)                                                            
      PRINT 906                                                                 
  906 FORMAT(1H ,8X,7HINITIAL,10X,5HFINAL,7X,8HAIR-SNOW,6X,9HSOIL-SNOW,         
     12X,13HPRECIPITATION,8X,7HOUTFLOW,8X,7HBALANCE)                            
      PRINT 907,WEI,WE,VAPOR,SUMSVT,SUMPX,WATER,BAL                             
  907 FORMAT(1H ,7F15.2)                                                        
C*******************************************************************************
C     ENERGY BALANCE COMPONENTS.                                                
      PRINT 908                                                                 
  908 FORMAT (1H0,46HENERGY BALANCE COMPONENTS--(UNITS ARE CAL/CM2))            
      PRINT 909                                                                 
  909 FORMAT (1H0,1X,19HSHORTWAVE RADIATION,2X,18HLONGWAVE RADIATION,           
     17X,13HSENSIBLE HEAT,9X, 11HLATENT HEAT,6X,14HSOIL-SNOW HEAT,              
     26X,14HHEAT FROM RAIN)                                                     
      PRINT 910,TQS,TQLW,TQH,TQE,TQG,TQPX                                       
  910 FORMAT (1H ,6F20.0)                                                       
C*******************************************************************************
C     FREQUENCY TABLES.                                                         
      DO 120 I=1,2                                                              
  120 PRINT 913                                                                 
      PRINT 915                                                                 
  915 FORMAT (1H0,28H***** FREQUENCY TABLES *****)                              
      PRINT 916                                                                 
  916 FORMAT (1H0,18HRICHARDSON NUMBER.)                                        
      PRINT 917,RL                                                              
  917 FORMAT (1H ,16HINTERVAL LIMITS.,9X,9F10.2)                                
      PRINT 918,IRF                                                             
  918 FORMAT (1H ,16HNUMBER OF CASES.,4X,10I10)                                 
      PRINT 919                                                                 
  919 FORMAT (1H0,50HAIR TEMPERATURE MINUS SNOW SURFACE TEMPERATURE (C),        
     15X,42HAIR TEMPERATURE .LE. ZERO DEGREES CELSIUS.)                         
      PRINT 917,TL                                                              
      PRINT 918,ITF                                                             
C*******************************************************************************
      RETURN                                                                    
      END                                                                       

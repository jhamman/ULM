      SUBROUTINE WINDF(DELTAT,UA,TA,TST,TSTDT,FUCOEF,PA,HEIGHT,IFU,FU,
     &RICRIT,ZO)          
C*******************************************************************************
C     COMPUTES THE WIND FUNCTION FOR THE PERIOD--(CM/MB)                        
C          UA IS AVERAGE WIND--METERS/SEC                                       
C          FUCOEF IS IN--(MM/MB/KM)                                             
C          HEIGHT IS IN METERS                                                  
C          ZO IS IN CENTIMETERS.                                                
C*******************************************************************************
C      COMMON/CRITRI/RICRIT,ZO                                                   
      DIMENSION RI(50),CWR(50),RATIO(50),RIBULK(50)                             
      DATA IFIRST/0/                                                            
      IF (IFIRST.GT.0) GO TO 100                                                
C     MIMIMUM WIND FUNCTION--MOLECULAR CONDUCTION.                              
C      FUMIN=(0.646*DELTAT)/(PA*HEIGHT*100.0)                                    
C Convert DELTAT to hours for emperical formula
      FUMIN=(0.646*(DELTAT/3600))/(PA*HEIGHT*100.0)              
      IFIRST=1                                                                  
      IF(IFU.LT.1)GO TO 100                                                     
C*******************************************************************************
C     GENERATE TABLE OF BULK TRANSFER COEFFICIENT RATIO FOR                     
C        WATER VAPOR AND HEAT TO ITS VALUE FOR NEUTRAL                          
C        CONDITIONS AS A FUNCTION OF THE BULK RICHARDSON NUMBER.                
C     COMPUTE NEUTRAL BULK TRANSFER COEFFICIENT                                 
      CDN=0.16/((ALOG(HEIGHT*100.0/ZO))**2)                                     
C     GENERATE TABLE                                                            
      DO 110 J=1,50                                                             
      ZL=-((J*0.1)**2)                                                          
      X=(1.0-16.0*ZL)**0.25                                                     
      TERM=1.0-(SQRT(CDN)/0.4)*(ALOG((1.0+X*X)/2.0)                             
     1     +2.0*ALOG((1.0+X)/2.0)-2.0*ATAN(X)+1.5708)                           
      CDR=1.0/(TERM*TERM)                                                       
      CHR=SQRT(CDR)/(1.0-5.0*SQRT(CDN)*ALOG((1+X*X)/2.0))                       
      RIB=(ZL*SQRT(CDN)*(CDR**1.5))/(0.4*CHR)                                   
      RI(J)=RIB                                                                 
      CWR(J)=CHR                                                                
  110 CONTINUE                                                                  
C     PUT TABLE INTO COMPUTABLE RI INCREMENTS.                                  
      DO 115 J=1,50                                                             
      RIB=-((J*0.03)**2)                                                        
      RIBULK(J)=RIB                                                             
      DO 116 JJ=1,50                                                            
      INC=JJ                                                                    
      IF(RIB.GT.RI(JJ))GO TO 120                                                
  116 CONTINUE                                                                  
      INC=51                                                                    
  120 IF(INC.GT.1)GO TO 122                                                     
      RATIO(J)=((RIB/RI(INC))*(CWR(INC)-1.0))+1.0                               
      GO TO 115                                                                 
  122 IF(INC.LE.50)GO TO 125                                                    
      RATIO(J)=(((RIB-RI(49))/(RI(50)-RI(49)))*                                 
     1         (CWR(50)-CWR(49)))+CWR(49)                                       
      GO TO 115                                                                 
  125 RATIO(J)=(((RIB-RI(INC-1))/(RI(INC)-RI(INC-1)))*                          
     1         (CWR(INC)-CWR(INC-1)))+CWR(INC-1)                                
  115 CONTINUE                                                                  
      CTN=FUCOEF                                                                
C*******************************************************************************
C     COMPUTE WIND TRAVEL IN KILOMETERS. 
C Convert kilometers to allow DELTAT to be in seconds                                       
C  100 UT=3.6*UA*DELTAT                                                          
  100 UT=(UA*DELTAT)/1000
      IF (IFU.EQ.1) GO TO 101                                                   
C*******************************************************************************
C     USE EMPIRICAL WIND FUNCTION                                               
      FU=0.1*FUCOEF*UT                                                          
      RETURN                                                                    
C*******************************************************************************
C     USE THEORETICAL WIND FUNCTION WITH A STABILITY CORRECTION.                
  101 IF (UA.GT.0.20) GO TO 102                                                 
C     CALM CONDITIONS                                                           
      FU=FUMIN                                                                  
      RETURN                                                                    
C     COMPUTE THE RICHARDSON NUMBER.                                            
  102 TS=0.5*(TST+TSTDT)                                                        
      RIB=(2.0*9.8*(TA-TS)*HEIGHT)/((TA+TS)*UA*UA)                              
      IF(RIB.GE.0.0)GO TO 103                                                   
C     UNSTABLE CONDITIONS.                                                      
      J=(SQRT(-RIB)/0.03)+1.0                                                   
      IF((J.EQ.1).OR.(J.GT.50))GO TO 106                                        
      R=(((RIB-RIBULK(J-1))/(RIBULK(J)-RIBULK(J-1)))                            
     1        *(RATIO(J)-RATIO(J-1)))+RATIO(J-1)                                
      GO TO 107                                                                 
  106 IF(J.GT.50)GO TO 108                                                      
      R=((RIB/RIBULK(1))*(RATIO(1)-1.0))+1.0                                    
      GO TO 107                                                                 
  108 R=(((RIB-RIBULK(49))/(RIBULK(50)-RIBULK(49)))                             
     1        *(RATIO(50)-RATIO(49)))+RATIO(49)                                 
  107 CT=R*CTN                                                                  
      GO TO 105                                                                 
C     NEUTRAL OR STABLE CONDITIONS.                                             
  103 IF(RIB.LT.RICRIT)GO TO 104                                                
      CT=0.0                                                                    
      GO TO 105                                                                 
  104 R=(1.0-(RIB/RICRIT))**2                                                   
      CT=R*CTN                                                                  
C     COMPUTE WIND FUNCTION.                                                    
  105 FU=0.1*CT*UT                                                              
      IF(FU.LT.FUMIN)FU=FUMIN                                                   
C*******************************************************************************
      RETURN                                                                    
      END                                                                       

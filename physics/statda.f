      SUBROUTINE STATDA(TST,TSTDT,TO,SCOUT,PO,WE,WESC,WEPW,VAPOUR,              
     1TDEPTH,DEPTH,STAKE,MONTH,IDAY,IYEAR,IHOUR,DELTAT,N,SOILVT,TA)             
C*******************************************************************************
C     THIS SUBROUTINE COMPARES VERIFICATION DATA BY.                            
C          1. STORING SURFACE TEMPERATURE,OUTFLOW AND VAPOR TRANSFER            
C                  FOR DISPLAY ONCE PER DAY.                                    
C          2. COMPUTING RUN STATISTICS BETWEEN OBSERVED AND COMPUTED            
C                  SURFACE TEMPERATURE AND OUTFLOW.                             
C          3. PRINTING A DAILY SUMMARY OF WATER-EQUIVALENT,                     
C                  AND DEPTH WHEN IHOUR=24.                                     
C*******************************************************************************
      COMMON/STATS/NCPO,OPO,SPO,SPO2,OSPO,OSPO2,PODIFF,NCTO,OTO,                
     1STO,STO2,OSTO,OSTO2,TODIFF,VAPOR,WATER,VSOIL,OPO2,OTO2,                   
     2 TQS,TQLW,TQH,TQE,TQG,TQPX,IRF(10),ITF(10)                                
      COMMON/EBAL/QS,QLW,QH,QE,QG,QPX,DQ                                        
      COMMON/IFU,UMSEC,FUCOEF,HEIGHT                                     
      COMMON/STEMP/TOSIM                                                        
C*******************************************************************************
      DIMENSION SPOA(24),POA(24),STOA(24),TOA(24),VT(24)                        
      DIMENSION QSA(24),QLWA(24),QHA(24),QEA(24),QGA(24),QPXA(24),              
     1DQA(24),RI(24),RL(9),TL(9)                                                
      DATA IFIRST/0/                                                            
      DATA RL/-0.2,-0.02,0.0,0.01,0.02,0.05,0.1,0.2,0.4/                        
      DATA TL/-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0/                          
      IF (IFIRST.GT.0) GO TO 100                                                
C     INITIALIZE STAT VALUES.                                                   
      NCPO=0                                                                    
      OPO=0.0                                                                   
      SPO=0.0                                                                   
      SPO2=0.0                                                                  
      OPO2=0.0                                                                  
      OSPO=0.0                                                                  
      OSPO2=0.0                                                                 
      PODIFF=0.0                                                                
      NCTO=0                                                                    
      OTO=0.0                                                                   
      STO=0.0                                                                   
      STO2=0.0                                                                  
      OTO2=0.0                                                                  
      OSTO=0.0                                                                  
      OSTO2=0.0                                                                 
      TODIFF=0.0                                                                
      VAPOR=0.0                                                                 
      VSOIL=0.0                                                                 
      WATER=0.0                                                                 
      TQS=0.0                                                                   
      TQLW=0.0                                                                  
      TQH=0.0                                                                   
      TQE=0.0                                                                   
      TQG=0.0                                                                   
      TQPX=0.0                                                                  
      NP=0                                                                      
      DO 116 J=1,10                                                             
      IRF(J)=0                                                                  
  116 ITF(J)=0                                                                  
      IFIRST=1                                                                  
C  100 IDT=DELTAT+0.01  
C Remove time control                                                         
  100 CONTINUE                                                           
C*******************************************************************************
C     CHECK FOR THE FIRST PERIOD OF A DAY -- INITIALIZE SUMS.                   
      IF (IHOUR.GT.IDT) GO TO 112                                               
      IF (N.EQ.0) RETURN                                                        
      SUMVT=0.0                                                                 
      SUMSPO=0.0                                                                
      SUMPO=0.0                                                                 
      SUMQS=0.0                                                                 
      SUMQLW=0.0                                                                
      SUMQH=0.0                                                                 
      SUMQE=0.0                                                                 
      SUMQG=0.0                                                                 
      SUMQPX=0.0                                                                
      SUMDQ=0.0                                                                 
      MPO=0                                                                     
      MTO=0                                                                     
      GO TO 104                                                                 
C     DO NOT COMPUTE VALUES IF A SNOW COVER DID NOT EXIST ON THE                
C         FIRST RERIOD OF THE DAY.  (RUN TOTALS MUST BE COMPUTED.)              
  112 IF ((NP.EQ.0).AND.(N.EQ.0)) RETURN                                        
      IF (NP.EQ.0) GO TO 105                                                    
      IF (N.EQ.0) GO TO 115                                                     
C*******************************************************************************
C     STORE VALUES FOR PRINTING AT END OF DAY.                                  
  104 I=IHOUR/IDT                                                               
      IF (N.LT.0) N=0                                                           
      NP=NP+1                                                                   
      IF (NP.EQ.I) GO TO 130                                                    
C     FILL IN INTERMEDIATE VALUES ON DAY WHEN SNOW DISAPPEARS AND REAPPEARS.    
      JJ=I-1                                                                    
      DO 131 J=NP,JJ                                                            
      SPOA(J)=999.9                                                             
      POA(J)=999.9                                                              
      STOA(J)=9999.                                                             
      TOA(J)=9999.                                                              
      VT(J)=99.99                                                               
      QSA(J)=999.9                                                              
      QLWA(J)=999.9                                                             
      QHA(J)=999.9                                                              
      QEA(J)=999.9                                                              
      QGA(J)=999.9                                                              
      QPXA(J)=999.9                                                             
      RI(J)=99.99                                                               
  131 DQA(J)=999.9                                                              
      NP=I                                                                      
  130 SPOA(I)=SCOUT                                                             
      POA(I)=PO                                                                 
      TS=(TST+TSTDT)*0.5                                                        
      IF (N.GT.1) TS=TOSIM                                                      
C     COMPUTE THE RICHARDSON NUMBER.                                            
      UA=UMSEC                                                                  
      IF (UA.LT.0.01) UA=0.01                                                   
      RIN=(2.0*9.8*((TA-TS)/HEIGHT))/((TA+TS)*((UA/HEIGHT)**2))                 
      IF (RIN.LT.-2.0) RIN=-2.0                                                 
      IF (RIN.GT.2.0) RIN=2.0                                                   
      RI(I)=RIN                                                                 
C     FILL FREQUENCY TABLES.                                                    
C          RICHARDSON NUMBER.                                                   
      DO 132 J=1,9                                                              
      IF (RIN.GT.RL(J)) GO TO 132                                               
      IRF(J)=IRF(J)+1                                                           
      GO TO 135                                                                 
  132 CONTINUE                                                                  
      IRF(10)=IRF(10)+1                                                         
C     AIR TEMP. MINUS SNOW SURFACE TEMP.                                        
  135 TDIFF=TA-TS                                                               
      IF (TA.GT.273.16) GO TO 140                                               
      DO 136 J=1,9                                                              
      IF(TDIFF.GT.TL(J)) GO TO 136                                              
      ITF(J)=ITF(J)+1                                                           
      GO TO 140                                                                 
  136 CONTINUE                                                                  
      ITF(10)=ITF(10)+1                                                         
  140 CONTINUE                                                                  
      TS=TS-273.16                                                              
      IF (N.EQ.0) TS=0.0                                                        
      STOA(I)=TS                                                                
      TOA(I)=TO                                                                 
      VT(I)=VAPOUR*10.0                                                         
      IF(N.GT.0)GO TO 120                                                       
      VAPOUR=0.0                                                                
      VT(I)=0.0                                                                 
      QSA(I)=0.0                                                                
      QLWA(I)=0.0                                                               
      QHA(I)=0.0                                                                
      QEA(I)=0.0                                                                
      QGA(I)=0.0                                                                
      QPXA(I)=0.0                                                               
      DQA(I)=0.0                                                                
      QS=0.0                                                                    
      QLW=0.0                                                                   
      QH=0.0                                                                    
      QE=0.0                                                                    
      QG=0.0                                                                    
      QPX=0.0                                                                   
      DQ=0.0                                                                    
      GO TO 121                                                                 
  120 QSA(I)=QS                                                                 
      QLWA(I)=QLW                                                               
      QHA(I)=QH                                                                 
      QEA(I)=QE                                                                 
      QGA(I)=QG                                                                 
      QPXA(I)=QPX                                                               
      DQA(I)=DQ                                                                 
  121 SUMQS=SUMQS+QS                                                            
      SUMQLW=SUMQLW+QLW                                                         
      SUMQH=SUMQH+QH                                                            
      SUMQE=SUMQE+QE                                                            
      SUMQG=SUMQG+QG                                                            
      SUMQPX=SUMQPX+QPX                                                         
      SUMDQ=SUMDQ+DQ                                                            
C*******************************************************************************
C     COMPUTE TOTALS AND ADD TO STATISTICS.                                     
      SUMVT=SUMVT+VT(I)                                                         
      SUMSPO=SUMSPO+SCOUT                                                       
      IF (PO.LT.9000.0) GO TO 101                                               
      MPO=MPO+1                                                                 
      GO TO 102                                                                 
C     SNOW COVER OUTFLOW STATISTICS                                             
  101 SUMPO=SUMPO+PO                                                            
      NCPO=NCPO+1                                                               
      OPO=OPO+PO                                                                
      SPO=SPO+SCOUT                                                             
      SPO2=SPO2+SCOUT*SCOUT                                                     
      OPO2=OPO2+PO*PO                                                           
      OSPO=OSPO+PO*SCOUT                                                        
      OSPO2=OSPO2+(PO-SCOUT)*(PO-SCOUT)                                         
      PODIFF=PODIFF+ABS(PO-SCOUT)                                               
  102 IF (TO.LT.9000.0) GO TO 103                                               
      MTO=MTO+1                                                                 
      GO TO 105                                                                 
C     SNOW SURFACE TEMPERATURE STATISTICS.                                      
  103 NCTO=NCTO+1                                                               
      OTO=OTO+TO                                                                
      STO=STO+TS                                                                
      STO2=STO2+TS*TS                                                           
      OTO2=OTO2+TO*TO                                                           
      OSTO=OSTO+TO*TS                                                           
      OSTO2=OSTO2+(TO-TS)*(TO-TS)                                               
      TODIFF=TODIFF+ABS(TO-TS)                                                  
C     RUN TOTALS                                                                
  105 VAPOR=VAPOR+VAPOUR                                                        
      WATER=WATER+SCOUT*0.1                                                     
      VSOIL=VSOIL+SOILVT                                                        
      TQS=TQS+QS                                                                
      TQLW=TQLW+QLW                                                             
      TQH=TQH+QH                                                                
      TQE=TQE+QE                                                                
      TQG=TQG+QG                                                                
      TQPX=TQPX+QPX                                                             
      IF (NP.EQ.0) RETURN                                                       
  115 IF (IHOUR.LT.24)RETURN                                                    
C*******************************************************************************
C     END OF DAY SUMMARY.                                                       
      IF (N.EQ.0) PRINT 926                                                     
  926 FORMAT (1H1)                                                              
      PRINT 900                                                                 
  900 FORMAT (1H0,120H**************************************************        
     1******************************************************************        
     2****)                                                                     
      PRINT 901,MONTH,IDAY,IYEAR                                                
  901 FORMAT (1H0,15HDAILY SUMMARY--,I2,1H/,I2,1H/,I4,10X,                      
     134H(ALL NINES INDICATES MISSING DATA))                                    
      PRINT 902,((I),I=IDT,IHOUR,IDT)                                           
  902 FORMAT (1H0,4HHOUR,1X,24I5)                                               
      PRINT 903                                                                 
  903 FORMAT (1H+,125X,2X,5HTOTAL)                                              
C     SNOW COVER OUTFLOW.                                                       
      PRINT 904                                                                 
  904 FORMAT (1H ,22HSNOW COVER OUTFLOW-MM.)                                    
      PRINT 905,(SPOA(I),I=1,NP)                                                
  905 FORMAT(1H ,4HSIM.,1X,24F5.1)                                              
      PRINT 906,SUMSPO                                                          
  906 FORMAT (1H+,125X,F7.1)                                                    
      IF (MPO.GE.NP) GO TO 106                                                  
      IF (MPO.GT.0)SUMPO=9999.9                                                 
      DO 113 I=1,NP                                                             
      IF (POA(I).LT.9000.0) GO TO 113                                           
      POA(I)=999.9                                                              
  113 CONTINUE                                                                  
      PRINT 907,(POA(I),I=1,NP)                                                 
  907 FORMAT (1H ,4HOBS.,1X,24F5.1)                                             
      PRINT 906,SUMPO                                                           
C     SNOW SURFACE TEMPERATURE.                                                 
  106 PRINT 908                                                                 
  908 FORMAT (1H ,31HSNOW SURFACE TEMPERATURE-DEG.C.)                           
      PRINT 915,(STOA(I),I=1,NP)                                                
  915 FORMAT(1H ,4HSIM.,1X,24F5.0)                                              
      IF (MTO.GE.NP) GO TO 107                                                  
      PRINT 916,(TOA(I),I=1,NP)                                                 
  916 FORMAT (1H ,4HOBS.,1X,24F5.0)                                             
C     VAPOR TRANSFER.                                                           
  107 PRINT 909                                                                 
  909 FORMAT (1H ,27HAIR-SNOW VAPOR TRANSFER-MM.)                               
      PRINT 910,(VT(I),I=1,NP)                                                  
  910 FORMAT (1H ,4HSIM.,1X,24F5.2)                                             
      PRINT 911,SUMVT                                                           
  911 FORMAT (1H+,125X,F7.2)                                                    
      PRINT 917                                                                 
  917 FORMAT(1H ,34HENERGY BALANCE COMPONENTS-CAL/CM2.)                         
      PRINT 918                                                                 
  918 FORMAT(1H ,2HQS)                                                          
      PRINT 919,(QSA(I),I=1,NP)                                                 
  919 FORMAT(1H+,5X,24F5.1)                                                     
      PRINT 906,SUMQS                                                           
      PRINT 920                                                                 
  920 FORMAT(1H ,3HQLW)                                                         
      PRINT 919,(QLWA(I),I=1,NP)                                                
      PRINT 906,SUMQLW                                                          
      PRINT 921                                                                 
  921 FORMAT(1H ,2HQH)                                                          
      PRINT 919,(QHA(I),I=1,NP)                                                 
      PRINT 906,SUMQH                                                           
      PRINT 922                                                                 
  922 FORMAT(1H ,2HQE)                                                          
      PRINT 919,(QEA(I),I=1,NP)                                                 
      PRINT 906,SUMQE                                                           
      PRINT 923                                                                 
  923 FORMAT(1H ,2HQG)                                                          
      PRINT 919,(QGA(I),I=1,NP)                                                 
      PRINT 906,SUMQG                                                           
      PRINT 924                                                                 
  924 FORMAT(1H ,3HQPX)                                                         
      PRINT 919,(QPXA(I),I=1,NP)                                                
      PRINT 906,SUMQPX                                                          
      PRINT 925                                                                 
  925 FORMAT(1H ,2HDQ)                                                          
      PRINT 919,(DQA(I),I=1,NP)                                                 
      PRINT 906,SUMDQ                                                           
      PRINT 927                                                                 
  927 FORMAT (1H ,2HRI)                                                         
      PRINT 928,(RI(I),I=1,NP)                                                  
  928 FORMAT (1H+,5X,24F5.2)                                                    
C     DEPTH AND WATER-EQUIVALENT.                                               
      IF (TDEPTH.EQ.0.0) GO TO 109                                              
      SDEN=(WE*0.1)/TDEPTH                                                      
      GO TO 111                                                                 
  109 SDEN=0.0                                                                  
  111 IF (WESC.GT.9000.0) GO TO 108                                             
      IF (DEPTH.LT.0.001) GO TO 108                                             
      ODEN=(WESC*0.1)/DEPTH                                                     
      GO TO 110                                                                 
  108 ODEN=99.99                                                                
  110 PRINT 912                                                                 
  912 FORMAT (1H0,7X,15HSIMULATED(2400),11X,24HSNOW COURSE(MID-MORNING),        
     16X,12HPILLOW(2400),3X,23HSNOW STAKE(MID-MORNING))                         
      PRINT 913                                                                 
  913 FORMAT (1H ,4X,6HWE(MM),1X,9HDEPTH(CM),3X,7HDENSITY,                      
     14X,6HWE(MM),1X,9HDEPTH(CM),3X,7HDENSITY,9X,6HWE(MM),                      
     26X,9HDEPTH(CM))                                                           
      PRINT 914,WE,TDEPTH,SDEN,WESC,DEPTH,ODEN,WEPW,STAKE                       
  914 FORMAT (1H ,2F10.0,F10.3,2F10.0,F10.3,2F15.0)                             
      NP=0                                                                      
C*******************************************************************************
      RETURN                                                                    
      END                                                                       

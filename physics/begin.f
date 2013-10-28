      SUBROUTINE BEGIN(N,D,P,TT,WT,NOKNOW,ITS,TI,WEI,THICK0,DEN,TEMP,
     &WATER0)
C*******************************************************************************
C     INITIALIZATION OF SNOW COVER                                              
C*******************************************************************************
C     INPUT CARDS ARE AS FOLLOWS                                                
C   FIRST CARD --
C     N             I5,TOTAL NUMBER OF LAYERS(100 IS THE MAXIMUM)               
C     MM            I5  =1, LIQUID-WATER IN MILLIMETERS.                        
C                       =0, LIQUID-WATER IN PERCENT(DECIMAL)                    
C     REMAINING CARDS                                                           
C  NN    I5    LAYER NUMBER (NN)                                                
C  THICK F5.0  THICKNESS OF THE LAYER IN CENTIMETERS                            
C  DEN   F5.2  DENSITY OF THE LAYER--DECIMAL. (ICE PLUS LIQUID)                 
C  TEMP  F5.0  MEAN TEMPERATURE IN DEGREES CELSIUS FOR THE LAYER                
C  WATER F5.2  LIQUID-WATER CONTENT OF THE LAYER--PERCENT OR MILLIMETERS.       
C     NOTE.... A CARD IS NEEDED WHENEVER INITIAL CONDITIONS CHANGE.             
C              THAT IS,EACH CARD DEFINES INITIAL VALUES FOR LAYER               
C              NP+1 TO NN, WHERE NP IS THE LAYER NUMBER ON THE                  
C              PREVIOUS CARD. THE LAYER NUMBER ON THE LAST CARD MUST            
C              MATCH WITH THE TOTAL NUMBER OF LAYERS.                           
C*******************************************************************************
C Initialize water equivalent and separate out ice and liquid water
C Remove cyclical calls (103) to read multiple snow layers as this
C subroutine will only be invoked for a new pack (1 layer)

      COMMON/DEBUG/IDEBUG,IDBUG,IHR,LHR,LISTF                                   
      DIMENSION D(100),P(100),TT(100),WT(100),NOKNOW(100)                       
      NP=0                                                                      
C     READ 900,N,MM,IDEBUG,IDBUG,IHR,LHR,LISTF                                  
C 900 FORMAT(2I5,45X,5I5)    
C Liquid water input in millimeters (should be zero), for single layer
      MM = 0.
      N = 1
      NN = 1
      IF (IDEBUG.LT.1) IDEBUG=0                                                 
C  103 READ 901,NN,THICK,DEN,TEMP,WATER                                          
C  901 FORMAT (I5,F5.0,F5.2,F5.0,F5.2)                                           
      IF (ITS.GT.1) TEMP=TI                                                     
      IF (NN.GT.N) STOP                                                         
      NS=NP+1                                                                   
      DO 100 I=NS,NN                                                            
      D(I)=THICK0      
C TWEL: TOTAL WATER EQUIVALENT LAYER: ICE AND LIQUID 
C SEPARATE LIQUID WATER FOR DENSITY COMPUTATION
      TWEL=D(I)*DEN                                                             
      IF (MM.EQ.1) GO TO 107                                                    
      WEL=TWEL/(1.0+WATER0)
      WT(I)=WATER0*WEL
      GO TO 108                                                                 
  107 WEL=TWEL-WATER0*0.1
      WT(I)=WATER0*0.1
  108 P(I)=WEL/D(I)                                                             
      TT(I)=TEMP+273.16  
      WRITE(*,*)'BEGIN TT(I)',TT(I)
      IF (TEMP.LT.0.0) GO TO 101                                                
      IF (WATER0.GT.0.0) GO TO 102
C     LAYER AT ZERO DEGREES CELSIUS AND NO LIQUID-WATER,                        
C        MAKE TEMPERATURE THE UNKNOWN.                                          
      NOKNOW(I)=0                                                               
      GO TO 100                                                                 
C     TEMPERATURE UNKNOWN                                                       
  101 NOKNOW(I)=0                                                               
      WT(I)=0.0                                                                 
      GO TO 100                                                                 
C     LIQUID-WATER IS THE UNKNOWN                                               
  102 NOKNOW(I)=1                                                               
      TT(I)=273.16                                                              
  100 CONTINUE                                                                  
      IF (NN.EQ.N) GO TO 105                                                    
      NP=NN                                                                     
C      GO TO 103                                         
C     SET NOKNOW VALUES TO NEGATIVE TO INDICATE NEW LAYER.                      
  105 DO 106 I=1,N                                                              
      NOKNOW(I)=NOKNOW(I)-2                                                     
  106 CONTINUE                                                                  
C     COMPUTE INITIAL SNOW COVER WATER-EQUIVALENT.                              
      WEI=0.0                                                                   
      DO 110 I=1,N                                                              
      WEL=P(I)*D(I)                                                             
      WEI=WEI+WEL+WT(I)                                                         
  110 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

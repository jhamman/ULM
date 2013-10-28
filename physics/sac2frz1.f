C  SUBROUTINE SPLITS FREE AND TENSION WATER STORAGES OF SAC-SMA
C  INTO TOTAL WATER CONTENTS OF FROZEN GROUND MODEL SOIL LAYERS
CVK 12/2005  Added parameter SWLT to the subroutine parameters
CBL Two major bug fixes: (1) when sh2o falls below swlt, corrected
CBL the accounting for the amount of water being added, and (2)
CBL do not let SH2O fall below SWLT (previously SMC only), since
CBL the algebra was not written to consider this case and becomes 
CBL backwards which, just like case (1) resulted in the model
CBL essentially *creating* water = IMBALANCE!

CVK  12/2005      SUBROUTINE SAC2FRZ1(DWT,DWF,SMC,SH2O,NUPL,NLOWL,ZSOIL,SMAX)
      SUBROUTINE SAC2FRZ1(DWT,DWF,SMC,SH2O,NUPL,NLOWL,ZSOIL,SMAX,SWLT,
     &     PRFLAG)
     
C   DWT - TENSION WATER CHANGE PER TIME INTERVAL, MM
C   DWF - FREE WATER CHANGE PER TIME INTERVAL, MM
C   SMC - TOTAL MOISTURE CONTENT OF FROZEN GROUND MODEL LAYERS
C   SH2O - UNFROZEN WATER CONTENT OF FROZEN GROUND MODEL LAYERS
C   NUPL - UPPER SOIL LAYER TO DISTRIBUTE SAC-SMA WATER
C   NLOWL - LOWER SOIL LAYER TO DISTRIBUTE SAC-SMA WATER
C   ZSOIL - SOIL LAYER DEPTHS, M
C   SMAX - SOIL POROSITY
CBL -- upper layer and wilting point bugs fixed - 3/2010 -- Ben Livneh

      real ZSOIL(*),SMC(*),SH2O(*)
      double precision SMCX(5),SH2OX(5)
c      double precision dsf,dz,z,dwx,dwf0,dwt0
      INTEGER MS(10)
      INTEGER PRFLAG
      double precision stot0,sliq0,s0,s1,frost,stot1,sliq1
      double precision  diftot,difliq,stot2,sliq2,d1,d2,dz,z,dsf

C  SPLIT FREE WATER BETWEEN SOIL LAYERS
C  FREE WATER IS SPLITTED EQUALLY RETWEEN LAYERS 

      DO I = 1,5
         SMCX(I) = 0.
         SH2OX(I) = 0.
      ENDDO
      if(prflag==1)write(*,*)'top SAC2FRZ'
      if(prflag==1)write(*,*)'i stot frost sliq s0 s1 smc sh2o'
      stot0 = 0.
      sliq0 = 0.
      frost = 0.
      s0 = 0.
      s1 = 0.
      if(prflag==1)write(*,*)'swlt smax',swlt,smax
      do i = nupl,nlowl
         SMCX(I) = SMC(I)
         SH2OX(I) = SH2O(I)
         if(i==1)then
            dz = 0.0 - zsoil(1)
         else
            dz = zsoil(i-1) - zsoil(i)
         endif
         stot0 = stot0+1000*(smcx(i) - swlt) * dz
         sliq0 = sliq0+1000*(sh2ox(i) - swlt) * dz
         frost = frost + 1000*(smcx(i) - sh2ox(i)) * dz
         s0 = 1000 * (smcx(i) - swlt) * dz
         s1 = 1000 * (sh2ox(i) - swlt) * dz
         if(prflag==1) write(*,"(IX,8(xf13.8))")i,stot0,frost,sliq0,
     +        s0,s1,smcx(i),sh2ox(i),dz
      enddo
      if(prflag==1)write(*,*)'DWF',dwf
      if(nupl .gt. 5) then
         write(*,*) 'nup',nupl,nlowl,dwt,dwf,(zsoil(i),i=1,5)
         stop
      endif 
      DO I=NUPL,NLOWL
         MS(I)=1
         if(prflag==1)write(*,*)'zsoil',i,zsoil(i)
      ENDDO 
cbl Correction so top-soil layer is not a null reference zsoil(0)
      if(nupl-1==0) then
         Z = 0.0 - ZSOIL(NLOWL)
      else
         Z = ZSOIL(NUPL-1) - ZSOIL(NLOWL)
      endif
      DSF=0.001*DWF/Z
      if(prflag==1)write(*,*)'z dsf',z,dsf
      M=0
 100  NN=0
      M=M+1
      IF(M .LE. NLOWL-NUPL+1) THEN
         DWX=0.
         IF(M .GT. 1) THEN
            Z=0.
            DO II=NUPL,NLOWL
               if(II-1==0) then
                  Z=Z+(0.0 - ZSOIL(II))*MS(II)
               else
                  Z=Z+(ZSOIL(II-1)-ZSOIL(II))*MS(II)
               endif
            ENDDO
         ENDIF  
         DO I=NUPL,NLOWL
            IF(MS(I) .NE. 0) THEN
               if(prflag==1)write(*,*)'smc dsf',i,smcx(i),dsf
               SMCX(I)=SMCX(I)+DSF          
               if(prflag==1)write(*,*)'smc dsf',i,smcx(i),dsf
               if(prflag==1)write(*,*)'sh2o dsf',i,sh2ox(i),dsf
               SH2OX(I)=SH2OX(I)+DSF
               if(prflag==1)write(*,*)'sh2o dsf',i,sh2ox(i),dsf
               DELTA=SMCX(I)-SMAX
CVK  12/2005       IF(DELTA .GT. 0.0001 .OR. SMC(I) .LT. 0.) THEN
cbl DELTA is excess over SMAX
               if(prflag==1)write(*,*)'DELTA',DELTA
cbl               IF(DELTA .GT. 0.0001 .OR. SMC(I) .LT. SWLT) THEN
cbl need to prevent SH2O from falling below SWLT here. 
cbl also prevent delta from getting > 0.0 (moisture above smax)
               IF(DELTA .GT. 0.0 .OR. SH2OX(I) .LT. SWLT) THEN
                  NN=I
                  MS(I)=0
                  if (i-1==0) then
                     DZ=0.0 - ZSOIL(I)
                  else
                      DZ=ZSOIL(I-1)-ZSOIL(I)
                  endif
                  if(prflag==1)write(*,*)'DZ',DZ
                  Z=Z-DZ
                  if(prflag==1)write(*,*)'new Z',z
cbl remove current layer from total depth of interest and spread
cbl the excess (delta)to remaining depth 'Z'
CVK  12/2005        IF(SMC(I) .GE. 0.) THEN
cbl                  IF(SMC(I) .GE. SWLT) THEN
                  IF(SH2OX(I) .GE. SWLT) THEN
                     DWX=DWX+DELTA*DZ*1000.
                     SH2OX(I)=SH2OX(I)-DELTA
                     SMCX(I)=SMAX
                     if(prflag==1)write(*,*)'DWX',DWX
                     if(prflag==1)write(*,*)'smc=smax',i,smcx(i),smax
                  ELSE
cbl wb error explained using example: smc(1) falls below swlt
cbl This appears to be an error in logic which is likely an artifact of
cbl the previous code where smc(i) was compared with zero and not swlt
cbl as it is done in the present code.
cbl In effect, what is happening is that when (eg) smc(1) falls below
cbl swlt, then it is set back to swlt. However rather than tallying the 
cbl remaining deficit (smc(1) - swlt) which is equal to the amount of water
cbl added, and applying it to smc(2), the amount smc(1) before being reset
cbl to swlt is added to the smc(2), which *creates* an amount of water 
cbl equal to swlt, because we add smc(1) to smc(2) and then add swlt-smc(1)
cbl to smc(1) when it is set to swlt.  DWX carries the deficit to smc(2)
CBL                     DWX=DWX+SMC(I)*DZ*1000.
cbl first fixt                     DWX=DWX+(SMC(I)-SWLT)*DZ*1000.
CVK  12/2005         SMC(I)=0.
CVK  12/2005         SH2O(I)=0.
cbl                     SMC(I)=SWLT
cbl                     SH2O(I)=SWLT
cbl In this secondary fix, the order is important, so that we can add back the
cbl same amount to SMC(I) that was removed, which brought SH2O(I) below SWLT
cbl as we bring SH2O(I) back to SWLT and carry the remainder in DWX
                     SMCX(I)=SMCX(I)+(SWLT-SH2OX(I))
                     DWX=DWX+(SH2OX(I)-SWLT)*DZ*1000.   
                     SH2OX(I)=SWLT
                     if(prflag==1)write(*,*)'smc=swlt',i,smc(i),swlt
                  ENDIF
                  if(prflag==1) write(*,*)'dwx',dwx
               ENDIF 
            ENDIF
         ENDDO

         IF(NN .NE. 0 .AND. Z .GT. 0.0001) THEN
            DSF=0.001*DWX/Z
            if(prflag==1)write(*,*)'dsf=dwx/z',dsf,dwx        
            GOTO 100
         ENDIF 
      ELSE
         WRITE(*,*) ' WARN: NO BALANCE IN SAC2FRZ',M,NUPL,NLOWL,DSF,
     +        DWF,NN,(MS(I),ZSOIL(I),SMCX(I),I=NUPL,NLOWL)
         WRITE(*,*) ' WARN: NO BALANCE IN SAC2FRZ',M,NLOWL,DSF,DWF,DWX
      ENDIF 

      if(prflag==1)write(*,*)'end of DWF',DWF
      stot1 = 0.
      sliq1 = 0.
      frost = 0.
      s0 = 0.
      s1 = 0.
      if(prflag==1)write(*,*)'swlt smax',swlt,smax
      do i = nupl,nlowl
         if(i==1)then
            dz = 0.0 - zsoil(1)
         else
            dz = zsoil(i-1) - zsoil(i)
         endif
         stot1 = stot1 + 1000 * (smcx(i) - swlt) * dz
         sliq1 = sliq1 + 1000 * (sh2ox(i) - swlt) * dz
         frost = frost + 1000 * (smcx(i) - sh2ox(i)) * dz
         s0 = 1000 * (smcx(i) - swlt) * dz
         s1 = 1000 * (sh2ox(i) - swlt) * dz
         if(prflag==1) write(*,"(IX,8(xf13.8))")i,stot1,frost,sliq1,
     +        s0,s1,smcx(i),sh2ox(i),dz
      enddo
      diftot = stot1 - stot0
      difliq = sliq1 - sliq0
      if(prflag==1)write(*,*)'diftot1 difliq1',diftot,difliq,DWF
      d1=diftot-dwf
      d2=difliq-dwf
      if(prflag==1)write(*,*)'# error',d1,d2,dwf

C  SPLIT TENSION WATER: IF TENSION WATER REDUCTION, SPLIT BY A RATIO
C  OF UNFROZEN WATER FROM PREVIOUS TIME STEP; IF TENSION WATER INCREASE,
C  SPLIT BY AN INVERSE RATIO OF TOTAL WATER DEFICIT.
C
C  CALCULATE A RATIO OF UNFROZEN WATER OR RATIO OF TOTAL WATER DEFICIT
      N=1
77    SAVG=0.
      if(prflag==1)write(*,*)'DWT',dwt
      DO I=NUPL,NLOWL
       IF(DWT .LT. 0.) THEN
CVK  12/2005        SAVG=SAVG+SH2OX(I)
        SAVG=SAVG+(SH2OX(I)-SWLT)
        if(prflag==1)write(*,*)'savg',i,sh2ox(i)
       ELSE
        SAVG=SAVG+(SMAX-SMCX(I))
        if(prflag==1)write(*,*)'savg',i,smcx(i)
       ENDIF   
      ENDDO
      if(prflag==1)write(*,*)'savg',savg
      SAVG=SAVG/(NLOWL-NUPL+1)
      if(prflag==1)write(*,*)'savg',savg
CVK  12/2005    IF(SAVG .LT. 1E-5) GOTO 7
      IF(SAVG .LT. 1E-6) GOTO 7
      
      ALP=0.
      DO I=NUPL,NLOWL
         if(i-1==0)then
            DZ=0.0 - ZSOIL(I)
         else
            DZ=ZSOIL(I-1)-ZSOIL(I)
         endif
         if(prflag==1)write(*,*)'DZ',DZ,I
       IF(DWT .LT. 0.) THEN
C  UNFROZEN WATER RATIO
CVK  12/2005         ALP=ALP+DZ*SH2O(I)/SAVG
        ALP=ALP+DZ*(SH2OX(I)-SWLT)/SAVG
        if(prflag==1)write(*,*)'ALP',ALP
       ELSE
C  TOTAL WATER DEFICIT RATIO
        ALP=ALP+DZ*(SMAX-SMCX(I))/SAVG
        if(prflag==1)write(*,*)'ALP',ALP
       ENDIF
      ENDDO
      ALP=1./ALP
      if(prflag==1)write(*,*)'1/=ALP',ALP

C  RUN REDISTRIBUTION OF WATER BETWEEN SOIL LAYERS
      DDTX=0.
      DO I=NUPL,NLOWL
         if(i-1==0)then
            DZ= 0.0 - ZSOIL(I)
         else
            DZ=ZSOIL(I-1)-ZSOIL(I)
         endif
         if(prflag==1)write(*,*)'DZ',DZ,I
       IF(DWT .LT. 0.) THEN
C  REDUCTION IN TENSION WATER. USE A RATIO OF UNFROZEN WATER
CVK  12/2005      DMAX=1000.*SH2O(I)*DZ
        DMAX=1000.*(SH2OX(I)-SWLT)*DZ
        if(prflag==1)write(*,*)'DMAX',DMAX,I
CVK  12/2005        DZTX=DZ*(DWT*(SH2O(I)/SAVG)*ALP+
        if(i-1==0)then
        DZTX=DZ*(DWT*((SH2OX(I)-SWLT)/SAVG)*ALP+ 
     +       DDTX/( 0.0 - ZSOIL(NLOWL)))
        else 
        DZTX=DZ*(DWT*((SH2OX(I)-SWLT)/SAVG)*ALP+ 
     +       DDTX/(ZSOIL(I-1)-ZSOIL(NLOWL)))
        endif
        if(prflag==1)write(*,*)'DZTX',DZTX
        xx1=ABS(DZTX)
        if(prflag==1)write(*,*)'XX1 DMAX',XX1,DMAX
        IF(xx1 .GT. DMAX) THEN
         DDTX=DZTX+DMAX
         if(prflag==1)write(*,*)'DDTX',DDTX
         DZTX=-DMAX
         if(prflag==1)write(*,*)'DZTX',DZTX
        ENDIF 
       ELSE
C  INCREASE IN TENSION WATER. USE AN INVERSE RATIO OF UNFROZEN WATER
        DMAX=1000.*(SMAX-SMCX(I))*DZ
        if(prflag==1)write(*,*)'dmax smc',i,dmax,smcx(i)
        if(i-1==0)then
        DZTX=DZ*(DWT*((SMAX-SMCX(I))/SAVG)*ALP+
     +       DDTX/( 0.0 - ZSOIL(NLOWL)))   
        else
        DZTX=DZ*(DWT*((SMAX-SMCX(I))/SAVG)*ALP+
     +       DDTX/(ZSOIL(I-1)-ZSOIL(NLOWL)))       
        endif
        if(prflag==1)write(*,*)'dztx',i,dztx
        IF(DZTX .GT. DMAX) THEN
         DDTX=DZTX-DMAX
         if(prflag==1)write(*,*)'ddtx',ddtx
         DZTX=DMAX
         if(prflag==1)write(*,*)'dztx>dmax',dztx,ddtx
        ENDIF 
       ENDIF
       
       A=0.001*DZTX/DZ
       if(prflag==1)write(*,*)'A DZTX',A,DZTX
       if(prflag==1)write(*,*)'smc sh2o',smcx(i),sh2ox(i)
       SH2OX(I)=SH2OX(I)+A
       SMCX(I)=SMCX(I)+A
       if(prflag==1)write(*,*)'smc+=A',smcx(i),A
       if(prflag==1)write(*,*)'sh2o+=A',sh2ox(i),A
CVK  12/2005     IF(SH2O(I) .LT. 0.) SH2O(I)=0.
CVK  12/2005     IF(SMC(I) .LT. 0.) SMC(I)=0.
       IF(SH2OX(I) .LT. SWLT) SH2OX(I)=SWLT
       IF(SMCX(I) .LT. SWLT) THEN
          if(prflag==1)write(*,*)'wilting breach',smcx(i),sh2ox(i),swlt
        SMCX(I)=SWLT
        SH2OX(I)=SWLT
       ENDIF 
      ENDDO

      xx1=ABS(DDTX)
      if(prflag==1)write(*,*)'xx1',i,xx1
CVK  12/2005    IF(xx1 .GT. 1E-2) THEN
      IF(xx1 .GT. 1E-4) THEN
       if(prflag==1)write(*,*)'xx1.gt.1e-4'
       DWT=DDTX
       if(prflag==1)write(*,*)'dwt-ddtx',dwt
       N=N+1
       if(prflag==1)write(*,*)'N',N
       IF(N .GT. 2) THEN
        WRITE(*,*) ' WARN: NO BALANCE IN SAC-SMA WATER CHANGE.'
       ELSE
        GOTO 77
       ENDIF 
      ENDIF 
      if(prflag==1)write(*,*)'bot1 SAC2FRZ'
      stot2 = 0.
      sliq2 = 0.
      frost = 0.
      s0 = 0.
      s1 = 0.
      if(prflag==1)write(*,*)'swlt smax',swlt,smax
      do i = nupl,nlowl
         if(i==1)then
            dz = 0.0 - zsoil(1)
         else
            dz = zsoil(i-1) - zsoil(i)
         endif
         stot2 = stot2 + 1000 * (smcx(i) - swlt) * dz
         sliq2 = sliq2 + 1000 * (sh2ox(i) - swlt) * dz
         frost = frost + 1000 * (smcx(i) - sh2ox(i)) * dz
         s0 = 1000 * (smcx(i) - swlt) * dz
         s1 = 1000 * (sh2ox(i) - swlt) * dz
         if(prflag==1) write(*,"(IX,8(xf13.8))")i,stot2,frost,sliq2,
     +        s0,s1,smcx(i),sh2ox(i),dz
      enddo
      diftot = stot2 - stot1
      difliq = sliq2 - sliq1
      if(prflag==1)write(*,*)'diftot2 difliq2',diftot,difliq,DWT
      d1=diftot-dwt
      d2=difliq-dwt
      if(prflag==1)write(*,*)'# error',d1,d2,dwt
CBL Convert precision back
      DO I = NUPL, NLOWL
         SMC(I) = SMCX(I)
         SH2O(I) = SH2OX(I)
      ENDDO
      RETURN

7     CONTINUE
      if(prflag==1)write(*,*)'bot2 SAC2FRZ'
      stot2 = 0.
      sliq2 = 0.
      frost = 0.
      s0 = 0.
      s1 = 0.
      if(prflag==1)write(*,*)'swlt smax',swlt,smax
      do i = nupl,nlowl
         if(i==1)then
            dz = 0.0 - zsoil(1)
         else
            dz = zsoil(i-1) - zsoil(i)
         endif
         stot2 = stot2 + 1000 * (smcx(i) - swlt) * dz
         sliq2 = sliq2 + 1000 * (sh2ox(i) - swlt) * dz
         frost = frost + 1000 * (smcx(i) - sh2ox(i)) * dz
         s0 = 1000 * (smcx(i) - swlt) * dz
         s1 = 1000 * (sh2ox(i) - swlt) * dz
         if(prflag==1) write(*,"(IX,8(xf13.8))")i,stot2,frost,sliq2,
     +        s0,s1,smcx(i),sh2ox(i),dz
      enddo
      diftot = stot2 - stot1
      difliq = sliq2 - sliq1
      if(prflag==1)write(*,*)'diftot2 difliq2',diftot,difliq,DWT
      d1=diftot-dwt
      d2=difliq-dwt
      if(prflag==1)write(*,*)'# error',d1,d2,dwt
CBL Convert precision back
      DO I = NUPL, NLOWL
         SMC(I) = SMCX(I)
         SH2O(I) = SH2OX(I)
      ENDDO

      RETURN      
      END

C  SUBROUTINE RECALCULATES SOIL MOISTURE STATES INTO SAC-SMA STORAGES
C
cbl Based on the differences between the newly updated liquid water
cbl states (SH2O) and the previous sac storages (e.g. LZTWH), the SAC
cbl storages are updated.
CBL Fixed major bug (4/29), wherein SH2O was being reset to SWLT any time
CBL it fell below SWLT.  This fact, combined with LZFPH (aka F1WH) not 
CBL being included in the below re-distribution of freezing water, was 
CBL causing a water balance to occur when: SH2O(i) states were at SWLT
CBL while LZFPH was > 0, meaning FLAND could still reduce LZFPH, but 
CBL the corresponding amount of moisture could not be removed from
CBL SH2O, since they cannot be pulled below SWLT.
CBL
CBL Solution: include F1WH in calculation; compare SWH to SLIQ (not SX)
CBL *note: secondary problem, by resetting SH2O to SWLT, could be 
CBL affecting the stability of heat flux solution...unclear.
      SUBROUTINE FST2SAC2(ZSOIL,SMC,SH2O,NUP,NLW,SWLT,SMAX,RT,TWH,FWH,
     +     F1WH,TWC,FWC,F1WC,FROST,PX,ZONE3,TWM,FWM,F1WM,PRFLAG,CELLID)
      
      REAL ZSOIL(*),SMC(*),SH2O(*),s0,s1,denom,twm,fwm,f1wm
      INTEGER PRFLAG,CELLID
      LOGICAL ZONE3
      DOUBLE PRECISION SWCRATIO,LIQRATIO,SWC,STOT,SLIQ

      if (prflag==1) write(*,*)'fst2sac swlt smax',swlt,smax
      FROST=0.         
      STOT=0.
      SLIQ=0.
      SX=0.
CBL Correct the layering to allow for top soil layer
cbl      IF(NUP==1) THEN
cbl         SWH=TWH+FWH+F1WH+1000.*SWLT*(0.0-ZSOIL(NLW))/RT
cbl      ELSE
cbl         SWH=TWH+FWH+F1WH+1000.*SWLT*(ZSOIL(NUP-1)-ZSOIL(NLW))/RT
cbl      ENDIF
CBL Do not need to add SWLT, since we're comparing SWH to SLIQ now
      SWH=TWH+FWH+F1WH
      SWC=TWC+FWC+F1WC
CBL This was a problem anytime we computing top soil layer depth
      if(prflag==1)write(*,*)'I STOT FROST SLIQ s0 s1 smc sh2o'
      DO I=NUP,NLW
c8  Correction liquid water to be not less than swlt
         if(sh2o(i) .lt. swlt) then
            if(prflag==1)write(*,*)'sh2o lt swlt',i,sh2o(i),swlt
            sh2o(i) = swlt
         endif
         IF(NUP==1.AND.I==1) THEN
            DZ=0.0-ZSOIL(I)
         ELSE
            DZ=ZSOIL(I-1)-ZSOIL(I)
         ENDIF
         SX=SX+1000.*SH2O(I)*DZ/RT
         FROST=FROST+1000.*(SMC(I)-SH2O(I))*DZ/RT
         STOT=STOT+1000.*(SMC(I)-SWLT)*DZ/RT
         SLIQ=SLIQ+1000.*(SH2O(I)-SWLT)*DZ/RT
         s0 = 1000*(smc(i) - swlt) * dz/rt
         s1 = 1000 * (sh2o(i) - swlt) * dz/rt
         if(prflag==1) write(*,"(IX,7(xf13.8))")i,stot,frost,sliq,s0,s1,
     +        smc(i),sh2o(i)
         
      ENDDO
      xxx = abs(stot-TWC-FWC-F1WC) 
      xx = stot-twc-fwc-f1wc
      if(prflag==1)write(*,*)'x stot twc fwc f1wc',xxx,stot,twc,fwc,f1wc
      if(xxx .gt. 0.5.and.prflag==1) then
         write(*,*) '** fst2sac initial imbalance:',cellid,xxx
         write(*,'(6f12.6)') xxx,stot,twc,fwc,f1wc
         write(*,*)'cellid swlt',CELLID,swlt
         write(*,'(5f12.6)') (zsoil(i),i=nup-1,nlw)
         write(*,'(5f12.6)') (smc(i),i=nup,nlw)
         write(*,'(5f12.6)') (sh2o(i),i=nup,nlw)
      endif 

CBL Redistribute the total water contents so that they align and ensures
CBL the removal of small but cumulatively problematic rounding errors?
CBL split the small discrepancy by their relative saturations.
CBL Add most water the the emptiest
      IF (STOT.GT.0.) THEN
        SWCRATIO = SWC / STOT
        LIQRATIO = MIN(1.0 , MAX(0.0 , (SLIQ / STOT)))
      ELSE
         SWCRATIO = 1.0
         LIQRATIO = 0.0
      ENDIF
      if(prflag==1)write(*,*)'swcratio',SWCRATIO,SWC,STOT
      if(prflag==1)write(*,*)'liqratio',LIQRATIO,SLIQ,STOT
      if(prflag==1)write(*,*)'SWLT SMAX',swlt,smax
CBL Adjust SMC states first -- update stot,sliq
      DO I=NUP,NLW
         if(prflag==1)write(*,*)'smc sh20',i,smc(i),sh2o(i)
         SMC(I) = SWLT + (SMC(I)-SWLT)*SWCRATIO
         SH2O(I) = SWLT + MAX(0.0,(SH2O(I)-SWLT))*SWCRATIO
         if(prflag==1)write(*,*)'smc sh20',i,smc(i),sh2o(i)
         IF(SMC(I).GT.SMAX) THEN
            if(prflag==1)write(*,*)'smc.gt.smax',i,smc(i),sh2o(i)
            SH2O(I) = SH2O(I) + SMAX - SMC(I)
            SMC(I) = SMAX
            if(prflag==1)write(*,*)'smc.gt.smax',i,smc(i),sh2o(i)
         ENDIF
      ENDDO
C  DH2O IS LIQUID WATER CHANGE DUE TO FREEZING/THAWING, since sac states
cbl were not yet updated to reflect sh2o changes
      if(prflag==1)write(*,*)'twh fwh f1wh',twh,fwh,f1wh
      if(prflag==1)write(*,*)'twc fwc f1wc',twc,fwc,f1wc
      TWH=TWC*LIQRATIO
      FWH=FWC*LIQRATIO
      F1WH=F1WC*LIQRATIO
  
C  CHECK CONSISTENCY OF THE WATER BALANCE BETWEEN SAC-SMA AND
      if(prflag==1)write(*,*)'bottom of fst2sac2'
      FROST=0.         
      STOT=0.
      SLIQ=0.
      SX=0.
      SWH=TWH+FWH+F1WH
      SWC=TWC+FWC+F1WC
      if(prflag==1)write(*,*)'I STOT FROST SLIQ s0 s1 smc sh2o'
      DO I=NUP,NLW
         if(sh2o(i) .lt. swlt) then
            if(prflag==1)write(*,*)'sh2o lt swlt',i,sh2o(i),swlt
            sh2o(i) = swlt
         endif
         IF(NUP==1.AND.I==1) THEN
            DZ=0.0-ZSOIL(I)
         ELSE
            DZ=ZSOIL(I-1)-ZSOIL(I)
         ENDIF
         SX=SX+1000.*SH2O(I)*DZ/RT
         FROST=FROST+1000.*(SMC(I)-SH2O(I))*DZ/RT
         STOT=STOT+1000.*(SMC(I)-SWLT)*DZ/RT
         SLIQ=SLIQ+1000.*(SH2O(I)-SWLT)*DZ/RT
         s0 = 1000*(smc(i) - swlt) * dz/rt
         s1 = 1000 * (sh2o(i) - swlt) * dz/rt
         if(prflag==1) write(*,"(IX,7(xf13.8))")i,stot,frost,sliq,s0,s1,
     +        smc(i),sh2o(i)
         
      ENDDO
C  FROZEN GROUND STATES 
      S=TWC+FWC+F1WC
      if(prflag==1)write(*,*)'s=twc+fwc+f1wc',s,twc,fwc,f1wc
      if (prflag==1) then
         write(*,*)'s twc fwc f1wc',s,twc,fwc,f1wc
         write(*,*)'stot',stot
      endif
      IF(ABS(STOT-S) .GT. 0.5) then
         WRITE(*,*)
     + ' ** WARN ** fst2sac2 WATER IMBALANCE:',abs(stot-s),NUP,STOT,TWC,
     +        FWC,F1WC,cellid
      endif       
      S=TWH+FWH+F1WH
      if (prflag==1) then
         write(*,*)'s twh fwh f1wh',s,twh,fwh,f1w
         write(*,*)'sliq',sliq
      endif
      IF(ABS(SLIQ-S) .GT. 0.5) then
         WRITE(*,*)
     + ' ** WARN **fst2sac2 LIQUID IMBALANCE:',abs(sliq-s),NUP,SLIQ,TWH,
     +        FWH,F1WH,cellid
      endif 
      
      RETURN
      END
            
C      if(prflag==1)write(*,*)'twm fwm f1wm'
C      if(prflag==1)write(*,*)twm,fwm,f1wm
C      if (xx.gt.0.) then
C         denom = (twm-twc)/twm + (fwm-fwc)/fwm
C         if(f1wm.gt.0.) denom = denom + (f1wm-f1wc)/f1wm
C         if (prflag==1)write(*,*)'X gt 0 DENOM',denom,xx
C         if (prflag==1)write(*,*)'add x',twc,fwc,f1wc
C         twc = twc + (((twm-twc)/twm)/denom)*xx
C         fwc = fwc + (((fwm-fwc)/fwm)/denom)*xx
C         if(f1wm.gt.0.) f1wc = f1wc + (((f1wm-f1wc)/f1wm)/denom)*xx
C         if (prflag==1)write(*,*)'add x',twc,fwc,f1wc
C      endif
CCBL Remove most water from the fullest
C      if (xx.lt.0.) then
C         denom = twc + fwc + f1wc
C         if (prflag==1)write(*,*)'X lt 0 DENOM',denom,xx
C         if (prflag==1)write(*,*)'sub x',twc,fwc,f1wc
C         twc = twc + (twc/denom)*xx
C         fwc = fwc + (fwc/denom)*xx
C         if(f1wm.gt.0.) f1wc = f1wc + (f1wc/denom)*xx
C         if (prflag==1)write(*,*)'sub x',twc,fwc,f1wc
C      endif

cbl update the SAC storages based on the newly calculated mositure 
cbl states via Richards equation (conventional Noah)
cbl 
      SUBROUTINE H2O2SAC1(ZSOIL,SMC,SH2O,NUP,NLW,SWLT,SMAX,RT,TWH,FWH,
     +     F1WH,TWC,FWC,F1WC,TWM,FWM,F1WM,PRFLAG,CELLID)
      
      REAL ZSOIL(*),SMC(*),SH2O(*),s0,s1,denom,twm,fwm,f1wm
      INTEGER PRFLAG,CELLID
      DOUBLE PRECISION SWCRATIO,LIQRATIO,SWC,SWH,STOT,SLIQ

      if (prflag==1) write(*,*)'H2O2SAC swlt smax',swlt,smax
      FROST=0.         
      STOT=0.
      SLIQ=0.
      SX=0.
      SWH=TWH+FWH+F1WH
      SWC=TWC+FWC+F1WC
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
         write(*,*) '** h2o2sac initial imbalance:',cellid,xxx
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
        if(prflag==1)write(*,*)'swcratio',SWCRATIO,SWC,STOT
        if(prflag==1)write(*,*)'liqratio',LIQRATIO,SLIQ,STOT
      ELSE
         SWCRATIO = 0.0
         TWC=0.
         FWC=0.
         F1WC=0.
         SWC=0.    
         STOT=0.
         if(prflag==1)write(*,*)'swcratio',SWCRATIO,SWC,STOT
         if(prflag==1)write(*,*)'liqratio',LIQRATIO,SLIQ,STOT
      ENDIF
      IF (SLIQ.GT.0) THEN
         SWHRATIO = SHW / SLIQ
         if(prflag==1)write(*,*)'swhratio',SWHRATIO,SWH,SLIQ
      ELSE
         SH2ORATIO = 0.
         SLIQ = 0.
         TWH=0.
         FWH=0.
         F1WH=0.
         SWH=0.
         if(prflag==1)write(*,*)'swhratio',SWHRATIO,SWH,SLIQ
      ENDIF
      DWC = STOT - SWC
      if(prflag==1)write(*,*)'DWC',DWC,STOT,SWC
      IF(DWC.GT.0) THEN
CBL Distribute surplus based on deficits
         DENOM = (TWM-TWC) + (FWM-FWC) + (F1WM-F1WC)
         if(prflag==1)write(*,*)'dwc.gt.0 denom',denom,twc,fwc,f1wc
         TWC = TWC + ((TWM-TWC)/DENOM)*DWC
         FWC = FWC + ((FWM-FWC)/DENOM)*DWC
         F1WC = F1WC + ((F1WM-F1WC)/DENOM)*DWC
         if(prflag==1)write(*,*)'dwc.gt.0 denom',denom,twc,fwc,f1wc
      ELSE
CBL Remove moisture based on saturations
         DENOM = TWC + FWC + F1WC
         if(prflag==1)write(*,*)'dwc.le.0 denom',denom,twc,fwc,f1wc
         TWC = TWC + (TWC/DENOM)*DWC
         FWC = FWC + (FWC/DENOM)*DWC
         F1WC = F1WC + (F1WC/DENOM)*DWC
         if(prflag==1)write(*,*)'dwc.le.0 denom',denom,twc,fwc,f1wc
      ENDIF
      if(prflag==1)write(*,*)'twh fwh f1wh',twh,fwh,f1wh
      TWH=TWC*LIQRATIO
      FWH=FWC*LIQRATIO
      F1WH=F1WC*LIQRATIO
      if(prflag==1)write(*,*)'twh fwh f1wh',twh,fwh,f1wh



C  CHECK CONSISTENCY OF THE WATER BALANCE BETWEEN SAC-SMA AND
C  FROZEN GROUND STATES 
      S=TWC+FWC+F1WC
      if(prflag==1)write(*,*)'s=twc+fwc+f1wc',s,twc,fwc,f1wc
      if (prflag==1) then
         write(*,*)'s twc fwc f1wc',s,twc,fwc,f1wc
         write(*,*)'stot',stot
      endif
      IF(ABS(STOT-S) .GT. 0.5) then
         WRITE(*,*)
     + ' ** WARN ** h2o2sac TOTAL IMBALANCE:',abs(stot-s),NUP,STOT,TWC,
     +        FWC,F1WC,cellid
      endif       
      S=TWH+FWH+F1WH
      if (prflag==1) then
         write(*,*)'s twh fwh f1wh',s,twh,fwh,f1w
         write(*,*)'sliq',sliq
      endif
      IF(ABS(SLIQ-S) .GT. 0.5) then
         WRITE(*,*)
     + ' ** WARN ** h2o2sac LIQUID IMBALANCE:',abs(sliq-s),NUP,SLIQ,TWH,
     +        FWH,F1WH,cellid
      endif 
      
      RETURN
      END


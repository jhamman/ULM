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
      SUBROUTINE FST2SAC1(ZSOIL,SMC,SH2O,NUP,NLW,SWLT,RT,TWH,FWH,
     +     F1WH,TWC,FWC,F1WC,FROST,PX,ZONE3,TWM,FWM,F1WM,PRFLAG,CELLID)
      
      REAL ZSOIL(*),SMC(*),SH2O(*),s0,s1,denom,twm,fwm,f1wm
      INTEGER PRFLAG,CELLID
      LOGICAL ZONE3
      DOUBLE PRECISION SWCRATIO,LIQRATIO,SWC,STOT,SLIQ

      if (prflag==1) write(*,*)'fst2sac swlt',swlt
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
      if(xxx .gt. 0.5) then
         write(*,*) '** fst2sac NO TOTAL BALANCE:'
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
        LIQRATIO = MIN(0.0 , (SLIQ / STOT))
      ELSE
         SWCRATIO = 1.0
         LIQRATIO = 0.0
      ENDIF
      if(prflag==1)write(*,*)'swcratio',SWCRATIO,SWC,STOT
      if(prflag==1)write(*,*)'liqratio',LIQRATIO,SLIQ,STOT
CBL Adjust SMC states first
      DO I=NUP,NLW
         if(prflag==1)write(*,*)'smc sh20',i,smc(i),sh2o(i)
         SMC(I) = SMC(I) + (SMC(I)-SWLT)*SWCRATIO
         SH2O(I) = SH2O(I) + MAX(0.0,(SH2O(I)-SWLT))*SWCRATIO
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
CBL Take the difference between new SWH and SLIQ
cbl      DH2O=SX-SWH
      DH2O=SLIQ-SWH
      if(prflag==1) then
         write(*,*)'dh2o sliq swh swlt'
         write(*,*)dh2o,sliq,swh,swlt
      endif
      if(frost .eq. 0. .and. twh. eq. twc .and. fwh .eq. fwc
     +     .and. f1wh.eq.f1wc) then
         continue
         if(prflag==1)write(*,*)'first case, do nothing'
      else 
         IF(FROST .EQ. 0. .AND. DH2O .EQ. 0.) THEN
C  NO FROZEN GROUND, NO LIQUID WATER CHANGE
            TWH=TWC
            FWH=FWC
            F1WH=F1WC
            if(prflag==1)write(*,*)'second case, set H=C'
         ELSE 
C  CHANGE LIQUID WATER STORAGES

CBL Include F1WH in the redistribution -- make it last to be affected
cbl            SWH=TWH+FWH        only important to include for freeze
            SWH=TWH+FWH+F1WH
            if(prflag==1)write(*,*)'SWH',swh,twh,fwh,f1wh
            IF(DH2O .LT. 0.) THEN
C  UNFROZEN WATER REDUCTION (FREEZING)      
               if(prflag==1)write(*,*)'freezing dh2o lt 0',dh2o
CBL replace: twh was always .ge. 0, because it was being reset in fland1
cbl               IF(TWH .GE. 0.) THEN
               IF(TWH .GT. 0.) THEN
                  if(prflag==1)write(*,*)'twh.gt.0',twh
                  DSH=SWH+DH2O
                  if(prflag==1)write(*,*)'dsh',dsh,swh,dh2o
                  IF(DSH .GE. 0.) THEN
                     IF(SWH .GT. 1E-04) THEN
c  new -unlikely case for LZ, spring conditions?
                        if(px .gt. 0. .and. sliq .lt. stot) then
                           if(prflag==1)write(*,*)'px gt 0'
                           if(prflag==1)write(*,*)'fwh+dh2o',fwh,dh2o
                           fwh=fwh+dh2o
                           if(prflag==1)write(*,*)'fwh+dh2o',fwh,dh2o
                        else 
CBL split between TWH,FWH based on their relative contents
                           ALP=DH2O*TWH/SWH
                           if(prflag==1)write(*,*)'alp',alp,twh,fwh,swh
                           TWH=TWH+ALP
                           if(prflag==1)write(*,*)'alp',alp,twh,fwh,swh
                           FWH = FWH+DH2O-ALP
                           if(prflag==1)write(*,*)'alp',alp,twh,fwh,swh
                        endif
                     ELSE
                        if(prflag==1)write(*,*)'twh+=d2ho',twh
                        TWH=TWH+DH2O
                        if(prflag==1)write(*,*)'twh+=d2ho',twh
                     ENDIF
                  ELSE
CBL This case shouldn't occur; its when sliq redux exceeds sac-h storages
CBL i.e. |DH2O|>SWH which will require negative SAC-H states and lead to
CBL small water balance error when they are set back to zero in FLAND.
CBL Thus, reduce DSH by other FWH when enumerating TWH.
cbl                     TWH=DSH
                     if(prflag==1)write(*,*)'twh=dsh+fwh',twh,dsh,fwh
                     TWH=DSH+FWH
                     FWH=0.
                     if(prflag==1)write(*,*)'twh=dsh+fwh',twh,dsh,fwh
                  ENDIF
               ELSE
CBL Here TWH=0 so split the reduction between FWH and F1WH
                  if(prflag==1)write(*,*)'twh.gt.0',twh
                  if(prflag==1)write(*,*)'fwh+=dh2o',fwh,dh2o
                  FWH=FWH+DH2O
                  if(prflag==1)write(*,*)'fwh',fwh
cbl                  IF(FWH .LT. 0.) THEN
cbl                     TWH=TWH+FWH
cbl                  ENDIF
               ENDIF
CBL Check for depletions < 0. in one place only.
CBL For all cases first try to split the deficits between FWH,TWH which
CBL would be the upper zone case.  If both of these are exhausted, take
CBL water from F1WH, which is the lower zone case for water balance.
               if(fwh .lt. 0.) then
                  if(prflag==1)write(*,*)'fwh<0',fwh,twh
                  twh=twh+fwh
                  if(prflag==1)write(*,*)'fwh<0',fwh,twh
                  fwh=0.
                  if(prflag==1)write(*,*)'fwh<0',fwh,twh
               endif
               if (twh.lt.0.) then
                  if(prflag==1)write(*,*)'twh<0',twh,fwh
                  fwh = fwh + twh 
                  if(prflag==1)write(*,*)'twh<0',twh,fwh
                  twh=0.
                  if(prflag==1)write(*,*)'twh<0',twh,fwh
               endif
               if (fwh.lt.0.) then
                  if(prflag==1)write(*,*)'REDUX fwh<0',fwh,f1wh
                  if(f1wm.gt.0) then
                     f1wh = f1wh + fwh 
                     if(prflag==1)write(*,*)'REDUX fwh<0',fwh,f1wh
                     fwh=0.
                     if(prflag==1)write(*,*)'REDUX fwh<0',fwh,f1wh
                  else
                     if(prflag==1)write(*,*)'warn UZ depleted',fwh,twh
                  endif
               endif
               if (f1wh.lt.0.)write(*,*)'warn F1WH redux:',cellid,f1wh
            ELSE
C  UNFROZEN WATER INCREASE (THAWING) --use old SWH definition,favors TWH      
               SWH = TWH + FWH
               if(prflag==1)write(*,*)'thawing dh2o ge 0',dh2o
               IF(TWH .LT. 0.) THEN
                  if(prflag==1)write(*,*)'twh lt 0',twh
                  TWH=TWH+DH2O
                  if(prflag==1)write(*,*)'twh+=dh2o',twh
               ELSE
                  IF(SWH .GT. 1E-04) THEN
                     if(prflag==1)write(*,*)'swh gt 1e-4',swh
                     ALP=DH2O*TWH/SWH
                     if(prflag==1)write(*,*)'alp',alp,twh,swh
                     TWH=TWH+ALP
                     if(prflag==1)write(*,*)'twh+=alp',twh
                     FWH=FWH+DH2O-ALP
                     if(prflag==1)write(*,*)'fwh',fwh,dh2o,alp
                  ELSE
                     TWH=TWH+DH2O
                     if(prflag==1)write(*,*)'twh+=dh2o',twh
                  ENDIF
               ENDIF
CBL Check for overflows in one place only.
               IF(TWH .GT. TWC) THEN
                  if(prflag==1)write(*,*)'twh gt twc',twh,twc,fwh,twm
                  FWH=FWH+TWH-TWC
                  TWH=TWC
                  if(prflag==1)write(*,*)'fwh+=twh-twc',twh,twc,fwh
               ENDIF
               IF(FWH .GT. FWC) THEN                       
                  if(prflag==1)write(*,*)'fwh gt fwc',fwh,fwc,twc,fwm
                  TWH=TWH+FWH-FWC
                  FWH=FWC
                  if(prflag==1)write(*,*)'twh+=fwh-fwc',fwh,fwc,twc
               ENDIF
               IF(TWH .GT. TWC) THEN
                  if(prflag==1)write(*,*)'twh gt twc',twh,twc,fwh,twm
                  FWH=FWH+TWH-TWC
                  TWH=TWC
                  if(prflag==1)write(*,*)'fwh+=twh-twc',twh,twc,fwh
               ENDIF
               IF (FWH.GT.FWC) THEN
                  if (F1WH.LT.F1WC) then
                     if(prflag==1)write(*,*)'f1wh+=fwh-fwc',fwh,fwc,f1wh
                     F1WH=F1WH+FWH-FWC
                     FWH=FWC
                     if(prflag==1)write(*,*)'f1wh+=fwh-fwc',fwh,fwc,f1wh
                  else
                     rem = fwh-fwc
                     write(*,*)'warn FWH over:',cellid,rem,fwh,fwc
                     write(*,*)'f1wh f1wc',f1wh,f1wc
                  endif
               ENDIF
            ENDIF
         ENDIF        
      endif
      
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
     + ' ** WARN ** NO TOTAL WATER BALANCE:',abs(stot-s),NUP,STOT,TWC,
     +        FWC,F1WC
      endif       
      S=TWH+FWH+F1WH
      if (prflag==1) then
         write(*,*)'s twh fwh f1wh',s,twh,fwh,f1w
         write(*,*)'sliq',sliq
      endif
      IF(ABS(SLIQ-S) .GT. 0.5) then
         WRITE(*,*)
     + ' ** WARN ** NO LIQUID WATER BALANCE:',abs(sliq-s),NUP,SLIQ,TWH,
     +        FWH,F1WH
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

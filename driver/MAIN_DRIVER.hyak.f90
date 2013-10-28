PROGRAM noah

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! Modifications:
  ! 2007-Nov-06 Fixed reporting of SnowT and added SnowSoilT.		Ben Livneh
  ! 2007-Nov-06 Fixed reporting of SoilMoist.				TJB
  ! 2007-Nov-12 Moved computation of SOLNET to inside SFLX.		TJB
  ! 2007-Nov-12 Added LSTSNW and day_of_year to SFLX call, for new
  !             snow albedo computation scheme.				Ben Livneh and TJB
  ! 2007-Nov-15 Fixed reporting of output record times.			TJB
  ! 2007-Dec-05 Improved looping over forcing, output files.		TJB
  ! 2008-May-05 Removed T12, LADJCH, and FFROZP, and removed day_of_year
  !             from SFLX call, for compatibility with NOAH 2.8.	Ben Livneh
  ! 2008-Jul-24 Added PackWater to track liquid water within snow       Ben Livneh
  ! 2008-Jul-24 Added SliqFrac and SnowTProf				Ben Livneh
  ! 2008-Aug-12 Added error_flag to alert bounding errors in SFLX
  !             SNWPAC iterative temperature solver			Ben Livneh
  ! 2008-Oct-07 Included SAC parameters for use as unified model         Ben Livneh
  ! 2010-Jan-06 Did a bunch of consistency checks between driver and sflx, units, wb components

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  CHARACTER*60 wbfile_debug

  ! Define local variables
  INTEGER :: I,J,K,L,p
  INTEGER :: MODEL_STEP
  INTEGER :: FORCING_STEP
  INTEGER :: OUTPUT_STEP
  INTEGER :: MODEL_STEP_COUNT
  INTEGER :: tsteps_per_day,FORCING_STEP_OF_DAY
  INTEGER :: days_per_year,day_of_year
  INTEGER :: year,month,day,hour,minute,sec
  INTEGER :: year_new,month_new,day_new,hour_new,minute_new,sec_new,day_of_year_new
  INTEGER :: timestp
  REAL    :: balance
  REAL    :: LAI_tmp
  REAL    :: DQSDT
  REAL    :: wb_error, eb_error
  INTEGER :: IARGC
  REAL    :: DT_TAIR_band,DT_PRCP_band,DT_SPFH_SAT_band,DT_SPFH_band
  REAL    :: DQSDT2_band,TH2_band,ALB_TOT_band,DT_SOLNET_band
  REAL    :: FFROZP,RAIN_ratio,SNOW_ratio
  REAL    :: ETA_band,H_band,EC1_band,EDIR1_band,ETT1_band,ESNOW_band,DRIP_band,DEW_band
  REAL    :: BETA_band,ETP_band,S_band,FLX1_band,FLX2_band,FLX3_band
  REAL    :: SNOMLT_band,RUNOFF1_band,RUNOFF2_band,RUNOFF3_band,EVAP_band
  REAL    :: RC_band,PC_band,RSMIN_band,RCS_band,RCT_band,RCQ_band,RCSOIL_band
  REAL    :: SOILW_band,SOILM_band,Fdepth_band,Tdepth_band
  REAL    :: RAIN_band,SNOW_band,SOILMOIST_total_band,SOILMOIST_save_band,SWE_band,SWE_save_band
  REAL    :: FX,RGL,HS,PackWater_band,CMC_total_band,CMC_save_band
  REAL, ALLOCATABLE:: ET_band(:), LWnet_model_step(:), Snowf_model_step(:), Rainf_model_step(:)
  LOGICAL :: forc_exist,lvrain
  INTEGER :: num_fsteps, tsteplen_save, FSTEP_COUNT, error_flag, prflag, refcellid
  REAL    :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13,snowcover
  REAL    :: wb_sum,eb_sum,balance_counter
  REAL    MAXSMC(30)
  REAL    WLTSMC(20)

  DATA MAXSMC/0.37308, 0.38568, 0.41592, 0.46758, 0.47766, 0.43482, 0.41592, 0.4764, 0.44868, 0.42348, 0.48144, 0.46128, 0.464, 0.000, 0.200, 0.421, 0.457, 0.200, 0.395, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
  DATA WLTSMC/0.03469064, 0.05199094, 0.08743051, 0.14637683, 0.10712489, 0.13941739, 0.15698002, 0.24386303, 0.21203782, 0.20755672, 0.28488226, 0.28290603, 0.069, 0.000, 0.012, 0.028, 0.135, 0.012, 0.023, 0.000/
  ! Initialize error status
  status = 0
  error_flag = 0

  ! READ CONTROL FILE
  IF(IARGC().NE.1)THEN
     PRINT*, 'USAGE:  noah <control_file>'
     STOP
  ENDIF
  CALL GETARG(1,CNTRFL)
!  CNTRFL = './controlfile'      
  WRITE(*,*) 'Reading control file ',TRIM(CNTRFL)
  CALL READCNTL
  ! Convert date/time information
  MODEL_DT_REAL = REAL(MODEL_DT)
  OUTPUT_DT_REAL = REAL(OUTPUT_DT)
  CALL GREG2JUL(YEAR0,MONTH0,DAY0,JULDAY0)
  CALL BUILD_DATE_STRING(YEAR0,MONTH0,DAY0,0,0,0,MODEL_START_TIME)

!!!!!! print the model parameters here...
  ! READ MODEL GEOMETRY
  IF (trim(PARAM_TYPE) == 'netcdf') THEN
    WRITE(*,*)'Reading model geometry ',TRIM(LSC)
  ELSE
    WRITE(*,*)'Reading model geometry ',TRIM(MASK_FILE)
  END IF
  CALL GET_GRID
  ! Allocate arrays based on the model geometry
  CALL ALLOCATE_ARRAYS
  ALLOCATE(ET_band(MAXNSOIL))
  ALLOCATE(LWnet_model_step(landlen))
  ALLOCATE(Snowf_model_step(landlen))
  ALLOCATE(Rainf_model_step(landlen))

  write(*,*)'Total active grid-cells',landlen
  ! READ LAND SURFACE CHARACTERISTICS
  IF (trim(PARAM_TYPE) == 'netcdf') THEN
    WRITE(*,*)'Reading model parameters ',TRIM(LSC)
  ELSE
    WRITE(*,*)'Reading model parameters, starting with ',TRIM(NSOIL_FILE)
  END IF
  CALL READ_LSC
  IF (nbands > 1) THEN
    WRITE(*,*)'Reading snowband file ',TRIM(SNOWBAND_FILE)
  ELSE
    WRITE(*,*)'Snowbands are turned OFF'
  END IF
  CALL READ_SNOWBANDS
!!!!!! print the main geometry/lsc parameters here...

  ! READ INITIAL CONDITIONS
  WRITE(*,*) 'Reading initial conditions'
  CALL READ_INITIAL
!!!!!! print the main initial parameters here...

  ! CHECK INITIAL CONDITIONS AND MODEL GEOMETRY
  WRITE(*,*) 'Checking initial conditions and model geometry'
  CALL CHECK_INITIAL

  ! BUILD CURRENT FILE SUFFIX
  WRITE(SUFFIX,'(".",i4.4,i2.2,".nc")')YEAR0,MONTH0

  ! OPEN ATMOSPHERIC FORCING DATA
  FORCFILE = TRIM(FORCING)//TRIM(SUFFIX)
  WRITE(*,*) 'Opening atmospheric data ',TRIM(FORCFILE)
  CALL OPEN_FORCING
  tsteplen_save = tsteplen
  ! OPEN OUTPUT FILES
  OUTFILES(1) = TRIM(RESULT)//"/eb"//TRIM(SUFFIX)
  OUTFILES(2) = TRIM(RESULT)//"/wb"//TRIM(SUFFIX)
  OUTFILES(3) = TRIM(RESULT)//"/sur"//TRIM(SUFFIX)
  OUTFILES(4) = TRIM(RESULT)//"/sub"//TRIM(SUFFIX)
  OUTFILES(5) = TRIM(RESULT)//"/eva"//TRIM(SUFFIX)
  OUTFILES(6) = TRIM(RESULT)//"/csp"//TRIM(SUFFIX)
  OUTFILES(7) = TRIM(RESULT)//"/pe"//TRIM(SUFFIX)
  WRITE(*,*) 'Opening output files:'
  CALL OPEN_OUTPUT
! ************ HACK for SW MONITOR **************
!  WRITE(*,*) 'EB file: ',TRIM(OUTFILES(1))
  WRITE(*,*) 'WB file: ',TRIM(OUTFILES(2))
!  WRITE(*,*) 'SUR file: ',TRIM(OUTFILES(3))
!  WRITE(*,*) 'SUB file: ',TRIM(OUTFILES(4))
!  WRITE(*,*) 'EVA file: ',TRIM(OUTFILES(5))
!  WRITE(*,*) 'CSP file: ',TRIM(OUTFILES(6))
  WRITE(*,*) 'PE file: ',TRIM(OUTFILES(7))
! ************ END HACK for SW MONITOR **************

  ! OPEN RESTART FILE
!  RESTFILE = TRIM(RESTART)//TRIM(SUFFIX)
!  WRITE(*,*) 'Opening restart file ',TRIM(RESTFILE)
!  CALL OPEN_RESTART
  ! Initialize time counters
  CALL COMPUTE_NUM_FSTEPS(FORCING_DT,YEAR0,MONTH0,DAY0,YEAR_FINAL,MONTH_FINAL,DAY_FINAL,num_fsteps)
  MADTT = FORCING_DT/MODEL_DT
  IF (MADTT < 1) THEN
    WRITE(*,*)'ERROR: MODEL_DT',MODEL_DT,'IS GREATER THAN FORCING_DT',FORCING_DT
    STOP
  END IF
  OUTPUT_STEP_RATIO = OUTPUT_DT/MODEL_DT
  IF (OUTPUT_STEP_RATIO < 1) THEN
    WRITE(*,*)'ERROR: MODEL_DT',MODEL_DT,'IS GREATER THAN OUTPUT_DT',OUTPUT_DT
    STOP
  END IF
  tsteps_per_day = 24*3600/FORCING_DT
  OUTPUT_STEP = 1
  MODEL_STEP_COUNT = 1
  year = YEAR0
  month = MONTH0
  day = DAY0
  day_of_year = JULDAY0
  hour = 0
  minute = 0
  sec = 0
  FORCING_STEP = 1
  FSTEP_COUNT = 1
  write(*,*)'tsteps_per_day ',tsteps_per_day
  write(*,*)'FORCING_DT ',FORCING_DT
  write(*,*)'MODEL_DT ',MODEL_DT
  write(*,*)'MADTT ',MADTT
  write(*,*)'OUTPUT_DT ',OUTPUT_DT
  write(*,*)'OUTPUT_STEP_RATIO ',OUTPUT_STEP_RATIO
  write(*,*)'num_fsteps ',num_fsteps
  ! Initialize "old" values
  DO I = 1,landlen
    SoilMoistTotal(I) = 0.0
    DO J = 1, NBANDS
      DO K = 1, NSOIL(I)
        SoilMoistTotal(I)=SoilMoistTotal(I)+SMC(I,J,K)*SOILDEPTH(I,K)*band_area(I,J)*1000.0
      END DO
    END DO
  END DO
  SoilMoistTotal_old = SoilMoistTotal
  DO I = 1, landlen
    SWE(I) = 0
    PackWater(I) = 0
    CanopInt(I) = 0
    DO J = 1, nbands
      SWE(I) = SWE(I) + SNEQV(I,J) * band_area(I,J) * 1000.0
      PackWater(I) = PackWater(I) + PACH20(I,J) * band_area(I,J) * 1000.0
      CanopInt(I) = CanopInt(I) + CMC(I,J) * band_area(I,J) * 1000.0
    END DO
  END DO
  SWE_old       = SWE
  CanopInt_old  = CanopInt

  ! READ ATMOSPHERIC FORCING DATA - BEGINNING OF 1ST FORCING STEP
  IF (FORCING_STEP .LE. tsteplen) THEN
    CALL READ_FORCING(FORCING_STEP,SWdownnow,LWdownnow,Tairnow,Qairnow,Rainfnow,Snowfnow,PSurfnow,Windnow)
  ELSE
    WRITE(*,*)'ERROR: FORCING_STEP',FORCING_STEP,'not found in forcing file',FORCFILE
    STOP
  ENDIF

  ! RUN TIME AND SPATIAL LOOP
  DO

!    write(*,*)
!    write(*,*)'FORCING_STEP ',FORCING_STEP
!    write(*,*)'FSTEP_COUNT ',FSTEP_COUNT

    ! READ ATMOSPHERIC FORCING DATA - BEGINNING OF NEXT FORCING STEP
    IF (FSTEP_COUNT .LT. num_fsteps) THEN
      ! At least one more forcing step remains in this simulation
      IF (FORCING_STEP .LT. tsteplen) THEN
        ! Next forcing step is contained in current forcing file
        CALL READ_FORCING(FORCING_STEP+1,SWdownplus,LWdownplus,Tairplus,Qairplus,Rainfplus,Snowfplus,PSurfplus,Windplus)
      ELSE
        ! We've reached the end of the current forcing file(s)
        tsteplen_save = tsteplen
        CALL CLOSE_FORCING
        ! Build new forcing file name(s)
        year_new = year
        month_new = month
        day_new = day
        hour_new = hour
        minute_new = minute
        sec_new = sec
        CALL INCREMENT_DATE(FORCING_DT,year_new,month_new,day_new,hour_new,minute_new,sec_new)
        WRITE(SUFFIX_NEW,'(".",i4.4,i2.2,".nc")') year_new,month_new
        FORCFILE = TRIM(FORCING)//TRIM(SUFFIX_NEW)
        ! Check whether next forcing files exist
        INQUIRE(FILE=FORCFILE,EXIST=forc_exist)
        IF (forc_exist) THEN
          ! Open new forcing file(s)
          WRITE(*,*) 'Opening atmospheric data ',TRIM(FORCFILE)
          CALL OPEN_FORCING
          CALL READ_FORCING(1,SWdownplus,LWdownplus,Tairplus,Qairplus,Rainfplus,Snowfplus,PSurfplus,Windplus)
        ELSE
          ! Forcings have run out before end of simulation
          WRITE(*,*) 'ERROR: Cannot find next atmospheric data file',TRIM(FORCFILE)
          STOP
        ENDIF
      ENDIF
    ELSE
      ! We are on the final forcing step of the simulation; set *plus vars equal to *now vars
      SWdownplus = SWdownnow
      LWdownplus = LWdownnow
      Tairplus   = Tairnow
      Qairplus   = Qairnow
      Rainfplus  = Rainfnow
      Snowfplus  = Snowfnow
      PSurfplus  = PSurfnow
      WINDplus   = WINDnow
    ENDIF

    FORCING_STEP_OF_DAY = mod((FORCING_STEP-1),tsteps_per_day) + 1

    ! INTERPOLATE MONTHLY ALBEDO, ETC TO VALUE FOR CURRENT DAY OF YEAR
    CALL INTERP_MONTHLY(day_of_year)


    IF (MODEL_STEP_COUNT == 1) THEN

      ! Initialize output variables that must be aggregated from model time step to output time step
      SWnet    = 0.0
      LWnet    = 0.0
      Qle      = 0.0
      Qh       = 0.0
      Qg       = 0.0
      Qf       = 0.0
      Qv       = 0.0
      Qa       = 0.0
      Snowf    = 0.0
      Rainf    = 0.0
      Evap     = 0.0
      Qs       = 0.0
      Qsb      = 0.0
      Qsm      = 0.0
      Qst      = 0.0
      PotEvap  = 0.0
      ECanop   = 0.0
      TVeg     = 0.0
      ESoil    = 0.0
      SubSnow  = 0.0

    END IF

    ! PotEvapPE is a special output variable - we aggregate it to the forcing step
    ! and do not aggregate across elevation bands
    PotEvapPE  = 0.0

    ! For each forcing step, there are MADTT model timesteps
    DO MODEL_STEP = 1,MADTT

!      write(*,*)'MODEL_STEP ', MODEL_STEP
!      write(*,*)'OUTPUT_STEP ', OUTPUT_STEP
!      write(*,*)'MODEL_STEP_COUNT ', MODEL_STEP_COUNT

      ! Reset timestep-specific grid-cell-total variables
      ETA = 0.0
      H = 0.0
      EVAP_total = 0.0
      EC1 = 0.0
      EDIR1 = 0.0
      ET = 0.0
      ETT1 = 0.0
      ESNOW = 0.0
      DRIP = 0.0
      DEW = 0.0
      BETA = 0.0
      ETP = 0.0
      S = 0.0
      FLX1 = 0.0
      FLX2 = 0.0
      FLX3 = 0.0
      SNOMLT = 0.0
      RUNOFF1 = 0.0
      RUNOFF2 = 0.0
      RUNOFF3 = 0.0
      Albedo_ALMA = 0.0
      ACond = 0.0
      PC = 0.0
      ACondMax = 0.0
      RCS = 0.0
      RCT = 0.0
      RCQ = 0.0
      RCSOIL = 0.0
      SOILW = 0.0
      SOILM = 0.0
      SWE = 0.0
      PackWater = 0.0
      CanopInt = 0.0
      SWEVeg = 0.0
      SoilMoist = 0.0
      SoilTemp = 0.0
      SMLiqFrac = 0.0
      SMFrozFrac = 0.0
      SoilWet = 0.0
      SoilMoistTotal = 0.0
      RootMoist = 0.0
      Fdepth = 0.0
      Tdepth = 0.0
      SnowFrac = 0.0
      SnowDepth = 0.0
      VegT       = 0.0
      BaresoilT  = 0.0
      AvgSurfT   = 0.0
      RadT       = 0.0
      SAlbedo    = 0.0
      band_area_with_snow    = 0.0
      LWnet_model_step    = 0.0
      Snowf_model_step    = 0.0
      Rainf_model_step    = 0.0
      SnowT = 0.0
 !Ben Livneh track effective snowpack temp TPACK = SnowTProf
      SnowTProf = 0.0
      wb_sum = 0.0
!      write(*,*)'day of year',day_of_year
!$OMP PARALLEL DO
      DO I = 1,landlen
!  Start debugging conditional
!         IF (I==1 .or. MOD(I,32) == 0) THEN
!         IF (I==1 .or. MOD(I,int(landlen*1/5)) == 0 .or. MOD(I,int(landlen*2/5)) == 0 .or. MOD(I,int(landlen*3/5)) == 0 .or. MOD(I,int(landlen*4/5)) == 0) THEN
! Major conus basins
         IF (I==1 .or. MOD(I,int(landlen*1/4)) == 0 .or. MOD(I,int(landlen*2/4)) == 0 .or. MOD(I,int(landlen*3/4)) == 0) THEN
! Minor conus basins
!         IF (I==1 .or. MOD(I,int(landlen*1/2)) == 0) THEN
! Apriori case or major basin 'all' case
!        IF (I >= 0) THEN
            refcellid = i
         ! INTERPOLATE FORCING DATA TO MODEL TIME STEP
         ! This produces DT_* variables and CZMODEL, given *now and *plus variables and time/location
         CALL INTERP_FORCING(day_of_year,FORCING_STEP_OF_DAY,MODEL_STEP,I)

!         if(i==1) write(*,*)'####msc',model_step_count,Rainf(i)
         if (I == -1) THEN
            write(*,*)'Before snowband loop'
            write(*,*)X(I),Y(I),LAT(X(I),Y(I)),LON(X(I),Y(I))
            write(*,*)
            write(*,*)NSOIL(I),SOILDEPTH(I,:),SOILTYP(I),SLOPETYP(I),TBOT(I)
            write(*,*)VEGTYP(I),SHDFAC_D(I),SHDMIN(I),ALBEDO_D(I),SNOALB(I),PTU(I)
            write(*,*)
            write(*,*)'soldn',DT_SOLDN(I),DT_LWDN(I),DT_TAIR(I),DT_PRCP(I),DT_PSFC(I),DT_WIND(I)
            write(*,*)
         END IF

         ! SAC
         UZTWM = UZTWM_2d(I)
         UZFWM = UZFWM_2d(I)
         UZK   = UZK_2d(I)
         PCTIM = PCTIM_2d(I)
         ADIMP = ADIMP_2d(I)
         RIVA  = RIVA_2d(I)
         ZPERC = ZPERC_2d(I)
         REXP  = REXP_2d(I)
         LZTWM = LZTWM_2d(I)
         LZFSM = LZFSM_2d(I)
         LZFPM = LZFPM_2d(I)
         LZSK  = LZSK_2d(I)
         LZPK  = LZPK_2d(I)
         PFREE = PFREE_2d(I)
         SIDE  = SIDE_2d(I)
         RSERV = RSERV_2d(I)
         ! Sensitivity parameters
         WCRIT = WCRIT_2d(I)
         PSNOW1 = PSNOW1_2d(I)
         PSNOW2 = PSNOW2_2d(I)
         RICHARDS = RICHARDS_2d(I)

         ! Loop over snow bands
         DO J = 1, nbands

            IF (band_area(I,J) > 0) THEN

               ! Adjust precip and temperature for the current snow band
               ! Note: we aren't taking into account changes in the other forcings with elevation
               ! e.g. air pressure, incoming long/short-wave radiation, wind speed
               if(i==9999)then
                  p=1
               else
                  p=0
               endif
               if(p==1) write(*,*)'###p',MODEL_STEP_COUNT,DT_PRCP(I),(dt_prcp(i)*model_dt)
               DT_PRCP_band = Pfactor(I,J)*DT_PRCP(I)
               DT_TAIR_band = Tfactor(I,J)+DT_TAIR(I)
               !            write(*,*)'prcp_init',DT_PRCP(I),Pfactor(I,J)
               ! Determine precip tyep
               FFROZP = 0.0
               IF (DT_PRCP_band .GT. 0.0) THEN
                  IF (DT_TAIR_band .LE. TFREEZ) THEN
                     FFROZP = 1.0
                  ENDIF
               ENDIF
               
               ! Aggregate Rainf and Snowf
               ! Note: these conditions must be exactly the same as those in SFLX to be correct
               ! Preserve ratios so that if total precip is adjusted in SFLX (dew,frost) then
               ! the relevant amounts can be updated after SFLX call.
               IF (FFROZP > 0.5 .OR. T1(I,J) < TFREEZ) THEN
                  SNOW_band = DT_PRCP_band
                  RAIN_band = 0
               ELSE
                  SNOW_band = 0
                  RAIN_band = DT_PRCP_band
               END IF
               if (DT_PRCP_band .gt. 0.0) then
                  RAIN_ratio = RAIN_band / (RAIN_band + SNOW_band)
                  SNOW_ratio = SNOW_band / (RAIN_band + SNOW_band)
               else
                  RAIN_ratio = 0.0
                  SNOW_ratio = 0.0
                  DT_PRCP_band = 0.
               endif
               ! CALCULATE SPECIFIC HUMIDITIES
               ! note - this doesn't take into accound pressure differences among elevation bands
               DT_SPFH_band = DT_SPFH(I)
               CALL QDATAP(DT_TAIR_band,DT_PSFC(I),DT_SPFH_SAT_band)
               IF (DT_SPFH_band .LE. 1.e-5) THEN
                  DT_SPFH_band = 1.e-5
               END IF
               IF (DT_SPFH_band .GE. DT_SPFH_SAT_band*0.99) THEN
                  DT_SPFH_band = DT_SPFH_SAT_band*0.99
               END IF
               
               ! CALCULATE SLOPE OF SAT SPECIFIC HUMIDITY CURVE FOR PENMAN: DQSDT2
               ! note - this doesn't take into accound pressure differences among elevation bands
               DQSDT2_band = DQSDT(DT_TAIR_band, DT_PSFC(I))
               
               ! CALCULATE OTHER VARS
               TH2_band = DT_TAIR_band + 0.0098 * Z(I)
               
               ! Set up miscellaneous water balance terms
               SOILMOIST_total_band = 0
               DO K = 1, NSOIL(I)
                  SOILMOIST_total_band = SOILMOIST_total_band + SMC(I,J,K)*SOILDEPTH(I,K)*1000.0
               END DO
               SOILMOIST_save_band = SOILMOIST_total_band
               SWE_band = SNEQV(I,J)*1000.0
               PackWater_band = PACH20(I,J)*1000
               SWE_save_band = SWE_band
               CMC_total_band = 0
               CMC_total_band = CMC_total_band + CMC(I,J) * 1000
               CMC_save_band = CMC_total_band
               
               if (I == 9999) then
                  write(*,*)'above sflx'
                  write(*,*)'cell',I,'band',J,'DT_PRCP_band',DT_PRCP_band,'DT_TAIR_band',DT_TAIR_band,'SNEQV',SNEQV(I,J)
                  write(*,*)'RAIN_band',RAIN_band,'SNOW_band',SNOW_band
                  write(*,*)'DT_SOLNET',DT_SOLNET_band,'DT_SPFH',DT_SPFH_band,'TH2',TH2_band,'DT_SPFH_SAT',DT_SPFH_SAT_band,'DQSDT2',DQSDT2_band
                  write(*,*)'State vars'
                  write(*,*)CMC(I,J),T1(I,J),STC(I,J,:),SMC(I,J,:),SH2O(I,J,:),SNOWH(I,J),SNEQV(I,J),ALB_TOT_band,CH(I,J),CM(I,J)
                  write(*,*)'STC',STC(I,J,:),'SMC',SMC(I,J,:),'SH20',SH2O(I,J,:)
                  write(*,*)
                  write(*,*)'MODEL TYPE',MODEL_TYPE
                  write(*,*)'SOILDEPTH',SOILDEPTH(I,:)
               end if
               
               ! ASSIGN CURRENT ELEMENTS OF 2D STATE VAR ARRAYS TO SINGLE VALUE VARS
               !            UZTWC = UZTWC_2d(I,J)
               !            UZFWC = UZFWC_2d(I,J)
               !            LZTWC = LZTWC_2d(I,J)
               !            LZFSC = LZFSC_2d(I,J)
               !            LZFPC = LZFPC_2d(I,J)
               !            ADIMC = ADIMC_2d(I,J)
               ! Replaced by SACST array
               
               !            SOILMOIST_total_band = UZTWC+UZFWC+LZTWC+LZFSC+LZFPC
               !            SOILMOIST_save_band = SOILMOIST_total_band
               
!               if (i==1) then
!                  prflag = 1
!                  write(*,*)'cell',i
!                                 write(*,*) '0.SMC:',(SMC(I,J,K),K=1,4)
!                                 write(*,*) '0.SH2O:',(SH2O(I,J,K),K=1,4)
!                                 write(*,*) '0.SACST:',(SACST(I,J,K),K=1,5)
!                                 write(*,*) '0.FRZST:',(FRZST(I,J,K+5),K=1,5)
!               else
                  prflag = 0
!               endif
               
               !           write(*,*)'t1before',T1(I,J),TPACK(I,J),DT_TAIR_band
               
               !########### Be sure to convert all moisture fluxes to M for SFLX
               
               ! CALL LAND-SURFACE PHYSICS
               
               ! Appended call for TPACK(I,J) - Ben Livneh Nov 2007
               !            write(*,*)'above sflx','NSNOW(I,J)',NSNOW(I,J),'MODEL_TYPE',MODEL_TYPE
               !            write(*,*)'prcp',DT_PRCP_band,'stc',STC(I,J,:)
               if (prflag == 1) then
               write(*,*)'cell id',I,'model step',MODEL_STEP_COUNT,'month',month,'day',day
               endif
               if(p==1) write(*,*)'##p',MODEL_STEP_COUNT,DT_PRCP_BAND,(dt_prcp_band*model_dt)
               CALL SFLX ( &
                    ICE(I),MODEL_DT_REAL,Z(I),NSOIL(I),SOILDEPTH(I,:), &
                    DT_LWDN(I),DT_SOLDN(I),DT_SOLNET_band, DT_PSFC(I),DT_PRCP_band, &
                    DT_TAIR_band,DT_SPFH_band,DT_WIND(I), &
                    TH2_band,DT_SPFH_SAT_band,DQSDT2_band, &
                    VEGTYP(I),SOILTYP(I),SLOPETYP(I),SHDFAC_D(I), &
                    SHDMAX(I),SHDMIN(I),PTU(I), &
                    ALBEDO_D(I),SNOALB(I),TBOT(I), &
                    CMC(I,J),T1(I,J),TPACK(I,J),PACH20(I,J),STC(I,J,:), &
                    SMC(I,J,:),SH2O(I,J,:),SNOWH(I,J),SNEQV(I,J),ALB_TOT_band,CH(I,J),CM(I,J), &
                    ! SAC PARAMETERS AND STATES
                    FRZST(I,J,:),FRZPAR(I,:),SACST(I,J,:),SACPAR(I,:), &
                    ETA_band,H_band, &
                    ! State and output variables for PEMB snow model
                    NSNOW(I,J),DSNOW(I,J,:),PSNOW(I,J,:),RTTSNOW(I,J,:),RTTDTSNOW(I,J,:),&
                    WTSNOW(I,J,:),WTDTSNOW(I,J,:),TTSNOW(I,J,:),TTDTSNOW(I,J,:),&
                    ! ----------------------------------------------------------------------
                    ! OUTPUTS, DIAGNOSTICS, PARAMETERS BELOW GENERALLY NOT NECESSARY WHEN
                    ! COUPLED WITH E.G. A NWP MODEL (SUCH AS THE NOAA/NWS/NCEP MESOSCALE ETA
                    ! MODEL).  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES.
                    ! ----------------------------------------------------------------------
                    EC1_band,EDIR1_band,ET_band(:),ETT1_band,ESNOW_band,DRIP_band,DEW_band, &
                    BETA_band,ETP_band,S_band, &
                    FLX1_band,FLX2_band,FLX3_band, &
                    SNOMLT_band,SNCOVR(I,J), &
                    RUNOFF1_band,RUNOFF2_band,RUNOFF3_band, &
                    RC_band,PC_band,RSMIN_band,LAI_tmp,RCS_band,RCT_band,RCQ_band,RCSOIL_band, &
                    FX,SOILW_band,SOILM_band,RGL,HS,error_flag,prflag,I,MODEL_STEP_COUNT, &
                    W_WILT(I),SMCDRY(I),SMCREF(I),W_SAT(I),NROOT(I),CZMODEL,LSTSNW(I,J),MODEL_TYPE,lvrain,WCRIT,PSNOW1,PSNOW2,RICHARDS)
               

               if (I == 9999) then
                  write(*,*)'below sflx'
                  write(*,*)'cell',I,'band',J,'DT_PRCP_band',DT_PRCP_band,'DT_TAIR_band',DT_TAIR_band,'SNEQV',SNEQV(I,J)
                  write(*,*)'RAIN_band',RAIN_band,'SNOW_band',SNOW_band
                  write(*,*)'DT_SOLNET',DT_SOLNET_band,'DT_SPFH',DT_SPFH_band,'TH2',TH2_band,'DT_SPFH_SAT',DT_SPFH_SAT_band,'DQSDT2',DQSDT2_band
                  write(*,*)'State vars'
                  write(*,*)CMC(I,J),T1(I,J),STC(I,J,:),SMC(I,J,:),SH2O(I,J,:),SNOWH(I,J),SNEQV(I,J),ALB_TOT_band,CH(I,J),CM(I,J)
                  write(*,*)'STC',STC(I,J,:),'SMC',SMC(I,J,:),'SH20',SH2O(I,J,:)
                  write(*,*)
                  write(*,*)'MODEL TYPE',MODEL_TYPE
                  write(*,*)'SOILDEPTH',SOILDEPTH(I,:)
               end if
               if(p==1) write(*,*)'#p',MODEL_STEP_COUNT,DT_PRCP_BAND,(dt_prcp_band*model_dt)
               ! If the snow temperature has exceeded its bounds in the iterative
               ! solver print an error message
               if (error_flag == 1) then
                  write(*,*)'snow temperature in SNWPAC beneath : cell',I,'band',J
                  write(*,*)'Model step',MODEL_STEP,'Output step',OUTPUT_STEP
                  error_flag = 0
               end if
               
               !            temp1 = DT_SOLDN(I)*(1-ALB_TOT_band)
               !            temp2 = DT_LWDN(I)
               !            temp3 = ETA_band
               !            temp4 = H_band
               !            temp5 = S_band
               !            temp6 = FLX3_band
               !            temp7 = FLX1_band + FLX2_band
               !            temp8 = DelSurfHeat
               !            temp9 = DelColdCont
               
               eb_error = DT_SOLDN(I)*(1-ALB_TOT_band) + DT_LWDN(I) - ETA_band - H_band - S_band - FLX3_band - FLX2_band - FLX1_band
               ! Error appears to be meaningless unless we consider heat changes of soil and snow over time
!               write(*,*)'#cell eb_error',I,eb_error,SWE_band,SNCOVR(I,J)
!               write(*,*)'SW',DT_SOLDN(I)*(1-ALB_TOT_band),'LW',DT_LWDN(I)
!               write(*,*)'ETA',ETA_band,'H_band',H_band
!               write(*,*)'S_band',S_band,'FLX3',FLX3_band
!               write(*,*)'FLX2',FLX2_band,'FLX1',FLX1_band
               
               ! Update miscellaneous water balance terms
               ! Convert units where appropriate
               SOILMOIST_total_band = 0
               DO K = 1, NSOIL(I)
                  SOILMOIST_total_band = SOILMOIST_total_band + SMC(I,J,K)*SOILDEPTH(I,K)*1000.0
               END DO
               CMC_total_band = 0
               CMC_total_band = CMC_total_band + CMC(I,J) * 1000
               SWE_band = SNEQV(I,J)*1000.0
               PackWater_band = PACH20(I,J)*1000
               EC1_band = EC1_band/LVH2O
               EDIR1_band = EDIR1_band/LVH2O
               ETT1_band = ETT1_band/LVH2O
               ! Be sure correct lsub is applied; rain on snow case
               if (lvrain .eq. .true.) then
                  ESNOW_band = ESNOW_band/LVH2O
               else
                  ESNOW_band = ESNOW_band/LSUBS
               endif
               EVAP_band = EC1_band + EDIR1_band + ETT1_band + ESNOW_band
               RUNOFF1_band = RUNOFF1_band*1000.0
               RUNOFF2_band = RUNOFF2_band*1000.0
               RUNOFF3_band = RUNOFF3_band*1000.0
               ETP_band = ETP_band/LVH2O
               SNOMLT_band = SNOMLT_band*1000.0/MODEL_DT_REAL
               if (prflag==1) write(*,*)'Evap_band ec edir ett esnow',EVAP_band,EC1_band,EDIR1_band,ETT1_band,ESNOW_band
               ! Water balance error check
               !            wb_error = (RAIN_band + SNOW_band - EVAP_band - RUNOFF1_band - RUNOFF2_band)*MODEL_DT_REAL + SWE_band - SWE_save_band + SOILMOIST_total_band - SOILMOIST_save_band + CMC_total_band - CMC_save_band
               !            temp1 = (RAIN_band + SNOW_band)*MODEL_DT_REAL
               temp1 = (DT_PRCP_band)*MODEL_DT_REAL
               temp2 = (EVAP_band)*MODEL_DT_REAL
               temp3 = (RUNOFF1_band + RUNOFF2_band)*MODEL_DT_REAL
               temp4 = SWE_band
               temp5 = SWE_save_band
               temp6 = SOILMOIST_total_band
               temp7 = SOILMOIST_save_band
               temp8 = SWE_band - SWE_save_band
               temp9 = SOILMOIST_total_band - SOILMOIST_save_band
               temp10 = CMC_total_band
               temp11 = CMC_save_band
               temp13 = CMC_total_band - CMC_save_band
               
               wb_error = temp1 - temp2 - temp3 - temp8*SNCOVR(I,J) - temp9 - temp13
               if (wb_error*wb_error > 0.1) then
                  write(*,*)'###cell wb_error',I,wb_error,SWE_band,SNCOVR(I,J)
                  write(*,*)'prcp',temp1,'evap',temp2,'runoff',temp3
                  write(*,*)'delswe',temp8*SNCOVR(I,J),'delsm',temp9,'delcmc',temp13
               endif
               
               if (abs(wb_error)>WB_ERROR_TOL*MODEL_DT_REAL/(24*3600)) then
                  !            if (I == 2) then
                  write(*,*)'below sflx'
                  write(*,*)
                  write(*,*)'cell',I,'band',J,'wb error',wb_error,'mm, over model time step of',MODEL_DT_REAL,'sec'
                  write(*,*)'FORCING STEP',FORCING_STEP,'MODEL STEP',MODEL_STEP,'MODEL STEP COUNT',MODEL_STEP_COUNT,'OUTPUT STEP',OUTPUT_STEP
                  write(*,*)
                  write(*,*)'WB terms (mm) for this band and model step:'
                  write(*,*)'RAIN',RAIN_band*MODEL_DT_REAL,'SNOW',SNOW_band*MODEL_DT_REAL,'EVAP',EVAP_band*MODEL_DT_REAL
                  write(*,*)'RUNOFF1',RUNOFF1_band*MODEL_DT_REAL,'RUNOFF2',RUNOFF2_band*MODEL_DT_REAL,'RUNOFF3',RUNOFF3_band*MODEL_DT_REAL
                  write(*,*)'SWE_band',SWE_band,'SWE_save_band',SWE_save_band,'SOILMOIST_total_band',SOILMOIST_total_band,'SOILMOIST_old',SOILMOIST_save_band
                  write(*,*)
                  write(*,*)'Soil Moisture'
                  DO K = 1, NSOIL(I)
!                     WRITE(*,*)'SMC(I,J,K)*SOILDEPTH(I,K)*1000.0',SMC(I,J,K),SOILDEPTH(I,K),(SMC(I,J,K)*SOILDEPTH(I,K)*1000)
                  END DO
               end if
               
               ! Calculate freeze/thaw depths
               Fdepth_band = 0.0
               DO K = 1, NSOIL(I)
                  IF (Fdepth_band > 0 .and. STC(I,J,K) > TFREEZ) exit
                  IF (STC(I,J,K) < TFREEZ) THEN
                     Fdepth_band = SOILDEPTH_ACCUM(I,K)
                  END IF
               END DO
               Tdepth_band = 0.0
               DO K = 1, NSOIL(I)
                  IF (STC(I,J,K) < TFREEZ) exit
                  Tdepth_band = SOILDEPTH_ACCUM(I,K)
               END DO
               
               ! Aggregate output variables across all snow bands
               ! Fluxes
               ETA(I) = ETA(I) + ETA_band*band_area(I,J)
               H(I) = H(I) + H_band*band_area(I,J)
               EVAP_total(I) = EVAP_total(I) + EVAP_band*band_area(I,J)
               EC1(I) = EC1(I) + EC1_band*band_area(I,J)
               EDIR1(I) = EDIR1(I) + EDIR1_band*band_area(I,J)
               ET(I,:) = ET(I,:) + ET_band(:)*band_area(I,J)
               ETT1(I) = ETT1(I) + ETT1_band*band_area(I,J)
               ESNOW(I) = ESNOW(I) + ESNOW_band*band_area(I,J)
               DRIP(I) = DRIP(I) + DRIP_band*band_area(I,J)
               DEW(I) = DEW(I) + DEW_band*band_area(I,J)
               BETA(I) = BETA(I) + BETA_band*band_area(I,J)
               ETP(I) = ETP(I) + ETP_band*band_area(I,J)
               S(I) = S(I) + S_band*band_area(I,J)
               FLX1(I) = FLX1(I) + FLX1_band*band_area(I,J)
               FLX2(I) = FLX2(I) + FLX2_band*band_area(I,J)
               FLX3(I) = FLX3(I) + FLX3_band*band_area(I,J)
               SNOMLT(I) = SNOMLT(I) + SNOMLT_band*band_area(I,J)
               RUNOFF1(I) = RUNOFF1(I) + RUNOFF1_band*band_area(I,J)
               RUNOFF2(I) = RUNOFF2(I) + RUNOFF2_band*band_area(I,J)
               RUNOFF3(I) = RUNOFF3(I) + RUNOFF3_band*band_area(I,J)
               ! State variables
               if (RC_band > 0 .AND. ACond(I) /= NODATA) then
                  ACond(I) = ACond(I) + (1/RC_band)*band_area(I,J)
               else
                  ACond(I) = NODATA
               end if
               PC(I) = PC(I) + PC_band*band_area(I,J)
               if (RSMIN_band > 0 .AND. ACondMax(I) /= NODATA) then
                  ACondMax(I) = ACondMax(I) + (1/RSMIN_band)*band_area(I,J)
               else
                  ACondMax(I) = NODATA
               end if
               RCS(I) = RCS(I) + RCS_band*band_area(I,J)
               RCT(I) = RCT(I) + RCT_band*band_area(I,J)
               RCQ(I) = RCQ(I) + RCQ_band*band_area(I,J)
               RCSOIL(I) = RCSOIL(I) + RCSOIL_band*band_area(I,J)
               SOILW(I) = SOILW(I) + SOILW_band*band_area(I,J)
               SOILM(I) = SOILM(I) + SOILM_band*band_area(I,J)
               Albedo_ALMA(I) = Albedo_ALMA(I) + ALB_TOT_band*band_area(I,J)
               SWE(I) = SWE(I) + SNEQV(I,J)*1000.0*band_area(I,J)
               PackWater(I) = PackWater(I) + PACH20(I,J) * 1000 * band_area(I,J)
               CanopInt(I) = CanopInt(I) + CMC(I,J)*1000.0*band_area(I,J)
               if (T1(I,J) < TFREEZ) then
                  SWEVeg(I) = SWEVeg(I) + CMC(I,J)*1000.0*band_area(I,J)
               end if
               DO K = 1, NSOIL(I)
                  SoilMoist(I,K) = SoilMoist(I,K) + SMC(I,J,K)*SOILDEPTH(I,K)*1000.0*band_area(I,J)
                  SoilTemp(I,K)  = SoilTemp(I,K) + STC(I,J,K)*band_area(I,J)
                  SMLiqFrac(I,K)  = SMLiqFrac(I,K) + SH2O(I,J,K)/SMC(I,J,K)*band_area(I,J)
                  SMFrozFrac(I,K)  = SMFrozFrac(I,K) + (SMC(I,J,K)-SH2O(I,J,K))/SMC(I,J,K)*band_area(I,J)
                  SoilWet(I) = SoilWet(I) + ( (SMC(I,J,K)-W_WILT(I)) / (W_SAT(I)-W_WILT(I)) )*(SOILDEPTH(I,K)/SOILDEPTH_ACCUM(I,MAXNSOIL))*band_area(I,J)
               END DO
               SoilMoistTotal(I)=SoilMoistTotal(I)+SOILMOIST_total_band*band_area(I,J)
               RootMoist(I) = RootMoist(I) + (SOILW_band - W_WILT(I))*SOILDEPTH_ACCUM(I,MAXNSOIL)*1000.0*band_area(I,J)
               SnowDepth(I) = SnowDepth(I) + SNOWH(I,J)*band_area(I,J)
               SnowFrac(I) = SnowFrac(I) + SNCOVR(I,J)*band_area(I,J)
               VegT(I) = VegT(I) + T1(I,J)*band_area(I,J)
               BaresoilT(I) = BaresoilT(I) + T1(I,J)*band_area(I,J)
               AvgSurfT(I) = AvgSurfT(I) + T1(I,J)*band_area(I,J)
               Fdepth(I) = Fdepth(I) + Fdepth_band*band_area(I,J)
               Tdepth(I) = Tdepth(I) + Tdepth_band*band_area(I,J)
               ! Misc
               LWnet_model_step(I) = LWnet_model_step(I)+(DT_LWDN(I)-SIGMA*T1(I,J)**4.0)*band_area(I,J)
               Snowf_model_step(I) = Snowf_model_step(I)+DT_PRCP_band*SNOW_ratio*band_area(I,J)
               Rainf_model_step(I) = Rainf_model_step(I)+DT_PRCP_band*RAIN_ratio*band_area(I,J)
               !if(p==1)write(*,*)'rat',snow_ratio,rain_ratio
               if(p==1)write(*,*)'##sn',snowf_model_step(i),(DT_PRCP_band*SNOW_ratio*band_area(I,J))
               if(p==1)write(*,*)'##rn',rainf_model_step(i),(DT_PRCP_band*RAIN_ratio*band_area(I,J))
               RadT(I) = RadT(I) + T1(I,J)*band_area(I,J)
               if (SNEQV(I,J) > 0) then
                  SAlbedo(I) = SAlbedo(I) + ALB_TOT_band*band_area(I,J)
                  SnowT(I) = SnowT(I) + T1(I,J)*band_area(I,J)
                  !Ben Livneh track effective snowpack temp TPACK = SnowTProf
                  SnowTProf(I) = SnowTProf(I) + TPACK(I,J)*band_area(I,J)
                  band_area_with_snow(I) = band_area_with_snow(I) + band_area(I,J)
               end if
               
               ! AGGREGATE POTEVAP TO THE FORCING TIME STEP FOR OUTPUT TO PE FILE
               PotEvapPE(I,J) = PotEvapPE(I,J) + ETP_band / MADTT
               
            END IF
            wb_sum = wb_sum + wb_error
         END DO
         
         ! Normalize variables that have NODATA over part of the grid cell
         ! (if appropriate)
         SAlbedo(I) = SAlbedo(I)/band_area_with_snow(I)
         SatSoil(I) = MAXSMC(SOILTYP(I))*SOILDEPTH_ACCUM(I,MAXNSOIL)*1000
         WltSoil(I) = WLTSMC(SOILTYP(I))*SOILDEPTH_ACCUM(I,MAXNSOIL)*1000
         SnowT(I) = SnowT(I)/band_area_with_snow(I)
         !Ben Livneh track effective snowpack temp TPACK = SnowTProf and liquid water PACH20 fraction SliqFrac
         SnowTProf(I) = SnowTProf(I)/band_area_with_snow(I)
         SliqFrac(I) = PackWater(I)/SWE(I)
         
         ! AGGREGATE OUTPUT FLUXES TO OUTPUT TIME STEP
         SWnet(I)= SWnet(I)+(1.0-Albedo_ALMA(I))*DT_SOLDN(I) / OUTPUT_STEP_RATIO
         LWnet(I)= LWnet(I)+LWnet_model_step(I) / OUTPUT_STEP_RATIO
         Qle(I)  = Qle(I) + ETA(I) / OUTPUT_STEP_RATIO
         Qh(I)   = Qh(I) + H(I) / OUTPUT_STEP_RATIO
         Qg(I)   = Qg(I) + S(I) / OUTPUT_STEP_RATIO
         Qf(I)   = Qf(I) + FLX3(I) / OUTPUT_STEP_RATIO
         Qv(I)   = Qv(I) + ESNOW(I) / OUTPUT_STEP_RATIO
         Qa(I)   = Qa(I) + (FLX1(I)+FLX2(I)) / OUTPUT_STEP_RATIO
         
         Snowf(I) = Snowf(I) + Snowf_model_step(I) / OUTPUT_STEP_RATIO
         Rainf(I) = Rainf(I) + Rainf_model_step(I) / OUTPUT_STEP_RATIO
         if(p==1)write(*,*)'#sn',snowf(i),(Snowf_model_step(I) / OUTPUT_STEP_RATIO)
         if(p==1)write(*,*)'#rn',rainf(i),(Rainf_model_step(I) / OUTPUT_STEP_RATIO)
         Evap(I) = Evap(I) + EVAP_total(I) / OUTPUT_STEP_RATIO
         Qs(I)   = Qs(I) + RUNOFF1(I) / OUTPUT_STEP_RATIO
         Qsb(I)  = Qsb(I) + RUNOFF2(I) / OUTPUT_STEP_RATIO
         Qsm(I)  = Qsm(I) + SNOMLT(I) / OUTPUT_STEP_RATIO
         
         write(*,*) "model_step",model_step, OUTPUT_STEP_RATIO
         write(*,*) "Qs i",Qs(I),I
         
         PotEvap(I) = PotEvap(I) + ETP(I) / OUTPUT_STEP_RATIO
         ECanop(I)  = ECanop(I) + EC1(I) / OUTPUT_STEP_RATIO
         TVeg(I)    = TVeg(I) + ETT1(I) / OUTPUT_STEP_RATIO
         ESoil(I)   = ESoil(I) + EDIR1(I) / OUTPUT_STEP_RATIO
         SubSnow(I) = SubSnow(I) + ESNOW(I) / OUTPUT_STEP_RATIO
         
         !IF (I == 2) THEN
         !      write(*,*)'end of landlen'
         !      write(*,*)'Rainf',Rainf(I),'Snowf',Snowf(I),'Evap',Evap(I),'Qs',Qs(I),'Qsb',Qsb(I),'CMC',CMC(I),'SMC',SMC(I,:),'SNEQV',SNEQV(I)
         !      write(*,*)
         !      write(*,*)'swnet',SWnet(I),'lwnet',LWnet(I),'qle',Qle(I),'qh',Qh(I),'qg',Qg(I),'qf',Qf(I),'qv',Qv(I),'qa',Qa(I),'delcc',DelColdCont(I)
         !      write(*,*)
         !END IF
         ! End debugging conditional
!      ENDIF
! set cell values to the reference cell values
      else
               Albedo_ALMA(I) = Albedo_ALMA(refcellid)
               SWE(I) = SWE(refcellid)
               PackWater(I) = PackWater(refcellid)
               CanopInt(I) = CanopInt(refcellid)
               SWEVeg(I) = SWEVeg(refcellid)
               DO K = 1, NSOIL(I)
                  SoilMoist(I,K) = SoilMoist(refcellid,K)
                  SoilTemp(I,K)  = SoilTemp(refcellid,K)
                  SMLiqFrac(I,K)  = SMLiqFrac(refcellid,K)
                  SMFrozFrac(I,K)  = SMFrozFrac(refcellid,K)
                  SoilWet(I) = SoilWet(refcellid)
               END DO
               SoilMoistTotal(I)=SoilMoistTotal(refcellid)
               RootMoist(I) = RootMoist(refcellid)
               SnowDepth(I) = SnowDepth(refcellid)
               SnowFrac(I) = SnowFrac(refcellid)
               VegT(I) = VegT(refcellid)
               BaresoilT(I) = BaresoilT(refcellid)
               AvgSurfT(I) = AvgSurfT(refcellid)
               Fdepth(I) = Fdepth(refcellid)
               Tdepth(I) = Tdepth(refcellid)
               ! Misc
               LWnet_model_step(I) = LWnet_model_step(refcellid)
               Snowf_model_step(I) = Snowf_model_step(refcellid)
               Rainf_model_step(I) = Rainf_model_step(refcellid)
               RadT(I) = RadT(refcellid)
               SAlbedo(I) = SAlbedo(refcellid)
               SatSoil(I) = SatSoil(refcellid)
               WltSoil(I) = WltSoil(refcellid)
               SnowT(I) = SnowT(refcellid)
               SnowTProf(I) = SnowTProf(refcellid)
               band_area_with_snow(I) = band_area_with_snow(refcellid)
               ! AGGREGATE POTEVAP TO THE FORCING TIME STEP FOR OUTPUT TO PE FILE
!               PotEvapPE(I,J) = PotEvapPE(,J)
               ! Normalize variables that have NODATA over part of the grid cell
               ! (if appropriate)
               SAlbedo(I) = SAlbedo(refcellid)
               SnowT(I) = SnowT(refcellid)
               SnowTProf(I) = SnowTProf(refcellid)
               SliqFrac(I) = SliqFrac(refcellid)
         
               ! AGGREGATE OUTPUT FLUXES TO OUTPUT TIME STEP
               SWnet(I)= SWnet(refcellid)
               LWnet(I)= LWnet(refcellid)
               Qle(I)  = Qle(refcellid)
               Qh(I)   = Qh(refcellid)
               Qg(I)   = Qg(refcellid)
               Qf(I)   = Qf(refcellid)
               Qv(I)   = Qv(refcellid)
               Qa(I)   = Qa(refcellid)
               
               Snowf(I) = Snowf(refcellid)
               Rainf(I) = Rainf(refcellid)
               Evap(I) = Evap(refcellid)
               Qs(I)   = Qs(refcellid)
               Qsb(I)  = Qsb(refcellid)
               Qsm(I)  = Qsm(refcellid)
               
               PotEvap(I) = PotEvap(refcellid)
               ECanop(I)  = ECanop(refcellid)
               TVeg(I)    = TVeg(refcellid)
               ESoil(I)   = ESoil(refcellid)
               SubSnow(I) = SubSnow(refcellid)

! end of the selective cell conditional         
        end if

     END DO
      !$OMP END PARALLEL DO
      
      IF (MODEL_STEP_COUNT == OUTPUT_STEP_RATIO) THEN

        ! FINISH CALCULATING OUTPUT VARIABLES

        ! Water balance variables
        Qfz          = 0.0
        Qst          = 0.0

        ! Surface state variables
        WHERE (SWE == 0.0)
          SnowT     = NODATA
 !Ben Livneh track effective snowpack temp TPACK = SnowTprof
          SnowTProf = NODATA
        END WHERE
        SurfStor    = 0.0

        ! Evap-related variables
        EWater   = 0.0
        EvapSnow = 0.0
        SubSurf  = 0.0

        ! Changes in storage terms
        DelSurfHeat  = 0.0
        DelColdCont  = 0.0
        DelSoilMoist = SoilMoistTotal - SoilMoistTotal_old
        DelSWE       = SWE - SWE_old
        DelSurfStor  = 0.0
        DelIntercept = CanopInt - CanopInt_old


!write(*,*)
!write(*,*)'Cell 1'
!write(*,*)'ALMA output terms (in ALMA units):'
!write(*,*)'Canopint',Canopint(1),'Canopint_old',Canopint_old(1),'SWE',SWE(1),'SWE_old',SWE_old(1),'SoilMoistTotal',SoilMoistTotal(1),'SoilMoistTotal_old',SoilMoistTotal_old(1)
!write(*,*)
!write(*,*)'Evap',Evap(1),'ECanop',ECanop(1),'TVeg',TVeg(1),'ESoil',ESoil(1),'SubSnow',SubSnow(1)
!write(*,*)
!write(*,*)'Rainf',Rainf(1),'Snowf',Snowf(1),'Evap',Evap(1),'Qs',Qs(1),'Qsb',Qsb(1),'DelIntercept',DelIntercept(1),'DelSurfStor',DelSurfStor(1),'DelSWE',DelSWE(1),'DelSoilMoist',DelSoilMoist(1)
!write(*,*)
!write(*,*)'SWdown',SWdownnow(1),'LWdown',LWdownnow(1),'swnet',SWnet(1),'lwnet',LWnet(1),'qle',Qle(1),'qh',Qh(1),'qg',Qg(1),'qf',Qf(1),'qv',Qv(1),'qa',Qa(1),'DelSurfHeat',DelSurfHeat(1),'DelColdCont',DelColdCont(1)
!write(*,*)

        ! Energy and Water Balance Check
        wb_error = wb_sum / landlen
!        wb_error = 0
!        eb_error = 0
!        DO I = 1, landlen
!          wb_error = wb_error + (Rainf(I) + Snowf(I) - Evap(I) - Qs(I) - Qsb(I)) * OUTPUT_DT - (DelIntercept(I) + DelSurfStor(I) + DelSWE(I) + DelSoilMoist(I))
!          eb_error = eb_error + (SWnet(I) + LWnet(I) - Qh(I) - Qle(I) - Qg(I) - Qf(I) - Qa(I)) - (DelSurfHeat(I) + DelColdCont(I))/OUTPUT_DT
!        END DO
!        wb_error = wb_error/landlen
!        eb_error = eb_error/landlen
        WRITE(*,*)
        WRITE(*,*)'OUTPUT STEP',OUTPUT_STEP,'FORCING STEP',FORCING_STEP
        WRITE(*,*)'avg wb_error per cell:',wb_error,'mm'
!        WRITE(*,*)'avg eb_error per cell:',eb_error,'W/m2'
!        wb_error = (Rainf(1) + Snowf(1) - Evap(1) - Qs(1) - Qsb(1)) * OUTPUT_DT - (DelIntercept(1) + DelSurfStor(1) + DelSWE(1) + DelSoilMoist(1))
!        eb_error = (SWnet(1) + LWnet(1) - Qh(1) - Qle(1) - Qg(1) - Qf(1) - Qa(1)) - (DelSurfHeat(1) + DelColdCont(1))/OUTPUT_DT
!        WRITE(*,*)'wb_error, cell 1:',wb_error,'mm'
!        WRITE(*,*)'eb_error, cell 1:',eb_error,'W/m2'
        WRITE(*,*)

        ! Update "old" values
        SoilMoistTotal_old= SoilMoistTotal
        SWE_old      = SWE
        CanopInt_old = CanopInt

        ! DRIVER STEP 9 
        
        ! Debug statement
        ! Write results for this interval to output files
        CALL WRITE_OUTPUT(OUTPUT_STEP)
        MODEL_STEP_COUNT = 1
        OUTPUT_STEP = OUTPUT_STEP + 1

      ELSE

        MODEL_STEP_COUNT = MODEL_STEP_COUNT + 1

      END IF

    END DO

    ! Write PotEvapPE to PE file
    CALL WRITE_PE(FORCING_STEP)

    ! Close output files if appropriate
    IF (FSTEP_COUNT == num_fsteps .OR. FORCING_STEP == tsteplen_save) THEN
      ! Close the output file(s)
      CALL CLOSE_OUTPUT
      ! Write restart file
!      WRITE(*,*) 'Writing restart file ',RESTFILE
!      CALL WRITE_RESTART
    ENDIF

    ! Check for end of simulation
    IF (FSTEP_COUNT == num_fsteps) THEN
      WRITE(*,*)'Reached end of simulation'
      EXIT
    ENDIF

    ! Save this forcing step's forcing values
    SWdownminus = SWdownnow
    LWdownminus = LWdownnow
    Tairminus   = Tairnow
    Qairminus   = Qairnow
    PSurfminus  = PSurfnow
    WINDminus   = WINDnow

    SWdownnow = SWdownplus
    LWdownnow = LWdownplus
    Tairnow   = Tairplus
    Qairnow   = Qairplus
    Rainfnow  = Rainfplus
    Snowfnow  = Snowfplus
    PSurfnow  = PSurfplus
    WINDnow   = WINDplus

    ! Increment forcing counters
    FORCING_STEP = FORCING_STEP + 1
    FSTEP_COUNT = FSTEP_COUNT + 1

    ! Calculate next forcing step's time counter values
    CALL INCREMENT_DATE(FORCING_DT,year,month,day,hour,minute,sec)
    CALL GREG2JUL(year,month,day,day_of_year)

    ! Open new output files if appropriate
    IF (FORCING_STEP .GT. tsteplen_save) THEN
      WRITE(SUFFIX,'(".",i4.4,i2.2,".nc")')year,month
      ! Check whether next forcing files exist
      FORCFILE = TRIM(FORCING)//TRIM(SUFFIX)
      INQUIRE(FILE=FORCFILE,EXIST=forc_exist)
      IF (forc_exist) THEN
        ! Open new restart file
!        RESTFILE = TRIM(RESTART)//TRIM(SUFFIX)
!        WRITE(*,*) 'Opening restart file ',TRIM(RESTFILE)
!        CALL OPEN_RESTART
        ! Open new output file(s)
        OUTFILES(1) = TRIM(RESULT)//"/eb"//TRIM(SUFFIX)
        OUTFILES(2) = TRIM(RESULT)//"/wb"//TRIM(SUFFIX)
        OUTFILES(3) = TRIM(RESULT)//"/sur"//TRIM(SUFFIX)
        OUTFILES(4) = TRIM(RESULT)//"/sub"//TRIM(SUFFIX)
        OUTFILES(5) = TRIM(RESULT)//"/eva"//TRIM(SUFFIX)
        OUTFILES(6) = TRIM(RESULT)//"/csp"//TRIM(SUFFIX)
        OUTFILES(7) = TRIM(RESULT)//"/pe"//TRIM(SUFFIX)
        WRITE(*,*) 'Opening output files:'
        CALL BUILD_DATE_STRING(YEAR,MONTH,DAY,0,0,0,MODEL_START_TIME)
        CALL OPEN_OUTPUT
! ************ HACK for SW MONITOR **************
!        WRITE(*,*) 'EB file: ',TRIM(OUTFILES(1))
        WRITE(*,*) 'WB file: ',TRIM(OUTFILES(2))
!        WRITE(*,*) 'SUR file: ',TRIM(OUTFILES(3))
!        WRITE(*,*) 'SUB file: ',TRIM(OUTFILES(4))
!        WRITE(*,*) 'EVA file: ',TRIM(OUTFILES(5))
!        WRITE(*,*) 'CSP file: ',TRIM(OUTFILES(6))
        WRITE(*,*) 'PE file: ',TRIM(OUTFILES(7))
! ************ END HACK for SW MONITOR **************
        OUTPUT_STEP = 1
        tsteplen_save = tsteplen
        FORCING_STEP = 1
      ELSE
        IF (FSTEP_COUNT .LT. num_fsteps) THEN
          WRITE(*,*)'ERROR: Reached end of forcing files before end of simulation'
        ENDIF
        STOP
      ENDIF

    ENDIF

  END DO

END PROGRAM noah

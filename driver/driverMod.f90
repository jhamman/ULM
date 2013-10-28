MODULE driverMod

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! DECLARATIONS OF GLOBAL VARIABLES & DEFINITIONS OF GLOBAL CONSTANTS

  ! Modifications:
  ! 2007-Nov-06 Added T12.							Ben Livneh
  ! 2007-Nov-12 Added LSTSNW.							TJB
  ! 2007-Dec-05 Added SUFFIX_NEW.						TJB
  ! 2008-Jul-15 Replaced T12 and SnowSoilT with TPACK and SnowTProf		Ben Livneh
  ! 2008-Jul-24 Added PACH20 and SliqFrac					Ben Livneh
  ! 2008-Oct-07 Included SAC parameters for use as unified model                Ben Livneh
  ! 2009-Feb-15 Added Anderson PEMB snow model parameters for unified model     Ben Livneh

  ! RCS ID string: $Id: driverMod.f90,v 1.8 2008/07/30 22:33:03 vicadmin Exp $

  IMPLICIT NONE

  SAVE

  ! NetCDF subroutines
!  include '/usr/local/i386/include/netcdf.inc'
!  include '$INC_NETCDF/netcdf.inc'
!  include '/usr/include/netcdf-3/netcdf.inc'
  include '/usr/include/netcdf.inc'


  ! Model constants
  REAL, PARAMETER :: TFREEZ=273.15
  REAL, PARAMETER :: LVH2O=2.501000E+6
  REAL, PARAMETER :: LSUBS=2.83E+6
  REAL, PARAMETER :: SIGMA=5.672E-8
  REAL, PARAMETER :: NODATA=1.E+20
  INTEGER, PARAMETER :: NODATA_INT=-9999
  INTEGER, PARAMETER :: MAXNSNOW=100 ! Max snow layers for Anderson PEMB snow model
  LOGICAL, PARAMETER :: LOCAL = .false.
  LOGICAL, PARAMETER :: LLANDUSE = .false.
  LOGICAL, PARAMETER :: LSOIL = .false.
  REAL, PARAMETER :: LAPSE_RATE = 6.5 ! deg C/km; for snow elevation bands
!  REAL, PARAMETER :: WB_ERROR_TOL = 10.0 ! water balance error tolerance for warning messages, in mm/day over the model time step
  REAL, PARAMETER :: WB_ERROR_TOL = 100.0 ! water balance error tolerance for warning messages, in mm/day over the model time step
  REAL, PARAMETER :: EB_ERROR_TOL = 1.0 ! energy balance error tolerance for warning messages, in W/m2 over the model time step
!  REAL, PARAMETER :: Rd=287.04 ! Gas constant for dry air
!  REAL, PARAMETER :: G=9.81 ! Gravitational acceleration at earth's surface

  ! Error-reporting
  INTEGER :: status

  ! Input/output filenames
  CHARACTER(len=200) :: CNTRFL
  CHARACTER(len=200) :: LSC
  CHARACTER(len=200) :: MASK_FILE
  CHARACTER(len=200) :: NSOIL_FILE
  CHARACTER(len=200) :: SOILDEPTH_FILE
  CHARACTER(len=200) :: SOILTYP_FILE
  CHARACTER(len=200) :: SLOPETYP_FILE
  CHARACTER(len=200) :: TBOT_FILE
  CHARACTER(len=200) :: VEGTYP_FILE
  CHARACTER(len=200) :: SHDFAC_FILE
  CHARACTER(len=200) :: ALBEDO_FILE
  CHARACTER(len=200) :: SNOALB_FILE
  CHARACTER(len=200) :: UZTWM_FILE
  CHARACTER(len=200) :: UZFWM_FILE
  CHARACTER(len=200) :: UZK_FILE
  CHARACTER(len=200) :: ZPERC_FILE
  CHARACTER(len=200) :: REXP_FILE
  CHARACTER(len=200) :: LZTWM_FILE
  CHARACTER(len=200) :: LZFSM_FILE
  CHARACTER(len=200) :: LZFPM_FILE
  CHARACTER(len=200) :: LZSK_FILE
  CHARACTER(len=200) :: LZPK_FILE
  CHARACTER(len=200) :: PFREE_FILE
  CHARACTER(len=200) :: ELEV_FILE
  CHARACTER(len=200) :: SNOWBAND_FILE
  CHARACTER(len=200) :: SAC_CONST
  CHARACTER(len=200) :: PE_SCALE
  CHARACTER(len=200) :: PE_ADJ
  CHARACTER(len=200) :: INITIAL  ! Name of the initial state file
  CHARACTER(len=200) :: FORCING  ! Path and prefix of the input forcing files
  CHARACTER(len=200) :: RESTART  ! Path and prefix of the restart files
  CHARACTER(len=200) :: RESULT   ! Path of the output result files
  CHARACTER(len=50)  :: SUFFIX   ! Suffix of current input/output files (".YYYYMM.nc")
  CHARACTER(len=50)  :: SUFFIX_NEW  ! Suffix of next input/output files (".YYYYMM.nc")
  CHARACTER(len=200) :: FORCFILE ! Path and full name of the current forcing file
  CHARACTER(len=200) :: PEFILE   ! Path and full name of the current pe file
  CHARACTER(len=200) :: RESTFILE ! Path and full name of the current restart file
  CHARACTER(len=200) :: OUTFILES(7) ! Paths and full names of the current output files

  ! Input/output netcdf file ids
  INTEGER :: LSC_NCID
  INTEGER :: INITIAL_NCID
  INTEGER :: FORCING_NCID
  INTEGER :: RESTART_NCID
  INTEGER :: OUT_NCIDS(7)
  INTEGER :: SAC_OUT_NCIDS(7)

  ! Other model options
  CHARACTER(len=10) :: PARAM_TYPE
  INTEGER :: MODEL_DT
  REAL    :: MODEL_DT_REAL
  INTEGER :: OUTPUT_DT
  REAL    :: OUTPUT_DT_REAL
  INTEGER :: MADTT
  INTEGER :: OUTPUT_STEP_RATIO
  INTEGER :: YEAR0
  INTEGER :: MONTH0
  INTEGER :: DAY0
  INTEGER :: JULDAY0
  INTEGER :: YEAR_FINAL
  INTEGER :: MONTH_FINAL
  INTEGER :: DAY_FINAL
  INTEGER :: ICE1
  REAL    :: Z1
  CHARACTER(len=19) :: MODEL_START_TIME

  ! Model data dimensions
  INTEGER :: MAXNSOIL         ! Number of soil layers
  INTEGER :: MODEL_TYPE       ! Flag specifying model extensions to use
  INTEGER :: ylen             ! Number of grid rows
  INTEGER :: xlen             ! Number of grid columns
  INTEGER :: landlen          ! Number of active (Land/InBasin) cells in grid
  INTEGER :: tsteplen         ! Number of timesteps in forcing data file
  INTEGER :: NMONTHS     ! Number of timesteps in land surface characteristics file
  INTEGER :: NBANDS     ! Number of snow elevation bands

  ! Grid cell geometric attributes
  REAL    :: XLLCORNER        ! Longitude of grid cell at lower left corner of grid
  REAL    :: YLLCORNER        ! Latitude of grid cell at lower left corner of grid
  REAL    :: cellsize         ! size of grid cell, in degrees
  INTEGER, ALLOCATABLE :: ROW(:,:)     ! Grid cell row number (indexed by y,x)
  INTEGER, ALLOCATABLE :: COL(:,:)     ! Grid cell column number (indexed by y,x)
  REAL, ALLOCATABLE    :: LAT(:,:)     ! Grid cell latitude (indexed by y,x)
  REAL, ALLOCATABLE    :: LON(:,:)     ! Grid cell longitude (indexed by y,x)
  INTEGER, ALLOCATABLE :: CELLID(:,:)  ! Cell ID number (indexed by y,x)
  INTEGER, ALLOCATABLE :: LANDMASK(:,:)! 1=Land/InBasin, 0=Ocean/OutOfBasin (indexed by y,x)
  INTEGER, ALLOCATABLE :: Y(:),X(:)    ! Y and X values of active (Land/InBasin) grid cells (indexed by cell id)
  INTEGER, ALLOCATABLE :: LAND(:)      ! List of (Y*ylen+X) values of active (Land/InBasin)
                                       ! grid cells (indexed by cell id) (used in compression)

  ! Grid cell soil/landcover attributes
  LOGICAL :: LREAL8
  INTEGER, ALLOCATABLE :: NSOIL(:)
  REAL, ALLOCATABLE    :: SOILDEPTH(:,:)
  REAL, ALLOCATABLE    :: SOILDEPTH_ACCUM(:,:)
  INTEGER, ALLOCATABLE :: SOILTYP(:)
  INTEGER, ALLOCATABLE :: SLOPETYP(:)
  REAL, ALLOCATABLE    :: TBOT(:)
  INTEGER, ALLOCATABLE :: VEGTYP(:)
  REAL, ALLOCATABLE    :: SHDFAC(:,:)
  REAL, ALLOCATABLE    :: SHDMIN(:)
  REAL, ALLOCATABLE    :: SHDMAX(:)
  REAL, ALLOCATABLE    :: ALBEDO(:,:)
  REAL, ALLOCATABLE    :: SNOALB(:)
  INTEGER, ALLOCATABLE :: ICE(:)
  REAL, ALLOCATABLE    :: W_WILT(:)
  REAL, ALLOCATABLE    :: W_SAT(:)
  REAL, ALLOCATABLE    :: W_BPOWER(:)
  REAL, ALLOCATABLE    :: W_SAT_HYDC(:)
  REAL, ALLOCATABLE    :: W_SAT_MATP(:)
  REAL, ALLOCATABLE    :: ELEV_2d(:)

  ! SAC parameters - 2D
  REAL, ALLOCATABLE    :: UZTWM_2d(:)
  REAL, ALLOCATABLE    :: UZFWM_2d(:)
  REAL, ALLOCATABLE    :: UZK_2d(:)
  REAL, ALLOCATABLE    :: PCTIM_2d(:)
  REAL, ALLOCATABLE    :: ADIMP_2d(:)
  REAL, ALLOCATABLE    :: RIVA_2d(:)
  REAL, ALLOCATABLE    :: ZPERC_2d(:)
  REAL, ALLOCATABLE    :: REXP_2d(:)
  REAL, ALLOCATABLE    :: LZTWM_2d(:)
  REAL, ALLOCATABLE    :: LZFSM_2d(:)
  REAL, ALLOCATABLE    :: LZFPM_2d(:)
  REAL, ALLOCATABLE    :: LZSK_2d(:)
  REAL, ALLOCATABLE    :: LZPK_2d(:)
  REAL, ALLOCATABLE    :: PFREE_2d(:)
  REAL, ALLOCATABLE    :: SIDE_2d(:)
  REAL, ALLOCATABLE    :: RSERV_2d(:)
! Sensitivity parameters
  REAL, ALLOCATABLE    :: WCRIT_2d(:)
  REAL, ALLOCATABLE    :: PSNOW1_2d(:)
  REAL, ALLOCATABLE    :: PSNOW2_2d(:)
  REAL, ALLOCATABLE    :: RICHARDS_2d(:)
! SAC-HT
  REAL, ALLOCATABLE    :: FRZPAR(:,:)
  REAL, ALLOCATABLE    :: FRZST(:,:,:)
  REAL, ALLOCATABLE    :: SACST(:,:,:)
  REAL, ALLOCATABLE    :: SACPAR(:,:)

  ! PE Scale/adjustment factors - 2D
!  REAL, ALLOCATABLE    :: PESCALE_2d(:)
!  REAL, ALLOCATABLE    :: PEADJ_2d(:,:)

  ! SAC State Variables - 2D
  REAL, ALLOCATABLE    :: UZTWC_2d(:,:)
  REAL, ALLOCATABLE    :: UZFWC_2d(:,:)
  REAL, ALLOCATABLE    :: LZTWC_2d(:,:)
  REAL, ALLOCATABLE    :: LZFSC_2d(:,:)
  REAL, ALLOCATABLE    :: LZFPC_2d(:,:)
  REAL, ALLOCATABLE    :: ADIMC_2d(:,:)

  ! SAC STATE VARS - 1D
!  REAL :: UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC

  ! SAC PARAMETERS - 1D
  REAL :: UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP, &
       LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,SIDE,RSERV,WCRIT, &
       PSNOW1,PSNOW2,RICHARDS

  ! Anderson PEMB snow model parameters 2D and 3D
  REAL, ALLOCATABLE    :: NSNOW(:,:)
  REAL, ALLOCATABLE    :: DSNOW(:,:,:)
  REAL, ALLOCATABLE    :: PSNOW(:,:,:)
  REAL, ALLOCATABLE    :: RTTSNOW(:,:,:)
  REAL, ALLOCATABLE    :: RTTDTSNOW(:,:,:)
  REAL, ALLOCATABLE    :: WTSNOW(:,:,:)
  REAL, ALLOCATABLE    :: WTDTSNOW(:,:,:)
  REAL, ALLOCATABLE    :: TTSNOW(:,:,:)
  REAL, ALLOCATABLE    :: TTDTSNOW(:,:,:)


  ! Interpolated landcover variables
  REAL, ALLOCATABLE    :: SHDFAC_D(:)
  REAL, ALLOCATABLE    :: ALBEDO_D(:)

  ! Snow elevation band parameters
  REAL, ALLOCATABLE :: band_elev(:,:)
  REAL, ALLOCATABLE :: band_area(:,:)
  REAL, ALLOCATABLE :: band_prec(:,:)
  REAL, ALLOCATABLE :: Tfactor(:,:)
  REAL, ALLOCATABLE :: Pfactor(:,:)

  ! Grid cell state variables
  REAL, ALLOCATABLE    :: T1(:,:)
  REAL, ALLOCATABLE    :: TPACK(:,:)
  REAL, ALLOCATABLE    :: STC(:,:,:)
  REAL, ALLOCATABLE    :: SMC(:,:,:)
  REAL, ALLOCATABLE    :: PACH20(:,:)
  REAL, ALLOCATABLE    :: SH2O(:,:,:)
  REAL, ALLOCATABLE    :: CMC(:,:)
  REAL, ALLOCATABLE    :: SNOWH(:,:)
  REAL, ALLOCATABLE    :: SNEQV(:,:)
  REAL, ALLOCATABLE    :: CH(:,:)
  REAL, ALLOCATABLE    :: CM(:,:)
  REAL, ALLOCATABLE    :: SNCOVR(:,:)
  INTEGER, ALLOCATABLE :: LSTSNW(:,:)

  ! Forcing variables
  INTEGER :: FORCING_DT
  REAL    :: FORCING_DT_REAL
  INTEGER :: FORCING_START_YEAR
  INTEGER :: FORCING_START_MONTH
  INTEGER :: FORCING_START_DAY
  INTEGER :: FORCING_START_HOUR
  INTEGER :: FORCING_START_MIN
  INTEGER :: FORCING_START_SEC
  CHARACTER(len=19)    :: FORCING_START_TIME
  ! If COMP_FORCING is FALSE, data from the FORCING file
  !   will be read into 2-d arrays and then converted to the
  !   1-d compressed arrays that are used by the model.
  LOGICAL              :: COMP_FORCING
  LOGICAL              :: LLSRAINF
  LOGICAL              :: LWIND2d
  REAL, ALLOCATABLE    :: SWdownminus(:)
  REAL, ALLOCATABLE    :: LWdownminus(:)
  REAL, ALLOCATABLE    :: Tairminus(:)
  REAL, ALLOCATABLE    :: Qairminus(:)
  REAL, ALLOCATABLE    :: Rainfminus(:)
  REAL, ALLOCATABLE    :: Snowfminus(:)
  REAL, ALLOCATABLE    :: PSurfminus(:)
  REAL, ALLOCATABLE    :: Windminus(:)
  REAL, ALLOCATABLE    :: SWdownnow(:)
  REAL, ALLOCATABLE    :: LWdownnow(:)
  REAL, ALLOCATABLE    :: Tairnow(:)
  REAL, ALLOCATABLE    :: Qairnow(:)
  REAL, ALLOCATABLE    :: Rainfnow(:)
  REAL, ALLOCATABLE    :: Snowfnow(:)
  REAL, ALLOCATABLE    :: PSurfnow(:)
  REAL, ALLOCATABLE    :: Windnow(:)
  REAL, ALLOCATABLE    :: SWdownplus(:)
  REAL, ALLOCATABLE    :: LWdownplus(:)
  REAL, ALLOCATABLE    :: Tairplus(:)
  REAL, ALLOCATABLE    :: Qairplus(:)
  REAL, ALLOCATABLE    :: Rainfplus(:)
  REAL, ALLOCATABLE    :: Snowfplus(:)
  REAL, ALLOCATABLE    :: PSurfplus(:)
  REAL, ALLOCATABLE    :: Windplus(:)
  REAL, ALLOCATABLE    :: DT_TAIR(:)
  REAL, ALLOCATABLE    :: DT_SPFH(:)
  REAL, ALLOCATABLE    :: DT_PSFC(:)
  REAL, ALLOCATABLE    :: DT_WIND(:)
  REAL, ALLOCATABLE    :: DT_LWDN(:)
  REAL, ALLOCATABLE    :: DT_SOLDN(:)
  REAL, ALLOCATABLE    :: DT_SOLDNB(:)
  REAL, ALLOCATABLE    :: DT_SOLNET(:)
  REAL, ALLOCATABLE    :: DT_PRCP(:)
  REAL, ALLOCATABLE    :: DT_RAIN(:)
  REAL, ALLOCATABLE    :: DT_SNOW(:)
  REAL, ALLOCATABLE    :: Z(:)

  ! Intermediate variables (used for i/o with SFLX subroutine)
!  LOGICAL :: LADJCH
  REAL, ALLOCATABLE    :: CZMODEL(:)
  REAL, ALLOCATABLE    :: DT_SPFH_SAT(:)
  REAL, ALLOCATABLE    :: DQSDT2(:)
  REAL, ALLOCATABLE    :: TH2(:)
  REAL, ALLOCATABLE    :: SOILM(:)
  REAL, ALLOCATABLE    :: PSISAT(:)
  REAL, ALLOCATABLE    :: LWDN(:)
  REAL, ALLOCATABLE    :: LWDN_MODEL(:)
  REAL, ALLOCATABLE    :: SoilMoistTotal(:)
  REAL, ALLOCATABLE    :: SoilMoistTotal_old(:)
  REAL, ALLOCATABLE    :: SWE_old(:)
  REAL, ALLOCATABLE    :: CanopInt_old(:)
  REAL, ALLOCATABLE    :: SFCSPD(:)
  REAL, ALLOCATABLE    :: PTU(:)
  REAL, ALLOCATABLE    :: ETA(:)
  REAL, ALLOCATABLE    :: H(:)
  REAL, ALLOCATABLE    :: EVAP_total(:)
  REAL, ALLOCATABLE    :: EC1(:)
  REAL, ALLOCATABLE    :: EDIR1(:)
  REAL, ALLOCATABLE    :: ET(:,:)
  REAL, ALLOCATABLE    :: ETT1(:)
  REAL, ALLOCATABLE    :: ESNOW(:)
  REAL, ALLOCATABLE    :: DRIP(:)
  REAL, ALLOCATABLE    :: DEW(:)
  REAL, ALLOCATABLE    :: BETA(:)
  REAL, ALLOCATABLE    :: ETP(:)
  REAL, ALLOCATABLE    :: S(:)
  REAL, ALLOCATABLE    :: FLX1(:)
  REAL, ALLOCATABLE    :: FLX2(:)
  REAL, ALLOCATABLE    :: FLX3(:)
  REAL, ALLOCATABLE    :: SNOMLT(:)
  REAL, ALLOCATABLE    :: RUNOFF1(:)
  REAL, ALLOCATABLE    :: RUNOFF2(:)
  REAL, ALLOCATABLE    :: RUNOFF3(:)
  REAL, ALLOCATABLE    :: RC(:)
  REAL, ALLOCATABLE    :: PC(:)
  REAL, ALLOCATABLE    :: RSMIN(:)
  REAL, ALLOCATABLE    :: RCS(:)
  REAL, ALLOCATABLE    :: RCT(:)
  REAL, ALLOCATABLE    :: RCQ(:)
  REAL, ALLOCATABLE    :: RCSOIL(:)
  REAL, ALLOCATABLE    :: SOILW(:)
  REAL, ALLOCATABLE    :: SMCDRY(:)
  REAL, ALLOCATABLE    :: SMCREF(:)
  INTEGER, ALLOCATABLE :: NROOT(:)
  REAL, ALLOCATABLE    :: RM(:)
  REAL, ALLOCATABLE    :: ACondMax(:)
  REAL, ALLOCATABLE    :: band_area_with_snow(:)

  ! Output variables
  INTEGER :: VARIDS(7,20)
  ! If COMP_OUTPUT is FALSE, output data will be
  !   converted from the 1- or 2-d compressed arrays that are used by the model
  !   into 2- or 3-d arrays that are written to the OUTPUT files
  LOGICAL              :: COMP_OUTPUT
  REAL, ALLOCATABLE    :: SWnet(:)
  REAL, ALLOCATABLE    :: LWnet(:)
  REAL, ALLOCATABLE    :: Qle(:)
  REAL, ALLOCATABLE    :: Qh(:)
  REAL, ALLOCATABLE    :: Qg(:)
  REAL, ALLOCATABLE    :: Qf(:)
  REAL, ALLOCATABLE    :: Qv(:)
  REAL, ALLOCATABLE    :: Qa(:)
  REAL, ALLOCATABLE    :: DelSurfHeat(:)
  REAL, ALLOCATABLE    :: DelColdCont(:)
  REAL, ALLOCATABLE    :: Snowf(:)
  REAL, ALLOCATABLE    :: Rainf(:)
  REAL, ALLOCATABLE    :: Evap(:)
  REAL, ALLOCATABLE    :: Qs(:)
  REAL, ALLOCATABLE    :: Qsb(:)
  REAL, ALLOCATABLE    :: Qsm(:)
  REAL, ALLOCATABLE    :: Qfz(:)
  REAL, ALLOCATABLE    :: Qst(:)
  REAL, ALLOCATABLE    :: DelSoilMoist(:)
  REAL, ALLOCATABLE    :: DelSWE(:)
  REAL, ALLOCATABLE    :: DelSurfStor(:)
  REAL, ALLOCATABLE    :: DelIntercept(:)
  REAL, ALLOCATABLE    :: SnowT(:)
!Ben Livneh track effective snowpack temp TPACK = SnowTProf
  REAL, ALLOCATABLE    :: SnowTProf(:)
  REAL, ALLOCATABLE    :: VegT(:)
  REAL, ALLOCATABLE    :: BareSoilT(:)
  REAL, ALLOCATABLE    :: AvgSurfT(:)
  REAL, ALLOCATABLE    :: RadT(:)
  REAL, ALLOCATABLE    :: Albedo_ALMA(:)
  REAL, ALLOCATABLE    :: SWE(:)
  REAL, ALLOCATABLE    :: UZTWC(:)
  REAL, ALLOCATABLE    :: UZFWC(:)
  REAL, ALLOCATABLE    :: LZTWC(:)
  REAL, ALLOCATABLE    :: LZFPC(:)
  REAL, ALLOCATABLE    :: LZFSC(:)
  REAL, ALLOCATABLE    :: PackWater(:)
  REAL, ALLOCATABLE    :: SliqFrac(:)
  REAL, ALLOCATABLE    :: SWEVeg(:)
  REAL, ALLOCATABLE    :: SurfStor(:)
  REAL, ALLOCATABLE    :: SoilMoist(:,:)
  REAL, ALLOCATABLE    :: SoilMoistSac(:,:)
  REAL, ALLOCATABLE    :: SoilTemp(:,:)
  REAL, ALLOCATABLE    :: SMLiqFrac(:,:)
  REAL, ALLOCATABLE    :: SMFrozFrac(:,:)
  REAL, ALLOCATABLE    :: SoilWet(:)
  REAL, ALLOCATABLE    :: PotEvap(:)
  REAL, ALLOCATABLE    :: ECanop(:)
  REAL, ALLOCATABLE    :: TVeg(:)
  REAL, ALLOCATABLE    :: ESoil(:)
  REAL, ALLOCATABLE    :: EWater(:)
  REAL, ALLOCATABLE    :: RootMoist(:)
  REAL, ALLOCATABLE    :: CanopInt(:)
  REAL, ALLOCATABLE    :: EvapSnow(:)
  REAL, ALLOCATABLE    :: SubSnow(:)
  REAL, ALLOCATABLE    :: SubSurf(:)
  REAL, ALLOCATABLE    :: ACond(:)
  REAL, ALLOCATABLE    :: SnowFrac(:)
  REAL, ALLOCATABLE    :: Fdepth(:)
  REAL, ALLOCATABLE    :: Tdepth(:)
  REAL, ALLOCATABLE    :: SAlbedo(:)
  REAL, ALLOCATABLE    :: SatSoil(:)
  REAL, ALLOCATABLE    :: WltSoil(:)
  REAL, ALLOCATABLE    :: SnowDepth(:)
  REAL, ALLOCATABLE    :: PotEvapPE(:,:)

END MODULE driverMod

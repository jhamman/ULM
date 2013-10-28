SUBROUTINE READCNTL()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS A CONTROL FILE

  ! Modifications:
  ! 2008-May-05 Removed LADJCH, for compatibility with NOAH 2.8.		Ben Livneh
  ! 2008-Jul-15 Added LADJCH again.						Ben Livneh
  ! 2008-Oct-07 Included SAC parameters for use as unified model         Ben Livneh

  ! driverMod contains definitions of all global variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READCNTL.f90,v 1.4 2008/07/30 22:33:03 vicadmin Exp $"/

  ! Define local variables
  NAMELIST /CONTROL/MODEL_DT,OUTPUT_DT,YEAR0,MONTH0,DAY0,YEAR_FINAL, &
                    MONTH_FINAL,DAY_FINAL,PARAM_TYPE,LSC,LREAL8, &
                    MASK_FILE,NSOIL_FILE,SOILDEPTH_FILE,SOILTYP_FILE, &
                    SLOPETYP_FILE,TBOT_FILE,VEGTYP_FILE,SHDFAC_FILE, &
                    ALBEDO_FILE,SNOALB_FILE,UZTWM_FILE,UZFWM_FILE, &
                    UZK_FILE,ZPERC_FILE,REXP_FILE,LZTWM_FILE, &
                    LZFSM_FILE,LZFPM_FILE,LZSK_FILE,LZPK_FILE, &
                    PFREE_FILE,ELEV_FILE,SNOWBAND_FILE,NBANDS, &
                    MODEL_TYPE,MAXNSOIL,NMONTHS,SAC_CONST,ICE1,Z1, &
                    INITIAL,FORCING,RESTART,RESULT,COMP_OUTPUT

  OPEN (UNIT=99, FILE=CNTRFL, STATUS='OLD')
  READ (99,CONTROL)
  CLOSE(99)

END

!  NAMELIST /CONTROL/MODEL_DT,OUTPUT_DT,YEAR0,MONTH0,DAY0,YEAR_FINAL, &
!                    MONTH_FINAL,DAY_FINAL,PARAM_TYPE,LSC,LREAL8, &
!                    MASK_FILE,NSOIL_FILE,SOILDEPTH_FILE,SOILTYP_FILE, &
!                    SLOPETYP_FILE,TBOT_FILE,VEGTYP_FILE,SHDFAC_FILE, &
!                    ALBEDO_FILE,SNOALB_FILE,UZTWM_FILE,UZFWM_FILE, &
!                    UZK_FILE,ZPERC_FILE,REXP_FILE,LZTWM_FILE, &
!                    LZFSM_FILE,LZFPM_FILE,LZSK_FILE,LZPK_FILE, &
!                    PFREE_FILE,ELEV_FILE,SNOWBAND_FILE,NBANDS, &
!                    MODEL_TYPE,MAXNSOIL,NMONTHS,SAC_CONST,PE_SCALE,PE_ADJ,ICE1,Z1, &
!                    LADJCH,INITIAL,FORCING,RESTART,RESULT,COMP_OUTPUT

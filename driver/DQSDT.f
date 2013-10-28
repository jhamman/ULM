      function dqsdt( t, p ) 

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     name: delta qs delta t
c
c     purpose:
c     =======
c     to retrieve the appropriate value of dqsdt (the change
c     of the saturation mixing ratio with respect to the 
c     change in temperature)
c
c     method:
c     ======
c     - if temperature is within an acceptable range...
c       - calculate dqsdt according to rogers and yau, chapter 2.
c     - else set dqsdt to zero.
c
c     process narrative:  flux3 - located in the flux3 sdf in dnxm
c     =================
c    
c     references
c     ==========
c     rogers and yau, 1989, a short course in clouds physics
c
c     called from:
c     ===========
c     flxcor
c
c     interface:
c     =========  
c
c       input variables  
c       ---------------  
c       t, p
c    
c       output variables 
c       ----------------
c       dqsdt
c
c     files accessed:
c     ============== 
c     filename/unit#           r/w               description 
c     ------------------------ --- -------------------------------------
c     none
c
c     remarks: none
c     =======
c
c     variables:
c     =========
c     label     .......................description......................
c
c     cpv       specific heat of water vapor (K kg-1 K-1)
c     cw        specific heat of liquid water (J kg-1 K-1)
c     desdt     change of saturation vapor pressure with respect to
c               temperature
c     dqsdt     change of saturation mixing ration with respect to
c               temperature
c     es        saturation vapor pressure (Pa)
c     eso       saturation vapor pressure at 0 Celsius (Pa)
c     eps       ratio of gas constant of dry air and the gas
c               constant for water vapor
c     lw        latent heat of vaporization  (J kg-1)
c     p         pressure (Pa)
c     rv        gas constant for water vapor (J K-1 kg-1)
c     t         temperature (K)
c     to        kelvin/celsius conversion
c
c     updates:
c     =======
c     30 jun 97   initial version.................grunmann/mitchell/ncep
c     08 sep 99   ported to ibm sp-2.  updated prolog.  added intent
c                 attributes to argument...................mr gayno/dnxm
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      implicit none


C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: DQSDT.f,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

      real,      parameter         :: cpv = 1870.0
      real,      parameter         :: cw  = 4187.0
      real,      parameter         :: eps = 0.622
      real,      parameter         :: eso = 611.2
      real,      parameter         :: rv  = 461.5
      real,      parameter         :: to  = 273.15

      real                         :: desdt
      real                         :: dqsdt
      real                         :: es
      real                         :: lw
      real,      intent(in)        :: p
      real,      intent(in)        :: t

c     ------------------------------------------------------------------
c     executable code starts here...evaluate clausius-clapeyron
c     equation.  see eq 2.10, rogers and yau.
c     ------------------------------------------------------------------
      
      if ( (t .ge. 173.0) .and. (t .le. 373.0) ) then

        lw    = 2.501e+06 - ( cw - cpv ) * ( t - to )
        es    = eso * exp (lw * (1 / to - 1 / t) / rv )  
        desdt = lw * es / (rv * t * t)

c     ------------------------------------------------------------------
c       calculate dqsdt:  dqsdt = dqs/p , where dqs = eps*desdt  
c     ------------------------------------------------------------------
c
        dqsdt = eps * desdt / p

      else

        dqsdt = 0.0

      end if

      return

      end


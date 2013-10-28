      function e (t) 

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     name: saturation vapor pressure
c     ====
c
c     purpose:
c     =======
c     to calculate values of saturation vapor pressure (e)
c
c     method:
c     ======
c     - read in temperature in kelvin
c     - calculate saturation vapor pressure in pascals according to
c       the clausius-clapyron equation.
c
c     process narrative:  flux3 - located in the flux3 sdf in dnxm
c     =================
c    
c     references
c     ==========
c     rogers and yau, 1989, a short course in clouds physics
c
c     called from:  svp
c     ===========
c
c     interface:
c     =========  
c
c       input variables  
c       ---------------  
c       t
c    
c       output variables 
c       ----------------
c       e
c
c     files accessed:
c     ============== 
c     filename/unit#           r/w               description 
c     ------------------------ --- --------------------------------------
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
c     e         saturation vapor pressure (pascals)
c     eso       saturation vapor pressure at 0 Celsius (pascals)
c     lw        latent heat of vaporization  (J kg-1)
c     rv        gas constant for water vapor (J K-1 kg-1)
c     t         temperature (kelvin)
c     to        kelvin/celsius conversion
c
c     updates:
c     =======
c     09 jul 1997  initial version................pablo j. grunmann/ncep
c     03 sep 1999  ported to ibm sp-2......................mr gayno/dnxm
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      implicit none

C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: E.f,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

      real,  parameter              :: cpv = 1870.0
      real,  parameter              :: cw  = 4187.0
      real,  parameter              :: eso = 611.2
      real,  parameter              :: rv  = 461.5
      real,  parameter              :: to  = 273.15      

      real                          :: e
      real                          :: t
      real                          :: lw

c     ------------------------------------------------------------------
c     executable code starts here ...
c     clausius-clapeyron: see rogers and yau page 14.
c     ------------------------------------------------------------------
     
      lw = 2.501e6 - ( cw - cpv ) * ( t - to )
      e  = eso * exp(lw *(1/to - 1/t)/rv)  
    
      return

      end


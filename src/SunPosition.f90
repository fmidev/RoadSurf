!Calculates sune levation and azimuth angle
Subroutine SunPosition(modelInput,localParam,i,elevation_angle,azimuth_angle)
      use RoadSurfVariables
      Implicit None
  
          type(inputArrays),intent(IN) ::modelInput !< model input arrays
          type(localParameters),intent(IN) :: localParam
          integer, intent(IN) ::i      !< model input index
          real(8) :: JDE                !< Julian Empheresis day
          real(8),intent(OUT) ::elevation_angle       !< solar elevation angle
          real(8),intent(OUT) :: azimuth_angle        !< solar azimuth angle
          call JulianEphemerisDay(modelInput,i,JDE)
          call calcElevationAzimuth(JDE,localParam%lat,localParam%lon,elevation_angle,&
               azimuth_angle)
end Subroutine


subroutine calcElevationAzimuth(JDE,lat,lon,elevation_angle,azimuth_angle)
!     This routine is based on the code written by Leif Backman, FMI
!
!     Purpose
!     -------
!
!     The purpose of this subroutine is to calculate the solar elevation 
!     angle and azimuth angle .
!     
!     Method
!     ------
!
!     Jean Meeus, Astronomical Algorithms
!     Chapters:
!       - 11, Sidereal Time at Greenwich
!       - 12, Transformation of Coordinates
!       - 21, Nutation and the Obliquity of the Ecliptic
!       - 24, Solar Coordinates
!
!     cos(solar zenith angle) = sin(latitude)*sin(solar declinationination)
!             + cos(latitude)*cos(solar declinationination)*cos(hour angle)
!     solar elevation angle = 90 - solar zenith angle
!
!     Azimuth angle calculated as in 
!     https://www.pveducation.org/pvcdrom/properties-of-sunlight/azimuth-angle
!
!     cos(solar azimuth angle)=(sin(solar declinationination)*cos(latitude)-
!                               cos(solar declinationination)*sin(latitude)*cos(hour_angle))
!                               /cos(sun elevation)
!
!     If elevation angle <=0, elevation angle and azimuth angle are set to
!     -9999.9
!
!
!-----------------------------------------------------------------------

      implicit none

  
      real*8, intent(IN) :: JDE   !< Julian Ephemeris Day
      real*8, intent(IN) :: lat,lon !< latitude and longtitude
      real*8, intent(INOUT) :: elevation_angle  !< solar elevation angle
      real*8, intent(INOUT) :: azimuth_angle !< solar azimuth angle
      real*8, parameter :: pi  = 4 * atan (1.0_8)
      real*8 :: Dyr, cos_lat,sin_lat,lat_radians
      real*8 :: cos_declination,sin_declination,cos_dec_lat,sin_dec_lat
      real*8 :: hour_angle_corr,cosah,cos_elev,chi,cosele,precos
      real*8 :: T,ml,ma,ecc,sunc,al,tilt,eps,ra,declination,stG
  
  !----constants
      Dyr = 365.25

!----T, Julian centuries since J2000.0
      T=(JDE-2451545.0)/(Dyr*100.)

!----ml, geometric mean longitude of the Sun
      ml=280.46645 +36000.76983*T +0.0003032*T*T 

      if ( ml.lt.0. ) ml=ml-360.*(AINT(ml/360.)-1.)
      if ( ml.gt.360. ) ml=ml-360.*AINT(ml/360.)

!----ma, mean anomaly of the Sun
      ma=357.52910 +35999.05030*T -0.0001559*T*T -0.00000048*T*T*T

      if ( ma.lt.0. ) ma=ma-360.*(AINT(ma/360.)-1.)
      if ( ma.gt.360. ) ma=ma-360.*AINT(ma/360.)

!----ecc, eccentricity of the earths orbit
      ecc=0.016708617 -0.000042037*T -0.0000001236*T*T

!----sunc, Sun's equation of center
      sunc=(1.913600 -0.004817*T -0.000014*T*T)*sin(ma*pi/180.) &
           +(0.019993 -0.000101*T)*sin(2.*ma*pi/180.)&
           +0.000290*sin(3.*ma*pi/180.)

!----al, apparent longitude of the Sun
      al=ml +sunc -0.00569 -0.00478*sin((125.04-1934.136*T)*pi/180.)
      al=al*pi/180.

!----tilt, mean obliquity of the ecliptic
      tilt= 23.43929111 -0.013004166*T -0.001638888*T*T &
           +0.005036111*T*T*T
      eps=tilt +0.00256*cos((125.04-1934.136*T)*pi/180.)
      eps=eps*pi/180.

!----ra, apparent right ascension of the Sun
      ra=atan2(cos(eps)*sin(al),cos(al))

      if ( ra.lt.0. ) ra=ra-2.*pi*(AINT(ra/(2.*pi))-1.)
      if ( ra.gt.2.*pi ) ra=ra-2.*pi*AINT(ra/(2.*pi))

!----declination, apparent declinationination of the Sun
      declination=asin(sin(eps)*sin(al))

!----stG, mean sidereal time at Greenwich
      stG=280.46061837+360.98564736629*(JDE-2451545.0)&
           +0.000387933*T*T-T*T*T/38710000.

      if ( stG.lt.0. ) stG=stG-360.*(AINT(stG/360.)-1.)
      if ( stG.gt.360. ) stG=stG-360.*AINT(stG/360.)

      stG=stG*pi/180.

!----trigonometry
      cos_declination = cos(declination)
      sin_declination = sin(declination)
      lat_radians=pi*lat/180.
      sin_lat=sin(lat_radians)
      cos_lat=cos(lat_radians)
      cos_dec_lat = cos_declination*cos_lat
      sin_dec_lat = sin_declination*sin_lat
  
!----hour_angle_corr, local hour angle
      hour_angle_corr = (stG +lon*pi/180. -ra)
      if ( ra.lt.0. ) hour_angle_corr=hour_angle_corr-2.*pi*(AINT(hour_angle_corr/(2.*pi))-1.)
      if ( ra.gt.2.*pi ) hour_angle_corr=hour_angle_corr-2.*pi*AINT(hour_angle_corr/(2.*pi))
  
!---- solar zenith angle
      cosah = cos(hour_angle_corr)
      cos_elev = sin_dec_lat + cos_dec_lat*cosah

      if (cos_elev.ge.1.0.and.cos_elev.lt.1.001) then
         cos_elev = 1.0000000000
         chi = 0.
      else if(cos_elev.ge.1.001) then
          write(*,*) 'Problem with zenith angle at lat/lon:',lat,lon
          stop
      else if(cos_elev.gt.-1.001.and.cos_elev.le.-1.0) then
          cos_elev = -1.0000000000
          chi = pi
      else
          chi = acos(cos_elev)
      endif

       elevation_angle = 90.0-chi*(180./pi)

!----hour angle correction       
       if (hour_angle_corr < 0.) then
           hour_angle_corr=2*pi+hour_angle_corr
       else if (hour_angle_corr > 2*pi) then
           hour_angle_corr=hour_angle_corr-2*pi
       end if

!----solar azimuth angle
       !Calculate ongly if elevation angle > 0
       if (elevation_angle>0) Then
          cosele=cos((pi/2.0)-chi)
          precos=0.0
          if (cosele .ge. -0.0001 .and. cosele .lt. 0.0001)Then
             azimuth_angle=-9999.9
          else 
             precos=(sin_declination*cos_lat-cos_declination*sin_lat*cosah)/cosele
             if (precos .ge. 1.0 .and. precos .lt. 1.001) then
                 precos=1.00000000
                 azimuth_angle=0.0
             
             else if(precos.ge.1.001) then
                 write(*,*) 'Problem with azimuth angle at lat/lon:',lat,lon
                 stop
             else if(precos.gt.-1.001.and.precos.le.-1.0) then
                 precos = -1.0000000000
                 azimuth_angle = pi
             else
                 azimuth_angle=acos(precos)
             end if
          end if
          if (hour_angle_corr<pi) then
             azimuth_angle=2*pi-azimuth_angle
          end if
          azimuth_angle=azimuth_angle*(180./pi)
       else
          azimuth_angle=-9999.9
          elevation_angle=-9999.9
       end if
end subroutine

subroutine JulianEphemerisDay(modelInput,idx,JDE)

!     Purpose
!     -------
!
!     The purpose of this subroutine is to calculate  the Julian Ephemeris Day.
!
!
!     Parameters  (i=integer, r=real)
!     ----------
!
!     juday  (i) = julian day
!     JDE    (r) = Julian Ephemeris Day (changed)
!
!
!     Results
!     -------
!     
!     Calculated result is stored in  JDE.
!
!
!     Method
!     ------
!
!     Jean Meeus, Astronomical Algorithms
!     Chapters:
!       - 7, Julian Day

      use RoadSurfVariables
      Implicit None
  
      type(inputArrays),intent(IN) ::modelInput !< model input arrays
      integer, intent(IN) ::idx     !< index in input data
      real*8, intent(OUT) :: JDE !< Julian Empheris Day	  
      integer :: mmsec,mmmin,mmhr,mmday,mmmon,mmyr
      real*8 :: yr,mo,A,B,day,Dyr
  
      mmyr=modelInput%year(idx)
      mmmon=modelInput%month(idx)
      mmday=modelInput%day(idx)
      mmhr=modelInput%hour(idx)
      mmmin=modelInput%minute(idx)
      mmsec=modelInput%second(idx)
  
!----Constants
      Dyr = 365.25

      !----Julian Ephemeris Day
      if ( mmmon.le.2 ) then
         yr=REAL(mmyr-1)
         mo=REAL(mmmon+12)
      else
         yr=REAL(mmyr)
         mo=REAL(mmmon)
      endif

      day =  REAL(mmday)+ REAL(mmhr)/24.+ REAL(mmmin)/(24.*60.)+ &
         REAL(mmsec)/(24.*60.*60.)

      A=DINT(yr/100.)
      B=2.-A+DINT(A/4.)

      JDE=DINT(Dyr*(yr+4716))+DINT(30.6001*(mo+1.))+day+B-1.5245D3

end subroutine

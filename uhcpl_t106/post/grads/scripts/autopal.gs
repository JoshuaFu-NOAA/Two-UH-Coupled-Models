function palette(cvar)

  say 'cvar = 'cvar

  if (cvar='code4' | cvar='precip')
    'set clevs 1 2 3 4 5 6 7 8 10 12 14 16 18 20'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
*   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31 33 35 37 39 41 43 45'
*   'set ccols 3 18  22  26  30  32  34   38  42'
*   'set ccols 0 25 41 40 39 38 37 36 35 34 33 32 31 30 29 28'
   'set ccols 0 28 26 24 22 20 18 16 45 42 40 38 35 33 31'

 return

  endif

  if (cvar='code77' | cvar='rel')
*    'set clevs 0 200 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800'
    'set clevs -500 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
*   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31 33 35 37 39 41 43 45'
*   'set ccols 3 18  22  26  30  32  34   38  42'
   'set ccols 25 41 40 39 38 37 36 35 34 33 32 31 30 29 28'

 return

  endif



  if (cvar='code173' | cvar='z0veg' | cvar='z0')
*    'set clevs 0 200 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800'
    'set clevs 1 5 30 100 250 1000 2000'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
*   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31 33 35 37 39 41 43 45'
*   'set ccols 3 18  22  26  30  32  34   38  42'
   'set ccols 18  22  26  30  32  34   38  42'

 return

  endif

*
* ********************* surface geopotential *************************************
*
if (cvar='code129' | cvar='GEOSP' | cvar='geosp')
   'set clevs -200 0 200 400 600 800 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000'
   colnum=subwrd(result,5)
   say 'Number of colors = 'colnum
   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'
   return
   'set clevs 0 200 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000 4200 4400 4600 4800 5000 5200 5400'
   colnum=subwrd(result,5)
   say 'Number of colors = 'colnum
   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'
   return
endif
*
* ********************* surface pressure ( hPa ) **********************************
*
  if (cvar='ps' | cvar='code134')
    'set clevs 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000 1025 1050'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38'

 return

  endif
*
* ********************** surface pressure ( Pa ) ***********************************
*
  if (cvar='APS' | cvar='aps')
    'set clevs 50000 52500 55000 57500 60000 62500 65000 67500 70000 72500 75000 77500 80000 82500 85000 87500 90000 92500 95000 97500 100000 102500 105000'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40'

 return

  endif
*
* *********************** surface temperature ( Celsius ) **************************
*
  if (cvar='ts' | cvar='code139'| cvar='code169' | cvar='TSURF' | cvar='TS' | cvar='tsurf')
    'set clevs -35 -30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30 35'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48'

 return

  endif
*
* *************************** soil wetness ( m ) *************************************
*
  if (cvar='code140' | cvar='WS' | cvar='ws')
    'set clevs 0 10 20 30 40 50 60 70 80 90'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
   'set ccols 17 20 22 26 28 32 36 38 40 42 44'

 return

  endif
*
* ********************************* snow depth ( cm ) *********************************
*
  if (cvar='code141' | cvar='SN' | cvar='sn')
    'set clevs 0 1 2 3 4 5 6 7 8 9 10 15 20 50 100'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37'

 return

  endif
*
* ************************ large scale precipitation ( m/s) **************************
*
  if (cvar='code142' | cvar='APRL' | cvar='aprl')
    'set clevs 0 0.5e-08 1e-08 1.5e-08 2e-08 2.5e-08 3e-08 3.5e-08 4e-08 4.5e-08 5e-08 5.5e-08 6e-08 6.5e-08 7e-08 7.5e-08 8e-08 8.5e-08 9e-08 9.5e-08 1e-07 1.05e-07 1.1e-07 1.15e-07 1.2e-07 1.25e-07 1.3e-07'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44'

 return

  endif
*
* ************************* convective precipitation ( m/s ) **************************
*
  if (cvar = 'code143' | cvar='APRC' | cvar='aprc')
    'set clevs 0 2e-08 4e-08 6e-08 8e-08 1e-07 1.1e-07 1.2e-07 1.3e-07 1.4e-07 1.5e-07 1.6e-07 1.7e-07 1.8e-07 1.9e-07 2e-07 2.1e-07 2.2e-07 2.3e-07 2.4e-07 2.5e-07 2.6e-07 2.7e-07 2.8e-07 2.9e-07 3e-07 3.1e-07'
  colnum=subwrd(result,5)
  say 'Number of colors = 'colnum
  'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44'

 return

  endif
*
* ******************************* snow fall ( m/s ) ************************************
*
  if (cvar = 'code144' | cvar='APRS' | cvar='aprs')
    'set clevs 0 0.5e-08 1e-08 1.5e-08 2e-08 2.5e-08 3e-08 3.5e-08 4e-08 4.5e-08 5e-08 5.5e-08 6e-08 6.5e-08 7e-08 7.5e-08 8e-08 8.5e-08 9e-08 9.5e-08 1e-07 1.05e-07 1.1e-07 1.15e-07 1.2e-07 1.25e-07 1.3e-07'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44'

 return

  endif
*
* *********************** surface sensible heat flux ( W/m**2 ) ************************
*
  if (cvar = 'code146' | cvar='AHFS' | cvar='ahfs')
    'set clevs -260 -240 -220 -200 -180 -160 -140 -120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120 140 160'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38'

 return

  endif
*
* *********************** surface latent heat flux ( W/m**2 ) **************************
*
  if (cvar = 'code147' | cvar='AHFL' | cvar='ahfl')
    'set clevs -300 -285 -270 -255 -240 -225 -210 -195 -180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 15 30'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39'

 return

  endif
*
* **************************** surface runoff ( m/s ) **********************************
*
  if (cvar = 'code160' | cvar='RUNOFF' | cvar='runoff')
    'set clevs 0 1e-08 2e-08 3e-08 4e-08 5e-08 6e-08 7e-08 8e-08 9e-08 1e-07 1.1e-07 1.2e-07 1.3e-07 1.4e-07 1.5e-07 1.6e-07 1.7e-07 1.8e-07 1.9e-07 2e-07'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38'

 return

  endif
*
* *************************** total cloud cover ( fract. ) *****************************
*
  if (cvar = 'code164' | cvar='ACLCOV' | cvar='aclcov')
    'set clevs 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 20 22 23 24 26 28 30 31 33 35 37 39 40 41 43 44'

 return

  endif
*
* **************************** 10m u-velocity ( m/s ) ************************************
*
  if (cvar = 'code165' | cvar='U10' | cvar='u10')
    'set clevs -8 -7 -6 -5 -4 -3 -2 -1 0- 1 2 3 4 5 6 7 8 9 10'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 20 22 23 24 26 27 28 30 31 33 35 36 37 39 40 41 43 44'

 return

  endif
*
* *************************** 10m v-velocity ( m/s ) **************************************
*
  if (cvar = 'code166' | cvar='V10' | cvar='v10')
    'set clevs -6 -5.5 -5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'

 return

  endif
*
* ****************************** 2m temperature ( Kelvin ) *******************************
*
  if (cvar = 'code167')
    'set clevs 230 235 237.5 240 242.5 245 247.5 250 252.5 255 257.5 260 262.5 265 267.5 270 272.5 275 277.5 280 282.5 285 287.5 290 292.5 300 302.5 305 310'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* **************************** 2m temperature ( Celsius ) *********************************
*
  if (cvar='TEMP2' | cvar='temp2')
    'set clevs -45 -40 -35 -30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 1. 12.5 15 17.5 20 22.5 25 27.5 30 35'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* **************************** deep soil temperature ( Celsius ) ***************************
*
  if (cvar='code170' | cvar='TD' | cvar='td')
  'set clevs -40 -35 -30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ***************************** 10m windspeed ( m/s ) ***************************************
*
  if (cvar='code171' | cvar='WIND10' | cvar='wind10')
    'set clevs 0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.011 0.012 0.013 0.014 0.015 0.016 0.017 0.018 0.019 0.02'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 19 20 21 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 41 42 43 44 45 46'

 return

  endif
*
* ******************************* land sea mask (0:sea ; 1:land ) ***************************
*
  if (cvar='code172' | cvar='SLM' | cvar='slm')
    'set clevs 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 19 20 21 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 41 42 43 44 45 46'

 return

  endif
*
* ******************************* surface roughness length ( m ) ****************************
* 
  if (cvar='code173' | cvar='AZ0' | cvar='az0')
    'set clevs 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 19 20 21 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 41 42 43 44 45 46'

 return

  endif
*
* *********************** surface background albedo ( fract. ) *******************************
*
  if (cvar='code174' | cvar='ALB' | cvar='alb')
    'set clevs 0 0.05 0.10.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 41 42 43 44 45 46'

 return

  endif
*
* ************************ net surface solar radiation ( W/m**2 ) ****************************
*
  if (cvar='code176' | cvar='SRADS' | cvar='srads')
    'set clevs 0 20 40 60 80 100 120 140 160 180 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ************************ net surface thermal radiation ( W/m*2 ) *************************
*
  if (cvar='code177' | cvar='TRADS' | cvar='trads')
    'set clevs -140 -135 -130 -125 -120 -115 -110 -105 -100 -95 -90 -85 -80 -75 -70 -65 -60 -55 -50 -45 -40 -35 -30 -25 -20 -15 -10 -5'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return 

  endif
*
* ************************ net top solar radiation ( W/m**2 ) ******************************
* 
  if (cvar='code178' | cvar='SRAD0' | cvar='srad0')
    'set clevs 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 19 21 22 23 24 25 26 27 29 30 31 33 34 35 37 38 39 41 42 43 44'

 return

  endif
*
* ************************* top thermal radiation (OLR) ( W/m**2 ) *************************
*    
  if (cvar='code179' | cvar='TRAD0' | cvar='trad0')
    'set clevs -300 -290 -280 -270 -260 -250 -240 -230 -220 -210 -200 -190 -180 -170 -160'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 24 26 28 30 32 34 36 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* *********************** surface u-stress ( Pa ) ******************************************
*
  if (cvar='code180' | cvar='USTR' | cvar='ustr')
    'set clevs -0.5 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 19 20 22 23 25 26 28 29 31 32 34 35 37 38 40 41 42 43 44 45 46'

 return

  endif
*
* ************************ surface v-stress ( Pa ) *****************************************
*
  if (cvar='code181' | cvar='VSTR' | cvar='vstr')
    'set clevs -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 19 21 23 24 26 28 29 31 33 34 36 38 40 41 42 43 44 45'

 return

  endif
*
* ***************************** soil temperature ( Celsius ) *******************************
*
  if (cvar='code183' | cvar='TDCL' | cvar='tdcl')
  'set clevs -40 -35 -30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ***************************** forest cover ( fract. ) ************************************
*
  if (cvar='code184' | cvar='FOREST' | cvar='forest')
    'set clevs 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ********************** east west orographic variance ( cm**2 ) ****************************
*
  if (cvar='code190' | cvar='EWOV' | cvar='ewov')
    'set clevs 2 4 6 8 10 15 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 160 170'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ********************** north south orographic variance ( cm**2 ) **************************
*
  if (cvar='code191' | cvar='NSOV' | cvar='nsov')
    'set clevs 2 4 6 8 10 15 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 160 170'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ****************** northwest southeast orographic variance ( cm**2 ) ***********************
*
  if (cvar='code192' | cvar='NWOV' | cvar='nwov')
    'set clevs 2 4  6 8 10 15 20 25 30 35 40 45 50 60 70 80 90 100 120 130 140 150 160 170'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ****************** northeast southwest orographic variance ( cm**2 ) **********************
*
  if (cvar='code193' | cvar='NEOV' | cvar='neov')
    'set clevs 2 4 6 8 10 15 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 160 170'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ********************** skin reservoir contents ( mm ) *************************************
*
  if (cvar='code198' | cvar='SRC' | cvar='src')
    'set clevs 0 0.00005 0.0001 0.00015 0.0002 0.00025 0.0003 0.00035 0.0004 0.00045 0.0005 0.00055 0.0006 0.00065 0.0007 0.00075 0.0008'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 20 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ************************** vegetation ratio ( *100=% ) ************************************
*
  if (cvar='code199' | cvar='VEG' | cvar='veg')
    'set clevs 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ******************** variance of subgrid scale ( m**2 ) ***********************************
*
  if (cvar='code200' | cvar='VSO' | cvar='vso')
*    'set clevs 1e05 2e05 3e05 4e05 5e05 6e05 7e05 8e05 9e05 1e06 1.1e06 1.2e06 1.3e06 1.4e06 1.5e06 1.6e06 1.7e06 1.8e06 1.9e06'
*    colnum=subwrd(result,5)
*    say 'Number of colors = 'colnum
*    'set ccols 16 20 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'
    'set clevs 2 4 6 8 10 15 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 160 170'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46'

 return

  endif
*
* ********************** top solar radiation upward ( W/m**2 ) *******************************
*
  if (cvar='code203' | cvar='STRAD0U' | cvar='strad0u')
    'set clevs -350 -330 -310 -290 -270 -250 -230 -210 -200 -190 -180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48'

 return

  endif
*
* ********************** surface solar radiation upward ( W/m**2 ) **************************
*
  if (cvar='code204' | cvar='SRADSU' | cvar='sradsu')
    'set clevs -350 -330 -310 -290 -270 -250 -230 -210 -190 -170 -150 -130 -110 -100 -90 -80 -70 -50 -40 -30 -20 -10'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48'

 return

  endif
*
* *********************** surface thermal radiation upward ( W/m**2 ) ************************
*
  if (cvar='code205' | cvar='TRADSU' | cvar='tradsu')
    'set clevs -500 -480 -460 -450 -440 -430 -420 -410 -400 -390 -380 -370 -360 -350 -340 -320 -300 -280 -260 -240 -220 -200'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'

 return

  endif
*
* **************************** sea ice cover ( fract. ) **************************************
*
  if (cvar='code210' | cvar='SEAICE' | cvar='seaice')
    'set clevs 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 16 18 20 22 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'

 return

  endif
*
* *************** vertically integrated specific humidity kg/m**2 ) **************************
*
  if (cvar='code230' | cvar='QVI' | cvar='qvi')
    'set clevs 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 60'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
   'set ccols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48'

 return

 endif
*
* **************** vertically integrated liquid water content ( kg/m**2 ) ********************
 if (cvar='code231' | cvar='ALWCVI' | cvar='alwcvi')
   'set clevs 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44'
*    colnum=subwrd(result,5)
*    say 'Number of colors = 'colnum
   'set ccols 16 18 20 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'

 return

 endif

return

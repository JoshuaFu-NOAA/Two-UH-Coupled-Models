function palette(cvar)

  say 'cvar = 'cvar
*
* ********************* surface geopotential *************************************
*
if (cvar='code129' | cvar='GEOSP' | cvar='geosp')
   'set clevs -200 0 200 400 600 800 1000 1200 1400 1600 1800 2000 2500 3000 4000 5000'
   colnum=subwrd(result,5)
   say 'Number of colors = 'colnum
   'set ccols 32 35 38 40 44 48 52 56 60 68 72 76 80 87 96 19 23 27'
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
* *********************** surface temperature ( Celsius ) **************************
*
  if (cvar='ts' | cvar='code139'| cvar='code169' | cvar='TSURF' | cvar='TS' | cvar='tsurf')
    'set clevs -40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30 35'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 68 64 60 56 52 48 44 40 36 33 28 23 18 96 85 78 73'
    return
  endif
*
* *************************** soil wetness ( m ) *************************************
*
  if (cvar='code140' | cvar='WS' | cvar='ws')
    'set clevs 0 5 10 15 20 30 40 50 60 70'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 36 39 42 45 48 51 54 57 60'
    return
  endif
*
* ********************************* snow depth ( cm ) *********************************
*
  if (cvar='code141' | cvar='SN' | cvar='sn')
    'set clevs 0 1 2 2 4 6 8 10 50 100 999'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 36 39 42 45 48 51 54 57 60 20'
    return
  endif
*
* ******************************* surface roughness length ( m ) ****************************
* 
  if (cvar='code173' | cvar='AZ0' | cvar='az0')
    'set clevs 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 38 43 48 53 58 63 68 73 82 87 22 27'

 return

  endif
*
* *********************** surface background albedo ( fract. ) *******************************
*
  if (cvar='code174' | cvar='ALB' | cvar='alb')
    'set clevs 10 15 20 25 30 35 40 45 50 60 70 80 90'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 38 43 48 53 58 63 68 73 81 95 22 27'
    return
  endif
*
* ********************** east west orographic variance ( km**2 ) ****************************
*
  if (cvar='code90' | cvar='EWOVAR' | cvar='ewovar')
    'set clevs 0 0.01 0.02 0.05 0.1 0.2 0.5 1 2'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 28 23 18 96 80 74 64 54'
    return
  endif
*
* ********************** north south orographic variance ( km**2 ) **************************
*
  if (cvar='code91' | cvar='NSOVAR' | cvar='nsovar')
    'set clevs 0 0.01 0.02 0.05 0.1 0.2 0.5 1 2'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 28 23 18 96 80 74 64 54'
    return
  endif
*
* ****************** northwest southeast orographic variance ( km**2 ) ***********************
*
  if (cvar='code92' | cvar='NWSEOV' | cvar='nwseov')
    'set clevs 0 0.01 0.02 0.05 0.1 0.2 0.5 1 2'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 28 23 18 96 80 74 64 54'
    return
  endif
*
* ****************** northeast southwest orographic variance ( km**2 ) **********************
*
  if (cvar='code93' | cvar='NESWOV' | cvar='neswov')
    'set clevs 0 0.01 0.02 0.05 0.1 0.2 0.5 1 2'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 28 23 18 96 80 74 64 54'
    return
  endif
*
* ******************  vegetation ratio [%] **********************
*
  if (cvar='code198' | cvar='VGRAT' | cvar='vgrat')
    'set clevs 0 10 20 30 40 50 60 70 80 90'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 38 43 48 53 58 63 68 75 83'
    'set ccols 0 28 31 34 36 38 40 42 45 48 53'
    return
  endif
*
* ****************** orographic variance ( km**2 ) **********************
*
  if (cvar='code199' | cvar='VAROR' | cvar='varor')
    'set clevs 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 33 28 23 18 96 80 74 64 54'
    return
  endif
*
* ****************** leaf area index **********************
*
  if (cvar='code200' | cvar='VLT' | cvar='vlt')
    'set clevs 0.1 1 2 3 4 5 6 7 8 9'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 28 31 34 36 38 40 42 45 48 53'
    return
  endif
*
* ******************  vegetation type [%] **********************
*
  if (cvar='code212' | cvar='FOREST' | cvar='forest')
    'set clevs 0 10 20 30 40 50 60 70 80 90'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 28 31 34 36 38 40 42 45 48 53'
    return
  endif
*
* ****************** FAO data set (soil data flags)        [0...5] **********************
*
  if (cvar='code226' | cvar='FAO' | cvar='fao')
    'set clevs 0.5 1.5 2.5 3.5 4.5'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 0 24 32 39 48 85'
    return
  endif
*
* ****************** field capacity of soil [cm] **********************
*
  if (cvar='code229' | cvar='WSMX' | cvar='wsmx')
    'set clevs 0 5 10 15 20 30 40 50 60 70 80'
    'set clevs 0 5 10 15 20 30 40 50 60 70'
    colnum=subwrd(result,5)
    say 'Number of colors = 'colnum
    'set ccols 33 38 43 48 53 58 63 68 0'
    'set ccols 33 36 39 42 45 48 51 54 57 60 0'
    'set ccols 0 33 36 39 42 45 48 51 54 57 60'
    return
  endif

return

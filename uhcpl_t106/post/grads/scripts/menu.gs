* *************************************************************************
* ******* Script to display model data with grads in an easy way **********
* *************************************************************************
* *******                      written by                        **********
* *******                                                        **********
* *******   Hans-Stefan Bauer( stefan.bauer@dkrz.d400.de )       **********
* *******       Max Planck Institut fuer Meteorologie            **********
* *******                                                        **********
* *************************************************************************
*
* Variables:
*
*      _autopal        - Logical variable which controls the color resolution
*      _btn, btn       - Buttonnumber
*      _butlength      - Length of button
*      _buthight       - Hight of button
*       cint._fno      - Contour intervall
*       cmax._fno      - Highest contour level
*       cmin._fno      - Lowest contour level
*      _cvar._fno, cvar._fno     - Stores the selected variable
*      _ctlfile._fno   - The name of the selected .ctl file
*      _ctllist        - List of the .ctl files in the working directory
*      _ctlno          - Number of .ctl files in the working directory
*      _gxout, gxout   - Stores the selected graphics type
*      _latfactor      - Latitudinal shift factor
*       latgrid        - Y-position in GRID-coordinates
*      _latmax         - Maximum latitude to display
*      _latmin         - Minimum latitude to display
*       latpos         - Y-position in world coordinates
*      _levs._fno      - Stores the levels
*      _logz._fno      - Logical variable which controls the scaling of the vertical axis
*      _lohigh.0       - Startcoordinates ( control file )
*      _lolow.0        - Startcoordinates ( control file )
*      _lonfactor      - Longitudinal shift factor
*       longrid        - X-position in GRID-coordinates
*      _lonmax         - Maximum longitude to display
*      _lonmin         - Minimum longitude to display
*       lonpos         - X-position in world coordinates
*      _lthigh.0       - Startcoordinates ( control file )
*      _ltlow.0        - Startcoordinates ( control file )
*       mainbut        - Number of the pull down menu
*       max._fno       - Maximum value of the variable in the displayed region
*       min._fno       - Minimum value of the variable in the displayed region
*      _palette._fno   - Obtains the selected color palette
*      _plot._fno, plot- Obtains the selected plot type
*      _print          - color or greyscale print
*      _proj._fno      - Stores the actual projrction
*      _reg._fno, reg  - Region to plot
*      _res._fno, res  - Stores the selected resolution
*       temp           - Stores the actual level before animating through the atmosphere
*      _time._fno, time- Stores the selected timestep
*      _time1._fno     - Number of the selected timestep
*      _times._fno     - Stores the timesteps from the control file
*      _title._fno     - Title of the plot
*      _tnum._fno      - Number of timesteps
*       totnum._fno    - Total Number of gridpoints
*      _trace          - Variable that enables/disables status messages
*      _update         - Logical variable which controls the plot update
*      _vno, vno       - Stores the selected variable number
*      _vars._fno      - Stores the variables
*      _vnum._fno      - Number of variables
*       xinch          - X-position in Inch coordinates
*      _xsize._fno, xnum._fno    - Number of values in x-direction
*       yinch          - Y-position in Inch coordinates
*      _ysize._fno, ynum._fno    - Number of values in y-direction
*       ypos           - Y-position of the submenu
*      _znum._fno      - Number of levels
*      _z              - Stores the number of zooms
*                
* *************************************************************************
* **************************** Main Loop **********************************
* *************************************************************************
*
  _yflip='off'
  _null='on'
  _machine=''
  _name=''
  if (_machine != 'mac')
    '!finger $USER | head -1 | cut -f 3 -d":" > grads_user'
    _name=sublin(read(grads_user),2)
    '!rm -f grads_user'
    "!date '+%m/%d/%y' > date"
    _date=sublin(read(date),2)
    '!rm -f date'
  endif
  _read='no'
  _push='off'
  _levelupdate='on'
  _timeupdate='on'
  _palettenupdate='yes'
  _projupdate='yes'
  _gxoutupdate='yes'
  _codeupdate='yes'
  _logzupdate='yes'
  _plotupdate='yes'
  _trace='off'
  _resflag='off'
  _filecount=1
  _fno=1
  _z=0
  _start='yes'
  'set grads off'
  _butlength=1.2
  _butheight=0.4
  rc = init()
  rc = mainmenu()
  while (1)
    'q pos'
    if (_trace='on'); say result; endif;
    _btn = subwrd(result,7)
    if (_trace='on'); say 'BUTTON :'_btn; endif;
    if (_btn>=1 & _btn<=9)
      btn=_btn
      rc = choice(btn)
    endif
    if (_btn=65)
      say 'Type in your title !'
      pull _title._fno
      rc = display()
    endif
    if (_btn=63); _proj._fno='mollweide'; rc = latlon(); endif
    if (_btn=64); _proj._fno='ortho'; rc = latlon(); endif
    if (_btn=66); _null='off'; rc = display(); endif;
    if (_btn=67); _null='on'; rc = display(); endif;
    if (_btn=68); _fno=_fno+1;_filecount=_filecount+1; _read='yes'; rc = datinput(); endif;
    if (_btn=69); rc = init(); endif;
    if (_btn=70); rc = stat(); endif;
    if (_btn=71); _logz._fno='on'; rc = display(); endif;
    if (_btn=72); _logz._fno='off'; rc = display(); endif;
    if (_btn=73); rc = colprn(); endif;
    if (_btn=74); rc = greyprn(); endif;
    if (_btn=75); _update='yes'; endif;
    if (_btn=76); _update='no'; endif;
    if (_btn=77); _autopal='on'; rc = display(); endif;
    if (_btn=78); _autopal='off'; rc = display(); endif
    if (_btn=79); _proj._fno='latlon'; _z=0; rc = latlon(); endif;
    if (_btn=80); _proj._fno='nps'; rc = latlon(); endif;
    if (_btn=81); _proj._fno='sps'; rc = latlon(); endif;
    if (_btn=82); _proj._fno='robinson'; rc = latlon(); endif;
    if (_btn=83 | _btn=84); rc = animate(); endif;
    if (_btn=85); rc = positio(); endif;
    if (_btn=86); rc = fastleft(); endif;
    if (_btn=87); rc = fastrig(); endif;
    if (_btn=88); rc = pushup(); endif;
    if (_btn=89); rc = pushleft(); endif;
    if (_btn=90); rc = pushrig(); endif;
    if (_btn=91); rc = pushdown(); endif;
    if (_btn=92); rc = command(); endif;
    if (_btn=93); rc = reset(); endif;
    if (_btn=94); rc = unzoom();endif;
    if (_btn=95); rc = zoom(); endif;
    if (_btn=96); _read='no'; rc = init(); endif;
    if (_btn=97); rc = clear(); endif;
    if (_btn=98); rc = display(); endif;
    if (_btn=99); break; endif;
  endwhile (1) 
*
* *************************************************************************
* ******************************** FUNCTIONS ******************************
* *************************************************************************
*
* *************************** Reads a data file ***************************
*
 function datinput()

  if (_start='yes')
    'reinit'  
    _start='no'
  endif
  _autopal='false'
  if (_machine != 'mac')
    '!rm -f liste'
    '!rm -f ctlno'
    '!ls *.ctl >liste'
    '!wc -w <liste | sed -e s/" "/""/g > ctlno'
    _ctlno=sublin(read(ctlno),2)
    if (_trace='on'); say 'Number of Control Files :';_ctlno; endif;
    rc = ctl()
  else
    say 'Which dataset do you want to open ?'
    pull _ctlfile._fno
    'open '_ctlfile._fno
  endif
  rc = getinfo()

 return
*
* *********** Gets Info on the File ( Scans the control file ) ************
*
 function getinfo()

  if (_trace='on'); say 'Scanning control file ...'; endif; 
  'set dfile '_fno
  'q file '
  res=result
  dum = sublin(res,1)
  _title._fno = subwrd(dum,4)' 'subwrd(dum,5)' 'subwrd(dum,6)' 'subwrd(dum,7)' 'subwrd(dum,8)
  _title._fno = _title._fno' 'subwrd(dum,9)' 'subwrd(dum,10)' 'subwrd(dum,11)' 'subwrd(dum,12)
  dum = sublin(res,5)
  _xsize._fno = subwrd(dum,3)
  _ysize._fno = subwrd(dum,6)
  _znum._fno = subwrd(dum,9)
  _tnum._fno = subwrd(dum,12)
  dum = sublin(res,6)
  _vnum._fno = subwrd(dum,5)
  'q dims'
  dum = sublin(result,2)
  _lolow._z = subwrd(dum,6)
  _lohigh._z = subwrd(dum,8)
  dum = sublin(result,3)
  _ltlow._z = subwrd(dum,6)
  _lthigh._z = subwrd(dum,8)

*
*** Get time values * 
* 
  i = 1
  _times._fno = ''
  while (i<=_tnum._fno)
    'set t 'i
    'q time'
    _times._fno = _times._fno % ' ' % subwrd(result,3)
    i = i + 1
  endwhile
*
*** Get level values * 
*
  if (_znum._fno > 0)
    i=1
    _levs._fno=''
    while (i<=_znum._fno) 
      'set z 'i
      _levs._fno = _levs._fno % ' ' % subwrd(result,4)
      i = i + 1
    endwhile
    else
      _levs._fno='surface'
  endif
*
*** Get variables *
*
  i = 0
  _vars._fno = ''
  while (i<_vnum._fno)
    dum = sublin(res,i+7)
    _vars._fno = _vars._fno % ' ' % subwrd(dum,1)
    if (subwrd(dum,1)='u'); _vec='u';endif
    if (subwrd(dum,1)='v' & _vec='u'); _vec='uv'; endif
    if (subwrd(dum,1)='u10'); _vec='u10'; endif
    if (subwrd(dum,1)='v10' & _vec='u10'); _vec='u10v10'; endif
    _vno = i + 1
    _vname._vno._fno = subwrd(dum,4)' 'subwrd(dum,5)' 'subwrd(dum,6)' 'subwrd(dum,7)' 'subwrd(dum,8)' 'subwrd(dum,9)' 'subwrd(dum,10)
    i = i + 1
  endwhile
  if (_trace='on')
    say ' FILENUMBER : '_fno
    say ' TITLE      : '_title._fno
    say ' NUMBER OF VALUES IN X-DIRECTION : '_xsize._fno
    say ' NUMBER OF VALUES IN Y-DIRECTION : '_ysize._fno
    say ' NUMBER OF LEVELS                : '_znum._fno
    say ' NUMBER OF TIMESTEPS             : '_tnum._fno
    say ' NUMBER OF VARIABLES             : '_vnum._fno
    say ' STARTCOORDINATES :'
    say '            XLOW  : '_lolow._z
    say '            XHIGH : '_lohigh._z
    say '            YLOW  : '_ltlow._z
    say '            YHIGH : '_lthigh._z
    say ' TIMES      : '_times._fno
    say ' LEVS       : '_levs._fno    
    say ' VARIABLES  : '_vars._fno
    say ' VARIABLENBESCHRIFTUNG           : '_vname._vno._fno
  endif;
  rc = startval()

 return
*
* ************************** Sets Startvalues *****************************
*
 function startval()

  if (_trace='on'); say 'Set start values ...'; endif;
  _file._fno=_ctlfile._fno
  if (_read='no')
    _z=0
    _zval=0
  endif 
  if (subwrd(_levs_fno,1) < 1)
    _zunit._fno='km'
  endif
  if (subwrd(_levs._fno,1) >= 1 & subwrd(_levs._fno,1) < 30)
    _zunit._fno='Modellevel'
  endif
  if (subwrd(_levs._fno,1) > 30)
    _zunit._fno='hPa'
  endif
  if (_trace='on'); say ' ZUNIT        : '_zunit._fno; endif
* Startvalues for the first dataset
  if (_read='no')
    _logz._fno='off'
    _autopal='off'
    _update='no'
    _proj._fno='latlon'
    _vno=1
    _reg._fno='world'
    _res._fno='lowres'
    _gxout='contour'
    _palette._fno='rainbow'
    _time._fno=subwrd(_times._fno,1)
    _time1._fno=1
    _clev._fno=subwrd(_levs._fno,1)
    _cvar._fno=subwrd(_vars._fno,1)
    _levno=1
    _templonmin=_lolow._z
    _templonmax=_lohigh._z
    _templatmin=_ltlow._z
    _templatmax=_lthigh._z
  else
*   If you read a new dataset, use the same startvalues
    _fold=_fno-1
    _logz._fno=_logz._fold
    _proj._fno=_proj._fold
    _reg._fno=_reg._fold
    _res._fno=_res._fold
    _palette._fno=_palette._fold
    _time._fno=_time._fold
    _time1._fno=_time1._fold
    _cvar._fno=_cvar._fold
    if (_znum._fno != 0)
      _highlev.0=subwrd(_levs._fold,_znum._fold)
    endif
    _clev._fno=_clev._fold
    _templonmin=_lolow._z
    _templonmax=_lohigh._z
    _templatmin=_ltlow._z
    _templatmax=_lthigh._z
  endif
*
  if (_xsize._fno > 1 & _ysize._fno > 1)
    say 'I choose normal XY-Plot !!'
    _plot._fno='latlon'
  endif
  if (_xsize._fno = 1 & _ysize._fno > 1 & _znum._fno > 1)
    say 'I choose meridional crossection !!'
    _plot._fno='meridional'
    'set lon '_lolow._z
    _lon._fno=_lolow._z
    'set lev 'subwrd(_levs._fno,1)' 'subwrd(_levs._fno,_znum._fno)
    'set t '1
  endif
  if (_xsize._fno > 1 & _ysize._fno = 1 & _znum._fno > 1)
    say 'I choose zonal crossection !!'
    _plot._fno='zonal'
    'set lat '_ltlow._z
    _lat._fno=_ltlow._z
    'set lev 'subwrd(_levs._fno,1)' 'subwrd(_levs._fno,_znum._fno)
    'set t '1
  endif
  if (_xsize._fno = 1 & _ysize._fno = 1 & _znum._fno > 1 & tnum._fno > 1)
    say 'I choose Hovmoeller !!'
    _plot._fno='hovm_t_z'
    'set lon '_lolow._z
    'set lat '_ltlow._z
    'set lev 'subwrd(_levs._fno,1)' 'subwrd(_levs._fno,_znum._fno)
    'set t '1' '_tnum._fno
  endif
  rc = mainmenu()

 return
*
* **************** Controls the choice in a submenu *************************
*
 function choice(btn)

  while (1)
    if (btn=1); rc = resmen(); break; endif
    if (btn=2); rc = plotmen(); break; endif
    if (btn=3); rc = timemen(); break; endif
    if (btn=4); rc = codmen(); break; endif
    if (btn=5); rc = levmen(); break; endif
    if (btn=6); rc = gxmen(); break; endif
    if (btn=7); rc = palette(); break; endif
    if (btn=8); rc = region(); break; endif
    if (btn=9); rc = filemenu(); break; endif
    if (btn>9); break; endif
  endwhile (1)

 return
*
* *************************** Clears the plot area ***************************
*
 function clear()

  if (_trace='on'); say 'Clear screen ...'; endif;
  rc = mainmenu() 

 return
*
* *************** Resets Startvalues and clears the plot area ****************
*
 function reset()

  if (_trace='on'); say 'Reset data ...'; endif;
  rc = startval()
  rc = display()

 return 
*
* ***** Resets all features, closes all datafiles, reads a new file *********
*
 function init()

  _start='yes'
  _filecount=1
  _fno=1
  rc = datinput()

 return
* ********************************* Info ************************************
*
 function info()

  'q gxinfo'
  dim=sublin(result,3)
  _xlow._fno=subwrd(dim,4)
  _xhigh._fno=subwrd(dim,6)
  dim=sublin(result,4)
  _ylow._fno=subwrd(dim,4)
  _yhigh._fno=subwrd(dim,6)

 return
*
* ************************ Gets numbers of levels ***************************
*
 function testlev()

 'q file'
 line = sublin(result,_vno+6)
* clev.tmp = _clev._fno 
 if (subwrd(line,2)=0)
*   znum.tmp = _znum._fno
   _znum._fno = 0
   _clev._fno = 'surface'
 endif
 if (subwrd(line,2)=19)
   _znum._fno=19
   _clev._fno=subwrd(_levs._fno,_levno)
 endif
 if (subwrd(line,2)=39)
   _znum._fno=39
   _clev._fno=subwrd(_levs._fno,_levno)
 endif

 return
*
* ****************** Switches back to latlon projection *********************
*
 function latlon()

  if (_trace='on'); say 'Projection = '_proj._fno; endif
  if (_plot._fno='latlon')
    _zval=_z
  else
   _z=_zval
  endif
  if (_xsize._fno>1 & _ysize._fno>1)
    'set lon '_lolow._z' '_lohigh._z
    'set lat '_ltlow._z' '_lthigh._z
    'set lev '_clev._fno
    _plot._fno='latlon'
  endif
  if (_xsize._fno=1 & _ysize._fno>1)
    'set lon '_lolow._z
    'set lat '_ltlow._z' '_lthigh._z
    say 'It is not possible to return to latlon, because the X-coordinate is fixed.'
    say 'I return to meridional crossection !'
    _plot._fno='meridional'
    if (_znum._fno>1)
      'set lev 'subwrd(_levs._fno,1)' 'subwrd(_levs._fno,_znum._fno)
    else
      'set lev '_clev._fno
    endif 
  endif
  if (_xsize._fno>1 & _ysize._fno=1)
    say '_zval = '_zval
    say 'Lohigh = '_lohigh._z
    'set lon '_lolow._z' '_lohigh._z
    'set lat '_ltlow._z
    say 'It is not possible to return to latlon, because Y-coordinate is fixed.'
    say 'Therefore I return to zonal crossection !!'
    _plot._fno='zonal'
    if (_znum._fno>1)
      'set lev 'subwrd(_levs._fno,1)' 'subwrd(_levs._fno,_znum._fno)
    else
      'set lev '_clev._fno
    endif
  endif 
  'set t '_time1._fno
  if (_proj._fno='nps' | _proj._fno='sps' | _proj._fno='robinson' | _proj._fno='mollweide')
    _z=0
  endif
  rc = display()

 return
*
* ******************* Draws the lines and chooses lat/lon *******************
*
 function choose()

 rc = info()
 'set line 1 1 6'
 while (1)
   'q pos'
   if (_trace='on'); say result; endif;
   _button=subwrd(result,5)
   if (_button=1)
     x1=subwrd(result,3)
     y1=subwrd(result,4)
     if (_xsize._fno>1 & _ysize._fno>1)
       if (_plot._fno != 'meridional' & _plot._fno != 'hovm_t_lat')
         'draw line '_xlow._fno' 'y1' '_xhigh._fno' 'y1
       endif
       if (_plot._fno != 'zonal' & _plot._fno != 'hovm_t_lon')
         'draw line 'x1' '_ylow._fno' 'x1' '_yhigh._fno
       endif
     else
       if (_plot._fno='time')
         'draw line 'x1' '_ylow._fno' 'x1' '_yhigh._fno
         'draw line '_xlow._fno' 'y1' '_xhigh._fno' 'y1
       else
         'draw line 'x1' '_ylow._fno' 'x1' '_yhigh._fno
       endif 
     endif
     'q xy2w 'x1' 'y1
     if (_xsize._fno>1 & _ysize._fno>1)
       _lon._fno=subwrd(result,3)
       _lat._fno=subwrd(result,6)
     endif
     if (_xsize._fno=1 & _ysize._fno>1)
       _lat._fno=subwrd(result,3)
       _lev._fno=subwrd(result,6)
     endif
     if (_xsize._fno>1 & _ysize._fno=1)
       _lon._fno=subwrd(result,3)
       _lev._fno=subwrd(result,6)
     endif
     if (_machine='mac')
       say 'Press c and then return to confirm the last selection !'
       say 'Press any other key and then return for a new selection !'
       pull variable
       if (variable='c')
         break
       endif
     endif
   else
     break
   endif
 endwhile (1)

 return
*
* ************************** Draws the main menu ****************************
* 
    function mainmenu()

    if (_trace='on'); say 'Set up main menu ...'; endif;
    if (_z=0); _reg._fno='world'; endif
    'set rgb 90 100 100 100'
    'set rgb 91  50 50 50 '
    'set rgb 92 200 200 200'
    'set button 1 90 91 92 0 90 91 92 6'
    'clear'
    'set string 1 c 6'
    'draw button 1 0.7 8.3 '_butlength' '_butheight' '_res._fno
    'draw button 2 1.9 8.3 '_butlength' '_butheight' '_plot._fno
    'draw button 3 3.1 8.3 '_butlength' '_butheight' '_time._fno
    'draw button 4 4.3 8.3 '_butlength' '_butheight' '_cvar._fno
    'draw button 5 5.5 8.3 '_butlength' '_butheight' '_clev._fno
    'draw button 6 6.7 8.3 '_butlength' '_butheight' '_gxout
    'draw button 7 7.9 8.3 '_butlength' '_butheight' '_palette._fno
    'draw button 8 9.1 8.3 '_butlength' '_butheight' '_reg._fno
    'draw button 9 10.3 8.3 '_butlength' '_butheight' '_file._fno
    'draw button 63 0.5 4.75 '%_butlength-0.3' '%_butheight-0.1' mollweide'
    'draw button 64 0.5 4.45 '%_butlength-0.3' '%_butheight-0.1' ortho'
    'draw button 65 0.5 1.65 '%_butlength-0.3' '%_butheight-0.1' Title'
    if (_null='on')
      'draw button 66 0.5 3.15 '%_butlength-0.3' '%_butheight-0.1' Zero off'
    else
      'draw button 67 0.5 3.15 '%_butlength-0.3' '%_butheight-0.1' Zero on'
    endif
    'draw button 68 0.7 0.9 '_butlength' '%_butheight-0.1' READ DATA'
    'draw button 69 5.0 0.2 '_butlength-0.2' '_butheight' REINIT'
    'draw button 70 0.5 1.95 '%_butlength-0.3' '%_butheight-0.1' statistics'
    if (_logz._fno='off')
      'draw button 71 0.5 2.25 '%_butlength-0.3' '%_butheight-0.1' logz on'
    else
      'draw button 72 0.5 2.25 '%_butlength-0.3' '%_butheight-0.1' logz off'
    endif
    'draw button 73 0.5 3.9 '%_butlength-0.3' '%_butheight-0.1' color print'
    'draw button 74 0.5 3.6 '%_butlength-0.3' '%_butheight-0.1' greyscale'
    if (_update='no')
      'draw button 75 0.5 2.85 '%_butlength-0.3' '%_butheight-0.1' update on'
    else  
      'draw button 76 0.5 2.85 '%_butlength-0.3' '%_butheight-0.1' update off'
    endif     
    if (_autopal='off')
      'draw button 77 0.5 2.55 '%_butlength-0.3' '%_butheight-0.1' autopal on'
    else
      'draw button 78 0.5 2.55 '%_butlength-0.3' '%_butheight-0.1' autopal off'
    endif
    'draw button 79 0.5 5.95 '%_butlength-0.3' '%_butheight-0.1' lat / lon' 
    'draw button 80 0.5 5.65 '%_butlength-0.3' '%_butheight-0.1' nps'
    'draw button 81 0.5 5.35 '%_butlength-0.3' '%_butheight-0.1' sps'
    'draw button 82 0.5 5.05 '%_butlength-0.3' '%_butheight-0.1' robinson'
    'draw button 83 0.5 6.45 '%_butlength-0.2' '%_butheight-0.1' animate z'
    'draw button 84 0.5 6.75 '%_butlength-0.2' '%_butheight-0.1' animate t'
    'draw button 85 2.0 0.2 '_butlength-0.2' '%_butheight' POSITION'
    'draw button 86 9.95 0.9 0.25 0.25 <<'
    'draw button 87 10.8 0.9 0.25 0.25 >>'
    'draw button 88 10.4 1.2 0.25 0.25 ^'
    'draw button 89 10.2 0.9 0.25 0.25 <'
    'draw button 90 10.55 0.9 0.25 0.25 >'
    'draw button 91 10.4 0.6 0.25 0.25 v'
    'draw button 92 3.0 0.2 '_butlength-0.2' '_butheight' COMMAND'
    'draw button 93 4.0 0.2 '_butlength-0.2' '_butheight' RESET'
    'draw button 94 6.0 0.2 '_butlength-0.2' '_butheight' UNZOOM'
    'draw button 95 7.0 0.2 '_butlength-0.2' '_butheight' ZOOM'
    'draw button 96 8.0 0.2 '_butlength-0.2' '_butheight' NEW  DATA'
    'draw button 97 9.0 0.2 '_butlength-0.2' '_butheight' CLEAR'
    'draw button 98 10.3 0.2 '_butlength' '_butheight' DISPLAY'
    'draw button 99 0.7 0.3 '%_butlength+0.1' '%_butheight+0.1' QUIT'

    return
*
* ********************** Reads a GrADS Command ********************************
*
 function command()

  say 'Type in one or more commands !'
  say 'Between each of them press the "enter" Button."'
  say 'If you want to stop, press "q" ant then "enter"!'
  while (1)
    say 'Command: '
    pull answer
    if (answer = "d" | answer = "display"); rc = display(); flag='on'; endif
    if (answer = "reset"); rc = reset(); flag='on'; endif
    if (answer = "clear" | answer = "c"); rc = clear(); flag='on'; endif
    if (answer = "q"); break; endif
    if (flag != 'on')
      comout="!echo "answer" > command"
      comout
      'exec command'
    endif
  endwhile
*  '!rm -f command'

 return
*
* ******** Draws the ctl-file menu and selects the desired file ***********************
*
 function ctl()

  if (_trace='on'); say 'Ctl menu'; endif;
  mainbut = 1
  'clear'
  'set string 1 c 6'
  'set strsiz 0.2 0.25'
  'draw string 5.5 8.0 Which dataset do you want ?'
  'set rgb 91 50 50 50'
  'set rgb 92 200 200 200'
  'set button 1 90 91 92 0 90 91 92 6'
  if (_ctlno <= 17); x=5.0; endif;
  if (_ctlno > 17 & _ctlno <= 34); x=3.5; endif;
  if (_ctlno > 34 & _ctlno <= 51); x=2.5; endif;
  if (_ctlno > 48); x=1.5; endif;
  x1=x+(_butlength+1.45)
  x2=x+2*(_butlength+1.45)
  x3=x+3*(_butlength+1.45)   
  y=7.5
  i=1
  while(i<=(_ctlno))
  _ctllist=sublin(read(liste),2)
    buttonnr=20+i
    if (i<=17)
      ypos=y-(i*_butheight)
      'draw button 'buttonnr' 'x' 'ypos' '%_butlength+1.2' '_butheight' '_ctllist''
    endif
    if (i>17 & i<=34)
      ypos=y-((i-17)*_butheight)
      'draw button 'buttonnr' 'x1' 'ypos' '%_butlength+1.2' '_butheight' '_ctllist''
    endif
    if (i>34 & i <=51)
      ypos=y-((i-34)*_butheight)
      'draw button 'buttonnr' 'x2' 'ypos' '%_butlength+1.2' '_butheight' '_ctllist''
    endif
    if (i>51 & i <=68)
       ypos=y-((i-51)*_butheight)
      'draw button 'buttonnr' 'x3' 'ypos' '%_butlength+1.2' '_butheight' '_ctllist''
    endif 
    if (i>68)
      say 'To many .ctl files in the directory !!!!'
      say 'I only use the first 68 !!!'
      break
    endif   
    i=i+1
  endwhile
  'q pos'
  if (_trace='on'); say result; endif;
  btn = subwrd(result,7)
  If (_trace='on'); say 'Button = 'btn; endif
  close(liste)
  if (btn=''); rc = ctl(); endif
  if (btn>20 & btn<=buttonnr)
    i=1
    while(i<=(btn-20))
     _ctlfile._fno=sublin(read(liste),2)
     i=i+1
    endwhile
    _ctlno=(btn-20)
    'open '_ctlfile._fno
    if (_trace='on'); say ' CTLFILE : '_ctlfile._fno; endif;
  endif
  '!rm -f liste'
  '!rm -f ctlno'
   close(liste)
   close(ctlno)

 return
*
* ************************* Prints out some Statistics **************************
*
 function stat()

  if (_trace='on'); say 'Statistic Function ...'; endif;
  'set gxout stat'
  'd '_cvar._fno
  sizes=sublin(result,5)
  minmax=sublin(result,8)
  cminmax=sublin(result,9)
  xnum._fno=subwrd(sizes,3)
  ynum._fno=subwrd(sizes,4)
  totnum._fno=subwrd(sizes,5)
  min._fno=subwrd(minmax,4)
  max._fno=subwrd(minmax,5)
  cmin._fno=subwrd(cminmax,5)
  cmax._fno=subwrd(cminmax,6)
  cint._fno=subwrd(cminmax,7)
  say 'Number of x-values:     'xnum._fno
  say 'Number of y-values:     'ynum._fno
  say 'Total number of values: 'totnum._fno
  say 'Minimum value of variable at this level: 'min._fno
  say 'Maximum value of variable at this level: 'max._fno
  say 'Contour minimum:        'cmin._fno
  say 'Contour maximum:        'cmax._fno
  say 'Contour intervall:      'cint._fno
  'set gxout '_gxout

 return
*
* ******* Draws the resolution menu and selects the desired resolution **********
*
 function resmen()

  if (_trace='on'); say 'Resolution menu'; endif;
  mainbut = 1
  i = 3
  x = 0.7
  y = 8.3
  'draw button 21 0.7 7.9 '_butlength' '_butheight' lowres'
  'draw button 22 0.7 7.5 '_butlength' '_butheight' mres'
  'draw button 23 0.7 7.1 '_butlength' '_butheight' hires'
  'q pos'
  if (_trace='on'); say result; endif;
  btn = subwrd(result,7)
  res = _res._fno
  if (btn=-999); rc = delmen(x,y,i,mainbut,res); endif
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,res)
    rc = choice(btn)
  endif
  _res._fno = res
  if (btn>=94 & btn<=99); rc = mainmenu(); endif
  if (btn=21); _res._fno='lowres'; endif
  if (btn=22); _res._fno='mres'; endif
  if (btn=23); _res._fno='hires'; endif
  if (_trace='on'); say _res._fno; endif;
  _resflag='on'
  if (_update='yes')
    rc = display()
    _resflag='off'
  else
    rc = mainmenu()
  endif;
      
 return
*
* ********* Draws the plotmenu and selects the desired kind of plot *************
*
 function plotmen()

  if (_trace='on'); say 'Plot menu'; endif;
  mainbut = 2
  i = 5
  x = 1.9
  y = 8.3
  'draw button 21 1.9 7.9 '_butlength' '_butheight' latlon'
  'draw button 22 1.9 7.5 '_butlength' '_butheight' meridional'
  'draw button 23 1.9 7.1 '_butlength' '_butheight' zonal'
  'draw button 24 1.9 6.7 '_butlength' '_butheight' time series'
  'draw button 25 1.9 6.3 '_butlength' '_butheight' hovm_t_z'
  'draw button 26 1.9 5.9 '_butlength' '_butheight' hovm_t_lat'
  'draw button 27 1.9 5.5 '_butlength' '_butheight' hovm_t_lon'
  'draw button 28 1.9 5.1 '_butlength' '_butheight' radiosonde'
  'q pos'
  if (_trace='on'); say result; endif;
  btn = subwrd(result,7)
  plot = _plot._fno
  if (btn=-999); rc = delmen(x,y,i,mainbut,plot); endif
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,plot)
    rc = choice(btn)
  endif
  _plot._fno = plot
  if (btn>=94 & btn<=99); rc = mainmenu(); endif
  if (btn=21)
    _proj._fno='latlon'
    rc = latlon()
    if (_z=0); 'set lat -90 90';endif
  endif
  if (btn=22)    
    if (_ysize._fno>1 & _znum._fno>1)
      if (_xsize._fno=1)
        _lon._fno=_lolow._z
      else
        rc = latlon()
        _plot._fno='meridional'
        say 'Left Mouse Button   : Select the longitude'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Draws the cross section with the last longitude'
        endif
        rc = choose()
      endif
      if (_trace='on'); say 'Longitude: '_lon._fno; endif
      'set lon '_lon._fno
      'set lat '_ltlow._z' '_lthigh._z
      _lolev._z=subwrd(_levs._fno,1)
      _highlev._z=subwrd(_levs._fno,_znum._fno)
      'set lev '_lolev._z' '_highlev._z
      'set t '_time1._fno
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (btn=23)
    if (_xsize._fno>1 & _znum._fno>1)
      if(_ysize._fno=1)
        _lat._fno=_ltlow._z
      else
        rc = latlon()
        _plot._fno='zonal'
        say 'Left Mouse Button   : Select the latitude'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Draws the cross section with the last latitude'
        endif
        rc = choose()
      endif
      if (_trace='on'); say 'Latitude: '_lat._fno; endif
      'set lon '_lolow._z' '_lohigh._z
      'set lat '_lat._fno
      _lolev._z=subwrd(_levs._fno,1)
      _highlev._z=subwrd(_levs._fno,_znum._fno)
      'set lev '_lolev._z' '_highlev._z
      'set t '_time1._fno
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (btn=24)
    if (_tnum._fno>1)
      if (_xsize._fno=1 & _ysize._fno=1)
        _lon._fno=_lolow._z
        _lat._fno=_ltlow._z
      else
        rc = latlon()
        _plot._fno='time'
        if (_xsize._fno=1 & _ysize._fno>1)
          _lon._fno=_lolow._z
          say 'Left Mouse Button   : Select the latitude and level you want'
          if (_machine != 'mac')
            say 'Middle Mouse Button : Confirms the last selection'
          endif
        endif
        if (_xsize._fno>1 & _ysize._fno=1)
          _lat._fno=_ltlow._z
          say 'Left Mouse Button   : Select the longitude and level you want'
          if (_machine != 'mac')
            say 'Middle Mouse Button : Confirms the last selection'
          endif
        endif
        if (_xsize._fno>1 & _ysize._fno>1)
          _lev._fno=_clev._fno
          say 'Left Mouse Button   : Select the location you want'
          if (_machine != 'mac')
            say 'Middle Mouse Button : Confirms the last selection'
          endif
        endif
      rc = choose()
      endif
      'set lon '_lon._fno
      'set lat '_lat._fno
      'set lev '_lev._fno
      'set t '1' '_tnum._fno
      _gxout='line'
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (btn=25)
    if (_tnum._fno>1 & _znum._fno>1)
      rc = latlon()
      _plot._fno='hovm_t_z'
      if (_xsize._fno=1 & _ysize._fno=1)
        _lon._fno=_lolow._z
        _lat._fno=_ltlow._z
      endif
      if (_xsize._fno=1 & _ysize._fno>1)
        _lon._fno=_lolow._z
        say 'Left Mouse Button   : Select the latitude you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last latitude'
        endif
        rc = choose()
      endif
      if (_xsize._fno>1 & _ysize._fno=1)
        _lat._fno=_ltlow._z
        say 'Left Mouse Button   : Select the longitude you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last longitude'
        endif
        rc = choose()
      endif
      if (_xsize._fno>1 & _ysize._fno>1)
        say 'Left Mouse Button   : Select the location you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last location'
        endif
        rc = choose()
      endif
      'set lon '_lon._fno
      'set lat '_lat._fno
      'set t '1' '_tnum._fno
      _lolev._z=subwrd(_levs._fno,1)
      _highlev._z=subwrd(_levs._fno,_znum._fno)
      'set lev '_lolev._z' '_highlev._z
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (btn=26)
    if (_tnum._fno>1 & _ysize._fno>1)
      if (_xsize._fno=1)
        _lon._fno=_lolow._z
        _lev._fno=_clev._fno
      else
        rc = latlon()
        _plot._fno='hovm_t_lat'
        say 'Left Mouse Button   : Select the longitude you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last longitude you select'
        endif
        rc = choose()
        _lev._fno=_clev._fno        
      endif 
      'set lon '_lon._fno
      'set lat '_ltlow._z' '_lthigh._z 
      'set t '1' '_tnum._fno
      'set lev '_lev._fno
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (btn=27)
    if (_tnum._fno>1 & _xsize._fno>1)
      _plot._fno='hovm_t_lon'
      if (_ysize._fno=1)
        _lat._fno=_ltlow._z
        _lev._fno=_clev._fno
      else
        rc = latlon()
        _plot._fno='hovm_t_lon'
        say 'Left Mouse Button   : Select the latitude you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last latitude you select'
        endif
        rc = choose()
        _lev._fno=_clev._fno           
      endif 
      'set lon '_lolow._z' '_lohigh._z
      'set lat '_lat._fno 
      'set t '1' '_tnum._fno
      'set lev '_lev._fno
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (btn=28)
    if (_znum._fno>1)
      rc = latlon()
      _plot._fno='radiosonde'
      if (_xsize._fno=1 & _ysize._fno>1)
        _lon._fno=_lolow._z
        say 'Left Mouse Button   : Select the latitude you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last latitude'
        endif
        rc = choose()
      endif
      if (_xsize._fno>1 & _ysize._fno=1)
        _lat._fno=_ltlow._z
        say 'Left Mouse Button   : Select the longitude you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last longitude'
        endif
        rc = choose()
      endif
      if (_xsize._fno>1 & _ysize._fno>1)
        say 'Left Mouse Button   : Select the location you want'
        if (_machine != 'mac')
          say 'Middle Mouse Button : Confirms the last location'
        endif
        rc = choose()
      endif
      'set lon '_lon._fno
      'set lat '_lat._fno
      'set lev 'subwrd(_levs._fno,1)' 'subwrd(_levs._fno,_znum._fno)
    else
      say 'NOT POSSIBLE !! Dimension environment not complete.'
    endif
  endif
  if (_plotupdate='yes')
    i=1
    while (i <= _filecount)
      _clev.i=_clev._fno
      _proj.i=_proj._fno
      _plot.i=_plot._fno
      if (btn=22)
         _lon.i=_lon._fno
      endif
      if (btn=23)
        _lat.i=_lat._fno
      endif
      i=i+1
    endwhile
  endif
  if (btn!=21)
    if (_update='yes')
      rc = display()
    else
      rc = mainmenu()
    endif;
  endif
      
 return
*
* ****** Draws the gxout menu and selects the desired Graphics output type ******
*
 function gxmen()
 
  if (_trace='on'); say 'Gxout menu'; endif
  mainbut = 6
  i = 4
  x = 6.7
  y = 8.3
  'draw button 21 6.7 7.9 '_butlength' '_butheight' contour'
  'draw button 22 6.7 7.5 '_butlength' '_butheight' shaded'
  'draw button 23 6.7 7.1 '_butlength' '_butheight' grfill'
  'draw button 24 6.7 6.7 '_butlength' '_butheight' grid'
  'draw button 25 6.7 6.3 '_butlength' '_butheight' vector'
  'draw button 26 6.7 5.9 '_butlength' '_butheight' stream'
  'q pos'
  if (_trace='on'); say result; endif;
  btn = subwrd(result,7)
  gxout=_gxout
  if (btn=-999); rc = delmen(x,y,i,mainbut,gxout); endif
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,gxout)
    rc = choice(btn)
  endif
  _gxout=gxout
  if (btn>=94 & btn<=99); rc = mainmenu(); endif
  if (btn=21); _gxout='contour'; endif
  if (btn=22); _gxout='shaded'; endif
  if (btn=23); _gxout='grfill'; endif
  if (btn=24); _gxout='grid'; endif
  if (btn=25); _gxout='vector'; endif
  if (btn=26); _gxout='stream'; endif
  if (_trace='on'); say 'gxout = '_gxout; endif;
  if (_gxout='vector' | _gxout='stream')
    rc = delmen(4.3,8.3,1,3,'u ; v')
  endif
  'set gxout '_gxout
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif;

 return
*
* ******** Draws the code menu and selects the desired code ***********************
*
 function codmen()

  if (_trace='on'); say 'Code menu'; endif;
  mainbut = 4
  x=4.3
  x1=x+(_butlength)
  x2=x1+(_butlength)
  x3=x2+(_butlength)   
  y=8.3
  i=1
  while(i<=_vnum._fno)
    buttonnr=20+i
    if (i<18)
      ypos=y-(i*_butheight)
      'draw button 'buttonnr' 'x' 'ypos' '_butlength' '_butheight' 'subwrd(_vars._fno,i)
    endif
    if (i>17 & i<34)
      ypos=y-((i-16)*_butheight)
      'draw button 'buttonnr' 'x1' 'ypos' '_butlength' '_butheight' 'subwrd(_vars._fno,i)
    endif
    if (i>33 & i<50)
      ypos=y-((i-32)*_butheight)
      'draw button 'buttonnr' 'x2' 'ypos' '_butlength' '_butheight' 'subwrd(_vars._fno,i)
    endif
    if (i>49 & i<66)
      ypos=y-((i-48)*_butheight)
      'draw button 'buttonnr' 'x3' 'ypos' '_butlength' '_butheight' 'subwrd(_vars._fno,i)
    endif
    i=i+1
  endwhile
  'q pos'
  if (_trace='on'); say result; endif;
  cvar._fno=_cvar._fno
  btn = subwrd(result,7)
  if (btn=-999); rc = delmen(x,y,i,mainbut,cvar._fno); endif;
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,cvar._fno)
    rc = choice(btn)
  endif
  _cvar._fno=cvar._fno
  if (btn>=94 & btn<=99); rc = mainmenu(); endif;
  if (btn>20 & btn<=buttonnr)
    _cvar._fno=subwrd(_vars._fno,btn-20)
    _vno=(btn-20)
  endif
  if (_trace='on'); say 'cvar._fno = '_cvar._fno; endif;
  if (_codeupdate='on')
    i=1
    while (i <= _filecount)
      _cvar.i=_cvar._fno
      i=i+1
    endwhile
  endif
*
  rc = testlev()
*
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif;

return
*
* ************** Draws the level menu and selects the desired level ***************
*
 function levmen()

  if (_trace='on'); say 'Level menu ...'; endif;
  mainbut=5
  x=5.5
  x1=x+(_butlength)
  x2=x+2*(_butlength)
  x3=x2+_butlength
  y=8.3
  i=1
  while(i<=_znum._fno)
    buttonnr=20+i
    if (i<18)
      ypos=y-(i*_butheight)
      'draw button 'buttonnr' 'x' 'ypos' '_butlength' '_butheight' 'subwrd(_levs._fno,i)
    endif
    if (i>17 & i<34)
      ypos=y-((i-16)*_butheight)
      'draw button 'buttonnr' 'x1' 'ypos' '_butlength' '_butheight' 'subwrd(_levs._fno,i)
    endif
    if (i>33 & i<50)
      ypos=y-((i-32)*_butheight)
      'draw button 'buttonnr' 'x2' 'ypos' '_butlength' '_butheight' 'subwrd(_levs._fno,i)
    endif
    if (i>49 & i<66)
      ypos=y-((i-48)*_butheight)
      'draw button 'buttonnr' 'x3' 'ypos' '_butlength' '_butheight' 'subwrd(_levs._fno,i)
    endif
    i=i+1
  endwhile
  'q pos'
  if (_trace='on'); say result; endif;
  btn=subwrd(result,7)
  clev._fno=_clev._fno
  if (btn=-999); rc = delmen(x,y,i,mainbut,clev._fno); endif
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,clev._fno)
    rc = choice(btn)
  endif
  _clev._fno=clev._fno
  if (btn>=94 & btn<=99); rc = mainmenu(); endif
  if (btn>20 & btn<=buttonnr)
    _clev._fno=subwrd(_levs._fno,btn-20)
    'set lev '_clev._fno
    _lev._fno=_clev._fno
    _levno = btn-20
     if (_trace='on'); say 'Level = '_clev._fno; endif;
  endif
  if (_levelupdate='on')
    i=1
    while (i <= _filecount)
      _clev.i=_clev._fno
      i=i+1
    endwhile
  endif
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif; 

 return
*
* ***************** Draws the time menu and selects the desired time ***********
* 
 function timemen()

  if (_trace='on'); say 'Time menu'; endif;
  mainbut=3
  x=3.1
  x1=x+(_butlength)
  x2=x+2*(_butlength)
  x3=x+3*(_butlength)
  x4=x+4*(_butlength)
  x5=x+5*(_butlength)
  x6=x+6*(_butlength)
  y=8.3
  i = 1
  while (i<=_tnum._fno)
    buttonnr=20+i
    if (i<18)
      ypos=y-(i*_butheight)
      'draw button 'buttonnr' 'x' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    if (i>17 & i<34)
      ypos=y-((i-16)*_butheight)
      'draw button 'buttonnr' 'x1' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    if (i>33 & i<50)
      ypos=y-((i-32)*_butheight)
      'draw button 'buttonnr' 'x2' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    if (i>49 & i<66)
      ypos=y-((i-48)*_butheight)
      'draw button 'buttonnr' 'x3' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    if (i>65 & i<82)
      ypos=y-((i-64)*_butheight)
      'draw button 'buttonnr' 'x4' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    if (i>81 & i<98)
      ypos=y-((i-80)*_butheight)
      'draw button 'buttonnr' 'x5' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    if (i>97 & i<114)
      ypos=y-((i-96)*_butheight)
      'draw button 'buttonnr' 'x6' 'ypos' '_butlength' '_butheight' 'subwrd(_times._fno,i)
    endif
    i=i+1
  endwhile
  'q pos'
  btn=subwrd(result,7)
  if (_trace='on'); say result; endif;
  time._fno=_time._fno
  if (btn=-999); rc = delmen(x,y,i,mainbut,time._fno); endif
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,time._fno)
    rc = choice(btn)
  endif
  _time._fno=time._fno
  if (btn>=94 & btn<=99); rc = mainmenu(); endif  
  if (btn>20 & btn<=buttonnr)
    _time._fno=subwrd(_times._fno,btn-20)
    _time1._fno=(btn-20)
    'set t '_time1._fno
    if (_trace='on'); say 'time = '_time._fno; endif;
  endif
  if (_timeupdate='on')
    i=1
    while (i <= _filecount)
      _time.i=_time._fno
      _time1.i=_time1._fno
      i=i+1
    endwhile
  endif
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif;

 return
*
* *********** Draws the palette-menu and selects the desired palette ****************
*
 function palette()

  if (_trace='on'); say 'Color palette menu ...'; endif
*
  palette=''
  'run palette.gs' palette
   _palette._fno=sublin(result,1)
  if (_palette._fno='')
    _palette._fno='rainbow'
    'set rbcols auto'
  endif
* 
  if (_palettenupdate='yes')
    i=1
    while (i <= _filecount)
      _palette.i=_palette._fno
      i=i+1
    endwhile
  endif
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif

 return
*
* ************ Draws the region menu and selects the desired region *****************
*
 function region()

  _proj._fno='latlon'
  _plot._fno='latlon'
  _z=1
  mainbut = 8
  i=9
  x=9.1
  y=8.3
  'draw button 21 9.1 7.9 '_butlength' '_butheight' Arctic Ocean'
  'draw button 22 9.1 7.5 '_butlength' '_butheight' Antarctica'
  'draw button 23 9.1 7.1 '_butlength' '_butheight' Europe'
  'draw button 24 9.1 6.7 '_butlength' '_butheight' Germany'
  'draw button 25 9.1 6.3 '_butlength' '_butheight' Asia'
  'draw button 26 9.1 5.9 '_butlength' '_butheight' Africa'
  'draw button 27 9.1 5.5 '_butlength' '_butheight' North America'
  'draw button 28 9.1 5.1 '_butlength' '_butheight' South America'
  'draw button 29 9.1 4.7 '_butlength' '_butheight' Australia'
  'draw button 30 9.1 4.3 '_butlength' '_butheight' Monsoon'
  'draw button 31 9.1 3.9 '_butlength' '_butheight' Sahelzone'
  'draw button 32 9.1 3.5 '_butlength' '_butheight' Indian Region'
  'draw button 33 9.1 3.1 '_butlength' '_butheight' Mediteranian'
  'q pos'
   if (_trace='on'); say result; endif;
   btn=subwrd(result,7)
   reg._fno=_reg._fno
   if (btn=-999); rc = delmen(x,y,i,mainbut,reg._fno); endif
   if (btn>=1 & btn<=9)
     rc = delmen(x,y,i,mainbut,reg._fno)
     rc = choice(btn)
  endif
  _reg._fno = reg._fno
  if (btn>=94 & btn<=99); rc = mainmenu(); endif 
  if (btn=21)
    _proj._fno='nps'
    _reg._fno='Arctic Ocean'
  endif
  if (btn=22)
    _proj._fno='sps' 
    _reg._fno='Antarctica'
  endif
  if (btn=23)
    _reg._fno='Europe'
    _lolow._z=-20
    _lohigh._z=50
    _ltlow._z=35
    _lthigh._z=75
  endif
  if (btn=24)
    _reg._fno='Germany'
    _lolow._z=5
    _lohigh._z=18
    _ltlow._z=46
    _lthigh._z=56
  endif
  if (btn=25)
    _reg._fno='Asia'
    _lolow._z=50
    _lohigh._z=190
    _ltlow._z=-20
    _lthigh._z=85
  endif
  if (btn=26)
    _reg._fno='Africa'
    _lolow._z=-30
    _lohigh._z=60
    _ltlow._z=-40
    _lthigh._z=40
  endif
  if (btn=27)
    _reg._fno='North America'
    _lolow._z=-170
    _lohigh._z=-50
    _ltlow._z=15
    _lthigh._z=80
  endif
  if (btn=28)
    _reg._fno='South America'
    _lolow._z=-110
    _lohigh._z=-20
    _ltlow._z=-60
    _lthigh._z=15
  endif
  if (btn=29)
    _reg._fno='Australia'
    _lolow._z=110
    _lohigh._z=160
    _ltlow._z=-45
    _lthigh._z=-10
  endif
  if (btn=30)
    _reg._fno='Monsoon'
    _lolow._z=-20
    _lohigh._z=180
    _ltlow._z=-30
    _lthigh._z=50
  endif
  if (btn=31)
    _reg._fno='Sahelzone'
    _lolow._z=-20
    _lohigh._z=55
    _ltlow._z=3
    _lthigh._z=25
  endif
  if (btn=32)
    _reg._fno='Indian Region'
    _lolow._z=65
    _lohigh._z=110
    _ltlow._z=0
    _lthigh._z=40
  endif
  if (btn=33)
    _reg._fno='Mediteranian'
    _lolow._z=-10
    _lohigh._z=40
    _ltlow._z=25
    _lthigh._z=50
  endif
  if (_trace='on'); say 'Region = '_reg._fno; endif;
  if (_regionupdate='yes')
    while (i <= _filecount)
      _proj.i=_proj._fno
      _reg.i=reg._fno
      i=i+1
    endwhile
  endif
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif;

 return
*
* ************ Draws print menu and selects color or greyscale print ****************
*
 function filemenu()

  if (_trace='on'); say 'File menu'; endif;
  mainbut=9
  i=1
  x=10.3
  y=8.3
  while(i <= _filecount)
    buttonnr=20+i
    ypos=y-(i*_butheight)
    'draw button 'buttonnr' 'x' 'ypos' '_butlength' '_butheight' '_ctlfile.i
    i=i+1
  endwhile
  'q pos'
  if (_trace='on'); say result; endif;
  btn=subwrd(result,7)
  file._fno=_file._fno
  if (btn=-999); rc = delmen(x,y,i,mainbut,file._fno); endif
  if (btn>=1 & btn<=9)
    rc = delmen(x,y,i,mainbut,file._fno)
    rc = choice(btn)
  endif
  _file._fno=file._fno
  if (btn>=94 & btn<=99); rc = mainmenu(); endif
  if (btn>20 & btn <= buttonnr)
    _fno=btn-20
    _file._fno=_ctlfile._fno
    'set dfile '_fno
  endif
  if (_trace='on'); say 'file = '_file._fno; endif;
  if (_update='yes')
    rc = display()
  else
    rc = mainmenu()
  endif

 return

*
* **************************** Deletes the Submenu **********************************
*
 function delmen(x,y,i,mainbut,choice)

  if (_trace='on'); say 'Delete a submenu ...'; endif;
  'set line 0'
  if (i<18); k=1; endif
  if (i>17 & i<34); k=2; endif
  if (i>33); k=3; endif
  if mainbut=6
    'draw recf '%(x-_butlength/2)' '%(y-4*_butheight)' '%(k*x+_butlength)' '%(y-_butheight/2)
    'draw button 'mainbut' 'x' 'y' '_butlength' '_butheight' 'choice
  endif  
  if mainbut=1
    'draw recf '%(x-_butlength/2)' '%(y-4*_butheight)' '%(k*x+_butlength)' '%(y-_butheight/2)
    'draw button 'mainbut' 'x' 'y' '_butlength' '_butheight' 'choice
  else
    'draw recf '%(x-_butlength/2)' '%(y-18*_butheight)' '%(k*x+_butlength)' '%(y-_butheight/2)
    'draw button 'mainbut' 'x' 'y' '_butlength' '_butheight' 'choice
  endif

 return
*
* ****************** Animate time/level of a desired variable ************************
*
 function animate()

  if (_trace='on'); say 'Animate ...'; endif;
  if (_btn=83)
    if (_plot._fno='meridional' | _plot._fno='zonal')
      say 'NOT POSSIBLE !! Vertical coordinate is already variable !!'
    else 
     'set dbuff on'
      count = 1
      temp._fno=_clev._fno
      while (count <= _znum._fno)
        'set lev 'subwrd(_levs._fno,count)
        if (_null='off')
          'set black -1.0E-10 1.0E-10'
        endif
        if (_gxout='vector' | _gxout='stream')
          'd u;v'
        else
          'd '_cvar._fno
        endif
        if (_gxout='shaded' | _gxout='grfill'); 'run cbarn.gs 0.8'; endif;
        'set clip 1.5 9.0 1.5 7.0'
        'swap'
        say subwrd(_levs._fno,count)
        count=count+1
        say 'Press enter to see next level (first q and then enter to break ) !'
        pull answer
        if (answer = 'q'); break; endif;
      endwhile
*      _clev._fno=temp._fno
      'set lev '_clev._fno
      'set dbuff off'
      _animate='yes'
      rc = display() 
    endif
  endif
  if (_btn=84)
    if (_plot._fno='time')
      say 'NOT POSSIBLE !! Time is already variable !!'
    else
      'set dbuff on'
      'set clip 1.5 9.0 1.5 7.0'
      'set t '1' '_tnum._fno
      'set loopdim t'
      _animate='no'
      rc = display()
      'set looping off'
     'set t '_time1._fno
      rc = mainmenu()
      'set dbuff off'
     _animate='yes'
      rc = display() 
    endif
  endif  

 return
*
* ********************************** Push Upward *************************************
*
 function pushup()

  if(_proj._fno = 'latlon')
    _push='on'
    if (_z<2)
      latfactor = 5
    else
      latfactor = 2
    endif
    if (_z>0 & (_lthigh._z+latfactor)<90)
      _ltlow._z = _ltlow._z+latfactor
      _lthigh._z = _lthigh._z+latfactor
      rc = display()
    else
      say 'Operation not possible !!'
    endif
  else
    say 'Operation not possible with '_proj._fno
  endif

 return
*
* ******************************** Push left ******************************************
*
 function pushleft()

  if (_proj._fno = 'latlon' | _proj._fno='ortho' | _proj._fno='mollweide')
    _push='on'
    lonfactor=(_lohigh._z-_lolow._z)/10
    _lolow._z = _lolow._z-lonfactor
    _lohigh._z = _lohigh._z-lonfactor
    rc = display()
  else
    say 'Operation not possible with '_proj._fno
  endif

 return
*
* ********************************** Fast left Push ********************************
*
 function fastleft()

  if (_proj._fno = 'latlon' | _proj._fno='ortho' | _proj._fno='mollweide')
    _push='on'
    lonfactor = (_lohigh._z-_lolow._z)/4
    _lolow._z = _lolow._z-lonfactor
    _lohigh._z = _lohigh._z-lonfactor
    rc = display()
  else
    say 'Operation not possible with '_proj._fno
  endif

 return
*
* ************************************ Push Right **********************************
*
 function pushrig()

  if (_proj._fno = 'latlon' | _proj._fno='ortho' | _proj._fno='mollweide')
    _push='on'
    lonfactor=(_lohigh._z-_lolow._z)/10
    _lolow._z = _lolow._z+lonfactor
    _lohigh._z = _lohigh._z+lonfactor
    rc = display()
*    _push='off'
  else
    say 'Operation not possible with '_proj._fno
  endif

 return
*
* ******************************* Fast right Push *********************************
*
 function fastrig()

  if (_proj._fno = 'latlon' | _proj._fno='ortho' | _proj._fno='mollweide')
    _push='on'
    lonfactor=(_lohigh._z-_lolow._z)/4
    _lolow._z = _lolow._z+lonfactor
    _lohigh._z = _lohigh._z+lonfactor
    rc = display()
  else
    say 'Operation not possible with '_proj._fno
  endif
 
 return
*
* ********************************** Push Downward *********************************
*
function pushdown()

  if(_proj._fno = 'latlon')
    _push='on'
    if (_z<2)
      latfactor = 5
    else
      latfactor = 2
    endif 
    if (_z>0 & (_ltlow._z-latfactor)>-90)
      _ltlow._z = _ltlow._z-latfactor
      _lthigh._z = _lthigh._z-latfactor
      rc = display()
    else
      say 'Operation not possible !!'
    endif
  else
    say 'Operation not possible with '_proj._fno
  endif

 return
*
* ********* Returns the position of in three coordinate-systems ********************
*
 function positio()

  if (_trace='on'); say 'Position-Function'; endif;
  if (_proj._fno = 'latlon')
    'q pos'
    if (_trace='on'); say result; endif;
    res = result
    if (_plot._fno='latlon')
      xinch._fno = subwrd(res,3)
      yinch._fno = subwrd(res,4)
      say 'INCH-Position:'
      say ' X='xinch._fno
      say ' Y='yinch._fno
      'q xy2w 'xinch._fno' 'yinch._fno
      lonpos._fno = subwrd(result,3)
      latpos._fno = subwrd(result,6)
      say 'WORLD-Position:'
      say ' LON='lonpos._fno
      say ' LAT='latpos._fno
      'q xy2gr 'xinch._fno' 'yinch._fno
      longrid._fno = subwrd(result,3)
      latgrid._fno = subwrd(result,6)
      say 'GRID-Position:'
      say ' XGR='longrid._fno
      say ' YGR='latgrid._fno
    endif
    if (_plot._fno='meridional')
      yinch._fno = subwrd(res,3)
      zinch._fno = subwrd(res,4)
      say 'INCH-Position:'
      say ' Y='yinch._fno
      say ' Z='zinch._fno
      'q xy2w 'yinch._fno' 'zinch._fno
      latpos._fno = subwrd(result,3)
      levpos._fno = subwrd(result,6)
      say 'WORLD-Position:'
      say ' LAT='latpos._fno
      say ' LEV='levpos._fno
      'q xy2gr 'yinch._fno' 'zinch._fno
      latgrid._fno = subwrd(result,3)
      levgrid._fno = subwrd(result,6)
      say 'GRID-Position:'
      say ' YGR='latgrid._fno
      say ' ZGR='levgrid._fno
    endif
    if (_plot._fno='zonal')
      xinch._fno = subwrd(res,3)
      zinch._fno = subwrd(res,4)
      say 'INCH-Position:'
      say ' X='xinch._fno
      say ' Z='zinch._fno
      'q xy2w 'xinch._fno' 'zinch._fno
      lonpos._fno = subwrd(result,3)
      levpos._fno = subwrd(result,6)
      say 'WORLD-Position:'
      say ' LON='lonpos._fno
      say ' LEV='levpos._fno
      'q xy2gr 'xinch._fno' 'zinch._fno
      longrid._fno = subwrd(result,3)
      levgrid._fno = subwrd(result,6)
      say 'GRID-Position:'
      say ' XGR='longrid._fno
      say ' ZGR='levgrid._fno
    endif    
  else
    say 'Function only works correctly with latlon projection !!!'
  endif
  
 return
*
* ********************************* Display data ***********************************
* 
 function display()

  if (_trace='on')
    say 'Display data ...'
    say 'Projection = '_proj._fno
    say 'Plot = '_plot._fno
    if (_plot._fno='zonal')
      say 'lonmin ='_lolow._z
      say 'lonmax ='_lohigh._z
      say 'latitude='_lat._fno
    endif
    if (_plot._fno='meridional')
      say 'latmin ='_ltlow._z
      say 'latmax ='_lthigh._z
      say 'longitude='_lon._fno
    endif
    if (_plot._fno='meridional' | _plot._fno='zonal' | _plot._fno='hovm_t_z')
      say 'Minlevel = '_lolev._z
      say 'Maxlevel = '_highlev._z
    else
      say 'Level  = '_clev._fno
    endif
    if (_plot._fno='time' | _plot._fno='hovm_t_z')
      say 'MinTime = 'subwrd(_times._fno,1)
      say 'MaxTime = 'subwrd(_times._fno,_tnum._fno)
    else
      say 'TIME = 'subwrd(_times._fno,_time1._fno)
    endif
  endif
*
* For global fields latitude should begin with -90 and end with 90
* 
  if (_z = 0 & _lolow._fno < -85 & _lohigh._fno > 85)
    if (_plot._fno='latlon' | _plot._fno='meridional')
      'set lat -90 90'
    endif
  endif
  if (_logz._fno='on')
    'set zlog on'
    if (_logzupdate='yes')
      i=1
      while (i <= _filecount)
        _logz.i=_logz._fno
        i=i+1
      endwhile
    endif
  else
    'set zlog off'
  endif
*
  rc = project()
*
  if (_plot._fno='latlon' & _resflag='off')
    if (_lohigh._z-_lolow._z<10 | _lthigh._z-_ltlow._z<10)
      _res._fno='hires'
    else
      if (_lohigh._z-_lolow._z<30 | _lthigh._z-_ltlow._z<30)
        _res._fno='mres'
      else
        _res._fno='lowres'
      endif
    endif
  endif 
  'set mpdset '_res._fno
*
  if (_btn != 84 ) 
    rc = mainmenu()
  endif
  'set mpvals off'
  'set grads off'
  'set parea 1.6 9.3 1.5 7.0'
*
  if (_autopal='on' & _plot._fno != 'time')
    cvar._fno=_cvar._fno
    'run autopal.gs 'cvar._fno
    'run stefpal.gs 'cvar._fno
    _cvar._fno=cvar._fno
  endif
  if (_trace='on')
    say 'lonmin='_lolow._z
    say 'lonmax='_lohigh._z
    say 'latmin='_ltlow._z
    say 'latmax='_lthigh._z
  endif
  if (_plot._fno='hovm_t_lon'); 'set xyrev on'; endif
  if (_plot._fno='hovm_t_lat' & _ysize._fno>1); 'set xyrev on'; endif
  if (_null = 'off')
    'set black -1.0E-10 1.0E-10'
  endif   
  if (_gxout='vector' | _gxout='stream')
    if (_vec='uv' | _vec='u10v10')
      if (_z = 0); 'set strmden 5'; endif
      if (_z >= 1); 'set strmden 3'; endif
      'set arrlab off'
      if (_vec='uv')
        'd u;v;mag(u,v)'
      else
        'd u10;v10;mag(u10,v10)'
      endif
    endif
  else
    'd '_cvar._fno
  endif
  rc = title()

 return
*
* **************************** Projection function *******************************
*
 function project()

  if (_xsize._fno>1 & _ysize._fno>1 & _plot._fno='latlon')
    if (_proj._fno='latlon')
      'set mproj latlon'
      if (_trace='on'); say 'PUSH = '_push; endif
      if (_z=0 & _push='off')
        _lolow._z=_templonmin
        _lohigh._z=_templonmax
        _ltlow._z=_templatmin
        _lthigh._z=_templatmax
        'set lon '_lolow._z' '_lohigh._z
        'set lat '_ltlow._z' '_lthigh._z
      endif
      if (z>0)
        'set lon '_lolow._z' '_lohigh._z
        'set lat '_ltlow._z' '_lthigh._z
      endif
      if (_znum._fno > 0)
        'set lev '_clev._fno
      endif
      if(_btn != 84)
        'set t '_time1._fno
      endif      
    endif
    if (_proj._fno='nps')
      'set mproj nps'
      _lolow._z=-180
      _lohigh._z=180
      _lthigh._z=90        
      if (_z=0)
        if (_xsize._fno=64) ; _ltlow._z=19.382; endif;
        if (_xsize._fno=96) ; _ltlow._z=20.411; endif;
        if (_xsize._fno=128); _ltlow._z=18.139; endif;
        if (_xsize._fno=192); _ltlow._z=19.585; endif;
        if (_xsize._fno=320); _ltlow._z=19.626; endif;
        if (_xsize._fno=360); _ltlow._z=19.500; endif;
        if (_xsize._fno=720); _ltlow._z=19.750; endif;
      endif
      if (_z=1)
        if (_xsize._fno=64) ; _ltlow._z=58.143; endif;
        if (_xsize._fno=96) ; _ltlow._z=62.232; endif;
        if (_xsize._fno=128); _ltlow._z=59.997; endif;
        if (_xsize._fno=192); _ltlow._z=60.620; endif;
        if (_xsize._fno=320); _ltlow._z=60.000; endif;
        if (_xsize._fno=360); _ltlow._z=59.500; endif;
        if (_xsize._fno=720); _ltlow._z=59.750; endif;
      endif
      'set lon '_lolow._z' '_lohigh._z
      'set lat '_ltlow._z' '_lthigh._z
      'set frame circle'
    endif
    if (_proj._fno='sps')
      'set mproj sps'
      _lolow._z=-180
      _lohigh._z=180
      _ltlow._z=-90
      if (_z=0)
        if (_xsize._fno=64) ; _lthigh._z=-19.382; endif;
        if (_xsize._fno=96) ; _lthigh._z=-20.411; endif;
        if (_xsize._fno=128); _lthigh._z=-18.139; endif;
        if (_xsize._fno=192); _lthigh._z=-19.585; endif;
        if (_xsize._fno=320); _lthigh._z=-19.626; endif;
        if (_xsize._fno=360); _lthigh._z=-19.500; endif;
        if (_xsize._fno=720); _lthigh._z=-19.750; endif;
      endif
      if (_z=1)
        if (_xsize._fno=64) ; _lthigh._z=-58.143; endif;
        if (_xsize._fno=96) ; _lthigh._z=-61.232; endif;
        if (_xsize._fno=128); _lthigh._z=-59.997; endif;
        if (_xsize._fno=192); _lthigh._z=-60.620; endif;
        if (_xsize._fno=320); _lthigh._z=-60.000; endif;
        if (_xsize._fno=360); _lthigh._z=-59.500; endif;
        if (_xsize._fno=720); _lthigh._z=-59.750; endif;
      endif
      'set lon '_lolow._z' '_lohigh._z
      'set lat '_ltlow._z' '_lthigh._z
      'set frame circle'
    endif
    if (_proj._fno='robinson')
      if (_z=0)
        'set mproj robinson'
        'set lat -90 90'
        'set lon -180 180'
        'set mpvals -180 180 -90 90'
      else
        say 'You have to unzoom before using robinson projection !!'
        _proj._fno='latlon'
      endif
    endif
    if (_proj._fno='mollweide')
      if (_z=0)
        'set mproj mollweide'
        'set lon '_lolow._z' '%_lolow._z+360
        'set lat -90 90'
      else
        say 'You have to unzoom before using mollweide projection !!'
        _proj._fno='latlon'
      endif 
    endif
    if (_proj._fno='ortho')
      'set mproj orthographic'
      'set lon '_lolow._z' '%_lolow._z+180
      'set lat -90 90'
    endif       
  endif
  if (_projupdate='yes')
    i=1
    while (i <= _filecount)
      _proj.i=_proj._fno
      i=i+1
    endwhile
  endif

 return 

*
* ************************ Set title,... on the script ***************************
*
 function title()
  if (_btn!=83 & _btn!=84)
    rc = info()
    'set string 1 l 6'
    'set strsiz 0.15 0.17'
    if (_gxout='vector')
      if (subwrd(_levs._fno,1) > 30)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' Vector Windfield'
        if (_znum._fno > 0)
          'set string 1 r 6'
          'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' '_clev._fno' hPa'
        endif
      else
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' Vector Windfield'
        if (_znum._fno > 0)
          'set string 1 r 6'
          'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Level: '_clev._fno
        endif
      endif
    endif
    if (_gxout='stream')
      if (subwrd(_levs._fno,1) > 30)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' Stream Windfield'
        if (_znum._fno > 0)
          'set string 1 r 6'
          'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' '_clev._fno' hPa'
        endif
      else
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' Stream Windfield'
        if (_znum._fno > 0)
          'set string 1 r 6'
          'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Level: '_clev._fno
        endif
      endif
    endif
    if (_gxout != 'vector' & _gxout != 'stream' & _plot._fno = 'latlon')
      if (subwrd(_levs._fno,1) > 30 & _clev._fno != 'surface' & _znum._fno > 1) 
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        if (_znum._fno > 0)
          'set string 1 r 6'
          'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' '_clev._fno' hPa'
        endif   
      else
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        if (_znum._fno > 0)
          'set string 1 r 6'
          'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Level: '_clev._fno
        endif
      endif
    endif
    if (_plot._fno='meridional')
      if (_lon._fno > 180)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        'set string 1 r 6'
        'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Longitude: '%(360-_lon._fno)%' W'
      endif
      if (_lon._fno < 0)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        'set string 1 r 6'
        'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Longitude: '_lon._fno' W'
      endif
      if (_lon._fno >= 0 & _lon._fno <= 180)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        'set string 1 r 6'
        'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Longitude: '_lon._fno' E'
      endif
      'set string 1 c 6 90'
      'draw string '%(_xlow._fno-0.5)%' '%(_ylow._fno+3)%' '_zunit._fno
      'set string 1 r 6 0'
    endif
    if (_plot._fno='zonal')
      if (_lat._fno > 0)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        'set string 1 r 6'
        'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Latitude: '_lat._fno' N'
      else
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno
        'set string 1 r 6'
        'draw string '_xhigh._fno' '%(_ylow._fno-0.5)%' Latitude: '%(0-_lat._fno)%' S'
      endif
      'set string 1 c 6 90'
      'draw string '%(_xlow._fno-0.5)%' '%(_ylow._fno+3)%' '_zunit._fno
      'set string 1 r 6 0'
    endif
    if (_plot._fno='time' | _plot._fno='radiosonde')
      if (_lon._fno > 180)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno' '%(360-_lon._fno)%' W'
      endif
      if (_lon._fno < 0)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno' '_lon._fno' W'
      endif
      if (_lon._fno >= 0 & _lon._fno <= 180)
        'draw string '_xlow._fno' '%(_ylow._fno-0.5)%' '%_vname._vno._fno' '_lon._fno' E'
      endif
       if (_lat._fno > 0)
        'draw string '%(_xlow._fno+4.3)%' '%(_ylow._fno-0.5)%' '_lat._fno' N'
      else
        'draw string '%(_xlow._fno+4.3)%' '%(_ylow._fno-0.5)%' '%(0-_lat._fno)%' S'
      endif
      if (_plot._fno='time')
        'draw string '%(_xlow._fno+6.2)%' '%(_ylow._fno-0.5)%' '_lev._fno' '_zunit._fno
      endif
      if (_plot._fno='radiosonde')
        'set string 1 c 6 90'
        'draw string '%(_xlow._fno-0.5)%' '%(_ylow._fno+3)%' '_zunit._fno
        'set string 1 r 6 0'
      endif
    endif
    if (_plot._fno='hovm_t_z') 
      if (_lon._fno > 180)
        'draw string '%(_xlow._fno+0.4)' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '%(360-_lon._fno)%' W'
      endif
      if (_lon._fno < 0)
        'draw string '%(_xlow._fno+0.4)%' '%(_ylow._fno-0.5)%' '%_vname._vno_fno'      '_lon._fno' W'
      endif
      if (_lon._fno >= 0 & _lon._fno <= 180)
        'draw string '%(_xlow._fno+0.4)%' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '_lon._fno' E'
      endif
       if (_lat._fno >= 0)
        'draw string '%(_xlow._fno+6)%' '%(_ylow._fno-0.5)%' '_lat._fno' N'
      else
        'draw string '%(_xlow._fno+6)%' '%(_ylow._fno-0.5)%' '%(0-_lat._fno)%' S'
      endif
      'set string 1 c 6 90'
      'draw string '%(_xlow._fno-0.5)%' '%(_ylow._fno+3)%' '_zunit._fno
      'set string 1 r 6 0'
    endif
    if (_plot._fno='hovm_t_lon')
      if (_lat._fno < 0)
        'draw string '%(_xlow._fno+0.4)%' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '_lat._fno' S'
      endif
      if (_lat._fno >= 0)
        'draw string '%(_xlow._fno+0.4)%' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '_lat._fno' N'
      endif
      if (_znum._fno > 0)
        'draw string '%(_xlow._fno+6)%' '%(_ylow._fno-0.5)%' '_lev._fno' '_zunit._fno
      endif
      'set string 1 c 6 90'
      'draw string '%(_xlow._fno-0.5)%' '%(_ylow._fno+3)%' Longitude'
      'set string 1 r 6 0'
    endif
    if (_plot._fno='hovm_t_lat')
      if (_lon._fno > 180)
        'draw string '%(_xlow._fno+0.4)' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '%(360-_lon._fno)%' W'
      endif
      if (_lon._fno < 0)
        'draw string '%(_xlow._fno+0.4)%' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '_lon._fno' W'
      endif
      if (_lon._fno >= 0 & _lon._fno <= 180)
        'draw string '%(_xlow._fno+0.4)%' '%(_ylow._fno-0.5)%' '%_vname._vno._fno'      '_lon._fno' E'
      endif
      if (_znum._fno > 0)
        'draw string '%(_xlow._fno+6)%' '%(_ylow._fno-0.5)%' '_lev._fno' '_zunit._fno
      endif
      'set string 1 c 6 90'
      'draw string '%(_xlow._fno-0.5)%' '%(_ylow._fno+3)%' Latitude'
      'set string 1 r 6 0'
    endif
  endif 
  if (_btn != 84 | _animate='yes')  
    'set string 1 c 6'
    if (_proj._fno='robinson')
      'set strsiz 0.2 0.25'
      'draw string 5.5 7.0 '%_title._fno
    else
      'draw title '_title._fno
    endif
    if (_gxout='shaded' | _gxout='grfill' | _gxout='vector' | _gxout='stream')
      'run cbarn.gs 0.8'
    endif
  endif
  'set strsiz 0.05 0.08'
  'set string 1 l 6'
  if (_machine != 'mac')
    'draw string '%_xlow._fno%' '%_ylow._fno-0.8' '_name' '_date
  else
    'draw string '%_xlow._fno%' '%_ylow._fno-0.8' '_name
  endif

 return
*
* ******************************* Zoom into Plot **********************
*
 function zoom()

  if (_trace='on')
    say 'Zoom function ...'
    say 'zoom = '_z
  endif
  rc = info()
  if (_proj._fno='latlon')
    _z=_z+1
    'set line 1 1 6'
    'set rband 1 box '_xlow._fno' '_ylow._fno' '_xhigh._fno' '_yhigh._fno
    'q pos'
    x1 = subwrd(result,3)
    y1 = subwrd(result,4)
    x2 = subwrd(result,8)
    y2 = subwrd(result,9)
*
*** Sort coordinates *
*
    if (x2<x1)
      temp = x1
      x1 = x2
      x2 = temp
    endif;
    if (y2<y1)
      temp = y1
      y1 = y2
      y2 = temp
    endif;
*
* ************** New minimum coordinate ******************
*
    'q xy2w 'x1' 'y1
    if (_plot._fno='latlon')
      _lolow._z=subwrd(result,3)
      _ltlow._z=subwrd(result,6)
    endif
    if (_plot._fno='meridional')
      _ltlow._z=subwrd(result,3)
      _lolev._z=subwrd(result,6)
    endif
    if (_plot._fno='zonal')
      _lolow._z=subwrd(result,3)
      _lolev._z=subwrd(result,6)
    endif
*
* ************** New maximum coordinate ********************
*
    'q xy2w 'x2' 'y2
    if (_plot._fno='latlon')
      _lohigh._z=subwrd(result,3)
      _lthigh._z=subwrd(result,6)
      'set lon '_lolow._z' '_lohigh._z
      'set lat '_ltlow._z' '_lthigh._z
    endif
    if (_plot._fno='meridional')
      _lthigh._z=subwrd(result,3)
      _highlev._z=subwrd(result,6)
      'set lat '_ltlow._z' '_lthigh._z
      'set lev '_lolev._z' '_highlev._z
    endif
    if (_plot._fno='zonal')
      _lohigh._z=subwrd(result,3)
      _highlev._z=subwrd(result,6)
      'set lon '_lolow._z' '_lohigh._z
      'set lev '_lolev._z' '_highlev._z
    endif
    rc = display()
  else
    say 'Zooming not possibe with '_proj._fno
  endif

 return
*
* ****************************** Zoom out of the Plot ******************************
*
 function unzoom()

  if (_trace='on'); say 'Unzoom function ...'; endif;
  if (_proj._fno='latlon')
    say 'ZOOM = '_z
    if (_z>0)
      _z=_z-1
      if (_plot._fno='latlon')
        'set lon '_lolow._z' '_lohigh._z
        'set lat '_ltlow._z' '_lthigh._z
      endif
      if (_plot._fno='meridional')
        'set lat '_ltlow._z' '_lthigh._z
        'set lev '_lolev._z' '_highlev._z
      endif
      if (_plot._fno='zonal')
        'set lon '_lolow._z' '_lohigh._z
        'set lev '_lolev._z' '_highlev._z
      endif
    else
      say 'Unzoom not possible !!!'
    endif
  else
    if (_proj._fno='nps' | _proj._fno='sps')
      if (_z=1)
        _z=_z-1
      else
        say 'Unzoom not possible !!!'
      endif
    endif
  endif
  rc = display()
 
 return
*
* **********************************************************************************
* ********************************* Print Data *************************************
* **********************************************************************************
*
* ********************************* Color Print ************************************
*
 function colprn()

  if (_trace='on'); say 'Color printing ...'; endif;
  rc = display()
  say 'Name for the print file ?'
  pull filename
  if (filename ='');filename='printdata';endif
  'enable print 'filename
  'print'
  'disable print 'filename
  if (_machine != 'mac')
    '!gxps -c -i 'filename' -o 'filename'.ps'
    '!rm 'filename
  else
    say 'It is not possible to start the postscript converter directly !'
    say 'You have to start gxeps seperately to produce a postscript file'
    say 'of your plot.'
  endif

 return
*
* ******************************** Greyscale Print *********************************
*
 function greyprn()

  if (_trace='on'); say 'Greyscale printing ...'; endif;
  rc = display()
  say 'Name for the print file ? '
  pull filename
  if (filename ='');filename='printdata';endif
  'enable print 'filename
  'print'
  'disable print 'filename
  if (_machine != 'mac')
    '!gxps -i 'filename' -o 'filename'.ps'
    '!rm 'filename
  else
    say 'It is not possible to start the postscript converter directly !'
    say 'You have to start gxeps seperately to produce a postscript file'
    say 'of your plot.'
  endif

 return
 


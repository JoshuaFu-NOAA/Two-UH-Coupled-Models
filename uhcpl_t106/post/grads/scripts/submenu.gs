**************************
function submenu(name)
**************************
*****************************************************************************
*  This GrADS function draws a submenu,  and then plots the selected data.
*  The name of the menu is given in the function call,  and this is then
*  searched for in the menudata file.  Data concerning the set-up of the
*  menu is then read.
*  The structure of the function is as follows:
*
*                        Initialize variables
*                                  |
*                                  |
*                                  V
*                            Draw title box<------------------------------
*                                  |                                      |
*                                  |                                      |
*                                  V                                      |
*      ---------------------->Draw buttons<---------------                |
*     |                            |                      |               |
*     |                            |                      |               |
*     |                            V                      |               |
*     |                     Wait for mouse click          |               |
*     |                            |                      |               |
*     |                            |                      |               |
*     |                            V                      |               |
*     |              Check which button is pressed        |               |
*     |              and set variables accordingly.       |               |
*     |                            |                      |               |
*     |                            |                      |               |
*     |                            V                      |               |
*   Neither<------Decide whether to plot or animate.---->Plot             |
*                                  |                                      |
*                                  |                                      |
*                                  V                                      |
*                              Clear screen                               |
*                              and do animation---------------------------
*
*
*****************************************************************************

****************************
*** Initialize variables ***

*  The path name for the  GrADS scripts and setup data.
_modules = '/afs/dkrz/pf/k/oracle/public/casino/grads/v1/bin'

*  The name and path of the menu data file.
*_menudata=_modules%'/menudata'
_menudata = 'menudata'

* This part of the program reads the information from the menudata file
* corresponding to the menu whose name is given in the function call.

* Find the correct entry *
* This loop is broken if the end of file is reached,  or if there is an error,
* since then the return code, 'code', will not be 0.

* 'data' contains the return code followed by one line of the file 'menudata'.
* 'entry' is the line of the file read.
* 'retcode' is the return code.

retcode=0
entry=0
while((entry!=name)&(retcode=0))
  data=read(_menudata)
  retcode=sublin(data,1)
  entry= sublin(data,2)
endwhile
* Error Handling *
if(retcode=2); say 'Menu data not found'; return; endif
if(retcode>0)
  say 'Error reading menu data. Return code: 'retcode
  return
endif

*** This routine reads all the neccessary information for the setup of
*** the menu from the menudata file.
*** This loop assumes all the variables are in a particular
*** format.

* Description file name  
data    =read(_menudata)
descript= sublin(data,2)

* Range file name.  (The range file contains infomation used to determine the
* default contour interval,  aswell as the units.)
data    =read(_menudata)
rangefil= sublin(data,2)

* The type of section (xy, xz or yz)
data    =read(_menudata)
section = sublin(data,2)

* The title of the menu
data    =read(_menudata)
menutitl =sublin(data,2)

* The number of buttons  
data    =read(_menudata)
numbuts = sublin(data,2)

*  This loop reads the title,  button label and variable name(s) corresponding
*  to each button.

button=1
while(button<=numbuts)
  data        =read(_menudata)
  title.button=sublin(data,2)
  data        =read(_menudata)
  label.button=sublin(data,2)
  data        =read(_menudata)
  line        =sublin(data,2)
  name.button =subwrd(line,1)
  vcomp.button=subwrd(line,2)
  button=button+1
endwhile
ret=close(_menudata)


*** standard variables ***
*** Use ************************************** Allowed values *************
animate =0
* Switching animation off and on               0,1
redraw  =1
* Switching title drawing off and on           0,1
print=0
* Switching prnting off and on                 0,1
buttonon=0
* Controlling which button is selected.        0(None selected) 1,2,..numbuts
side=1
* Controlling which side is selected.          1,2,3,4,5
* side=1:left side=2:right side=3:bottom side=4:top side=5:contour scaling
time    =1
* The index number of the time in the dataset  1,2,3,...lasttime
z       =1
* The index number of the z level              1,2,3,...lastz
width   =(4.43/numbuts)
* The width of the variable buttons            4.43/numbuts (always +ve)
mapproj ='Global'
* The map projection                           Global, North, South
projlab ='G->N'
* The projection button label                  G->N,  N->S,  S->G
graphics='shaded'
* The graphics type                            shaded, contour
gxlab   ='C->W'
* The graphics button label                    C->W,  W->C
btn=1000
* The number of the button pressed             All button numbers, and -999.
* The value 1000 is a dummy value,  to plot the buttons initially.
axisscl=2
* map scale increment/adjustment increment     1,2,3.... 
rotate=0
* Controls rotation of polar projections       0...180

** This section of code defines the default positions of the sides of the
** map,  and the increments by which they can be changed.  The third dimension
** (that which is altered with the buttons) is also defined,  along with its
** name and units.
** N.B. The defaults given must be multiples of their respective increments.
if(section='xy')
  deflt.1=0; deflt.2=360; deflt.3=-90; deflt.4=90
  xincr=10; yincr=10 
  thirddim='z';thirddnm='depth'; thirddun='m'
endif
if(section='xz')
  deflt.1=0   ;deflt.2=360 ;deflt.3=0   ;deflt.4=6000
  xincr=10; yincr=-200
  thirddim='y';thirddnm='latitude' ; thirddun=''
endif
if(section='yz')
  deflt.1=-90 ;deflt.2=90  ;deflt.3=0   ;deflt.4=6000
  xincr=10; yincr=-200
  thirddim='x';thirddnm='longitude'; thirddun=''
endif
* positinc is the change in the selected variable when the plus or minus
* button is clicked.
positinc=xincr


** This section of code sets the positions of the sides to their defaults.
sideno=1
while(sideno<5)
  position.sideno=deflt.sideno
  sideno=sideno+1
endwhile

** This section of code switches all the buttons off initially,  so that they
** can be turned on when each is found in the range file. 
button=1
while(button<=numbuts)
  switch.button=0
  button=button+1
endwhile

*** reinitialize grads ***

'reinit'

*** open descriptor file ***

'open 'descript

*** get number of x,y, and z levels and times from description file ***

'q file'; all=result
line=sublin(all,5)
lastx   =subwrd(line,3)
lasty   =subwrd(line,6)
lastz   =subwrd(line,9)
lasttime=subwrd(line,12)

*** get the number of variables (treating vectors as two scalars)
line  =sublin(all,6)
numvar=subwrd(line,5)

********************************************************************
*** This routine gets information about each variable from the   ***
*   'RAN' file.  'min', 'max' and 'nkol' are used to determine the *
*   contour increment.  'type' classifies the variable as a vector *
*   or a scalar.  'unit' specifies the unit.                       *
********************************************************************

record=1
while(record<=numvar)
  data=read (rangefil)
  line=sublin(data,2)
  retcode=sublin(data,1)
** Read variable name from file.
  variable=subwrd(line,1)
  button=1
** Search through variable names from menudata until the name is found.
  while(button<=numbuts)
** If the variable is found,  switch that button on,  and assign min.button
** and max.button to the values in the range file.  Read nkol,  the number
** of colours, aswell, and assign del.button,  the default contour interval,
** accordingly.  colscal.button,  the factor by which the contour interval
** is multiplied,   is initialized to 100%,  and type.button is set to 1 or
** 2 according to whether the variable name.button is a scalar or vector.
** unit.button is the unit.
    if(name.button=variable)
      switch.button=1
      min.button=subwrd(line,2)
      max.button=subwrd(line,3)
      nkol=subwrd(line,4)
      del.button=((max.button-min.button)/(nkol-1))
      colscal.button=100
      type.button=subwrd(line,5)
      unit.button=subwrd(line,6)
    endif
    button=button+1
  endwhile
  record=record+1
endwhile
ret=close (rangefil)

*** This routine reads the z levels,  and time information from the
*** description file,  so that these can be displayed on the screen.

*  levelcou is a counter for the number of z levels.
*  continue is a control variable,  either 1, or 0.
*  wordno is the word number counting from the left.
levelcou=1
continue=1
wordno=4
* This loop reads each line of the description file until ZDEF is reached,
* so that the z levels can then be read.
* The z levels are read,  starting at the fourth word along, and continuing
* until there are no more words on that line.  The next line is then read,
* and if the first word is not TDEF,  more zlevels are read,  until the end
* of the line is reached.  The next line is then read,  and so on.
while(firstwrd!='ZDEF')
  data=read (descript)
  retcode=sublin(data,1)
  line=sublin(data,2)
  firstwrd=subwrd(line,1)
*  Error handling:
  if(retcode=2); firstwrd='ZDEF'; continue=0; endif;
endwhile
*ijotta - merke: im DES-file bleibt weiterhin anzahl levels=0 ->
*                anzahl levels bei VARS musz auch 0 bleiben.
* es geht hier nur darum, lastz=1 zu setzen und anzahl worte in line auf 4.
if(continue!=0)
if(subwrd(line,2)='0'); line='ZDEF 1 LEVELS 0'; lastz=1; endif;
endif
*ijotte
while(continue)
  zlevel.levelcou=subwrd(line,wordno)
  if(zlevel.levelcou='')
    data=read (descript)
    line=sublin(data,2)
    wordno=0
    levelcou=levelcou-1
  endif
  if(zlevel.levelcou='TDEF'); continue=0; endif;
  wordno=  wordno  +1
  levelcou=levelcou+1
endwhile
* The last line read above is the TDEF line,  so the start time,  and
* time increment can be found below.
*ijotta - contiZ und continue kommen jeweils +1 h"oher aus Schleife!
strttime=subwrd(line,4)
incrword=subwrd(line,5)
incremen=substr(incrword,1,1)
interval=substr(incrword,2,2)
rtsymmud=substr(incrword,4,1)
if(rtsymmud!='')
  incremen=substr(incrword,1,2)
  interval=substr(incrword,3,2)
endif
strthour=0
strtmin=0
if(substr(strttime,3,1)=':')
  strthour=substr(strttime,1,2)
  strtmin=substr(strttime,4,2)
endif
contiZ=1
while(strtday!='Z')
  strtday=substr(strttime,contiZ,1)
  contiZ=contiZ+1
  if(strtday='')
    contiZ=1
    strtday='Z'
  endif
endwhile
continue=contiZ
while(strtmon!='jan' & strtmon!='feb' & strtmon!='mar' & strtmon!='apr' & strtmon!='may' & strtmon!='jun' & strtmon!='jul' & strtmon!='aug' & strtmon!='sep' & strtmon!='oct' & strtmon!='nov' & strtmon!='dec')
  strtmon=substr(strttime,continue,3)
  continue=continue+1
endwhile
continue=continue-1
strtday=0
if(continue-contiZ=1); strtday=substr(strttime,contiZ,1); endif;
if(continue-contiZ=2); strtday=substr(strttime,contiZ,2); endif;
continue=continue+3
if(substr(strttime,continue,3)='')
  strtyear=substr(strttime,continue,2)
else
  strtyear=substr(strttime,continue,4)
endif
if(strtmon='jan'); month=1; endif;
if(strtmon='feb'); month=2; endif;
if(strtmon='mar'); month=3; endif;
if(strtmon='apr'); month=4; endif;
if(strtmon='may'); month=5; endif;
if(strtmon='jun'); month=6; endif;
if(strtmon='jul'); month=7; endif;
if(strtmon='aug'); month=8; endif;
if(strtmon='sep'); month=9; endif;
if(strtmon='oct'); month=10; endif;
if(strtmon='nov'); month=11; endif;
if(strtmon='dec'); month=12; endif;
monstrt=month
*ijotte
ret=close (descript)  

*** Set some GrADs settings if the section is not xy ***
if(section!='xy')
  'set mproj scaled'
  'set frame on'
endif

************************************************************
********************** Main Loop ***************************
* The loop is repeated indefinitely,  until exit is selected.

 
while (1)
* plot is a switch for controlling plotting and animating.
  plot=1

  if(redraw)

*ijotta - um submenu.gs direkt (i.e. nicht "uber topmenu.gs und level_menu.gs)
*         aufrufen zu k"onnen
*** load color palette C32a

  'run '_modules'/C32a'

*** define greyscale rgb color indices

  'set rgb 88 150 150 150'
*ijotte

*** draw main frame in black

    'set line 0 1 6'
    'draw recf 0.0 0.0 11.0 8.3'

*** draw top section like a title box

    'set line 22 1 6'
    'draw recf 2.0 7.5 11.0 8.3'
    'set line 0 1 6'
    'draw recf 2.05 7.54 10.95 8.25'
    'set line 21 1 6'
    'draw recf 2.1 7.58 10.9 8.2'

    'set string 1 c 10'
    'set strsiz 0.43 0.45'
    'draw string 6.5 7.9 CERA Graphics'

*** draw section like a title box

    'set line 22 1 6'
    'draw recf 2.0 6.7 11.0 7.5'
    'set line 0 1 6'
    'draw recf 2.05 6.75 10.95 7.45'
    'set line 21 1 6'
    'draw recf 2.1 6.8 10.9 7.4'

    'set string 32 c 10'
    'set strsiz 0.36 0.4'
    'draw string 6.5 7.1 'menutitl

*** Set the viewport

    'set vpage off'
    'set line 7 1 6'
    'draw line  2.0  0.1 10.9 0.1'
    'draw line 10.9  0.1 10.9 6.6'
    'draw line 10.9  6.6  2.0 6.6'
    'draw line  2.0  6.6  2.0 0.1'
    'set vpage  2.0 10.9  0.1 6.6'
    'set parea  1.0 10.0  1.2 6.5'

  endif


*** Control which sides can be moved                                ***
**  Reinitialize plus and minus stops to zero.

   sideno=1
   while(sideno<6)
     minusstp.sideno=0
     plusstp.sideno =0
     sideno=sideno+1
   endwhile


**  Set plus and minus stops according to which section is selected.
*   If plusstp.side=1 then that side cannot be moved in the positive
*   direction,  and if minusstp.side=1 the side cannot be moved in the
*   negative direction.
   if(section='xy')
     if(position.1=deflt.1); minusstp.1=1; endif
     if(position.2=deflt.2); plusstp.2 =1; endif
     if(position.3=deflt.3); minusstp.3=1; endif
     if(position.4=deflt.4); plusstp.4 =1; endif
     if(position.2-position.1=xincr); plusstp.1=1; minusstp.2=1; endif
     if(position.4-position.3=yincr); plusstp.3=1; minusstp.4=1; endif
     if((mapproj='North')&(position.3=0)); minusstp.3=1;endif
     if((mapproj='South')&(position.4=0)); plusstp.4 =1;endif
   endif
   if(section='xz')
     if(position.1=deflt.1); minusstp.1=1; endif
     if(position.2=deflt.2); plusstp.2 =1; endif
     if(position.3=deflt.3); plusstp.3 =1; endif
     if(position.4=deflt.4); minusstp.4=1; endif
     if(position.2-position.1=xincr); plusstp.1 =1; minusstp.2=1; endif
     if(position.3-position.4=yincr); minusstp.3=1; plusstp.4 =1; endif
   endif
   if(section='yz')
     if(position.1=deflt.1); minusstp.1=1; endif
     if(position.2=deflt.2); plusstp.2 =1; endif
     if(position.3=deflt.3); plusstp.3 =1; endif
     if(position.4=deflt.4); minusstp.4=1; endif
     if(position.2-position.1=xincr); plusstp.1 =1; minusstp.2=1; endif
     if(position.3-position.4=yincr); minusstp.3=1; plusstp.4 =1; endif
   endif

*** Set position.5 for first time ***
   if(buttonon=0); position.5=100; endif
   if(position.5=10 ); minusstp.5=1; endif


********************
*** Draw buttons ***

*** Draw variable buttons with labels ***********************************
*   The buttons are drawn red if they are not available, blue if they are
*   available but not selected,  and green if selected.
*   The buttons are scaled to fit in the available space.
*************************************************************************
  button=1
  while(button<=numbuts)
    if(switch.button)
      if button=buttonon
        'set button 1 28 24 23 1 28 24 23 6'
      else
        'set button 1 17 20 23 0 17 20 23 6'
      endif
    else
      'set button 1 37 36 34 0 37 36 34 6'    
    endif
    position=(6.62-(width*(button-0.5)))
    'draw button 'button' 0.98 'position' 1.9 'width' 'label.button  
    button=button+1
  endwhile

*** Set and draw quit, level, animate, map projection and boundary buttons.
*** Colour the buttons according to whether they have been selected,
*** and whether that option is available.


  'set button 1 11 20 23 0 28 24 23 6'
  'draw button 99 0.98 0.35 1.9 0.5 Exit'

  if(time=1)
    'set button 1 37 36 34 0 37 36 34 6'
  else
    'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 91 0.3466 1.35 0.6333 0.5  t-1'
  if(time=lasttime)
    'set button 1 37 36 34 0 37 36 34 6'
  else
    'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 92  1.6133 1.35 0.6333 0.5 t+1'
  if(lasttime=1)
    'set button 1 37 36 34 0 37 36 34 6'
  else
    'set button 1 17 20 23 0 28 24 23 6'
  endif
  if buttonon=0
    'set button 1 37 36 34 0 37 36 34 6'
  endif
  'draw button 93  0.9800 1.35 0.6333 0.5 >>'

  if(lastz=1)
    'set button 17 17 20 23 17 17 20 23 6'
    'draw button 94 0.9800 0.85 0.6333 0.5 -'
    'draw button 95 0.9800 1.85 0.6333 0.5 -'
  else
    if(z=1)
      'set button 1 37 36 34 0 37 36 34 6'
    else
      'set button 1 17 20 23 0 28 24 23 6'
    endif
  'draw button 94 0.9800 0.85 0.6333 0.5 'thirddim'-1'
    if(z=lastz)
      'set button 1 37 36 34 0 37 36 34 6'
    else
      'set button 1 17 20 23 0 28 24 23 6'
    endif
  'draw button 95 0.9800 1.85 0.6333 0.5 'thirddim'+1'
  endif

  if((time=1)&(z=1))
    'set button 1 37 36 34 0 37 36 34 6'
  else
    'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 96 0.3466 0.85 0.6333 0.5  |<' 

  if(section='xy')
    'set button 1 17 20 23 0 28 24 23 6'
    'draw button 81  0.3466 1.85 0.6333 0.5 'projlab
  else
    'set button 17 17 20 23 17 17 20 23 6'
    'draw button 81  0.3466 1.85 0.6333 0.5 -'
  endif
  if(type.buttonon=2)
    'set button 1 37 36 34 0 37 36 34 6'
    graphics='shaded'
    gxlab='C->W'
  else
  'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 82  1.6133 1.85 0.6333 0.5 'gxlab

  if buttonon=0
    'set button 1 37 36 34 0 37 36 34 6'
  else
  'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 83  1.6133 0.85 0.6333 0.5 .ps'

 
  
  if side=1
    'set button 1 28 24 23 1 28 24 23 6'
  else
    'set button 1 17 20 23 0 17 20 23 6'
  endif 
  'draw button 61  0.3466 7.5000 0.6333 0.5333 Left'
  if side=2
    'set button 1 28 24 23 1 28 24 23 6'
  else
    'set button 1 17 20 23 0 17 20 23 6'
  endif 
  'draw button 62  1.6133 7.5000 0.6333 0.5333 Right'
  if side=3
    'set button 1 28 24 23 1 28 24 23 6'
  else
    'set button 1 17 20 23 0 17 20 23 6'
  endif 
  'draw button 63  0.980 6.9667 0.6333 0.5333 Bottm'
  if side=4
    'set button 1 28 24 23 1 28 24 23 6'
  else
    'set button 1 17 20 23 0 17 20 23 6'
  endif 
  'draw button 64  0.980 8.0333 0.6333 0.5333 Top'
  if side=5
    'set button 1 28 24 23 1 28 24 23 6'
  else
    'set button 1 17 20 23 0 17 20 23 6'
  endif
  if buttonon=0
    'set button 1 37 36 34 0 37 36 34 6'
  endif
  'draw button 65  0.980 7.5000 0.6333 0.5333 Incr'
  if(minusstp.side)
    'set button 1 37 36 34 0 37 36 34 6'
  else
    'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 66  0.3466 8.0333 0.6333 0.5333 -'
  if(plusstp.side)
    'set button 1 37 36 34 0 37 36 34 6'
  else
    'set button 1 17 20 23 0 28 24 23 6'
  endif
  'draw button 67  1.6133 8.0333 0.6333 0.5333 +'
  'set button 1 17 20 23 0 28 24 23 6'
  'draw button 68  0.3466 6.9667 0.6333 0.5333 |<'
  'set button 1 17 20 23 0 17 20 23 6'
* This section of code determines the data that should be displayed in the
* button at the top right.  If the side selected is the left or right, then
* the data displayed is the position of that side. If the side selected is
* the top or bottom,  then the information displayed is determined by the
* coordinate mapping between the GrADS coordinates and those diplayed. If
* contour interval is selected,  the scaling factor is displayed as a
* percentage.
  if(side<5)
    if(positinc=xincr)
      value=position.side
    else
* Coordinate mapping
      if(section='xy'); value=position.side          ;endif
      if(section='xz'); value=deflt.4-position.side;endif
      if(section='yz'); value=deflt.4-position.side;endif
    endif
  else
    value=position.side'%'
  endif
  'draw button 69  1.613 6.9667 0.6333 0.5333 'value
 

*** get information on the pressed button,  unless the screen is being
*** redrawn after an animation.

  if(redraw=0)

    'query pos'

*** get the number of pressed button to get a control criteria
    btn = subwrd(result,7)    

  else

* Reset redraw and plot if the screen is being redrawn after an animation.
    redraw=0
    plot=1
  endif

***************************************************
*********** Check buttons ************************* 

* Click not on any button.
  if(btn = -999); plot=0; endif;

* Exit button
*  if(btn = 99); 'set vpage off'; return; endif;
*ijotta
  if(btn = 99); 'quit'; endif;
*ijotte

* Variable buttons
  if((btn<20)&(btn>0))
    if(switch.btn)
      buttonon=btn
      position.5=colscal.buttonon
    else
      plot=0
    endif
  endif

*** Check animation buttons ***
* Animate cannot be selected if there is only one frame,  so the map is not
* redrawn in that case.
* N.B. The first if statement is included to make the program run faster. 
* The vpage command is included here,  so that other settings made later in
* the program are retained.

  if(btn>80)
* If print is selected
    if(btn=83)
      'set vpage off'
      'c'
    endif 
* t-1
    if(btn = 91);time=time-1;endif;
* t+1
    if(btn = 92);time=time+1;endif;
* >>
    if(btn = 93)
      if((lasttime>1)&(buttonon!=0))
        animate=1
        'set vpage off'
      else
        plot=0
      endif
    endif
* |<
    if(btn = 96)
      if((time=1)&(z=1))
        plot=0
      else
        time=1
        z   =1
      endif
    endif
* Prevent time from being changed if it as the end of its range.
    if(time<1);       time=1;       plot=0;endif;
    if(time>lasttime);time=lasttime;plot=0;endif;
*ijotta
* Set time,  and assign time.
    'set t 'time
if(interval='yr')
    year=(strtyear-incremen)+(time*incremen)
    month=monstrt
    day=strtday
endif
if(interval='mo')
    year=strtyear
    month=(monstrt-incremen)+(time*incremen)
    agecount=0
    while(month>12)
      month=month-12
      agecount=agecount+1
    endwhile
    year=year+agecount
    day=strtday
endif
*ijotte

*** Check Level Buttons ***
    if(btn = 94);z=z-1;endif;
    if(btn = 95);z=z+1;endif;
    if(z<1);    z=1;    plot=0;endif;
    if(z>lastz);z=lastz;plot=0;endif;
    'set z 'z
  endif

*** Check Boundary Buttons ***

*** Set plotting off if any of the boundary buttons is pressed      ***
  if((btn>60)&(btn<70))
    plot=0
  endif

* left
  if(btn=61); side=1;  positinc=xincr;endif
* right
  if(btn=62); side=2;  positinc=xincr;endif
* bottom
  if(btn=63); side=3;  positinc=yincr;endif
* top
  if(btn=64); side=4;  positinc=yincr;endif
* incr
  if((buttonon>0)&(btn=65)); side=5; positinc=10   ;endif
* -
  if((btn=66)&(minusstp.side=0))
    position.side=position.side-positinc
  endif
* +
  if((btn=67)&(plusstp.side=0))
    position.side=position.side+positinc
  endif

* If reset is selected,  return to original map projection and contour scaling.
* For xy sections,  reset map by simulating a press of the corresponding
* projection button.  For other sections,  return to defaults.
* |<
  if(btn=68)
    position.5=100
    if(section='xy')
      btn=81
      if(mapproj='North')
        mapproj='Global'
      else
        if(mapproj='South')
          mapproj='North'
        else
          if(mapproj='Global')
            mapproj='South'
          endif
        endif
      endif
    else
      counter=1
      while(counter<5)
        position.counter=deflt.counter
        counter=counter+1
      endwhile
    endif
  endif
  colscal.buttonon=position.5
  
*** Check projection buttons ***

* G->N, N->S, S->G
  if(btn=81)
    if(section='xy')
      if(mapproj='South')
        mapproj='Global'
        projlab='G->N'
        'set mproj scaled'
        'set frame on'
        position.1=0  ;position.2=360
        position.3=-90;position.4=90
      else
        if(mapproj='Global')
          mapproj='North'
          projlab='N->S'
          position.1=0; position.2=360
          position.3=50;position.4=90
          'set mproj nps'
          'set frame off'
        else
          if(mapproj='North')
            mapproj='South'
            projlab='S->G'
            position.1=0;  position.2=360
            position.3=-90;position.4=-50
            'set mproj sps'
            'set frame off'
          endif
        endif
      endif
    else
      plot=0
    endif
  endif

*** Check graphics buttons ***

* C->W, W->C
  if(btn=82)
    if(type.buttonon=2)
      plot=0
    else
      if(graphics='shaded')
        graphics='contour'
        gxlab   ='W->C'
      else
        graphics='shaded'
        gxlab   ='C->W'
      endif
    endif
  endif

*** Check print button ***
* Print
  if(btn=83)
    '!rm -f  print.tmp'
    'enable print print.tmp'
    print=1
  endif

*** Set plotting area **
  'set lon 'position.1' 'position.2
  'set lat 'position.3' 'position.4

  if(mapproj='Global')
    mpval.1=position.1;mpval.2=position.2;
    mpval.3=position.3;mpval.4=position.4
  endif
  if(mapproj='North')
    mpval.1=0-180-rotate;mpval.2=180-rotate;
    mpval.3=position.3-5;mpval.4=position.4
  endif
  if(mapproj='South')
    mpval.1=0-180-rotate;mpval.2=180-rotate;
    mpval.3=position.3  ;mpval.4=position.4+5
  endif


  'set mpvals 'mpval.1' 'mpval.2' 'mpval.3' 'mpval.4 



*** Determine axis labels (mapping GrADS coordinates to displayed coordinates)
  xaxisinc=xincr*axisscl
  yaxisinc=yincr*axisscl
  if(section='xy')
    bottmscl=position.3
    topscl=  position.4
  endif
  if(section='xz')
    bottmscl=deflt.4-position.3
    topscl=  deflt.4-position.4
  endif
  if(section='yz')
    bottmscl=deflt.4-position.3
    topscl=  deflt.4-position.4
  endif

*** Set axis labelling ***

  'set xaxis 'position.1' 'position.2' 'xaxisinc  
  'set yaxis 'bottmscl' 'topscl' 'yaxisinc


*************************
*** Plotting Routines ***

* Animate routine
  
* set contour colour to white
  'set ccolor 1'
  if(animate)
    'set looping on'
    'set loopdim t'
    'set t 1 'lasttime
    if(type.buttonon=2)
      'set gxout vector'
      arrscl=max.buttonon*colscal.buttonon/100
      'd 'name.buttonon' ; 'vcomp.buttonon
    else
      'set gxout 'graphics
      cint=del.buttonon*colscal.buttonon/100
      'set cint 'cint
      'd 'name.buttonon
    endif
    'set looping off'
    'set t 1'
    redraw=1
    animate=0
    time=1;
  endif

* This routine draws a static plot:
* 1) if the screen is being re-drawn after an animation,
* 2) if an animate or variable button has been clicked,
* but not immediately after an animation. 

  if((buttonon=0)|(redraw=1)|(btn>96));plot=0;endif;
  
  if(plot)
* Clear screen by drawing a black filled rectangle over graphics window.
    'set line 0'
    'draw recf 0.03 0.02 10.97 8'
    'set mpdraw off'
    if(type.buttonon=2)
      'set gxout vector'
      arrscl=max.buttonon*colscal.buttonon/100
      'set arrscl 0.2 'arrscl
      'd 'name.buttonon' ; 'vcomp.buttonon
    else
      'set gxout 'graphics
      cint=del.buttonon*colscal.buttonon/100
      'set cint 'cint
      'd 'name.buttonon
* Call the label bar routine if shaded graphics are selected.
      if(graphics='shaded'); 'run '_modules'/labbar.gs 'unit.buttonon; endif
    endif
    if(lastz=1)
      title=title.buttonon'  year:'year '  month:'month '  day:'day
    else
      title=title.buttonon'  year:'year' 'thirddnm': 'zlevel.z' 'thirddun
    endif
    'draw title 'title
  endif
  btn=0
*********************************************************
*** Print Routine ***************************************
  if(print)
* Print displayed graphics to metafile.
    'print'
    'disable print'
* Convert to postscript file.
    '!gxpscw -i print.tmp -o print.ps'
     say 'file print.ps created'
* Print file.
*    '!lpr print.ps'
*    say 'printing...'
*    say 'Printer queue:'
*    '!lpq'
    redraw=1
    print=0
    '!rm -f  print.tmp'
*    '!rm -f  print.ps'   
  endif
endwhile

return









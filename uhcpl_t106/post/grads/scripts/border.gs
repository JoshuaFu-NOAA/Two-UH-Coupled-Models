*****************************************************************
function boundary()
*****************************************************************
*------------------------------------------------------------
*------------------------------------------------------------
* get the map project info
'q gxinfo'
 map=sublin(result,6)
 yline=sublin(result,4)
 project=subwrd(map,3)
 ybottom=subwrd(yline,4)
 ytop=subwrd(yline,6)
*
* get the dimension boundary
'q dims'
  card=sublin(result,2)
  east=subwrd(card,6)
  west=subwrd(card,8)

  card=sublin(result,3)
  lowlat=subwrd(card,6)
  highlat=subwrd(card,8)
*------------------------------------------------------------
east1=east
west1=west
delta=west1-east1
if(delta >= 180)
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(delta = 360)
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
east=east1
'q w2xy 'east' 'lowlat
xold=subwrd(result,3)
yold=subwrd(result,6)
lowlat2=lowlat-15
if(project=4 & lowlat<=-75);lowlat2=-90;endif
'q w2xy 'east' 'lowlat2
xold2=subwrd(result,3)
yold2=subwrd(result,6)
*********
while(east<=west)
'q w2xy 'east' 'lowlat
x=subwrd(result,3)
y=subwrd(result,6)
'q w2xy 'east' 'lowlat2
xnew=subwrd(result,3)
ynew=subwrd(result,6)
'set line 0'
'draw polyf 'xold' 'yold' 'xold2' 'yold2' 'xnew' 'ynew' 'x' 'y
'set line 1 1 4'
'draw line 'xold' 'yold' 'x' 'y
xold=x
yold=y
xold2=xnew
yold2=ynew
east=east+2.5
endwhile
*****************************
* wipe out the highlat circle
east=east1
'q w2xy 'east' 'highlat
xold=subwrd(result,3)
yold=subwrd(result,6)
highlat2=highlat+15
if(project = 3 & highlat >= 75) ; highlat2=90 ; endif
xold2=subwrd(result,3)
yold2=subwrd(result,6)
*****************
while(east<=west)
'q w2xy 'east' 'highlat
x=subwrd(result,3)
y=subwrd(result,6)
'q w2xy 'east' 'highlat2
xnew=subwrd(result,3)
ynew=subwrd(result,6)
'set line 0'
'draw polyf 'xold' 'yold' 'xold2' 'yold2' 'xnew' 'ynew' 'x' 'y
'set line 1 1 4'
'draw line 'xold' 'yold' 'x' 'y
xold=x
yold=y
xold2=xnew
yold2=ynew
east=east+2.5
endwhile
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
else
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
east=east1
'q w2xy 'east' 'lowlat
xold=subwrd(result,3)
yold=subwrd(result,6)
lowlat2=lowlat-15
if(project=4 & lowlat<=-75);lowlat2=-90;endif
'q w2xy 'east' 'lowlat2
xold2=subwrd(result,3)
yold2=subwrd(result,6)
*********
while(east<=west)
'q w2xy 'east' 'lowlat
x=subwrd(result,3)
y=subwrd(result,6)
'q w2xy 'east' 'lowlat2
xnew=subwrd(result,3)
ynew=subwrd(result,6)
'set line 0'
'draw polyf 'xold' 'yold' 'xold2' 'yold2' 'xnew' 'ynew' 'x' 'y
'set line 1 1 4'
'draw line 'xold' 'yold' 'x' 'y
xold=x
yold=y
xold2=xnew
yold2=ynew
east=east+2.5
endwhile
*****************************
* wipe out the highlat circle
east=east1
'q w2xy 'east' 'highlat
xold=subwrd(result,3)
yold=subwrd(result,6)
highlat2=highlat+15
if(project = 3 & highlat >= 75) ; highlat2=90 ; endif
xold2=subwrd(result,3)
yold2=subwrd(result,6)
*****************
while(east<=west)
'q w2xy 'east' 'highlat
x=subwrd(result,3)
y=subwrd(result,6)
'q w2xy 'east' 'highlat2
xnew=subwrd(result,3)
ynew=subwrd(result,6)
'set line 0'
'draw polyf 'xold' 'yold' 'xold2' 'yold2' 'xnew' 'ynew' 'x' 'y
'set line 1 1 4'
'draw line 'xold' 'yold' 'x' 'y
xold=x
yold=y
xold2=xnew
yold2=ynew
east=east+2.5
endwhile
**********************
*wipe out the west part
'q w2xy 'west' 'lowlat2
x1=subwrd(result,3)
y1=subwrd(result,6)
west2=west+360-delta
'q w2xy 'west2' 'lowlat2
xlf=subwrd(result,3)
ylf=subwrd(result,6)
'q w2xy 'west' 'highlat2
x2=subwrd(result,3)
y2=subwrd(result,6)
'q w2xy 'west2' 'highlat2
xlf1=subwrd(result,3)
ylf1=subwrd(result,6)
'set line 0'
'draw polyf 'x1' 'y1' 'xlf' 'ylf' 'xlf1' 'ylf1' 'x2' 'y2
'q w2xy 'west' 'lowlat
x1=subwrd(result,3)
y1=subwrd(result,6)
'q w2xy 'west' 'highlat
x2=subwrd(result,3)
y2=subwrd(result,6)
'set line 1 1 4'
'draw line 'x1' 'y1' 'x2' 'y2
****************
*wipe out the east part
'q w2xy 'east1' 'lowlat
x1=subwrd(result,3)
y1=subwrd(result,6)
'q w2xy 'east1' 'highlat
x2=subwrd(result,3)
y2=subwrd(result,6)
'set line 1 1 4'
'draw line 'x1' 'y1' 'x2' 'y2
*++++++++++++++++++++++++++++++
endif
*++++++++++++++++++++++++++++++
*----------------------------------------------------------------------
* for normal nps
*----------------------------------------------------------------------
else
*wipe out the lowlat circle
east=east1
'q w2xy 'east' 'lowlat
xold=subwrd(result,3)
yold=subwrd(result,6)
ybt=ybottom-0.6
*****************
while(east<=west)
'q w2xy 'east' 'lowlat
x=subwrd(result,3)
y=subwrd(result,6)
'set line 0'
'draw polyf 'xold' 'yold' 'xold' 'ybt' 'x' 'ybt' 'x' 'y
'set line 1 1 4'
'draw line 'xold' 'yold' 'x' 'y
xold=x
yold=y
east=east+2.5
endwhile
*****************************
* wipe out the highlat circle
east=east1
'q w2xy 'east' 'highlat
xold=subwrd(result,3)
yold=subwrd(result,6)
ytp=ytop+0.6
*****************
while(east<=west)
'q w2xy 'east' 'highlat
x=subwrd(result,3)
y=subwrd(result,6)
'set line 0'
'draw polyf 'xold' 'yold' 'xold' 'ytp' 'x' 'ytp' 'x' 'y
'set line 1 1 4'
'draw line 'xold' 'yold' 'x' 'y
xold=x
yold=y
east=east+2.5
endwhile
*******************************************
*wipe out the east part
'q w2xy 'east1' 'lowlat
x1=subwrd(result,3)
y1=subwrd(result,6)
'q w2xy 'east1' 'highlat
x2=subwrd(result,3)
y2=subwrd(result,6)
if(project = 3)
'set line 0'
'draw polyf 'x1' 'y1' 'x1' 'y2' 'x2' 'y2
'set line 1 1 4'
'draw line 'x1' 'y1' 'x2' 'y2
else
'set line 0'
'draw polyf 'x1' 'y1' 'x2' 'y1' 'x2' 'y2
'set line 1 1 4'
'draw line 'x1' 'y1' 'x2' 'y2
endif
**************************************************************
*wipe out the west part
'q w2xy 'west1' 'lowlat
x1=subwrd(result,3)
y1=subwrd(result,6)
'q w2xy 'west1' 'highlat
x2=subwrd(result,3)
y2=subwrd(result,6)
if(project = 3)
'set line 0'
'draw polyf 'x1' 'y1' 'x1' 'y2' 'x2' 'y2
'set line 1 1 4'
'draw line 'x1' 'y1' 'x2' 'y2
else
'set line 0'
'draw polyf 'x1' 'y1' 'x2' 'y1' 'x2' 'y2
'set line 1 1 4'
'draw line 'x1' 'y1' 'x2' 'y2
endif
endif
return

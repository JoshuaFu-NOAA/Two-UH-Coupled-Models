function page (args)
*
* run page rows cols gnum
* or
* run page rows cols row col
*
*
rows=subwrd(args,1)
cols=subwrd(args,2)
gnum=subwrd(args,3)
col=subwrd(args,4)
row=gnum
*
if(rows='');say "syntax: run pagem rows cols gnum";return;endif
maxg=rows*cols
if (col='')
  while (gnum > maxg) ; gnum = gnum - maxg ; endwhile
  row=1
  dum=gnum
  while (dum > cols);dum=dum-cols;row=row+1;endwhile
  col=gnum-cols*(row-1)
  say 'row='%row%' col='%col
endif
*
if (rows=5 & cols=4)
xoffl=0.05
xoffr=0.8
yoffb=0.15
yofft=0.8
xsc=0.05
ysc=0.05
endif
*
if (rows=4 & cols=4)
xoffl=0.05
xoffr=0.05
yoffb=0.65
yofft=0.8
xsc=0.05
ysc=0.01
endif
*
if (rows=3 & cols=4)
xoffl=0.12
xoffr=0.12
yoffb=1.0
yofft=1.0
xsc=0.12
ysc=-0.17
endif
*
if (rows=3 & cols=3)
xoffl=0.1
xoffr=0.99
yoffb=0.4
yofft=0.8
xsc=0.1
ysc=-0.13
endif
*
if (rows=3 & cols=2)
xoffl=0.42
*xoffr=0.12
xoffr=0.3
yoffb=0.8
yofft=0.8
*xsc=0.12
xsc=0.3
*ysc=-0.17
ysc=-0.14
endif
*
if (rows=2 & cols=2)
xoffl=0.12
xoffr=0.12
yoffb=0.1
yofft=0.7
xsc=0.1
ysc=0.
endif
*
*xoff=0.1
*yoff=0.99
*xysc=0.03
if (rows=1 & cols=1) ; xoff=0 ; yoff=0 ; endif
*
'set vpage off'
'query gxinfo'
rec2 = sublin(result,2)
xlo = xoffl
xhi = subwrd(rec2,4)
xhi = xhi - xoffr
xde = xhi - xlo
ylo = yoffb
yhi = subwrd(rec2,6)
yhi = yhi - yofft
yde = yhi - ylo
*say xhi
*say yhi
*say xde
*say yde
xleft   = xlo + xde / cols*(col-1)
xright  = xlo + xde / cols*(col)
ybottom = yhi - yde/rows*(row)
ytop    = yhi - yde/rows*(row-1)
xleft   = xleft   - xsc
xright  = xright  + xsc
ybottom = ybottom - ysc
ytop    = ytop    + ysc
*if (xleft   < 0)   ; xleft=0    ; endif
*if (xright  > xhi) ; xright=xhi ; endif
*if (ybottom < 0)   ; ybottom=0  ; endif
*if (ytop    > yhi) ; ytop=yhi   ; endif
*say xleft
*say xright
*say ybottom
*say ytop
say 'set vpage '%xleft%' '%xright%' '%ybottom%' '%ytop
*'draw rec ' xleft ' ' ybottom ' ' xright ' ' ytop
'set vpage ' xleft ' ' xright ' ' ybottom ' ' ytop
'set grads off'
return

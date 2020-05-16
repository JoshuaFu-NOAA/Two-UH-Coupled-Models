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
if(rows='');say "syntax: run page rows cols gnum";return;endif
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
xoff=0
yoff=0.5
if (rows=1 & cols=1) ; xoff=0 ; yoff=0 ; endif
*
'set vpage off'
'query gxinfo'
rec2 = sublin(result,2)
xlo = xoff
xhi = subwrd(rec2,4)
xhi = xhi - xoff
xde = xhi - xlo
ylo = yoff
yhi = subwrd(rec2,6)
yhi = yhi - yoff
yde = yhi - ylo
*say xhi
*say yhi
*say xde
*say yde
xleft   = xde / cols*(col-1)
xright  = xde / cols*(col)
ybottom = yhi - yde/rows*(row)
ytop    = yhi - yde/rows*(row-1)
if (xleft   < 0)   ; xleft=0    ; endif
if (xright  > xhi) ; xright=xhi ; endif
if (ybottom < 0)   ; ybottom=0  ; endif
if (ytop    > yhi) ; ytop=yhi   ; endif
*say xleft
*say xright
*say ybottom
*say ytop
say 'set vpage '%xleft%' '%xright%' '%ybottom%' '%ytop
'set vpage ' xleft ' ' xright ' ' ybottom ' ' ytop
'set grads off'
return

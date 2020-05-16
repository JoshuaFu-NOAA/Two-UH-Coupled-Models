function DrawLine (args)

quiet=0

filename=subwrd(args,1)

if(filename='')
  say 'syntax: run DrawLine filename'
  return
endif

'q gxinfo'
fstat=sublin(result,1)
mproj=sublin(result,6)
*say mproj
*
xline=sublin(result,3)
yline=sublin(result,4)
*
xlo=subwrd(xline,4)
xhi=subwrd(xline,6)
*
ylo=subwrd(yline,4)
yhi=subwrd(yline,6)
*
*say xlo%'  '%xhi%'  '%ylo%'  '%yhi
'set clip 'xlo' 'xhi' 'ylo' 'yhi


line=read(filename)
fstat=sublin(line,1)
if (quiet > 0) ; say fstat ; endif
if (fstat = 0)
ianz=sublin(line,2)
ianz=subwrd(ianz,1)
if (quiet > 0) ; say ianz ; endif
endif

while(fstat = 0)
ip=1
if (quiet > 0) ; say ip; say ianz ; endif
while(ip<=ianz)
  line=read(filename)
  fstat=sublin(line,1)
  lin2=sublin(line,2)
  xval=subwrd(lin2,1)
  yval=subwrd(lin2,2)
  x.ip=xval
  y.ip=yval
*  say ip;
  ip=ip+1
endwhile
 
'q w2xy 'x.1' 'y.1
x1=subwrd(result,3)
y1=subwrd(result,6)
ip=2
while(ip<=ianz)
  'q w2xy 'x.ip' 'y.ip
  x2=subwrd(result,3)
  y2=subwrd(result,6)
*  dst='draw line '%x1%' '%y1%' '%x2%' '%y2
* if (ianz = 25) ;   say dst; endif
* say x1; say y1; say  x2 ; say y2
  'draw line 'x1' 'y1' 'x2' 'y2
*  say ip;
  x1=x2
  y1=y2
  ip=ip+1
endwhile

line=read(filename)
fstat=sublin(line,1)
ianz=sublin(line,2)
ianz=subwrd(ianz,1)
if (quiet > 0) ; say fstat ; endif
if (quiet > 0) ; say ianz ; endif

endwhile

return

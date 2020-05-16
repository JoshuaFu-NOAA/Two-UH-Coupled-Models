function DrawBox (args)
*
* run DrawBox xlo ylo xhi yhi
*
bxlo=subwrd(args,1)
bylo=subwrd(args,2)
bxhi=subwrd(args,3)
byhi=subwrd(args,4)
*
if (bxlo='' | bylo='' | bxhi='' | byhi='')
  say "syntax: run DrawBox xlo ylo xhi yhi"
  return
endif
*
quiet=0
*
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

 
'q w2xy 'bxlo' 'bylo
x1=subwrd(result,3)
y1=subwrd(result,6)
'q w2xy 'bxhi' 'byhi
x2=subwrd(result,3)
y2=subwrd(result,6)


*say x1; say y1; say  x2 ; say y2
'draw rec 'x1' 'y1' 'x2' 'y2

return

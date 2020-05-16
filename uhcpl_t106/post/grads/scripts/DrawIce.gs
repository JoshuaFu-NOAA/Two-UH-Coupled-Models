function DrawIce (args)
quiet=0

code=subwrd(args,1)

if(code='')
  say 'syntax: run DrawIce expr'
  return
endif (args)

'q dims'
fstat=sublin(result,1)
* check fstat No files open !!!

xline=sublin(result,2)
yline=sublin(result,3)
xva=subwrd(yline,13)

'set mpdraw off'
'set grads off'

'define tmean=aave('code',x=1,x=1,y='xva',y='xva')'
'd tmean'
xmean=subwrd(result,4)
say code%'(1,'xva')='xmean

sfact=1
xval=xmean
while (xval < 10000)
 xval=xmean*sfact
 sfact=sfact*10
endwhile

'define SFACT='sfact
say 'fact='%sfact

xmean=xmean*sfact
* say xmean

'define EPS=0.005'

'define tmax=(-tmean-EPS)*SFACT'
'define tmin=(-tmean+EPS)*SFACT'
'd tmin'
x=sublin(result,1)
* say x
'd tmax'
x=sublin(result,1)
* say x
'define temp=('code')*SFACT'
'define c=maskout(maskout(temp,temp+tmin),-temp-tmax)'
'set gxout fgrid'
'set fgvals 'xmean' 4'
'd c'

return

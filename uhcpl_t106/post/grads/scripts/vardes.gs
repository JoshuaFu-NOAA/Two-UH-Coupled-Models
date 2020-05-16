function vardes (args)
*
* get the VARiable DESription for var
*
* run vardes var <part>
*
var=subwrd(args,1)
*
if (var='')
  say "usages: run vardes var <part> <filen>"
  return
endif
*
part=subwrd(args,2)
if (part='')
  part=0
endif
filen=subwrd(args,3)
if (filen='')
  filen=1
endif
*
'q file 'filen
res   = result
dum   = sublin(res,1)
title = subwrd(dum,4)' 'subwrd(dum,5)' 'subwrd(dum,6)' 'subwrd(dum,7)' 'subwrd(dum,8)
title = title' 'subwrd(dum,9)' 'subwrd(dum,10)' 'subwrd(dum,11)' 'subwrd(dum,12)
dum   = sublin(res,5)
xsize = subwrd(dum,3)
ysize = subwrd(dum,6)
znum  = subwrd(dum,9)
tnum  = subwrd(dum,12)
dum   = sublin(res,6)
vnum  = subwrd(dum,5)
*
*say ' TITLE      : 'title
*say ' NUMBER OF VALUES IN X-DIRECTION : 'xsize
*say ' NUMBER OF VALUES IN Y-DIRECTION : 'ysize
*say ' NUMBER OF LEVELS                : 'znum
*say ' NUMBER OF TIMESTEPS             : 'tnum
*say ' NUMBER OF VARIABLES             : 'vnum
*
vdes='???'
i = 0
while (i<vnum)
  dum  = sublin(res,i+7)
  varn = subwrd(dum,1)
  varc = subwrd(dum,3)
  if ( var = varn | var = code%varc | var = varc )
    vdes=''
    j=0
    while (j<10)
      vh = subwrd(dum,4+j)
      if (part = 1)
        fchar = substr(vh,1,1)
        if (fchar = '[') ; break ; endif
      endif
      if (vh = '') ; break ; endif
      vdes = vdes' 'vh
      j = j + 1
    endwhile
    break
  endif
  i = i + 1
endwhile
*
return vdes

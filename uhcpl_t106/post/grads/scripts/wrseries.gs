function wrseries(arg)
* write out a time series in ASCII.
* USAGE: run wrseries filename fld ts te [-q]
* where 
*   fld:   grads expressions
*    ts:   first time step
*    te:   last time step
*    -q:   option supresses output of the data to screen.
*
* NOTE: The space units must be set to one data point!
*
* Example: run wrseries datafile.txt zeta 1 100
file=subwrd(arg,1)
fld=subwrd(arg,2)
ts=subwrd(arg,3)
te=subwrd(arg,4)
vrbos=subwrd(arg,5)
if(substr(vrbos,1,2)!= '-q')
   vrbos=v
endif
if(vrbos='')
   vrbos=v
endif
say '     writing 'fld' for 'ts'<=t<='te' to 'file'.'

tt=ts
'!rm -f 'file
while (tt <= te)
   'set t 'tt
   'd 'fld
   yy=subwrd(result,4)
   if(vrbos='v')
      say '  y('tt')='yy
   endif
   good=write(file,tt' 'yy,append)
   tt=tt + 1
endwhile
'set t 'ts
return






function saveshd (args)
*
* save shading information to a file
*
* run saveshd filename
*
* externals: none
*
shdfile=subwrd(args,1)
*
if (shdfile='')
  say 'usages: run saveshd filename'
  return
endif
*say shdfile
*
'query shades'
shdinfo = result
if (subwrd(shdinfo,1)='None') 
   say 'Cannot save to file: No shading information'
   return
endif
*
ret=write(shdfile,shdinfo)
*
sret=sublin(ret,1)
if (sret=1)
  say 'Cannot open '%shdfile
  return
endif
*
return

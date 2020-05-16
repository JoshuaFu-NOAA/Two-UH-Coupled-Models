function signatur (args)
*
* make a signature on the bottom left corner of the page
*
* run signatur <scale> <rotation> <username> <date>
*
* externals: user.gs date.gs page.gs
*
scale=subwrd(args,1)
rotation=subwrd(args,2)
user=subwrd(args,3)
date=subwrd(args,4)
*
if (user='')
  'run user.gs'
  user=sublin(result,1)
endif
*
if (date='')
  'run date.gs'
  date=sublin(result,1)
endif
*
if (scale='')
  scale=1.0
endif
*
if (rotation='')
  rotation=0
endif
*
*say "user="%user
*say "date="date
*say "scale="scale
*say "rotation="rotation
size=0.07*scale
*
'run page.gs 1 1 1 1'
'set font 1'
'set string 1 bl 1 0'
if (rotation>0)
  'set string 1 tl 1 90'
endif
'set strsiz 'size
'draw string 0.0 0.0 'user'  'date
*
return

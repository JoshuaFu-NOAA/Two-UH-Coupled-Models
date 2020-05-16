function pannel (args)
pmax=subwrd(args,1)
pcur=subwrd(args,2)
*
*say pmax
*say pcur
if (pcur = 1)
'set line 1 1 7'
'draw line 0.5 8 10.5 8'
'draw line 0.5 8.3 10.5 8.3'
'draw line 0.5 8 0.5 8.3'
'draw line 10.5 8 10.5 8.3'
endif
'set line 3 1 1'
xpos = 0.5 + 10. * pcur / pmax
xpos0=0.5
*xpos=1.5
*say xpos
*'draw line 'xpos0' 8 'xpos' 8.3'
'draw recf 'xpos0' 8 'xpos' 8.3'
return

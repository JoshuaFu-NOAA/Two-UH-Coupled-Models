spawn, 'gv idl.ps &'
;------------------------------
; default plotting options

ystyle=1
xstyle=0
thick=1.
xthick=1.
ythick=1.
title=' ' 
xtitle=' '
ytitle=' '
charthick=1. 
charsize=1.
index=3
;------------------------------

max=9999
nitems=5
yy=fltarr(nitems,max)
x=fltarr(nitems,max)

xr=0.0
yr=fltarr(nitems)
i=0
ans='1'

repeat begin
  if (ans eq '1') then spawn, 'post_1d | gnuplot' $
                  else spawn, 'post_1d'
  openr, 1, 'gnu.plot'
  j=0
  print, 'reading data'
  while (eof(1) ne 1) do begin
    readf, 1, xr,yr 
;    print, xr, yr
    x(j)=xr
    yy(*,j)=yr
    j=j+1
  endwhile
  close, 1
  i=i+1

  y=yy(index,*)

; cut arrays off to correct size
  x=x(0:j-1)
  y=y(0:j-1)

; find max/min
  ymin=MIN(y)
  ymax=MAX(y)
  xmin=MIN(x)
  xmax=MAX(x)

  exit=0

  set_plot, 'PS'
  !p.font=0
  device, xsize=18, ysize=18, xoffset=2, yoffset=5, /portrait, /helvetica, $
  filename='idl.ps'

  plot, x, y, $
     thick=thick, ystyle=ystyle, xstyle=xstyle, $
     xthick=xthick, ythick=ythick, $
     xtitle=xtitle, ytitle=ytitle, $
     yrange=[ymin,ymax],xrange=[xmin,xmax], $
     charthick=charthick, charsize=charsize 

  device, /close
  set_plot, 'X'

  repeat begin
    print
    print, ' Note some options require axis specification (marked *)'
    print, ' enter x or y in lower case, EXAMPLE: "3x" '
    print
    print, ' option ?              current value'
    print
    print, ' 0: exit '
    print, ' 1: swap x and y axis'
    print, ' 2: *reverse axis x/y'
    print, ' 3: *scale x/y'
    print
    print, ' 5: thick            ',thick
    print, ' 6: *axis thick x/y  ',xthick,ythick
    print, ' 8: character thick  ',charthick
    print, ' 9: character size   ',charsize
    print
    print, '10: Title             ',title
    print, '11: *axis Title   x/y ',xtitle,ytitle
    print
    print, '13: specify xmin/max ',xmin,xmax
    print, '14: specify ymin/max ',ymin,ymax
    print, '15: xy min/max < file  '
    print, '16: *axis style   x/y ',xstyle,ystyle
    print
    print, '99: print to psm9'

    option = 'xxx'
    read, option
    case option of 
    '0': print, '0: exit, 1: read with GNUPLOT, 2:simple read' 
    '1': begin
         store=x
         x=y
         y=store
         store=xmin
         xmin=ymin
         ymin=store
         store=xmax
         xmax=ymax
         ymax=store
         end
   '2x': begin
         store = xmax
         xmax  = xmin 
         xmin  = store
         end
   '2y': begin
         store = ymax
         ymax  = ymin 
         ymin  = store
         end
   '3x': begin
         read, 'scale factor (negative = 1./scale) ',scale
         if (scale lt 0.) then scale = 1./ABS(scale)
         x=x*scale
         xmin=xmin*scale
         xmax=xmax*scale
         end
   '3y': begin
         read, 'scale factor (negative = 1./scale) ',scale
         if (scale lt 0.) then scale = 1./ABS(scale)
         y=y*scale
         ymin=ymin*scale
         ymax=ymax*scale
         end
    '5': read, 'thick ',thick
   '6x': read, 'xthick ',xthick
   '6y': read, 'ythick ',ythick
    '8': read, 'char thick ',charthick
    '9': read, 'charsize ',charsize
   '10': read, 'title string (!U superscript !D subscript !N normal) ',title
  '11x': read, 'xtitle string (!U superscript !D subscript !N normal) ',xtitle
  '11y': read, 'ytitle string (!U superscript !D subscript !N normal) ',ytitle
   '13': begin
         read, 'x min ',xmin
         read, 'x max ',xmax
         end
   '14': begin
         read, 'y min ',ymin
         read, 'y max ',ymax
         end
   '15': begin
         xmin=MIN(x)
         xmax=MAX(x)
         ymin=MIN(y)
         ymax=MAX(y)
         end
  '16x': read, 'x style (0 round, 1 exact range) ',xstyle
  '16y': read, 'y style (0 round, 1 exact range) ',ystyle 
   '99': spawn, 'lp -d psm9 idl.ps'
    else: PRINT, '*** illegal value ***'
    endcase

    set_plot, 'PS'
    !p.font=0
    device, xsize=18, ysize=18, xoffset=2, yoffset=5, $
                  /portrait, /helvetica, filename='idl.ps'
    plot, x, y, $
     thick=thick, ystyle=ystyle, xstyle=xstyle, $
     xthick=xthick, ythick=ythick, $
     xtitle=xtitle, ytitle=ytitle, title=title, $
     yrange=[ymin,ymax],xrange=[xmin,xmax], $
     charthick=charthick, charsize=charsize 

    device, /close
    set_plot, 'X'
  endrep until option eq '0'

  read, ans
endrep until ans eq '0'
END
function qsaturation, t, p

; replace function for qsat with tetons formula
; t=temp in K
; p=pressure in hPa

qsat=(t gt 273.15)*3.8/(p*EXP(-17.2693882*(t-273.15)/(t-35.86))-6.109) + $
     (t le 273.15)*3.8/(p*EXP(-21.8745584*(t-273.15)/(t-7.66))-6.109) 
return, qsat
end

pro tplot,dat,ilen,ctit,plotmoi,sd_plot
forward_function readdat,f87,svp,dew,rmix,potwet,entrop,wet,alatnt
;
;  tplot plots the temp and dew point lines and the 
;  titles on ht etephigram
;  dat: the data to be plotted ( p(mb), t(c), td(c) )
;  ilen: number of points in dat
;  ctit: title to be plotted at top of tephigram
;
; plotmoi - moisture plot yes/no
; sd_plot - shading for stan dev

common corner, t0, h0, t1, h1, t2, h2

x=fltarr(1000)
y=fltarr(1000)
;
akappa = 0.286
;
;  the tephigram is drawn in a 9x9 square
      xw=9.0 & yw=9.0
;
      ct=t0
      bt=(t2-ct)/yw
      at=(t1-ct)/xw
      ch=h0
      bh=(h2-ch)/yw
      ah=(h1-ch)/xw
      divisx=at*bh-ah*bt 
      divisy=bt*ah - bh*at
      ax=bh/divisx
      bx=-bt/divisx
      cx=(ch*bt-ct*bh)/divisx
      ay=ah/divisy
      by=-at/divisy
      cy=(ch*at-ct*ah)/divisy

      p0=1000.*(((273.15+t0)/(273.15+h0))^(1./akappa))
      p1=1000.*(((273.15+t1)/(273.15+h1))^(1./akappa))
      p2=1000.*(((273.15+t2)/(273.15+h2))^(1./akappa))
;
;  set character height to 0.20, align to bottom and write title
;
lab1='' & for n=0,34 do lab1=lab1+ctit(1,n)
lab2='' & for n=0,34 do lab2=lab2+ctit(2,n)
lab3='' & for n=0,34 do lab3=lab3+ctit(0,n)
xyouts, CHARTHICK=3,0.1,yw+0.1,lab1,/data
xyouts, CHARTHICK=3,0.1,yw+0.4,lab2,/data
xyouts, CHARTHICK=3,0.1,yw+0.7,lab3,/data

;*AMT* 
; ADD LABELS

if (1 eq 1) then begin

device, /HELVETICA, /BOLD
!p.noclip=1

;print,'window=',!d.x_size,!d.y_size
xs=!d.x_size
ys=!d.y_size

xyouts, charthick=3,.36*xs,.08*ys,'Water Vapour Mixing Ratio (g kg!U-1!n)', /device,alignment=0.5, size=1.
xyouts, charthick=3,.08*xs,!d.y_size/2,'Pressure (hPa)', /device,alignment=0.5, orientation=90
xyouts, charthick=3,.92*xs,!d.y_size/2,'Temperature  (!Uo!nC)', /device,alignment=0.5,orientation=270

xpi=xs*0.71 ;***TOP RIGHT***
xpi=xs*0.77 ;***TOP RIGHT***
;xpi=xs*0.36 ;***TOP LEFT***
ypi=.87*ys

ww=4000
polyfill,[xpi-ww/2,xpi+ww/2,xpi+ww/2,xpi-ww/2],[ypi,ypi,ypi+400,ypi+400],color=-1, /device
xyouts, charthick=3,xpi,ypi+100,'Temperature (!Uo!nC)', /device,alignment=0.5

ypi=.818*ys
ww=6000
polyfill,[xpi-ww/2,xpi+ww/2,xpi+ww/2,xpi-ww/2],[ypi,ypi,ypi+400,ypi+400],color=-1, /device
xyouts, charthick=3,xpi-700*1,ypi+00,'Potential Temperature (!Uo!nC)', /device,alignment=0.5

xpi=0.185*xs
ypi=!d.y_size/2
ww=5500
polyfill,[xpi-300,xpi+80,xpi+80,xpi-300],[ypi-ww/2,ypi-ww/2,ypi+ww/2,ypi+ww/2],color=-1, /device
xyouts, charthick=3,xpi,ypi+100,'Potential Temperature (!Uo!nC)', /device,alignment=0.5, orientation=90

DEVICE, /HELVETICA
!p.font=0
endif ;*AMT*

; turn clipping on now title have been plotted
!p.noclip=0
;
;  create t points
;
      for i=0,ilen-1 do begin
       t=dat(1,i)
       p=dat(0,i)
;print, '*',i,t,p
       h=(273.15+t)*(1000./p)^akappa - 273.15
       x(i)=ax*t+bx*h+cx
       y(i)=ay*t+by*h+cy
      endfor
!p.thick=10.0
;print
;print,'TEMPERATURE:'
;read,'Colour index (0-255)(0=black,2=red,3=green,4=blue,5=yellow): ',col
;read,'linestyle (0=solid,1=dot,2=dash,3=d/d,5=longdash): ',lin
;read,'line thickness: ',thick
col=0
lin=0
thick=10

if (sd_plot eq 1) then polyfill,x(0:ilen-1),y(0:ilen-1), noclip=0,color=col else $
plots,x(0:ilen-1),y(0:ilen-1),$
    color=col,noclip=0,linestyle=lin,thick=thick,/data

;
;  repeat for dew points (but some - imiss - are missing)
;
      imiss=0
      for i=0,ilen-1 do begin
       t=dat(2,i)
       if(t eq -999.) then begin
        imiss=imiss+1
       endif else begin
        p=dat(0,i)
        h=(273.15+t)*(1000./p)^akappa - 273.15
        x(i-imiss)=ax*t+bx*h+cx
        y(i-imiss)=ay*t+by*h+cy
       endelse
      endfor
;!p.linestyle=5

if (plotmoi ne 'n') then begin
;print
;print,'MOISTURE:'
;read,'Colour index (0-255)(0=black,2=red,3=green,4=blue,5=yellow): ',col
;read,'linestyle (0=solid,1=dot,2=dash,3=d/d,5=longdash): ',lin
;read,'line thickness: ',thick
col=2
lin=5
thick=10

if (sd_plot eq 1) then polyfill,x(0:ilen-1),y(0:ilen-1), noclip=0,color=col else $
plots,x(0:ilen-imiss-1),y(0:ilen-imiss-1),color=col,thick=thick,linestyle=lin,$
noclip=0,/data
endif

;!p.linestyle=0
!p.thick=1.0
;
end
;
pro teph
; teph plots the tephigram structure, i.e
; 1) Pressure, Temp and Hmr lines
; 2) Wet and Dry adiabats
; 3) zero isotherm

common corner, t0, h0, t1, h1, t2, h2

akappa = 0.286
xw=9.0 & yw=9.0
iprec =100
;nhmr = 18
;nhmr=14


xp=fltarr(iprec)
yp=fltarr(iprec)
; following are dimensioned 1 more than in corresponing
; fortran so don't have to change the subscripts
x2=fltarr(3)
y2=fltarr(3)
xs2=fltarr(3)
ys2=fltarr(3)
xout=fltarr(6)
yout=fltarr(6)

       ct=t0
       bt=(t2-ct)/yw
       at=(t1-ct)/xw
       ch=h0
       bh=(h2-ch)/yw
       ah=(h1-ch)/xw
       divisx=at*bh-ah*bt
       divisy=bt*ah-bh*at
       ax=bh/divisx
       bx=-bt/divisx
       cx=(ch*bt-ct*bh)/divisx
       ay=ah/divisy
       by=-at/divisy
       cy=(ch*at-ct*ah)/divisy

      xout=[0.0,0.0,xw,xw,0.0,0.0]
      yout=[0.0,0.0,0.0,yw,yw,0.0]
;
;  the following data lines define the beginning and end of hmr lines
;  since they are very difficult to calculate analytically
;
;  alterations should be made as follows:  alter nhmr to be the
;  number of lines to be plotted : hms holds the values used for
;  labelling; tt1s and hh1s are the (t,theta) co-ordinates of
;  the bottom left of the line (n.b. they don't have to lie
;  within the actual area of the plotted tephigram) and tt2s,hh2s
;  the top right
;
; **** original hmr values **** set nhmr to 14 (AMT) ****
      hms=[0.02,0.05,0.15,0.3,0.6,0.8,1.5,2.,4.,6.,8.,10.,14.0,20.] 
;      else hms= [0.02,0.05,0.2,0.4,0.6,0.8,1.0,1.5,2.,3.,4.,6.,8., $
;              10.,14.,16.,18.,20.]
      nhmr=n_elements(hms)

; - NEW CODE TO INTERATIVELY FIND THE MOISTURE LINES USING SAME 
; - QSAT FUNCTION AS CRM (function for qsat in qsaturation.pro)
; - USES SIMPLE BISECTION METHOD (Numerical Recipes FORTRAN p.346)
  tt1s=fltarr(nhmr)
  hh1s=fltarr(nhmr)
  tt2s=fltarr(nhmr)
  hh2s=fltarr(nhmr)

  press=1050.
  for i=0,nhmr-1 do begin
    tit1=150.
    tit2=350.
    qg=0.
    dt=tit2-tit1
    while (abs(qg-hms(i)) gt 0.0001*abs(hms(i))) do begin
       dt=dt*0.5
       tmid=tit1+dt
       qg=1000.*qsaturation(tmid,press)
       if (qg lt hms[i]) then tit1=tmid
       ct=ct+1
       if (ct eq 1000) then stop
    endwhile
    tt1s[i]=tmid-273.15
    hh1s[i]=(tmid*(1000./press)^akappa)-273.15
  endfor

; have to repeat due to method of using 2 arrays in
; original program (yuk)

  press=300.
  for i=0,nhmr-1 do begin
    tit1=150.
    tit2=350.
    qg=0.
    dt=tit2-tit1
    while (abs(qg-hms(i)) gt 0.0001*abs(hms(i))) do begin
       dt=dt*0.5
       tmid=tit1+dt
       qg=1000.*qsaturation(tmid,press)
       if (qg lt hms[i]) then tit1=tmid
       ct=ct+1
       if (ct eq 1000) then stop
    endwhile
    tt2s[i]=tmid-273.15
    hh2s[i]=(tmid*(1000./press)^0.285)-273.15
  endfor

;      tt1s=[ -58.,-49.,-42.,-32.5,-22.5,-19.,-11.5,-8.,1.2, $
;            7.,11.2,14.5,20.0,25.5]
;      hh1s=[ 1.3,-8.2,-10.,-20.,-26.7,-22.7,-15.2,-12., $
;            -2.7,3.1,7.1,10.2,16.30,21.3]
;      tt2s=[ -64.,-57.,-49.5,-46.,-38.,-35.2,-25.5,-20., $
;            -10.2,-5.,-1.,2.5,10.0,13.0]
;      hh2s=[ 62.5,71.3,48.6,84.6,80.,86.1,70.,58.0,60.,66., $
;            70.,70.,60.0,80.]

; **** new hmr values **** set nhmr to 18 ****
;     hms= [0.02,0.05,0.2,0.4,0.6,0.8,1.0,1.5,2.,3.,4.,6.,8., $
;              10.,14.,16.,18.,20.]

      p0=1000.0*(((273.15+t0)/(273.15+h0))^(1.0/akappa))
      p1=1000.0*(((273.15+t1)/(273.15+h1))^(1.0/akappa))
      p2=1000.0*(((273.15+t2)/(273.15+h2))^(1.0/akappa))
;
charlen=0.2  ; approximate length of a character

; draw border of tephigram
;
plots,xout(1:5),yout(1:5),/data
!p.thick=1.0
;
;**********************************************************************
; *AMT* first lot of lines up to 400mb
;**********************************************************************
; draw pressure lines with precision iprec from 1050mb to 300mb
;
      for ip=1050,450.,-100 do begin
;      for ip=1000,p2,-100 do begin
       press=float(ip)
       if(ip mod 50 eq 0) then begin
       endif else begin
       endelse
       tmin=t0+(t2-t0)*(press-p0)/(p2-p0)
       for ilp=0,iprec-1 do begin
        t=tmin+float(ilp)*(t1+15.0-t0)/float(iprec-1)
        h=(273.15+t)*(1000.0/press)^akappa - 273.15
        xp(ilp)=ax*t+bx*h+cx
        yp(ilp)=ay*t+by*h+cy
       endfor
       if(ip mod 50 eq 0) then begin
;
clab=string(ip,format='(i4)')
;
;  find first real y point
;
        for ilp=0,iprec-1 do begin
         if(xp(ilp) gt 0.0) then begin
          ytxt=yp(ilp)
          goto,lab201
         endif
        endfor
lab201:
!p.noclip=0
xyouts, CHARTHICK=3,0.03,ytxt,clab,noclip=0,/data
;
;  don't overwrite label
;
        xnew=0.03+4.0*charlen
        for ilp=0,iprec-1 do begin
         if(xp(ilp) gt xnew) then begin
          if(ilp ne 0) then begin
           xp(ilp-1)=xnew
           yp(ilp-1)=yp(ilp)
           istart=ilp-2
          endif else begin
           istart=ilp-1
          endelse
          goto,lab101
         endif
        endfor
       endif
;      give labels a bit more room by increasing starting index by 2
;

lab101:
plots,xp(istart:iprec-1),yp(istart:iprec-1),noclip=0,/data
      endfor
;**********************************************************************
; *AMT* sdcond lot of pressure lines!!! from 400mb
;**********************************************************************
;
; draw pressure lines with precision iprec from 1050mb to 300mb
;
      for ip=400,p2,-50 do begin
;      for ip=1000,p2,-100 do begin
       press=float(ip)
       if(ip mod 50 eq 0) then begin
       endif else begin
       endelse
       tmin=t0+(t2-t0)*(press-p0)/(p2-p0)
       for ilp=0,iprec-1 do begin
        t=tmin+float(ilp)*(t1+15.0-t0)/float(iprec-1)
        h=(273.15+t)*(1000.0/press)^akappa - 273.15
        xp(ilp)=ax*t+bx*h+cx
        yp(ilp)=ay*t+by*h+cy
       endfor
       if(ip mod 50 eq 0) then begin
;

clab=string(ip,format='(i4)')
;
;  find first real y point
;
        for ilp=0,iprec-1 do begin
         if(xp(ilp) gt 0.0) then begin
          ytxt=yp(ilp)
          goto,lab202
         endif
        endfor
lab202:
!p.noclip=0
xyouts, CHARTHICK=3,0.03,ytxt,clab,noclip=0,/data
;
;  don't overwrite label
;
        xnew=0.03+4.0*charlen
        for ilp=0,iprec-1 do begin
         if(xp(ilp) gt xnew) then begin
          if(ilp ne 0) then begin
           xp(ilp-1)=xnew
           yp(ilp-1)=yp(ilp)
           istart=ilp-2
          endif else begin
           istart=ilp-1
          endelse
          goto,lab102
         endif
        endfor
       endif
;      give labels a bit more room by increasing starting index by 2
;
lab102:
plots,xp(istart:iprec-1),yp(istart:iprec-1),noclip=0,/data
      endfor
;
; draw temperature lines
;
!p.thick=2.0
!p.linestyle=0
lower=0.05
forward=0.05
      tmin=float(round(t2/10.))
      tmin=10.*tmin
      for it=tmin,t1,5 do begin
;      for it=tmin,t1,10 do begin
       t=float(it)
       h=(273.15+t)*(1000.0/1100.0)^akappa - 273.15
       x2(1)=ax*t+bx*h+cx
       y2(1)=ay*t+by*h+cy
       h=(273.15+t)*(1000.0/100.0)^akappa - 273.15
       x2(2)=ax*t+bx*h+cx
       y2(2)=ay*t+by*h+cy
;
;  put in label at 0.92*yw or 0.1*xw at 10 degree intervals
;
       if(it mod 5 eq 0) then  begin
clab=string(it,format='(i4)')
       agrd=(y2(2)-y2(1))/(x2(2)-x2(1))
       ac=y2(1)-agrd*x2(1)
;  character up vector will be along -agrd,1
       theta=atan(agrd)*180.0/3.1415
       if(agrd lt 0.0) then theta=180.0-theta
       ytx=0.92*yw
       xtx=(ytx-ac)/agrd
       xs2(1)=x2(1)
       ys2(1)=y2(1)
       if(xtx lt xw) then begin
        xs2(2)=xtx
        ys2(2)=ytx
       endif else begin
        xs2(2)=0.93*xw
        ys2(2)=agrd*xs2(2)+ac
       endelse
       angle=atan(agrd)
       x2(1)=xs2(2)+4.0*charlen*cos(angle)
       y2(1)=ys2(2)+4.0*charlen*sin(angle)
!p.noclip=0
; decrease height of labels as can't be vertically centred in idl
; and move forward slightly to improve appearance
if (xs2(2)+forward gt 1. and xs2(2)+forward lt xw and ys2(2)-lower gt 1. and ys2(2)-lower lt yw) then $
xyouts, CHARTHICK=3,xs2(2)+forward,ys2(2)-lower,clab,noclip=1,/data,orientation=theta

;stop

!p.linestyle=0
!p.thick=2.0
       if(t eq 0.0) then begin
; plot dashed line for zero isotherm
!p.linestyle=2
       endif
plots,xs2(1:2),ys2(1:2),noclip=0,/data
       endif
       if(t eq 0.0) then begin
       endif
plots,x2(1:2),y2(1:2),noclip=0,/data
!p.linestyle=0
!p.thick=1.0
      endfor

; draw pot. temperature lines (Dry Adiabats)
;
!p.linestyle=0
!p.thick=1.0
      hmin=float(round(h0/10.0))
      hmin=10.0*hmin
;      hmin=float(round(h0/5.0))
;      hmin=5.0*hmin
;      for ih=hmin,h2+h1-h0,5 do begin
      for ih=hmin,h2+h1-h0,10 do begin
       h=float(ih)
       t=(h+273.15)/(1000.0/100.0)^akappa - 273.15
       x2(1)=ax*t+bx*h+cx
       y2(1)=ay*t+by*h+cy
       t=(h+273.15)/(1000.0/1100.0)^akappa - 273.15
       x2(2)=ax*t+bx*h+cx
       y2(2)=ay*t+by*h+cy
;
;  put in label at 0.97*yw or 0.1*xw at 10 deg intervals
;
       if(ih mod 10 eq 0) then begin
clab=string(ih,format='(i4)')
       agrd=(y2(2)-y2(1))/(x2(2)-x2(1))
       ac=y2(1)-agrd*x2(1)
;  character up vector will be along -agrd,1
       theta=atan(agrd)*180.0/3.1415
       if(theta lt 0.0) then theta=360.0+theta
; ytx was .97*yw
       ytx=0.90*yw
       xtx=(ytx-ac)/agrd
       xs2(1)=x2(1)
       ys2(1)=y2(1)
       if(xtx gt 0.0) then begin
        xs2(2)=xtx
        ys2(2)=ytx
       endif else begin
        xs2(2)=0.10*xw
        ys2(2)=agrd*xs2(2)+ac
       endelse
       angle=atan(agrd)
       x2(1)=xs2(2)+4.*charlen*cos(angle)
       y2(1)=ys2(2)+4.*charlen*sin(angle)
       if(y2(2) gt 0.0) then begin
        y2(2)=0.
        x2(2)=-ac/agrd
       endif
!p.linestyle=0
!p.noclip=0
if (xs2(2)+forward gt .05 and xs2(2)+forward lt xw-.5 and ys2(2)-lower gt 1. and ys2(2)-lower lt yw) then $
xyouts, CHARTHICK=3,xs2(2)-forward+.1,ys2(2)-lower-.1,clab,noclip=0,/data,orientation=theta
!p.linestyle=2
plots,xs2(1:2),ys2(1:2),noclip=0,/data
       endif
plots,x2(1:2),y2(1:2),noclip=0,/data
      endfor
;
;  draw wet adiabats to -50 degrees
;
!p.linestyle=0
!p.thick=1.0
      tmin=float(round(t0/5.))
      tmin=5*tmin
      for it=tmin,t1,5 do begin
       t=float(it)
       told=-20.
       p=1000.
       pinc=(100.0-1000.0)/float(iprec-1)
       h=(273.15+t)*(1000.0/p)^akappa - 273.15
       xp(0)=ax*t+bx*h+cx
       yp(0)=ay*t+by*h+cy
       for ilp=1,iprec-1 do begin
        pav=0.5*(p+pinc+p)
        p=p+pinc
        tk=t+273.15
        e=svp(t)
        al=alatnt(t)
        v1=al*e/(pav*461.51)
        v2=1005.0*pav/287.05
        v3=al*al*0.622*e/(tk*tk*287.05*461.51)
        t=t+pinc*(tk+v1)/(v2+v3)
        h=(273.15+t)*(1000.0/p)^akappa - 273.15
        xp(ilp)=ax*t+bx*h+cx
        yp(ilp)=ay*t+by*h+cy
; don't draw portion of line where t < -50
        if(t lt -50.0) then begin
         scal=(-50.0-told)/(t-told)
         xp(ilp)=xp(ilp-1)+scal*(xp(ilp)-xp(ilp-1))
         yp(ilp)=yp(ilp-1)+scal*(yp(ilp)-yp(ilp-1))
         ilen=ilp
         goto,lab8
        endif
        told=t
       endfor
       ilen=iprec
lab8:
plots,xp(0:ilen-1),yp(0:ilen-1),noclip=0,/data
      endfor

!p.linestyle=0
!p.thick=1.0
;
;  draw hmr's
;
!p.linestyle=1
!p.thick=1.0
      for ih=2,nhmr-1 do begin
       iflag=0
       hm=hms(ih)/1000.0
       tt1=tt1s(ih)
       hh1=hh1s(ih)
       tt2=tt2s(ih)
       hh2=hh2s(ih)
       x2(1)=ax*tt1+bx*hh1+cx
       y2(1)=ay*tt1+by*hh1+cy
       x2(2)=ax*tt2+bx*hh2+cx
       y2(2)=ay*tt2+by*hh2+cy
;
;  put in label at 0.80*yw or 0.85*xw
;
       if(hms(ih) lt 1.0) then begin
clab=string(hms(ih),format='(f4.2)')
       endif else begin
clab=string(hms(ih),format='(f4.1)')
       endelse
;      find equation of line i.e y=mx+c where m=agrd,c=ac
 	
        agrd=(y2(2)-y2(1))/(x2(2)-x2(1))
	
        ac=y2(1)-agrd*x2(1)
;  character up vector will be along -agrd,1
       theta=atan(agrd)*180.0/3.1415
       if(agrd lt 0.0) then theta=180.0-theta
;
; check lines reach boundaries
;
       x2(1)=-1.0
       y2(1)=agrd*x2(1)+ac
       x2(2)=xw+1.0
       y2(2)=agrd*x2(2)+ac
       ytx=0.02*yw
       xtx=(ytx-ac)/agrd
       xs2(1)=x2(1)
       ys2(1)=y2(1)
       if(xtx lt xw) then begin
         xs2(2)=xtx
         ys2(2)=ytx
       endif else begin
         xs2(2)=0.85*xw
         ys2(2)=agrd*xs2(2)+ac
         iflag=1
       endelse
       angle=atan(agrd)
       x2(1)=xs2(2)+4.0*charlen*cos(angle)
       y2(1)=ys2(2)+4.0*charlen*sin(angle)
!p.linestyle=0
!p.thick=1.0
!p.noclip=0
; decrease height of labels as can't be vertically centred in idl
if (xs2(2)+forward-.1 gt .6 and xs2(2)+forward lt xw ) then $
xyouts, CHARTHICK=3,xs2(2)+forward-.1,ys2(2)-lower-.15,clab,noclip=0,/data,orientation=theta else print,'missing label is ',clab

;.....modify end points so on boundaries of plot
;
       if(ac lt 0.0) then begin
         ys2(1)=0.0
         xs2(1)=-ac/agrd
       endif else begin
         xs2(1)=0.0
         ys2(1)=ac
       endelse
!p.linestyle=1
!p.thick=1.0
plots,xs2(1:2),ys2(1:2),noclip=0,/data
	   
       if(iflag eq 1) then begin
         x2(2)=xw
         y2(2)=agrd*xw+ac
       endif
plots,x2(1:2),y2(1:2),noclip=0,/data
	  
!p.linestyle=0
!p.thick=1.0
     endfor
end

function readdat,inunit,header,np,dat
;
; to read standard data format
;
;.....initialise variables
;
maxpt = 1000
;
;.....read in header 
;
;;;count=-1 *amt*
count=np

head=strarr(35)
on_ioerror,endpro
;readf,inunit,format='(35a1)',  head
;header(0,0:34)=head(0:34)
;readf,inunit,format='(35a1)',  head
;header(1,0:34)=head(0:34)
;readf,inunit,format='(35a1)',  head
;header(2,0:34)=head(0:34)

char=strarr(1)
readf,inunit,format='(a1)',  char
readf,inunit,format='(a1)',  char  ;*AMT*
readf,inunit,format='(a1)',  char  ;*AMT*
readf,inunit,format='(a1)',  char  ;*AMT*
;
;.....read in data
;
vals=fltarr(4)
for i=1,maxpt do begin
    readf,inunit,format='(4(f16.8,x))',vals     ;*AMT* changed format
    a=vals(0) & b=vals(1) & c=vals(2) & d=vals(3)
    if(a eq 9999.0 or a eq -9999.0) then goto,endpro
    if(c ne 999.0 and d ne 999.0) then begin
       count=count+1
       dat(0,count)=a
       dat(1,count)=c
       dat(2,count)=d
 print, a,b,c,d
    endif
endfor
endpro:
np=count+1
return,1
end

function f87,inunit,cline,m,dat
;
; to read fronts 87 format data for dropsonde descents.
;
;.....initialise variables
;
maxpt = 1000
m = -1
;
;.....read in header
;
line=strarr(76)
cline=strarr(160)
readf,inunit,format='(a76)',  line
cline(0:75)    = line(0:75)
readf,inunit,format='(a76)',  line
cline(76:151)  = line(0:75)
readf,inunit,format='(a76)',  line
cline(152:159) = line(0:8)
;
;.....read in data
;
on_ioerror,endpro
first:
readf,inunit,format='(i4,i6,6(i8),3(i4),i1,i2,3i1)', $
      ifile,ilat,ilon,itime,igeoht,iradht, $
      ipp,itt,iuu,idd,iff,iwm,iwq,ip,it,ih
m = m + 1
;
;.....check that data is available
;
if (ipp ne 99999999 and itt ne 99999999 and iuu ne 9999) then begin
;
;.....convert data to real values
;
;     lat(m)   = ilat   * 0.001
;     lon(m)   = ilon   * 0.001
;     time(m)  = itime  * 0.01
;     geoht(m) = igeoht * 0.1
;
      dat(1,m) = ipp    * 0.01
      dat(2,m) = itt    * 0.01
      dat(3,m) = iuu    * 0.1
;
;     dd(m)    = idd    * 0.1
;     ff(m)    = iff    * 0.1
;
;.....if data is missing decrement counter.
;
endif else begin
      m = m - 1
endelse
;
;.....loop back to read in the next line.
;
goto, first
;
;.....finish
;
endpro:
m=m-1
return,1
end

function svp,temp
;
; *AMT* changed assumption to above or below freezing
; change to use TETON's formula
; t=temperature in Kelvin

  t=temp+273.15       
  svp=(t gt 273.15)*6.109/EXP(-17.2693882*(t-273.15)/(t-35.86)) + $
      (t le 273.15)*6.209/EXP(-21.8745584*(t-273.15)/(t-7.66))

return, svp
end

function dew,t,rh
;
; uses rh=svp(dew)/svp(t)
;
;      if (rh gt 1.0) then rh=rh/100.0
      edew=rh*svp(t)
      test=t
      for i=0,-3,-1 do begin
       tinc=5.0*10^float(i)
testlab:
       test=test-tinc
       etest=svp(test)
       if(etest gt edew) then goto,testlab
       test=test+tinc
      endfor
      dew=test
return, dew
end
 
function rmix,temp,press,rhumy
;
; calculates the mixing ratio of water vapour
;
      e=rhumy * svp(temp)
      rmix=( e / (press - e)) * 0.622
return, rmix
end

function potwet,temp,rhumy,press
;
;  integrates along a moist adiabat to find wet bulb potential temp.
;
      cp=1005.0 & r=8.31 & rv=461.51 & ra=287.05 & inc=20
      tw=wet(temp,rhumy,press)
      tw=tw+273.15
      pinc=(1000.-press)/inc
      for i=1,inc do begin
       pav=0.5*(press+pinc+press)
       press=press+pinc
       e=svp(tw-273.15)
       al=alatnt(tw-273.15)
       v1=al*e/(press*rv)
       v2=cp*press/ra
       v3=al*al*0.622*e/(tw*tw*ra*rv)
       tw=tw+pinc*(tw+v1)/(v2+v3)
      endfor
      potwet=tw-273.15
return,potwet
end
 
function entrop,t,p
;
; calculates the entropy of a parcel of air just saturated
;
     r=287.05 & cp=1005.0 & c=4218.0
     e=svp(t)
     entrop=(cp+0.622*e/p)*alog(t+273.15) - r*alog(p-e) + alatnt(t) $
       *e*0.622/(t+273.15)/p
return, entrop
end
 
function wet,temp,rhumy,press
;
;   iterates to find the wet bulb temperature of a mass of air
;
      cp=1005.0 & cw=4218.0 & cpv=0.441*cw
      if(rhumy gt 1.0) then rhumy=rhumy/100.0
      rhs=temp+alatnt(temp)*rmix(temp,press,rhumy)/cp
      tbot=temp-40.*(1.-rhumy)-1.
      ttop=temp
      its=0
twlab:
      tw=0.5*(ttop+tbot)
      its=its+1
      alhs=tw+alatnt(tw)*rmix(tw,press,1.0)/cp
      if(abs(alhs-rhs) lt 0.0005*abs(rhs)) then begin
       wet=tw
       return, wet
      endif
      if(alhs lt rhs) then begin
       tbot=tw
      endif else begin
       ttop=tw
      endelse
      if(its lt 100) then goto,twlab
      wet=1000.+tw
return,wet
end
function alatnt,t
;
;  approximates a linear relationship between latent heat of
;  vaporization of water and temperature.
;
     alatnt=2500590.-2350.*t
return,alatnt
end
;****************** Start of Main Program ***********************
pro tephi_auto

common corner, t0, h0, t1, h1, t2, h2
data_type=3
;  data_type - format of cfile
;    1 for Fronts87
;    2 for standard
;    3 for standard but reading in Q instead of Dew point Temp

; specify corners of tephigram plot with t and h values
; n.b 0 - lower left 1 - lower right 2 - upper left

;; THIS IS THE ORIGINAL TEPHI FRAME
;;t0=-23.0 & h0=-27.0 & t1=45.0 & h1=40.0 & t2=-86.0 & h2=30.0

;print, 'Standard Frame (y/n)?'
;frame=get_kbrd(1)
frame='y'

if (frame eq 'y' or frame eq 'Y') then begin
;; THIS IS NEW STANDARD FRAME FOR THE TROPICS *AMT*
t0=-23.0 & h0=-27.0 & t1=45.0 & h1=40.0 & t2=-110.0 & h2=60.0

endif else begin
;; THIS IS A BLOW UP FRAME FOR THE TROPICS *AMT*
t0=-5.0 & h0=-10.0 & t1=35.0 & h1=30.0 & t2=-95.0 & h2=80.0

endelse

;; THIS IS STANDARD FRAME NUMBER 2
;;t0=-23.0 & h0=-27.0 & t1=45.0 & h1=40.0 & t2=-130.0 & h2=40.0


port_size=[16.0,22.0]
window=[0.1,0.1,0.9,0.9]
offsets=[3.0,4.0]

;; OLD DEFAULT FILE
cfile='DATA'

;****************** END of User Options *************************
set_plot,'PS'
gfile='TEPHI_AUTO.ps'
!p.font=0
device, portrait=1,filename=gfile,xsize=port_size(0),ysize=port_size(1), $
       xoffset=offsets(0),yoffset=offsets(1),/color,/HELVETICA,/BOLD   ;*AMT* add color
;RED = [0, 1, 1, 0, 0, 1]     ;Specify the red component of each color.
;GREEN = [0, 1, 0, 1, 0, 1]  ;Specify the green component of each color.
;BLUE = [0, 1, 0, 0, 1, 0]   ;Specify the blue component of each color.

;TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE;Load the first six elements of the color table.

loadct, 0                                             ;*AMT*

;device, landscape=1,filename=gfile,xsize=port_size(0),ysize=port_size(1), $
; establish coordinate system for further plotting
; n.b lenx-1 and leny-1 must be same as xw and yw in tplot and tephi
lenx=10.0
leny=10.0
plot,[0,lenx-1,lenx-1,0,0],[0,0,leny-1,leny-1,0], $
  xstyle=5,ystyle=5,/norm,position=window

maxpt=1000
;
dat=fltarr(3,maxpt)
header=strarr(3,35)


;print, 'Moisture?'
;plotmoi=get_kbrd(1)
plotmoi='y'

;
;  read standard or fronts87 dataset
;
;--------------
ans='n'
repeat begin
;--------------

;BADFILE: read,'enter file name (without .DAT extension): ',cfile

openr,inunit,cfile+'.DAT',/get_lun,ERROR=err
if (err ne 0) then begin
 print, 'ERROR OPENING FILE'
; goto, badfile
 stop
endif

if(data_type eq 1) then begin
    dummy=f87(inunit,header,np,dat)
endif else begin
    np=-1
    dummy=readdat(inunit,header,np,dat)

;    print, 'Standard Deviation?'
;    ans=get_kbrd(1)
ans='n'
    if (ans eq 'y' or ans eq 'Y') then sd_plot=1 else sd_plot=0
;    print, sd_plot

    if (sd_plot eq 1) then begin
      ; reverse the first half
      d3=fltarr(3)
      for count=0,(np-1)/2 do begin
        d3=dat(*,count)
        dat(*,count)=dat(*,np-1-count)
        dat(*,np-1-count)=d3
      endfor

      ; read in second file
BADFILE2: read,'enter 2nd file name (without .DAT extension): ',cfile

      openr,inunit,cfile+'.DAT',/get_lun,ERROR=err
      if (err ne 0) then begin
         print, 'ERROR OPENING FILE'
         goto, badfile2
      endif
      np=np-1
      dummy=readdat(inunit,header,np,dat)      
    endif

endelse
;
;  find dew point
;
case data_type of
1: begin
;  FRONTS87 - calc dew point from T and RH
    for i=0,np-1 do begin ;*amt* was np-1
       if (dat(2,i) gt 100.0) then dat(2,i) = 100.0
       dat(2,i) = dew(dat(1,i),dat(2,i))
    endfor
   end
2:         
;  standard data format dew point is in file
3: begin
;  standard data containing Q instead of dew point
;  calculate Td from Q,T,P
;
;  formula below is rh=(0.1*P*Q)/(epsilon*SVP(T))
;  factor of 0.1 is 100/1000
;  the 100 scales rh and the 1000 converts Q in g/Kg to Kg/Kg
;  n.b need temp in celcius
     for i=0,np-1 do begin ;*AMT* was np-1
        t=dat(1,i)-273.15

;       dat(1,i)=dat(1,i)-273.15
;       rh=(0.1*dat(0,i)*dat(2,i))/(0.622*svp(t))
;       dat(2,i) = dew(t,rh)

; *AMT* Change for TETON's formula from LEM

;       qsat=(t gt 0)*3.8/(dat(0,i)*EXP(-17.2693882* $
;             (dat(1,i)-273.15)/(dat(1,i)-35.86))-6.109) + $
;            (t le 0)*3.8/(dat(0,i)*EXP(-21.8745584* $
;             (dat(1,i)-273.15)/(dat(1,i)-7.66))-6.109) 

       qsat=qsaturation(dat(1,i),dat(0,i))

       dat(1,i) = dat(1,i)-273.15      
       rh=0.001*dat(2,i)/qsat
       dat(2,i) = dew(t,rh)
     endfor
   end
endcase
;
;  plot T and Td lines and border and titles
;
teph
;
;   plot data
;
tplot,dat,np,header,plotmoi,sd_plot


;----------------------
;  print, 'Again?'
;  ans=get_kbrd(1)
ans='n'
endrep until ans ne 'y'
;----------------------

device,/close
set_plot,'X'
free_lun,inunit
; reset plotting parameters
!p.linestyle=0
!p.thick=1.0
!p.noclip=1
!p.charsize=1.0
!p.font=-1
      
end

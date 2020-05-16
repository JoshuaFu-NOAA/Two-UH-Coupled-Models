!c===================================================================c 
!c  modification for pgi compiler
!c  Qing Bao . Oct. 20, 2005 
!c===================================================================c

subroutine restart_ocn
  
   USE mo_couple, ONLY: restart_u1, restart_u2, restart_v1, restart_v2, &
       restart_h1, restart_h2, restart_tmp, restart_itt
       parameter(im=720, jm=121)
        common /var0/ u2(im,jm,2),v2(im,jm,2),h2(im,jm,2),&
                     u1(im,jm,2),v1(im,jm,2),h1(im,jm,2),tmp(im,jm,2)
        common /scalar/ grav,omega,itt,ndays,tfr,tfnz,tfnm,epsln,dt,&
                       degrad,nergy,ansamp,alpha,nrstrt,maps,nsmth,nf,&
                       netr,tmn,tmpmax,thmexp,rhoa,cd,cw,wm,tn,r,&
                       gamma1,gamma2,nbuoy,nbounx,nbouny,forcef,forcev,&
                       deltat,nseason,nsmthstr,hepsln,nbeta,beta,&
                       ittend,tspd,txyk1,txyk,cdfv,ufcoef,&
                       q0unit,earthr,daysave,ittwrt,ittatm,phlag,&
                       lcktim,imthlck,day,nday,iyear,tcount,tadd0,&
                       nladv,cdfv1,wdh,theta,timep,nqin,ncouple,&
                       iloop
      
   !c******************************************************************** 
!  include 'pram1.inc'
!  include 'com1.inc'
   !c******************************************************************** 

     restart_itt = itt
     restart_u1  = u1
     restart_u2  = u2
     restart_v1  = v1
     restart_v2  = v2
     restart_h1  = h1
     restart_h2  = h2
     restart_tmp = tmp
end subroutine restart_ocn

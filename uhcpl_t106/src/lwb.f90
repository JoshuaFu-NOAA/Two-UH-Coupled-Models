!+ computes black-body functions for longwave calculations
!+ $Id: lwb.f90,v 1.5 1998/12/11 13:23:34 m214003 Exp $

SUBROUTINE lwb(kdlon,kflev,ptave,ptl,pb,pbint,pbsuin,pbtop,pdbsl,pga,pgb, &
&      pgasur,pgbsur,pgatop,pgbtop,pgc,pgd,pgcsur,pgdsur,pgctop,pgdtop)

  ! Description:
  !
  ! Computes planck functions.
  !
  ! Method:
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! ptave  : (kdlon,kflev)        ; temperature
  ! ptl    : (kdlon,0:kflev)      ; half level temperature
  ! ==== outputs ===
  ! pbint  : (kdlon,0:kflev)      ; half level Planck function
  ! pbsuin : (kdlon)              ; surface Planck function
  ! pbtop  : (kdlon,nsint)        ; top spectral Planck function
  ! pdbsl  : (kdlon,nsint,kflev*2); sub-layer Planck function gradient
  ! pga    : (kdlon,8,2,kflev)    ; db/dt-weighted layer pade approximants
  ! pgb    : (kdlon,8,2,kflev)    ; db/dt-weighted layer pade approximants
  ! pgc    :(kdlon,5,2,kflev)     ; db/dt-weighted layer pade approximants
  ! pgd    :(kdlon,5,2,kflev)     ; db/dt-weighted layer pade approximants
  !
  ! 1. Computes the Planck function on all levels and half levels
  ! from a polynomial development of Planck function.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model".
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_longwave

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdlon, kflev

  !  Array arguments 
  REAL :: pb(kdlon,nsint,kflev+1), pbint(kdlon,kflev+1), pbsuin(kdlon), &
&      pbtop(kdlon,nsint), pdbsl(kdlon,nsint,kflev*2), pga(kdlon,8,2,kflev), &
&      pgasur(kdlon,8,2), pgatop(kdlon,8,2), pgb(kdlon,8,2,kflev), &
&      pgbsur(kdlon,8,2), pgbtop(kdlon,8,2), pgc(kdlon,5,2,kflev), &
&      pgcsur(kdlon,5,2), pgctop(kdlon,5,2), pgd(kdlon,5,2,kflev), &
&      pgdsur(kdlon,5,2), pgdtop(kdlon,5,2), ptave(kdlon,kflev), &
&      ptl(kdlon,kflev+1)

  !  Local scalars: 
  REAL :: zdst1, zdsto1, zdstox, zdstx, zres, zres2, zti, zti2
  INTEGER :: indsu, indt, indto, indtp, ixtox, ixtx, jf, jg, jk, jk1, jk2, &
&      jl, jnu

  !  Local arrays: 
  REAL :: zblay(kdlon,kflev), zblev(kdlon,kflev+1)
  INTEGER :: indb(kdlon), inds(kdlon)

  !  Intrinsic functions 
  INTRINSIC INT, MAX, MERGE, MIN


  !  Executable statements 

!-- 1.0 Planck functions and gradients

  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pbint(jl,jk) = 0.
    END DO
  END DO
  DO jl = 1, kdlon
    pbsuin(jl) = 0.
  END DO

  DO jnu = 1, nsint

!-- 1.1 Levels from surface to kflev

    DO jk = 1, kflev
      DO jl = 1, kdlon
        zti = (ptl(jl,jk)-tstand)/tstand
        zres = xp(1,jnu) + zti*(xp(2,jnu)+zti*(xp(3, &
&            jnu)+zti*(xp(4,jnu)+zti*(xp(5,jnu)+zti*(xp(6,jnu))))))
        pbint(jl,jk) = pbint(jl,jk) + zres
        pb(jl,jnu,jk) = zres
        zblev(jl,jk) = zres
        zti2 = (ptave(jl,jk)-tstand)/tstand
        zres2 = xp(1,jnu) + zti2*(xp(2,jnu)+zti2*(xp(3, &
&            jnu)+zti2*(xp(4,jnu)+zti2*(xp(5,jnu)+zti2*(xp(6,jnu))))))
        zblay(jl,jk) = zres2
      END DO
    END DO

!-- 1.2 Top of the atmosphere and surface

    DO jl = 1, kdlon
      zti = (ptl(jl,kflev+1)-tstand)/tstand
      zti2 = (ptl(jl,1)-tstand)/tstand
      zres = xp(1,jnu) + zti*(xp(2,jnu)+zti*(xp(3,jnu)+zti*(xp(4,jnu)+zti*(xp &
&          (5,jnu)+zti*(xp(6,jnu))))))
      zres2 = xp(1,jnu) + zti2*(xp(2,jnu)+zti2*(xp(3, &
&          jnu)+zti2*(xp(4,jnu)+zti2*(xp(5,jnu)+zti2*(xp(6,jnu))))))
      pbint(jl,kflev+1) = pbint(jl,kflev+1) + zres
      pb(jl,jnu,kflev+1) = zres
      zblev(jl,kflev+1) = zres
      pbtop(jl,jnu) = zres
      pbsuin(jl) = pbsuin(jl) + zres2
    END DO

!-- 1.3 Gradients in sub-layers

    DO jk = 1, kflev
      jk2 = 2*jk
      jk1 = jk2 - 1
      DO jl = 1, kdlon
        pdbsl(jl,jnu,jk1) = zblay(jl,jk) - zblev(jl,jk)
        pdbsl(jl,jnu,jk2) = zblev(jl,jk+1) - zblay(jl,jk)
      END DO

    END DO
  END DO

!-- 2.0 Choose the relevant sets of pade approximants

  DO jl = 1, kdlon
    zdsto1 = (ptl(jl,kflev+1)-tintp(1))/tstp
    ixtox = MAX(1,MIN(mxixt,INT(zdsto1+1.)))
    zdstox = (ptl(jl,kflev+1)-tintp(ixtox))/tstp
    indto = MERGE(ixtox,ixtox+1,zdstox<0.5)
    indb(jl) = indto
    zdst1 = (ptl(jl,1)-tintp(1))/tstp
    ixtx = MAX(1,MIN(mxixt,INT(zdst1+1.)))
    zdstx = (ptl(jl,1)-tintp(ixtx))/tstp
    indt = MERGE(ixtx,ixtx+1,zdstx<0.5)
    inds(jl) = indt
  END DO

  DO jf = 1, 2
    DO jg = 1, 8
      DO jl = 1, kdlon
        indsu = inds(jl)
        pgasur(jl,jg,jf) = ga(indsu,2*jg-1,jf)
        pgbsur(jl,jg,jf) = gb(indsu,2*jg-1,jf)
        indtp = indb(jl)
        pgatop(jl,jg,jf) = ga(indtp,2*jg-1,jf)
        pgbtop(jl,jg,jf) = gb(indtp,2*jg-1,jf)
      END DO
    END DO
  END DO

  DO jf = 1, 2
    DO jg = 1, 5
      DO jl = 1, kdlon
        indsu = inds(jl)
        pgcsur(jl,jg,jf) = gc(indsu,2*jg-1,jf)
        pgdsur(jl,jg,jf) = gd(indsu,2*jg-1,jf)
        indtp = indb(jl)
        pgctop(jl,jg,jf) = gc(indtp,2*jg-1,jf)
        pgdtop(jl,jg,jf) = gd(indtp,2*jg-1,jf)
      END DO
    END DO

  END DO

  DO jk = 1, kflev
    DO jl = 1, kdlon
      zdst1 = (ptave(jl,jk)-tintp(1))/tstp
      ixtx = MAX(1,MIN(mxixt,INT(zdst1+1.)))
      zdstx = (ptave(jl,jk)-tintp(ixtx))/tstp
      indt = MERGE(ixtx,ixtx+1,zdstx<0.5)
      indb(jl) = indt
    END DO

    DO jf = 1, 2
      DO jg = 1, 8
        DO jl = 1, kdlon
          indt = indb(jl)
          pga(jl,jg,jf,jk) = ga(indt,2*jg,jf)
          pgb(jl,jg,jf,jk) = gb(indt,2*jg,jf)
        END DO
      END DO
    END DO
    DO jf = 1, 2
      DO jg = 1, 5
        DO jl = 1, kdlon
          indt = indb(jl)
          pgc(jl,jg,jf,jk) = gc(indt,2*jg,jf)
          pgd(jl,jg,jf,jk) = gd(indt,2*jg,jf)
        END DO
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE lwb

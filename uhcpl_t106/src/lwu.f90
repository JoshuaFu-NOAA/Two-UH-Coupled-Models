!+ longwave effective absorber amounts
!+ $Id: lwu.f90,v 1.10 1999/07/20 14:21:37 m214003 Exp $

SUBROUTINE lwu(kdlon,kflev,kewaer,kaerh,paer,pcco2,pdp,ppmb,pqof,ptave,pwv, &
&      pabcu,pcfcabs)

  ! Description:
  !
  ! computes absorber amounts including pressure and
  ! temperature effects
  !
  ! Method:
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! paer   : (kdlon,kflev,5+kewaer); aerosol optical thickness (1,,5)
  ! tanre et al., 1984
  ! aerosol mass mixing ratio (kg/kg)
  ! (6,.,5+kewaer) computed in echam4
  ! pcco2  :                       ; concentration in co2 (pa/pa)
  ! pdp    : (kdlon,kflev)         ; layer pressure thickness (pa)
  ! ppmb   : (kdlon,0:kflev)       ; half level pressure
  ! pqof   : (kdlon,kflev)         ; concentration in ozone (pa/pa)
  ! ptave  : (kdlon,kflev)         ; temperature
  ! pwv    : (kdlsur,kflev)        ; specific humidity pa/pa
  ! ==== outputs ===
  ! kx.    : (kdlon)               ; temperature indices
  ! pabcu  :(kdlon,nua,3*kflev+1)  ; effective absorber amounts
  !
  ! 1. Computes the pressure and temperature weighted amounts of
  !    absorbers.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! R. Van Dorland, KNMI, May 1992, changed
  ! R. Van Dorland, KNMI, May 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants
  USE mo_aerosols
  USE mo_longwave
  USE mo_cfc
  USE mo_radint
  USE mo_radiation

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: pcco2
  INTEGER :: kdlon, kewaer, kflev

  !  Array arguments 
  REAL :: pabcu(kdlon,nua,3*kflev+1), paer(kdlon,kflev,5+kewaer), pcfcabs(4), &
&      pdp(kdlon,kflev), ppmb(kdlon,kflev+1), pqof(kdlon,kflev), &
&      ptave(kdlon,kflev), pwv(kdlon,kflev)
  INTEGER :: kaerh(kdlon,kflev)

  !  Local scalars: 
  REAL :: zalup, zcac8, zcah1, zcah2, zcah3, zcah4, zcah5, zcah6, zcam4, &
&      zcam5, zcan1, zcan2, zcan3, zcbc8, zcbh1, zcbh2, zcbh3, zcbh4, zcbh5, &
&      zcbh6, zcbm4, zcbm5, zcbn1, zcbn2, zcbn3, zch4up, zdiff, zdpm, zfppw, &
&      zn2oup, ztavi, ztx, ztx2, zu6, zup, zupm, zupmco2, zupmh2o, zupmo3, &
&      zzably, zzabme, zzabni
  INTEGER :: iae, icae, ig1, ih, jae1, jae2, jae3, jaer, jcp1, jj, jjpn, &
&      jk, jki, jkip1, jkj, jkjp, jkjpn, jkjr, jkk, jkl, jl

  !  Local arrays: 
  REAL :: zably(kdlon,nua,3*kflev+1), zduc(kdlon,3*kflev+1), zphio(kdlon), &
&      zphio2(kdlon), zpsc2(kdlon), zpsc3(kdlon), zpsh1(kdlon), zpsh2(kdlon), &
&      zpsh3(kdlon), zpsh4(kdlon), zpsh5(kdlon), zpsh6(kdlon), zpsio(kdlon), &
&      zpsio2(kdlon), zpsm3(kdlon), zpsm6(kdlon), zpsn2(kdlon), zpsn3(kdlon), &
&      zpsn6(kdlon), zssig(kdlon,3*kflev+1), ztcon(kdlon), zuaer(kdlon,nsint), &
&      zxoz(kdlon), zxwv(kdlon)

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF,ALOGHF
!DIR$ VFUNCTION EXPHF
!DIR$ VFUNCTION ALOGHF
#define EXP(x)  EXPHF(x)
#define LOG(x)  ALOGHF(x)
#else
  INTRINSIC EXP, LOG
#endif
  INTRINSIC MAX, MIN

!DIR$ NOBOUNDS PABCU

  !  Executable statements 

!-- 1. Initialization

  zdiff = diff

!-- 2. Pressure over gauss sub-levels

  DO jl = 1, kdlon
    zssig(jl,1) = ppmb(jl,1)*100.
  END DO

  DO jk = 1, kflev
    jkj = (jk-1)*ng1p1 + 1
    jkjr = jkj
    jkjp = jkj + ng1p1
    DO jl = 1, kdlon
      zssig(jl,jkjp) = ppmb(jl,jk+1)*100.
    END DO
    DO ig1 = 1, ng1
      jkj = jkj + 1
      DO jl = 1, kdlon
        zssig(jl,jkj) = (zssig(jl,jkjr)+zssig(jl,jkjp))*0.5 + &
&            rt1(ig1)*(zssig(jl,jkjp)-zssig(jl,jkjr))*0.5
      END DO
    END DO
  END DO

!-- 3. Pressure thickness and mean pressure of sub-layers

  DO jki = 1, 3*kflev
    jkip1 = jki + 1
    DO jl = 1, kdlon
      zably(jl,5,jki) = (zssig(jl,jki)+zssig(jl,jkip1))*0.5
      zably(jl,3,jki) = (zssig(jl,jki)-zssig(jl,jkip1))/(10.*g)
    END DO
  END DO

  DO jk = 1, kflev
    jkl = kflev + 1 - jk
    DO jl = 1, kdlon
      zxwv(jl) = MAX(pwv(jl,jkl),zepscq)
      zxoz(jl) = MAX(pqof(jl,jkl)/pdp(jl,jkl),zepsco)
    END DO
    jkj = (jk-1)*ng1p1 + 1
    jkjpn = jkj + ng1
    DO jkk = jkj, jkjpn
      DO jl = 1, kdlon
        zdpm = zably(jl,3,jkk)
        zupm = zably(jl,5,jkk)*zdpm/101325.
        zupmh2o = (zably(jl,5,jkk)+pvgh2o)*zdpm/101325.
        zupmco2 = (zably(jl,5,jkk)+pvgco2)*zdpm/101325.
        zupmo3 = (zably(jl,5,jkk)+pvgo3)*zdpm/101325.
        zduc(jl,jkk) = zdpm
        zably(jl,12,jkk) = zxoz(jl)*zdpm
        zably(jl,13,jkk) = zxoz(jl)*zupmo3
        zably(jl,28,jkk) = zxoz(jl)*zdpm
        zably(jl,29,jkk) = zxoz(jl)*zupmo3
        zu6 = zxwv(jl)*zupm
        zfppw = 1.6078*zxwv(jl)/(1.+0.608*zxwv(jl))
        zably(jl,6,jkk) = zxwv(jl)*zupmh2o
        zably(jl,11,jkk) = zu6*zfppw
        zably(jl,10,jkk) = zu6*(1.-zfppw)
        zably(jl,9,jkk) = pcco2*zupmco2
        zably(jl,8,jkk) = pcco2*zdpm
        zably(jl,16,jkk) = znitox*zupm
        zably(jl,18,jkk) = zmetha*zupm
        zably(jl,22,jkk) = zdpm
      END DO
    END DO
  END DO

!-- 4. Cumulative absorber amounts from top of atmosphere

  pabcu(:,:,3*kflev+1) = 0.

  DO jk = 1, kflev
    jj = (jk-1)*ng1p1 + 1
    jjpn = jj + ng1
    jkl = kflev + 1 - jk

!-- 4.1 Cumulative aerosol amounts from top of atmosphere

    jae1 = 3*kflev + 1 - jj
    jae2 = 3*kflev + 1 - (jj+1)
    jae3 = 3*kflev + 1 - jjpn
    DO iae = 1, 5
      DO jl = 1, kdlon
        zuaer(jl,iae) = (caer(iae,1)*paer(jl,jkl,1)+caer(iae,2)*paer(jl,jkl,2 &
&            )+caer(iae,3)*paer(jl,jkl,3)+caer(iae,4)*paer(jl,jkl,4)+ &
&            caer(iae,5)*paer(jl,jkl,5))/(zduc(jl,jae1)+zduc(jl,jae2)+zduc(jl, &
&            jae3))
      END DO
      ! Contribution gads aerosols
      DO jaer = 1, kewaer
        DO jl = 1, kdlon
          icae = ndfaer(jaer)
          ih = kaerh(jl,jk)
          zuaer(jl,iae) = zuaer(jl,iae) + caern(ih,iae,icae)*fcvaer(icae)* &
&              paer(jl,jkl,5+jaer)
        END DO
      END DO
    END DO

!-- 4.2 Introduces temperature effects on absorber amounts

    DO jl = 1, kdlon
      ztavi = ptave(jl,jkl)
      ztcon(jl) = EXP(6.08*(296./ztavi-1.))
      ztx = ztavi - tref
      ztx2 = ztx*ztx
      zzably = zably(jl,6,jae1) + zably(jl,6,jae2) + zably(jl,6,jae3)
      zup = MIN(MAX(0.5*c10e*LOG(zzably)+5.,0.),6.0)
      zcah1 = at(1,1) + zup*(at(1,2)+zup*(at(1,3)))
      zcbh1 = bt(1,1) + zup*(bt(1,2)+zup*(bt(1,3)))
      zpsh1(jl) = EXP(zcah1*ztx+zcbh1*ztx2)
      zcah2 = at(2,1) + zup*(at(2,2)+zup*(at(2,3)))
      zcbh2 = bt(2,1) + zup*(bt(2,2)+zup*(bt(2,3)))
      zpsh2(jl) = EXP(zcah2*ztx+zcbh2*ztx2)
      zcah3 = at(3,1) + zup*(at(3,2)+zup*(at(3,3)))
      zcbh3 = bt(3,1) + zup*(bt(3,2)+zup*(bt(3,3)))
      zpsh3(jl) = EXP(zcah3*ztx+zcbh3*ztx2)
      zcah4 = at(4,1) + zup*(at(4,2)+zup*(at(4,3)))
      zcbh4 = bt(4,1) + zup*(bt(4,2)+zup*(bt(4,3)))
      zpsh4(jl) = EXP(zcah4*ztx+zcbh4*ztx2)
      zcah5 = at(5,1) + zup*(at(5,2)+zup*(at(5,3)))
      zcbh5 = bt(5,1) + zup*(bt(5,2)+zup*(bt(5,3)))
      zpsh5(jl) = EXP(zcah5*ztx+zcbh5*ztx2)
      zcah6 = at(6,1) + zup*(at(6,2)+zup*(at(6,3)))
      zcbh6 = bt(6,1) + zup*(bt(6,2)+zup*(bt(6,3)))
      zpsh6(jl) = EXP(zcah6*ztx+zcbh6*ztx2)

      zzably = zably(jl,9,jae1) + zably(jl,9,jae2) + zably(jl,9,jae3)
      zalup = c10e*LOG(zzably)
      zup = MAX(0.0,5.0+0.5*zalup)
      zpsc2(jl) = (ztavi/tref)**zup
      zcac8 = at(8,1) + zup*(at(8,2)+zup*(at(8,3)))
      zcbc8 = bt(8,1) + zup*(bt(8,2)+zup*(bt(8,3)))
      zpsc3(jl) = EXP(zcac8*ztx+zcbc8*ztx2)
      zphio(jl) = EXP(oct(1)*ztx+oct(2)*ztx2)
      zpsio(jl) = EXP(2.*(oct(3)*ztx+oct(4)*ztx2))
      zphio2(jl) = EXP(odt(1)*ztx+odt(2)*ztx2)
      zpsio2(jl) = EXP(2.*(odt(3)*ztx+odt(4)*ztx2))
      zzabni = zably(jl,16,jae1) + zably(jl,16,jae2) + zably(jl,16,jae3)
      zn2oup = MAX(-1.5,5.0+0.5*c10e*LOG(zzabni))
      zcan1 = ct(1,1) + zn2oup*(ct(1,2)+zn2oup*ct(1,3))
      zcbn1 = dt(1,1) + zn2oup*(dt(1,2)+zn2oup*dt(1,3))
      zpsn2(jl) = EXP(zcan1*ztx+zcbn1*ztx2)
      zcan2 = ct(2,1) + zn2oup*(ct(2,2)+zn2oup*ct(2,3))
      zcbn2 = dt(2,1) + zn2oup*(dt(2,2)+zn2oup*dt(2,3))
      zpsn3(jl) = EXP(zcan2*ztx+zcbn2*ztx2)
      zcan3 = ct(3,1) + zn2oup*(ct(3,2)+zn2oup*ct(3,3))
      zcbn3 = dt(3,1) + zn2oup*(dt(3,2)+zn2oup*dt(3,3))
      zpsn6(jl) = EXP(zcan3*ztx+zcbn3*ztx2)
      zzabme = zably(jl,18,jae1) + zably(jl,18,jae2) + zably(jl,18,jae3)
      zch4up = MAX(-1.5,5.0+0.5*c10e*LOG(zzabme))
      zcam4 = ct(4,1) + zch4up*(ct(4,2)+zch4up*ct(4,3))
      zcbm4 = dt(4,1) + zch4up*(dt(4,2)+zch4up*dt(4,3))
      zpsm3(jl) = EXP(zcam4*ztx+zcbm4*ztx2)
      zcam5 = ct(5,1) + zch4up*(ct(5,2)+zch4up*ct(5,3))
      zcbm5 = dt(5,1) + zch4up*(dt(5,2)+zch4up*dt(5,3))
      zpsm6(jl) = EXP(zcam5*ztx+zcbm5*ztx2)
    END DO

    DO jkk = jj, jjpn
      jc = 3*kflev + 1 - jkk
      jcp1 = jc + 1
      DO jl = 1, kdlon
        pabcu(jl,10,jc) = pabcu(jl,10,jcp1) + zably(jl,10,jc)*zdiff
        pabcu(jl,11,jc) = pabcu(jl,11,jcp1) + zably(jl,11,jc)*ztcon(jl)*zdiff

        pabcu(jl,12,jc) = pabcu(jl,12,jcp1) + zably(jl,12,jc)*zphio(jl)*zdiff
        pabcu(jl,13,jc) = pabcu(jl,13,jcp1) + zably(jl,13,jc)*zpsio(jl)*zdiff

        pabcu(jl,28,jc) = pabcu(jl,28,jcp1) + zably(jl,28,jc)*zphio2(jl)* &
&            zdiff
        pabcu(jl,29,jc) = pabcu(jl,29,jcp1) + zably(jl,29,jc)*zpsio2(jl)* &
&            zdiff

        pabcu(jl,7,jc) = pabcu(jl,7,jcp1) + zably(jl,9,jc)*zpsc2(jl)*zdiff
        pabcu(jl,8,jc) = pabcu(jl,8,jcp1) + zably(jl,9,jc)*zpsc3(jl)*zdiff
        pabcu(jl,9,jc) = pabcu(jl,9,jcp1) + zably(jl,9,jc)*zpsc3(jl)*zdiff

        pabcu(jl,1,jc) = pabcu(jl,1,jcp1) + zably(jl,6,jc)*zpsh1(jl)*zdiff
        pabcu(jl,2,jc) = pabcu(jl,2,jcp1) + zably(jl,6,jc)*zpsh2(jl)*zdiff
        pabcu(jl,3,jc) = pabcu(jl,3,jcp1) + zably(jl,6,jc)*zpsh5(jl)*zdiff
        pabcu(jl,4,jc) = pabcu(jl,4,jcp1) + zably(jl,6,jc)*zpsh3(jl)*zdiff
        pabcu(jl,5,jc) = pabcu(jl,5,jcp1) + zably(jl,6,jc)*zpsh4(jl)*zdiff
        pabcu(jl,6,jc) = pabcu(jl,6,jcp1) + zably(jl,6,jc)*zpsh6(jl)*zdiff

        pabcu(jl,14,jc) = pabcu(jl,14,jcp1) + zably(jl,16,jc)*zpsn2(jl)*zdiff
        pabcu(jl,15,jc) = pabcu(jl,15,jcp1) + zably(jl,16,jc)*zpsn3(jl)*zdiff
        pabcu(jl,16,jc) = pabcu(jl,16,jcp1) + zably(jl,16,jc)*zpsn6(jl)*zdiff

        pabcu(jl,17,jc) = pabcu(jl,17,jcp1) + zably(jl,18,jc)*zpsm3(jl)*zdiff
        pabcu(jl,18,jc) = pabcu(jl,18,jcp1) + zably(jl,18,jc)*zpsm6(jl)*zdiff

        pabcu(jl,19,jc) = pabcu(jl,19,jcp1) + zably(jl,22,jc)*pcfcabs(1)* &
&            zdicfc
        pabcu(jl,20,jc) = pabcu(jl,20,jcp1) + zably(jl,22,jc)*pcfcabs(2)* &
&            zdicfc
        pabcu(jl,21,jc) = pabcu(jl,21,jcp1) + zably(jl,22,jc)*pcfcabs(3)* &
&            zdicfc
        pabcu(jl,22,jc) = pabcu(jl,22,jcp1) + zably(jl,22,jc)*pcfcabs(4)* &
&            zdicfc

        pabcu(jl,23,jc) = pabcu(jl,23,jcp1) + zuaer(jl,1)*zduc(jl,jc)*zdiff
        pabcu(jl,24,jc) = pabcu(jl,24,jcp1) + zuaer(jl,2)*zduc(jl,jc)*zdiff
        pabcu(jl,25,jc) = pabcu(jl,25,jcp1) + zuaer(jl,3)*zduc(jl,jc)*zdiff
        pabcu(jl,26,jc) = pabcu(jl,26,jcp1) + zuaer(jl,4)*zduc(jl,jc)*zdiff
        pabcu(jl,27,jc) = pabcu(jl,27,jcp1) + zuaer(jl,5)*zduc(jl,jc)*zdiff
      END DO

    END DO
  END DO

  RETURN

END SUBROUTINE lwu

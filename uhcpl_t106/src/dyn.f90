!+ computes adiabatic tendencies and auxilliary hybrid variables.
!+ $Id: dyn.f90,v 1.21 2000/06/05 12:35:22 m214003 Exp $

!OCL NOALIAS

SUBROUTINE dyn(petadot,ulz,vmaxz)

  ! Description:
  !
  ! Computes adiabatic tendencies and auxilliary hybrid variables.
  !
  ! Method:
  !
  ! The primary purpose is to compute adiabatic tendencies.
  ! Auxilliary variables connected with the vertical difference
  ! scheme are saved as they are subsequently required to compute
  ! input to physical parameterizations.
  !
  ! *dyn* is called from *gpc*. 
  ! Input is from long-term storage and modules *mo_hyb* and *mo_gaussgrid*
  ! Output is to long-term storage.
  !
  ! Externals:
  ! *pres* and *auxhyb* are called to calculate auxilliary variables.
  ! *geopot* is called to calculate full-level geopotentials.
  !
  ! *External documentation of the model equations, of
  ! the organization of the vertical calculation, and of the
  ! organization of the spectral horizontal calculation.
  !
  ! Authors:
  !
  ! A. J. Simmons, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! I. Kirchner, MPI, May 2000, tendency diagnostics revision
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_scan_buffer,   ONLY: alnpr_scb, alpha_scb, alps_scb, alpse_scb,   &
                              d_scb, qe_scb, rh_scb, t_scb, te_scb, u_scb, &
                              ul_scb, v_scb, vo_scb, vervel_scb, vol_scb,  &
                              vom_scb, xe_scb, dalpsl_scb, dalpsm_scb,     &
                              dtl_scb, dtm_scb
  USE mo_memory_gl,     ONLY: q, x
  USE mo_memory_g3a,    ONLY: apsm, geospm
  USE mo_memory_g3b,    ONLY: aps
  USE mo_control,       ONLY: ltdiag, nlev, nlevp1, nvclev, vct, ngl
  USE mo_gaussgrid,     ONLY: coriol, cst, rcsth, sqcst
  USE mo_hyb,           ONLY: alpham, alrrdic, ardprc, cetah, cpg, delb,    &
                              delpr, nlevm1, nlmsgl, nlmslp, nplev, nplvp1, &
                              ralpha, rddelb, rdlnp0i, rdt0ral, rlnpr, t0icao
  USE mo_constants,     ONLY: cpd, rcpd, rd, vtmpc1, vtmpc2
  USE mo_diag_tendency, ONLY: pdiga, pdsga
  USE mo_global_op,     ONLY: maxval_zonal, minval_zonal
  USE mo_decomposition, ONLY: dc => local_decomposition

  IMPLICIT NONE

  !  Array arguments (now indexed by continuous latitude index)

  REAL ,INTENT(OUT) :: petadot(dc% nglon ,nlevp1, dc% nglat)
  REAL ,INTENT(OUT) :: ulz    (           nlev  , dc% nglat)
  REAL ,INTENT(OUT) :: vmaxz  (           nlev  , dc% nglat)

  !  Local array bounds
  INTEGER :: nglon, nglpx, nglat

  !  Local scalars: 
  REAL :: dpde, sdoth, zalpha, zalphm, zardrc, zcorio, zcpdrrd, zcpg, zcst, &
&      zdelb, zdlp, zdt, zetdpde, zln, zlnpr, zrcsth, zrddlb, zrdwrp, &
&      zrgr, zrrd, zrtgr, ztq, zvcikp
  INTEGER :: ikp, irow, jrow, jglat, irow1, jk, jl

  !  Local arrays: 
  REAL :: zdpsl (dc%nglpx),         zdpsm  (dc%nglpx),       &
          zgeos (dc%nglpx),         zrdelp (dc%nglpx, nlev), &
          zsdiv (dc%nglpx, nlevp1), ztbar  (dc%nglpx, nlev), &
          ztv   (dc%nglpx, nlev),   zvgrps (dc%nglpx, nlev), &
          delp  (dc%nglpx,nlev),    aph   (dc%nglpx,nlevp1), &
          zrcst (dc%nglat)

  !  External subroutines 
  EXTERNAL auxhyb, geopot, pres

  !  Intrinsic functions 
  INTRINSIC EXP, SQRT


  !  Executable statements 

  !  Local array bounds
  nglon = dc%nglon ! number of longitudes
  nglpx = dc%nglpx ! number of longitudes allocated
  nglat = dc%nglat ! number of latitudes

  ! Latitude indices

  DO jrow = 1, nglat
    jglat = dc%glat(jrow)                 ! global continuous north -> south
    irow  = min(2*jglat-1,2*(ngl+1-jglat)) ! global ping pong index

!-- 1. Preliminary calculations

!-- 1.1 Set local values

    zcorio = coriol(irow)
    irow1 = (irow+1)/2
    zcst = cst(irow1)
    zrcsth = rcsth(irow1)
    zcpdrrd = cpd/rd
    zrrd = 1.0/rd

!-- 1.2 Compute surface pressure and its gradient

    aph(:,nlevp1) = EXP(alps_scb(:,jrow))
    zdpsl(:) = aph(:,nlevp1)*dalpsl_scb(:,jrow)
    zdpsm(:) = aph(:,nlevp1)*dalpsm_scb(:,jrow)

!-- 1.3 Compute half-level pressures and auxilliary variables.

    CALL pres(aph,nglpx,aph(1,nlevp1),nglon)
    CALL auxhyb(delp,zrdelp,alnpr_scb(:,:,jrow),alpha_scb(:,:,jrow),&
              aph,nglpx,nglon)

!-- 1.4 Compute v.grad(ps)

    DO jk = nplvp1, nlev

      DO jl = 1, nglon
        zvgrps(jl,jk) = u_scb(jl,jk,jrow)*zdpsl(jl) + &
                        v_scb(jl,jk,jrow)*zdpsm(jl)
      END DO
    END DO

!-- 2. Sum divergence and compute surface pressure tendency

!-- 2.1 Compute pressure-level sums

    zsdiv (:nglon,1) = 0.

    DO jk = 1, nplev
      ikp = jk + 1
      zdlp = delpr(jk)

      DO jl = 1, nglon
        zsdiv(jl,ikp) = d_scb(jl,jk,jrow)*zdlp + zsdiv(jl,jk)
      END DO
    END DO

!-- 2.2 Compute hybrid-level sums

    DO jk = nplvp1, nlev
      ikp = jk + 1
      zdelb = delb(jk)

      DO jl = 1, nglon
        zsdiv(jl,ikp) = d_scb(jl,jk,jrow)*delp(jl,jk) + &
                        zdelb*zvgrps(jl,jk) + zsdiv(jl,jk)
      END DO
    END DO

    IF (ltdiag) THEN
       ! store pressure tendency term for each layer
       DO jk = 1,nlev
          pdiga(:,jk,20,jrow) = (zsdiv(:,jk+1) - zsdiv(:,jk))/aph(:,nlevp1)
       ENDDO
    END IF

!-- 2.3 Tendency of logarithm of surface pressure

    DO jl = 1, nglon
      alpse_scb(jl,jrow) = -zsdiv(jl,nlevp1)/aph(jl,nlevp1)
    END DO
    IF (ltdiag) &
      pdsga(1:nglon,1,jrow) = pdsga(1:nglon,1,jrow) + alpse_scb(1:nglon,jrow)

!-- 3. Compute reference temperature and deviation
!      of virtual temprature

    DO jl = 1, nglon
      zgeos(jl) = rd*alps_scb(jl,jrow)
    END DO

    DO jk = nlev, 1, -1
      DO jl = 1, nglon
        ztbar(jl,jk) = t0icao*EXP(alrrdic*(zgeos(jl)-&
                       alpha_scb(jl,jk,jrow)-rdlnp0i))
        zgeos(jl) = zgeos(jl) - alnpr_scb(jl,jk,jrow)
        ztv(jl,jk) = t_scb(jl,jk,jrow)*&
          (1.+vtmpc1*q(jl,jk,jrow)-x(jl,jk,jrow)) - ztbar(jl,jk)
      END DO
    END DO

!-- 4. Compute vertical advection
!      and do zonal mean and box diagnostics.

!-- 4.1 Compute vertical advection

    vom_scb    (:,1,jrow) = 0.
    vol_scb    (:,1,jrow) = 0.
    te_scb     (:,1,jrow) = 0.
    qe_scb     (:,1,jrow) = 0.
    xe_scb     (:,1,jrow) = 0.
    vervel_scb (:,1,jrow) = 0.

    DO jk = 1, nlevm1
      ikp = jk + 1
      zvcikp = vct(nvclev+ikp)

      DO jl = 1, nglon
        sdoth = 0.5*(zvcikp*zsdiv(jl,nlevp1)-zsdiv(jl,ikp))
        vom_scb(jl,jk ,jrow) = (u_scb(jl,jk,jrow)-u_scb(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + vom_scb(jl,jk,jrow)
        vom_scb(jl,ikp,jrow) = (u_scb(jl,jk,jrow)-u_scb(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))
        vol_scb(jl,jk ,jrow) = (v_scb(jl,jk,jrow)-v_scb(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + vol_scb(jl,jk,jrow)
        vol_scb(jl,ikp,jrow) = (v_scb(jl,jk,jrow)-v_scb(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))
        te_scb (jl,jk ,jrow) = (t_scb(jl,jk,jrow)-t_scb(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + te_scb(jl,jk,jrow)
        te_scb (jl,ikp,jrow) = (t_scb(jl,jk,jrow)-t_scb(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))

        ! compute eta-dot for semi Lagrangian scheme

        zetdpde = 2.*sdoth
        dpde = (aph(jl,jk+2)-aph(jl,jk))/(cetah(jk+2)-cetah(jk))
        petadot(jl,jk+1,jrow) = zetdpde/dpde

      END DO
    END DO

    IF (ltdiag) THEN
       ! store vertical advection increment
       pdiga(1:nglon,:, 2,jrow) = pdiga(1:nglon,:, 2,jrow) + vom_scb(1:nglon,:,jrow)
       pdiga(1:nglon,:, 7,jrow) = pdiga(1:nglon,:, 7,jrow) + vol_scb(1:nglon,:,jrow)
       pdiga(1:nglon,:,13,jrow) = pdiga(1:nglon,:,13,jrow) + te_scb (1:nglon,:,jrow)
       ! prepare next parts
       pdiga(1:nglon,:, 1,jrow) = pdiga(1:nglon,:, 1,jrow) - vom_scb(1:nglon,:,jrow)
       pdiga(1:nglon,:, 6,jrow) = pdiga(1:nglon,:, 6,jrow) - vol_scb(1:nglon,:,jrow)
       pdiga(1:nglon,:,14,jrow) = pdiga(1:nglon,:,14,jrow) - te_scb (1:nglon,:,jrow)
    END IF

    DO jl = 1, nglon
      petadot(jl,1,jrow) = 0.
      petadot(jl,nlevp1,jrow) = 0.
    END DO

!-- 5. Compute energy conversion term for pressure levels
!      and do zonal mean and box diagnostics.

!-- 5.1 Compute energy conversion term for pressure levels

    DO jk = 1, nplev
      zardrc = ardprc(jk)
      zalphm = alpham(jk)

      DO jl = 1, nglon
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1.+vtmpc2*q(jl,jk,jrow))
        zdt = -ztq*(zsdiv(jl,jk)*zardrc+zalphm*d_scb(jl,jk,jrow))
        te_scb(jl,jk,jrow) = te_scb(jl,jk,jrow) + zdt

      END DO
    END DO

!-- 6. Compute pressure-gradient terms,complete calculation
!      of energy conversion term

!-- 6.1 Hybrid levels

    DO jk = nplvp1, nlmsgl
      zdelb = delb(jk)
      zrddlb = rddelb(jk)
      zcpg = cpg(jk)

!OCL NOVREC,NOALIAS
      DO jl = 1, nglon
        zrgr = (zrddlb+zcpg*alnpr_scb(jl,jk,jrow)*zrdelp(jl,jk))*zrdelp(jl,jk)
        zrtgr = zrgr*ztv(jl,jk)
        vom_scb(jl,jk,jrow) = vom_scb(jl,jk,jrow) - zcst*zrtgr*zdpsl(jl)
        vol_scb(jl,jk,jrow) = vol_scb(jl,jk,jrow) - zcst*zrtgr*zdpsm(jl)
        zrdwrp = (zrgr*zvgrps(jl,jk)-(zrdelp(jl,jk)*(zsdiv(jl,jk)*&
                 alnpr_scb(jl,jk,jrow)+alpha_scb(jl,jk,jrow)*zdelb*&
                 zvgrps(jl,jk))+alpha_scb(jl,jk,jrow)*d_scb(jl,jk,jrow)))*rcpd
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1.+vtmpc2*q(jl,jk,jrow))
        zdt = ztq*zrdwrp
        te_scb(jl,jk,jrow) = te_scb(jl,jk,jrow) + zdt
      END DO
    END DO

!-- 6.2 Sigma levels

    DO jk = nlmslp, nlev
      zalphm = alpham(jk)
      zlnpr  = rlnpr(jk)
      zalpha = ralpha(jk)
      zrddlb = rddelb(jk)
      DO jl = 1, nglon
        zrgr  = zrddlb*zrdelp(jl,jk)
        zrtgr = zrgr*ztv(jl,jk)
        vom_scb(jl,jk,jrow) = vom_scb(jl,jk,jrow) - zcst*zrtgr*zdpsl(jl)
        vol_scb(jl,jk,jrow) = vol_scb(jl,jk,jrow) - zcst*zrtgr*zdpsm(jl)
        zrdwrp = (zrgr*zvgrps(jl,jk)-(zrdelp(jl,jk)*(zsdiv(jl,jk)*zlnpr+zalphm* &
            zvgrps(jl,jk))+zalpha*d_scb(jl,jk,jrow)))*rcpd
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1.+vtmpc2*q(jl,jk,jrow))
        zdt = ztq*zrdwrp
        te_scb(jl,jk,jrow) = te_scb(jl,jk,jrow) + zdt
      END DO
    END DO

    IF (ltdiag) THEN
       ! store energy conversion increment
       pdiga(1:nglon,:,14,jrow) = pdiga(1:nglon,:,14,jrow) + te_scb (1:nglon,:,jrow)
       ! prepare next part
       pdiga(1:nglon,:,12,jrow) = pdiga(1:nglon,:,12,jrow) - te_scb (1:nglon,:,jrow)
    END IF

!-- 6.3 Compute vertical velocity for mass-flux scheme

    DO jk = 1, nplev
      zardrc = ardprc(jk)
      zalphm = alpham(jk)
      DO jl = 1, nglon
        vervel_scb(jl,jk,jrow) = &
          -(zsdiv(jl,jk)*zardrc+zalphm*d_scb(jl,jk,jrow))*zcpdrrd
      END DO
    END DO

    DO jk = nplvp1, nlmsgl
      zdelb  = delb(jk)
      zrddlb = rddelb(jk)
      zcpg   = cpg(jk)
      DO jl = 1, nglon
        zrgr = (zrddlb+zcpg*alnpr_scb(jl,jk,jrow)*zrdelp(jl,jk))*zrdelp(jl,jk)
        vervel_scb(jl,jk,jrow) = (zrgr*zvgrps(jl,jk)-zrdelp(jl,jk)*(zsdiv(jl,jk) &
&            *alnpr_scb(jl,jk,jrow)+alpha_scb(jl,jk,jrow)*zdelb*zvgrps(jl,jk)) &
&            -alpha_scb(jl,jk,jrow)*d_scb(jl,jk,jrow))*zrrd
      END DO
    END DO

    DO jk = nlmslp, nlev
      zalphm = alpham(jk)
      zlnpr  = rlnpr(jk)
      zalpha = ralpha(jk)
      zrddlb = rddelb(jk)
      DO jl = 1, nglon
        zrgr = zrddlb*zrdelp(jl,jk)
        vervel_scb(jl,jk,jrow) = (zrgr*zvgrps(jl,jk)-zrdelp(jl,jk)*(zsdiv(jl,jk) &
&            *zlnpr+zalphm*zvgrps(jl,jk))-zalpha*d_scb(jl,jk,jrow))*zrrd
      END DO
    END DO

    DO jk = 1, nlev
      DO jl = 1, nglon
        vervel_scb(jl,jk,jrow) = vervel_scb(jl,jk,jrow) &
                               * 0.5*(aph(jl,jk)+aph(jl,jk+1))
      END DO
    END DO

!-- 7. Compute geopotential

!-- 7.1 Compute deviation of geopotential height at surface

    DO jl = 1, nglon
      zln = rd*alps_scb(jl,jrow) - rdlnp0i
      zgeos(jl) = geospm(jl,jrow) + rdt0ral*EXP(alrrdic*zln)
    END DO

!-- 7.2 Compute deviation of geopotential height

    CALL geopot(rh_scb(:,:,jrow),ztv,alnpr_scb(:,:,jrow),alpha_scb(:,:,jrow),&
                zgeos,nglpx,nglon)

!-- 8. Compute horizontal advection terms

    DO jk = 1, nlev
      DO jl = 1, nglon
        rh_scb(jl,jk,jrow) = zrcsth*(u_scb(jl,jk,jrow)*u_scb(jl,jk,jrow)+&
          v_scb(jl,jk,jrow)*v_scb(jl,jk,jrow)) + rh_scb(jl,jk,jrow)
        vom_scb(jl,jk,jrow) = (vo_scb(jl,jk,jrow)+zcorio)*v_scb(jl,jk,jrow) +&
          vom_scb(jl,jk,jrow)
        vol_scb(jl,jk,jrow) = -(vo_scb(jl,jk,jrow)+zcorio)*u_scb(jl,jk,jrow) +&
           vol_scb(jl,jk,jrow)
        zdt = -u_scb(jl,jk,jrow)*dtl_scb(jl,jk,jrow) &
              - v_scb(jl,jk,jrow)*dtm_scb(jl,jk,jrow)
        te_scb (jl,jk,jrow) = te_scb(jl,jk,jrow) + zdt
      END DO
    END DO

    IF (ltdiag) THEN
       ! pressure gradient term, horizontal advection and coriolisterm
       pdiga(1:nglon,:, 1,jrow) = pdiga(1:nglon,:, 1,jrow) + vom_scb(1:nglon,:,jrow)
       pdiga(1:nglon,:, 6,jrow) = pdiga(1:nglon,:, 6,jrow) + vol_scb(1:nglon,:,jrow)
       ! horizontal advection
       pdiga(1:nglon,:,12,jrow) = pdiga(1:nglon,:,12,jrow) + te_scb (1:nglon,:,jrow)
       ! G-term, potential energy and kinetic energy
       pdiga(1:nglon,:,11,jrow) = pdiga(1:nglon,:,11,jrow) - rh_scb (1:nglon,:,jrow)
    END IF

!-- 10. Duplicate ps

    aps(:,jrow)  = aph(:,nlevp1)
    apsm(:,jrow) = aph(:,nlevp1)

!-- 9. Compute maximum !u!+!v!

    zrcst(jrow) = 1./sqcst(irow)

  END DO

  ul_scb(:,:) = (maxval_zonal(u_scb(:nglon,:,:))  + &
                 minval_zonal(u_scb(:nglon,:,:))) * 0.5

  ulz(:,:) = ul_scb(:,:)

  vmaxz(:,:) = SQRT(maxval_zonal(                           &
                      u_scb(:nglon,:,:)*u_scb(:nglon,:,:)   &
                    + v_scb(:nglon,:,:)*v_scb(:nglon,:,:))) &
                 * spread (zrcst,dim=1,ncopies=nlev)

END SUBROUTINE dyn

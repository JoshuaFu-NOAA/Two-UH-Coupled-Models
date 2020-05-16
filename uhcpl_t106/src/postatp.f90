!+ complete budgets for physics.
!+ $Id: postatp.f90,v 1.11 1999/11/10 08:21:43 m214003 Exp $

SUBROUTINE postatp

  ! Description:
  !
  ! Complete budgets for physics.
  !
  ! Method:
  !
  ! *postatp* is called from *scan1*.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, June 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tmp_buffer
  USE mo_parameters
  USE mo_control
  USE mo_diagnostics
  USE mo_param_switches
  USE mo_constants
  USE mo_start_dataset
  USE mo_stat_global
  USE mo_stat_zonal
  USE mo_doctor
  USE mo_year

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zgqm, zgqm0, zgsn, zgsn0, zgwd, zgwd0, zgws, zgws0, zinter, znorm, &
&      zpsgrho, zqday, zqsec, zsec, ztimdia
  INTEGER :: icd, id, iday, idiapfr, ihour, im, isec, ishour, iy, iyy

  !  Local arrays: 
  REAL :: zdef(8)

  !  Intrinsic functions 
  INTRINSIC NINT, SUM


  !  Executable Statements 

!-- 1. Set up some constants

  zpsgrho = gps/g/rhoh2o
  znorm   = 1000.
  ishour  = 3600

!-- 2. Complete diagnostics for first time step

  ! Complete global sums

  gts = SUM(gtsz(1:ngl))
  gtd = SUM(gtdz(1:ngl))
  gws = SUM(gwsz(1:ngl))
  gwd = SUM(gwdz(1:ngl))
  gsn = SUM(gsnz(1:ngl))

  IF (nstep==nstart) THEN

    gke0 = gke
    gpe0 = gpe
    gqm0 = gq*zpsgrho
    gts0 = gts
    gtd0 = gtd
    gws0 = gws
    gwd0 = gwd
    gsn0 = gsn

!-- 3. Complete diagnostics for other time steps

  ELSE

    gqm = gq*zpsgrho

    zgqm0 = gqm0*znorm
    zgsn0 = gsn0*znorm
    zgws0 = gws0*znorm
    zgwd0 = gwd0*znorm

    zgqm = gqm*znorm
    zgsn = gsn*znorm
    zgws = gws*znorm
    zgwd = gwd*znorm

    zdef(1) = gpe - gpe0 - (dia(1)+dia(2)+dia(3)+dia(5)+dia(7)+dia(9)+dia(15) &
&        +dia(17)+dia(22)+dia(24)+dia(26)+dia(27)+dia(29)+dia(33)+dia(35)+ &
&        dia(37)+dia(38)+dia(40)+dia(53))
    zdef(2) = gke - gke0 - (dia(8)+dia(16)+dia(18)+dia(54))
    zdef(3) = zgqm - zgqm0 - (dia(11)+dia(19)+dia(20)+dia(21)+dia(23)+dia(25) &
&        +dia(28)+dia(30)+dia(34)+dia(36)+dia(39)+dia(41))
    zdef(4) = gts - gts0 - (dia(4)+dia(6)+dia(10)+dia(12)+dia(44)+dia(48))
    zdef(5) = gtd - gtd0 - (dia(45)+dia(56))
    zdef(6) = zgsn - zgsn0 - (dia(13)+dia(31)+dia(42)+dia(49))
    zdef(7) = zgws - zgws0 - (dia(14)+dia(32)+dia(43)+dia(46)+dia(50)+dia(51) &
&        )
    zdef(8) = zgwd - zgwd0 - (dia(55)+dia(47)+dia(52))

    ztimdia = nstep*0.5*twodt
    zsec  = nstep*0.5*twodt
    iday  = NINT(zsec/dayl)
    ihour = NINT((zsec-(iday*dayl))/ishour)
    isec  = NINT(zsec-ihour*ishour-iday*dayl)

    icd = ((nstep)*dtime+ntbase)/dayl + ncbase
    CALL cd2dat(icd,id,im,iy)
    iyy = iy - 0
    IF (nstep==nstart) THEN
      idiapfr = 1
    ELSE
      idiapfr = ndiapfr
    END IF

!-- 4. Print diagnostics

    zinter = idiapfr*0.5*twodt/ishour

    zqsec = 1./(ndiapfr*dtime)
    zqday = zqsec*dayl

    WRITE (nout,'(a,a,f5.1,a)') &
         ' Global diagnostics for energy and water cycles.', &
         ' Diagnostic interval: ', zinter , ' hours.'
    WRITE (nout, '(a,i2.2,a,i2.2,a,i4.4,a,i4,a,i7)') &
       ' Date: ', id, '.', im, '.',iy, &
       '  Forecast year: ', iyy, ', nstep = ', nstep
    WRITE (nout,'(1x,135("-"))')
    WRITE (nout,'(a,i4,a,i2,7x,a,14x,a,15x,a,26x,a,27x,a)') &
         ' |  Day',iday,' Hour ',ihour,'||','Atmosphere','| |', &
         'Land surface','|'
    WRITE (nout,'(a,f8.1,10x,110("-"))') &
      ' | sec. ', zsec
    WRITE (nout,'(a,13x,a,a,a)') &
      ' |', 'sub.=2.8345|| potential  |   kinetic  ', &
      '||water vapour| | top layer  | deep layer || equivalent ', &
      '| top layer  | deep layer |'
    WRITE (nout,'(a,a,a)') &
         ' |latent heats vap.=2.5008||   energy   |   energy   ||', &
         ' equivalent | |   energy   |   energy   || snow depth |', &
         '   water    |   water    |'
    WRITE (nout,'(a,13x,a,a,a)') &
         ' |', 'fus.=0.3337||  [W/m**2]  |  [W/m**2]  ||  [mm/day]  ', &
         '| |  [W/m**2]  |  [W/m**2]  ||  [mm/day]  |  [mm/day]  ', &
         '|  [mm/day]  |'
    WRITE (nout,'(1x,135("-"))')
    WRITE (nout,'(1x,135("-"))')
    WRITE (nout,'(a,8(f12.3,a))') &
         ' |     initial values     ||', gpe0*zqsec, '|', &
         gke0*zqsec, '||', zgqm0*zqday, '| |', gts0*zqsec, '|', &
         gtd0*zqsec, '||', zgsn0*zqday, '|', zgws0*zqday, '|', &
         zgwd0*zqday , '|'

    WRITE (nout,'(1x,135("-"))')
    WRITE (nout,'(1x,135("-"))')
    
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,12x,a,12x,a,12x,a,3(12x,a))') & 
         ' |','|top short wave  ||',dia(1)*zqsec, &
         '|', '||', '| |', '|', '||','|','|','|'
    WRITE (nout,'(a,3x,a,f12.3,a,12x,a,12x,a,12x,a,12x,a,3(12x,a))') &
         ' |radiat.|top long wave','||',dia(2)*zqsec, &
         '|','||','| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,12x,a,f12.3,a,12x,a,3(12x,a))') & 
         ' |','|surf. short wave||',dia(3)*zqsec,'|','||','| |', &
         dia(4)*zqsec,'|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,12x,a,f12.3,a,12x,a,3(12x,a))') &
         ' |','|surf. long wave ||',dia(5)*zqsec,'|','||', &
         '| |',dia(6)*zqsec,'|','||','|','|','|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,7x,a,f12.3,a,f12.3,a,12x,a,12x,a,12x,a,3(12x,a))') &
         ' |','|friction        ||',dia(7)*zqsec,'|',dia(8)*zqsec,'||', &
         '| |','|','||','|','|','|'
    WRITE (nout,'(a,f12.3,a,12x,a,12x,a,f12.3,a,12x,a,3(12x,a))') &
         ' |turbul.|sensible heat   ||',dia(9)*zqsec,'|','||', &
         '| |',dia(10)*zqsec,'|','||','|','|','|'
    WRITE (nout,'(a,7x,a,12x,a,12x,a,f12.3,a,f12.3,a,12x,a,3(f12.3,a))') &
         ' |','|latent heat     ||','|','||',dia(11)*zqday,'| |', &
         dia(12)*zqsec,'|','||',dia(13)*zqday,'|',dia(14)*zqday,'|', &
         dia(55)*zqday,'|'
    WRITE (nout,'(a,7x,a,f12.3,a,f12.3,a,12x,a,12x,a,12x,a,3(12x,a))') &
         ' |','|gwd dissipation ||',dia(15)*zqsec,'|',dia(16)*zqsec, &
         '||','| |','|','||','|','|','|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,7x,a,f12.3,a,f12.3,a,12x,a,12x,a,12x,a,3(12x,a))') &
         ' |','|friction        ||',dia(17)*zqsec,'|',dia(18)*zqsec,'||','| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,12x,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))') &    
         ' |','|moisture supply ||','|','||',dia(19)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,12x,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))') &
         ' |','|env. moistening ||','|','||',dia(20)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,12x,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))') &
         ' |','|prec. formation ||','|','||',dia(21)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,f12.3,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))')  &
         ' |convec.|rain generation ||',dia(22)*zqsec,'|','||',dia(23)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))') &
         ' |','|snow generation ||',dia(24)*zqsec,'|','||',dia(25)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,12x,a,12x,a,12x,a,3(12x,a))') &
         ' |','|snow/rain melt  ||',dia(26)*zqsec,'|','||','| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))') &
         ' |','|rain evaporation||',dia(27)*zqsec,'|','||',dia(28)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,12x,a,12x,a,3(12x,a))') &
         ' |','|snow evaporation||',dia(29)*zqsec,'|','||',dia(30)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,12x,a,12x,a,12x,a,12x,a,12x,a,2(f12.3,a),12x,a)') &
         ' |','|surf. precipit. ||','|','||','| |','|','||',dia(31)*zqday,'|',dia(32)*zqday,'|','|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,5(12x,a))') & 
         ' |','|rain generation ||',dia(33)*zqsec,'|','||', &
         dia(34)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,5(12x,a))') & 
          ' |','|snow generation ||',dia(35)*zqsec,'|','||', &
         dia(36)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,f12.3,a,7(12x,a))') & 
         ' |l.s.pr.|snow/rain melt  ||',dia(37)*zqsec,'|', &
         '||','| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,5(12x,a))') & 
         ' |','|rain evaporation||',dia(38)*zqsec,'|','||', &
         dia(39)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,f12.3,a,12x,a,f12.3,a,5(12x,a))') & 
         ' |','|snow evaporation||',dia(40)*zqsec,'|','||', &
         dia(41)*zqday,'| |','|','||','|','|','|'
    WRITE (nout,'(a,7x,a,5(12x,a),2(f12.3,a),12x,a)') & 
         ' |','|surf. precipit. ||','|','||','| |','|','||', &
         dia(42)*zqday,'|',dia(43)*zqday,'|','|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,7x,a,3(12x,a),2(f12.3,a),12x,a,2(f12.3,a))') & 
         ' |','|exch. in soil   ||','|','||','| |',dia(44)*zqsec, &
         '|',dia(45)*zqsec,'||','|',dia(46)*zqday,'|',dia(47)*zqday,'|'
    WRITE (nout,'(a,3(12x,a),2(f12.3,a),12x,a,2(f12.3,a),12x,a)') & 
         ' |soil p.|snow melt       ||','|','||','| |',dia(48)*zqsec, &
         '|',dia(56)*zqsec,'||','|',dia(49)*zqday,'|',dia(50)*zqday,'|'
    WRITE (nout,'(a,7x,a,6(12x,a),2(f12.3,a))') & 
         ' |','|run-off         ||','|','||','| |','|','||','|', &
         dia(51)*zqday,'|',dia(52)*zqday,'|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,f12.3,a,f12.3,a,6(12x,a))') & 
         ' |adiabatic conversion    ||',dia(53)*zqsec,'|', &
         dia(54)*zqsec,'||','| |','|','||','|','|','|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,8(f12.3,a))') & 
         ' |residual terms          ||',zdef(1)*zqsec,'|',zdef(2)*zqsec, &
         '||',zdef(3)*zqday,'| |',zdef(4)*zqsec,'|',zdef(5)*zqsec,'||', &
         zdef(6)*zqday,'|',zdef(7)*zqday,'|',zdef(8)*zqday,'|'

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(1x,135("-"))')

    WRITE (nout,'(a,8(f12.3,a))') & 
         ' |     present values     ||',gpe*zqsec,'|',gke*zqsec, &
         '||',zgqm*zqday,'| |',gts*zqsec,'|',gtd*zqsec,'||', &
         zgsn*zqday,'|',zgws*zqday,'|',zgwd*zqday,'|'

    WRITE (nout,'(1x,135("-"),/)')

    ! New initial values

    gpe0 = gpe
    gke0 = gke
    gqm0 = gqm
    gts0 = gts
    gtd0 = gtd
    gws0 = gws
    gwd0 = gwd
    gsn0 = gsn

  END IF

  RETURN
END SUBROUTINE postatp

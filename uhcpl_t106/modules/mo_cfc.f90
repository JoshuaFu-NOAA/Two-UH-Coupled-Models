MODULE mo_cfc

  IMPLICIT NONE

  !      ------------------------------------------------------------------
  ! *    lw-absorption coefficients of cfc's
  !      ------------------------------------------------------------------

  INTEGER :: cfcwmo(4,16) !  lw-absorption coefficients in 4 lw-intervals
  INTEGER :: zdicfc       !  diffusivity factor for cfc's
  INTEGER :: ncfc         !  number of cfc's (=2 -> cfc11 and cfc12 )

  INTEGER :: JC
  !        ABSORPTION COEFFICIENTS CFC'S, HCFC'S AND HFC'S IN 4 LW
  !        FREQUENCY INTERVALS (IN CM-1.ATM-1)
  !        DATA EXTRACTED FROM INTEGRATED BAND INTENSITIES,
  !        AFEAS REPORT, WMO 1989

  !
  !*    1. SET VALUES
  !     -------------
  DATA NCFC, ZDICFC /16, 2./
  !
  !
  !*    2. ABSORPTION COEFFICIENTS (CM-1.ATM-1) IN 4 LW INTERVALS
  !     ---------------------------------------------------------
  !
  !---- INTERVAL 2 --- 500-800 CM-1 --- FROM WMO 1989 ---
  DATA (CFCWMO(1,JC),JC=1,16)/  &
       & 0.121004E+01,   0.000000E+00,   0.521869E+00,   0.417986E+00, &
       & 0.357761E+00,   0.409101E+00,   0.844864E+00,   0.848554E+00, &
       & 0.677562E+00,   0.497803E+00,   0.216235E+01,   0.616038E+00, &
       & 0.302616E+00,   0.727887E-01,   0.210073E+01,   0.246229E+01/
  !
  !---- INTERVAL 3 --- 800-970 + 1110-1250 CM-1 --- FROM WMO 1989 ---
  DATA (CFCWMO(2,JC),JC=1,16)/  &
       & 0.539411E+01,   0.858383E+01,   0.803155E+01,   0.887493E+01, &
       & 0.941552E+01,   0.525779E+01,   0.542247E+01,   0.787479E+01, &
       & 0.932012E+01,   0.472236E+01,   0.280806E+01,   0.595839E+01, &
       & 0.618597E+01,   0.326485E+01,   0.178601E+01,   0.817861E+00/
  !
  !---- INTERVAL 4 --- 970-1110 CM-1 --- FROM WMO 1989 ---
  DATA (CFCWMO(3,JC),JC=1,16)/  &
       & 0.102980E+01,   0.359796E+01,   0.488178E+01,   0.701477E+01, &
       & 0.742985E+01,   0.338980E+01,   0.189000E+01,   0.284286E+01, &
       & 0.225714E+01,   0.385562E+01,   0.211261E+01,   0.294679E+01, &
       & 0.243261E+01,   0.270085E+01,   0.000000E+00,   0.122632E+01/
  !
  !--- INTERVAL 6 --- 1250-1450 + 1880-2820 CM-1 --- FROM WMO 1989 ---
  DATA (CFCWMO(4,JC),JC=1,16)/  &
       & 0.000000E+00,   0.000000E+00,   0.155309E-01,   0.176863E+00, &
       &  0.486060E+00,   0.276579E+00,   0.571606E+00,   0.800380E+00, &
       & 0.436065E+00,   0.983304E+00,   0.641228E-01,   0.117235E+00, &
       & 0.795934E+00,   0.201075E+00,   0.000000E+00,   0.280702E-01/

END MODULE mo_cfc

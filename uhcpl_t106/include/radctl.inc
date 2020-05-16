!*    namelist *RADCTL* to control options having a bearing on the      
!                       radiation scheme.                               
!                                                                       
!      name     type      purpose                                       
!      ----     ----      -------                                       
!                                                                       
!     *lrad*    LOGICAL   *true for radiation.                          
!     *ldiur*   LOGICAL   *true for diurnal cycle.                      
!     *lsolc*   LOGICAL   *true for solar clear sky diagnostic.         
!     *laer*    LOGICAL   *true for aerosols                            
!     *lgadsrh* LOGICAL   *true for rel.hum. dependency                 
!                          of *gads* aerosol optical parameters         
!     *nradfr*  INTEGER   *frequency of full radiation computations.    
!                         if(nradfr.gt.0):rad every *nradfr* time step. 
!                         if(nradfr.lt.0):rad every *-nradfr* hours.    
!     *nrpart*  INTEGER   *partition factor for radiation computations. 
!     *nradpfr* INTEGER   *print frequency for radiation statistics     
!                         (in number of radiation steps).               
!     *nradpla* INTEGER   *print radiation statistics every *nradpla*   
!                         latitude line.                                
!     *ndfaer*  INTEGER   aerosol definition array (see *suaerx*)       
!                                                                       

NAMELIST /radctl/                                                      &
& lrad, ldiur, lsolc, laer, lcfc, lgadsrh, co2fac, nmonth, nradfr,     &
& nradpfr, nradpla ,ndfaer

!     ------------------------------------------------------------------

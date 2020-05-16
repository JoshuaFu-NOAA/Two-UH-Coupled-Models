!*    namelist *PHYSCTL* to control options having a bearing on the     
!                        parameterisations of diabatic processes.       
!                                                                       
!      name     type      purpose                                       
!      ----     ----      -------                                       
!                                                                       
!     *lphys*   LOGICAL   *true for parameterisation of diabatic        
!                          processes.                                    
!     *lvdiff*  LOGICAL   *true for vertical diffusion.                 
!     *lconv*   LOGICAL   *true* to allow convection                   
!     *lmfpen*  LOGICAL   *true if penetrative convection is switched on 
!     *lcond*   LOGICAL   *true for large scale condensation scheme.    
!     *lgwdrag* LOGICAL   *true for gravity wave drag scheme  
!     *lsurf*   LOGICAL   *true for surface exchanges.                  
!     *lice*    LOGICAL   *true* for sea-ice temperature calculation    
!                                                                       
!     *ndiapfr* INTEGER   *frequency of physics budgets.                
!                          if(nradfr.gt.0):rad every *nradfr* time step. 
!                          if(nradfr.lt.0):rad every *-nradfr* hours.    
!                                                                       

NAMELIST /physctl/                                       &
& lphys, lvdiff, lcond, lsurf, lconv, lmfpen, lgwdrag, lice, ndiapfr

!-----------------------------------------------------------------------------

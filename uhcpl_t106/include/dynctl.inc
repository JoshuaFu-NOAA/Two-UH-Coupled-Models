!*    namelist *DYNCTL* to control options having a bearing on the dy-  
!                       namics and the general control of the forecast. 
!                                                                       
!      name     type      purpose                                       
!      ----     ----      -------                                       
!                                                                       
!     *ndiadfr* INTEGER   *frequency of dynamical diagnostics.          
!                         if(ndiadfr.gt.0):every ndiadfr time step.     
!                         if(ndiadfr.eq.0):no diagnostic.               
!                         if(ndiadfr.lt.0):every -ndiadfr hour.         
!     *ndiavfr* INTEGER   *frequency of vertical dynamical diagnostics. 
!                         if(ndiavfr.gt.0):every ndiavfr time step.     
!                         if(ndiavfr.eq.0):no diagnostic.               
!                         if(ndiavfr.lt.0):every -ndiavfr hour.         
!     *ntrn*    INTEGER   *critical wavenumbers for strong              
!                         stratospheric damping.                        
!     *nlvstd1* INTEGER   *LAST (UPPERMOST) LAYER AT WHICH              
!                         STRATOSPHERIC HORIZONTAL DIFFUSION IS         
!                         ENHANCED.                                     
!     *NLVSTD2* INTEGER   *first (lowest) layer at which                
!                         stratospheric horizontal diffusion is         
!                         enhanced.                                     
!     *ldrag*   LOGICAL   .true. to apply drag to upper 2 levels        
!     *vcrit*   REAL      *critical velocity above which                
!                         horizontal diffusion is enhanced for          
!                         t63 with dt=20min.                            
!     *vcheck*  REAL      threshold value for check of high windspeed   
!     *hdamp*   REAL      *damping factor for strong                    
!                         stratospheric damping.                        
!     *enstdif* REAL      *factor by which stratospheric                
!                         horizontal diffusion is increased from one    
!                         level to next level above.                    
!     *cdrag*   REAL      drag coefficient in seconds                   
!                         or -(drag coefficent) in days                 
!     *lumax*   LOGICAL   *true to compute and print information on     
!                         maximum wind.                                 
!     *eps*     REAL      *coefficient for the time filtering.          
!     *dampth*  REAL      *damping time (n=nk) in hours for the         
!                         horizontal                                    
!                         diffusion of vorticity (linear square         
!                         laplacian), divergence, temperature.          
!     *apsurf*  REAL      fixed global mean of surface pressure         
!         
!     *damhih*  REAL                                                    
!     *spdrag*  REAL      upper sponge layer coefficient (sec)-1
!     *enspodi* REAL      factor by which   upper sponge layer 
!                         coefficient is increased from one        
!                         level to next level above. 
!     *nlvspd1* INTEGER   last (uppermost) layer of upper sponge  
!     *nlvspd2* INTEGER   first (lowest) layer of upper sponge  

NAMELIST /dynctl/                                                      &
& ndiadfr, ndiavfr, ntrn, nlvstd1, nlvstd2, lumax, lzondia, ldrag,     &
& ldiahdf, vcrit, hdamp, enstdif, apsurf, vcheck, eps, dampth, cdrag,  &
& damhih, spdrag, enspodi, nlvspd1, nlvspd2

!     ------------------------------------------------------------------

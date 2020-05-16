!**** *namelist* *SDSCTL* - modify start data set default values.       
!                                                                       
!     J.K. Gibson, ECMWF, December 1982, original source
!                                                                       
!      name       type      purpose                                       
!      ----       ----      -------                                       
!     *nstart*    INTEGER   time step for start/restart.                  
!     *nstep*     INTEGER   current time step.                            
!     *nf1a*                                                            
!      *to*       INTEGER   indices for fourier buffers.                  
!     *nf4a*                                                            
!     *ng1a*                                                            
!      *to*       INTEGER   indices for grid point buffers.               
!     *ng4a*                                                            
!     *ngl1a*     INTEGER   index for grid point slt buffer.              
!     *nl1a*      INTEGER   index for legendre coefficients buffer        
!     *ndiahdf*   INTEGER   logical unit for hdiff diagnostics.           
!     *nhf1*                                                            
!      *to*       INTEGER   logical units for fourier history files.      
!     *nhf4*                                                            
!     *nhg1*                                                            
!      *to*       INTEGER   logical units for grid point history files.   
!     *nhg4*                                                            
!     *nhgl1*     INTEGER   logical unit for grid point slt work file     
!     *nhg4*                                                            
!     *nisp*      INTEGER   logical unit  for initial spectral fields.    
!     *nigp*      INTEGER   logical unit  for initial grid point fields.  
!     *nist*      INTEGER   logical unit for surf.temp. file              
!     *nemi*      INTEGER   logical unit for surface emission file        
!     *nini*      INTEGER   logical unit for tracer initial file          
!     *nfl1*      INTEGER   logical unit for optional file read at nstep=0
!     *nfl2*      INTEGER   logical unit for optional file read at nresum 
!     *ngribs*    INTEGER   logical unit for spectral grib file           
!     *ngribg*    INTEGER   logical unit for gridpoint grib file          
!     *ngribx*    INTEGER   logical unit for g4x grib file                
!     *njin*      INTEGER   logical unit for "jobn" input file            
!     *njout*     INTEGER   logical unit for "subjobn" output file        
!      ndstart    INTEGER   start date of a run
!      ntstart    INTEGER   start time of a run
!     *ntimeadj*  INTEGER   time offset adjustment in seconds
!     *lres*      LOGICAL   .true. if forecast is restarted.              
!     *ldebugio*  LOGICAL   .true. to debug io           
!     *ldebugmem* LOGICAL   .true. to debug memory
!     *ly365*     LOGICAL   .true. for runs with 365 days per year
!     *lg3force*  LOGICAL   .true. to force read from g3b if g3a not in restart
!     *lg3setzero LOGICAL   .true. to set G3 vars to 0. if not in restart

NAMELIST /sdsctl/                                                      &
& nstart, nstep, nf1a, ng1a, ngl1a, nl1a, nisp, nigp, ndiahdf,         &
& nemi, nist, nhf1, nhf3, nhf4, nhg1, nhg2, nhg3, ntimeadj,            &
& nhg4, nhgl1, nfl1, nfl2, nini, ngribs, ngribg, ngribx, njin, njout,  &
& lres, ldebugio, ldebugmem, ly365, lg3force, lg3setzero,              &
& ndstart, ntstart 

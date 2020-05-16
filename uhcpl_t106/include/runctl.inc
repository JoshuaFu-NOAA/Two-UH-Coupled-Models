!*    namelist *RUNCTL* controls general options to run the model.      
!                                                                       
!      name       type      purpose                                       
!      ----       ----      -------                                       
!                                                                       
!     *nresum*    INTEGER   time step at which the run started or was     
!                           reumed after interruption.                    
!     *ncbase*    INTEGER   century date of the initial data.             
!     *ntbase*    INTEGER   time of the initial data.                     
!     *nstop*     INTEGER   time step to stop the forecast.               
!                           a negative value indicates number of          
!                           forecast hours required, and will be          
!                           converted to time steps after the namelist    
!                           has been read.                                
!     *nwrite*    INTEGER   time step of next write-up time.              
!     *nwtime*    INTEGER   *frequency of write up times.                 
!     *nwlag*     INTEGER   frequency of rerun-files saving               
!     *nsub*      INTEGER   number of jobs to be submitted at end of run  
!     *nsubint*   INTEGER   frequency of submits in months                
!     *nptime*    INTEGER   frequency of write up times                   
!     *nctime*    INTEGER   frequency of coupling with ocean
!     *n4ptime*   INTEGER   frequency of write up times of g4x fields     
!     *numfl1*    INTEGER   number of optional fields read at nstep=0     
!     *numfl2*    INTEGER   number of optional fields read at nstep=nresum
!     *dtime*     REAL      time step in seconds.                         
!     *ldebug*    LOGICAL   .true. for mass fixer diagnostics.            
!     *lamip*     LOGICAL   .true. for *amip* type experiment             
!     *ng3xp*     INTEGER   number of g3x fields                          
!     *ng4xp*     INTEGER   number of g4x fields                          
!     *ng3xl*     INTEGER   number of g3x levels (array)                  
!     *ng4xl*     INTEGER   number of g4x levels (array)                  
!     *nxpbits*   INTEGER   number of bits of g3x for grib format (array) 
!     *nx4pbits*  INTEGER   number of bits of g4x for grib format (array) 
!
!     *lxaccu*    LOGICAL   .true. for zero reset g3x after write up (array) 
!     *l4xaccu*   LOGICAL   .true. for zero reset g4x after write up (array) 
!     *lhg3x*     LOGICAL   .true. to use g3x info from history files     
!     *lsstadj*   LOGICAL   .true. for orographic adjustment of sst       
!     *lsub*      LOGICAL   .true. to submit *nsub* jobs                  
!     *lrepro *   LOGICAL   .true. for reproducable results in multi task 
!                                                                       
!     *lnwp*      LOGICAL   .false. for climate mode, .true. for nwp mode
!     *lanalysis* LOGICAL   .false. for analysis mode, .false. usually
!     *lnudge*    LOGICAL   .true.  for nudging mode
!      lmidatm    LOGICAL   .true.  for middle atmosphere model version
!
!     *lnmi*      LOGICAL   .true.  use normal mode initialisation
!     *ltdiag*    LOGICAL   .true.  additional diagnostics of tendencies
!
!     *lpci*      LOGICAL   .true.  for pci version
!
!     *lcolumn    LOGICAL   .true.  to run the column model
!
!     *nproca*    INTEGER   number of processors for set A
!     *nprocb*    INTEGER   number of processors for set B	
!     *lamip2*    LOGICAL   .false.  true for *amip2* type experiments
!
!     *lalai*     LOGICAL   .true.  for annual lai
!     *lavgrat*   LOGICAL   .true.  for annual vegetation
!     *imtt*      (1,1000) fu++
!
!  add options for ensemble generation
!     *lens*      LOGICAL   .true.  generate rerun file ensemble
!     *nens*      INTEGER   number of ensemble size
!fu++ *lonudg*    LOGICAL   .true. nudging ocean sst from nudging data
!fu++ *ldsst*     LOGICAL   .true. turn on SST diurnal cycle
!fu++ *lstratiform*     LOGICAL   .true. turn on stratiform forcing 
 
NAMELIST /runctl/                                                      &
& nresum, ncbase, ntbase, ng3xp, ng3xl, nxpbits, ng4xp, ng4xl,         &
& n4xpbits, n4ptime, nstop, nwtime, nwlag,                             &
& nptime, nctime, nsub, nsubint, numfl1, numfl2, dtime,                &
& lamip2,lonudg,ldsst,lstratiform,                                     &
& lnwp, lanalysis, lnudge, lmidatm,                                    &
& lxaccu, l4xaccu, lcouple, lamip, lsstadj, labort, lsub,              & 
& lhg3x, ldebug, lrepro, lnmi, ltdiag, lpci, nproca, nprocb,           &
& lalai, lavgrat, lcolumn,imtt,                                        &
& lens, nens

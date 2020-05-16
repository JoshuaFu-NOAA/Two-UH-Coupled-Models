export

SHELL = /bin/sh

SYSTEM = $(shell uname -msr|sed -e 's= =_=g')

 
FC       = f90
CC       = cc


ifeq ($(SYSTEM),SunOS_5.6_sun4u)
CC       = gcc 
CFLAGS   = -O
FFLAGS   = 
endif

ifeq ($(SYSTEM),SUPER-UX_8.1_SX-4)
CFLAGS   = -copp,transform -Wo,"-IExtraZ,-IExtraT" -hfloat0,ansi,vector \
           -Onooverlap -pvctl,vworksz=4M,fullmsg,noassume,loopcnt=1000000 -D_BUILTIN_
FFLAGS   = -float0
endif

ifeq ($(SYSTEM),sn9114_9.0.2.5_CRAY_J90)
CFLAGS   = -O2 -h ivdep -h report=v
FFLAGS   = 
endif

ifeq ($(SYSTEM),sn4001_9.0.2.5_CRAY_C90)
CFLAGS   = -O2 -h ivdep -h report=v
FFLAGS   = 
endif

ifeq ($(SYSTEM),UNIX_System_V_4.1_ES_F700)
CC       = vcc
CFLAGS   = -O -K4,lib -Wv,-m3 -xExtraZ,ExtraT
FFLAGS   = 
endif

ifeq ($(SYSTEM),IRIX64_6.2_IP19)
CFLAGS   = -O -OPT:Olimit=9999
FFLAGS   = 
endif

ifeq ($(SYSTEM),IRIX_6.2_IP22)
CC       = gcc
CFLAGS   = -O -Olimit 9999 -sopt
CFLAGS   = -O2
FFLAGS   = 
endif

SRCDIRS = burn table/src service grads/gradsdes

all:
	@echo `pwd`
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done       

clean:
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done


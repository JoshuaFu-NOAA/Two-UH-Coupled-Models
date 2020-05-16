include Make.config 

SRCDIRS = src initial sst ozon fld restart
INCLUDEDIR = include
.PHONY : src initial restart sst ozon fld

all:
	@echo `pwd`
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done 

#%:	src    # doesn't function !!
#	@echo "testing if $(%) is a target .." ; \
#	@if [ -d $% ] ; then echo "making $%" ; \
#	cd $% &&  $(MAKE) ; fi

clean:
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
# Rules for singel target
src:
	cd src &&  $(MAKE)
initial: src
	cd initial &&  $(MAKE)
restart: src
	cd restart &&  $(MAKE)
ozon: src
	cd ozon &&  $(MAKE)
sst: src
	cd sst &&  $(MAKE)
fld: src
	cd fld &&  $(MAKE)

# Rules for cleaning singel target

srcclean :
	cd src && $(MAKE) clean
iniclean : initialclean
initialclean :
	cd initial && $(MAKE) clean
resclean : restartclean
restartclean :
	cd restart && $(MAKE) clean
sstclean :
	cd sst && $(MAKE) clean
ozonclean :
	cd ozon && $(MAKE) clean
fldclean :
	cd fld $(MAKE) clean
indent: indent_include
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) indent ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make indent  was $$?" ; exit 1 ; fi ;) ;\
	  done 

indent_include:
	$(INDENT) $(INCLUDEDIR)/*.h

editclean:
	@for DIR in $(SRCDIRS) ;
	  do \
	  (cd $$DIR ;\
	  $(MAKE) editclean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
distclean:
	@for DIR in $(SRCDIRS) ;\
	do \
	  (cd $$DIR ;\
	  $(MAKE) distclean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

include  ../Make.config
PROG =	../bin/res2rnc

SRCS =	header.c netcdf_interface.c  read-data.c \
	restart-convert.c ncheader.c

OBJS =	header.o  netcdf_interface.o read-data.o  \
	restart-convert.o ncheader.o

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(PROG) $(OBJS)

editclean:
	-rm -f *.c~ *.h~

indent:
	$(INDENT) *.c 

distclean: indent editclean 

header.o: ../include/prep.h 
netcdf_interface.o: ../include/prep.h 
ncheader.o: ../include/prep.h 
read-data.o: ../include/prep.h 
restart-convert.o: ../include/prep.h

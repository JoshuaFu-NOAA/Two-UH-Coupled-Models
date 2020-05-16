include ../Make.config

PROG1 =	../bin/year2ync 
PROG2 = ../bin/ync2srv 
PROG3 = ../bin/srv2ync
PROG4 = ../bin/sst2srv
PROG5 = ../bin/srv2sst

SRCS1 =	year_table.c year_io.c year_ncio.c year2ync.c
SRCS2 =	year_table.c year_io.c year_ncio.c ync2srv.c
SRCS3 =	year_table.c year_io.c year_ncio.c srv2ync.c
SRCS4 =	year_io.c sst2srv.c
SRCS5 =	year_io.c srv2sst.c

OBJS1 := $(SRCS1:.c=.o)
OBJS2 := $(SRCS2:.c=.o)
OBJS3 := $(SRCS3:.c=.o)
OBJS4 := $(SRCS4:.c=.o)
OBJS5 := $(SRCS5:.c=.o)

all: $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5)

$(PROG1): $(OBJS1)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS1) $(LIBS)

$(PROG2): $(OBJS2) 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS2) $(LIBS)

$(PROG3): $(OBJS3)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS3) $(LIBS)

$(PROG4): $(OBJS4)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS4) $(LIBS)

$(PROG5): $(OBJS5)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS5) $(LIBS)

clean:
	-rm -f $(PROG1) $(OBJS1) $(PROG2) $(OBJS2)  $(PROG3) $(OBJS3) \
	 $(PROG4) $(OBJS4)  $(PROG5) $(OBJS5)

editclean:
	-rm -f *.c~ *.h~ 

indent:
	$(INDENT) *.c 

distclean: indent editclean 



year_io.o: ../include/prep.h ../include/io_nc.h
ync2srv.o: ../include/prep.h
srv2ync.o: ../include/prep.h
sst2ync.o: ../include/prep.h

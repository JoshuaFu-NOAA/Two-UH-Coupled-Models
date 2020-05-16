include ../Make.config

PROG1 = ../bin/inispec2anc
PROG2 = ../bin/inisurf2snc
PROG3 = ../bin/snc2srv
PROG4 = ../bin/srv2snc

SRCS1 =	ini_table.c ini_io.c header.c inispec2anc.c
SRCS2 =	ini_table.c ini_io.c header.c inisurf2snc.c
SRCS3 = ini_table.c ini_io.c snc2srv.c
SRCS4 = ini_table.c ini_io.c srv2snc.c

OBJS1 := $(SRCS1:.c=.o)
OBJS2 := $(SRCS2:.c=.o)
OBJS3 := $(SRCS3:.c=.o)
OBJS4 := $(SRCS4:.c=.o)


all:  $(PROG1) $(PROG2) $(PROG3)  $(PROG4)

$(PROG1): $(OBJS1)
	$(CC)  $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS1) $(LIBS)
$(PROG2): $(OBJS2)
	$(CC)  $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS2) $(LIBS)
$(PROG3): $(OBJS3)
	$(CC)  $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS3) $(LIBS)
$(PROG4): $(OBJS4)
	$(CC)  $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS4) $(LIBS)


clean:
	-rm -f $(PROG1) $(OBJS1) $(PROG2) $(OBJS2) $(PROG3) $(OBJS3) $(PROG4) $(OBJS4)

editclean:
	-rm -f *.c~ *.h~

indent:
	$(INDENT) *.c 

distclean: indent editclean 


header.o: ../include/prep.h  ../include/pp.h
ini_io.o: ../include/prep.h ../include/io_nc.h ../include/pp.h
inispec2anc.o: ../include/pp.h ../include/prep.h
inisurf2snc.o: ../include/pp.h ../include/prep.h
snc2srv.o: ../include/pp.h ../include/prep.h
srv2snc.o: ../include/pp.h ../include/prep.h

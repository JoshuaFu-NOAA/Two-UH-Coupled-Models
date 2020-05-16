include ../Make.config

PROG1 =	../bin/fld2srv 
PROG2 = ../bin/fnc2srv 
PROG3 = ../bin/srv2fnc

SRCS1 =	fld_table.c fld_io.c fld2srv.c
SRCS2 =	fld_table.c fld_io.c fnc2srv.c
SRCS3 =	fld_table.c fld_io.c srv2fnc.c

OBJS1 =	$(SRCS1:.c=.o)
OBJS2 =	$(SRCS2:.c=.o)
OBJS3 =	$(SRCS3:.c=.o)

all: $(PROG1)   $(PROG2) $(PROG3)

$(PROG1): $(OBJS1)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS1) $(LIBS)

$(PROG2): $(OBJS2) 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS2) $(LIBS)

$(PROG3): $(OBJS3)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS3) $(LIBS)

clean:
	-rm -f $(PROG1) $(OBJS1) $(PROG2) $(OBJS2)  $(PROG3) $(OBJS3)

editclean:
	-rm -f *.c~ *.h~

indent:
	$(INDENT) *.c 

distclean: indent editclean 


fld_io.o: ../include/prep.h ../include/io_nc.h ../include/pp.h
fnc2srv.o: ../include/prep.h ../include/pp.h ../lib/libpp.a
srv2fnc.o: ../include/prep.h ../include/pp.h ../lib/libpp.a
fld2srv.o: ../include/prep.h ../include/pp.h ../lib/libpp.a


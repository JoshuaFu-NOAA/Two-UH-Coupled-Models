include ../Make.config

PROG1 =	../bin/ozon2srv 
PROG2 = ../bin/onc2srv 
PROG3 = ../bin/srv2onc

SRCS1 = ozon_table.c ozon_io.c ozon2srv.c
SRCS2 =	ozon_table.c ozon_io.c onc2srv.c
SRCS3 =	ozon_table.c ozon_io.c srv2onc.c

OBJS1 =	$(SRCS1:.c=.o)
OBJS2 =	$(SRCS2:.c=.o)
OBJS3 =	$(SRCS3:.c=.o) 

all: $(PROG1)  $(PROG2)   $(PROG3)

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


ozon_io.o: ../include/prep.h ../include/io_nc.h
onc2srv.o: ../include/prep.h
srv2onc.o: ../include/prep.h 
ozon2srv.o: ../include/prep.h

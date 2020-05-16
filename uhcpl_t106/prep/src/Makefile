include ../Make.config

LIBRARY = ../lib/libpp.a

PROG = ../bin/srv

SRCS =	dim_table.c io_nc.c pp.c

OBJS := $(SRCS:.c=.o)

all: $(LIBRARY) $(PROG)

$(LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS)

$(PROG): srv.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ srv.o $(LIBS)

pp.o: ../include/pp.h

clean:
	-rm -f $(OBJS) $(LIBRARY) srv.o $(PROG)

editclean:
	-rm -f *.c~ *.h~ ;
	-rm -f ../include/*.h~ 

indent:
	$(INDENT) $(SRCS)
	$(INDENT) ../include/*.h

distclean: indent  editclean

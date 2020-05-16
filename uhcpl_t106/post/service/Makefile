PROG1 = ../bin/read_srv
PROG2 = ../bin/read_srv8

SRCS1 = read_srv.f90
SRCS2 = read_srv8.f90

OBJS1 =	read_srv.o
OBJS2 =	read_srv8.o

all: $(PROG1) $(PROG2)

$(PROG1): $(OBJS1)
	$(FC) $(FFLAGS) -o $@ $(OBJS1)

$(PROG2): $(OBJS2)
	$(FC) $(FFLAGS) -o $@ $(OBJS2)

clean:
	rm -f $(PROG1) $(PROG2) $(OBJS1) $(OBJS2)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(FC) $(FFLAGS) -c $<


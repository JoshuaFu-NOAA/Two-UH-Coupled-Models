PROG1 = ../../bin/createtab
PROG2 = ../../bin/anatab
PROG3 = ../../bin/tabinfo

SRCS1 = createtab.f90
SRCS2 = anatab.f90
SRCS3 = tabinfo.f90

OBJS1 =	createtab.o
OBJS2 =	anatab.o
OBJS3 =	tabinfo.o

all: $(PROG1) $(PROG2) $(PROG3)

$(PROG1): $(OBJS1)
	$(FC) $(FFLAGS) -o $@ $(OBJS1)

$(PROG2): $(OBJS2)
	$(FC) $(FFLAGS) -o $@ $(OBJS2)

$(PROG3): $(OBJS3)
	$(FC) $(FFLAGS) -o $@ $(OBJS3)

clean:
	rm -f $(PROG1) $(PROG2) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS3) 

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

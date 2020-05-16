# Generated automatically by fu.pl
 
SRCS =	

OBJS := $(SRCS:.f90=.o)

all: $(OBJS)

clean:
	rm -f $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<


ifeq ($(ARCH), CRAY_PVP)
mo_buffer_fft.o:
	$(F90) $(F90FLAGS) -Ovector1 -c mo_buffer_fft.f90
mo_grib.o:
	$(F90) $(F90FLAGS) -Ovector1 -c mo_grib.f90
endif

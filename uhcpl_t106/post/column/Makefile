
F90=frt -Am -Cpp

# F90=/net/gardiken/export/home/NAG/bin/f95 -C=all -g -gline

all : ../../data/post_1d ../../data/chvct ../../data/chua ../../ab_ecmwf_32

../../data/post_1d : post_1d.o
	$(F90) -o ../../data/post_1d post_1d.o ask_module.o mo_system.o \
                                     unit_number_module.o mo_1d_io.o

../../data/chvct: chvct.o
	$(F90) -o ../../data/chvct chvct.o mo_system.o mo_1d_io.o \
                                   ask_module.o unit_number_module.o

../../ab_ecmwf_32: ab_ecmwf_32
	cp -p ab_ecmwf_32 ../../ab_ecmwf_32

../../data/chua: chua.o
	$(F90) -o ../../data/chua chua.o mo_system.o mo_1d_io.o \
                                   ask_module.o unit_number_module.o

mo_1d_io.o : mo_1d_io.f90 ask_module.o mo_system.o
	$(F90) -c mo_1d_io.f90

post_1d.o : post_1d.f90 mo_system.o ask_module.o mo_1d_io.o
	$(F90) -c post_1d.f90

chvct.o : chvct.f90 mo_system.o mo_1d_io.o
	$(F90) -c chvct.f90

chua.o : chua.f90 mo_system.o mo_1d_io.o
	$(F90) -c chua.f90

mo_system.o : mo_system.f90
	$(F90) -c mo_system.f90

ask_module.o : ask_module.f90 unit_number_module.o mo_system.o
	$(F90) -c ask_module.f90

unit_number_module.o : unit_number_module.f90 mo_system.o
	$(F90) -c unit_number_module.f90

clean:
	\rm -f unit_number_module.o   unit_number_module.mod \
               ask_module.o           ask_module.mod         \
               post_1d.o              post_1d.mod            \
               mo_system.o            mo_system.mod          \
               mo_1d_io.o             mo_1d_io.mod           \
               chvct.o                html/[a-z]*

index:
	-rm -f html/[a-z]*
	../../util/f2html.pl -f ../../util/fgenrc -d html .



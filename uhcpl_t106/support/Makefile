LIB =	support

SRCS =	util_reshape.c codegb5.c util_convert.c util_gwpack.c util_pbio.c util_sysinfo.c util_system.c

OBJS := $(SRCS:.c=.o)

.SUFFIXES:
.SUFFIXES: .c .o
%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: $(LIB)

$(LIB): ../lib/lib$(LIB).a 

../lib/lib$(LIB).a: $(OBJS)
	$(AR) $(ARFLAGS) ../lib/lib$(LIB).a $(OBJS) 

clean:
	rm -f $(OBJS)

ifeq ($(ARCH), SX)
codegb5.o:
	$(CC) $(CFLAGS) -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c codegb5.c
endif

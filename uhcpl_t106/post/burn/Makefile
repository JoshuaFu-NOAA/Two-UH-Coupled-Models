PROG = ../bin/after	

SRCS = after.c

OBJS = 	after.o

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) -lm $(LDFLAGS)

clean:
	rm -f $(PROG) $(OBJS)

.SUFFIXES: $(SUFFIXES) .c

%.o: %.c
	$(CC) $(CFLAGS) -c $<

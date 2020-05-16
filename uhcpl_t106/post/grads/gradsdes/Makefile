PROG =	../../bin/gradsdes 

SRCS =	gradsdes.c gamach.c gauss.c 

OBJS =	gradsdes.o gamach.o gauss.o 

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) -lm

clean:
	rm -f $(PROG) $(OBJS)

gradsdes.o: gradsdes.h

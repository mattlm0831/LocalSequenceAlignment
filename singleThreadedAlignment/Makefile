CC=gcc
CFLAGS= -g -O0 -Wall
OBJS= $(patsubst %.o,../headers/%.o, max_data.o cigar.o matrix.o) main.o

all: alignST

alignST: $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) -o $@

.c:.o
	$(CC) -c $(CFLAGS) $<

clean:
	rm -rf *.o alignMT


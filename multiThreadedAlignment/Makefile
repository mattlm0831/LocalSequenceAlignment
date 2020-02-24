CC=nvcc
CFLAGS= -g -O0
OBJS= $(patsubst %.o,../headers/%.o, max_data.o cigar.o matrix.o) unixtimer.o main.o

all: alignMT

unixtimer.o: unixtimer.c
	$(CC) $(CFLAGS) -c unixtimer.c

alignMT: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

.cu:.o
	$(CC) -c $(CFLAGS) $<

clean:
	rm -rf *.o alignMT


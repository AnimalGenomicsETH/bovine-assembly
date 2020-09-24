CC=gcc
CFLAGS=-g -Wall -O3 -march=native -Wno-unused-function

all: fasterqc

fasterqc: fasterqc.c kseq.h
		$(CC) $(CFLAGS) fasterqc.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out *~ *.a *.dSYM session*

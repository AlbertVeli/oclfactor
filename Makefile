CC = gcc
#CFLAGS = -W -Wall -O0 -g -fprofile-arcs -ftest-coverage
CFLAGS = -W -Wall -O0 -g
#LIBS = -L. -lmini-gmp
LIBS = -lgmp

all: factor

#mini-gmp.o: mini-gmp.c mini-gmp.h
#	$(CC) -c $(CFLAGS) -o $@ mini-gmp.c

#mini-gmp.a: mini-gmp.o
#	$(AR) rcs libmini-gmp.a mini-gmp.o

factor.o: factor.c
	$(CC) -c $(CFLAGS) -o $@ $<

#factor: factor.o mini-gmp.a
#	$(CC) $(CFLAGS) -o $@ factor.o $(LIBS)

factor: factor.o
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

.PHONY: clean

clean:
	rm -f factor *.o *.a *.gcov *.gcda *.gcno

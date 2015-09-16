# GNU make Makefile to compile simple C++ project

CFLAGS=-Wall -fPIC -O2

LDFLAGS= 
LDLIBS= -lc

.PHONY: all

all: hadamard.so

hadamard.o: hadamard.c
	gcc $(CFLAGS) $(LDLIBS) $(LDFLAGS) -c $<

hadamard.so: hadamard.o
	ld -shared $< -o $@ -fPIC $(LDLIBS) $(LDFLAGS)

clean:
	rm -f *.o
	rm -f hadamard.so

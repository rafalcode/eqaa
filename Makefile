CC=gcc
CFLAGS0=-g -Wall
CFLAGS=-g -Wall -std=c99
DBGCFLAGS=-g -Wall -DDBG

EXECUTABLES=eqaa0 eqaa subeqa

#Note if you're looking for the occurence coutner in here, forget it, it has its own repo in "juegooca"

# the original convoluted eqaa.c ifrom summer 2017 renamed to eqaa0.c
eqaa0: eqaa0.c
	${CC} ${CFLAGS} -o $@ $^

# the main one
eqaa: eqaa.c
	${CC} ${CFLAGS} -o $@ $^

# this is a sanity checker / utility, run without arguments to see how to use.
subeqa: subeqa.c
	${CC} ${CFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXECUTABLES}

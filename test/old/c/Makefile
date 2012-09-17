OS := $(shell uname)

CFLAGS=

ifeq ($(OS), Linux)
	CFLAGS=-fPIC -O3
endif

CC=gcc
DEPFILE=.depend
PROGS=mk_fbin rd_seq_fbin iterate_fbin
OBJS=libreria_gz.o
LIBS=-lz
LIB_NAME=liblibreria_gz

all: $(OBJS) $(PROGS) lib

$(PROGS):
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) $@.c -o $@

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

lib:
ifeq ($(OS), Linux)
	$(CC) -shared -Wl,-soname,$(LIB_NAME).so.1 -o $(LIB_NAME).so $(OBJS) $(LIBS) 
endif

ifeq ($(OS), Darwin)
	$(CC) $(LIBS) -dynamiclib -o $(LIB_NAME).dylib -dylib $(OBJS)
endif
 
clean:
	rm -f *.o $(PROGS) *.dylib *.so*

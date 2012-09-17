prefix=/usr/local
bin_prefix=$(prefix)/bin
lib_prefix=$(prefix)/lib


OS := $(shell uname)

CFLAGS=

ifeq ($(OS), Linux)
	CFLAGS=-fPIC -O3
	lib_name=$(LIB_NAME).so
endif

ifeq ($(OS), Darwin)
	lib_name=$(LIB_NAME).dylib
endif


CC=gcc
DEPFILE=.depend
PROGS=mk_fbin read_fbin iterate_fbin hash_fbin fq 
OBJS=libfbin.o
LIBS=-lz
LIB_NAME=libfbin

all: $(OBJS) $(PROGS) lib

$(PROGS):
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) $@.c -o $@

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

lib:
ifeq ($(OS), Linux)
	lib_name=$(LIB_NAME).so
	$(CC) -shared -Wl,-soname,$(LIB_NAME).so.1 -o $(LIB_NAME).so $(OBJS) $(LIBS)
endif

ifeq ($(OS), Darwin)
	lib_name=$(LIB_NAME).dylib
	$(CC) $(LIBS) -dynamiclib -o $(LIB_NAME).dylib -dylib $(OBJS)
endif
 
clean:
	rm -f *.o $(PROGS) *.dylib *.so*

install:
	
	install -m 0755 $(lib_name) $(lib_prefix)
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(bin_prefix); \
	done; \
	install -m 0755 sort_index.sh $(bin_prefix); 

uninstall:
	rm  $(lib_prefix)/$(lib_name)
	for prog in $(PROGS); do \
		rm $(bin_prefix)/$$prog; \
	done; \
	rm $(bin_prefix)/sort_index.sh
	
	

CC = g++
AR = ar
PREFIX=/usr/local

# Uncomment one of the following to switch between debug and opt mode
OPT = -DNDEBUG -Ofast
#OPT = -g -ggdb -pg

CFLAGS += --std=c++11 -fno-strict-aliasing -Wall -c -I. -I./include -I./src/ $(OPT)

LDFLAGS+= -Wall -lpthread -lssl -lcrypto -lstxxl

LIBOBJECTS = \
	./src/hashutil.o \

HEADERS = $(wildcard src/*.h)
ALIB = libcuckoofilter.a

TEST = test test_hash_accesses test_ext_throughput

all: $(TEST)

clean:
	rm -f $(TEST) */*.o

test: example/test.o $(LIBOBJECTS) 
	$(CC) example/test.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_hash_accesses: example/test_hash_accesses.o $(LIBOBJECTS) 
	$(CC) example/test_hash_accesses.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_ext_throughput: example/test_ext_throughput.o $(LIBOBJECTS) 
	$(CC) example/test_ext_throughput.o $(LIBOBJECTS) $(LDFLAGS) -o $@

%.o: %.cc ${HEADERS} Makefile
	$(CC) $(CFLAGS) $< -o $@

$(ALIB): $(LIBOBJECTS)
	$(AR) rcs $@ $(LIBOBJECTS)

.PHONY: install
install: $(ALIB)
	install -D -m 0755 $(HEADERS) -t $(DESTDIR)$(PREFIX)/include/cuckoofilter
	install -D -m 0755 $< -t $(DESTDIR)$(PREFIX)/lib

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/$(ALIB)
	rm -rf $(DESTDIR)$(PREFIX)/include/cuckoofilter

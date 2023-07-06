CC = g++
AR = ar
PREFIX=/usr/local

# Uncomment one of the following to switch between debug and opt mode
OPT = -DNDEBUG -Ofast #-pg
#OPT = -g -ggdb

CFLAGS += --std=c++11 -fno-strict-aliasing -Wall -c -I. -I./include -I./src/ $(OPT)

LDFLAGS+= -Wall -lpthread -lssl -lcrypto -lstxxl -lsplinterdb

LIBOBJECTS = \
	./src/hashutil.o \

HEADERS = $(wildcard src/*.h)
ALIB = libcuckoofilter.a

TEST = test_ext_throughput test_churn #test test_fp test_bits test_micro test_fill_varied_throughput

all: $(TEST)

clean:
	rm -f $(TEST) */*.o

test: example/test.o $(LIBOBJECTS) 
	$(CC) example/test.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_fill_varied_throughput: example/test_fill_varied_throughput.o $(LIBOBJECTS) 
	$(CC) example/test_fill_varied_throughput.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_fp: example/test_fp.o $(LIBOBJECTS) 
	$(CC) example/test_fp.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_bits: example/test_bits.o $(LIBOBJECTS)
	$(CC) example/test_bits.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_micro: example/test_micro.o $(LIBOBJECTS)
	$(CC) example/test_micro.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_ext_throughput: example/test_ext_throughput.o $(LIBOBJECTS)
	$(CC) example/test_ext_throughput.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_churn: example/test_churn.o $(LIBOBJECTS)
	$(CC) example/test_churn.o $(LIBOBJECTS) $(LDFLAGS) -o $@

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

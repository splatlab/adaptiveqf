CC = g++
AR = ar
PREFIX=/usr/local

# Uncomment one of the following to switch between debug and opt mode
ifdef D
	OPT = -g -ggdb
	#SPLINTERPATH=../../external/splinterdb/build/debug/lib
	SPLINTERPATH=../../external/splinterdb/btree
else
	OPT = -DNDEBUG -Ofast
	#SPLINTERPATH=../../external/splinterdb/build/release/lib
	SPLINTERPATH=../../external/splinterdb/btree
endif

CFLAGS += --std=c++11 -fno-strict-aliasing -Wall -c -I. -I./include -I./src/ -I../../external/splinterdb/include -DSPLINTERDB_PLATFORM_DIR=platform_linux $(OPT)

LDFLAGS+= -Wall -lpthread -lssl -lcrypto -L$(SPLINTERPATH) -lsplinterdb -Wl,-rpath=$(SPLINTERPATH)

LIBOBJECTS = \
	./src/hashutil.o \

HEADERS = $(wildcard src/*.h)
ALIB = libcuckoofilter.a

TEST = test_splinter_ops test_splinter_throughput test_splinter_zipfian_histogram test_splinter_adversarial test_micro_throughput #test_acf_ops #test_churn test_ext_throughput #test test_fp test_bits test_micro test_fill_varied_throughput

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

test_splinter_ops: example/test_splinter_ops.o $(LIBOBJECTS)
	$(CC) example/test_splinter_ops.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_splinter_zipfian_histogram: example/test_splinter_zipfian_histogram.o $(LIBOBJECTS)
	$(CC) example/test_splinter_zipfian_histogram.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_splinter_throughput: example/test_splinter_throughput.o $(LIBOBJECTS)
	$(CC) example/test_splinter_throughput.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_splinter_adversarial: example/test_splinter_adversarial.o $(LIBOBJECTS)
	$(CC) example/test_splinter_adversarial.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_acf_ops: example/test_acf_ops.o $(LIBOBJECTS)
	$(CC) example/test_acf_ops.o $(LIBOBJECTS) $(LDFLAGS) -o $@

test_micro_throughput: example/test_micro_throughput.o $(LIBOBJECTS)
	$(CC) example/test_micro_throughput.o $(LIBOBJECTS) $(LDFLAGS) -o $@

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

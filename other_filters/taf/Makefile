CTARGETS=test test_threadsafe test_pc bm test_throughput test_fill_varied_throughput test_near_full test_deletions test_merge test_hash_accesses test_bulk test_whitelist
CXXTARGETS=taf
# test_progress

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-Ofast -DNDEBUG
endif

ifndef NH
	ARCH=
else
	ARCH=-msse4.2 -D__SSE4_2_
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

LOC_INCLUDE=src
LOC_SRC=src
LOC_TEST=src
OBJDIR=src

CC = gcc -std=gnu11
CXX = g++ -std=c++11
LD= gcc -std=gnu11

CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -lpthread -lssl -lcrypto -m64 -I. -Iinclude

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lm

#
# declaration of dependencies
#

all: $(CXXTARGETS)

taf:									$(OBJDIR)/arcd.o $(OBJDIR)/murmur3.o $(OBJDIR)/bit_util.o $(OBJDIR)/set.o

# dependencies between .o files and .h files

# dependencies between .o files and .cc (or .c) files

$(OBJDIR)/arcd.o:							$(LOC_SRC)/arcd.c $(LOC_INCLUDE)/arcd.h
$(OBJDIR)/murmur3.o:					$(LOC_SRC)/murmur3.c $(LOC_INCLUDE)/murmur3.h
$(OBJDIR)/bit_util.o:					$(LOC_SRC)/bit_util.c $(LOC_INCLUDE)/bit_util.h
$(OBJDIR)/set.o:						$(LOC_SRC)/set.c $(LOC_INCLUDE)/set.h

#
# generic build rules
#

$(CTARGETS):
	$(LD) $^ -o $@ $(LDFLAGS)

$(CXXTARGETS):
	$(CXX) $^ -o $@ $(CXXFLAGS)

$(OBJDIR)/%.o: $(LOC_SRC)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

$(OBJDIR)/%.o: $(LOC_SRC)/%.c | $(OBJDIR)
	$(CC) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

$(OBJDIR)/%.o: $(LOC_TEST)/%.c | $(OBJDIR)
	$(CC) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(CTARGETS) $(CXXTARGETS) core

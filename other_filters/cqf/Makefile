CTARGETS=test_splinter_ops test_splinter_throughput test_splinter_zipfian_histogram test_splinter_adversarial test_micro_throughput test_micro_super test_parallel test_threadsafe #test_pc bm test_hash_accesses #test_fill_varied_throughput #test
CXXTARGETS=test_ext_throughput test_adversarial

ifdef D
	DEBUG=-g
	OPT=
	#SPLINTERPATH=../../external/splinterdb/build/debug/lib
	SPLINTERPATH=../../external/splinterdb/btree
else
	DEBUG=
	OPT=-O3 -g
	#SPLINTERPATH=../../external/splinterdb/build/release/lib
	SPLINTERPATH=../../external/splinterdb/btree
endif

ifdef NH
	ARCH=
else
	ARCH=-msse4.2 -D__SSE4_2_
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

LOC_INCLUDE=include
LOC_SRC=src
LOC_TEST=test
OBJDIR=obj

CC = gcc -std=gnu11
CXX = g++ -std=c++11
LD= gcc -std=gnu11

#CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -Iinclude
CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -DSPLINTERDB_PLATFORM_DIR=platform_linux -DSKIP_BOOL_DEF -D_GNU_SOURCE -m64 -I. -Iinclude -I../../external/splinterdb/include

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lm -L$(SPLINTERPATH) -lsplinterdb -Wl,-rpath=$(SPLINTERPATH)

#
# declaration of dependencies
#

all: $(CTARGETS) #$(CXXTARGETS)

# dependencies between programs and .o files

test:								$(OBJDIR)/test.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_micro_throughput:                                                    $(OBJDIR)/test_micro_throughput.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o 

test_micro_super:                                                    $(OBJDIR)/test_micro_super.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_parallel:                                                    $(OBJDIR)/test_parallel.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_fill_varied_throughput:								$(OBJDIR)/test_fill_varied_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o


test_hash_accesses:								$(OBJDIR)/test_hash_accesses.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_threadsafe:		$(OBJDIR)/test_threadsafe.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/gqf_file.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_pc:						$(OBJDIR)/test_partitioned_counter.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/gqf_file.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

bm:									$(OBJDIR)/bm.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_ext_throughput:                                                    $(OBJDIR)/test_ext_throughput.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_adversarial:                                                    $(OBJDIR)/test_adversarial.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_splinter_ops:                                                    $(OBJDIR)/test_splinter_ops.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_splinter_throughput:                                                    $(OBJDIR)/test_splinter_throughput.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_splinter_zipfian_histogram:                                                    $(OBJDIR)/test_splinter_zipfian_histogram.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

test_splinter_adversarial:                                                    $(OBJDIR)/test_splinter_adversarial.o $(OBJDIR)/gqf.o \
                                                                                $(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
                                                                                $(OBJDIR)/partitioned_counter.o

# dependencies between .o files and .h files

$(OBJDIR)/test.o: 						$(LOC_INCLUDE)/gqf.h $(LOC_INCLUDE)/gqf_file.h \
															$(LOC_INCLUDE)/hashutil.h \
															$(LOC_INCLUDE)/partitioned_counter.h

$(OBJDIR)/test_threadsafe.o: 	$(LOC_INCLUDE)/gqf.h $(LOC_INCLUDE)/gqf_file.h \
															$(LOC_INCLUDE)/hashutil.h \
															$(LOC_INCLUDE)/partitioned_counter.h

$(OBJDIR)/bm.o:								$(LOC_INCLUDE)/gqf_wrapper.h \
															$(LOC_INCLUDE)/partitioned_counter.h


# dependencies between .o files and .cc (or .c) files

$(OBJDIR)/gqf.o:							$(LOC_SRC)/gqf.c $(LOC_INCLUDE)/gqf.h
$(OBJDIR)/gqf_file.o:					$(LOC_SRC)/gqf_file.c $(LOC_INCLUDE)/gqf_file.h
$(OBJDIR)/hashutil.o:					$(LOC_SRC)/hashutil.c $(LOC_INCLUDE)/hashutil.h
$(OBJDIR)/partitioned_counter.o:	$(LOC_INCLUDE)/partitioned_counter.h

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
	rm -rf $(OBJDIR) $(TARGETS) core


CTARGETS=test test_threadsafe test_pc bm test_throughput test_fill_varied_throughput test_near_full test_deletions test_merge test_hash_accesses test_bulk test_whitelist test_resize
CXXTARGETS=test_ext_throughput test_zipf_throughput test_ext_churn test_adversarial taf
# test_progress

ifndef D
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

LOC_INCLUDE=include
LOC_SRC=src
LOC_TEST=test
OBJDIR=obj

CC = gcc -std=gnu11
CXX = g++ -std=c++11
LD= gcc -std=gnu11

CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -lpthread -lssl -lcrypto -lstxxl -m64 -I. -Iinclude

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lm

#
# declaration of dependencies
#

all: $(CTARGETS) $(CXXTARGETS)

# dependencies between programs and .o files

test:									$(OBJDIR)/test.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_progress:								$(OBJDIR)/test_progress.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_throughput:							$(OBJDIR)/test_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_fill_varied_throughput:						$(OBJDIR)/test_fill_varied_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_near_full:								$(OBJDIR)/test_near_full.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_deletions:								$(OBJDIR)/test_deletions.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_merge:								$(OBJDIR)/test_merge.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_bulk:								$(OBJDIR)/test_bulk.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_whitelist:								$(OBJDIR)/test_whitelist.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_hash_accesses:							$(OBJDIR)/test_hash_accesses.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_resize:							$(OBJDIR)/test_resize.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_threadsafe:							$(OBJDIR)/test_threadsafe.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/gqf_file.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_pc:								$(OBJDIR)/test_partitioned_counter.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/gqf_file.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

bm:									$(OBJDIR)/bm.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_ext_throughput:							$(OBJDIR)/test_ext_throughput.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_zipf_throughput:							$(OBJDIR)/test_zipf_throughput.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_ext_churn:							$(OBJDIR)/test_ext_churn.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

test_adversarial:							$(OBJDIR)/test_adversarial.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o \
										$(OBJDIR)/partitioned_counter.o

taf:									$(OBJDIR)/taf.o $(OBJDIR)/hashutil.o

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

$(OBJDIR)/gqf.o:					$(LOC_SRC)/gqf.c $(LOC_INCLUDE)/gqf.h
$(OBJDIR)/gqf_file.o:					$(LOC_SRC)/gqf_file.c $(LOC_INCLUDE)/gqf_file.h
$(OBJDIR)/hashutil.o:					$(LOC_SRC)/hashutil.c $(LOC_INCLUDE)/hashutil.h
$(OBJDIR)/partitioned_counter.o:			$(LOC_INCLUDE)/partitioned_counter.h
$(OBJDIR)/ll_table.o:					$(LOC_SRC)/ll_table.c $(LOC_INCLUDE)/ll_table.h

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

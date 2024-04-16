CTARGETS=test_micro test_throughput test_split_throughput test_adversarial test_splinter_lltable_throughput test_splinter_parallel test_parallel# test_fill_varied_throughput test_near_full test_deletions test_merge test_hash_accesses test_bulk test_whitelist test_resize test_micro_throughput test_micro_write test_micro_read test_lltable_throughput
CXXTARGETS=test_ext_throughput test_ext_inc_throughput test_zipf_throughput test_ext_churn taf
SPLTARGETS=test_splinter_ops test_splinter_inserts test_splinter_inserts_2 test_splinter_throughput test_splinter_zipfian_histogram test_splinter_adversarial
# test_progress

ifndef D
	DEBUG=
	OPT=-O3 -DNDEBUG
	SPLINTERPATH=external/splinterdb/build/release/lib
	#SPLINTERPATH=external/splinterdb/btree
else
	DEBUG=-g
	OPT=-O0
	SPLINTERPATH=external/splinterdb/build/debug/lib
	#SPLINTERPATH=external/splinterdb/btree
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

CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -Iinclude -Iexternal/splinterdb/include -DSPLINTERDB_PLATFORM_DIR=platform_linux -DSKIP_BOOL_DEF -D_GNU_SOURCE

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lm -L$(SPLINTERPATH) -lsplinterdb -Wl,-rpath=$(SPLINTERPATH)
#LDFLAGS += -L/usr/lib/ -lstxxl

#
# declaration of dependencies
#

all: $(CTARGETS) #$(SPLTARGETS) #$(CXXTARGETS)

# dependencies between programs and .o files

test_unit:								$(OBJDIR)/test_unit.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_throughput:						$(OBJDIR)/test_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/splinter_util.o $(OBJDIR)/test_driver.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/ll_table.o $(OBJDIR)/rand_util.o

test_split_throughput:						$(OBJDIR)/test_split_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/splinter_util.o $(OBJDIR)/test_driver.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/ll_table.o $(OBJDIR)/rand_util.o

test_adversarial:						$(OBJDIR)/test_adversarial.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/splinter_util.o $(OBJDIR)/test_driver.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/ll_table.o $(OBJDIR)/rand_util.o

test_micro:								$(OBJDIR)/test_micro.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/splinter_util.o $(OBJDIR)/test_driver.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/ll_table.o $(OBJDIR)/rand_util.o

test_splinter_parallel:						$(OBJDIR)/test_splinter_parallel.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/splinter_util.o $(OBJDIR)/test_driver.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/ll_table.o $(OBJDIR)/rand_util.o

test_parallel:						$(OBJDIR)/test_parallel.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/splinter_util.o $(OBJDIR)/test_driver.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/ll_table.o $(OBJDIR)/rand_util.o

test_progress:								$(OBJDIR)/test_progress.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_throughput_old:							$(OBJDIR)/test_throughput_old.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_micro_throughput:							$(OBJDIR)/test_micro_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_micro_write:							$(OBJDIR)/test_micro_write.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_micro_read:							$(OBJDIR)/test_micro_read.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_lltable_throughput:							$(OBJDIR)/test_lltable_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_fill_varied_throughput:						$(OBJDIR)/test_fill_varied_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_near_full:								$(OBJDIR)/test_near_full.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_deletions:								$(OBJDIR)/test_deletions.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_merge:								$(OBJDIR)/test_merge.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_bulk:								$(OBJDIR)/test_bulk.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_whitelist:								$(OBJDIR)/test_whitelist.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_hash_accesses:							$(OBJDIR)/test_hash_accesses.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_resize:							$(OBJDIR)/test_resize.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_threadsafe:							$(OBJDIR)/test_threadsafe.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/gqf_file.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_pc:								$(OBJDIR)/test_partitioned_counter.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/gqf_file.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

bm:										$(OBJDIR)/bm.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_ext_throughput:							$(OBJDIR)/test_ext_throughput.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_ext_inc_throughput:							$(OBJDIR)/test_ext_inc_throughput.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_zipf_throughput:							$(OBJDIR)/test_zipf_throughput.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_ext_churn:							$(OBJDIR)/test_ext_churn.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

taf:									$(OBJDIR)/taf.o $(OBJDIR)/hashutil.o

test_splinter_ops:							$(OBJDIR)/test_splinter_ops.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_splinter_inserts:							$(OBJDIR)/test_splinter_inserts.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_splinter_inserts_2:							$(OBJDIR)/test_splinter_inserts_2.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_splinter_throughput:							$(OBJDIR)/test_splinter_throughput.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_splinter_zipfian_histogram:							$(OBJDIR)/test_splinter_zipfian_histogram.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_splinter_adversarial:							$(OBJDIR)/test_splinter_adversarial.o $(OBJDIR)/gqf.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o \
										$(OBJDIR)/partitioned_counter.o

test_splinter_lltable_throughput:		$(OBJDIR)/test_splinter_lltable_throughput.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
										$(OBJDIR)/zipf.o $(OBJDIR)/hashutil.o $(OBJDIR)/ll_table.o $(OBJDIR)/splinter_util.o \
										$(OBJDIR)/partitioned_counter.o $(OBJDIR)/rand_util.o

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

$(OBJDIR)/gqf.o:						$(LOC_SRC)/gqf.c $(LOC_INCLUDE)/gqf.h
$(OBJDIR)/gqf_file.o:					$(LOC_SRC)/gqf_file.c $(LOC_INCLUDE)/gqf_file.h
$(OBJDIR)/hashutil.o:					$(LOC_SRC)/hashutil.c $(LOC_INCLUDE)/hashutil.h
$(OBJDIR)/partitioned_counter.o:		$(LOC_INCLUDE)/partitioned_counter.h
$(OBJDIR)/ll_table.o:					$(LOC_SRC)/ll_table.c $(LOC_INCLUDE)/ll_table.h
$(OBJDIR)/splinter_util.o:				$(LOC_SRC)/splinter_util.c $(LOC_INCLUDE)/splinter_util.h# $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o
$(OBJDIR)/test_driver.o:				$(LOC_SRC)/test_driver.c $(LOC_INCLUDE)/test_driver.h# $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o $(OBJDIR)/splinter_util.o

#
# generic build rules
#

$(CTARGETS):
	$(LD) $^ -o $@ $(LDFLAGS)

$(SPLTARGETS):
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
	rm -rf $(OBJDIR) $(CTARGETS) $(CXXTARGETS) $(SPLTARGETS) core

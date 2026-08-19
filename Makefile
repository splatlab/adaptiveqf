# [TODO-P5] This Makefile is being cleaned up. Current: ~57 targets → target: 4-5.
# [TODO-P5] Migrate to sponge/build/ output.
# [TODO-P2] Add bench_ filters and workload_gen targets when ready.
# See adaptiveqf/TODO.md Phase 5 for the full plan.

ifndef D
	DEBUG=
	OPT=-O3 -DNDEBUG
	SPLINTERPATH=external/splinterdb/build/release/lib
else
	DEBUG=-g
	OPT=-O0
	SPLINTERPATH=external/splinterdb/build/debug/lib
endif

ifdef NH
	ARCH=
else
	ARCH=-msse4.2 -D__SSE4_2_
endif

ifdef P
	PROFILE=-pg -no-pie
endif

LOC_INCLUDE=include
LOC_SRC=src
OBJDIR=obj

CC = gcc -std=gnu11
CXX = g++ -std=c++17
LD = gcc -std=gnu11

CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 \
	-I. \
	-Iinclude \
	-Ireverse_maps/ll_table \
	-Ireverse_maps/splinterdb \
	-Ireverse_maps/mock \
	-Iexternal \
	-Iexternal/splinterdb/include \
	-DSPLINTERDB_PLATFORM_DIR=platform_linux \
	-DSKIP_BOOL_DEF \
	-D_GNU_SOURCE

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) \
	-lpthread -lssl -lcrypto -lm \
	-L$(SPLINTERPATH) -lsplinterdb -Wl,-rpath=$(SPLINTERPATH)

# Core source files
CORE_OBJS = \
	$(OBJDIR)/gqf.o \
	$(OBJDIR)/hashutil.o \
	$(OBJDIR)/partitioned_counter.o \
	$(OBJDIR)/aqf.o

# Reverse map source files
LLTABLE_OBJS = $(OBJDIR)/ll_table.o

SPLINTER_OBJS = $(OBJDIR)/splinter_util.o

all: libaqf.a liblltable.a

libaqf.a: $(CORE_OBJS)
	ar rcs $@ $^

liblltable.a: $(LLTABLE_OBJS)
	ar rcs $@ $^

ifneq ($(wildcard $(SPLINTERPATH)/libsplinterdb.a),)
all: libsplinterutil.a
libsplinterutil.a: $(SPLINTER_OBJS)
	ar rcs $@ $^
endif

# Core source compilation
$(OBJDIR)/gqf.o: src/gqf.c include/gqf.h include/gqf_int.h reverse_maps/ll_table/ll_table.h | $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/hashutil.o: src/hashutil.c include/hashutil.h | $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/partitioned_counter.o: include/partitioned_counter.h | $(OBJDIR)
	$(CC) $(CXXFLAGS) -c src/partitioned_counter.c -o $@

$(OBJDIR)/aqf.o: src/aqf.c include/aqf.h include/gqf.h include/gqf_int.h include/hashutil.h | $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@

# Reverse map compilation
$(OBJDIR)/ll_table.o: reverse_maps/ll_table/ll_table.c reverse_maps/ll_table/ll_table.h | $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@

ifneq ($(wildcard $(SPLINTERPATH)/libsplinterdb.a),)
$(OBJDIR)/splinter_util.o: reverse_maps/splinterdb/splinter_util.c reverse_maps/splinterdb/splinter_util.h | $(OBJDIR)
	$(CC) $(CXXFLAGS) -c $< -o $@
endif

$(OBJDIR):
	@mkdir -p $(OBJDIR)

.PHONY: clean
clean:
	rm -rf $(OBJDIR) libaqf.a liblltable.a libsplinterutil.a

# LIBRARY SETTINGS - SET AS NECESSARY
# 
# For recent versions of BamTools (>= 2.5):
#
BOOST_PATH ?= /your/boost/dir
BAMTOOLS_PATH ?= /your/bamtools/dir
BOOST_INCLUDE = $(BOOST_PATH)/include
BOOST_LIB = $(BOOST_PATH)/lib
BAMTOOLS_INCLUDE = $(BAMTOOLS_PATH)/include/bamtools
BAMTOOLS_LIB = $(BAMTOOLS_PATH)/lib64

INCS = -I$(BOOST_INCLUDE) -I$(BAMTOOLS_INCLUDE) -I.
LIBS = -L$(BOOST_LIB) -L$(BAMTOOLS_LIB) -lboost_random -lboost_filesystem -lboost_system  -lbamtools -lz -lboost_serialization -lbam -pthread

# use the following for older versions of BamTools (e.g. commit https://github.com/pezmaster31/bamtools/commit/2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e)
# 
# BOOST_PATH ?= /data/projects/phillippy/software/boost_1_60_0/
# BAMTOOLS_PATH ?= /data/projects/phillippy/software/bamtools
# BOOST_INCLUDE = $(BOOST_PATH)/include
# BOOST_LIB = $(BOOST_PATH)/lib
# BAMTOOLS_INCLUDE = $(BAMTOOLS_PATH)/include
# BAMTOOLS_LIB = $(BAMTOOLS_PATH)/lib

# INCS = -I$(BOOST_INCLUDE) -I$(BAMTOOLS_INCLUDE) -I$(BAMTOOLS_SRC)
# LIBS = -L$(BOOST_LIB) -L$(BAMTOOLS_LIB) -lboost_random -lboost_filesystem -lboost_system  -lbamtools -lbamtools-utils -lz -lboost_serialization


MKDIR_P = mkdir -p

.PHONY: directories
	
# END LIBRARY SETTINGS

#
# object and binary dirs  
#

DIR_OBJ = ./obj1
DIR_OBJ2 = ./obj2
DIR_OBJ3 = ./obj3
DIR_OBJ4 = ./obj4
DIR_OBJ5 = ./obj5
DIR_BIN = ./bin
OUT_DIR = ./bin ./obj1 ./obj2 ./obj3 ./obj4 ./obj5
directories: ${OUT_DIR}

#-- for now these MUST point to the included "samtools-0.x.x" and "gclib" sub-directories
BAM  := ./samtools-0.1.18
GDIR := ./gclib
#--

INCDIRS := -I. -I${GDIR} -I${BAM}

CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++0x -fno-strict-aliasing -fno-exceptions -fno-rtti
#for gcc 8+ add: -Wno-class-memaccess
GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

LDFLAGS += -L${BAM}

# LIBS    := -lbam -lz

ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
 NOTHREADS=1
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
ifdef WINDOWS
 RM = del /Q
else
 RM = rm -f
endif

# File endings
ifdef WINDOWS
 EXE = .exe
else
 EXE =
endif

# Non-windows systems need pthread
#ifndef WINDOWS
# ifndef NOTHREADS
#   LIBS := -pthread ${LIBS}
#   BASEFLAGS += -pthread
# endif
#endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
  else
     #just plain debug build
     DEBUG_BUILD=1
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     #CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-ggdb -g3 -O0 -fvar-tracking-assignments -fno-omit-frame-pointer)
     ifneq (, $(findstring darwin, $(DMACH)))
        CXXFLAGS += -gdwarf-3
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  endif
endif

ifdef RELEASE_BUILD
 ifneq (,$(findstring static, $(MAKECMDGOALS)))
    STATIC_CLIB=1
 endif
endif

ifdef STATIC_CLIB
 LDFLAGS += -static-libgcc -static-libstdc++
endif

ifdef DEBUG_BUILD
  #$(warning Building DEBUG version of gingko_abundance.. )
  DBG_WARN=@echo
#  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make clean release" for a faster, optimized version of the program.'
endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GBam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CXXFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

OBJS += rlink.o tablemaker.o tmerge.o

all release static debug: $(DIR_BIN)/gingko_abundance${EXE}
memcheck memdebug tsan tcheck thrcheck: $(DIR_BIN)/gingko_abundance${EXE}
memuse memusage memtrace: $(DIR_BIN)/gingko_abundance${EXE}
nothreads: $(DIR_BIN)/gingko_abundance${EXE}

${GDIR}/GBam.o : $(GDIR)/GBam.h
gingko_abundance.o : $(GDIR)/GBitVec.h $(GDIR)/GHash.hh $(GDIR)/GBam.h
rlink.o : rlink.h tablemaker.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h
tmerge.o : rlink.h tmerge.h
tablemaker.o : tablemaker.h rlink.h
${BAM}/libbam.a: 
	cd ${BAM} && make lib
$(DIR_BIN)/gingko_abundance: ${BAM}/libbam.a $(OBJS) gingko_abundance.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
test demo tests: gingko_abundance${EXE}
	@./run_tests.sh



COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS = 
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS) ${LDFLAGS}
VPATH = .:metapathsearch
        
OBJS1 = \
        $(DIR_OBJ)/assembler.o \
        $(DIR_OBJ)/simplify_graph.o \
        $(DIR_OBJ)/junction_graph.o \
        $(DIR_OBJ)/QuadProg++.o \

OBJS2 = \
        $(DIR_OBJ2)/meta_simplify_graph.o \
        $(DIR_OBJ2)/meta_junction_graph.o \
        $(DIR_OBJ2)/meta_QuadProg++.o \
# list executable file names
#
EXECS1 = gingko_graph
EXECS2 = gingko_path_search
EXECS3 = gingko_cover
EXECS4 = gingko_merge
EXECS5 = gingko_individual
#
# compile and link
#
default:
	@echo
	@echo " to build:"
	@echo "    make all"
	@echo
	@echo " to clean:"
	@echo "    make clean"
	@echo "    make realclean"
	@echo

all: directories $(EXECS1) 

$(EXECS1): $(OBJS1)
	$(foreach EX, $(EXECS1), $(COMPILE) $(EX).cpp  -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, $(EXECS1), $(COMPILE) $(OBJS1) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)


$(DIR_OBJ)/%.o: %.cpp %.h
	$(COMPILE) $< -c -o $@

all: directories $(EXECS2)
$(EXECS2): $(OBJS2)
	$(foreach EX, $(EXECS2), $(COMPILE) $(EX).cpp  -c -o $(DIR_OBJ2)/$(EX).o;)
	$(foreach EX, $(EXECS2), $(COMPILE) $(OBJS2) $(DIR_OBJ2)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)

$(DIR_OBJ2)/%.o: %.cpp %.h
	$(COMPILE) $< -c -o $@

all: directories $(EXECS3)

$(EXECS3):
	$(COMPILE) $(EXECS3).cpp  -c -o $(DIR_OBJ3)/$(EXECS3).o
	$(COMPILE) $(DIR_OBJ3)/$(EXECS3).o -o $(DIR_BIN)/$(EXECS3) $(LIBS);

all: directories $(EXECS4)

$(EXECS4):
	$(COMPILE) $(EXECS4).cpp  -c -o $(DIR_OBJ4)/$(EXECS4).o
	$(COMPILE) $(DIR_OBJ4)/$(EXECS4).o -o $(DIR_BIN)/$(EXECS4) $(LIBS);

all: directories $(EXECS5)
$(EXECS5):
	$(COMPILE) $(EXECS5).cpp  -c -o $(DIR_OBJ5)/$(EXECS5).o
	$(COMPILE) $(DIR_OBJ5)/$(EXECS5).o -o $(DIR_BIN)/$(EXECS5) $(LIBS);
#
# odds and ends
#
.PHONY : clean cleanall cleanAll allclean

# target for removing all object files

#	echo $(PATH)
clean:
	${RM} $(DIR_BIN)/gingko_abundance${EXE} gingko_abundance.o*  $(OBJS)
	${RM} core.*
	/bin/rm $(DIR_OBJ)/* $(DIR_OBJ2)/* $(DIR_OBJ3)/* $(DIR_OBJ4)/* $(DIR_OBJ5)/*

allclean cleanAll cleanall:
	cd ${BAM} && make clean
	${RM} $(DIR_BIN)/gingko_abundance${EXE} gingko_abundance.o* $(OBJS)
	${RM} core.*
	/bin/rm $(DIR_BIN)/*

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

#OBJ_DIR = ./

HEADERS = ../src/Variant.h \
		  ../src/split.h \
		  ../src/join.h \
		  pdflib.hpp \
		  var.hpp \
		  cdflib.hpp \
		  rnglib.hpp \

SOURCES = ../src/Variant.cpp \
		  ../src/split.cpp \
		  rnglib.cpp \
		  var.cpp \
		  pdflib.cpp \
		  cdflib.cpp \

BIN_SOURCES = dumpContigsFromHeader.cpp \
			  iHS.cpp \
			  bFst.cpp \
			  hapLrt.cpp \
			  popStats.cpp \
			  wcFst.cpp \
			  sequenceDiversity.cpp \
			  pFst.cpp \
			  xpEHH.cpp \
			  smoother.cpp \
			  LD.cpp \
			  plotHaps.cpp \
			  abba-baba.cpp \
			  permuteGPAT++.cpp \

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
#CXXFLAGS = -O2
#CXXFLAGS = -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual

INCLUDES = -I. -I.. -I../src -I../tabixpp/htslib/
LDINCLUDES = -L. -L.. -L../tabixpp/ -L../tabixpp/htslib/
LDFLAGS = -lvcflib -ltabix -lhts -lpthread -lz -lm

OBJECTS= $(SOURCES:.cpp=.o)
BIN_OBJECTS= $(BIN_SOURCES:.cpp=.o)

BINS = $(addprefix ../bin/,$(notdir $(BIN_SOURCES:.cpp=)))
SHORTBINS = $(notdir $(BIN_SOURCES:.cpp=))

all: $(OBJECTS) $(BIN_OBJECTS) $(BINS)

openmp:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fopenmp -D HAS_OPENMP"

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS) $(BIN_OBJECTS): $(SOURCES) $(BIN_SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(@:.o=.cpp) $(INCLUDES) $(CXXFLAGS)

$(BINS): $(BIN_OBJECTS) ../libvcflib.a $(OBJECTS) $(SMITHWATERMAN) $(FASTAHACK) $(DISORDER) $(LEFTALIGN) $(INDELALLELE) $(SSW) $(FILEVERCMP) $(GPAT_BINS)
	$(CXX) $(notdir $@).o $(OBJECTS) -o $@ $(INCLUDES) $(LDINCLUDES) $(LDFLAGS) $(CXXFLAGS)

#test: $(BINS)
#	@prove -Itests/lib -w tests/*.t

clean:
	rm -f $(BINS) $(BIN_OBJECTS) $(OBJECTS)

.PHONY: clean all test

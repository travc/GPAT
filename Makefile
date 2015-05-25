VCFLIB_PATH = ..

EXTRA_HEADERS = ../src/join.h

SOURCES = $(VCFLIB_PATH)/src/Variant.cpp \
		  $(VCFLIB_PATH)/src/split.cpp \
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

INCLUDES = -I. -I$(VCFLIB_PATH) -I$(VCFLIB_PATH)/src -I$(VCFLIB_PATH)/tabixpp/htslib/
LDINCLUDES = -L. -L$(VCFLIB_PATH) -L$(VCFLIB_PATH)/tabixpp/ -L$(VCFLIB_PATH)/tabixpp/htslib/
LDFLAGS = -lvcflib -ltabix -lhts -lpthread -lz -lm

OBJECTS= $(SOURCES:.cpp=.o)
HEADERS= $(SOURCES:.cpp=.h) $(EXTRA_HEADERS)

BINS = $(addprefix $(VCFLIB_PATH)/bin/,$(notdir $(BIN_SOURCES:.cpp=)))
SHORTBINS = $(notdir $(BIN_SOURCES:.cpp=))

all: $(OBJECTS) $(BINS)

openmp:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fopenmp -D HAS_OPENMP"

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS): $(SOURCES) $(BIN_SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(@:.o=.cpp) $(INCLUDES) $(CXXFLAGS)

$(BINS): $(VCFLIB_PATH)/libvcflib.a $(OBJECTS)
	$(CXX) $(notdir $@).cpp $(OBJECTS) -o $@ $(INCLUDES) $(LDINCLUDES) $(LDFLAGS) $(CXXFLAGS)

#test: $(BINS)
#	@prove -Itests/lib -w tests/*.t

clean:
	rm -f $(BINS) $(OBJECTS)

.PHONY: clean all test

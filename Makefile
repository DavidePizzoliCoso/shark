CFLAGS	= -DNDEBUG -march=native -Wno-char-subscripts -Wall -O3 -std=c++14 -I. -g3
CXXFLAGS= ${CFLAGS}
LIBS = -L./lib -lz -ltbb

.PHONY: all clean

all: shark

shark: main.o
	@echo "* Linking shark"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main.o: common.hpp argument_parser.hpp simpleBF.hpp bloomtree.hpp BloomfilterFiller.hpp KmerBuilder.hpp FastaSplitter.hpp FastqSplitter.hpp ReadAnalyzer.hpp ReadOutput.hpp kmer_utils.hpp

clean:
	rm -rf *.o

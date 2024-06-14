CXX = g++
CFLAGS = -pthread -std=c++17 -Wall

SEQUENCE_ALIGNMENT_FILES = alignment.hpp alignment.cpp alignment_parallel.cpp alignment_sequential.cpp 

all: main.o $(SEQUENCE_ALIGNMENT_FILES)
	$(CXX) $(CFLAGS) -o main main.o

main.o: main.cpp $(SEQUENCE_ALIGNMENT_FILES)
	$(CXX) -c $(CFLAGS) -o main.o main.cpp

BENCHMARKER_OBJECTS = benchmarker.o

benchmarker: $(BENCHMARKER_OBJECTS) $(SEQUENCE_ALIGNMENT_FILES)
	$(CXX) $(CFLAGS) -o benchmarker $(BENCHMARKER_OBJECTS) -latomic

benchmarker.o: benchmarker.cpp $(SEQUENCE_ALIGNMENT_FILES)
	$(CXX) -c $(CFLAGS) -o benchmarker.o benchmarker.cpp -I.

clean:
	rm -f *.o
	rm -f grader
	rm main test benchmarker

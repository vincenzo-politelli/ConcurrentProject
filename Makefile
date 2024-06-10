CXX = g++
CFLAGS = -pthread -std=c++17 -Wall

SOURCES = main.cpp sequence_alignment_sequential.cpp  test.cpp sequence_alignment_parallel.cpp  #sequence_alignment.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# The default target
all: main

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o main

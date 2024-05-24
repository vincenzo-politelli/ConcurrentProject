CXX = g++
CFLAGS = -pthread -std=c++17 -Wall

SOURCES = main.cpp sequence_alignment_sequential.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# The default target
all: main

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o main

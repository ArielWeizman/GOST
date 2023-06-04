#makefile

CXX=g++
CXXFLAGS=-std=c++11 -O2

all: main

main: main.cpp
	$(CXX) $(CXXFLAGS) -o main.out main.cpp

run: main
	./main.out

clean:
	rm main.out

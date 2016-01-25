CC = g++
CC_FLAGS = -O2 -std=c++11 -Wall

all: basic-test advanced-test

basic-test: basic-test.cc
	$(CC) basic-test.cc -o basic-test $(CC_FLAGS)
advanced-test: advanced-test.cc
	$(CC) advanced-test.cc -o advanced-test $(CC_FLAGS)

clean:
	-rm basic-test
	-rm advanced-test
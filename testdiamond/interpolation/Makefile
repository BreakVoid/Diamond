CC = g++
CC_FLAGS = -O2 -std=c++11 -Wall

all: natural-spline-shooting-method natural-spline-b-spline

natural-spline-shooting-method: natural-spline-shooting-method.cc
	$(CC) natural-spline-shooting-method.cc -o natural-spline-shooting-method $(CC_FLAGS)
natural-spline-b-spline: natural-spline-b-spline.cc
	$(CC) natural-spline-b-spline.cc -o natural-spline-b-spline $(CC_FLAGS)

clean:
	-rm natural-spline-b-spline
	-rm natural-spline-shooting-method
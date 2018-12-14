output: main.o methods.o
	c++ -O3 -std=c++11 main.o methods.o -larmadillo -o 2DHeatEQ.x
 main.o: main.cpp
	c++ -O3 -c main.cpp
methods.o: methods.cpp
	c++ -O3 -c methods.cpp
 clean:
	rm *.o 2DHeatEQ.x

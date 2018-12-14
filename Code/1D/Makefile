output: main.o methods.o
	c++ -O3 -std=c++11 main.o methods.o -larmadillo -o HeatEQ.x
 main.o: main.cpp
	g++ -O3 -c main.cpp
methods.o: methods.cpp
	g++ -O3 -c methods.cpp
 clean:
	rm *.o HeatEQ.x

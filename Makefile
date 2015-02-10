all:
	g++ -o test.o test.cpp
	./test.o
doc:
	doxygen
clean:
	rm -rf test.o simulation.dat

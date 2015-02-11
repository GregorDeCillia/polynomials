all:
	g++ polynomial_test.cpp -otest -lboost_unit_test_framework -std=c++11
	./test
doc:
	doxygen
clean:
	rm -rf test

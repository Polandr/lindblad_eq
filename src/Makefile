CXX= mpicxx
HDRS= complex_matrix.h complex_matrix.hpp complex_matrix_io.hpp scalapack.h

first: my_test

my_test: test.cpp $(HDRS)
	$(CXX) test.cpp -o my_test -lscalapack -llapack -lgfortran -lblas

clean:
	rm -f my_test

.PHONY: first clean

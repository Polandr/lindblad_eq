CXX= mpicxx
HDRS= solver.h solver.hpp processors.h processors.hpp complex_matrix.h complex_matrix.hpp complex_matrix_io.hpp scalapack.h quantum_system.h

first: my_test

my_test: test.cpp $(HDRS)
	$(CXX) test.cpp -o my_test -lscalapack -llapack -lblas -lgfortran

clean:
	rm -f my_test

.PHONY: first clean

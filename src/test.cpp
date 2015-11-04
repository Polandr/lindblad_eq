#include "complex_matrix.h"

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	Matrix A(8,8);
	A.generate(diagonal_natural_sequence);
	A.writef(1,"Matrix_1");

	Matrix B;
	B.readf("Matrix_1");
	B.writef(2,"Matrix_2");

	Matrix C = B;
	C.writef(2,"Matrix_3");

	Matrix D(C);
	D.writef(2,"Matrix_4");

	cout << D;

	ProcessorGrid::exit();
}
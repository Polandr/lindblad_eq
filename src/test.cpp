#include "complex_matrix.h"

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	Matrix A(8,8);
	A.generate(index_indicator);
	A.writef(2,"matrices/Matrix_1");

	Matrix B;
	B.readf("matrices/Matrix_1");
	B.writef(2,"matrices/Matrix_2");

	Matrix C = B;
	C.writef(2,"matrices/Matrix_3");

	Matrix D(C);
	D.writef(2,"matrices/Matrix_4");

	cout << D;

	ProcessorGrid::exit();
}
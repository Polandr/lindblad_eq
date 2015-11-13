#include "complex_matrix.h"

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	Matrix A(16,16);
	//A.generate(diagonal_natural_sequence);
	A.generate(identity_matrix);
	A.writef(1,"Matrix_1");

	Matrix B(16,16);
	//B.readf("Matrix_1");
	//B.generate(identity_matrix);
	B.generate(diagonal_natural_sequence);
	B.writef(2,"Matrix_2");
	Matrix C = A*B;
	C.writef(2,"Matrix_3");
	cout<<B<<endl;
	cout<<C<<endl;

	/*Matrix C = B;
	C.writef(2,"Matrix_3");

	Matrix D(C);
	D.writef(2,"Matrix_4");

	cout << D;*/

	//cout<<A<<endl;
	//cout<<B<<endl;

	/*Matrix C1 = A*B;
	cout<<C1<<endl;
	C1.writef(2,"Matrix_3");

	Matrix C2  = ~C1;
	cout<<C2<<endl;*/
	ProcessorGrid::exit();
}

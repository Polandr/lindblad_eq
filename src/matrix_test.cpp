#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	Matrix A(8,8), B(8,8);

	A.generate(identity_matrix);
	B.generate(index_indicator);

	complexd imag_unit(0,1);
	Matrix C;

	/*C = A;
	C *= 2;
	cout << C;
	ProcessorGrid::endline();
	C =  A * 2;
	cout << C;
	ProcessorGrid::endline();
	C =  2 * A;
	cout << C;
	ProcessorGrid::endline();

	C = A;
	C *= imag_unit;
	cout << C;
	ProcessorGrid::endline();
	C = A * imag_unit;
	cout << C;
	ProcessorGrid::endline();
	C = imag_unit * A;
	cout << C;
	ProcessorGrid::endline();

	C = A;
	C += A;
	cout << C;
	ProcessorGrid::endline();
	C = A + A;
	cout << C;
	ProcessorGrid::endline();

	C = B;
	C -= A;
	cout << C;
	ProcessorGrid::endline();
	C = B - A;
	cout << C;
	ProcessorGrid::endline();

	C = A * B;
	cout << C;
	ProcessorGrid::endline();
	C = B * A;
	cout << C;
	ProcessorGrid::endline();*/

	Matrix D(8,8);
	D.generate(identity_matrix);

	C = exp(D,1.0);
	cout << C;
	ProcessorGrid::endline();
	C = exp(D,imag_unit);
	cout << C;
	ProcessorGrid::endline();
	

	ProcessorGrid::exit();
}

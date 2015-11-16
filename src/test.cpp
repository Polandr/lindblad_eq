#include "complex_matrix.h"

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	//Matrix A(8,8);
	//A.generate(index_indicator);
	//A.writef(2,"Matrix_1");

	Matrix B(8,8);
	B.generate(diagonal_natural_sequence);

	/*vector<complexd> eigenvalues;
	Matrix A = B.diagonalize(eigenvalues);

	for (int i=0; i<A.global_n_rows(); i++)
	{
		cout<<eigenvalues[i]<<endl;
	}*/

	Matrix A = exp(B);
	A.writef(2,"MatrixD");

	ProcessorGrid::exit();
}

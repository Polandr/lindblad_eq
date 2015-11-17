#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	Matrix H(8,8), R0(8,8);
	H.generate(identity_matrix);
	R0.generate(identity_matrix);

	Solver solver;
	solver.init_H(H);
	solver.init_R0(R0);
	solver.init_dT();
	solver.init_step_num();
	cout << solver << endl << flush;

	//A.writef(2,"matrices/eigenvectors");

	//cout << A;

	ProcessorGrid::exit();
}
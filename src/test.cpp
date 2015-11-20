#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	/*Matrix H(8,8), R0(8,8);
	H.generate(identity_matrix);
	R0.generate(identity_matrix);

	Solver solver;
	solver.init_hamiltonian(H);
	solver.init_density_matrix(R0);
	solver.init_time_step(0.2);
	solver.init_step_num(4);

	solver.get_hamiltonian().writef(2,"matrices/hamiltonian");
	solver.get_density_matrix().writef(2,"matrices/init_density");

	solver.solve(NULL);*/

	Matrix A(8,8);

	A.generate(index_indicator);

	cout << A << flush;

	//cout << "A(3,4) = " << A(3,4) << endl << flush;

	ProcessorGrid::exit();
}
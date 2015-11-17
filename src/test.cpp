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
	solver.init_dT(0.2);
	solver.init_step_num(4);

	solver.get_H().writef(2,"matrices/hamiltonian");
	solver.get_R0().writef(2,"matrices/init_density");

	solver.solve(NULL);

	ProcessorGrid::exit();
}
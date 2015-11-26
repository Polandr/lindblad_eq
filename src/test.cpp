#include "complex_matrix.h"
#include "solver.h"

using namespace std;
complexd myH(int i, int j)
{
	if ((i==0)&&(j==1)) return 1;
	if ((i==1)&&(j==0)) return 1;
	if (i==j) return 1;
	return 0;
}

complexd myR(int i, int j)
{
	return ((i == 0)&&(j==0)) ? 1 : 0;
}

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	Matrix H(4,4), R0(4,4);
	H.generate(myH);
	R0.generate(myR);

	Solver solver;
	solver.init_hamiltonian(H);
	solver.init_density_matrix(R0);
	solver.init_time_step(0.2);
	solver.init_step_num(1);

	solver.get_hamiltonian().writef(2,"matrices/hamiltonian");
	solver.get_density_matrix().writef(2,"matrices/init_density");

	solver.solve(NULL);

	ProcessorGrid::exit();
}

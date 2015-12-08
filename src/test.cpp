#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	Matrix H(4,4), R0(8,8);
	H.generate(identity_matrix);
	R0.generate(identity_matrix);

	Lindblad_part lindblad;
	std::vector<complexd> di;
	di.push_back(0);
	di.push_back(0);
	//lindblad.init(1,di);

	//H.generate(test_H);
	//R0.generate(test_R);

	Solver solver;

	//solver.init_hamiltonian(H);
	vector<complexd> a(1), w(2);
	a[0] = 1;
	w[0] = 2; w[1] = 3;
	solver.init_hamiltonian(2,2,0,2,a,w);
	solver.init_density_matrix(R0);

	solver.init_time_step(0.2);
	solver.init_step_num(1);

	solver.init_lindblad(1,di);

	solver.solve(NULL);


	// Hamiltonian constructing test---------------------------------------------------------------

	/*vector<complexd> a(1), w(2);

	a[0] = 1;

	w[0] = 2;
	w[1] = 3;

	solver.init_hamiltonian(2,2,0,2,a,w);

	//solver.get_hamiltonian().writef(2,"matrices/hamiltonian");

	cout << solver.get_hamiltonian();*/

	ProcessorGrid::exit();
}

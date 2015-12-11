#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();


	std::vector<complexd> di;
	di.push_back(1);
	di.push_back(2);
	//lindblad.init(1,di);

	//H.generate(test_H);
	//R0.generate(test_R);

	Solver solver;

	//solver.init_hamiltonian(H);
	vector<complexd> a(1), w(2);
	a[0] = 1;
	w[0] = 2; w[1] = 3;
	solver.init_hamiltonian(2,2,0,2,a,w);
	solver.init_density_matrix(0);

	solver.init_time_step(1);
	solver.init_step_num(1);

	//solver.init_lindblad(3,di);

	cout << solver;

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

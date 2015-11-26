#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	//ProcessorGrid::square_init();
	ProcessorGrid::default_init();

	Matrix H(8,8), R0(8,8);
	H.generate(test_H);
	R0.generate(test_R);

	Solver solver;
	solver.init_hamiltonian(H);
	solver.init_density_matrix(R0);
	solver.init_time_step(0.5);
	solver.init_step_num(4);

	solver.get_hamiltonian().writef(2,"matrices/hamiltonian");
	solver.get_density_matrix().writef(2,"matrices/init_density");

	solver.solve(NULL);

	vector<complexd> eigenvalues;
	Matrix base_H = H.diagonalize(eigenvalues);

	//cout << H;
	//if (ProcessorGrid::is_root())
	//	cout << endl;

	/*cout << base_H;
	if (ProcessorGrid::is_root())
		cout << endl;

	if (ProcessorGrid::is_root())
	{
		for (int i = 0; i < eigenvalues.size(); i++)
			cout << eigenvalues[i] << endl;
		cout << endl;
	}*/


	ProcessorGrid::exit();
}

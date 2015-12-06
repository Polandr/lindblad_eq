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

	Matrix H(8,8), R0(8,8);

	//Lindblad_part lindblad;
	//std::vector<complexd> di;
	//complexd a(1.0,0);
	//complexd a1(3.0,0);
	//di.push_back(a);
	//di.push_back(a1);
	//lindblad.init(0,di);

	//H.generate(test_H);
	//R0.generate(test_R);

	Solver solver;

	/*solver.init_hamiltonian(H);
	solver.init_density_matrix(R0);

	solver.init_time_step(0.2);
	solver.init_step_num(1);

	solver.init_lindblad(0,di);

	solver.solve(NULL);


	vector<complexd> eigenvalues;
	Matrix base_H = H.diagonalize(eigenvalues);

	cout << H;
	if (ProcessorGrid::is_root())
		cout << endl;

	cout << base_H;
	if (ProcessorGrid::is_root())
		cout << endl;

	if (ProcessorGrid::is_root())
	{
		for (int i = 0; i < eigenvalues.size(); i++)
			cout << eigenvalues[i] << endl;
		cout << endl;
	}

	Matrix exp_H = exp(H,1.0);
	cout << exp_H;*/

	// Hamiltonian constructing test---------------------------------------------------------------

	vector<complexd> a(1), w(2);

	a[0] = 1;

	w[0] = 2;
	w[1] = 3;

	solver.init_hamiltonian(2,2,0,2,a,w);

	//solver.get_hamiltonian().writef(2,"matrices/hamiltonian");

	cout << solver.get_hamiltonian();

	ProcessorGrid::exit();
}

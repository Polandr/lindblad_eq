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

<<<<<<< HEAD
	Matrix H(4,4), R0(4,4);
	H.generate(myH);
	R0.generate(myR);
=======
	Solver solver;

	/*Matrix H(8,8), R0(8,8);
	H.generate(test_H);
	R0.generate(test_R);
>>>>>>> 3cee25457b2161fc6e3608b22938163077142391

	Solver solver;
	solver.init_hamiltonian(H);
	solver.init_density_matrix(R0);
<<<<<<< HEAD
	solver.init_time_step(0.2);
	solver.init_step_num(1);
=======
	solver.init_time_step(1);
	solver.init_step_num(4);
>>>>>>> 3cee25457b2161fc6e3608b22938163077142391

	//solver.get_hamiltonian().writef(2,"matrices/hamiltonian");
	//solver.get_density_matrix().writef(2,"matrices/init_density");

	solver.solve(NULL);

<<<<<<< HEAD
=======
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

	vector<complexd> a(2), w(3);

	a[0] = 1;
	a[1] = 2;

	w[0] = 1;
	w[1] = 2;
	w[2] = 3;

	solver.init_hamiltonian(3,0,0,3,a,w);

	//solver.get_hamiltonian().writef(2,"matrices/hamiltonian");

	cout << solver.get_hamiltonian();

>>>>>>> 3cee25457b2161fc6e3608b22938163077142391
	ProcessorGrid::exit();
}

/* 
Experiment about relation between:
	- Stock coefficient (Stock)
	- Number of corpuscles (N)
	- Evolution time (E_time)

Experiment model:

w_N       w_2       w_1      w_0
 o---...---o---a_1---o---a_0--o
                              |
                            Stock
                              |
                              o

Other parameters:
	- Energy level is 2 (both min and max)
	- Maximum stock energy is 2
	- Amplitudes between corpuscles are 1
	- Amplitudes on corpuscles are 1

	- Dephase coefficients are 0

	- Time step is 0.05

Initial density matrix:

Stop conditions:
	- Stock amplitude is more than 0.99
	or
	- Evolution time is more than 1000

*/

#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	/*vector<double> q_probs(2), s_probs(2);
	q_probs[0] = 1; q_probs[1] = 0;
	s_probs[0] = 1; s_probs[1] = 0;
	solver.init_density_matrix(q_probs,s_probs);*/

	int N_min = 7, N_max = 7, N_step = 1;
	double Stock_min = 1.0, Stock_max = 1.0, Stock_step = 0.02;


	for (int N = N_min; N <= N_max; N+=N_step)
		for (double Stock = Stock_min; Stock <= Stock_max; Stock+=Stock_step)
	{
		Solver solver;

		vector<complexd> a, w, d;
		for (int i = 0; i < N-1; i++)
			a.push_back(1);
		for (int i = 0; i < N; i++)
			w.push_back(1);
		for (int i = 0; i < N; i++)
			d.push_back(0);

		solver.init_hamiltonian(N,2,2,2,a,w);
		solver.init_density_matrix(0);
		solver.init_time_step(0.05);
		solver.init_step_num(1000.0/0.05);
		solver.init_lindblad(Stock,d);

		//cout << solver;
		ProcessorGrid::endline();

		double E_time = 0;
		E_time = solver.solve_to_max_stock();
		//solver.solve(NULL);

		if (ProcessorGrid::is_root())
			cout << "Evolution time is " << E_time << endl;

		ProcessorGrid::endline();

	}
	ProcessorGrid::exit();
}

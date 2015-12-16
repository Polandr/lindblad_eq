#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();


	

	//H.generate(test_H);
	//R0.generate(test_R);

	Solver solver;

	//solver.init_hamiltonian(H);
	vector<complexd> a(1), w(2);
	a[0] = 1;
	w[0] = 2; w[1] = 3;
	solver.init_hamiltonian(2,2,0,2,a,w);

	vector<double> q_probs(2), s_probs(3);
	q_probs[0] = 1; q_probs[1] = 0;
	s_probs[0] = 1; s_probs[1] = 0; s_probs[2] = 0;
	solver.init_density_matrix(q_probs,s_probs);

	//cout << solver.get_density_matrix();

	solver.init_time_step(0.1);
	solver.init_step_num(100);

	std::vector<complexd> d;
	d.push_back(1);
	d.push_back(1);
	solver.init_lindblad(1,d);

	//cout << solver;

	solver.solve(NULL);


	ProcessorGrid::exit();
}

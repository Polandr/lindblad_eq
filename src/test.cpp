#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	//H.generate(test_H);
	//R0.generate(test_R);

	Solver solver;

	int step = 10;
	double st = 0.001;

	//solver.init_hamiltonian(H);
	/*vector<complexd> a(1), w(2);
	a[0] = 1;
	w[0] = 1; w[1] = 1;
	solver.init_hamiltonian(2,1,0,1,a,w);*/

	vector<complexd> a(1), w(2);
	for (int i=0; i<1; i++)
	{
		a[i] = 1;
	}

	for (int i=0; i<2; i++)
	{
		w[i] = 1;
	}

	solver.init_hamiltonian(2,1,0,1,a,w);

	/*vector<double> q_probs(2), s_probs(2);
	q_probs[0] = 1; q_probs[1] = 0;
	s_probs[0] = 1; s_probs[1] = 0;// s_probs[2] = 0;
	solver.init_density_matrix(q_probs,s_probs);*/

	solver.init_density_matrix(1);

	//cout << solver.get_density_matrix();

	solver.init_time_step(st);
	solver.init_step_num(step);

	std::vector<complexd> d;
	for (int i=0; i<2; i++)
		d.push_back(1);
	solver.init_lindblad(1,d);

	//cout << solver;
	//solver.print_base_states(cout);

	solver.solve(NULL);
	
	//cout<<step*st<<endl;


	ProcessorGrid::exit();
}

#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	int step = 10000000;
	double st = 0.003;

	vector<complexd> a(6), w(7);
		for (int i=0; i<6; i++)
		{
			a[i] = 1;
		}
	
	w[0] = 1; w[1] = 2;
	w[2] = 1; w[3] = 2;
	w[4] = 1; w[5] = 2;
	w[6] = 1;

	double op = 0.05;
	double now = 0.0;
	double stockStep = 0.2;
	double stockKoeff = 2.2;
	int steps = 2/op + 1;

/********************************get time*******************************************************/

	for (int j=0; j<10; j++)
	{
		for (int i=0; i<steps; i++)
		{
			Solver solver;
			solver.init_hamiltonian(7,1,1,1,a,w);
			solver.init_density_matrix(1);
			solver.init_time_step(st);
			solver.init_step_num(step);

			std::vector<complexd> d;
			for (int i=0; i<7; i++)
				d.push_back(now); 
			solver.init_lindblad(stockKoeff,d);

			now+=op; 

			solver.solve(NULL);
		}
	cout<<"-------------------------------------------------"<<endl;
	stockKoeff +=stockStep;
	now = 0.0;
	}
/***************************end**************************************************************/
	ProcessorGrid::exit();
}

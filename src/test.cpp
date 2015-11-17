#include "complex_matrix.h"
#include "solver.h"

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	int N = 3;

	/*Matrix H(8,8);
	H.generate(random_lower_triangle);
	H.writef(2,"Matrix_H");

	Matrix R(8,8);
	R.generate(diagonal_natural_sequence);
	R.writef(2,"Matrix_R");*/

	Solver sol;
	sol.init_H("Matrix_H");
	sol.init_R0("Matrix_R");
	sol.init_dT(0.1);
	sol.init_step_num(N);

	Matrix H = sol.get_H();
	Matrix R = sol.get_R0();
	double dT = sol.get_dT();
	int step = sol.get_step_num();

	sol.solve("Matrix_out");

	ProcessorGrid::exit();
}

#include "complex_matrix.h"
#include "solver.h"

using namespace std;

int main(int argc, char** argv)
{
	ProcessorGrid::default_init();

	Matrix R(8,8);
	R.generate(test_R);

	cout << R;
	ProcessorGrid::endline();

	herm_normalization(R);

	cout << R;	

	ProcessorGrid::exit();
}

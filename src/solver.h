#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <iostream>
#include <cstdio>
#include <utility>
#include <fstream>
#include <vector>

#include "complex_matrix.h"
#include "quantum_system.h"


class Solver
{
	Matrix H;
	Matrix R;

	std::vector<int> base_states;
	std::vector<int> state_nums;
	int N;
	
	Lindblad_part L;
	double dT;
	int step_n;

public:

	Matrix& get_hamiltonian () { return H; }
	Matrix& get_density_matrix () { return R; }
	double& get_time_step () { return dT; }
	int& get_step_num () { return step_n; }

	void init_dimension(int dim) { N = dim; }

	void init_hamiltonian (const char*);
	void init_hamiltonian (const Matrix&);
	void init_hamiltonian (int, int, int, int, std::vector<complexd>, std::vector<complexd>);

	void init_density_matrix (const char*);
	void init_density_matrix (const Matrix&);
	void init_density_matrix (vector<complexd>);
	void init_density_matrix (int);

	void init_lindblad (complexd, std::vector<complexd>);
	void init_time_step (double);
	void init_step_num (int);

	void init_system ();

	void solve (const char*);

	void print_base_states(std::ostream&);
	void operator >> (std::ostream&);
	void operator << (std::istream&);
};

std::ostream& operator << (std::ostream&, Solver&);
std::istream& operator >> (std::istream&, Solver&);


#include "solver.hpp"

#endif

// __SOLVER_H__ 

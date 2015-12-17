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

	void init_hamiltonian (const char* filename);
	void init_hamiltonian (const Matrix& matrix_H);
	void init_hamiltonian (int sys_dim, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w);
	// Make hamiltonian for definite system with parameters:
	// sys_dim - dimension of system
	// s - maximum stock level
	// E_min, E_max - minimum and maximum energy levels
	// a - probabilities between atoms
	// w - probabilities on atoms

	void init_density_matrix (const char* filename);
	void init_density_matrix (const Matrix& matrix_R);
	void init_density_matrix (vector<complexd> state);
	// Make density matrix for quantum state
	void init_density_matrix (int pos);
	// Make density matrix with single 1 on main diagonal on position number <pos> and other elemants set as 0
	// Requires base_states initialization
	void init_density_matrix (vector<double> qbit_probs, vector<double> stock_probs);
	// Make density matrix with probabilities of {qbit = 1} event as <qbit_probs> and {stock = i} event as <stock_probs>
	// If <qbit_probs> designate set of possible states (at least 1): states which were set during hamiltonian construction
	// then set of possible states are chosen with certain ratio
	// If <qbit_probs> designate set of impossible states
	// then all states are set as equiprobable
	// Requires base_states initialization

	void init_time_step (double time);
	void init_step_num (int steps);
	void init_lindblad (complexd stock, std::vector<complexd> l);

	void init_system ();
	// Unsafe!
	// Some default initialization

	void solve (const char* filename);

	void print_base_states(std::ostream&);
	void operator >> (std::ostream&);
	void operator << (std::istream&);
};

std::ostream& operator << (std::ostream&, Solver&);
std::istream& operator >> (std::istream&, Solver&);


#include "solver.hpp"

#endif

// __SOLVER_H__ 

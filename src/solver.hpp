// Solver class realization

#include <cmath>
#include <cstdio>

#define DEFAULT_H_FILE "Matrix_H"
#define DEFAULT_R_FILE "Matrix_R"
#define DEFAULT_DT 0.1
#define DEFAULT_STEP_NUM 1

#ifndef EPS
#define EPS 0.000001
#endif

using namespace std;

#include "exceptions.hpp"

// Service functions:

void normalize_probabilities_with_warning(vector<double>& probs)
{
	double probs_sum = 0;
	for (int i = 0; i < probs.size(); ++i)
		probs_sum += probs[i];
	if (fabs(probs_sum - 1) >= EPS)
	// If sum of probabilities is not equal to 1 then normalize probabilities
	{
		if (ProcessorGrid::is_root())
			cout << "Incorrect probability sequence, normalizing it.\n";
		for (int i = 0; i < probs.size(); ++i)
			if (probs_sum != 0)
				probs[i] /= probs_sum;
			else
				probs[i] = 1/probs.size();

	}
}

void print_header (FILE* file)
{
	if (ProcessorGrid::is_root())
		fprintf(file, "Magnitudes of diagonal elements are:\n");
}

// Solver functions realization------------------------------------------------------------

// Hamiltonian initialization:

void Solver::init_hamiltonian (const char* filename = DEFAULT_H_FILE)
{
	H.readf(filename);
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian initialization");
}

void Solver::init_hamiltonian (const Matrix& matrix_H)
{
	H = matrix_H;
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian initialization");
}

void Solver::init_hamiltonian (int sys_dim, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w)
// sys_dim - dimension of system
// s - maximum stock level
// E_min, E_max - minimum and maximum energy levels
// a - probabilities between atoms
// w - probabilities on atoms
{
	N = sys_dim;
	if (a.size() != N-1 || w.size() != N)
		throw Solver_exception("incorrect parameters in hamiltonian initialization");
	H = hamiltonian(N, s, E_min, E_max, a, w, base_states, state_nums);
}

// Initial density matrix initialization:

void Solver::init_density_matrix (const char* filename = DEFAULT_R_FILE)
{
	R.readf(filename);
	if (!(R.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix initialization");
}

void Solver::init_density_matrix (const Matrix& matrix_R)
{
	R = matrix_R;
	if (!(R.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix initialization");
}

void Solver::init_density_matrix (vector<complexd> state)
{
	R = density_matrix(state);
}

void Solver::init_density_matrix (int pos)
{
	if (pos >= base_states.size())
		throw  Solver_exception("invalid density matrix initialization");
	R = density_matrix(base_states.size(), pos);
}

void Solver::init_density_matrix (vector<double> qbit_probs, vector<double> stock_probs)
{
	if (qbit_probs.size() != N || stock_probs.size() != state_nums.size())
		throw Solver_exception("incorrect parameters in initital density matrix initialization");

	for (int i = 0; i < qbit_probs.size(); ++i)
		if (qbit_probs[i] > 1)
			throw Solver_exception("incorrect parameters in initital density matrix initialization");
	normalize_probabilities_with_warning(stock_probs);	

	R = density_matrix(qbit_probs, stock_probs, base_states, state_nums);
}

// Other parameters initialization:

void Solver::init_lindblad (complexd out, std::vector<complexd> l)
{
	if (l.size() != N)
		throw Solver_exception("incorrect parameters in lindblad initialization");
	L.init(out, l, base_states, state_nums);
}

void Solver::init_base_states ()
{
	base_states.resize(0);
	state_nums.resize(0);
	for (int i = 0; i < H.global_n_rows(); i++)
		base_states.push_back(i);
	state_nums.push_back(base_states.size());
	int log_size = static_cast<int>(ceil(log(base_states.size())/log(2)));
	N = log_size;
}

void Solver::init_time_step (double dt = DEFAULT_DT)
{
	dT = dt;
}

void Solver::init_step_num (int steps = DEFAULT_STEP_NUM)
{
	step_n = steps;
}

void Solver::init_system ()
{
	init_hamiltonian(DEFAULT_H_FILE);
	init_base_states();
	init_density_matrix(DEFAULT_R_FILE);
	init_time_step(DEFAULT_DT);
	init_step_num(DEFAULT_STEP_NUM);
}

void Solver::clear_system ()
{
	//H.~Matrix();
	//R.~Matrix();
	H.init();
	R.init();
	base_states.clear();
	base_states.resize(0);
	state_nums.clear();
	state_nums.resize(0);
	N = 0;
	dT = 0;
	step_n = 0;
	L.destroy();
}

// Main functions-------------------------------------------------------------------------

void Solver::solve (const char* filename)
{
	FILE* file;
	if (filename != NULL)
	{
		file = fopen(filename,"w");
		print_header(file);
	}
	else
		print_header(stdout);

	complexd imag_unit(0,1);
	Matrix U = exp(H,(-imag_unit)*dT/Plank_const);
	Matrix conj_U = U.herm_conj();

	for (int i = 0; i < step_n; i++)
	{
		R = U*R;
		R = R*conj_U;

		if (L.active)
			R += dT*L(R);

		if (filename != NULL)
			R.print_diagonal_abs(file);
		else
			R.print_diagonal_abs(stdout);

		/*if (ProcessorGrid::is_root())
			cout << "Time of evolution is " << (i*dT) << endl << flush;
		int last_elem = R.global_n_rows()-1;
		if (abs(R(last_elem,last_elem)) >= 0.99)
			break;*/
	}

	if (filename != NULL)
		fclose(file);
}

// Solver provides evolution till all energy flow down to stock
// Returns time of evolution
const double MAX_STOCK_LVL = 0.99;
const double MAX_TIME = 1000.0;
double Solver::solve_to_max_stock ()
{
	double evolution_time = 0;

	complexd imag_unit(0,1);
	Matrix U = exp(H,(-imag_unit)*dT/Plank_const);
	Matrix conj_U = U.herm_conj();

	for (bool done = false; (!done && evolution_time < MAX_TIME);)
	{
		R = U*R;
		R = R*conj_U;

		if (L.active)
			R += dT*L(R);

		evolution_time += dT;

		// State describing that stock has maximum energy is always in last block of density matrix
		// So to get amplitude of maximum energy in stock last block should be summarized
		int last_block = state_nums.size() - 1;
		double amplitude = summarize_amplitudes_in_stock_block(R, last_block, base_states, state_nums);
		if (amplitude >= MAX_STOCK_LVL)
			done = true;
	}

	return evolution_time;
}

// I/O-------------------------------------------------------------------------

void Solver::print_base_states(ostream& out)
{
	if (state_nums.size() == 0)
		ProcessorGrid::root_print("<No base states>\n");
	if (ProcessorGrid::is_root())
		for (int i = 0, ofs = 0; i < state_nums.size(); ++i)
		{
			out << "Stock is " << i << ":\n";
			for (int j = 0; j < state_nums[i]; ++j, ++ofs)
				print_ketbra_stream(base_states[ofs], N, out);
			out << endl;
		}
	out << flush;
}

void Solver::operator >> (ostream& out)
{
	ProcessorGrid::root_print("System configuration is:\n");

	ProcessorGrid::root_print("Base states:\n");
	print_base_states(out);

	ProcessorGrid::root_print("Matrix H:\n");
	out << get_hamiltonian();

	ProcessorGrid::root_print("Matrix R:\n");
	out << get_density_matrix();

	if (ProcessorGrid::is_root())
		out << endl << "dT: " << get_time_step();
	if (ProcessorGrid::is_root())
		out << endl << "step number: " << get_step_num() << endl << endl;

	if (ProcessorGrid::is_root() && L.active)
		out << "Lindblad:\n";
	if (L.active)
		L.print_matrices(out);

	out << flush;
}
void Solver::operator << (istream& in)
{
	in >> get_hamiltonian();
	in >> get_density_matrix();
	in >> get_time_step();
	in >> get_step_num();
	init_base_states();
}

ostream& operator << (ostream& out, Solver& src)
{
	src >> out;
	return out;
}

istream& operator >> (istream& in, Solver& trg)
{
	trg << in;
	return in;
}

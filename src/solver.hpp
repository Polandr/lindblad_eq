#include <complex>
#include <cmath>

#define DEFAULT_H_FILE "Matrix_H"
#define DEFAULT_R_FILE "Matrix_R"
#define DEFAULT_DT 0.1
#define DEFAULT_STEP_NUM 1

#ifndef EPS
#define EPS 0.000001
#endif

using namespace std;

// Solver exception class-----------------------------------------------------------------

class Solver_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Solver_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Solver_exception() throw()
	{
		delete [] errstr;
	}
	virtual const char* what() const throw()
	{
		char* tmp = errstr;
		char* prefix = const_cast <char*> ("Solver class error: ");
		try
		{
			errstr = new char [strlen(prefix)+strlen(errstr)+2];
		}
		catch (std::exception& smth)
		{
			return "Couldn't generate an error message (there is no memory)\n";
		}
		sprintf(errstr, "%s%s.\n", prefix, tmp);
		return errstr;
	}	
};

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
	base_states.resize(0);
	for (int i = 0; i < H.global_n_rows(); i++)
		base_states.push_back(i);
	state_nums.push_back(base_states.size());
}

void Solver::init_hamiltonian (const Matrix& matrix_H)
{
	H = matrix_H;
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian initialization");
	base_states.resize(0);
	for (int i = 0; i < H.global_n_rows(); i++)
		base_states.push_back(i);
	state_nums.push_back(base_states.size());
}

void Solver::init_hamiltonian (int sys_dim, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w)
// sys_dim - dimension of system
// s - maximum stock level
// E_min, E_max - minimum and maximum energy levels
// a - probabilities between atoms
// w - probabilities on atoms
{
	init_dimension(sys_dim);
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
	L.active = true;
	L.init(out, l, base_states, state_nums);
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
	init_density_matrix(DEFAULT_R_FILE);
	init_time_step(DEFAULT_DT);
	init_step_num(DEFAULT_STEP_NUM);
}

// Main function-------------------------------------------------------------------------

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

	}
	if (filename != NULL)
		fclose(file);
}

// I/O :

void Solver::print_base_states(ostream& out)
{
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
	if (ProcessorGrid::is_root())
		out << "System configuration is:\n";

	if (ProcessorGrid::is_root())
		out << "Base states:\n";
	print_base_states(out);

	if (ProcessorGrid::is_root())
		out << "Matrix H:\n";
	out << get_hamiltonian();
	if (ProcessorGrid::is_root())
		out << endl << "Matrix R:\n";
	out << get_density_matrix();
	if (ProcessorGrid::is_root())
		out << endl << "dT: " << get_time_step();
	if (ProcessorGrid::is_root())
		out << endl << "step number: " << get_step_num() << endl;

	if (ProcessorGrid::is_root())
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

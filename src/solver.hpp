#include <complex>
#include <cmath>

#define DEFAULT_H_FILE "Matrix_H"
#define DEFAULT_R_FILE "Matrix_R"
#define DEFAULT_DT 0.1
#define DEFAULT_STEP_NUM 1

const double Plank_const = 1.0;

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

// Solver functions realization------------------------------------------------------------

// Hamiltonian initialization:

void Solver::init_hamiltonian (const char* filename = DEFAULT_H_FILE)
{
	H.readf(filename);
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian");
	base_states.resize(0);
	for (int i = 0; i < H.global_n_rows(); i++)
		base_states.push_back(i);
}

void Solver::init_hamiltonian (const Matrix& matrix_H)
{
	H = matrix_H;
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian");
	base_states.resize(0);
	for (int i = 0; i < H.global_n_rows(); i++)
		base_states.push_back(i);
}

void Solver::init_hamiltonian (int N, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w)
// N - dimension of system
// s - maximum stock level
// E_min, E_max - minimum and maximum energy levels
// a - probabilities between atoms
// w - probabilities on atoms
{
	if (a.size() != N-1 || w.size() != N)
		throw Solver_exception("incorrect parameters in hamiltonian initialization");
	H = hamiltonian(N, s, E_min, E_max, a, w, base_states);
}

// Initial density matrix initialization:

void Solver::init_density_matrix (const char* filename = DEFAULT_R_FILE)
{
	R.readf(filename);
	if (!(R.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix");
}

void Solver::init_density_matrix (const Matrix& matrix_R)
{
	R = matrix_R;
	if (!(R.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix");
}

void Solver::init_density_matrix (vector<complexd> state)
{
	R = density_matrix(state);
}

// Other parameters initialization:

void Solver::init_lindblad (complexd out, std::vector<complexd> ls, std::vector<complexd> Ls)
{
	L.active = true;
	L.init(out, ls, Ls);
}

void Solver::init_time_step (double dt = DEFAULT_DT)
{
	dT = dt;
}

void Solver::init_step_num (int steps = DEFAULT_STEP_NUM)
{
	step_num = steps;
}

void Solver::init_system ()
{
	init_hamiltonian(DEFAULT_H_FILE);
	init_density_matrix(DEFAULT_R_FILE);
	init_time_step(DEFAULT_DT);
	init_step_num(DEFAULT_STEP_NUM);
}

void print_header (FILE* file)
{
	if (ProcessorGrid::is_root())
		fprintf(file, "Magnitudes of diagonal elements are:\n");
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

	for (int i = 0; i < step_num; i++)
	{
		Matrix dL;
		if (L.active)
			dL = dT*L(R,base_states);
		R = U*R;
		R = R*conj_U;
		if (L.active)
			R += dL;
		if (filename != NULL)
			R.print_diagonal_abs(file);
		else
			R.print_diagonal_abs(stdout);
	}
	if (filename != NULL)
		fclose(file);
}

ostream& operator << (ostream& out, Solver& src)
{
	out << "System configuration is:\n";
	out << "Matrix H:\n" << src.get_hamiltonian() << endl;
	out << "Matrix R:\n" << src.get_density_matrix() << endl;
	out << "dT: " << src.get_time_step() << endl;
	out << "step number: " << src.get_step_num() << endl;

	return out;
}

istream& operator >> (istream& in, Solver& trg)
{
	in >> trg.get_hamiltonian();
	in >> trg.get_density_matrix();
	in >> trg.get_time_step();
	in >> trg.get_step_num();

	return in;
}

#include <complex>

#define DEFAULT_H_FILE "Matrix_H"
#define DEFAULT_R0_FILE "Matrix_R0"
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

void Solver::init_hamiltonian(const char* filename = DEFAULT_H_FILE)
{
	H.readf(filename);
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian");
}

void Solver::init_density_matrix(const char* filename = DEFAULT_R0_FILE)
{
	R0.readf(filename);
	if (!(R0.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix");
}

void Solver::init_hamiltonian(const Matrix& matrix_H)
{
	H = matrix_H;
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian");
}

void Solver::init_density_matrix(const Matrix& matrix_R0)
{
	R0 = matrix_R0;
	if (!(R0.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix");
}

void Solver::init_time_step(double dt = DEFAULT_DT)
{
	dT = dt;
}

void Solver::init_step_num(int steps = DEFAULT_STEP_NUM)
{
	step_num = steps;
}

void Solver::init_system()
{
	init_hamiltonian(DEFAULT_H_FILE);
	init_density_matrix(DEFAULT_R0_FILE);
	init_time_step(DEFAULT_DT);
	init_step_num(DEFAULT_STEP_NUM);
}

void print_header(FILE* file)
{
	if (ProcessorGrid::is_root())
		fprintf(file, "Magnitudes of diagonal elements are:\n");
}

void Solver::solve(const char* filename)
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
	Matrix U = exp(H*((-imag_unit)*dT/Plank_const));
	Matrix conj_U = ~U;
	Matrix Rt = R0;

	for (int i = 0; i < step_num; i++)
	{
		Rt = U*Rt;
		Rt = Rt*conj_U;
		if (filename != NULL)
			Rt.print_diagonal_abs(file);
		else
			Rt.print_diagonal_abs(stdout);
	}
	if (filename != NULL)
		fclose(file);
}

ostream& operator << (ostream& out, Solver& src)
{
	out << "System configuration is:\n";
	out << "Matrix H:\n" << src.get_hamiltonian() << endl;
	out << "Matrix R0:\n" << src.get_density_matrix() << endl;
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

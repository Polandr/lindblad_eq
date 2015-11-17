

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

void Solver::init_H(const char* filename)
{
	H.readf(filename);
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian");
}

void Solver::init_R0(const char* filename)
{
	R0.readf(filename);
	if (!(R0.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix");
}

void Solver::init_dT(double dt)
{
	dT=dt;
}

void Solver::init_step_num(int steps)
{
	step_num=steps;
}

void Solver::init_system()
{
	dT = DEFAULT_DT;
	step_num = DEFAULT_STEP_NUM;
	H.readf(DEFAULT_H_FILE);
	if (!(H.is_square()))
		throw Solver_exception("incorrect matrix dimensions in hamiltonian");

	R0.readf(DEFAULT_R0_FILE);
	if (!(R0.is_square()))
		throw Solver_exception("incorrect matrix dimensions in initital density matrix");
}

void Solver::solve(const char* filename)
{
	ofstream file;
	double i = sqrt(-1);
	if (filename != NULL)
		file.open(filename, ios::out);
	/*else
		file = cout;*/

	//using namespace std::literals;
	Matrix H0 = H*(-i*dT/Plank_const);
	Matrix U = exp(H0); // This function need to be overloaded
	Matrix conj_U = ~U; // This function need to be realized
	Matrix Rt = R0;

	for (int i = 0; i < step_num; i++)
	{
		Rt = U*Rt;
		Rt = Rt*conj_U; /* ScaLapack multiplication U*Rt*conj_U */


		//TODO: Fix printing abs in file


		//Rt.print_on_condition(file,diagonal);
	}
}

/*ostream& operator << (ostream& out, const Solver& src)
{
	out << "System configuration is:\n"
	out << "Matrix H:\n" << src.get_H() << endl;
	out << "Matrix R0:\n" << src.get_R0() << endl;
	out << "dT: " << src.get_dT() << endl;
	out << "step number: " << src.get_step_num() << endl;

	return out;
}

istream& operator >> (istream& in, const Solver& trg)
{
	in >> trg.get_H();
	in >> trg.get_R0();
	in >> trg.get_dT();
	in >> trg.get_step_num();

	return in;
}*/



// Need:
//
// Matrix& exp(Matrix&)

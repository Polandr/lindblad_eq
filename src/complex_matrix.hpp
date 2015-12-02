// Matrix class realization

using namespace std;

const double Plank_const = 1.0;
#define ONE_VAL_TAG 512

class Matrix_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Matrix_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Matrix_exception() throw()
	{
		delete [] errstr;
	}
	virtual const char* what() const throw()
	{
		char* tmp = errstr;
		char* prefix = const_cast <char*> ("Complex matrix error: ");
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

// Matrix grid initialization

void Matrix::init_distribution (int rows = 0, int cols = 0,
	int row_block = R_BLOCK_SIZE, int col_block = C_BLOCK_SIZE)
{
	if (!ProcessorGrid::initialized)
	{
		ProcessorGrid::default_init();
	}

	info.set_matrix_dims(rows,cols);
	info.set_block_sizes(row_block,col_block);

	int lld = max(1,info.local_row_num());
	int out = 0;
	descinit_((info.descriptor),
			  &(info.matrix_global_rows), &(info.matrix_global_cols), &row_block, &col_block,
			  &(ProcessorGrid::start_proc_row), &(ProcessorGrid::start_proc_col),
			  &(ProcessorGrid::context), &lld, &out);
	if (out != 0)
		throw Matrix_exception("incorrect descriptor initialization");
}

// Service functions-----------------------------------------------------------------

void Matrix::create ()
{
	data = new complexd [n_rows*n_cols];
}

void Matrix::destroy()
{
	delete [] data;
}

bool Matrix::in_block(int i, int j) const
{
	return ((i >= info.row_offset()) &&
			(i < info.row_offset() + n_rows) &&
			(j >= info.col_offset()) &&
			(j < info.col_offset() + n_cols));
}
		
void Matrix::set (int i, int j, complexd val)
// Global set
{
	if (in_block(i,j))
		data[(i - info.row_offset()) + n_rows * (j - info.col_offset())] = val;
}

void Matrix::set (int i, complexd val)
{
	data[i] = val;
}

const complexd Matrix::get (int i) const
{
	return data[i];
}

double* Matrix::get_data () const
{
	double* array = (double*)malloc(2*n_rows*n_cols*sizeof(double));
	for (int i = 0; i < n_rows*n_cols; i++)
	{
		array[2*i] = data[i].real();
		array[2*i+1] = data[i].imag();
	}
	return array;
}

void Matrix::get_row (double* array, int row) const
// The row is got according to column-major order
// so elements are taken fragmently
{
	for (int i = 0; i < n_cols; i++)
	{
		array[2*i] = data[row+i*n_rows].real();
		array[2*i+1] = data[row+i*n_rows].imag();
	}
}

void Matrix::set_data (double* array)
{
	for (int i = 0; i < n_rows*n_cols; i++)
	{
		data[i].real() = array[2*i];
		data[i].imag() = array[2*i+1];
	}
}

void Matrix::set_row (double* array, int row)
// The row is set according to column-major order
// so elements are set fragmently
{
	for (int i = 0; i < n_cols; i++)
	{
		data[row+i*n_rows].real() = array[2*i];
		data[row+i*n_rows].imag() = array[2*i+1];
	}
}

// Constructors and destructor---------------------------------------------------------

Matrix::Matrix ()
{
	init_distribution();
	n_rows = 0;
	n_cols = 0;
}
   
Matrix::Matrix (int rows, int cols,
	int row_block = R_BLOCK_SIZE, int col_block = C_BLOCK_SIZE)
{
	if (rows < 0 || cols < 0)
		throw Matrix_exception("invalid matrix size");
	init_distribution(rows,cols,row_block,col_block);
	n_rows = info.local_row_num();
	n_cols = info.local_col_num();
	create();
	for (int i = 0; i < n_rows*n_cols; i++)
		data[i] = 0;
}

void Matrix::init (int rows, int cols,
	int row_block = R_BLOCK_SIZE, int col_block = C_BLOCK_SIZE)
{
	this->~Matrix();
	if (rows < 0 || cols < 0)
		throw Matrix_exception("invalid matrix size");
	init_distribution(rows,cols,row_block,col_block);
	n_rows = info.local_row_num();
	n_cols = info.local_col_num();
	create();
	for (int i = 0; i < n_rows*n_cols; i++)
		data[i] = 0;
}

Matrix::Matrix (const Matrix& other)
{
	info = other.get_distribution();
	n_rows = other.n_rows;
	n_cols = other.n_cols;
	create();
	memcpy(data, other.data, n_rows*n_cols*sizeof(complexd));
}
		
Matrix::~Matrix ()
{
	if (n_rows != 0 && n_cols != 0)
		destroy();
}

// Some necessary stuff-----------------------------------------------------------------

Matrix& Matrix::operator = (const Matrix& other)
{
	this->~Matrix();
	info = other.get_distribution();
	n_rows = other.n_rows;
	n_cols = other.n_cols;
	create();
	//memcpy(data, other.data, n_rows*n_cols*sizeof(complexd));
	for (int i = 0; i < n_rows*n_cols; i++)
		data[i] = other.data[i];

	return *this;
}

void Matrix::in_place_transposition ()
{
	for (int i = 1; i < (n_rows*n_cols-1); i++)
	{
		vector<int> idx_group(1,i);
		bool swapped = false;
		for (int idx = i; idx*n_cols % (n_rows*n_cols-1) != i;)
		{
			idx = idx*n_cols % (n_rows*n_cols-1);
			idx_group.push_back(idx);
			if (idx < i)
			{
				swapped = true;
				break;
			}
		}
		if (!swapped)
		{
			complexd tmp = data[idx_group.back()];
			for (int j = 0; j < idx_group.size()-1; j++)
				data[idx_group[j+1]] = data[idx_group[j]];
			data[idx_group[0]] = tmp;
		}
	}
}

const complexd Matrix::operator () (int row, int col) const
// Global read-only access to the element
{
	if (row >= global_n_rows() || row < 0 || col >= global_n_cols() || col < 0) 
		throw Matrix_exception("out of bounds");
	else
	{
		double value[2];
		if (in_block(row,col))
		{
			value[0] = data[(row - info.row_offset()) + n_rows * (col - info.col_offset())].real();
			value[1] = data[(row - info.row_offset()) + n_rows * (col - info.col_offset())].imag();
			MPI_Send(value,2,MPI_DOUBLE,ProcessorGrid::root,ONE_VAL_TAG,MPI_COMM_WORLD);
		}
		if (ProcessorGrid::is_root())
		{
			MPI_Status status;
			MPI_Recv(value,2,MPI_DOUBLE,MPI_ANY_SOURCE,ONE_VAL_TAG,MPI_COMM_WORLD,&status);
		}
		MPI_Bcast(value,2,MPI_DOUBLE,ProcessorGrid::root,MPI_COMM_WORLD);
		return complexd(value[0],value[1]);
	}
}

// Arithmetics----------------------------------------------------------------------------

Matrix& Matrix::operator *= (const complexd val)
{
	for (int i = 0; i < n_rows*n_cols; i++)
		data[i] *= val;
	return *this;
}

Matrix Matrix::operator * (const complexd val) const
{
	Matrix out(*this);
	out *= val;
	return out;
}

Matrix operator * (const complexd val, const Matrix& mat)
{
	return  mat*val;
}

Matrix& Matrix::operator += (const Matrix& A)
{
	if (global_n_rows() != A.global_n_rows() || global_n_cols() != A.global_n_cols())
		throw Matrix_exception("incomatible sizes of matrices in summation");
	for (int i = 0; i < n_rows*n_cols; i++)
		data[i] += A.get(i);
	return *this;
}

Matrix Matrix::operator + (const Matrix& A) const
{
	Matrix out(*this);
	out += A;
	return out;
}

Matrix& Matrix::operator -= (const Matrix& A)
{
	if (global_n_rows() != A.global_n_rows() || global_n_cols() != A.global_n_cols())
		throw Matrix_exception("incomatible sizes of matrices in summation");
	for (int i = 0; i < n_rows*n_cols; i++)
		data[i] -= A.get(i);
	return *this;
}

Matrix Matrix::operator - (const Matrix& A) const
{
	Matrix out(*this);
	out -= A;
	return out;
}


Matrix Matrix::operator * (Matrix& b) const
{
	Matrix c_matr(global_n_rows(), b.global_n_cols());

	int m = global_n_rows();
	int n = b.global_n_cols();
	int k = global_n_cols();

	Distribution distr = get_distribution();
	Distribution distrB = b.get_distribution();
	Distribution distrC = c_matr.get_distribution();

	int row_offset = 1;
	int col_offset = 1;

	double* a_data = get_data();
	double* b_data = b.get_data();
	double* c_data = c_matr.get_data();
	double alpha[2];
	double beta[2];
	alpha[0] = 1.0;
	alpha[1] = 0.0;
	beta[0] = 0.0; 
	beta[1] = 0.0;

	pzgemm_((char*) "N", (char*) "N", &m, &n, &k, 
		alpha, a_data, &row_offset, &col_offset, distr.descriptor, 
		b_data, &row_offset, &col_offset, distrB.descriptor, 
		beta, c_data, &row_offset, &col_offset, distrC.descriptor);

	c_matr.set_data(c_data);

	return c_matr;
}

Matrix Matrix::operator ~ () const
{
	Matrix c_matr(global_n_cols(), global_n_rows());

	Distribution distrA = get_distribution();
	Distribution distrC = c_matr.get_distribution();

	int row_offset = 1;
	int col_offset = 1;

	int m = global_n_rows();
	int n = global_n_cols();

	double alpha[2];
	double beta[2];
	alpha[0] = 1.0;
	alpha[1] = 0.0;
	beta[0] = 0.0; 
	beta[1] = 0.0;

	double* a_data = get_data();
	double* c_data = c_matr.get_data();
	
	pzgeadd_((char*) "C", &m, &n, 
		alpha, a_data, &row_offset, &col_offset, distrA.descriptor, 
		beta, c_data, &row_offset, &col_offset, distrC.descriptor);

	c_matr.set_data(c_data);

	return c_matr;
}

Matrix Matrix::conj () const
{
	Matrix out(*this);
	for (int i = 0; i < n_rows*n_cols; i++)
		out.set(i,std::conj(data[i]));
	return out;
}

Matrix Matrix::herm_conj () const
{
	Matrix out = this->conj();
	return ~out;
}

Matrix Matrix::diagonalize (vector<complexd>& eigenvalues) const
{
	double *w = (double*)malloc(global_n_rows()*sizeof(double));
	Matrix Z(global_n_rows(), global_n_rows());

	Distribution distrA = get_distribution();
	Distribution distrZ = Z.get_distribution();

	int n = global_n_rows();

	double* a = get_data();
	double* z = Z.get_data();
	int row_offset = 1;
	int col_offset = 1;

	int lrwork = 2*n + 2*n-2;
	int lwork = -1;

	double* work = (double*)malloc(2*sizeof(double));
	double* rwork  = (double*)malloc(lrwork*sizeof(double));
	int ret_info;

	pzheev_((char*) "V", (char*) "L", &n, 
		a, &row_offset, &col_offset, distrA.descriptor, 
		w, 
		z, &row_offset, &col_offset, distrZ.descriptor, 
		work, &lwork, rwork, &lrwork, &ret_info);

	lwork = work[0];

	free(work);
	free(rwork);

	work = (double*)malloc(lwork*2*sizeof(double));
	rwork  = (double*)malloc(lrwork*sizeof(double));

	pzheev_((char*) "V", (char*) "L", &n, 
		a, &row_offset, &col_offset, distrA.descriptor, 
		w, 
		z, &row_offset, &col_offset, distrZ.descriptor, 
		work, &lwork, rwork, &lrwork, &ret_info);

	free(work);
	free(rwork);

	for (int i=0; i<global_n_rows(); i++)
	{
		eigenvalues.push_back(w[i]);
	}

	Z.set_data(z);

	Z = ~Z;

	return Z;
}

Matrix exp (Matrix& A, complexd c)

{
	vector<complexd> eigenvalues;

	Matrix Out(A.global_n_rows(), A.global_n_cols());
	Matrix U(A.global_n_rows(), A.global_n_cols());
	Matrix U_c(A.global_n_cols(), A.global_n_rows());
	Matrix D(A.global_n_rows(), A.global_n_cols()); 

	Distribution distr = D.get_distribution();

	U = A.diagonalize(eigenvalues);
	U_c = U.herm_conj();

	for (int i=distr.row_offset()+D.n_rows; i>=distr.row_offset(); i--)
		for (int j=distr.col_offset()+D.n_cols; j>=distr.col_offset(); j--)
			if (i==j)
				D.set(i,j,exp(eigenvalues[i]*c));

	Out = U_c * D; 
	Out = Out * U;
			
	return Out;
}

Matrix commutator (Matrix& A, Matrix& B)
{
	return A*B-B*A;
}

Matrix diagonal_matrix(vector<complexd> values)
{
	Matrix out(values.size(),values.size());
	for (int i = 0; i < values.size(); i++)
		out.set(i,i,values[i]);
	return out;
}

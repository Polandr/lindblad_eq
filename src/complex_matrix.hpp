using namespace std;

//typedef complex<double> complexd;

class My_exception: public std::exception
{
	mutable char* errstr; 

	public:

	My_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~My_exception() throw()
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

// Parallel distribution-------------------------------------------------------------

// Distribution struct:

void Distribution::set_matrix_dims (int rows = 0, int cols = 0)
{
	matrix_global_rows = rows;
	matrix_global_cols = cols;
}

void Distribution::set_block_sizes (int row_size = R_BLOCK_SIZE, int col_size = C_BLOCK_SIZE)
{
	row_block_size = row_size;
	col_block_size = col_size;		
}

// Size and offset of local matrix:

int Distribution::local_row_num ()
{
	return numroc_(&matrix_global_rows,&row_block_size,
		&(ProcessorGrid::my_proc_row),&(ProcessorGrid::start_proc_row),&(ProcessorGrid::proc_row_num));
}

int Distribution::local_col_num ()
{
	return numroc_(&matrix_global_cols,&col_block_size,
		&(ProcessorGrid::my_proc_col),&(ProcessorGrid::start_proc_col),&(ProcessorGrid::proc_col_num));
}

int Distribution::row_offset ()
{
	return npreroc_(&matrix_global_rows,&row_block_size,
		&(ProcessorGrid::my_proc_row),&(ProcessorGrid::start_proc_row),&(ProcessorGrid::proc_row_num));
}

int Distribution::col_offset ()
{
	return npreroc_(&matrix_global_cols,&col_block_size,
		&(ProcessorGrid::my_proc_col),&(ProcessorGrid::start_proc_col),&(ProcessorGrid::proc_col_num));
}

void Distribution::operator = (const Distribution& other)
{
	row_block_size = other.row_block_size;
	col_block_size = other.col_block_size;
	matrix_global_rows = other.matrix_global_rows;
	matrix_global_cols = other.matrix_global_cols;
	for (int i = 0; i < DESC_LEN; i++)
		descriptor[i] = other.descriptor[i];
}


// ProcessorGrid struct:

void ProcessorGrid::init (int proc_rows = PRC_R, int proc_cols = PRC_C,
	int start_row = STRT_R, int start_col = STRT_C, int rt = ROOT)
{
	int grid_size;
	Cblacs_pinfo(&my_proc,&grid_size);
	if (proc_rows*proc_cols > grid_size ||
		my_proc < 0 || my_proc >= grid_size)
		throw My_exception("incorrect grid initialization");

	root = rt;
	start_proc_row = start_row;
	start_proc_col = start_col;

	Cblacs_get(0,0,&context);
	Cblacs_gridinit(&context, (char *)"Column", proc_rows, proc_cols);
	Cblacs_gridinfo(context,&proc_row_num,&proc_col_num,&my_proc_row,&my_proc_col);

	if (my_proc == rt)
		printf("Grid initialized!\n");

	initialized = true;
}

void ProcessorGrid::default_init ()
{
	int proc_num, grid_size;
	Cblacs_pinfo(&proc_num,&grid_size);
	int proc_dim = (int)sqrt(grid_size);
	if (proc_num == ROOT)
	{
		cout << flush;
		cout << "WARNING" << endl;
		cout << "Using default process grid initialization:" << endl;
		cout << "Root process has number: " << ROOT << endl;
		cout << "Process grid is " << proc_dim << "x" << proc_dim << 
		" starting from (" << STRT_R << "," << STRT_C << ") process." << endl;
		cout << flush;
	}
	ProcessorGrid::init(proc_dim, proc_dim, STRT_R, STRT_C, ROOT);
}

void ProcessorGrid::square_init (int start_row = STRT_R, int start_col = STRT_C, int rt = ROOT)
{
	int proc_num, grid_size;
	Cblacs_pinfo(&proc_num,&grid_size);
	int proc_dim = (int)sqrt(grid_size);
	ProcessorGrid::init(proc_dim, proc_dim, start_row, start_col, rt);
}

void ProcessorGrid::exit()
{
	Cblacs_gridexit(context);
	Cblacs_exit(0);
}


// Matrix class:

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
		throw My_exception("incorrect descriptor initialization");

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
		
void Matrix::set (int i, int j, complexd val)
{
	data[i*n_cols+j] = val;
}

void Matrix::set (int i, complexd val)
{
	data[i] = val;
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

void Matrix::get_row(double* array, int row) const
{
	for (int i = 0; i < n_cols; i++)
	{
		array[2*i] = data[row*n_cols+i].real();
		array[2*i+1] = data[row*n_cols+i].imag();
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
{
	for (int i = 0; i < n_cols; i++)
	{
		data[row*n_cols+i].real() = array[2*i];
		data[row*n_cols+i].imag() = array[2*i+1];
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
	if (rows <= 0 || cols <= 0)
		throw My_exception("invalid matrix size");
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
	if (rows <= 0 || cols <= 0)
		throw My_exception("invalid matrix size");
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

void Matrix::local_data_transpose ()
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

complexd& Matrix::operator () (int row, int col)
// Local indexing
{
	if (row >= n_rows || row < 0 || col >= n_cols || col < 0) 
		throw My_exception("out of bounds");
	else
		return data[row*n_cols+col];
}

// IT DOESN'T WORK :(
/*complexd& Matrix::operator [] (int row, int col)
// Global indexing
{
	if (row >= info.global_n_rows() || row < 0 || col >= info.global_n_cols() || col < 0) 
		throw My_exception("out of bounds");
	else
		if ((row >= info.row_offset()) &&
			(row < info.row_offset() + n_rows) &&
			(col >= info.col_offset()) &&
			(col < info.col_offset() + n_cols))
			return data[(row - info.row_offset()) + n_rows * (col - info.col_offset())];
}*/

// Arithmetics----------------------------------------------------------------------------

Matrix Matrix::operator * (Matrix& b) const
{
	Matrix c_matr(global_n_rows(), b.global_n_cols());

	int m = n_rows;
	int n = b.n_cols;
	int k = n_cols;

	Distribution distr = get_distribution();
	Distribution distrB = b.get_distribution();
	Distribution distrC = c_matr.get_distribution();

	int row_offset = 1;
	int col_offset = 1;

	double* a_data = get_data();
	double* b_data = b.get_data();
	double* c_data = c_matr.get_data();
	double* alpha = (double*)malloc(2*sizeof(double));
	double* beta = (double*)malloc(2*sizeof(double));
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

	int m = n_cols;
	int n = n_rows;

	double* alpha = (double*)malloc(2*sizeof(double));
	double* beta = (double*)malloc(2*sizeof(double));
	alpha[0] = 1.0;
	alpha[1] = 0.0;
	beta[0] = 1.0; 
	beta[1] = 0.0;


	int row_offset2 = distrA.row_offset()+1;
	int col_offset2 = distrA.col_offset()+1;


	double* a_data = get_data();
	double* c_data = c_matr.get_data();
	
	pzgeadd_((char*) "N", &m, &n, 
		alpha, a_data, &row_offset2, &col_offset2, distrA.descriptor, 
		beta, c_data, &row_offset2, &col_offset2, distrC.descriptor);

	c_matr.set_data(c_data);

	return c_matr;
}

Matrix Matrix::diagonalize (vector<complexd>& eigenvalues) const
{
	double *w = (double*)malloc(global_n_rows()*sizeof(double));
	Matrix Z(global_n_rows(), global_n_rows());

	Distribution distrA = get_distribution();
	Distribution distrZ = Z.get_distribution();

	int n = n_rows;

	double* a = get_data();
	double* z = Z.get_data();
	int row_offset = 1;
	int col_offset = 1;

	int lrwork = 2*n_rows + 2*n_rows-2;
	int lwork = -1;
	double *work = (double*)malloc(1*sizeof(double));
	double* rwork  = (double*)malloc(lwork*sizeof(double));
	int ret_info;

	pzheev_((char*) "V", (char*) "L", &n, 
		a, &row_offset, &col_offset, distrA.descriptor, 
		w, 
		z, &row_offset, &col_offset, distrZ.descriptor, 
		work, &lwork, rwork, &lrwork, &ret_info);

	for (int i=0; i<global_n_rows(); i++)
	{
		eigenvalues.push_back(w[i]);
		cout<<w[i]<<endl;
	}

	Z.set_data(z);
	
	cout<<Z<<endl;

	return Z;
}

Matrix exp (Matrix& A)
{
	vector<complexd> eigenvalues;
	Matrix Out(A.global_n_rows(), A.global_n_cols());
	Matrix U(A.global_n_rows(), A.global_n_cols());
	Matrix U_c(A.global_n_cols(), A.global_n_rows());
	Matrix D(A.global_n_rows(), A.global_n_cols());

	U = A.diagonalize(eigenvalues);
	U_c = ~U;
	
	for (int i=0; i<D.n_rows; i++)
		for (int j=0; j<D.n_cols; j++)
			if (i==j)
				D(i,j) = exp(eigenvalues[i]);

	Out = U_c * D;
	Out = Out * U;
	
	return Out;
}



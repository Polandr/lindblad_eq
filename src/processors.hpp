// ProcessorGrid and Distribution structures realization

using namespace std;

class Processors_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Processors_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Processors_exception() throw()
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

// ProcessorGrid struct:

void ProcessorGrid::init (int proc_rows = PRC_R, int proc_cols = PRC_C,
	int start_row = STRT_R, int start_col = STRT_C, int rt = ROOT)
{
	int grid_size;
	Cblacs_pinfo(&my_proc,&grid_size);
	if (proc_rows*proc_cols > grid_size ||
		my_proc < 0 || my_proc >= grid_size)
		throw Processors_exception("incorrect grid initialization");

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

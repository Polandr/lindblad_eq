
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <utility>

#include <mpi/mpi.h>
#include "scalapack.h"

typedef std::complex<double> complexd;

#define ROOT 0
#define STRT_R 0
#define STRT_C 0
#define PRC_R 1
#define PRC_C 1

struct ProcessorGrid
{
	static bool initialized;

	static int context;

	static int root;

	static int start_proc_row;
	static int start_proc_col;
	static int proc_row_num;
	static int proc_col_num;
	static int my_proc;
	static int my_proc_row;
	static int my_proc_col;

// Methods:

	static void init (int, int, int, int, int);
	static void default_init ();
	static void square_init (int, int, int);

	static int size () { return proc_row_num*proc_col_num; }
	static int is_root () { return my_proc == root; }
	static void barrier () { MPI_Barrier(MPI_COMM_WORLD); }

	static void exit ();
};

bool ProcessorGrid::initialized = false;
int ProcessorGrid::context = 0;
int ProcessorGrid::root = ROOT;
int ProcessorGrid::start_proc_row = STRT_R;
int ProcessorGrid::start_proc_col = STRT_C;
int ProcessorGrid::proc_row_num = PRC_R;
int ProcessorGrid::proc_col_num = PRC_C;
int ProcessorGrid::my_proc = 0;
int ProcessorGrid::my_proc_row = 0;
int ProcessorGrid::my_proc_col = 0;

#define DESC_LEN 9

#define R_BLOCK_SIZE 2
#define C_BLOCK_SIZE 2

struct Distribution
{
	// We can not to use representative of this class
	// cause this is purely static class:
	//ProcessorGrid proc_grid;

	int descriptor[DESC_LEN];

	int row_block_size;
	int col_block_size;

	int matrix_global_rows;
	int matrix_global_cols;

// Methods:

	void set_matrix_dims (int, int);
	void set_block_sizes (int, int);

	int proc_grid_size () const { return ProcessorGrid::size(); }

	int local_row_num ();
	int local_col_num ();

	int row_offset ();
	int col_offset ();

	void operator = (const Distribution&);
};


class Matrix
{
	complexd *data;
	Distribution info;

	void create ();
	void destroy ();

	bool in_block (int, int);

	int* gather_info ();

public:
	int n_rows, n_cols;

	int get_n_rows () const { return n_rows; };
	int get_n_cols () const { return n_cols; };
	int global_n_rows () const { return info.matrix_global_rows; }
	int global_n_cols () const { return info.matrix_global_cols; }
	const Distribution& get_distribution () const { return info; }
	const int* get_descriptor() const { return info.descriptor; }
	double* get_data () const;
	void get_row (double*, int) const;

	void set (int, int, complexd);
	void set (int, complexd);
	complexd& operator () (int, int);
	void set_data (double*);
	void set_row (double*, int);

	Matrix ();
	Matrix (int, int, int, int);
	void init (int, int, int, int);
	Matrix (const Matrix&);
	~Matrix ();

	bool is_square() const { return (global_n_rows() == global_n_cols()); }
	void init_distribution (int, int, int, int);

	void local_data_transpose();

	Matrix& operator = (const Matrix&);
	Matrix operator * (Matrix&) const;
	Matrix operator ~ () const;
	Matrix diagonalize(std::vector<complexd>&) const;

	int readf (const char*, int, int);
	void writef (int, const char*);
	void print_on_condition (std::ofstream&, bool (*cond)(int, int));
	void generate (complexd (*func)(int, int));

	void operator >> (std::ostream&);
	void operator << (std::istream&);
};

std::ostream& operator << (std::ostream&, Matrix&);
std::istream& operator >> (std::istream&, Matrix&);

#include "complex_matrix.hpp"
#include "complex_matrix_io.hpp"

#ifndef __PROCESSORS_H__
#define __PROCESSORS_H__

#define ROOT 0
#define STRT_R 0
#define STRT_C 0
#define PRC_R 1
#define PRC_C 1

// Purely static class which describes grid of processors:

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

#include "processors.hpp"

#endif

// __PROCESSORS_H__

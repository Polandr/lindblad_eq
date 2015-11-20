#ifndef __MATRIX_H__
#define __MATRIX_H__

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
#include "processors.h"

typedef std::complex<double> complexd;



class Matrix
{
	complexd *data;
	Distribution info;

	void create ();
	void destroy ();

	bool in_block (int, int) const;

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
	void set_data (double*);
	void set_row (double*, int);

	const complexd operator () (int, int) const;

	Matrix ();
	Matrix (int, int, int, int);
	void init (int, int, int, int);
	Matrix (const Matrix&);
	~Matrix ();

	bool is_square() const { return (global_n_rows() == global_n_cols()); }
	void init_distribution (int, int, int, int);

	void in_place_transposition();

	Matrix& operator = (const Matrix&);
	Matrix& operator *= (complexd);
	Matrix operator * (complexd) const;
	Matrix operator * (Matrix&) const;
	Matrix operator ~ () const;
	Matrix diagonalize(std::vector<complexd>&) const;

	int readf (const char*, int, int);
	void writef (int, const char*);
	void print_on_condition (std::ostream&, bool (*cond)(int, int));
	void print_diagonal_abs(FILE*);
	void generate (complexd (*func)(int, int));

	void operator >> (std::ostream&);
	void operator << (std::istream&);
	void stream_output (std::ostream&);
	void stream_input (std::istream&);
};

std::ostream& operator << (std::ostream&, Matrix&);
std::istream& operator >> (std::istream&, Matrix&);

Matrix exp (Matrix);

#include "complex_matrix.hpp"
#include "complex_matrix_io.hpp"

#endif

// __MATRIX_H__

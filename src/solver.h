#include <iostream>
#include <cstdio>
#include <utility>
#include <fstream>


#ifndef __SOLVER_H__
#define __SOLVER_H__

class Solver
{
	Matrix H;
	Matrix R0;
	double dT;
	int step_num;

public:

	Matrix& get_H() { return H; }
	Matrix& get_R0() { return R0; }
	double& get_dT() { return dT; }
	int& get_step_num() { return step_num; }

	void init_H(const char*);
	void init_R0(const char*);
	void init_dT(double);
	void init_step_num(int);
	void init_system();

	void solve(const char*);
};

std::ostream& operator << (std::ostream&, const Solver&);
std::istream& operator >> (std::istream&, const Solver&);


#include "solver.hpp"

#endif 

/*  __SOLVER_H__  */
#ifndef __QUANTUM_SYSTEM_H__
#define __QUANTUM_SYSTEM_H__

#include <vector>
#include <cstdio>
#include <iostream>

#include "complex_matrix.h"

#ifndef EPS
#define EPS 0.000001
#endif

#include "exceptions.hpp"

using namespace std;

// Normalization of hermitian matrix:

void herm_normalization (Matrix& H);

// Hamiltonian:

Matrix hamiltonian (int N, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w, vector<int>& base_states);
// N - dimension of system
// s - maximum stock level
// E_min, E_max - minimum and maximum energy levels
// a - probabilities between atoms
// w - probabilities on atoms
// base_states - base states of the system

// Density matrix:

Matrix density_matrix (vector<complexd> state);
// Make density matrix for quantum state
Matrix density_matrix (int N, int i);
// Make density matrix with single 1 on main diagonal on position number <pos> and other elemants set as 0
Matrix density_matrix (vector<double> qbit_probs, vector<double> stock_probs, vector<int> base_states, vector<int> state_nums);
// Make density matrix with probabilities of {qbit = 1} event as <qbit_probs> and {stock = i} event as <stock_probs>
// If <qbit_probs> designate set of possible states (at least 1): states which were set during hamiltonian construction
// then set of possible states are chosen with certain ratio
// If <qbit_probs> designate set of impossible states
// then all states are set as equiprobable

double summarize_amplitudes_in_stock_block (Matrix& density, int block_num, vector<int> base_states, vector<int> state_nums);
// Summarize amplitudes of block corresponded to <block_num> stock state in density matrix <density>

// Lindblad:

class Lindblad_part
{
	complexd stock;
	std::vector<complexd> l;

	Matrix L0;
	vector<complexd> L;
	int L_num () { return l.size(); }
	int L_dim () { return L.size()/L_num(); }

public:

	bool active;

	Lindblad_part () { active = false; }
	void init (complexd, vector<complexd>, vector<int> base_states, vector<int> state_nums);
	void destroy ();

	Matrix operator () (Matrix&);

	void print_matrices (ostream&);
};


// Implementation:-------------------------------------------------------------------------------------------------------

using namespace std;

// Service functions:

int combination_num(int k, int n)
{
	int out = 1;
	for (int i = k+1; i <= n; i++)
		out *= i;
	for (int i = 1; i <= n-k; i++)
		out /= i;
	return out;
}

void print_ketbra(int state, int N)
{
	printf("|");
	for (int i = 0, mask = 1 << N; i < N; i++)
	{
		mask = mask >> 1;
		if ((mask & state) != 0)
			printf("1");
		else
			printf("0");			
	}
	printf(">");
}

void print_ketbra_stream(int state, int N, ostream& out = cout)
{
	out << '|';
	for (int i = 0, mask = 1 << N; i < N; i++)
	{
		mask = mask >> 1;
		if ((mask & state) != 0)
			out << '1';
		else
			out << '0';			
	}
	out << '>';
}

int unit_num(int state, int N)
{
	int out = 0;
	for (int i = 0, mask = 1; i < N; i++)
	{
		if ((mask & state) != 0)
			out++;
		mask = mask << 1;
	}
	return out;
}

// For hamiltonian:

vector<int> collect_base_states (int E_low, int E_high, int N)
{
	vector<int> cur_base_states;
	int pow_N = pow(2,N);
	for (int i = 0; i < pow_N; i++)
	{
		int E_lvl = unit_num(i,N);
		if (E_lvl >= E_low && E_lvl <= E_high)
			cur_base_states.push_back(i);
	}

	return cur_base_states;
}

int simple_transition(int state_1, int state_2, int N)
{
	int out = -1;
	int state_sum = state_1 ^ state_2;

	if (unit_num(state_sum,N) != 2)
		return -1;
	
	for (int i = 0, mask = 1; i < N; i++)
	{
		if ((mask & state_sum) != 0)
		{
			if ((mask & state_1) != 0)
			{
				mask = mask << 1;
				if (((mask & state_sum) != 0) && ((mask & state_2) != 0))
					return i;
				else
					return -1;
			}
			if ((mask & state_2) != 0)
			{
				mask = mask << 1;
				if (((mask & state_sum) != 0) && ((mask & state_1) != 0))
					return i;
				else
					return -1;
			}
		}
		mask = mask << 1;
	}
}

complexd hamiltonian_element(int row, int col, int N, vector<complexd> a, vector<complexd> w, vector<int> states)
{
	if (row == col)
	{
		complexd out(0,0);
		for (int i = 0, mask = 1; i < N; i++)
		{
			if ((mask & states[row]) != 0)
			{
				out += w[i];
			}
			mask = mask << 1;
		}
		return out;
	}
	else
	{
		int pos = simple_transition(states[row],states[col],N);
		if (pos >= 0)
		{
			return a[pos];
		}
		else
		{
			return complexd(0,0);
		}
	}
}

// For density matrix:

double state_probability(int state, vector<double> qbit_probs)
{
	double out = 1;
	for (int i = 0, mask = 1; i < qbit_probs.size(); i++)
	{
		if ((mask & state) != 0)
			out *= qbit_probs[i];
		else
			out *= (1 - qbit_probs[i]);
		mask = mask << 1;
	}
	return out;
}

// For Lindblad:

int makeMask(int i)
{
	int mask = pow(2,i);
	return mask;
}

int getPrevState(int state, bool& found)
{
	if ((state & 1) == 0)
	{
		found = true;
		return (state | 1);
	}
	else
	{
		found = false;
		return -1;
	}
}

Matrix createStockMatrix(vector<int> base_states, vector<int> state_nums)
{
	Matrix stockMatrix(base_states.size(),base_states.size());

	if (state_nums.size() == 1)
		return stockMatrix;

	int ofs = state_nums[0];
	for (int i = 1; i < state_nums.size(); ++i)
	{
		for (int row = ofs; row < ofs + state_nums[i]; ++row)
		{
			bool found;
			int prevState = getPrevState(base_states[row], found);
			if (found)
			{
				int col = -1;
				for (int k = ofs-state_nums[i-1]; k < ofs; ++k)
					if (base_states[k] == prevState)
						col = k;
				if (col != -1)
				{
					stockMatrix.set(row,col,1);
					//stockMatrix.set(col,row,1);
				}
			}
		}
		ofs += state_nums[i];
	}

	return stockMatrix;
}

void addDephaseMatrix(vector<complexd>& L, int pos, vector<int> base_states)
{
	int mask = makeMask(pos);

	for (int i = 0; i < base_states.size(); ++i)
	{
		if ((mask & base_states[i]) == mask)
			L.push_back(1);
		else
			L.push_back(0);
	}
}

// Normalization of hermitian matrix:------------------------------------------------------------------------------------

void herm_normalization (Matrix& H)
{
	H = (H + H.herm_conj()) * 0.5 / trace(H);
}

// Hamiltonian constructing:---------------------------------------------------------------------------------------------

Matrix hamiltonian (int N, int s, int E_min, int E_max, 
	vector<complexd> a, vector<complexd> w, vector<int>& base_states, vector<int>& state_nums)
// N - dimension of system
// s - maximum stock level
// E_min, E_max - minimum and maximum energy levels
// a - probabilities between atoms
// w - probabilities on atoms
// base_states - base states of the system
{
	Matrix H;

	s = min(s,E_max);

	for (int i = 0; i <= s; i++)
	{
		int low = max(0,E_min-i);
		int high = max(0,E_max-i);
		vector<int> cur_base_states = collect_base_states(low,high,N);
		state_nums.push_back(cur_base_states.size());
		base_states.insert(base_states.end(),cur_base_states.begin(),cur_base_states.end());
	}

	H.init(base_states.size(),base_states.size());

	for (int block_num = 0, ofs = 0; block_num < state_nums.size(); block_num++)
	{
		for (int i = ofs; i < ofs + state_nums[block_num]; i++)
			for (int j = i; j < ofs + state_nums[block_num]; j++)
			{
				complexd val = hamiltonian_element(i,j,N,a,w,base_states);
				if (val != complexd(0,0))
				{
					H.set(i,j,val);
					if (i != j)
						H.set(j,i,conj(val));
				}
			}
		ofs += state_nums[block_num];
	}

	return H;	
}

// Density matrix constructing:------------------------------------------------------------------------------------------

Matrix density_matrix (vector<complexd> state)
{
	Matrix out(state.size(),state.size());
	for (int i = 0; i < out.global_n_rows(); i++)
		for (int j = 0; j < out.global_n_cols(); j++)
			out.set(i,j,state[i]*std::conj(state[j]));
	return out;
}

Matrix density_matrix (int N, int i)
{
	vector<complexd> state(N,0);
	state[i] = 1;
	return density_matrix(state);
}

Matrix density_matrix (vector<double> qbit_probs, vector<double> stock_probs, vector<int> base_states, vector<int> state_nums)
{
	vector<complexd> diagonal(base_states.size(),0);

	double probs_sum = 0;
	for (int block_num = 0, ofs = 0; block_num < state_nums.size(); block_num++)
	{
		for (int i = ofs; i < ofs + state_nums[block_num]; i++)
		{
			diagonal[i] = state_probability(base_states[i],qbit_probs) * stock_probs[block_num];
			probs_sum += diagonal[i].real();
		}
		ofs += state_nums[block_num];
	}

	if (fabs(probs_sum - 1) >= EPS)
	// If sum of probabilities is not equal to 1 then normalize probabilities
	{
		for (int i = 0; i < diagonal.size(); ++i)
			if (probs_sum != 0)
				diagonal[i] /= probs_sum;
			else
				diagonal[i] = 1.0/diagonal.size();
	}

	return diagonal_matrix(diagonal);
}

double summarize_amplitudes_in_stock_block (Matrix& density, int block_num, vector<int> base_states, vector<int> state_nums)
{
	if (block_num >= state_nums.size())
		throw  Solver_exception("invalid number of block to summarize in density matrix");

	int ofs = 0;
	for (int i = 0; i < block_num; ++i)
		ofs += state_nums[i];

	double sum = 0;
	for (int i = ofs; i < ofs+state_nums[block_num]; ++i)
		sum += abs(density(i,i));

	return sum;
}

// Lindblad part:--------------------------------------------------------------------------------------------------------

void Lindblad_part::init (complexd out, vector<complexd> ls, vector<int> base_states, vector<int> state_nums)
{
	if (active)
		destroy();

	active = true;
	stock = out;

	for (int i = 0; i < ls.size(); ++i)
		l.push_back(ls[i]);

	L0 = createStockMatrix(base_states, state_nums);

	for (int i = 0; i < L_num(); ++i)
		addDephaseMatrix(L,i,base_states);
}

void Lindblad_part::destroy ()
{
	active = false;
	l.clear();
	l.resize(0);
	L.clear();
	L.resize(0);
	//L0.~Matrix();
	L0.init();
}

Matrix Lindblad_part::operator () (Matrix& R)
{
	Matrix Out;

	// Stock matrix
	Out = stock*(L0*R*(L0.herm_conj()) - 0.5*((L0.herm_conj())*L0*R + R*(L0.herm_conj())*L0));

	// Dephase matrices
	int dim = L_dim();
	for (int i = 0; i < L_num(); i++)
	{
		vector<complexd> diagonal(L.begin() + i*dim, L.begin() + (i+1)*dim);
		Matrix Li = diagonal_matrix(diagonal);
		Out += l[i]*(Li*R*(Li.herm_conj()) - 0.5*((Li.herm_conj())*Li*R + R*(Li.herm_conj())*Li));
	}

	return Out;
}

void Lindblad_part::print_matrices (ostream& out = cout)
{
	// Stock matrix
	if (ProcessorGrid::is_root())
	{
		out << "L0:\n";
	}
	out << L0;

	// Dephase matrices
	int dim = L_dim();
	for (int i = 0, ofs = 0; i < L_num(); ++i)
	{
		if (ProcessorGrid::is_root())
		{
			out << "L" << (i+1) << ":\n";
			for (int j = 0; j < dim; ++j)
			{
				out << L[ofs+j] << ' ';
			}
			out << endl;
		}
		/*vector<complexd> diagonal(L.begin() + ofs, L.begin() + ofs+dim);
		Matrix Li = diagonal_matrix(diagonal);
		out << Li;*/
		ofs += dim;
	}

	out << flush;
}

#endif

// __QUANTUM_SYSTEM_H__ 

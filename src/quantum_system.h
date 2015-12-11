#ifndef __QUANTUM_SYSTEM_H__
#define __QUANTUM_SYSTEM_H__

#include <vector>
#include <cstdio>
#include <iostream>

#include "complex_matrix.h"

using namespace std;

// Hamiltonian:

Matrix hamiltonian (int N, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w, vector<int>& base_states);
// N - dimension of system
// s - maximum stock level
// E_min, E_max - minimum and maximum energy levels
// a - probabilities between atoms
// w - probabilities on atoms
// base_states - base states of the system

// Density matrix:

Matrix density_matrix(vector<complexd> state);
Matrix density_matrix(int N, int i);

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
					stockMatrix.set(col,row,1);
				}
			}
		}
		ofs += state_nums[i];
	}

	return stockMatrix;
}

void addDephaseMatrix(vector<complexd>& L, int pos, complexd coeff, vector<int> base_states)
{
	int mask = makeMask(pos);

	for (int i = 0; i < base_states.size(); ++i)
	{
		if ((mask & base_states[i]) == mask)
			L.push_back(coeff);
		else
			L.push_back(0);
	}
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

Matrix density_matrix(vector<complexd> state)
{
	Matrix out(state.size(),state.size());
	for (int i = 0; i < out.global_n_rows(); i++)
		for (int j = 0; j < out.global_n_cols(); j++)
			out.set(i,j,state[i]*std::conj(state[j]));
	return out;
}

Matrix density_matrix(int N, int i)
{
	if (i >= N)
		throw  Matrix_exception("invalid density matrix initialization");
	vector<complexd> state(N,0);
	state[i] = 1;
	return density_matrix(state);
}

// Lindblad part:--------------------------------------------------------------------------------------------------------

void Lindblad_part::init (complexd out, vector<complexd> ls, vector<int> base_states, vector<int> state_nums)
{
	stock = out;

	for (int i = 0; i < ls.size(); ++i)
		l.push_back(ls[i]);

	L0 = createStockMatrix(base_states, state_nums);

	for (int i = 0; i < L_num(); ++i)
		addDephaseMatrix(L,i,l[i],base_states);
}

Matrix Lindblad_part::operator () (Matrix& R)
{
	Matrix Out;
	complexd imag_unit(0,1);

	// Stock matrix
	Out = stock*(L0*R*(L0.herm_conj()) - (1/2)*((L0.herm_conj())*L0*R + R*(L0.herm_conj())*L0));

	// Dephase matrices
	int dim = L_dim();
	for (int i = 0; i < L_num(); i++)
	{
		vector<complexd> diagonal(L.begin() + i*dim, L.begin() + (i+1)*dim);
		Matrix Li = diagonal_matrix(diagonal);
		Out += l[i]*(Li*R*(Li.herm_conj()) - (1/2)*((Li.herm_conj())*Li*R + R*(Li.herm_conj())*Li));
	}
	
	Out = Out*imag_unit;

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

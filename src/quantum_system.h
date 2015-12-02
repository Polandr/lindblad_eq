#ifndef __QUANTUM_SYSTEM_H__
#define __QUANTUM_SYSTEM_H__

#include <vector>

#include "complex_matrix.h"

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
	int stock;
	std::vector<complexd> ls;

public:

	bool active;

	Lindblad_part() { active = false; }
	void init (int, std::vector<complexd>);

	Matrix operator () (Matrix&, vector<int>);
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

void collect_base_states (int E_low, int E_high, int N, vector<int>& base_states)
{
	int mask = 1;
	int full_energy_state = 0;
	for (int i = 0; i < N; i++)
	{
		full_energy_state = full_energy_state | mask;
		mask = mask << 1;
	}
	vector<int> cur_base_states;
	int pow_N = pow(2,N);
	for (int i = 0; i < pow_N; i++)
	{
		int E_lvl = unit_num(i,N);
		if (E_lvl >= E_low && E_lvl <= E_high)
			cur_base_states.push_back(i);
	}
	base_states.insert(base_states.end(),cur_base_states.begin(),cur_base_states.end());
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

// Hamiltonian constructing:---------------------------------------------------------------------------------------------

Matrix hamiltonian (int N, int s, int E_min, int E_max, vector<complexd> a, vector<complexd> w, vector<int>& base_states)
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
		int low = max(0,E_min-s);
		int high = max(0,E_max-s);
		collect_base_states(low,high,N,base_states);
	}

	H.init(base_states.size(),base_states.size());
	for (int i = 0; i < base_states.size(); i++)
		for (int j = i; j < base_states.size(); j++)
		{
			complexd val = hamiltonian_element(i,j,N,a,w,base_states);
			if (val != complexd(0,0))
			{
				H.set(i,j,val);
				if (i != j)
					H.set(j,i,conj(val));
			}
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

int makeMask(int i)
{
	int x = i-1;
	int digit = pow(2,x);
	return digit;
}

Matrix createStockMatrix()
{
	;
}

Matrix createDiffaseMatrix(int i, vector<int> base_states)
{
	Matrix diffaseMatr(base_states.size(),base_states.size());
	i++;
	int mask = makeMask(i);

	for (int j=0; j<base_states.size(); j++)
	{
		if ((mask & base_states[j]) == mask)
			diffaseMatr.set(j,j,1);
	}
	return diffaseMatr;
}

// Lindblad part:--------------------------------------------------------------------------------------------------------

void Lindblad_part::init (int out, std::vector<complexd> v1)
{

	stock = out;
	for (int i = 0; i < v1.size(); ++i)
		ls.push_back(v1[i]);
}

Matrix Lindblad_part::operator () (Matrix& R, vector<int> base_states)
{
	vector<Matrix> Li;
	Matrix Out(base_states.size(),base_states.size());
	Matrix Out2, Out3, Li_conj, Li_conj_Li;
	complexd imag_unit(0,1.0);
	int N = ls.size();

//zero step - stock
if (stock != 0)
{
	Li.push_back(createStockMatrix());
	Li_conj = ~Li[0];
	Li_conj_Li = Li_conj*Li[0];

	Out = Li[0]*R;

	Out = Out*Li_conj;

	Out2 = Li_conj_Li*R;
	Out3 = R*Li_conj_Li;

	Out2 += Out3;
	Out2 = Out2*(1/2);

	Out -= Out2;
	Out = Out*stock;
}
// end stock

	for (int i=0; i<N; i++)
	{
		Li.push_back(createDiffaseMatrix(i, base_states));
		Li_conj = ~Li[i];
		Li_conj_Li = Li_conj*Li[i];
		
		Out += Li[i]*R;

		Out = Out*Li_conj;

		Out2 = Li_conj_Li*R;
		Out3 = R*Li_conj_Li;

		Out2 += Out3;
		Out2 = Out2 *(1/2);

		Out -= Out2;
		Out = Out*ls[i];
	}
	
	Out = Out*imag_unit;

	return Out;
}

#endif

// __QUANTUM_SYSTEM_H__ 

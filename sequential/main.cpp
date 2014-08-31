#include <iostream>
#include <vector>
#include <list>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "Magic_State.h"

//Magic square:
//Sum = n(1+n^2)/2

using namespace std;

int n = 0;
int Sum = 0;
bool stopSearch = false;
clock_t t1 = 0, t2 = 0;

State state;

void print_state(const State& s);
void solve_state(State& s);
void get_possible_transforms(State& s, list<Transform>& transfs);
bool is_correct_state(const State& s);

int get_column_sum(const State& s, short col);
int get_row_sum(const State& s, short row);
int get_mdiag_sum(const State& s);
int get_sdiag_sum(const State& s);

int main(int argc, char *argv[])
{
	n = 6;

	Sum = n*(1+n*n)/2;
	state.init(n);
	
	state.put(0, 0, 1);
	state.put(0, 1, 35);
	state.put(0, 2, 34);
	state.put(0, 3, 3);
	state.put(0, 4, 32);
	state.put(0, 5, 6);

	state.put(1, 0, 30);
	state.put(1, 1, 8);
	state.put(1, 2, 28);
	state.put(1, 3, 27);
	state.put(1, 4, 11);
	state.put(1, 5, 7);
	
	state.put(2, 0, 24);
	//state.put(2, 1, 23);
	//state.put(2, 2, 15);
	//state.put(2, 3, 16);

	t1 = clock();
	solve_state(state);
	t2 = clock();

	print_state(state);
	printf("Time elapsed: %.2f s \n", (float)(t2-t1)/CLOCKS_PER_SEC);
	return 0;
}

bool is_correct_state(const State& s)
{
	short vert_sum = 0;		short vert_k = 0;
	short hor_sum = 0;		short hor_k = 0;
	short main_diag_sum = 0;	short main_k = 0;
	short sec_diag_sum = 0;		short sec_k = 0;

	for(int i=0; i<n; i++)
	{
		if( s.m[i][i] != 0 )
		{
			main_diag_sum += s.m[i][i];
			main_k++;
		}
		if( s.m[i][n-i-1] != 0 )
		{
			sec_diag_sum += s.m[i][n-i-1];
			sec_k++;
		}
		vert_sum = 0; vert_k = 0;
		hor_sum = 0; hor_k = 0;
		for(int j=0; j<n; j++)
		{
			if( s.m[i][j] != 0 )
			{
				hor_sum += s.m[i][j];
				hor_k++;
			}
			if( s.m[j][i] != 0 )
			{
				vert_sum += s.m[j][i];
				vert_k++;
			}
		}
		if( hor_sum != Sum && hor_k == n ) return false;
		if( vert_sum != Sum && vert_k == n ) return false;
	}
	if( main_diag_sum != Sum && main_k == n ) return false;
	if( sec_diag_sum != Sum && sec_k == n ) return false;
    
	return true;
}

void get_possible_transforms(State& state, list<Transform>& transfs)
{
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			if( state.m[i][j]==0 )
			{
				for(int k=1; k<=n*n; k++)
				{
					if( !state.used[k] )
					{
                        Transform s;
                        s.t.push_back(OneStepTransform(i, j, k));
						state.put(i, j, k);
						if( n-i==2 )
						{
							short lastVertValue = Sum - get_column_sum(state, j);
							if( lastVertValue <= n*n && lastVertValue >= 1 && state.used[lastVertValue]==0 )
							{
								state.put(i+1, j, lastVertValue);
								s.t.push_back(OneStepTransform(i+1, j, lastVertValue));
							}
							else
							{
								goto rollback;
							}
						}
						if( n-j==2 )
						{
							short lastVertValue = Sum - get_row_sum(state, i);
							if( lastVertValue <= n*n && lastVertValue >= 1 && state.used[lastVertValue]==0 )
							{
								state.put(i, j+1, lastVertValue);
								s.t.push_back(OneStepTransform(i, j+1, lastVertValue));
							}
							else
							{
								goto rollback;
							}
						}
						if( is_correct_state(state) )
						{
							transfs.push_back(s);
						}
                        
						rollback:
                        state.revert(s);
					}
				}
				goto exit;
			}
		}
	}

exit:
	return;
}

void solve_state(State& state)
{
	if( !stopSearch )
	{
		if( state.isFinal == n*n )
		{
			stopSearch = true;
		}
		else
		{
			list<Transform> transforms;
			get_possible_transforms(state, transforms);
			for(list<Transform>::iterator iTrans = transforms.begin(); iTrans != transforms.end(); iTrans++)
			{
                state.apply(*iTrans);

				solve_state(state);
                
				if( stopSearch ) { return; }
				else
                {
                    state.revert(*iTrans);
				}
			}
		}
	}
}

void print_state(const State& s)
{
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			cout<<s.m[i][j]<<" ";
		}
		cout<<endl;
	}
}

int get_column_sum(const State& state, short col)
{
	short sum = 0;
	for(int i=0; i<n; i++)
	{
		sum += state.m[i][col];
	}
	return sum;
}
int get_row_sum(const State& state, short row)
{
	short sum = 0;
	for(int i=0; i<n; i++)
	{
		sum += state.m[row][i];
	}
	return sum;
}
int get_mdiag_sum(const State& state)
{
	short sum = 0;
	for(int i=0; i<n; i++)
	{
		sum += state.m[i][i];
	}
	return sum;
}
int get_sdiag_sum(const State& state)
{
	short sum = 0;
	for(int i=0; i<n; i++)
	{
		sum += state.m[i][n-i-1];
	}
	return sum;
}

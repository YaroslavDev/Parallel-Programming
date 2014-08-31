#include "mpi.h"
#include <iostream>
#include <vector>
#include <list>
#include <stdlib.h>
#include <time.h>
#include "Magic_State.h"

//Magic square:
//Sum = n(1+n^2)/2

using namespace std;

int n = 0;
int Sum = 0;
bool stopSearch = false;

int msg_length = 0;
int my_rank = 0;
int p = 0;
int tag = 1;
MPI_Status status;
char *buffer = NULL;
int position = 0;
int flag = 0;

clock_t t1 = 0, t2 = 0;

State g_state;

//Master functions
void generate_jobs(State& s, list<Transform>& jobs);
void send_job(State& s, int recipient);
void receive_response(State& s, int sender);
void shutdown(State& s);
//Slave functions
void send_response(State& s);
void receive_job(State& s);
void solve_state(State& s);

//Auxiliary functions
void initialize_state(State& s);
void print_state(const State& s);
void get_possible_transforms(State& s, list<Transform>& transfs);
bool is_correct_state(const State& s);

int get_column_sum(const State& s, short col);
int get_row_sum(const State& s, short row);
int get_mdiag_sum(const State& s);
int get_sdiag_sum(const State& s);

void initialize_state(State& s)
{
    s.clear();
    s.put(0, 0, 1);
    s.put(0, 1, 35);
    s.put(0, 2, 34);
    s.put(0, 3, 3);
    s.put(0, 4, 32);
    s.put(0, 5, 6);

    s.put(1, 0, 30);
    s.put(1, 1, 8);
    s.put(1, 2, 28);
    s.put(1, 3, 27);
    s.put(1, 4, 11);
    s.put(1, 5, 7);

    s.put(2, 0, 24);
    s.put(2, 1, 23);
    s.put(2, 2, 15);
    //s.put(2, 3, 16);
    //s.put(2, 4, 14);
    //s.put(2, 5, 19);
}

int main(int argc, char **argv)
{
    n = 6;

    msg_length = 4*n*n;
    buffer = new char[msg_length];
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    
    Sum = n*(1+n*n)/2;
    g_state.init(n);
    initialize_state(g_state);

    if( my_rank == 0 )  //master branch
    {
    	t1 = clock();
        list<Transform> jobs;
        generate_jobs(g_state, jobs);
        printf("Master: Total # of processors - %d\n", p);
        printf("Master: Total # of jobs - %lu\n", jobs.size());
        printf("Master: N = %d\n", n);
        printf("Master: Msg Length = %d\n", msg_length);
        for(int slave=1; slave < p; slave++)
        {
            //Send jobs to all slaves
            Transform job = *jobs.begin();
            jobs.pop_front();
            
            g_state.apply(job);
            
            send_job(g_state, slave);
            
            g_state.revert(job);
        }
        int slave = 1;
        while( !jobs.empty() && !stopSearch )
        {
            flag = 0;
            MPI_Iprobe(slave, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if( flag )
            {
                receive_response(g_state, slave);
                if( !stopSearch )
                {
		    //initialize_state(g_state);
                    Transform job = *jobs.begin();
                    jobs.pop_front();
                
                    g_state.apply(job);
                
                    send_job(g_state, slave);
                
                    g_state.revert(job);
                }
            }
            slave++;
            if( slave >= p ) slave=1;
        }
        print_state(g_state); //print solution to screen
	t2 = clock();
        shutdown(g_state);
        printf("Master finished work in %.2f s \n", (float)(t2-t1)/CLOCKS_PER_SEC);
    }
    else                //slave branch
    {
        while( !stopSearch )
        {
            receive_job(g_state);

	    solve_state(g_state);
            if( g_state.isFinal != n*n && !stopSearch )
            {
                initialize_state(g_state);
            }
            send_response(g_state);
        }
        
        //solve_state(g_state);
        printf("Slave #%d finished work\n", my_rank);
    }
    
    MPI_Finalize();
    
    delete[] buffer;
    buffer = NULL;

	return 0;
}

bool is_correct_state(const State& s)
{
	short vert_sum = 0;			short vert_k = 0;
	short hor_sum = 0;			short hor_k = 0;
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

void solve_state(State& s)
{
	if( !stopSearch )
	{
		if( s.isFinal == n*n )
		{
			stopSearch = true;
		}
		else
		{
			list<Transform> transforms;
			get_possible_transforms(s, transforms);
			for(list<Transform>::iterator iTrans = transforms.begin(); iTrans != transforms.end(); iTrans++)
			{
                s.apply(*iTrans);

				solve_state(s);
                
				if( stopSearch ) { return; }
				else
                {
                    s.revert(*iTrans);
				}
			}
		}
	}
}

//Master functions
void generate_jobs(State& state, list<Transform>& jobs)
{
    get_possible_transforms(state, jobs);
    int emptyCells = n - state.isFinal;
    int nIters = emptyCells*emptyCells;
    for(int k=0; k<nIters; k++)
    {
        Transform job = *jobs.begin();
        jobs.pop_front();
        state.apply(job);
        list<Transform> new_jobs;
        get_possible_transforms(state, new_jobs);
        for(list<Transform>::iterator it=new_jobs.begin(); it != new_jobs.end(); it++)
        {
            Transform new_job = job;
            new_job.merge(*it);
            jobs.push_back(new_job);
        }
        state.revert(job);
    }
}

void send_job(State& s, int recipient)
{
    position = 0;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            short v = s.m[i][j];
            MPI_Pack(&v, 1, MPI_SHORT, buffer, msg_length, &position, MPI_COMM_WORLD);
        }
    }
    MPI_Send(buffer, position, MPI_PACKED, recipient, 1, MPI_COMM_WORLD);
}
void receive_response(State& s, int sender)
{
    position = 0;
    MPI_Recv(buffer, msg_length, MPI_PACKED, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            MPI_Unpack(buffer, msg_length, &position, &s.m[i][j], 1, MPI_SHORT, MPI_COMM_WORLD);
        }
    }
    s.analyse();
    stopSearch = (s.isFinal == n*n);
    //stopSearch = (s.m[0][0] != 0);
}
void shutdown(State& s)
{
    s.clear();  //send clear state - means no more jobs, shutdown
    position = 0;
    for(int slave=1; slave < p; slave++)
    {
        send_job(s, slave);
    }
}
//Slave functions
void send_response(State& s)
{
    position = 0;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            short v = s.m[i][j];
            MPI_Pack(&v, 1, MPI_SHORT, buffer, msg_length, &position, MPI_COMM_WORLD);
        }
    }
    MPI_Send(buffer, position, MPI_PACKED, 0, 1, MPI_COMM_WORLD);
}
void receive_job(State& s)
{
    position = 0;
    MPI_Recv(buffer, msg_length, MPI_PACKED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            MPI_Unpack(buffer, msg_length, &position, &s.m[i][j], 1, MPI_SHORT, MPI_COMM_WORLD);
        }
    }
    s.analyse();
    if( s.m[0][0] == 0 )
    {
        stopSearch = true;
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

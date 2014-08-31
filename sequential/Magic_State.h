//
//  Magic_State.h
//  Magic_Square
//
//  Created by Ciupin Iaroslav on 10/11/13.
//  Copyright (c) 2013 Ciupin Iaroslav. All rights reserved.
//

#ifndef Magic_Square_Magic_State_h
#define Magic_Square_Magic_State_h

#include <vector>
#include <list>

using namespace std;

struct OneStepTransform
{
	OneStepTransform(int _i = 0, int _j = 0, int _val = 0)
    : i(_i), j(_j), val(_val)
	{}
	short i;
	short j;
	short val;
};

struct Transform
{
    void merge(Transform& tr)
    {
        for(list<OneStepTransform>::iterator it=tr.t.begin(); it != tr.t.end(); it++ )
        {
            t.push_back(*it);
        }
    }
    list<OneStepTransform> t;
};

struct State
{
	State()
	{
        n = 0;
	}
	void init(int _n)
	{
        n = _n;
		isFinal = 0;
		m.resize(n);
		for(int i=0; i<n; i++)
		{
			m[i].resize(n);
		}
		used.resize(n*n+1);
	}
	void put(short i, short j, short val)
	{
		m[i][j] = val;
		used[val] = 1;
		isFinal++;
	}
	void remove(short i, short j)
	{
		used[m[i][j]] = 0;
		m[i][j] = 0;
		isFinal--;
	}
    void apply(const Transform& tr)
    {
        for(list<OneStepTransform>::const_iterator it=tr.t.begin(); it != tr.t.end(); it++)
        {
            put(it->i, it->j, it->val);
        }
    }
    void revert(const Transform& tr)
    {
        for(list<OneStepTransform>::const_iterator it=tr.t.begin(); it != tr.t.end(); it++)
        {
            remove(it->i, it->j);
        }
    }
    int n;
	vector< vector<short> > m;
	vector<short> used;
	int isFinal;
};


#endif

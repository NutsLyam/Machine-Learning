#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ctime> 
#include <cstdlib>

using namespace std;
struct Edge
{
	int from;
	int to;
	float s = 1;
	float a = 0;
	float r = 0;
	Edge(int from, int to, float s, float a, float r) :
		from(from), to(to), s(s), a(a), r(r)
	{}
};

float alfa = 0.6;

using Vector_of_edges = vector<Edge*>;
//split line 
template <class Container>
void split1(const std::string& str, Container& cont, char delim = '\t')
{
	std::stringstream ss(str);
	std::string token;
	cont.clear();
	while (std::getline(ss, token, delim)) {
		cont.push_back(token);
	}
}


float exponential_smoothing(float prev_value, float next_value, float alfa)
{
	return alfa * next_value + (1 - alfa) * prev_value;
}

void find_max12(Vector_of_edges row, float * max1, float * max2, int * max_index)
{
	float value;
	for (int j = 0; j < row.size(); j++)
	{
		value = row[j]->a + row[j]->s;
		if (value > *max1)
		{
			*max1 = value;
			*max_index = j;
		}
		else if (value > *max2)
		{
			*max2 = value;
		}
	}
}
//compare with diaginal elements
void cheak_diag(Edge * element, float*max1, float * max2, int*max_index, bool*flag, int _size)
{
	float value = element->a + element->s;
	if (value > *max1)
	{
		*max2 = *max1;
		*max1 = value;
		*max_index = _size;
		element->r = exponential_smoothing(element->r, element->s - *max2, alfa);
		*flag = true;  // true if diagonal element is changed
	}
	else if (value > *max2)
	{
		*max2 = value;
		//element->r = element->s - *max1;
		//*flag = true;
	}
}

void update_R(vector<Vector_of_edges> from_edge, vector<Edge*> diagonal)
{
	for (int i = 0; i < from_edge.size(); i++)
	{
		Vector_of_edges  row = from_edge[i];
		float max1 = -10000;
		float max2 = -10000;
		float value, new_r;
		int max_index = 0;
		bool diag_chanded = false;
		//
		find_max12(row, &max1, &max2, &max_index);
		cheak_diag(diagonal[i], &max1, &max2, &max_index, &diag_chanded, from_edge[i].size());

		//cout << from_edge[i].size() << endl;
		for (int k = 0; k < max_index; k++)
		{
			new_r = from_edge[i][k]->s - max1;
			from_edge[i][k]->r = exponential_smoothing(from_edge[i][k]->r, new_r, alfa);

		}
		for (int k = max_index + 1; k <from_edge[i].size(); k++)
		{
			new_r = from_edge[i][k]->s - max1;
			from_edge[i][k]->r = exponential_smoothing(from_edge[i][k]->r, new_r, alfa);
		}


		if (!diag_chanded)
		{
			diagonal[i]->r = exponential_smoothing(diagonal[i]->r, diagonal[i]->s - max1, alfa);
			new_r = from_edge[i][max_index]->s - max2;
			from_edge[i][max_index]->r = exponential_smoothing(from_edge[i][max_index]->r, new_r, alfa);
		}
	}
}
void update_A(vector<Vector_of_edges> to_edge, vector<Edge*> diagonal)
{
	for (int k = 0; k < to_edge.size(); k++)
	{
		{
			auto col = to_edge[k];
			float new_a, sum = 0;
			for (int j = 0; j < col.size(); j++)
			{
				sum += (col[j]->r > 0) * col[j]->r;
			}
			for (int j = 0; j < to_edge[k].size(); j++)
			{
				new_a = diagonal[k]->r + sum - col[j]->r * (col[j]->r > 0); // sum excluding element with j = i
				new_a = (new_a < 0) * new_a;
				to_edge[k][j]->a = exponential_smoothing(to_edge[k][j]->a, new_a, alfa);
			}
			//diagonal
			diagonal[k]->a = exponential_smoothing(diagonal[k]->a, sum, alfa);
		}
	}
}

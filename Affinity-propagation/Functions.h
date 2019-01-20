#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ctime> 
#include <cstdlib>
#include <limits>


using namespace std;
struct Edge
{
	int from;
	int to;
	double s = 1;
	double a = 0;
	double r = 0;
	Edge(int from, int to, double s, double a, double r) :
		from(from), to(to), s(s), a(a), r(r)
	{}
};

double alfa = 0.6, N = 196591;
clock_t start, finish;

using Vector_of_edges = vector<Edge*>;
//split line 
template <class Container>
void split1(const std::string& str, Container& cont, char delim = '\t')
{
	stringstream ss(str);
	string token;
	cont.clear();
	while (getline(ss, token, delim)) {
		cont.push_back(token);
	}
}

void read(vector<Vector_of_edges> from_edge, vector<Vector_of_edges> to_edge, vector<Edge*> diagonal)
{
	srand(time(NULL));
	vector<string> coordinates; // [from, to]
	string line;
	ifstream in;
	int x, y, N = 196591;
	int counter = 0;
	in.open("C:\\Users\\Nuts\\Documents\\Visual Studio 2015\\Projects\\AffinityPropagation\\Gowalla_edges.txt");
	if (in.is_open())
	{
		cout << "Begin!!" << endl;
		//counter = 0;
		start = clock();
		while (getline(in, line)) //&& counter < 1000) 
		{
			split1(line, coordinates); // split
			x = stoi(coordinates[0]);
			y = stoi(coordinates[1]);
			//cout << x << "  "  << y << endl;
			double press = 0.00001 * (rand() % 101);
			Edge *edge = new Edge(x, y, 1 + press, 0, 0); //pointer at the new edge, указатель  на новую считанную вершину
			from_edge[x].push_back(edge);
			to_edge[y].push_back(edge);
			//counter++;
		}
		finish = clock();
		cout << "Time of reading : " << (double)(finish - start) / ((double)CLOCKS_PER_SEC) << endl;
		in.close();
	}
 //diagonal
	for (int i = 0; i < N; i++)
	{
		double press = 0.00001 * (rand() % 101);
		Edge *edge = new Edge(i, i, -3 - press, 0, 0);
		diagonal.push_back(edge);
	}
}
double exponential_smoothing(double prev_value, double next_value, double alfa)
{
	return alfa * next_value + (1 - alfa) * prev_value;
}

void find_max12(Vector_of_edges row, double * max1, double * max2, int * max_index)
{
	//cout << row.size() << endl;
	double value;
	if (row.size() > 1) {
		for (int j = 0; j < row.size(); j++)
		{
			value = row[j]->a + row[j]->s;
			if (value > *max1)
			{
				*max2 = *max1;
				*max1 = value;
				*max_index = j;
			}
			else if (value > *max2)
			{
				*max2 = value;
			}
		}
	}
	//
	else {
		value = row[0]->a + row[0]->s;
		*max2 = *max1 = value;
		*max_index = 0;
	}
}
void update_R(vector<Vector_of_edges> from_edge, vector<Edge*> diagonal)
{
	for (int i = 0; i < from_edge.size(); i++)
	{
		Vector_of_edges  row = from_edge[i];
		
		double max1 = -numeric_limits<double>::infinity();
		double max2 = -numeric_limits<double>::infinity();
		double value, new_r;
		int max_index = 0;
		//

		find_max12(row, &max1, &max2, &max_index);
	
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

		diagonal[i]->r = exponential_smoothing(diagonal[i]->r, diagonal[i]->s - max1, alfa);
		new_r = from_edge[i][max_index]->s - max2;
		from_edge[i][max_index]->r = exponential_smoothing(from_edge[i][max_index]->r, new_r, alfa);

	}
}
void update_A(vector<Vector_of_edges> to_edge, vector<Edge*> diagonal)
{
	for (int k = 0; k < to_edge.size(); k++)
	{
		{
			auto col = to_edge[k];
			double new_a, sum = 0;
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
void choose_cluster(vector<Vector_of_edges> from_edge, vector<Edge*> diagonal)
{
	ofstream fout;
	vector<int> c(N);
	for (size_t i = 0; i < from_edge.size(); i++)
	{
		auto row = from_edge[i];
		double value, max = -numeric_limits<float>::infinity();;
		int index = 0;
		for (size_t j = 0; j < row.size(); j++)
		{
			value = row[j]->a + row[j]->r;
			if (value > max)
			{
				max = value;
				index = row[j]->to;
			}
		}

		// проверка диагональных элементов
		if (diagonal[i]->a + diagonal[i]->r > max)
		{
			index = diagonal[i]->to;
		}
		c[i] = index;

	}
	//write labels
	fout.open("fout_exemplars.txt");
	for (int i = 0; i < N; i++)
	{
		fout << i << " \t" << c[i] << endl;
	}
	fout.close();
}
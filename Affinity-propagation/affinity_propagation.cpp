#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ctime> 
#include <cstdlib>

#include "Functions.h"

using namespace std;

int main() {
	clock_t start, end;

	ifstream in;          
	in.open("C:\\Users\\Nuts\\Documents\\Visual Studio 2015\\Projects\\AffinityPropagation\\Gowalla_edges.txt"); 
	
	vector<string> coordinates; // [from, to]
	string line;
	int x, y, N = 196591; 
	int iteration = 0, Maxiteration = 500, counter = 0;
	vector<int> c(N);
	
	using Vector_of_edges = vector<Edge*>; //вектор из указателей на вершины
	vector<Vector_of_edges> from_edge(N); // row, вектор из векторов
	vector<Vector_of_edges> to_edge(N); // column
	vector<Edge*> diagonal;

	//вывод в файл
	ofstream fout;
	srand(time(NULL));
	
	// reading and filling 
	if (in.is_open())
	{
		cout << "Begin!!" << endl;
		//counter = 0;
		start = clock();
		while (getline(in, line)) //counter < 10) 
		{
			split1(line, coordinates); // split
			x = stoi(coordinates[0]) ;
			y = stoi(coordinates[1]);
			//cout << x << "  "  << y << endl;
			float press = 0.0001 * (rand() % 101);
			Edge *edge = new Edge(x,y,1 + press,0,0); //pointer at the new edge, указатель  на новую считанную вершину
			from_edge[x].push_back(edge);
			to_edge[y].push_back(edge);
			//counter++;
		}
		end = clock();
		cout << "Time of reading : " << (double)(end - start) / ((double)CLOCKS_PER_SEC) << endl;
		in.close();
	}
	for (int i = 0; i < N; i++)
	{
		float press = 0.0001 * (rand() % 101);
		Edge *edge = new Edge(i, i, -1+press, 0, 0);
		diagonal.push_back(edge);
	}
	while (iteration < Maxiteration) 
	{
		iteration++;
	// update matrix R
	 start = clock();
	update_R(from_edge, diagonal);
	
	//matrix A
	update_A(to_edge, diagonal);

	end = clock();
	cout << iteration<< " iteration : " << (double)(end - start) / ((double)CLOCKS_PER_SEC) << " s" << endl;
	}

//choose cluster labels
	for (size_t i = 0; i < from_edge.size(); i++)
	{
		auto row = from_edge[i];
		float value, max = -10000;
		int index = 1;
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

	system("pause");
	return 0;
}

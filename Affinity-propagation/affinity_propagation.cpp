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
//choose cluster labels


int main() {
	int iteration = 0, Maxiteration = 600;
	

	using Vector_of_edges = vector<Edge*>; //вектор из указателей на вершины
	vector<Vector_of_edges> from_edge(N); // row, вектор из векторов
	vector<Vector_of_edges> to_edge(N); // column
	vector<Edge*> diagonal;

	/*x
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
			double press = 0.0001 * (rand() % 101);
			Edge *edge = new Edge(x, y, 1 + press, 0, 0); //pointer at the new edge, указатель  на новую считанную вершину
			from_edge[x].push_back(edge);
			to_edge[y].push_back(edge);
			//counter++;
		}
		finish = clock();
		cout << "Time of reading : " << (double)(finish - start) / ((double)CLOCKS_PER_SEC) << endl;
		in.close();

		for (int i = 0; i < N; i++)
		{
			double press = 0.0001 * (rand() % 101);
			Edge *edge = new Edge(i, i, -1 - press, 0, 0);
			diagonal.push_back(edge);
		}
	}
*/
read(from_edge, to_edge, diagonal);

	// reading and filling 

	while (iteration < Maxiteration)
	{
		iteration++;
		// update matrix R
		start = clock();
		update_R(from_edge, diagonal);

		//matrix A
		update_A(to_edge, diagonal);

		finish = clock();
		cout << iteration << " iteration : " << (double)(finish - start) / ((double)CLOCKS_PER_SEC) << " s" << endl;

		if (iteration%50==0)
		{
			choose_cluster(from_edge, diagonal);
		}
	}
	choose_cluster(from_edge, diagonal);
	system("pause");
	return 0;
}

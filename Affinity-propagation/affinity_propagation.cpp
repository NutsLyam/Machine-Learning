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

void cheak_diag(Edge * element, float*max1, float * max2, int*max_index, bool*flag, int _size)
{
	float value = element->a + element->s;
	if (value > *max1)
	{
		*max2 = *max1;
		*max1 = value;
		*max_index = _size;
		element->r = exponential_smoothing(element->r , element->s - *max2, alfa);
		*flag = true;  // изменили диагональный элемент
	}
	else if (value > *max2)
	{
		*max2 = value;
		//element->r = element->s - *max1;
		//*flag = true;
	}
}

void update_R(  vector<Vector_of_edges> from_edge, vector<Edge*> diagonal)
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
				new_a  = diagonal[k]->r + sum - col[j]->r * (col[j]->r > 0); // сумма без элемeнта где j = i
				new_a =  (new_a < 0) * new_a;
				to_edge[k][j]->a = exponential_smoothing(to_edge[k][j]->a, new_a, alfa);
			}
			//diagonal
			diagonal[k]->a = exponential_smoothing(diagonal[k]->a,sum, alfa);
		}
	}
}

int main() {
	clock_t start, end;

	ifstream in;          // поток для записи
	in.open("C:\\Users\\Nuts\\Documents\\Visual Studio 2015\\Projects\\AffinityPropagation\\Gowalla_edges.txt"); // открываем файл
	
	vector<string> coordinates; // [from, to]
	string line;
	int x, y, N = 196591; 
	int iteration = 0, Maxiteration = 500, counter = 0;
	vector<int> c(N);
	
	using Vector_of_edges = vector<Edge*>; //вектор из указателей на вершины
	vector<Vector_of_edges> from_edge(N); // вектор из векторов
	vector<Vector_of_edges> to_edge(N);
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
			split1(line, coordinates); // разделение строки
			x = stoi(coordinates[0]) ;
			y = stoi(coordinates[1]);
			//cout << x << "  "  << y << endl;
			float press = 0.0001 * (rand() % 101);
			Edge *edge = new Edge(x,y,1 + press,0,0); //указатель  на новую считанную вершину
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

	for (size_t i = 0; i < from_edge.size(); i++)
	{
		auto row = from_edge[i];
		float value, max = -10000;
		int index = 1;
		for (int j = 0; j < row.size(); j++)
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

	fout.open("fout_exemplars.txt");
	for (int i = 0; i < N; i++)
	{
		fout << i << " \t" << c[i] << endl;
	}
	fout.close();

		
fout.open("sra_matrix.txt");
		for (int i = 0; i < from_edge.size(); i++){
			//запись диагонального элемента
			fout << diagonal[i]->from << " \t" << diagonal[i]->to << " \t" << diagonal[i]->s << " \t" <<
				diagonal[i]->r << " \t" << diagonal[i]->a << endl;
			for (int j = 0; j < from_edge[i].size(); j++)
			{
				fout << from_edge[i][j]->from << " \t" << from_edge[i][j]->to << " \t" <<
					from_edge[i][j]->s << " \t" << from_edge[i][j]->r << " \t" << from_edge[i][j]->a << endl;
			}
		}
		fout.close();

	system("pause");
	return 0;
}

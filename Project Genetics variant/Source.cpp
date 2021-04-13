#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <string>
#include <unordered_map>
#include <algorithm>
using namespace std;
int n, ncarro;
int localPercentage = 40;
vector < int > caps;
vector <double> Demands;
vector < vector <double> > matrizC;
int Crom;
int contador = 1, cuanti;
double testingtime =0;


struct solution {
	vector <int> route;
	double latency;
	//Default constructor:
	solution(){}
	//Constructor with route as first param:
	solution(vector <int> route) : solution() {
		this->route = route;
	}
	//Constructor with route and latency as params:
	solution(vector <int> route, double latency) : solution(route) {
		this->latency = latency;
	}
};

bool compareByLatency(const solution& a, const solution& b)
{
	return a.latency < b.latency;
}

void reset()
{
	caps.clear();
	Demands.clear();
	matrizC.clear();
}
bool busqueda(int valor, vector <int> buscar)
{
	bool por = false;

	for (int i = 0; i < buscar.size(); i++)
	{
		if (buscar[i] == valor)
		{
			por = true;
			break;
		}
	}

	return por;
}

double calcLatency(vector<int> soli)
{
	double latency = 0;

	int id = 1;
	for (int j = 1; j < soli.size(); j++)
	{
		if (soli[j] != 0)
		{
			double time = 0;

			for (int k = id; k <= j; k++)
			{
				time += matrizC[soli[k - 1]][soli[k]];

			}
			latency += time;
		}
		else
		{
			id = j + 1;
		}
	}

	return latency;

}

void theFixer(vector <int>& tofix, vector <int> &capas) {

	vector <int> copy = tofix;
	tofix = copy;
	vector <int> copyCaps = capas;
	capas = copyCaps;

	if (tofix[tofix.size() - 1] != 0) {
		tofix.emplace_back(0);
	}
	//Obtain ones that are missing:
	vector <int> faltan;
	unordered_map <int, double> clientCapacities;
	for (int i = 1; i <= n; i++) {
		if (!busqueda(i, tofix)) {
			faltan.emplace_back(i);
			clientCapacities[i] = Demands[i];
		}
	}

	//Obtain all car indexes
	vector <int> indcars;
	indcars.emplace_back(0);
	for (int i = 1; i < tofix.size(); i++) {
		if (tofix[i] == 0) {
			indcars.emplace_back(i);
		}
	}
	//The good stuff:
	int dond;

	while (faltan.size() != 0) {
		//Pick car:
		int indixcarro = rand() % ((indcars.size() - 2) - 0 + 1) + 0;
		int indixcarrofix = indcars[indixcarro];
		if (capas[indixcarro] > 0) {
			double toSubstract = capas[indixcarro] < clientCapacities[faltan[0]] ? capas[indixcarro] : clientCapacities[faltan[0]];
			if (indixcarrofix + 1 == indcars[indixcarro + 1]) {
				dond = 1 + indixcarrofix;
			}
			else
			{
				dond = rand() % ((indcars[indixcarro + 1] - 1) - (indixcarrofix + 1) + 1) + (indixcarrofix + 1);
			}
			capas[indixcarro] -= toSubstract;
			clientCapacities[faltan[0]] -= toSubstract;
			tofix.insert(tofix.begin() + dond, faltan[0]);
			//Push indexes 1 space to right:
			for (int i = indixcarro + 1; i < indcars.size(); i++) indcars[i]++;
			if(!clientCapacities[faltan[0]])
				faltan.erase(faltan.begin() + 0);
		}
	}

	
}

vector <int> Mutate(vector<int> cromida)
{
	vector <int> mutated, route, capasi = caps;

	int verificado = n + 1, vamos = 0;

	for (int i = 0; i < cromida.size(); i++)
	{
		if (cromida[i] != 0)
		{
			route.emplace_back(verificado - cromida[i]);
		}
	}

	mutated.emplace_back(0);
	bool banderaza = false;
	vector <int> indices;
	while (route.size() != 0)
	{

		if (vamos == capasi.size())
		{
			theFixer(mutated, capasi);
			break;
		}

		if (banderaza == true)
		{
			indices.clear();
			indices.emplace_back(0);
			for (int i = 1; i < mutated.size(); i++)
			{
				if (mutated[i] == 0)
				{
					indices.emplace_back(i);
				}
			}
		}

		if (capasi[vamos] - Demands[route[0]] >= 0)
		{
			capasi[vamos] -= Demands[route[0]];
			if (banderaza == false)
			{
				mutated.emplace_back(route[0]);
			}
			else
			{
				mutated.insert(mutated.begin() + indices[vamos] + 1, route[0]);
			}
			route.erase(route.begin() + 0);

			if (route.size() == 0 && banderaza == false)
			{
				mutated.emplace_back(0);
			}
		}
		else
		{
			bool bandera = true;
			while (bandera == true)
			{
				bandera = false;
				int ran = rand() % (2 - 1 + 1) + 1;
				if (ran == 1)
				{
					for (int j = 0; j < route.size(); j++)
					{
						if (capasi[vamos] - Demands[route[j]] >= 0)
						{
							capasi[vamos] -= Demands[route[j]];
							if (banderaza == false)
							{
								mutated.emplace_back(route[j]);
							}
							else
							{
								mutated.insert(mutated.begin() + indices[vamos] + 1, route[j]);
							}
							route.erase(route.begin() + j);
							bandera = true;
							break;
						}
					}
				}
				else
				{
					for (int j = route.size() - 1; j >= 0; j--)
					{
						if (capasi[vamos] - Demands[route[j]] >= 0)
						{
							capasi[vamos] -= Demands[route[j]];
							if (banderaza == false)
							{
								mutated.emplace_back(route[j]);
							}
							else
							{
								mutated.insert(mutated.begin() + indices[vamos] + 1, route[j]);
							}
							route.erase(route.begin() + j);
							bandera = true;
							break;
						}
					}

				}
			}
			vamos++;
			if (banderaza == false)
			{
				mutated.emplace_back(0);
			}
		}
	}


	return mutated;
}
vector <int >  cruzar(vector <int> padre1, vector <int> padre2)
{
	
	vector <int > hijo, route;
	route.reserve(n);

	vector <int> acom(n +1,0);
	route.clear();
	hijo.emplace_back(0);
	acom[0] = 1;
	int enturno = 0;
	vector <int> capas = caps;
	int debuga = padre1.size() < padre2.size() ? padre1.size() : padre2.size();

	for (int i = 1; i < debuga; i++)
	{
		//Armar ruta:
		int randi = rand() % 2;
		int pick1 = padre1[i], pick2 = padre2[i];
		if (randi && !acom[pick1]) {
			route.emplace_back(pick1);
			acom[pick1] = 1;
		}
		else {
			if (!acom[pick2]) {
				route.emplace_back(pick2);
				acom[pick2] = 1;
			}
		}
	}

	while (route.size() != n)
	{
		//Faltaron de colocar
		int ran = rand() % 2;
		if (ran)
		{
			for (int i = 1; i <= n; i++)
			{
				if (acom[i] == 0)
				{
					route.emplace_back(i);
					acom[i] = 1;
					break;
				}
			}
		}
		else
		{
			for (int i = n; i >= 1; i--)
			{
				if (acom[i] == 0)
				{
					route.emplace_back(i);
					acom[i] = 1;
					break;
				}
			}
		}
	}

	enturno = 0;
	while (route.size() != 0)
	{
		if (enturno > ncarro - 1)
		{
			//The real good stuff:
			theFixer(hijo, capas);
			break;

		}
		if (capas[enturno] - Demands[route[0]] < 0)
		{
			enturno++;
			hijo.emplace_back(0);
		}
		else
		{
			hijo.emplace_back(route[0]);
			capas[enturno] -= Demands[route[0]];
			route.erase(route.begin() + 0);
		}

	}
	int ran = rand() % (100 - 1 + 1) + 1;
	if (ran <= -1)
	{
		hijo = Mutate(hijo);
	}
	

	return hijo;
}
vector<solution> GenerarHijos(vector <solution> Population)
{
	int vamos = 0;
	int frent = 0;

	vector <solution> NewG;
	for (int k = 0; k < Crom; k++)
	{
		//Select fathers:
		vector <int> parents[2];
		parents[0] = Population[rand() % (Crom) + 0].route;
		parents[1] = Population[rand() % (Crom) + 0].route;

		//Cruzamiento
		vector <int> newChrom = cruzar(parents[0], parents[1]);
		NewG.emplace_back(newChrom, calcLatency(newChrom));
	}

	//See the better one
	vector <solution> dieBesten;
	dieBesten.reserve(Population.size() + NewG.size());
	dieBesten.insert(dieBesten.end(), Population.begin(), Population.end());
	dieBesten.insert(dieBesten.end(), NewG.begin(), NewG.end());
	//sort based on latencies:
	sort(dieBesten.begin(), dieBesten.end(), compareByLatency);
	//Keep the better Ones:
	dieBesten.resize(Crom);
	dieBesten.shrink_to_fit();

	return dieBesten;
}



void Local(solution &hijo)
{
	vector<int> ho, &orig = hijo.route;
	double currlat = hijo.latency;
	vector<int> indices;
	vector<double> capis;
	double capa = 0;
	indices.emplace_back(0);
	for (int i = 1; i < orig.size(); i++)
	{
		if (orig[i] == 0)
		{
			indices.emplace_back(i);
			capis.emplace_back(capa);
			capa = 0;
		}
		else
		{
			capa += Demands[orig[i]];
		}

	}

	//Busqueda intra:
	for (int wann = 0; wann < indices.size() - 1; wann++)
	{
		for (int idi = indices[wann] + 1; idi != indices[wann + 1] - 1; idi++)
		{
			for (int i = idi + 1; i < indices[wann + 1]; i++)
			{
				ho = orig;
				ho[idi] = orig[i];
				ho[i] = orig[idi];
				double nuevlat = calcLatency(ho);
				if (currlat > nuevlat)
				{
					orig = ho;
					currlat = nuevlat;
					goto endofintra;
				}
			}
		}
	}
endofintra:
	//Busqueda interna 1:
	for (int wann = 0; wann != indices.size() - 1; wann++)
	{
		for (int ou = wann + 1; ou < indices.size(); ou++)
		{
			if (ou + 1 > indices.size() - 1)
			{
				break;
			}
			for (int idi = indices[wann] + 1; idi != indices[wann + 1] - 1; idi++)
			{
				for (int doda = indices[ou] + 1; doda != indices[ou + 1] - 1; doda++)
				{
					if ((capis[wann] - Demands[orig[idi]] + Demands[orig[doda]]) < caps[wann] && (capis[ou] - Demands[orig[doda]] + Demands[orig[idi]]) < caps[ou])
					{
						ho = orig;
						ho[idi] = orig[doda];
						ho[doda] = orig[idi];
						double nuevlat = calcLatency(ho);
						if (currlat > nuevlat)
						{
							capis[wann] = capis[wann] - Demands[orig[idi]] + Demands[orig[doda]];
							capis[ou] = capis[ou] - Demands[orig[doda]] + Demands[orig[idi]];
							orig = ho;
							currlat = nuevlat;
							goto endofinter1;
						}
					}
				}
			}
		}
	}
endofinter1:
	//Busqueda interna 2:

	for (int wann = 0; wann != indices.size() - 1; wann++)
	{
		for (int ou = wann + 1; ou < indices.size(); ou++)
		{
			if (ou + 1 > indices.size() - 1)
			{
				break;
			}
			for (int idi = indices[wann] + 1; idi != indices[wann + 1] - 1; idi++)
			{
				for (int doda = indices[ou] + 1; doda != indices[ou + 1] - 1; doda++)
				{

					if (capis[ou] + Demands[orig[idi]] < caps[ou])
					{
						ho = orig;
						ho.insert(ho.begin() + doda, ho[idi]);
						ho.erase(ho.begin() + idi);
						double nuevlat = calcLatency(ho);
						if (currlat > nuevlat)
						{
							capis[wann] -= Demands[orig[idi]];
							capis[ou] += Demands[orig[idi]];
							orig = ho;
							currlat = nuevlat;
							goto endofinter2;
						}
					}

				}
			}
		}
	}
endofinter2:
	hijo.latency = currlat;
}



vector <solution> FirstGen(int Crom)
{
	vector <solution> Population;
	Population.reserve(Crom);
	vector <int> Chromosome;
	for (int i = 0; i < Crom; i++)
	{
		while (true)
		{
			vector <int> Chromosome;
			Chromosome.emplace_back(0);
			vector <int> capturn = caps;
			vector <bool> placed(n + 1, 0);
			int cuantos = 0;
			placed[0] = 1;
			int carro = 0;
			int randi = 0;
			for (int j = 0; j <= n + ncarro; j++)
			{
				while (placed[randi])
					randi = rand() % (n - 1 + 1) + 1;

				if (Demands[randi] > capturn[carro])
				{	
					Chromosome.emplace_back(0);
					if (carro + 1 == ncarro) {
						theFixer(Chromosome, capturn);
						goto fixed;
					}
					carro++;
				}
				else
				{
					placed[randi] = 1;
					cuantos++;

					Chromosome.emplace_back(randi);
					capturn[carro] -= Demands[randi];
					if (cuantos == n)
					{
						Chromosome.emplace_back(0);
						break;
					}
				}

			}
			if (cuantos == n )
			{
				fixed:
				Population.emplace_back(Chromosome, calcLatency(Chromosome));
				break;
			}
		}

	}
	return Population;
}


solution AlgoritmoBueno(vector<solution> Population)
{
	solution best;
	best.latency = INFINITY;
	//double localtime = 0, paretotime = 0, generarhijostime = 0, inicio = 0;
	while (contador <= cuanti)
	{
		//system("cls");
		for (int sola = 0; sola < Population.size(); sola++)
		{
			//Busqueda local:
			int randi = rand() % (100 - 1 + 1) + 1;
			//inicio = clock();
			
			if (randi <= localPercentage)
			{
				Local(Population[sola]);
			}
			//localtime += (clock() - inicio) / (double)CLOCKS_PER_SEC;
			//Catalogar el cromosoma/hijo en los frentes
			//inicio = clock();
			//paretotime += (clock() - inicio) / (double)CLOCKS_PER_SEC;

			if (Population[sola].latency < best.latency) {
				best.route = Population[sola].route;
				best.latency = Population[sola].latency;
			}

		}

		//Generar nueva población:

		//inicio = clock();
		Population = GenerarHijos(Population);
		//generarhijostime += (clock() - inicio) / (double)CLOCKS_PER_SEC;
		
		contador++;

	}
	//cout << generarhijostime << " local: " << localtime;
	contador = 1;
	return best;
}

int main()
{
	
	srand(time(NULL));
	vector <string> inst;

	/*inst.emplace_back("CMT1.txt");
	inst.emplace_back("CMT2.txt");
	inst.emplace_back("CMT3.txt");
	inst.emplace_back("CMT4.txt");*/
	inst.emplace_back("CMT5.txt");
	inst.emplace_back("CMT11.txt");
	inst.emplace_back("CMT12.txt");

	vector <string> nombres;

	/*nombres.emplace_back("CMT1");
	nombres.emplace_back("CMT2");
	nombres.emplace_back("CMT3");
	nombres.emplace_back("CMT4");*/
	nombres.emplace_back("CMT5");
	nombres.emplace_back("CMT11");
	nombres.emplace_back("CMT12");

	for (int itar = 0; itar < inst.size(); itar++)
	{
		reset();
		cout.flush();
		cout << "Problema: " << nombres[itar];
		cout << "\nChromosomes: ", cin >> Crom;
		cout << "\nGenerations: ", cin >> cuanti;
		ifstream data(inst[itar]);
		double cap = 0,demanda = 0;
		data >> n >> cap;
		vector <vector <double>> matrizD;
		vector <double> xs;
		vector <double> ys;
		xs.reserve(n + 1);
		ys.reserve(n + 1);
		double x, y, dem;
		data >> x >> y;
		xs.emplace_back(x);
		ys.emplace_back(y);
		Demands.reserve(n + 1);
		Demands.emplace_back(0);

		for (int i = 0; i < n; i++)
		{
			data >> x >> y >> dem;
			xs.emplace_back(x);
			ys.emplace_back(y);
			demanda += dem;
			Demands.emplace_back(dem);
		}

		ncarro = ceil(demanda / cap);
		caps.resize(ncarro, cap);


		for (int i = 0; i <= n; i++)
		{
			vector <double> dists;
			for (int j = 0; j <= n; j++)
			{

				double dist = sqrt(pow(xs[j] - xs[i], 2) + pow(ys[j] - ys[i], 2));
				dists.emplace_back(dist);
			}

			matrizD.emplace_back(dists);
		}

		for (int i = 0; i <= n; i++)
		{
			vector <double> costs;
			for (int j = 0; j <= n; j++)
			{
				double costo = matrizD[i][j];
				costs.emplace_back(costo);
			}
			matrizC.emplace_back(costs);
		}

		/*vector <int> test{ 0, 12, 5, 49, 10, 30, 34, 21, 29, 3 ,0 };

		cout << calcLatency(test);
		cin >> n;*/

		solution best = AlgoritmoBueno(FirstGen(Crom));
		cout << endl;
		for (int& i : best.route)
			cout << i << " ";
		cout << endl << best.latency;
		cin >> Crom;
	}





}
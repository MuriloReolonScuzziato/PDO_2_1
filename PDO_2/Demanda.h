#pragma once
#include <string>
using std::string;

#include <vector>
using std::vector;
typedef vector<double> vetorfloat;
typedef vector<vetorfloat> vetorfloat2;

#include <sstream> 
using std::ostringstream;

class CDemanda
{
	string nome_demanda;
	//vetorfloat D;
	vetorfloat2 cen_demanda;
	//vetorfloat D_1, D_2;
public:
	CDemanda(int n_demanda_a);
	~CDemanda(void);
	string GetNome() {return nome_demanda;}
	vetorfloat GetCenD(int cen) {return cen_demanda[cen];}
	//double GetD(int t_a) {return cen_demanda[0][t_a];}
	double GetD(int cen, int t_a) {return cen_demanda[cen][t_a];}
	//vetorfloat GetD1() {return D_1;}
	//vetorfloat GetD2() {return D_2;}
	void SetDemanda(vetorfloat D_a) {cen_demanda.push_back(D_a);}
	
};

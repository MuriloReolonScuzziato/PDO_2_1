#pragma once

#include "Barra.h"

class CLinha
{
	string nome_linha;
	int ind_de_barra, ind_para_barra;
	double reatancia, capacidade;
	
public:
	CBarra * de_barra, * para_barra;
	
	CLinha(string, CBarra * , CBarra * , double, double, int, int);
	~CLinha(void);

	string GetNome() {return nome_linha;}
	string GetDeBarra() {return de_barra->GetNome();}
	string GetParaBarra() {return para_barra->GetNome();}
	//int GetDeBarra() {return de_barra;}
	//int GetParaBarra() {return para_barra;}
	double GetReatancia() {return reatancia;}
	double GetCapacidade() {return capacidade;}
	int GetIndDeBarra() {return ind_de_barra;}
	int GetIndParaBarra() {return ind_para_barra;}
};


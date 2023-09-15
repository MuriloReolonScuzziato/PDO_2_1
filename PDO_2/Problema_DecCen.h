#pragma once

#include "SubProblemaCen.h"
#include <time.h>
#include <windows.h>

#include <eigen>
using namespace Eigen;

//#include "BundleClean_1709.h"

typedef vector<LARGE_INTEGER> vetorTempo;

namespace met_DecCen
{

typedef vector<CSubProblemaCen *> vetorCenarios;

class CProblema_DecCen
{
	GRBEnv ambGRB;     // Ambiente Gurobi
	vetorCenarios problemaCenario;		// Vetor de ponteiros para os subproblemas

	//CBundle * bundle;

	double fdual_RL, norma_sg, fdual_RP, norma_g;
	int n_var_duplicadas, n_var_est1;
	vetorfloat2 x, Lambda, subgrad, grad;		// os vetores Lambda, subgrad e grad tem o tamanho do numero de variáveis duplicadas
	vetorfloat x_med, mi, Lambda_med;		// x_med: o mesmo valor para todos cenários (um valor para cada tipo de variável em cada período)
	CSistema * sistema_a;
	int flag1, flag4, JJ;

	LARGE_INTEGER clockPerSec;
	LARGE_INTEGER inicio_RL, fim_RL, inicio_Modelo, fim_Modelo, inicio_Resol, fim_Resol, inicio_Gravar_Escrever, fim_Gravar_Escrever;
	vetorTempo inicio_cenarios, fim_cenarios;
	vetorfloat tempo_cenarios;
	CResultados * resultadosGurobi;
	int it_RL, it_RP;
	double tempoRL, tempoModelo, tempoResolucao, tempoGravar_Escrever;

	void oracle(const vetorfloat2 &lambda, double &fdual, vetorfloat2 &subgradiente);
	void subgradiente(int &itmax);
	void oracle(const VectorXd &lambda_e, double &fdual, VectorXd &subgradiente_e);

public:
	CProblema_DecCen(CSistema * const sistema_end, CResultados * const resultadosGurobi_end);
	~CProblema_DecCen(void);

	void ResolverRL();
};

};
#pragma once

// Subproblema para otimização de cada cenário

#include "Sistema.h"
#include "MatrizEsparsa.h"
#include "ResultadosConj.h"

#include "gurobi_c++.h"
using namespace std;

namespace met_DecCen
{

class CSubProblemaCen
{
	GRBModel modeloGRB;  // Modelo Gurobi
	GRBVar * vars;

	CSistema * sistema_a;
	int n_restricoes, n, cenario;

	vetorfloat x, L;
	double fo;

	void CriarVariaveis();
	void CriarFuncaoObjetivoRL(const int cenario, const vetorfloat * const Lambda);
	void CriarFuncaoObjetivoRP(const int, const vetorfloat * const, const vetorfloat * const, const vetorfloat * const);
	void CriarRestricoes();
	void AtualizarRestricoes();

public:
	// 2 constructores, um para criar modelo e outro para copiar modelo
	CSubProblemaCen(CSistema * const sistema_end, int cenario_a, const GRBEnv &ambiente_gurobi);
	CSubProblemaCen(CSistema * const sistema_end, int cenario_a, CSubProblemaCen &cenario0);
	~CSubProblemaCen(void);

	int GetN_Restricoes() {return n_restricoes;}
	int GetN() {return n;}
	
	int ResolverProblemaRL(CResultados * resultadoGurobi, const vetorfloat * const lambda, const int iter);
	int ResolverProblemaRP(CResultados * resultadoGurobi, const vetorfloat * const lambda, const vetorfloat * const mi, const vetorfloat * const x_med, const int iter);
};

};
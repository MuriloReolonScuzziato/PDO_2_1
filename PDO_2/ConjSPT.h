#pragma once

#include "ResultadosConj.h"

class ConjSPT
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada usina)
	vetorGRBVar vars;

	CSistema * sistema_a;
	int n_usinas, flag4, flag7;
	vetorint n;
	double fo_i;
	vetorfloat x_i;

	void CriarVariaveis(int nu);
	void CriarRestricoes(int nu);
	void CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int nu);
	void CriarFuncaoObjetivoRL2(const vetorfloat * const lambda, int nu);
	void FixarCondIniciais(int nu);

	CMatrizEsparsa MatrizTup(int n_a, int nu);
	CMatrizEsparsa MatrizTdown(int n_a, int nu);
	CMatrizEsparsa MatrizTupDown(int n_a, int nu);
	CMatrizEsparsa MatrizRampaUp(int n_a, int nu);
	CMatrizEsparsa MatrizRampaDown(int n_a, int nu);
	CMatrizEsparsa MatrizLimPtMin(int n_a, int nu);
	CMatrizEsparsa MatrizLimPtMax(int n_a, int nu);
	CMatrizEsparsa MatrizRestCP(int n_a, int nu);
	CMatrizEsparsa MatrizCortesF(int n_a, int nu);
	vetorfloat2 LimTup(int n_a, int nu);
	vetorfloat2 LimTdown(int n_a, int nu);
	vetorfloat2 LimTupDown(int n_a, int nu);
	vetorfloat2 LimRampaUp(int n_a, int nu);
	vetorfloat2 LimRampaDown(int n_a, int nu);
	vetorfloat2 LimPtMin(int n_a, int nu);
	vetorfloat2 LimPtMax(int n_a, int nu);
	vetorfloat2 LimRestCP(int n_a, int nu);
	vetorfloat2 LimCortesF(int n_a, int nu);
	void MatrizRestricoesLineares(int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int nu, int &n_restricoes);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int nu);

public:
	ConjSPT(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	~ConjSPT(void);

	int ResolverProblemaRL(ResultadosConj * resultadoGurobi, const vetorfloat * const lambda, const int iter, int usina);
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
};
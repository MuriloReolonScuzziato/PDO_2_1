#pragma once

#include "ResultadosConj.h"

class ConjSPHE
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada usina)
	vetorGRBVar vars;
	
	CSistema * sistema_a;
	int n_usinas;
	vetorint n;
	double fo_rn;
	vetorfloat x_rn;

	void CriarVariaveis(int nu);
	void CriarRestricoes(int nu);
	void CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int nu);

	CMatrizEsparsa MatrizLimPhgMin(int n_a, int nu);
	CMatrizEsparsa MatrizLimPhgMax(int n_a, int nu);
	CMatrizEsparsa MatrizLimqMin(int n_a, int nu);
	CMatrizEsparsa MatrizLimqMax(int n_a, int nu);
	CMatrizEsparsa MatrizBalPotencia(int n_a, int nu);
	CMatrizEsparsa MatrizBalVazao(int n_a, int nu);
	CMatrizEsparsa MatrizFuncProd(int n_a, int nu);
	CMatrizEsparsa MatrizPhMax(int n_a, int nu);
	vetorfloat2 LimPhgMin(int n_a, int nu);
	vetorfloat2 LimPhgMax(int n_a, int nu);
	vetorfloat2 LimQMin(int n_a, int nu);
	vetorfloat2 LimQMax(int n_a, int nu);
	vetorfloat2 LimBalPotenciaL(int n_a, int nu);
	vetorfloat2 LimBalVazaoL(int n_a, int nu);
	vetorfloat2 LimFuncProdL(int n_a, int nu);
	vetorfloat2 LimPhMax(int n_a, int nu);
	void MatrizRestricoesLineares(int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int nu, int &n_restricoes);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int nu);
	
public:
	ConjSPHE(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	~ConjSPHE(void);

	int ResolverProblemaRL(ResultadosConj * resultadoGurobi, const vetorfloat * const lambda, const int iter, int usina, int no);
	double GetFO() { return fo_rn;}
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
};
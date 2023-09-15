#pragma once

#include "Spcdec2Results.h"

class Spcdec2SPH
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada cascata)
	vetorGRBVar vars;
	
	CSistema * sistema_a;
	int n_cascatas, T, flag2, flag3;
	vetorint n;
	double fo_ch;
	vetorfloat x_ch;
	vetorint2 cascata, cascata_jusante;		// vetor com os indices das hidros de cada cascata

	void CriarVariaveis(int ch);
	void CriarRestricoes(int ch);
	void CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int ch);

	CMatrizEsparsa MatrizBalHid(int n_a, int ch);
	CMatrizEsparsa MatrizVmeta(int n_a, int ch);
	CMatrizEsparsa MatrizLimPhgMin(int n_a, int ch);
	CMatrizEsparsa MatrizLimPhgMax(int n_a, int ch);
	CMatrizEsparsa MatrizLimqMin(int n_a, int ch);
	CMatrizEsparsa MatrizLimqMax(int n_a, int ch);
	CMatrizEsparsa MatrizBalPotencia(int n_a, int ch);
	CMatrizEsparsa MatrizBalVazao(int n_a, int ch);
	CMatrizEsparsa MatrizFuncProd(int n_a, int ch);
	CMatrizEsparsa MatrizPhMax(int n_a, int ch);
	vetorfloat2 LimBalHid(int n_a, int ch);
	vetorfloat2 LimVmeta(int n_a, int ch);
	vetorfloat2 LimPhgMin(int n_a, int ch);
	vetorfloat2 LimPhgMax(int n_a, int ch);
	vetorfloat2 LimQMin(int n_a, int ch);
	vetorfloat2 LimQMax(int n_a, int ch);
	vetorfloat2 LimBalPotenciaL(int n_a, int ch);
	vetorfloat2 LimBalVazaoL(int n_a, int ch);
	vetorfloat2 LimFuncProdL(int n_a, int ch);
	vetorfloat2 LimPhMax(int n_a, int ch);

	void MatrizRestricoesLineares(int n_a, CMatrizEsparsa &MM, int ch);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int ch);
	
public:
	Spcdec2SPH(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	~Spcdec2SPH(void);

	int ResolverProblemaRL(Spcdec2Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int ch);
	int GetNCascatas() {return n_cascatas;}
	double GetFO() { return fo_ch;}
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
};
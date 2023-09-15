#pragma once

#include "Spcdec1Results.h"

class Spcdec1SPH
{
	GRBModel * modelosGRB;  // Modelo do Gurobi (um modelo somente)
	GRBVar * vars;
	
	CSistema * sistema_a;
	int n, T, flag2;
	double fo;
	vetorfloat x;

	void CriarVariaveis();
	void CriarRestricoes();
	void CriarFuncaoObjetivoRL(const vetorfloat * const lambda);

	CMatrizEsparsa MatrizBalHid();
	CMatrizEsparsa MatrizVmeta();
	CMatrizEsparsa MatrizFuncProd();
	CMatrizEsparsa MatrizBalVazao();
	CMatrizEsparsa MatrizBalPotencia();
	CMatrizEsparsa MatrizLimPhgMin();
	CMatrizEsparsa MatrizLimPhgMax();
	CMatrizEsparsa MatrizLimqMin();
	CMatrizEsparsa MatrizLimqMax();
	CMatrizEsparsa MatrizReserva();

	vetorfloat2 LimBalHid();
	vetorfloat2 LimVmeta();
	vetorfloat2 LimFuncProdL();
	vetorfloat2 LimBalVazaoL();
	vetorfloat2 LimBalPotenciaL();
	vetorfloat2 LimPhgMin();
	vetorfloat2 LimPhgMax();
	vetorfloat2 LimQMin();
	vetorfloat2 LimQMax();
	vetorfloat2 LimReserva();
	void MatrizRestricoesLineares(CMatrizEsparsa &MM);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);
	
public:
	Spcdec1SPH(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	~Spcdec1SPH(void);

	int ResolverProblemaRL(Spcdec1Results * resultadoGurobi, const vetorfloat * const lambda, const int iter);
	double GetFO() { return fo;}
	int GetNVarPorComp() {return n;}		// comp de cada tipo de subproblema
};
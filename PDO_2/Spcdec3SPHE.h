#pragma once

#include "Spcdec3Results.h"

class Spcdec3SPHE
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
	void MatrizRestricoesLineares(int n_a,  CMatrizEsparsa &MM, int nu);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int nu);
	
public:
	Spcdec3SPHE(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	~Spcdec3SPHE(void);

	int ResolverProblemaRL(Spcdec3Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int usina, int no);
	double GetFO() { return fo_rn;}
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
	void SetPrecision( double precision);
	double GetLBvar(int comp, int var) {return vars[comp][var].get(GRB_DoubleAttr_LB);}

};
#pragma once

#include "Spcdec3Results.h"

#include "OPTUtils.h"

#include <algorithm>
using std::find;

class Spcdec3SPHA
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada cascata)
	vetorGRBVar vars;
	
	CSistema * sistema_a;
	int n_cascatas, N;
	vetorint n;
	double fo_ch;
	vetorfloat x_ch;
	vetorint2 cascata, cascata_jusante;		// vetor com os indices das hidros de cada cascata

	// Atributos somente usados se o subproblema for Easy Component
	vetorint2 Mbeg, Mind, dtipo;
	vetorfloat2 Mval, dvalor;
	vetorint Nrestricoes;
	//vetorfloat2 var_ub, var_lb;

	void CriarVariaveis(int ch);
	void CriarRestricoes(int ch);
	void CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int ch);
	void GerarRestricoes(int ch);
	void GerarVarBounds(int ch, double *lbd , double *ubd);
	void GerarCoefFuncaoObjetivo(int ch, double *cst);

	CMatrizEsparsa MatrizBalHid(int n_a, int ch);
	CMatrizEsparsa MatrizVmeta(int n_a, int ch);
	vetorfloat2 LimBalHid(int n_a, int ch);
	vetorfloat2 LimVmeta(int n_a, int ch);
	void MatrizRestricoesLineares(int n_a, CMatrizEsparsa &MM, int ch);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int ch);

public:
	Spcdec3SPHA(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	Spcdec3SPHA(CSistema * const sistema_end);
	~Spcdec3SPHA(void);

	int ResolverProblemaRL(Spcdec3Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int ch);
	int GetNCascatas() {return n_cascatas;}
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
	int GetNRestricoes(int comp) {return Nrestricoes[comp];}
	int GetNNZsRest(int comp);
	double GetLBvar(int comp, int var) {return vars[comp][var].get(GRB_DoubleAttr_LB);}
	//int GetNRestricoes(int comp) {return modelosGRB[comp]->get(GRB_IntAttr_NumConstrs);}
	//int GetNNZsRest(int comp) {return modelosGRB[comp]->get(GRB_IntAttr_NumNZs);}
	void GetBDesc( int comp , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd );
	int GetANZ( int comp);
	void GetADesc( int comp , int *Abeg , int *Aind , double *Aval );
	void SetPrecision( double precision );
};
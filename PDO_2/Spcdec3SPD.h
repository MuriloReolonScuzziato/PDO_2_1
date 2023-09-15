#pragma once

#include "Spcdec3Results.h"

#include "OPTUtils.h"

#include <time.h>

class Spcdec3SPD
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada nó)
	vetorGRBVar vars;

	CSistema * sistema_a;
	int n_nos, cont;
	vetorint n;
	double fo_n;
	vetorfloat x_n;

	// Atributos somente usados se o subproblema for Easy Component
	vetorint2 dtipo;
	vetorfloat2 dvalor;
	vetorint Mbeg, Mind;
	vetorfloat Mval;
	int Nrestricoes;

	void CriarVariaveis(int no);
	void CriarRestricoes(int no);
	void CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int no);
	void GerarRestricoes(int no);
	void GerarVarBounds(int no, double *lbd , double *ubd);
	void GerarCoefFuncaoObjetivo(int no, double *cst);

	CMatrizEsparsa MatrizRestDemanda(int n_a, int no);
	CMatrizEsparsa MatrizRestDemandaBarraUnica(int n_a, int no);
	CMatrizEsparsa MatrizLimFluxo(int n_a, int no);
	CMatrizEsparsa MatrizReserva(int n_a, int no);
	vetorfloat2 LimRestDemanda(int n_a, int no);
	vetorfloat2 LimRestDemandaBarraUnica(int n_a, int no);
	vetorfloat2 LimFluxo0(int n_a, int no);
	vetorfloat2 LimFluxo2(int n_a, int no);
	vetorfloat2 LimReserva(int n_a, int no);
	void MatrizRestricoesLineares(int n_a, CMatrizEsparsa &MM, int no);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int no);

public:
	Spcdec3SPD(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	Spcdec3SPD(CSistema * const sistema_end);
	~Spcdec3SPD(void);

	int ResolverProblemaRL(Spcdec3Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int no);
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
	int GetNRestricoes(int comp);
	int GetNNZsRest(int comp);
	void GetBDesc( int comp , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd );
	int GetANZ( int comp);
	void GetADesc( int comp , int *Abeg , int *Aind , double *Aval );
	void SetPrecision( double precision );
	double GetLBvar(int comp, int var) {return vars[comp][var].get(GRB_DoubleAttr_LB);}

};

#pragma once

#include "ResultadosConj.h"

class ConjSPD
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada nó)
	vetorGRBVar vars;

	CSistema * sistema_a;
	int n_nos;
	vetorint n;
	double fo_n;
	vetorfloat x_n;

	// Atributos somente usados se o subproblema for Easy Component
	vetorint2 Blinha, Bcoluna, Bnnz, dtipo;
	vetorfloat2 Bvalor, dvalor;
	vetorint Nrestricoes;

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
	void MatrizRestricoesLineares(int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int no, int &n_restricoes);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int no);

public:
	ConjSPD(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	ConjSPD(CSistema * const sistema_end);
	~ConjSPD(void);

	int ResolverProblemaRL(ResultadosConj * resultadoGurobi, const vetorfloat * const lambda, const int iter, int no);
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
	int GetNRestricoes(int comp) {return Nrestricoes[comp];}
	int GetNNZsRest(int comp) {int nnz = 0;for (size_t i = 1; i < Bnnz[comp].size(); i++) nnz += Bnnz[comp][i]; return nnz;}
	void GetBDesc( int comp , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd );
	int GetANZ( int comp);
	void GetADesc( int comp , int *Abeg , int *Aind , double *Aval );
};

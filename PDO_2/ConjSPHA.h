#pragma once

#include "ResultadosConj.h"

#include <algorithm>
using std::find;

class ConjSPHA
{
	vetorGRBModel modelosGRB;  // Modelos do Gurobi (um modelo para cada usina)
	vetorGRBVar vars;
	
	CSistema * sistema_a;
	int n_cascatas, N;
	vetorint n;
	double fo_ch;
	vetorfloat x_ch;
	vetorint2 cascata, cascata_jusante;		// vetor com os indices das hidros de cada cascata

	// Atributos somente usados se o subproblema for Easy Component
	vetorint2 Blinha, Bcoluna, Bnnz, dtipo;
	vetorfloat2 Bvalor, dvalor;
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
	void MatrizRestricoesLineares(int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int ch, int &n_restricoes);
	void MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int ch);

public:
	ConjSPHA(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	ConjSPHA(CSistema * const sistema_end);
	~ConjSPHA(void);

	int ResolverProblemaRL(ResultadosConj * resultadoGurobi, const vetorfloat * const lambda, const int iter, int ch);
	int GetNCascatas() {return n_cascatas;}
	int GetNVarPorComp(int comp) {return n[comp];}		// comp de cada tipo de subproblema
	int GetNRestricoes(int comp) {return Nrestricoes[comp];}
	int GetNNZsRest(int comp) {int nnz = 0;for (size_t i = 1; i < Bnnz[comp].size(); i++) nnz += Bnnz[comp][i]; return nnz;}
	//int GetNRestricoes(int comp) {return modelosGRB[comp]->get(GRB_IntAttr_NumConstrs);}
	//int GetNNZsRest(int comp) {return modelosGRB[comp]->get(GRB_IntAttr_NumNZs);}
	void GetBDesc( int comp , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd );
	int GetANZ( int comp);
	void GetADesc( int comp , int *Abeg , int *Aind , double *Aval );
};
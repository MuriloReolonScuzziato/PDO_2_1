#pragma once

#include "Sistema.h"

#include "MatrizEsparsa.h"

#include "gurobi_c++.h"
using namespace std;

typedef vector<string> vetorstring;
typedef vector<GRBModel *> vetorGRBModel;
typedef vector<GRBVar *> vetorGRBVar;
typedef vector<GRBConstr *> vetorGRBCons;

class ResultadosConj
{
	CSistema * sistema_a;
	int flag1, flag2, flag3, flag4, flag7;
	vetorfloat x, xa, obj_subp, x_med, lambda;
	int n, na, CH, R, N, I, status;
	vetorint2 cascata;			// vetor com os indices das hidros de cada cascata
	double obj, CFO, CIO, DEF, CTL, CTQ, VFOL, mipgap;
	bool Aggrgtd;
	vetorstring stringVariaveis;
	ofstream * inFile;

	vetorfloat2 ptr_x_til;			// Convexified primal solution for each component (numero de componentes x variáveis de cada componente)
	// estrutura dos vetores acima depende se o problema é agregado ou n
	// x_hat vai ser um dos vetores abaixo:
	vetorfloat2 ptr_x_subp;			// Primal solution for each subproblem (existe para modelo agregado e desagregado)
	vetorfloat2 ptr_x_subp_agg;		// Primal solution for the aggregated subproblem (1 x n + na)

public:
	ResultadosConj(CSistema * const sistema_end);
	~ResultadosConj(void);

	CMatrizEsparsa MatrizCalcFluxo(int n_a);
	CMatrizEsparsa MatrizCalcFluxo_cen(int n_a);

	void SetComponents( int * comp_inf, bool aggr );
	void GravarSolucao(double obj_a, vetorfloat x_a, int status_a, int comp);
	void GetSubGradiente(int wFi_a, double * SubG);
	void GetSubGradiente(double * SubG);
	double GetFobj();
	vetorfloat GetCompX(int comp);		// Retorna a solução de cada componente
	vetorfloat GetX_med();
	vetorfloat GetX_med(bool hat_or_til);
	vetorfloat GetX() {return x;}
	vetorfloat GetXa() {return xa;}

	double GetobjComp(int comp) {return obj_subp[comp];}
	int GetNA () {return na;}
	int Getn () {return n;}
	int GetCH () {return CH;}
	int GetR () {return R;}
	int GetI () {return I;}
	int GetN () {return N;}
	
	//vetorfloat * GetPtrXtil() { return ptr_x_til; }
	vetorfloat * GetPtrXtil() { return &ptr_x_til[0]; }
	vetorfloat * GetPtrXhat() { if(Aggrgtd) return &ptr_x_subp_agg[0];else return &ptr_x_subp[0];}
	void SetAgrModel( bool agg ) { Aggrgtd = agg; }

	void AlocarX_med();
	void AlocarPtrXHat();
	void AlocarXeXa();
	
	void ExportarX(string nome_arquivo);
	void ExportarXA(string nome_arquivo);
	void ExportarXmed(string nome_arquivo);

	void ZerarSolucoes();
};

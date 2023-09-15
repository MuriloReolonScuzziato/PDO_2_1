#pragma once

#include "Results.h"

class Spcdec3Results : public Results
{
	//vetorfloat x, xa, obj_subp, x_med, lambda;
	//int n, na, nd, CH, R, N, I, status;
	vetorint2 cascata;			// vetor com os indices das hidros de cada cascata
	//double obj, CFO, CIO, DEF, CTL, CTQ, VFOL, mipgap;
	bool Aggrgtd;
	//vetorstring stringVariaveis;
	//ofstream * inFile;
	//CMatrizEsparsa * X_hrst;		// solução da heuristica

	vetorfloat2 ptr_x_til;			// Convexified primal solution for each component (numero de componentes x variáveis de cada componente)
	// estrutura do vetor acima depende se o problema é agregado ou n
	// x_hat vai ser um dos vetores abaixo:
	vetorfloat2 ptr_x_subp;			// Primal solution for each subproblem (existe para modelo agregado e desagregado)
	vetorfloat2 ptr_x_subp_agg;		// Primal solution for the aggregated subproblem (1 x n + na)

public:
	Spcdec3Results(CSistema * const sistema_end);
	~Spcdec3Results(void);

	//CMatrizEsparsa MatrizCalcFluxo(int n_a);
	//CMatrizEsparsa MatrizCalcFluxo_cen(int n_a);

	void SetComponents( int * comp_inf, bool aggr );
	void GravarSolucao(double obj_a, vetorfloat &x_a, int status_a, int comp);
	void GetSubGradiente(int wFi_a, double * SubG);
	void GetSubGradiente(double * SubG);
	double GetFobj();
	vetorfloat GetCompX(int comp);		// Retorna a solução de cada componente
	vetorfloat GetX_med();
	vetorfloat GetX_med(bool til = true);
	vetorfloat GetX() {return x;}
	vetorfloat GetXa() {return xa;}
	//double GetLambda(int posicao) {return lambda[posicao];}

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
	void CalcularNormaSubgXtil(vetorfloat &norma1, vetorfloat &norma2, vetorfloat &normaInf);
	void CalcularNormaSubgXtil(vetorfloat &normas);

	//void EscreverArquivo(string nome_arquivo, double tempo_modelo, double tempo_resol, string saida);
	//void GravarSolHeuristica(double obj_a, CMatrizEsparsa * x_spr, double LB);
	//void GravarSolExtDE(GRBVar * vars, GRBConstr * constrs, int nvar, int nvar_a, int nconstr);
	//void GravarSolExtDE(int nvar, int nvar_a, int nconstr);
	void AlocarXmed(CMatrizEsparsa * x_spr);
	void AdicionarLBXtil(CMatrizEsparsa * x_spr, double mult);
	//void ArmazenarXtilhat();
	//void ImpromirArmazenarXtilhat();
};

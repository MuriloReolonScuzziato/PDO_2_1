#pragma once

#include "Sistema.h"

#include "MatrizEsparsa.h"

#include "gurobi_c++.h"
using namespace std;

typedef vector<string> vetorstring;
typedef vector<int> vetorint;
typedef vector<vetorint> vetorint2;
typedef vector<vetorfloat2> vetorfloat3;
//typedef vector<CMatrizEsparsa> vatorMatrizEsparsa;

typedef vector<GRBModel *> vetorGRBModel;		// deve ser um vetor de ponteiros??? como apagar?
typedef vector<GRBVar *> vetorGRBVar;
typedef vector<GRBConstr *> vetorGRBCons;

class ResultadosConj
{
	CSistema * sistema_a;
	int flag1, flag2, flag3, flag4, flag7, decomp_str;
	vetorfloat x, xa, obj_subp, x_med, lambda;
	int n, na, CH, R, N, I, status;
	vetorint2 cascata;			// vetor com os indices das hidros de cada cascata
	double obj, CFO, CIO, DEF, CTL, CTQ, VFOL, mipgap;
	bool Aggrgtd;
	vetorstring stringVariaveis;
	ofstream * inFile;
	//vetorfloat * ptr_x_hat;		// Primal solution for each component (depende do numero de componentes)
	//vetorfloat * ptr_x_til;		// Convexified primal solution for each component (depende do numero de componentes)
	//vetorfloat * ptr_x_subp;		// Primal solution for each subproblem ( = ptr_x_hat in the disaggregated model)
	

	//??? diferença entre a estrutura de dados acima e abaixo, qual usar?


	//vetorfloat2 * ptr_x_hat;		// ponteiro q aponta para ptr_x_subp ou ptr_x_subp_agg dependendo do modelo desag. ou agreg.
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
	void CarregarResultadosED(double, vetorfloat, vetorfloat, int, double);
	void EscreverArquivo(string, double, double);
	void CarregarResultadosDecEsp(double obj_a, vetorfloat x_a, int status_a, int tipoSubproblema, int id1, int id2);
	void GravarSolucao(double obj_a, vetorfloat x_a, int status_a, int comp);
	void AlocarPtrXHat();
	void AlocarXeXa();
	void AlocarX_med();

	vetorfloat GetX() {return x;}
	vetorfloat GetXa() {return xa;}
	vetorfloat GetX_med();
	vetorfloat GetX_med(bool hat_or_til);
	vetorfloat GetCompX(int comp);		// Retorna a solução de cada componente
	double GetFobj();
	double GetobjComp(int comp) {return obj_subp[comp];}
	int GetNA () {return na;}
	int Getn () {return n;}
	int GetCH () {return CH;}
	int GetR () {return R;}
	int GetI () {return I;}
	int GetN () {return N;}
	void GetSubGradiente(int comp, double * SubG);
	void GetSubGradiente(double * SubG);
	void SetDecompositionStrategy( int str ) { decomp_str = str; }
	void SetComponents( int * comp_inf, bool aggr );
	//vetorfloat * GetPtrXtil() { return ptr_x_til; }
	vetorfloat * GetPtrXtil() { return &ptr_x_til[0]; }
	vetorfloat * GetPtrXhat() { if(Aggrgtd) return &ptr_x_subp_agg[0];else return &ptr_x_subp[0];}
	void SetAgrModel( bool agg ) { Aggrgtd = agg; }

	void GetSubGradiente2(int wFi_a, double * SubG);
	void GetSubGradiente2(double * SubG);
	
	void ExportarX(string nome_arquivo);
	void ExportarXA(string nome_arquivo);
	void ExportarXmed(string nome_arquivo);
	void ZerarSolucoes();
};

// Guarda ultimos resultados dos compontentes do problema dual... 
// tb possui funçoes para calcular subgradiente...
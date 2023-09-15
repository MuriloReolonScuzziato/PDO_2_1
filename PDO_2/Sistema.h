#pragma once
#include <fstream>
using std::ifstream;
using std::ios;
using std::ofstream;

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;
using std::setw;
using std::left;
using std::internal;
using std::right;

#include <eigen>
#include <Sparse>
using namespace Eigen;

#include <math.h>

#define CNM 5
#define CNJ 5
#define CRH 6
#define COT 3
#define CPT 3

#include "Linha.h"

#include "MatrizEsparsa.h"

typedef vector<CBarra> vetorbarra;
typedef vector<CLinha> vetorlinha;
typedef vector<CHidreletrica> vetorhidreletrica;
typedef vector<CTermeletrica> vetortermeletrica;
typedef vector<CDemanda> vetordemanda;

#define FPH_APP 0		// método de aproximação dos planos da função de produção, 0 = Taylor de primeira ordem, 1 = mínimos quadrados
#define FPH_REM 1		// remover cortes que anulam outros (em q), 0 = não remover planos; 1 = remover planos que a derivada em 'q' diminui

class CSistema
{
	string diretorio;
	string arq_hidreletricas;
	string arq_unidades;
	string arq_termeletricas;
	string arq_barras;
	string arq_linhas;
	string arq_parametros;
	string arq_demanda;
	string arq_afluencia;
	string arq_cond_inic_h;
	string arq_cond_inic_t;

	int Tt1, Tt2, fator_reserva, barra_ref, n_cenarios;
	double delta_t1, delta_t2, custo_def, custo_vfol;
	vetorfloat precond;			// vetor de valores dos preconditioners (1 x na), mesmo tamanho que o vetor de multiplicadores
	vetorfloat res_gir;

	int flag_modelo_rede;		// 0: barra unica; 1: rede completa; 2: modelo compacto; 3: modelo compacto com lazy constraints
	bool flag_vfol;				// 1 se for considerar a variável de vfol
	bool flag_phmax;			// 1 se for considerar no ED a variável de phmax
	int flag_init_custoT;		// Numero de partes para linearizar a função de custo de produção das termicas (0:quadratica, 1:uma reta, 2:duas retas,...)
	bool flag_var_bin;			// 1 se for considerar variáveis binárias e 0 se forem todas contínuas
	int flag_appr_custo;		// 0 use least-square approximantion, 1 use perspective-cut formulation, 2 perspective-cut formulation dynamically (entr4 is the initial number of cuts (no min. = 2), and entr6 is the maximum number of cuts that can be added)
	int flag_Tbinary_model;		// the binary variables of the thermal units: 0 old model; 1 modern model.
	int flag_preconditioner;	// Kind of Preconditioner: 0 none of them; 1 probability of the node; 2 square root of the probability of the node; 3 maximum value of the relaxed constraint; 4 combinaion between 1 and 3; 5 combination between 2 and 3
	
	void CriarHidreletricas();
	void CriarGrupos();
	void GerarFpgh();
	void AproximarFCT();
	void CriarTermeletricas();
	void GerarPhg(double, double, double, double);	// Criar a função de produção hidreletrica do grupo
	vetorfloat2 GerarCoeficientes(CHidreletrica * hidreletricaPtr_a, CUnidades * grupoVtr_a);	// Gerar coneficientes da função de produção hidreletrica do grupo
	vetorfloat GerarCoeficientes2(CHidreletrica * hidreletricaPtr_a, CUnidades * grupoVtr_a, CSistema * sistema_a);	// 
	void CarregarParametros();
	void CriarDemandas();
	void CriarBarras();
	void CriarLinhas();
	void CarregarAfluencias();
	void CarregarCIH();
	void CarregarCIT();
	void CalcularReserva();
	void CalcularMatrizBeta();
	void CalcularMatrizBetaS();

public:
	vetorhidreletrica hidreletricasVtr;
	vetortermeletrica termeletricasVtr;
	vetordemanda demandasVtr;
	vetorbarra barrasVtr;
	vetorlinha linhasVtr;
	MatrixXd Beta;
	CMatrizEsparsa BetaS;
	
	CSistema(string, string, string, string, string, string, string, string, string, string, string);
	~CSistema(void);
	int GetTt1() {return Tt1;}
	int GetTt2() {return Tt2;}
	int GetFatorReserva() {return fator_reserva;}
	double GetDeltaT1() {return delta_t1;}
	double GetDeltaT2() {return delta_t2;}
	int GetBarraRef() {return barra_ref;}
	int GetNCenarios() {return n_cenarios;}
	double GetCustoDeficit() {return custo_def;}
	double GetCustoVfol() {return custo_vfol;}
	double GetReserva(int t) {return res_gir[t];}
	void SetFlags(int mod_rede, bool vfol, bool phmax, int init_custoT, bool var_bin, int appr_custo, int Tbinary_model, int preconditioner);
	void SetFlagVarBin(bool var_bin) {flag_var_bin = var_bin;}
	// substituir nas demais classes e remover 2 funções abaixo
	//void SetFlagBarraUnica(bool bar_uni) {cout << "funcao substituida!" << endl;flag_modelo_rede = bar_uni;}
	//bool GetFlagBarraUnica() {cout << "funcao substituida!" << endl; return flag_modelo_rede;}
	//
	/// 0: barra unica; 1: rede completa; 2: modelo compacto; 3: modelo compacto com lazy constraints
	void SetFlagModeloRede(int mod_rede) {flag_modelo_rede = mod_rede;}		
	int GetFlagModeloRede() {return flag_modelo_rede;}
	bool GetFlagVfol() {return flag_vfol;}
	bool GetFlagPhmax() {return flag_phmax;}
	int GetFlagInitAproxCT() {return flag_init_custoT;}
	bool GetFlagVarBin() {return flag_var_bin;}
	int GetFlagMaxAproxCT() {return flag_appr_custo;}
	int GetFlagTbinaryModel() {return flag_Tbinary_model;}
	int GetFlagPrecond() {return flag_preconditioner;}
	double GetPrecondidioner(int indice) {return precond[indice];}
	void SetPreconditioners(vetorfloat precond_a);
	void DeterminarLimitesV();
};
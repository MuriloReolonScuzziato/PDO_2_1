/*--------------------------------------------------------------------------*/
/*---------------------------- File Hrstc.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para aplicar heuristicas na solução da RL (considerando a decom-
 * posiçao SpcDec_3
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#pragma once

#define DEBUG 0

//#include "Spcdec3Results.h"
#include "Results.h"

#include <algorithm>
using std::sort;
using std::reverse;

#include "OPTUtils.h"

typedef vector<CMatrizEsparsa> vetorMtrzEsparsa;

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Hrstc --------------------------*/
/*--------------------------------------------------------------------------*/

class Hrstc 
{
	// Parametros ( carregados do arquivo)
	// parametros de ponderaçao da funçao objetivo
	double ALFA;		// fator de ponderaçao entre a funçao objetivo original e os termos proximais. = [0,1]
	double BETA;		// fator de ponderaçao entre os termos proximais (soluçao da RL e convexificada). = [0,1]
	double GAMA;		// fator de ponderaçao dos termos proximais entre variáveis continnuas e binárias. = [0,1]
	int NORM_P;			// Norma 1 (linear); 2 (quadratica) ou 0 norma infinita (não implementada ainda)
	int APPROACH; 		// 0 without benders, 1 with benders
	int INITCUT;		// Usar cortes do modelo Estendido (continuo, linear) na iteração 0 da heuristica
	int itBenders;		// numero de iterações da heuristica a cada iteração do bundle
	int itBendersF;		// numero de iterações da heuristica ao final do bundle
	// Numero de períodos que compõem cada problema da Heuristica
	int WIN_SIZE;		// Numero de periodos que tem a soluçao salva (deve ser proporcional ao numero de períodos do problema)
	int TIME_R;			// means that I give to each subproblem (  TIME_R x total running time ) / ( number of subproblems )
	// ajustes ( constructor )
	int WIN_FUT_SIZE;	// Numero de periodos utilizados para ter se uma ideia do futuro (essa janela perde o sentido se forem utilizados cortes de Benders)
	// natureza das variáveis "binárias"
	int	BIN_MW; 		// Define as var. binárias p/ a janela principal: 0 = fixed in x_hat, 1 = free and continuous, 2 = free and binary
	int	BIN_FW;			// Define as var. binárias p/ a janela futuro: 0 = fixed in x_hat, 1 = free and continuous, 2 = free and binary

	CSistema * sistema_a;
	//Spcdec3Results * resultadosGurobi;
	Results * resultadosGurobi;
	ofstream * log_auxiliar;

	GRBEnv ambGRB;				// Ambiente Gurobi
	vetorGRBModel modelosGRB;   // Modelos do Gurobi (um modelo para cada problema-janela)
	vetorGRBVar vars;
	vetorGRBCons constr;
	// PL
	GRBModel modeloPL;
	GRBVar * varPL;

	double alfa, beta, gama;			// parametros de ponderaçao da funçao objetivo
	double Cvfol, PenVF;
	int n_modelos;
	int n, flag1, flag2, flag3, flag4, flag7, T;
	int flagH1, flagH2, flagH3;			// Flags para as heuristicas
	int ws, wlas, nw;					// ws largura da janela em periodos, wlas largura da janela-futuro em periodos, nw numero de janelas(ou problemas)
	vetorint nStatus;					// vetor com status de resoluçao de cada subproblema
	int Status;							// 0 se o resultado for ótimo ou por tempo, e soma-se 1 para subproblema que for inviável;
	vetorfloat fo_subp, ff_subp, foRL;	// vetor com a funçao objetivo de cada subproblema (considerando somente janela principal) e vetor com os valores da aprox. do custo futuro
	double fo, fo_ub_a;					// funçao objetivo total
	//vetorfloat fo_media;				// media movel da função objetivo
	// O problema de cada janela (conjunto de periodos) é resolvido considerando ws e wlas, porem o resultado gravada é só o de ws
	CMatrizEsparsa * X;					// x em formato esparso
	// passar desse formato para x apos resolver td o problema!!!
	CMatrizEsparsa * A;					// partes da matriz de restriçoes (com acoplamento temporal: Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida)
	vetorfloat3 lim_iniciais;			// vetor com os limites das iniciais das restriçoes, similiar à matriz A porém com os RHS
	vetorint2 posicao_rest_din;			// posiçao inicial das restriçoes que sao alteradas para cada subproblema, relacionado com a matriz A
	vetorint posicao_rest_din_Ext;		// posiçao inicial das restriçoes que sao alteradas para cada subproblema, correspondente ao modelo Ext
	vetorfloat x, x_hat, x_til;			// x = [pt u cp F teta ph v d s phmax phg q z def], x_hat out of LR solution, x_til convexified solution
	// Considerando que x é o valor médio das variáveis que foram duplicadas, x = (xo + xa) / 2
	vetorint2 vtr_nos_pri;				// vetor com os numeros dos nos principais de cada subproblema
	vetorint2 vtr_nos_fut;				// vetor com os numeros dos nos futuro de cada subproblema
	vetorint numero_var_subp, numero_var_subp_principal;	// numero de variáveis (originais) do modelo
	vetorint numero_var_ad, numero_var_ad2;					// numero de variáveis adicionais para o termo proximal linear, x_til e x_hat
	vetorint restr_mod;					// numero de restrições (originais) por modelo
	vetorint nvar_acum;				// numero de var. acumuladas por subproblema
	vetorint nrest_orig;			// numero de restrições de cada subproblema (sem considerar os cortes)
	vetorint subpB1etapa, subpB2etapa, subpB3etapa;	// identificador dos subproblemas de cada etapa do Backward
	vetorfloat3 acomplamentos;		// vetorint2 para cada modelo com informações de seus acoplamentos
	//CMatrizEsparsa PI;				// coeficiente acumulado das variáveis nos cortes: 1 x n
	vetorMtrzEsparsa L;				// coeficiente das variáveis nos cortes: n_iterações x n, uma matriz para cada subproblema!!
	int itB;						// numero de iterações totais Backward e Forward
	vetorint2 index_var_bin;		// indice das var. binarias por subproblema
	vetorfloat lhs_var;				// lhs variavel dos cortes
	vetorfloat2 rhs_orig;			// rhs original dos cortes: itB x n_subproblemas, um valor para cada subproblema (usado para atualizar os cortes ao final do forward)
	vetorint n_var_folga;			// variáveis de folga para cada subproblema (para usar cortes de otimalidade ao invés de cortes de viabilidade)
	//vetorint ind_restr_folga;		// indice inicial das restrições com var. de folga mesmo para todos os subproblemas (já está gravado no vetor posicao_rest_din)
	//vetorfloat2 rhs_orig_bin;		// rhs original dos cortes considerando as var. bin. como complicadas (backward): itB x n_subproblemas, um valor para cada subproblema (usado para atualizar os cortes ao final do forward)
	//vetorMtrzEsparsa L_bin;			// coeficiente das variáveis bin. nos cortes: n_iterações x n, uma matriz para cada subproblema!!
	vetorint2 ind_var_fol;			// indice das variáveis de folga adicionadas em cada modelo
	vetorint ident_var;				// identificador de variável que é considerada como termo proximal
	vetorstring stringVariaveis;
	double coeftol;					// tolerance for the coefficients of the cuts 

	vetorfloat debug_alfa_B2;

	CMatrizEsparsa MatrizRestDemanda();
	CMatrizEsparsa MatrizRestDemandaBarraUnica();
	CMatrizEsparsa MatrizLimFluxo();
	CMatrizEsparsa MatrizLimPhgMin();
	CMatrizEsparsa MatrizLimPhgMax();
	CMatrizEsparsa MatrizLimqMin();
	CMatrizEsparsa MatrizLimqMax();
	CMatrizEsparsa MatrizBalHid();
	CMatrizEsparsa MatrizTup();
	CMatrizEsparsa MatrizTdown();
	CMatrizEsparsa MatrizTupDown();
	CMatrizEsparsa MatrizRampaUp();
	CMatrizEsparsa MatrizRampaDown();
	CMatrizEsparsa MatrizLimPtMax();
	CMatrizEsparsa MatrizLimPtMin();
	CMatrizEsparsa MatrizRestCP();
	CMatrizEsparsa MatrizVmeta();
	CMatrizEsparsa MatrizBalPotencia();
	CMatrizEsparsa MatrizBalVazao();
	CMatrizEsparsa MatrizFuncProd();
	CMatrizEsparsa MatrizPhMax();
	CMatrizEsparsa MatrizReserva();
	CMatrizEsparsa MatrizCortesF();
	vetorfloat2 LimRestDemanda();
	vetorfloat2 LimRestDemandaBarraUnica();
	vetorfloat2 LimFluxo0();
	vetorfloat2 LimFluxo2();
	vetorfloat2 LimPhgMin();
	vetorfloat2 LimPhgMax();
	vetorfloat2 LimQMin();
	vetorfloat2 LimQMax();
	vetorfloat2 LimBalHid();
	vetorfloat2 LimTup();
	vetorfloat2 LimTdown();
	vetorfloat2 LimTupDown();
	vetorfloat2 LimRampaUp();
	vetorfloat2 LimRampaDown();
	vetorfloat2 LimRestCP();
	vetorfloat2 LimPtMax();
	vetorfloat2 LimPtMin();
	vetorfloat2 LimVmeta();
	vetorfloat2 LimBalPotenciaL();
	vetorfloat2 LimBalVazaoL();
	vetorfloat2 LimFuncProdL();
	vetorfloat2 LimPhMax();
	vetorfloat2 LimReserva();
	vetorfloat2 LimCortesF();
	void MatrizRestricoesLineares(CMatrizEsparsa  * MM, CMatrizEsparsa &MMpl);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor, vetorint * LimTipo_pl, vetorfloat * LimValor_pl);

	void CriarModelos(double max_time);
	
	void CriarVariaveis(int modelo);
	void CriarVariaveisPL();
	
	void CriarRestricoes();
	void AlocarRestricoes(CMatrizEsparsa * matriz, CMatrizEsparsa * MM);
	void AlocarRestricoesVmeta(CMatrizEsparsa * matriz, CMatrizEsparsa * MM);
	void AlocarLimitesSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor);
	void AlocarLimitesVmetaSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor);
	void CriarRestriçõesProximais(int modelo);
	void AtualizarRHSproximais(int modelo);

	void CriarFuncaoObjetivo(int modelo, bool backward);
	void CriarFuncaoObjetivoQ(int modelo, bool backward);
	void CriarFuncaoObjetivoL(int modelo, bool backward);
	void CriarFuncaoObjetivoPL();
	
	void FixarCondIniciais(int modelo);
	void FixarCondIniciaisPL();
	void AtualizarVariaveis(int janela);
	void AtualizarRestricoes(int modelo, CMatrizEsparsa * x_spr);
	
	double ResolverPL(CMatrizEsparsa * x_spr);
	double ResolverSubpLinear(int modelo, CMatrizEsparsa * x_spr);
	void CalcularFuncaoObjetivoRL(vetorfloat & fo);
	void AlocarSolucao(int modelo, CMatrizEsparsa * x_spr);

	int ResolverForward(CMatrizEsparsa * x_spr, double * fo_a);
	int ResolverBackward(CMatrizEsparsa * x_spr);
	int ResolverBackward0(CMatrizEsparsa * x_spr);
	int ResolverSubp(int modelo, bool backward = false);
	
	void IdentAcoplamentos(int modelo);
	void SubpMILP2LP(CMatrizEsparsa * x_spr, int modelo, vetorint * limites);
	void SubpMILP2LPFixing(CMatrizEsparsa * x_spr, int modelo, vetorint * limites);
	void SubpLP2MILP(CMatrizEsparsa * x_spr, int modelo, vetorint * limites);
	void DeterminarIndiceVarBin();
	void AtualizarCortes(CMatrizEsparsa & x_spr, int modelo);
	void AtualizarCortesBin(CMatrizEsparsa * x_spr, int modelo);
	double CalcularFuncaoObjetivo(int modelo, CMatrizEsparsa * x_spr);
	void RemoverVarFol(int modelo);
	void AdicionarVarFol(int modelo);
	void SelecionarVarProx();
	void ImprimirSolSub(int modelo);
	void ImprimirSolSub(int modelo, int nvar);

public:
	Hrstc(CSistema * const sistema_end, Results * const resultadosGurobi_end, double max_time, istream *iStrm = NULL);
	Hrstc(CSistema * const sistema_end, Results * const resultadosGurobi_end, int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a, double max_time);		// Usa parametros da entrada
	~Hrstc(void);

	double ResolverHeuristica();
	double ResolverHeuristica(bool so_heuristica);
	void EscreverX();
	int GetStatus() {return Status;}
	double GetFO() {return fo;}
	CMatrizEsparsa * GetX() {return X;}
	void SetLog(ofstream * log);
	int GetInitCut() {return INITCUT;}
	int GetBenders() {return APPROACH;}
	//double ResolverHeuristica(vetorfloat &x_a);		// gravar soluçao em x externo... ver o que é gravado e apagado a cada chamada da funçao ResolverHeuristica (tempo de processamento!!)
};
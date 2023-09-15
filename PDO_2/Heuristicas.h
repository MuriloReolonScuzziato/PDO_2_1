/*--------------------------------------------------------------------------*/
/*---------------------------- File Heuristicas.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para aplicar heuristicas na solução da RL (considerando a decom-
 * posiçao SpcDec_2
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

// Heuristic strategies
// => Original objective function
// => Original o.f. + quadratic proximal term
// => Original o.f. + linear proximal term (Gurobi and FP)

#define APPROACH 0			// 0 forward approach, 1 backward approach and 2 mixed approach

// Numero de períodos que compõem cada problema da Heuristica
#define	WIN_SIZE 1			// Numero de periodos que tem a soluçao salva (deve ser proporcional ao numero de períodos do problema)
#define WIN_FUT_SIZE 0		// Numero de periodos utilizados para ter se uma ideia do futuro
// parametros de ponderaçao da funçao objetivo
#define ALFA 0			// fator de ponderaçao entre a funçao objetivo original e os termos proximais. = [0,1]
#define BETA 1			// fator de ponderaçao entre os termos proximais (soluçao da RL e convexificada). = [0,1]
#define GAMA 0.5			// fator de ponderaçao dos termos proximais entre variáveis continnuas e binárias. = [0,1]
// natureza das variáveis "binárias"
#define	BIN_MW 2			// Define as var. binárias p/ a janela principal: 0 = fixed in x_hat, 1 = free and continuous, 2 = free and binary
#define	BIN_FW 2			// Define as var. binárias p/ a janela futuro: 0 = fixed in x_hat, 1 = free and continuous, 2 = free and binary
// norma para os termos proximais (só aplicada para os termos proximais das var. binárias)
#define NORM_P 2			// Norma 1 (linear) ou 2 (quadratica)

#define TIME_R 4			// means that I give to each subproblem (  TIME_R x total running time ) / ( number of subproblems )

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#pragma once

#include "ResultadosConj.h"

#include <algorithm>
using std::sort;
using std::reverse;

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Heuristicas --------------------------*/
/*--------------------------------------------------------------------------*/

class Heuristicas
{
	CSistema * sistema_a;
	ResultadosConj * resultadosGurobi;

	GRBEnv ambGRB;			// Ambiente Gurobi
	vetorGRBModel modelosGRB;   // Modelos do Gurobi (um modelo para cada problema-janela)
	vetorGRBVar vars;
	vetorGRBCons constr;
	double alfa, beta, gama;			// parametros de ponderaçao da funçao objetivo
	int n_modelos;
	int n, flag1, flag2, flag3, flag4, T;
	int flagH1, flagH2, flagH3;			// Flags para as heuristicas
	int ws, wlas, nw;					// ws largura da janela em periodos, wlas largura da janela-futuro em periodos, nw numero de janelas(ou problemas)
	vetorint nStatus;					// vetor com status de resoluçao de cada subproblema
	int Status;							// 0 se o resultado for ótimo ou por tempo, e soma-se 1 para subproblema que for inviável;
	vetorfloat fo_subp;				// vetor com a funçao objetivo de cada subproblema (considerando somente janela principal)
	double fo;							// funçao objetivo total
	// O problema de cada janela (conjunto de periodos) é resolvido considerando ws e wlas, porem o resultado gravada é só o de ws
	CMatrizEsparsa * X;					// x em formato esparso
	// passar desse formato para x apos resolver td o problema!!!
	CMatrizEsparsa * A;					// partes da matriz de restriçoes (com acoplamento temporal: Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida)
	vetorfloat3 lim_iniciais;			// vetor com os limites das iniciais das restriçoes, similiar à matriz A porém com os RHS
	vetorint2 posicao_rest_din;			// posiçao inicial das restriçoes que sao alteradas para cada subproblema, relacionado com a matriz A
	vetorfloat x, x_hat, x_til;			// x = [pt u cp F teta ph v d s phmax phg q z def], x_hat out of LR solution, x_til convexified solution
	// Considerando que x é o valor médio das variáveis que foram duplicadas, x = (xo + xa) / 2
	vetorint2 vtr_nos_pri;				// vetor com os numeros dos nos principais de cada subproblema
	vetorint2 vtr_nos_fut;				// vetor com os numeros dos nos futuro de cada subproblema
	vetorint numero_var_subp, numero_var_subp_principal;

	CMatrizEsparsa MatrizRestDemanda();
	CMatrizEsparsa MatrizRestDemandaBarraUnica();
	CMatrizEsparsa MatrizLimFluxo();
	CMatrizEsparsa MatrizLimPhgMin();
	CMatrizEsparsa MatrizLimPhgMax();
	CMatrizEsparsa MatrizBalHid();
	CMatrizEsparsa MatrizLimqMin();
	CMatrizEsparsa MatrizLimqMax();
	CMatrizEsparsa MatrizTup();
	CMatrizEsparsa MatrizTdown();
	CMatrizEsparsa MatrizRampaUp();
	CMatrizEsparsa MatrizRampaDown();
	CMatrizEsparsa MatrizLimPtMin();
	CMatrizEsparsa MatrizLimPtMax();
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
	vetorfloat2 LimBalHid();
	vetorfloat2 LimQMin();
	vetorfloat2 LimQMax();
	vetorfloat2 LimTup();
	vetorfloat2 LimTdown();
	vetorfloat2 LimRampaUp();
	vetorfloat2 LimRampaDown();
	vetorfloat2 LimPtMin();
	vetorfloat2 LimPtMax();
	vetorfloat2 LimRestCP();
	vetorfloat2 LimVmeta();
	vetorfloat2 LimBalPotenciaL();
	vetorfloat2 LimBalVazaoL();
	vetorfloat2 LimFuncProdL();
	vetorfloat2 LimPhMax();
	vetorfloat2 LimReserva();
	vetorfloat2 LimCortesF();
	void MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes, int n_modelos);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);

	void CriarModelos(double max_time);
	void CriarVariaveis(int modelo);
	void CriarRestricoes();
	void AlocarRestricoes(CMatrizEsparsa * matrizcompleta, vetorint * indexLinha, vetorint * indexColuna, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count, int n_modelos);
	void AlocarRestricoesVmeta(CMatrizEsparsa * matrizcompleta, vetorint * indexLinha, vetorint * indexColuna, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count, int n_modelos);
	void AlocarLimitesSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor);
	void AlocarLimitesVmetaSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor);
	void CriarFuncaoObjetivo(int modelo);
	void CriarFuncaoObjetivoQ(int modelo);
	void CriarFuncaoObjetivoL(int modelo);


	void AtualizarVariaveis(int janela);
	void AtualizarRestricoes(int modelo, CMatrizEsparsa * x_spr);
	double CalcularFuncaoObjetivo(int modelo, CMatrizEsparsa * x_spr);
	double CalcularFuncaoObjetivoRL(int modelo);
	void ResolverSubps(CMatrizEsparsa * x_spr, double * fo_a);
	void AlocarSolucao(int modelo, CMatrizEsparsa * x_spr);
	void AlterarEstadoVarBin(int modelo, int posicao, int estado);
	void ResolverInviabilidadeHidro(int modelo);


	void ConferirDeficit(int modelo);
	void ResolverInviabilidadeTermo(int modelo);
	vetorint ListaPrioridades(int modelo, int no, int tipo);

public:
	Heuristicas(CSistema * const sistema_end, ResultadosConj * const resultadosGurobi_end, double max_time);
	Heuristicas(CSistema * const sistema_end, ResultadosConj * const resultadosGurobi_end, int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a, double max_time);		// Usa parametros da entrada
	~Heuristicas(void);

	double ResolverHeuristica();
	double ResolverHeuristica(vetorfloat x_a, vetorfloat x_a2);
	int GetStatus() {return Status;}
	double GetFO() {return fo;}
	//double ResolverHeuristica(vetorfloat &x_a);		// gravar soluçao em x externo... ver o que é gravado e apagado a cada chamada da funçao ResolverHeuristica (tempo de processamento!!)
};
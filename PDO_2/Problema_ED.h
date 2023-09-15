#pragma once

#include "Resultados_ED.h"

#include "MatrizEsparsa.h"

#include "CallbackED.h"

#include "gurobi_c++.h"
using namespace std;

class CProblema_ED
{
	GRBEnv   ambGRB;     // Ambiente Gurobi
	GRBModel modeloGRB;  // Modelo Gurobi
	GRBVar * vars;
	size_t n;
	int flag1, flag2, flag3, flag4, flag7, T;
	vetorfloat x, L;
	double fo, Cvfol;
	CSistema * sistema_a;
	Resultados_ED * resultados;
	//MatrixXd Beta;

	void CriarVariaveis();
	void CriarFuncaoObjetivo();
	void CriarRestricoes();
	void CriarFuncaoObjetivo2();
	void FixarCondIniciais();

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
	CMatrizEsparsa MatrizTupDown();
	CMatrizEsparsa MatrizRampaUp();
	CMatrizEsparsa MatrizRampaDown();
	CMatrizEsparsa MatrizLimPtMin();
	CMatrizEsparsa MatrizLimPtMax();
	CMatrizEsparsa MatrizLimPtMinMax();
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
	vetorfloat2 LimTupDown();
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
	void MatrizRestricoesLineares(CMatrizEsparsa &matriz);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);

public:
	CProblema_ED(CSistema * const sistema_ptr);
	~CProblema_ED(void);
	int ResolverProblema();
	void EscreverResultados(string nome_arquivo, double tempo_modelo, double tempo_resol, string saida) {resultados->EscreverArquivo(nome_arquivo, tempo_modelo, tempo_resol, saida);}
	void ConferirCoerenciaDados();
};

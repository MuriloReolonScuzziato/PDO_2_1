#pragma once

#include "ResultadosConj.h"

#include "MatrizEsparsa.h"

#include "CallbackED.h"

#include "gurobi_c++.h"
using namespace std;

typedef vector<int> vetorint;

namespace met_ED
{

class CProblema_ED
{
	GRBEnv   ambGRB;     // Ambiente Gurobi
	GRBModel modeloGRB;  // Modelo Gurobi
	GRBVar * vars;
	int n, flag1, flag2, flag3, flag4, flag7, T;
	vetorfloat x, L;
	double fo;
	CSistema * sistema_a;

	void CriarVariaveis();
	void CriarFuncaoObjetivo();
	void CriarRestricoes();
	void CriarFuncaoObjetivo2();
	void FixarCondIniciais();

	CMatrizEsparsa CProblema_ED::MatrizRestDemanda();
	CMatrizEsparsa CProblema_ED::MatrizRestDemandaBarraUnica();
	CMatrizEsparsa CProblema_ED::MatrizLimFluxo();
	CMatrizEsparsa CProblema_ED::MatrizLimPhgMin();
	CMatrizEsparsa CProblema_ED::MatrizLimPhgMax();
	CMatrizEsparsa CProblema_ED::MatrizBalHid();
	CMatrizEsparsa CProblema_ED::MatrizLimqMin();
	CMatrizEsparsa CProblema_ED::MatrizLimqMax();
	CMatrizEsparsa CProblema_ED::MatrizTup();
	CMatrizEsparsa CProblema_ED::MatrizTdown();
	CMatrizEsparsa CProblema_ED::MatrizTupDown();
	CMatrizEsparsa CProblema_ED::MatrizRampaUp();
	CMatrizEsparsa CProblema_ED::MatrizRampaDown();
	CMatrizEsparsa CProblema_ED::MatrizLimPtMin();
	CMatrizEsparsa CProblema_ED::MatrizLimPtMax();
	CMatrizEsparsa CProblema_ED::MatrizLimPtMinMax();
	CMatrizEsparsa CProblema_ED::MatrizRestCP();
	CMatrizEsparsa CProblema_ED::MatrizVmeta();
	CMatrizEsparsa CProblema_ED::MatrizBalPotencia();
	CMatrizEsparsa CProblema_ED::MatrizBalVazao();
	CMatrizEsparsa CProblema_ED::MatrizFuncProd();
	CMatrizEsparsa CProblema_ED::MatrizPhMax();
	CMatrizEsparsa CProblema_ED::MatrizReserva();
	CMatrizEsparsa CProblema_ED::MatrizCortesF();
	vetorfloat2 CProblema_ED::LimRestDemanda();
	vetorfloat2 CProblema_ED::LimRestDemandaBarraUnica();
	vetorfloat2 CProblema_ED::LimFluxo0();
	vetorfloat2 CProblema_ED::LimFluxo2();
	vetorfloat2 CProblema_ED::LimPhgMin();
	vetorfloat2 CProblema_ED::LimPhgMax();
	vetorfloat2 CProblema_ED::LimBalHid();
	vetorfloat2 CProblema_ED::LimQMin();
	vetorfloat2 CProblema_ED::LimQMax();
	vetorfloat2 CProblema_ED::LimTup();
	vetorfloat2 CProblema_ED::LimTdown();
	vetorfloat2 CProblema_ED::LimTupDown();
	vetorfloat2 CProblema_ED::LimRampaUp();
	vetorfloat2 CProblema_ED::LimRampaDown();
	vetorfloat2 CProblema_ED::LimPtMin();
	vetorfloat2 CProblema_ED::LimPtMax();
	vetorfloat2 CProblema_ED::LimRestCP();
	vetorfloat2 CProblema_ED::LimVmeta();
	vetorfloat2 CProblema_ED::LimBalPotenciaL();
	vetorfloat2 CProblema_ED::LimBalVazaoL();
	vetorfloat2 CProblema_ED::LimFuncProdL();
	vetorfloat2 CProblema_ED::LimPhMax();
	vetorfloat2 CProblema_ED::LimReserva();
	vetorfloat2 CProblema_ED::LimCortesF();
	void CProblema_ED::MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes);
	void CProblema_ED::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);

public:
	CProblema_ED(CSistema * const);
	~CProblema_ED(void);
	int ResolverProblema(ResultadosConj * resultadoGurobi);
};

};
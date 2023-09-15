#pragma once

//#include "Spcdec2Results.h"
#include "Results.h"

#include "MatrizEsparsa.h"

#include "gurobi_c++.h"
using namespace std;

#include "OPTtypes.h"
using namespace OPTtypes_di_unipi_it;

typedef vector<int> vetorint;

class Spcdec2ExtDE
{
	GRBEnv   ambGRB;     // Ambiente Gurobi
	GRBModel modeloGRB;  // Modelo Gurobi
	GRBVar * vars;
	GRBConstr * constr;
	int n, no, na, nd, flag1, flag2, flag3, flag4, flag7, T, rest_o, rest_tot;
	vetorfloat x, L;
	double fo;
	CSistema * sistema_a;

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
	CMatrizEsparsa MatrizRestCP();
	CMatrizEsparsa MatrizVmeta();
	CMatrizEsparsa MatrizBalPotencia();
	CMatrizEsparsa MatrizBalVazao();
	CMatrizEsparsa MatrizFuncProd();
	CMatrizEsparsa MatrizPhMax();
	CMatrizEsparsa MatrizReserva();
	CMatrizEsparsa MatrizCortesF();
	CMatrizEsparsa MatrizRestAux();
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
	vetorfloat2 LimRestAux();
	void MatrizRestricoesLineares(CMatrizEsparsa &MM);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);

public:
	Spcdec2ExtDE(CSistema * const);
	~Spcdec2ExtDE(void);
	double ResolverProblema(LMRow lambda, Results * results = NULL);
};
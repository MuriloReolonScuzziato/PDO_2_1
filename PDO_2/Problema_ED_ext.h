#pragma once

#include "ResultadosConj.h"

#include "MatrizEsparsa.h"

#include "gurobi_c++.h"
using namespace std;

#include "OPTtypes.h"
using namespace OPTtypes_di_unipi_it;

typedef vector<int> vetorint;

class CProblema_ED_ext
{
	GRBEnv   ambGRB;     // Ambiente Gurobi
	GRBModel modeloGRB;  // Modelo Gurobi
	GRBVar * vars;
	GRBConstr * constr;
	int n, no, na, flag1, flag2, flag3, flag4, flag7, T, rest_o, rest_tot;
	vetorfloat x, L;
	double fo;
	CSistema * sistema_a;

	void CriarVariaveis();
	void CriarFuncaoObjetivo();
	void CriarRestricoes();
	void CriarFuncaoObjetivo2();
	void FixarCondIniciais();

	CMatrizEsparsa CProblema_ED_ext::MatrizRestDemanda();
	CMatrizEsparsa CProblema_ED_ext::MatrizRestDemandaBarraUnica();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimFluxo();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimPhgMin();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimPhgMax();
	CMatrizEsparsa CProblema_ED_ext::MatrizBalHid();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimqMin();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimqMax();
	CMatrizEsparsa CProblema_ED_ext::MatrizTup();
	CMatrizEsparsa CProblema_ED_ext::MatrizTdown();
	CMatrizEsparsa CProblema_ED_ext::MatrizTupDown();
	CMatrizEsparsa CProblema_ED_ext::MatrizRampaUp();
	CMatrizEsparsa CProblema_ED_ext::MatrizRampaDown();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimPtMin();
	CMatrizEsparsa CProblema_ED_ext::MatrizLimPtMax();
	CMatrizEsparsa CProblema_ED_ext::MatrizRestCP();
	CMatrizEsparsa CProblema_ED_ext::MatrizVmeta();
	CMatrizEsparsa CProblema_ED_ext::MatrizBalPotencia();
	CMatrizEsparsa CProblema_ED_ext::MatrizBalVazao();
	CMatrizEsparsa CProblema_ED_ext::MatrizFuncProd();
	CMatrizEsparsa CProblema_ED_ext::MatrizPhMax();
	CMatrizEsparsa CProblema_ED_ext::MatrizReserva();
	CMatrizEsparsa CProblema_ED_ext::MatrizCortesF();
	CMatrizEsparsa CProblema_ED_ext::MatrizRestAux();
	vetorfloat2 CProblema_ED_ext::LimRestDemanda();
	vetorfloat2 CProblema_ED_ext::LimRestDemandaBarraUnica();
	vetorfloat2 CProblema_ED_ext::LimFluxo0();
	vetorfloat2 CProblema_ED_ext::LimFluxo2();
	vetorfloat2 CProblema_ED_ext::LimPhgMin();
	vetorfloat2 CProblema_ED_ext::LimPhgMax();
	vetorfloat2 CProblema_ED_ext::LimBalHid();
	vetorfloat2 CProblema_ED_ext::LimQMin();
	vetorfloat2 CProblema_ED_ext::LimQMax();
	vetorfloat2 CProblema_ED_ext::LimTup();
	vetorfloat2 CProblema_ED_ext::LimTdown();
	vetorfloat2 CProblema_ED_ext::LimTupDown();
	vetorfloat2 CProblema_ED_ext::LimRampaUp();
	vetorfloat2 CProblema_ED_ext::LimRampaDown();
	vetorfloat2 CProblema_ED_ext::LimPtMin();
	vetorfloat2 CProblema_ED_ext::LimPtMax();
	vetorfloat2 CProblema_ED_ext::LimRestCP();
	vetorfloat2 CProblema_ED_ext::LimVmeta();
	vetorfloat2 CProblema_ED_ext::LimBalPotenciaL();
	vetorfloat2 CProblema_ED_ext::LimBalVazaoL();
	vetorfloat2 CProblema_ED_ext::LimFuncProdL();
	vetorfloat2 CProblema_ED_ext::LimPhMax();
	vetorfloat2 CProblema_ED_ext::LimReserva();
	vetorfloat2 CProblema_ED_ext::LimCortesF();
	vetorfloat2 CProblema_ED_ext::LimRestAux();
	void CProblema_ED_ext::MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes);
	void CProblema_ED_ext::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);

public:
	CProblema_ED_ext(CSistema * const);
	~CProblema_ED_ext(void);
	double ResolverProblema(LMRow lambda);
};
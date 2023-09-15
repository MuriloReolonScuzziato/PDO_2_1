#pragma once

#include "Scndec2Results.h"

#include "MatrizEsparsa.h"

#include "gurobi_c++.h"
using namespace std;

#include "OPTtypes.h"
using namespace OPTtypes_di_unipi_it;

typedef vector<int> vetorint;

class Scndec2ExtDE
{
	GRBEnv   ambGRB;     // Ambiente Gurobi
	GRBModel modeloGRB;  // Modelo Gurobi
	GRBVar * vars;
	GRBConstr * constr;
	int n, no, flag1, flag2, flag3, flag4, flag7, T, rest_o, rest_tot;
	// n é o numero de variáveis de cada cenário
	// no é o número total de variáveis = n*n_cenarios
	vetorfloat L; //, x;
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
	CMatrizEsparsa MatrizLimPtMinMax();
	CMatrizEsparsa MatrizRestCP();
	CMatrizEsparsa MatrizVmeta();
	CMatrizEsparsa MatrizBalPotencia();
	CMatrizEsparsa MatrizBalVazao();
	CMatrizEsparsa MatrizFuncProd();
	CMatrizEsparsa MatrizPhMax();
	CMatrizEsparsa MatrizReserva();
	CMatrizEsparsa MatrizCortesF();
	CMatrizEsparsa MatrizRestNAnt();
	vetorfloat2 LimRestDemanda(int &cenario);
	vetorfloat2 LimRestDemandaBarraUnica(int &cenario);
	vetorfloat2 LimFluxo0(int &cenario);
	vetorfloat2 LimFluxo2(int &cenario);
	vetorfloat2 LimPhgMin();
	vetorfloat2 LimPhgMax();
	vetorfloat2 LimBalHid(int &cenario);
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
	vetorfloat2 LimReserva(int &cenario);
	vetorfloat2 LimCortesF();
	vetorfloat2 LimRestNAnt();
	void MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes, int &cenario);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor, int &cenario);

public:
	Scndec2ExtDE(CSistema * const);
	~Scndec2ExtDE(void);
	double ResolverProblema(LMRow lambda);
};
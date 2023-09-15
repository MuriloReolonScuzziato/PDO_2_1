#pragma once

#include "Scndec3Results.h"

class Scndec3SPC
{
	GRBModel * modeloGRB;  // Modelos do Gurobi (um modelo para cada cenário)
	GRBVar * vars;
	
	CSistema * sistema_a;
	int n, T, flag1, flag2, flag3, flag4, flag7, cenario, nStatus;
	double fo, sol_mipgap, sol_mipgap_curr, lb;
	vetorfloat x;
	vetorint i_rest_cen;
	int reset_stt;

	void CriarVariaveis();
	void CriarRestricoes();
	void AtualizarRestricoes();
	void FixarCondIniciais();
	void CriarFuncaoObjetivoRL(const vetorfloat &coefRL);
	void AtualizarLUB();

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

	void MatrizRestricoesLineares(CMatrizEsparsa &MM);
	void MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor);

public:
	// 2 construtores, um para criar modelo e outro para copiar o modelo (cenário 0 é criado e demais são somente copiados, com as devidas restrições atualizadas)
	Scndec3SPC(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi);
	Scndec3SPC(CSistema * const sistema_end, int cenario_a, Scndec3SPC &cenario0);
	~Scndec3SPC(void);

	int ResolverProblemaRL(Scndec3Results * resultadoGurobi, const vetorfloat2 &lambda, const int iter);
	double GetFO() { return fo;}
	int GetNVarPorComp() {return n;}		// comp de cada tipo de subproblema
	double GetMipgap(){return sol_mipgap_curr;}
	int GetStatus() {return nStatus;}
	//void SetMaxMipgap(double mipgap_a) {modeloGRB->getEnv().set(GRB_DoubleParam_MIPGap, mipgap_a); modeloGRB->update(); if (cenario == 0) cout << "Precisao setada para: " << mipgap_a << endl;}
	void SetMaxMipgap(double mipgap_a) {modeloGRB->getEnv().set(GRB_DoubleParam_MIPGap, mipgap_a); modeloGRB->update(); cout << "Precisao subp" << cenario << " setada para: " << mipgap_a << endl;}
	int GetResetStt() {return reset_stt;}
	void SetResetStt(int status) {reset_stt = status;}
};
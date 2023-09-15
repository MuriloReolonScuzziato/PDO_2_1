#pragma once

#include "Sistema.h"

#include "MatrizEsparsa.h"

typedef vector<string> vetorstring;

class Resultados_ED
{
	CSistema * sistema_a;
	int flag1, flag2, flag3, flag4, flag7;
	vetorfloat x, lambda;
	int n, CH, R, N, I, status;
	vetorint2 cascata;			// vetor com os indices das hidros de cada cascata
	double obj, CFO, CIO, DEF, CTL, CTQ, VFOL, mipgap;
	vetorstring stringVariaveis;
	ofstream * inFile;

	CMatrizEsparsa MatrizCalcFluxo(int n_a);
	CMatrizEsparsa MatrizCalcFluxo_cen(int n_a);
	void CriarStringVars();

public:
	Resultados_ED(CSistema * const sistema_end);
	~Resultados_ED(void);
	
	void CarregarResultadosED(double, vetorfloat, vetorfloat, int, double);
	void EscreverArquivo(string, double, double, string saida);

	vetorfloat GetX() {return x;}
	int Getn () {return n;}
	int GetCH () {return CH;}
	int GetR () {return R;}
	int GetI () {return I;}
	int GetN () {return N;}
};
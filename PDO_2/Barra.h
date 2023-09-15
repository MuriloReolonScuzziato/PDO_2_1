#pragma once

#include "Hidreletrica.h"
#include "Termeletrica.h"
#include "Demanda.h"

typedef vector<CHidreletrica* > vetorhidreletricaPtr;
typedef vector<CTermeletrica* > vetortermeletricaPtr;
typedef vector<CDemanda* > vetordemandaPtr;
//typedef vector<CHidreletrica> vetorhidreletrica;
//typedef vector<CTermeletrica> vetortermeletrica;
//typedef vector<CDemanda> vetordemanda;

class CBarra
{
	string nome_barra;
	int ident_barra;

public:
	vetordemandaPtr demandasPtr;
	vetorhidreletricaPtr hidrosPtr;
	vetortermeletricaPtr termosPtr;
	
	CBarra(string nome_barra_a, int ident);
	~CBarra(void);

	string GetNome() {return nome_barra;}
	int GetIdentBarra() {return ident_barra;}
	void AddDemanda(CDemanda * );
	void AddHidro(CHidreletrica * );
	void AddTermo(CTermeletrica * );
	int NDemandas() {return int (demandasPtr.size());}
	int NHidros() {return int (hidrosPtr.size());}
	int NTermos() {return int (termosPtr.size());}
};


#include "ConjSPD.h"

// Criar restrições
// ------------------------------------------------

// Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa ConjSPD::MatrizRestDemanda(int n_a, int no)	// Monta matriz da restrição de atendimento a demanda
{
	CMatrizEsparsa Agh(sistema_a->barrasVtr.size(), sistema_a->hidreletricasVtr.size());
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
				if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
					Agh.InserirElemento(b, r, 1);
	CMatrizEsparsa Agt(sistema_a->barrasVtr.size(), sistema_a->termeletricasVtr.size());
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
				if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
					Agt.InserirElemento(b, i, 1);
	// algoritmo para formar a Ybarra-> matriz B
	int ifr, ito;
	CMatrizEsparsa B(int(sistema_a->barrasVtr.size()), int(sistema_a->barrasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
	{
		ifr = sistema_a->linhasVtr[l].GetIndDeBarra() - 1;
		ito = sistema_a->linhasVtr[l].GetIndParaBarra() - 1;
		B.SubstituirElemento(ifr, ifr, B.GetElemento(ifr, ifr) + 100/(sistema_a->linhasVtr[l].GetReatancia()));
		B.SubstituirElemento(ifr, ito, B.GetElemento(ifr, ito) - 100/(sistema_a->linhasVtr[l].GetReatancia()));
		B.SubstituirElemento(ito, ifr, B.GetElemento(ito, ifr) - 100/(sistema_a->linhasVtr[l].GetReatancia()));
		B.SubstituirElemento(ito, ito, B.GetElemento(ito, ito) + 100/(sistema_a->linhasVtr[l].GetReatancia()));
	}
	B.RemoverColuna(sistema_a->GetBarraRef() - 1);
	B.MultiplicarPorEscalar( -1);

	CMatrizEsparsa Alin(1 * sistema_a->barrasVtr.size(), n_a);
	Alin.InserirMatriz(0, 0, 0 + sistema_a->barrasVtr.size() - 1, 0 + sistema_a->termeletricasVtr.size() - 1, &Agt, 0, 0);
	Alin.InserirMatriz(0, 0 + sistema_a->termeletricasVtr.size(),0 + sistema_a->barrasVtr.size() - 1, 0 + (sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &B, 0, 0);
	Alin.InserirMatriz(0, 0 + (sistema_a->barrasVtr.size() - 1 + sistema_a->termeletricasVtr.size()), 0 + sistema_a->barrasVtr.size() - 1, 0 + (sistema_a->barrasVtr.size() - 1 + sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size()) - 1, &Agh, 0, 0);
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	int c = sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 2*(sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa Adef(sistema_a->barrasVtr.size());
	Alin.InserirMatriz(0, c, 0 + sistema_a->barrasVtr.size() - 1,c + sistema_a->barrasVtr.size() - 1, &Adef, 0, 0);

	return Alin;
}
CMatrizEsparsa ConjSPD::MatrizRestDemandaBarraUnica(int n_a, int no)	// Monta matriz da restrição de atendimento a demanda
{
	size_t I = sistema_a->termeletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	a.RemoverTodosElementos();
	for (size_t i = 0; i < I; i++)		// pt
		a.InserirElemento(0, i, 1);
	for (size_t r = 0; r < R; r++)		// ph
		a.InserirElemento(0, r + sistema_a->termeletricasVtr.size(), 1);
	a.InserirElemento(0, sistema_a->termeletricasVtr.size() + 2*sistema_a->hidreletricasVtr.size(), 1);		// def
	Alin.JuntarColuna(&a);
	
	return Alin;
}
CMatrizEsparsa ConjSPD::MatrizLimFluxo(int n_a, int no)
{
	CMatrizEsparsa Alb(sistema_a->linhasVtr.size(),sistema_a->barrasVtr.size());
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			if (sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[b])
				Alb.InserirElemento(l, b, 1);
			else if (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[b])
				Alb.InserirElemento(l, b, -1);			
	Alb.RemoverColuna(sistema_a->GetBarraRef() - 1);
	CMatrizEsparsa TT(sistema_a->linhasVtr.size(),sistema_a->linhasVtr.size());
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
	{
		for ( size_t ll = 0; ll < sistema_a->linhasVtr.size(); ll++)
		{
			if (l == ll)
			{
				TT.InserirElemento(l, ll, 100 / (sistema_a->linhasVtr[l].GetReatancia()));
			}
		}
	}
	CMatrizEsparsa Alin(sistema_a->linhasVtr.size(), n_a);
	TT.MultiplicarPorMatriz(&Alb);
	Alin.InserirMatriz(0, 0 + (sistema_a->termeletricasVtr.size()), 0 + sistema_a->linhasVtr.size() - 1, 0 + (sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
	
	return Alin;
}
CMatrizEsparsa ConjSPD::MatrizReserva(int n_a, int no)
{
	int flag1 = int (1 - sistema_a->GetFlagBarraUnica());
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t B = flag1*(sistema_a->barrasVtr.size() - 1);
	int nt;
	CMatrizEsparsa Alin(0,n_a);
	CMatrizEsparsa a(1,n_a);
	int dd = I + B + R;
	a.RemoverTodosElementos();
	for (size_t r = 0; r < R; r++)
	{
		nt = (n_a / T);
		a.InserirElemento(0, dd + r, 1);
		a.InserirElemento(0, dd - R + r, - 1);
	}
	Alin.JuntarColuna(&a);
	
	return Alin;
}
// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 ConjSPD::LimRestDemanda(int n_a, int no)
{
	int cen = 0;
	int t = 0;
	if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
		t = no;
	else
	{
		cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
	}

	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, B, 2);
	IniciaMatriz(&Lim, 0);

	for (size_t b = 0; b < B; b++)
	{
		double D = 0;
		for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
			D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t);		// cen = 0, ou qq outro cenários, pois os T1 primeiros periodos são iguais (deterministicos)
		Lim[b][0] = 1;
		Lim[b][1] = D;
	}
	return Lim;
}
vetorfloat2 ConjSPD::LimRestDemandaBarraUnica(int n_a, int no)
{
	int cen = 0;
	int t = 0;
	if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
		t = no;
	else
	{
		cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
	}

	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1, 2);
	IniciaMatriz(&Lim, 0);
	double D = 0;
	for (size_t b = 0; b < B; b++)
	{
		for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
			D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t);
	}
	Lim[0][0] = 1;
	Lim[0][1] = D;

	return Lim;
}
vetorfloat2 ConjSPD::LimFluxo0(int n_a, int no)
{
	int L = sistema_a->linhasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, L, 2);
	IniciaMatriz(&Lim, 0);
	for (int l = 0; l < L; l++)
	{
		Lim[l][0] = 0;
		Lim[l][1] = sistema_a->linhasVtr[l].GetCapacidade();
	}
	return Lim;
}
vetorfloat2 ConjSPD::LimFluxo2(int n_a, int no)
{
	int L = sistema_a->linhasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, L, 2);
	IniciaMatriz(&Lim, 0);
	for (int l = 0; l < L; l++)
	{
		Lim[l][0] = 2;
		Lim[l][1] = - sistema_a->linhasVtr[l].GetCapacidade();
	}
	return Lim;
}
vetorfloat2 ConjSPD::LimReserva(int n_a, int no)
{
	int cen = 0;
	int t = 0;
	if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
		t = no;
	else
	{
		cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
	}

	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1, 2);
	IniciaMatriz(&Lim, 0);
	double D = 0;
	for (size_t b = 0; b < B; b++)
	{
		for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
			D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t);
	}
	Lim[0][0] = 2;
	Lim[0][1] = D*sistema_a->GetFatorReserva()/100;
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void ConjSPD::MatrizRestricoesLineares(int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int no, int &n_restricoes)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	int count = 0;
	if ( sistema_a->GetFlagBarraUnica() == false )
	{
		M = MatrizRestDemanda(n_a, no);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
		n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo(n_a, no);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
		n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo(n_a, no);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
		n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	else
	{
		M = MatrizRestDemandaBarraUnica(n_a, no);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
		n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	M = MatrizReserva(n_a, no);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
}
void ConjSPD::MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int no)
{
	vetorfloat2 L;
	if ( sistema_a->GetFlagBarraUnica() == false )
	{
		L = LimRestDemanda(n_a, no);
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo0(n_a, no);
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo2(n_a, no);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	else
	{
		L = LimRestDemandaBarraUnica(n_a, no);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimReserva(n_a, no);
	AlocarLimites(&L, LimTipo, LimValor);
}
// ------------------------------------------------

ConjSPD::ConjSPD(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;

	n_nos = sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios();
	modelosGRB.resize(n_nos);
	vars.resize(n_nos);
	for (int nno = 0; nno < n_nos; nno++)
		modelosGRB[nno] = new GRBModel(ambiente_gurobi);		// um modelo de otimização para cada no
	n.resize(n_nos);
	for (int nno = 0; nno < n_nos; nno++)
	{
		n[nno] = sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 2*sistema_a->hidreletricasVtr.size() + sistema_a->barrasVtr.size();
		if ( sistema_a->GetFlagBarraUnica() == true )		// Remover teta, def e deixar uma variável de deficit por período
			n[nno] -= (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	}

	// Criar variáveis
	for (int nno = 0; nno < n_nos; nno++)
	{
		CriarVariaveis(nno);
		modelosGRB[nno]->update();	// Atualiza o modelo Gurobi.
		vars[nno] = modelosGRB[nno]->getVars();
	}

	// Adicionar restrições
	for (int nno = 0; nno < n_nos; nno++)
	{
		CriarRestricoes(nno);
		modelosGRB[nno]->update();	// Atualiza o modelo Gurobi.
	}

	// Definir ajustes do solver
	for (int nno = 0; nno < n_nos; nno++)
	{
		//modelosGRB[nno]->write("probD.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modelosGRB[nno]->getEnv().set(GRB_IntParam_OutputFlag, 0);
		
		//modelosGRB[nno]->getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[nno]->getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		//modelosGRB[nno]->getEnv().set(GRB_DoubleParam_MIPGap, 0.0001);	// Define o gap de tolerancia
		//modelosGRB[nno]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[nno]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modelosGRB[nno]->getEnv().set(GRB_IntParam_Threads, 4);
		//modelosGRB[nno].getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[nno].getEnv().set(GRB_IntParam_PreSparsify, 1);
		//modelosGRB[nno]->getEnv().set(GRB_IntParam_Method, 0);
		//modelosGRB[nno]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);

		//modelosGRB[nno].getEnv().set(GRB_IntParam_MIPFocus, 2);
		//modelosGRB[nno].getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

		modelosGRB[nno]->getEnv().set(GRB_DoubleParam_TimeLimit, 60);		// Limita o tempo de resolução do problema
		//modelosGRB[nno]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-2);	// Define tolerancia de otimalidade
	}
}
ConjSPD::~ConjSPD(void)
{
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
}

void ConjSPD::CriarVariaveis(int no)
{
	try 
	{
		// variáveis em cada problema somente para um nó
		int cen = 0;
		int t = 0;
		if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
			t = no;
		else
		{
			cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
		}
		size_t I = sistema_a->termeletricasVtr.size();
		size_t R = sistema_a->hidreletricasVtr.size();
		size_t B = sistema_a->barrasVtr.size() - 1;
		//x = [pt teta ph phmax def]
		for (size_t i = 0; i < I; i++)	//pt
			modelosGRB[no]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
		if ( sistema_a->GetFlagBarraUnica() == false)
			for (size_t b = 0; b < B; b++)	//teta
				modelosGRB[no]->addVar(-6.2832, 6.2832, 0.0, GRB_CONTINUOUS, "");
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
		{
			double phmax = 0;
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
			modelosGRB[no]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
		}
		// Esse tipo de decomposição sempre terá a variável phmax
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
		{
			double phmax = 0;
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
			modelosGRB[no]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
		}
		if ( sistema_a->GetFlagBarraUnica() == false)
		{
			for (size_t b = 0; b < B + 1; b++)	//def
			{
				double cap_d = 0;
				for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
					cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
				modelosGRB[no]->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
			}
		}
		else
		{
			double cap_d = 0;			//def barra unica
			for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
					cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
			modelosGRB[no]->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
void ConjSPD::CriarRestricoes(int no)
{
	try
	{
		vetorint indexL, indexC, nnZ, LimTipo;
		vetorfloat indexV, LimValor;
		int n_restricoes = 0;
		indexL.resize(0);indexC.resize(0);indexV.resize(0);nnZ.resize(0);LimTipo.resize(0);LimValor.resize(0);
		nnZ.push_back(0);
		MatrizRestricoesLineares(n[no], &indexL, &indexC, &indexV, &nnZ , no, n_restricoes);
		MatrizLimitesLineares(n[no], &LimTipo, &LimValor, no);

		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		for (int i = 0; i < n_restricoes; i++)
		{
			for (int j = nnZ[i]; j < nnZ[i + 1]; j++)
			{
				coeficiente = indexV[j];
				variavel = vars[no][indexC[j]];
				restricao.addTerms( &coeficiente, &variavel, 1);
			}
			switch (LimTipo[i]) {
			case 0:
				modelosGRB[no]->addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modelosGRB[no]->addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modelosGRB[no]->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
				break;
			default:
				cout << "Tipo inválido de restrição adicionada" << endl;
			}
			restricao.clear();
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
void ConjSPD::CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int no)
{
	try 
	{
		int cen = 0;
		int t = 0;
		if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
			t = no;
		else
		{
			cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
		}

		//x = [pt teta ph phmax def]
		int flag1 = int (1 - sistema_a->GetFlagBarraUnica());

		// Termos lineares
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int nt_dual = lambda->size() / T;
		int delta_dual = nt_dual * no;
		double deltaT;
		int delta = 0;
		
		// Deficit
		delta = sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size();
		if (t < sistema_a->GetTt1())
		{
			deltaT = sistema_a->GetDeltaT1();
			if (sistema_a->GetFlagBarraUnica() == false)
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					vars[no][delta + b].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT);		//def
			else
				vars[no][delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT);		//def
		}
		else
		{
			deltaT = sistema_a->GetDeltaT2();
			if (sistema_a->GetFlagBarraUnica() == false)
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					vars[no][delta + b].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//def
			else
				vars[no][delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//def
		}

		// Termos da RL
		// lambda [pt ph v d phmax]
		//x = [pt teta ph phmax def]
		// Constantes

		// Lineares
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
		{
			vars[no][i].set(GRB_DoubleAttr_Obj, lambda->at(i + delta_dual));		//pt
		}
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			vars[no][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1)].set(GRB_DoubleAttr_Obj, lambda->at(r + sistema_a->termeletricasVtr.size() + delta_dual));		//ph
			vars[no][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()].set(GRB_DoubleAttr_Obj, lambda->at(r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + delta_dual));		//phmax
		}

	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
int ConjSPD::ResolverProblemaRL(ResultadosConj * resultadoGurobi, const vetorfloat * const lambda, const int iter, int no)
{
	// Resolve o subproblema de um nó especifico
	try
	{
		// Criar função objetivo, a seleçao dos lambdas de interesse é feita dentro da funçao abaixo
		CriarFuncaoObjetivoRL(lambda, no);
		modelosGRB[no]->update();	// Atualiza o modelo Gurobi.
				
		modelosGRB[no]->reset();

		// Otimizar
		modelosGRB[no]->optimize();

		// Salvar resultados
		x_n.clear();
		x_n.resize(n[no]);
		int nStatus = modelosGRB[no]->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The D subproblem was not solved until optimality! " << nStatus << endl;
		if ((nStatus == 2) || (nStatus == 9) || (nStatus == 13))
		{
			fo_n = double(modelosGRB[no]->get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n[no]; i++)
			{	
				x_n[i] = double(vars[no][i].get(GRB_DoubleAttr_X));
			}
		}
		else
		{
			for (int i = 0; i < n[no]; i++)
			{	
				x_n[i] = 0;
			}
			fo_n = 0;
		}
		resultadoGurobi->GravarSolucao(fo_n, x_n, nStatus, resultadoGurobi->GetCH() + resultadoGurobi->GetN() * resultadoGurobi->GetR() + resultadoGurobi->GetI() + no);
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->GravarSolucao(fo_n, x_n, e.getErrorCode(), resultadoGurobi->GetCH() + resultadoGurobi->GetN() * resultadoGurobi->GetR() + resultadoGurobi->GetI() + no);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
	return 0;
}

// Caso os subproblemas sejam Easy Components
ConjSPD::ConjSPD(CSistema * const sistema_end)
{
	sistema_a = sistema_end;
	n_nos = sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios();
	n.resize(n_nos);
	for (int nno = 0; nno < n_nos; nno++)
	{
		n[nno] = sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 2*sistema_a->hidreletricasVtr.size() + sistema_a->barrasVtr.size();
		if ( sistema_a->GetFlagBarraUnica() == true )		// Remover teta, def e deixar uma variável de deficit por período
			n[nno] -= (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	}

	Blinha.resize(n_nos);
	Bcoluna.resize(n_nos);
	Bnnz.resize(n_nos);
	dtipo.resize(n_nos);
	Bvalor.resize(n_nos);
	dvalor.resize(n_nos);
	Nrestricoes.resize(n_nos);

	for (int nno = 0; nno < n_nos; nno++)
		GerarRestricoes(nno);
}
void ConjSPD::GerarRestricoes(int no)
{
	int n_restricoes = 0;
	Blinha[no].resize(0);Bcoluna[no].resize(0);Bvalor[no].resize(0);Bnnz[no].resize(0);dtipo[no].resize(0);dvalor[no].resize(0);
	Bnnz[no].push_back(0);
	MatrizRestricoesLineares(n[no], &Blinha[no], &Bcoluna[no], &Bvalor[no], &Bnnz[no] , no, n_restricoes);
	MatrizLimitesLineares(n[no], &dtipo[no], &dvalor[no], no);
	Nrestricoes[no] = n_restricoes;

	//// Ordenar elementos da matriz de restriçao (para serem lidas conforme GetBDesc())
	//vector<pair<size_t, myiter> > order(Bcoluna[no].size());
	//size_t nn = 0;
	//for (myiter it = Bcoluna[no].begin(); it != Bcoluna[no].end(); ++it, ++nn)
	//	order[nn] = make_pair(nn, it);
	//sort(order.begin(), order.end(), ordering());
	//sort_from_ref(Bcoluna[no], order);
	//sort_from_ref(Blinha[no], order);
	//sort_from_ref(Bvalor[no], order);

	//// considerando a ordem do vetor de linhas tb
	vector<pair<size_t, pair<myiter, myiter>> > order(Bcoluna[no].size());
	size_t nn = 0;
	myiter it2 = Blinha[no].begin();
	for (myiter it = Bcoluna[no].begin(); it != Bcoluna[no].end(); ++it, ++nn, ++it2)
		order[nn] = make_pair(nn, make_pair(it, it2));
	sort(order.begin(), order.end(), ordering2());
	sort_from_ref2(Bcoluna[no], order);
	sort_from_ref2(Blinha[no], order);
	sort_from_ref2(Bvalor[no], order);
}
void ConjSPD::GerarVarBounds(int no, double *lbd , double *ubd)
{
	// variáveis em cada problema somente para um nó
	int cen = 0;
	int t = 0;
	if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
		t = no;
	else
	{
		cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
	}
	size_t I = sistema_a->termeletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t B = sistema_a->barrasVtr.size() - 1;
	//x = [pt teta ph phmax def]
	int iub = 0;
	int ilb = 0;
	for (size_t i = 0; i < I; i++)	//pt
	{
		ubd[iub++] = double (sistema_a->termeletricasVtr[i].GetPmax());
		lbd[ilb++] = 0;
	}
	if ( sistema_a->GetFlagBarraUnica() == false)
		for (size_t b = 0; b < B; b++)	//teta
		{
			ubd[iub++] = 6.2832;
			lbd[ilb++] = -6.2832;
		}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
	{
		double phmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();

		ubd[iub++] = double (phmax);
		lbd[ilb++] = 0;
	}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
	{
		double phmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
		ubd[iub++] = double (phmax);
		lbd[ilb++] = 0;
	}
	if ( sistema_a->GetFlagBarraUnica() == false)
	{
		for (size_t b = 0; b < B + 1; b++)	//def
		{
			double cap_d = 0;
			for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
				cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
			ubd[iub++] = cap_d;
			lbd[ilb++] = 0;
		}
	}
	else
	{
		double cap_d = 0;			//def barra unica
		for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
				cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
		ubd[iub++] = cap_d;
		lbd[ilb++] = 0;
	}
}
void ConjSPD::GerarCoefFuncaoObjetivo(int no, double *cst)
{
	int cen = 0;
	int t = 0;
	if (no < sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
		t = no;
	else
	{
		cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		t = (no - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1();
	}

	//x = [pt teta ph phmax def]
	int flag1 = int (1 - sistema_a->GetFlagBarraUnica());
	// Termos lineares
	double deltaT;
	int delta = 0;
	for (int i = 0; i < n[no]; i++)		// zerar todos elementos de cst
		cst[i] = 0;
		
	// Deficit
	delta = sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size();
	if (t < sistema_a->GetTt1())
	{
		deltaT = sistema_a->GetDeltaT1();
		if (sistema_a->GetFlagBarraUnica() == false)
			for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				cst[delta + b] = - sistema_a->GetCustoDeficit()*deltaT;		//def
		else
			cst[delta] = - sistema_a->GetCustoDeficit()*deltaT;		//def
	}
	else
	{
		deltaT = sistema_a->GetDeltaT2();
		if (sistema_a->GetFlagBarraUnica() == false)
			for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				cst[delta + b] = - sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
		else
			cst[delta] = - sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
	}
}
void ConjSPD::GetBDesc( int comp , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd )
{
	// fazer uma funçao GetBDesc nessa classe, que passa os valores da matriz jah criada no constructor para os respectivos vetores!!
	int col_ant = 0;
	int beg = 0;
	int ind = 0;
	Bbeg[ind++] = beg;
	for (int i = 0; i < Bcoluna[comp].size(); i++)		// loop para todos elementos não zeros
	{
		Bind[i] = Blinha[comp][i];
		Bval[i] = Bvalor[comp][i];
		if (col_ant != Bcoluna[comp][i])
		{
			for (int j = col_ant; j < Bcoluna[comp][i]; j++)
				Bbeg[ind++] = beg;
		}
		beg++;
		col_ant = Bcoluna[comp][i];
	}
	for (int j = col_ant; j < n[comp]; j++)			// loop para completar o vetor com o numero de colunas(= num. variáveis)
			Bbeg[ind++] = beg;

	for (int i = 0; i < Nrestricoes[comp]; i++)			// loop no numero de restriçoes (=num. linhas da matriz)
	{
		switch (dtipo[comp][i]) {
		case 0:
			rhs[i] = dvalor[comp][i];
			lhs[i] = - GRB_INFINITY;
			break;
		case 1:
			rhs[i] = dvalor[comp][i];
			lhs[i] = dvalor[comp][i];
			break;
		case 2:
			rhs[i] = GRB_INFINITY;
			lhs[i] = dvalor[comp][i];
			break;
		default:
			cout << "Tipo inválido de restrição adicionada" << endl;
		}
	}

	GerarVarBounds(comp, lbd, ubd);				// definir os limites das variaveis
	
	GerarCoefFuncaoObjetivo(comp, cst);			// definir a funçao objetivo (inverter sinal pois o prob. no MP é de max.
}
int ConjSPD::GetANZ(int comp)
{
	// Retorna o numero de multiplicadores de lagrange considerado em cada subproblema
	return ( sistema_a->termeletricasVtr.size() + 2 * sistema_a->hidreletricasVtr.size() );
}
void ConjSPD::GetADesc( int comp , int *Abeg , int *Aind , double *Aval )
{
	// Termos da RL
	// lambda [pt ph v d phmax]
	// x = [pt teta ph phmax def]
	int nt_dual = sistema_a->termeletricasVtr.size() + 4 * sistema_a->hidreletricasVtr.size();
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	int ind = 0;
	int iAbeg = 0;
	Abeg[iAbeg++] = 0;
	int delta_dual = nt_dual * comp;
	// loop para todas as variáris do subproblema (BNC), deve ser na ordem das colunas, variáveis...
	for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)		// pt
	{
		Aind[ind] = i + delta_dual;		// numero da linha correspondente ao lambda
		//Aval[ind++] = 1;
		Aval[ind++] = sistema_a->GetPrecondidioner(i + delta_dual);
		Abeg[iAbeg++] = ind;
	}
	if ( sistema_a->GetFlagBarraUnica() == false)		// teta (n tem lambda na funçao objetivo)
		for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
			Abeg[iAbeg++] = ind;		// adiciona-se um elemento em Abeg porem com o valor dos anteriores, indicando que n teve-se adição de elementos n nulos
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)		// ph
	{
		Aind[ind] = r + sistema_a->termeletricasVtr.size() + delta_dual;
		//Aval[ind++] = 1;
		Aval[ind++] = sistema_a->GetPrecondidioner(r + sistema_a->termeletricasVtr.size() + delta_dual);
		Abeg[iAbeg++] = ind;
	}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)		// phmax
	{
		Aind[ind] = r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + delta_dual;
		//Aval[ind++] = 1;
		Aval[ind++] = sistema_a->GetPrecondidioner(r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + delta_dual);
		Abeg[iAbeg++] = ind;
	}
	if ( sistema_a->GetFlagBarraUnica() == false)		// def	(n tem lambda na funçao objetivo)
		for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			Abeg[iAbeg++] = ind;
	else
		Abeg[iAbeg++] = ind;
}
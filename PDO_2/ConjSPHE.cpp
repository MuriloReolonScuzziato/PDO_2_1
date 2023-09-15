#include "ConjSPHE.h"

// Criar restrições
// ------------------------------------------------
// Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa ConjSPHE::MatrizLimPhgMin(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ,JJ);
	int cont = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmin());
		cont++;
	}
	CMatrizEsparsa Alin(JJ, n_a);
	Alin.InserirMatriz(0,5, JJ - 1,5 + JJ - 1, &A1, 0, 0);
	Alin.InserirMatriz(0,5 + 2*JJ, JJ - 1,5 + 3*JJ - 1, &A2, 0, 0);
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizLimPhgMax(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades());
		cont++;
	}
	CMatrizEsparsa Alin(JJ,n_a);
	Alin.InserirMatriz(0, 5, JJ - 1, 5 + JJ - 1, &A1, 0, 0);
	Alin.InserirMatriz(0, 5 + 2*JJ, JJ - 1, 5 + 3*JJ - 1, &A2, 0, 0);
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizLimqMin(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetQmin());
		cont++;
	}
	CMatrizEsparsa Alin(JJ, n_a);
	Alin.InserirMatriz(0, 5 + JJ, JJ - 1, 5 + 2*JJ - 1, &A1, 0, 0);
	Alin.InserirMatriz(0, 5 + 2*JJ, JJ - 1, 5 + 3*JJ - 1, &A2, 0, 0);
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizLimqMax(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades());
		cont++;
	}
	CMatrizEsparsa Alin(JJ,n_a);
	Alin.InserirMatriz(0, 5 + JJ, JJ - 1, 5 + 2*JJ - 1, &A1, 0, 0);
	Alin.InserirMatriz(0, 5 + 2*JJ, JJ - 1, 5 + 3*JJ - 1, &A2, 0, 0);
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizBalPotencia(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa Alin(1, n_a);
	int jj = 5;
	Alin.RemoverTodosElementos();
	Alin.InserirElemento(0, 0, 1);
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		Alin.InserirElemento(0, jj + j, -1);
	}
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizBalVazao(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa Alin(1, n_a);
	int dd = 2;
	int jj = dd + JJ + 3;
	Alin.RemoverTodosElementos();
	Alin.InserirElemento(0, dd, 1);
	Alin.InserirElemento(0, dd + 1, - 1);
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		Alin.InserirElemento(0, jj + j, -1);
	}
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizFuncProd(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int vv = 1;
	int jj = vv + 4 + JJ;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		a.RemoverTodosElementos();
		a.InserirElemento(0, jj - JJ + j, 1);															// phg
		a.InserirElemento(0, vv, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].CoefPhgV());				// v
		a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].CoefPhgQ());			// q
		a.InserirElemento(0, vv + 1, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].CoefPhgD());			// d
		a.InserirElemento(0, jj + j + JJ, sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmin());		// z -> adiciona um termo (1 - z)*Pmin na função, para phg <= pmin qdo z = 0; o rhs estava negativo, devido a aproximação linear, fazendo com q a hidro ficasse obrigatoriamente ligada
		if (sistema_a->hidreletricasVtr[nu].GetInflueVert() == 1 )		// Se o s influencia já tenho contabilizado isso no d, caso contrário só subtrair o valor de s de d
		{
		}
		else		// vertimento n influencia no canal de fuga, entao subtraio o valor de s da defluencia
		{
			a.InserirElemento(0, vv + 2, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].CoefPhgS());		// s
		}
		Alin.JuntarColuna(&a);
	}
	return Alin;
}
CMatrizEsparsa ConjSPHE::MatrizPhMax(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	CMatrizEsparsa Alin(1, n_a);
	int dd = 4;
	int jj = dd + 1 + 2*JJ;
	Alin.RemoverTodosElementos();
	Alin.InserirElemento(0, dd, 1);
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		Alin.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades());
	}
	return Alin;
}

// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 ConjSPHE::LimPhgMin(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ, 2);
	IniciaMatriz(&Lim, 0);
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		Lim[j][0] = 2;
	return Lim;
}
vetorfloat2 ConjSPHE::LimPhgMax(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 ConjSPHE::LimQMin(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ, 2);
	IniciaMatriz(&Lim, 0);
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		Lim[j][0] = 2;
	return Lim;
}
vetorfloat2 ConjSPHE::LimQMax(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 ConjSPHE::LimBalPotenciaL(int n_a, int nu)
{
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 ConjSPHE::LimBalVazaoL(int n_a, int nu)
{
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 ConjSPHE::LimFuncProdL(int n_a, int nu)
{
	int JJ = 0;
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
		JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ, 2);
	IniciaMatriz(&Lim, 0);
	for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
	{
		Lim[j][0] = 0;
		//Lim[j][1] = sistema_a->hidreletricasVtr[nu].grupoVtr[j].CoefPhg();
		Lim[j][1] = sistema_a->hidreletricasVtr[nu].grupoVtr[j].CoefPhg() + sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmin();
	}
	return Lim;
}
vetorfloat2 ConjSPHE::LimPhMax(int n_a, int nu)
{
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void ConjSPHE::MatrizRestricoesLineares(int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int nu, int &n_restricoes)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	int count = 0;
	M = MatrizLimPhgMin(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPhgMax(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimqMin(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimqMax(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalPotencia(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalVazao(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizFuncProd(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizPhMax(n_a, nu);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
}
void ConjSPHE::MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int nu)
{
	vetorfloat2 L;
	L = LimPhgMin(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhgMax(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMin(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMax(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalPotenciaL(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalVazaoL(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimFuncProdL(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhMax(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
}
// ------------------------------------------------

ConjSPHE::ConjSPHE(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;

	n_usinas = sistema_a->hidreletricasVtr.size();
	modelosGRB.resize(n_usinas);
	vars.resize(n_usinas);
	for (int nu = 0; nu < n_usinas; nu++)
		modelosGRB[nu] = new GRBModel(ambiente_gurobi);		// um modelo de otimização para cada usina
	n.resize(n_usinas);
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		n[r] = 5 + 3*sistema_a->hidreletricasVtr[r].GetNGrupos();

	// Criar variáveis
	for (int nu = 0; nu < n_usinas; nu++)
	{
		CriarVariaveis(nu);
		modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.
		vars[nu] = modelosGRB[nu]->getVars();
	}
	
	// Adicionar restrições
	for (int nu = 0; nu < n_usinas; nu++)
	{
		CriarRestricoes(nu);
		modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.
	}

	// Definir ajustes do solver
	for (int nu = 0; nu < n_usinas; nu++)
	{
		//modelosGRB[nu]->write("probHE.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modelosGRB[nu]->getEnv().set(GRB_IntParam_OutputFlag, 0);
		
		//modelosGRB[nu].getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[nu].getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		modelosGRB[nu]->getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);	// Define o gap de tolerancia
		//modelosGRB[nu]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[nu]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modelosGRB[nu]->getEnv().set(GRB_IntParam_Threads, 1);
		//modelosGRB[nu]->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[nu].getEnv().set(GRB_IntParam_PreSparsify, 1);
		//modelosGRB[nu].getEnv().set(GRB_IntParam_Method, 0);

		//modelosGRB[nu].getEnv().set(GRB_IntParam_MIPFocus, 2);
		//modelosGRB[nu].getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

		modelosGRB[nu]->getEnv().set(GRB_DoubleParam_TimeLimit, 60);		// Limita o tempo de resolução do problema
	}
}
ConjSPHE::~ConjSPHE(void)
{
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
}

void ConjSPHE::CriarVariaveis(int nu)
{
	try 
	{
		// variáveis em cada problema somente para um período
		size_t R = sistema_a->hidreletricasVtr.size();
		// Problema de 1 período para cada usina
		//x = [ph v d s phmax phg q z]
		double phmax;
		//ph
		phmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
			phmax = phmax + sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades();
		modelosGRB[nu]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
		//v
		modelosGRB[nu]->addVar(double(sistema_a->hidreletricasVtr[nu].GetVmin()), double(sistema_a->hidreletricasVtr[nu].GetVmax()), 0.0, GRB_CONTINUOUS, "");
		//d
		double qhmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
			qhmax = qhmax + double (sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades());
		modelosGRB[nu]->addVar(0, qhmax + sistema_a->hidreletricasVtr[nu].GetSmax(), 0.0, GRB_CONTINUOUS, "");
		//s
		modelosGRB[nu]->addVar(0, sistema_a->hidreletricasVtr[nu].GetSmax(), 0.0, GRB_CONTINUOUS, "");
		//phmax
		phmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)
			phmax = phmax + sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades();
		modelosGRB[nu]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
		//phg
		for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)	
			modelosGRB[nu]->addVar(0, double (sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
		//q
		for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)	
			modelosGRB[nu]->addVar(0, double (sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[nu].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
		//z
		if (sistema_a->GetFlagVarBin() == true)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)	
				modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");
		}
		else
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[nu].GetNGrupos(); j++)	
				modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
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
void ConjSPHE::CriarRestricoes(int nu)
{
	try
	{
		vetorint indexL, indexC, nnZ, LimTipo;
		vetorfloat indexV, LimValor;
		int n_restricoes = 0;
		indexL.resize(0);indexC.resize(0);indexV.resize(0);nnZ.resize(0);LimTipo.resize(0);LimValor.resize(0);
		nnZ.push_back(0);
		MatrizRestricoesLineares(n[nu], &indexL, &indexC, &indexV, &nnZ, nu, n_restricoes);
		MatrizLimitesLineares(n[nu], &LimTipo, &LimValor, nu);

		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		for (int i = 0; i < n_restricoes; i++)
		{
			for (int j = nnZ[i]; j < nnZ[i + 1]; j++)
			{
				coeficiente = indexV[j];
				variavel = vars[nu][indexC[j]];
				restricao.addTerms( &coeficiente, &variavel, 1);
			}
			switch (LimTipo[i]) {
			case 0:
				modelosGRB[nu]->addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modelosGRB[nu]->addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modelosGRB[nu]->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
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
void ConjSPHE::CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int nu)
{
	try 
	{
		//GRBLinExpr fo;
		////x = [ph v d s phmax phg q z]
		//double * coeffs;
		//coeffs = new double[n[nu]];

		//// Termos da RL
		//// lambda [ph v d phmax]
		//// Constantes
		//
		//// Lineares
		//for (int i = 0; i < n[nu]; i++)
		//	coeffs[i] = 0;
		//coeffs[0] = - lambda->at(0);		//ph
		//coeffs[1] = - lambda->at(1);		//v
		//coeffs[2] = - lambda->at(2);		//d
		//coeffs[4] = - lambda->at(3);		//phmax
		//fo.addTerms(coeffs, vars[nu], n[nu]);

		//delete coeffs;
		//modelosGRB[nu]->setObjective(fo, GRB_MINIMIZE);

		//modelosGRB[nu]->set(GRB_IntAttr_ModelSense, 1);			// minimize obj. f.

		vars[nu][0].set(GRB_DoubleAttr_Obj, - lambda->at(0));		//ph
		vars[nu][1].set(GRB_DoubleAttr_Obj, - lambda->at(1));		//v
		vars[nu][2].set(GRB_DoubleAttr_Obj, - lambda->at(2));		//d
		vars[nu][4].set(GRB_DoubleAttr_Obj, - lambda->at(3));		//phmax

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
int ConjSPHE::ResolverProblemaRL(ResultadosConj * resultadoGurobi, const vetorfloat * const lambda, const int iter, int usina, int no)
{
	// Resolve o subproblema de uma usina em um nó especifico
	try
	{
		vetorfloat lambda_rn;		// vetores para um período
		lambda_rn.resize(4);
		int nStatus;

		// Criar função objetivo ("pesca" lambdas de interesse)
		for (int c = 0; c < 4; c++)
			lambda_rn[c] = lambda->at(no*(sistema_a->termeletricasVtr.size() + 4*sistema_a->hidreletricasVtr.size()) + sistema_a->termeletricasVtr.size() + c * sistema_a->hidreletricasVtr.size() + usina);
		
		CriarFuncaoObjetivoRL(&lambda_rn, usina);
		modelosGRB[usina]->update();	// Atualiza o modelo Gurobi.

		// ajustes do solver

		modelosGRB[usina]->reset();

		// Otimizar
		modelosGRB[usina]->optimize();

		// Salvar resultados
		x_rn.clear();
		x_rn.resize(n[usina]);
		nStatus = modelosGRB[usina]->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The HE subproblem was not solved until optimality! " << nStatus << endl;
		if ((nStatus == 2) || (nStatus == 9))
		{
			fo_rn = double(modelosGRB[usina]->get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n[usina]; i++)
			{	
				x_rn[i] = double(vars[usina][i].get(GRB_DoubleAttr_X));
			}
		}
		else
		{
			for (int i = 0; i < n[usina]; i++)
			{	
				x_rn[i] = 0;
			}
			fo_rn = 0;
		}
		resultadoGurobi->GravarSolucao(fo_rn, x_rn, nStatus, resultadoGurobi->GetCH() + no*resultadoGurobi->GetR() + usina);

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->GravarSolucao(fo_rn, x_rn, e.getErrorCode(), resultadoGurobi->GetCH() + no*resultadoGurobi->GetR() + usina);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}
#include "Spcdec3SPHA.h"

// Criar restrições
// ------------------------------------------------
// Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa Spcdec3SPHA::MatrizBalHid(int n_a, int ch)
{
	int flag2 = int (sistema_a->GetFlagVfol());
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	double deltaT;
	int tempo_viagem;
	size_t R = cascata[ch].size();
	CMatrizEsparsa Addd(R * T, R);
	CMatrizEsparsa Add(R, R);
	std::vector<int>::iterator it;
	for (int t = 0; t < T; t++)
	{
		deltaT = sistema_a->GetDeltaT1();
		Add.RemoverTodosElementos();
		for (int r = 0; r < R; r++)
		{
			if (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetTempoViagem() < deltaT)
				tempo_viagem = 0;
			else
				tempo_viagem = int (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetTempoViagem() / deltaT);
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetUsinaJusante() != 0) && (tempo_viagem < sistema_a->GetTt2()))
			{
				it = find(cascata[ch].begin(), cascata[ch].end(), sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetUsinaJusante() - 1);
				Add.InserirElemento(int(it-cascata[ch].begin()), r, double (-0.0036 * deltaT));
			}
		}
		Addd.InserirMatriz(t * R,0,(t + 1) * R - 1,R - 1, &Add, 0, 0);
	}
	CMatrizEsparsa Addd2(R * T,R);
	for (int t = 0; t < T; t++)
	{
		deltaT = sistema_a->GetDeltaT2();
		Add.RemoverTodosElementos();
		for (int r = 0; r < R; r++)
		{
			if (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetTempoViagem() < deltaT)
				tempo_viagem = 0;
			else
				tempo_viagem = int (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetTempoViagem() / deltaT);
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetUsinaJusante() != 0) && (tempo_viagem < sistema_a->GetTt2()))
			{
				it = find(cascata[ch].begin(), cascata[ch].end(), sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetUsinaJusante() - 1);
				Add.InserirElemento(int(it-cascata[ch].begin()), r, double (-0.0036 * deltaT));
			}
		}
		Addd2.InserirMatriz(t * R,0,(t + 1) * R - 1,R - 1, &Add, 0, 0);
	}
	CMatrizEsparsa Ad(R * T);
	CMatrizEsparsa MSoma1(R, R);
	CMatrizEsparsa MSoma2(R, R);
	int cen;
	for (int t = 0; t < T; t++)
	{
		MSoma1.RemoverTodosElementos();
		MSoma2.RemoverTodosElementos();
		if (t < sistema_a->GetTt1())
		{
			deltaT = sistema_a->GetDeltaT1();
			for (int tt = 0; tt <= t; tt++)
			{
				MSoma1.InserirMatriz(0, 0, R - 1, R - 1, &Addd, (t - tt)*R, 0);
				MSoma2.InserirMatriz( 0, 0, R - 1, R - 1, &Ad, t*R, tt*R);
				MSoma2.MultiplicarPorEscalar(double (0.0036 * deltaT));
				MSoma1.SomarComMatriz(&MSoma2);
				Ad.InserirMatriz(t*R,tt*R,(t + 1)*R - 1,(tt + 1)*R - 1, &MSoma1, 0, 0);
			}
		}
		else
		{
			deltaT = sistema_a->GetDeltaT2();
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			if (t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) == sistema_a->GetTt1())
			{
				for (int tt = 0; tt <= t; tt++)
				{
					int tt_ = (t - tt) - cen*(sistema_a->GetTt2() - sistema_a->GetTt1());
					if (tt_ >= 0 )
					{
						MSoma1.InserirMatriz(0, 0, R - 1, R - 1, &Addd2, (tt_)*R, 0);
						MSoma2.InserirMatriz(0, 0, R - 1, R - 1, &Ad, t*R, tt*R);
						MSoma2.MultiplicarPorEscalar(double (0.0036 * deltaT));
						MSoma1.SomarComMatriz(&MSoma2);
						Ad.InserirMatriz(t*R,tt*R,(t + 1)*R - 1,(tt + 1)*R - 1, &MSoma1, 0, 0);
					}
					else
					{
						MSoma2.InserirMatriz(0, 0, R - 1, R - 1, &Ad, t*R, tt*R);
						MSoma2.MultiplicarPorEscalar(double (0.0036 * deltaT));
						Ad.InserirMatriz(t*R,tt*R,(t + 1)*R - 1,(tt + 1)*R - 1, &MSoma2, 0, 0);
					}
				}
			}
			else
			{
				for (int tt = 0; tt <= t; tt++)
				{
					MSoma1.InserirMatriz(0, 0, R - 1, R - 1, &Addd2, (t - tt)*R, 0);
					MSoma2.InserirMatriz(0, 0, R - 1, R - 1, &Ad, t*R, tt*R);
					MSoma2.MultiplicarPorEscalar(double (0.0036 * deltaT));
					MSoma1.SomarComMatriz(&MSoma2);
					Ad.InserirMatriz(t*R,tt*R,(t + 1)*R - 1,(tt + 1)*R - 1, &MSoma1, 0, 0);
				}
			}
		}
	}
	CMatrizEsparsa Av(R);
	CMatrizEsparsa Alin(R * T,n_a);
	CMatrizEsparsa AvNeg(R);
	AvNeg.MultiplicarPorEscalar( -1);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*R;
	int l = 0;
	int c = 0;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			Alin.InserirMatriz(l, c, l + R - 1,c + R - 1, &Av, 0, 0);
			for (int tt = 0; tt <= t; tt++)
				Alin.InserirMatriz(l, tt*(n_a / T) + R,l + R - 1, tt*(n_a / T) + 2*R - 1, &Ad, t*R, tt*R);
			if (t > 0)
				Alin.InserirMatriz(l, c - (n_a / T), l + R - 1, c - (n_a / T) + R - 1, &AvNeg, 0, 0);
		}
		else
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			//if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
			Alin.InserirMatriz(l, c, l + R - 1,c + R - 1, &Av, 0, 0);
			cc = 0;
			for (int tt = 0; tt <= t; tt++)
			{
				Alin.InserirMatriz(l, cc + R, l + R - 1, cc + 2*R - 1, &Ad, t*R, tt*R);
				if (((tt + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((tt + 1 - sistema_a->GetTt1()) > 0))
					cc += (n_a / T) + flag2*R;
				else
					cc += (n_a / T);
			}
			if (t > 0)
				if (t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) == sistema_a->GetTt1())
					Alin.InserirMatriz(l, (n_a / T)*sistema_a->GetTt1() - (n_a / T), l + R - 1,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + R - 1, &AvNeg, 0, 0);
				else
					Alin.InserirMatriz(l, c - (n_a / T), l + R - 1, c - (n_a / T) + R - 1, &AvNeg, 0, 0);
		}
		l = l + R;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*R;
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec3SPHA::MatrizVmeta(int n_a, int ch)
{
	int flag2 = int (sistema_a->GetFlagVfol());
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	size_t R = cascata[ch].size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n_a = n_a - flag2*sistema_a->GetNCenarios()*R;
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*R);
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*R);
	int c;
	//for (int cen = sistema_a->GetNCenarios() - 1 ; cen >= 0; cen--)
	//{
	//	for (size_t r = 0; r < R; r++)
	//	{	
	//		c = r + (sistema_a->GetNCenarios() - 1 - cen)*R;
	//		c = c + (T - 1)*(n_a / T) - cen * (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T);
	//		a.RemoverTodosElementos();
	//		a.InserirElemento(0, c, 1);
	//		if (sistema_a->GetFlagVfol() == true)
	//			a.InserirElemento(0, c + 2*R, 1);
	//		Alin.JuntarColuna(&a);
	//	}
	//}
	c = (sistema_a->GetTt2() - 1)*(n_a / T);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (size_t r = 0; r < R; r++)
		{	
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + r, 1);
			if (sistema_a->GetFlagVfol() == true)
				a.InserirElemento(0, c + r + 2*R, 1);
			Alin.JuntarColuna(&a);
		}
		c += (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T) + flag2*R;
	}
	return Alin;
}

// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Spcdec3SPHA::LimBalHid(int n_a, int ch)
{
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	double deltaT;
	int cenario;
	int periodo;
	size_t R = cascata[ch].size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0); 
	for (int t = 0; t < T; t++)
	{
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			deltaT = sistema_a->GetDeltaT1();
			periodo = t;
			cenario = 0;
		}
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				{
					deltaT = sistema_a->GetDeltaT2();
					periodo = t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1());
					cenario = cen;
				}
		for (int r = 0; r < R; r++)
		{
			if (t == 0)
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetV0() + 0.0036*deltaT*sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetAfluencia(cenario,periodo));
			}
			else
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (0.0036*deltaT*sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetAfluencia(cenario,periodo));
			}
		}
	}
	return Lim;
}
vetorfloat2 Spcdec3SPHA::LimVmeta(int n_a, int ch)
{
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	size_t R = cascata[ch].size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * sistema_a->GetNCenarios(), 2);
	IniciaMatriz(&Lim, 0);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (int r = 0; r < R; r++)
		{
			Lim[r + R*cen][1] = sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVMeta();
			//Lim(r,1) = sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax();
			Lim[r + R*cen][0] = 2;
		}
	}
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void Spcdec3SPHA::MatrizRestricoesLineares(int n_a, CMatrizEsparsa &MM, int ch)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	//int count = 0;
	M = MatrizBalHid(n_a, ch);
	//SparseMatriz(&M, indexL, indexC, indexV, nnZ, &n_restricoes, &count);
	//n_restricoes = n_restricoes + M.GetNlin(); M.ZerarMatriz();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizVmeta(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
}
void Spcdec3SPHA::MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int ch)
{
	vetorfloat2 L;
	L = LimBalHid(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimVmeta(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
}
// ------------------------------------------------

Spcdec3SPHA::Spcdec3SPHA(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;
	int flag2 = int (sistema_a->GetFlagVfol());
	n_cascatas = 0;
	vetorint usinas_fim;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		if ( sistema_a->hidreletricasVtr[r].GetUsinaJusante() == 0)
		{
			usinas_fim.push_back(r);
			n_cascatas++;
		}
	modelosGRB.resize(n_cascatas);
	vars.resize(n_cascatas);
	for (int ch = 0; ch < n_cascatas; ch++)
		modelosGRB[ch] = new GRBModel(ambiente_gurobi);		// um modelo de otimização para cada cascata
	n.resize(n_cascatas);
	cascata.resize(n_cascatas);
	N = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );

	int cc, ult_usina;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
	{
		cc = sistema_a->hidreletricasVtr[r].GetUsinaJusante();
		ult_usina = r;
		while ( cc != 0 )
		{
			ult_usina = cc - 1;
			cc = sistema_a->hidreletricasVtr[ult_usina].GetUsinaJusante();
		}
		for (int ch = 0; ch < n_cascatas; ch++)
		{
			if ( ult_usina == usinas_fim[ch] )
				cascata[ch].push_back(r);			// cada vetor cascata[ch] pode ter tamanho diferente !!
		}	
	}

	for (int ch = 0; ch < n_cascatas; ch++)
		n[ch] = 2 * N * cascata[ch].size() + flag2 * sistema_a->GetNCenarios() * cascata[ch].size();

	// Criar variáveis
	for (int ch = 0; ch < n_cascatas; ch++)
	{
		CriarVariaveis(ch);
		modelosGRB[ch]->update();	// Atualiza o modelo Gurobi.
		vars[ch] = modelosGRB[ch]->getVars();
	}
	
	// Adicionar restrições
	for (int ch = 0; ch < n_cascatas; ch++)
	{
		CriarRestricoes(ch);
		modelosGRB[ch]->update();	// Atualiza o modelo Gurobi.
	}

	// Definir ajustes do solver
	for (int ch = 0; ch < n_cascatas; ch++)
	{
		//modeloGRB.write("probHA.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modelosGRB[ch]->getEnv().set(GRB_IntParam_OutputFlag, 0);

		//modelosGRB[ch]->getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[ch]->getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		//modelosGRB[ch]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_Threads, 4);
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_PreSparsify, 1);
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_Method, 0);
		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);

		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_NodefileStart, 0.5);
		//modelosGRB[ch]->getEnv().set(GRB_StringParam_NodefileDir, "D:\Nodefile");
		
		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_IterationLimit, 100000);	// Limita o numero de iterações do método simplex
		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_NodefileStart, 0.5);
		
		modelosGRB[ch]->getEnv().set(GRB_DoubleParam_TimeLimit, 5);		// Limita o tempo de resolução do problema

		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-2);	// Define tolerancia de otimalidade
	}
}
Spcdec3SPHA::~Spcdec3SPHA(void)
{
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
}
void Spcdec3SPHA::SetPrecision(double precision)
{
	for (int ch = 0; ch < n_cascatas; ch++)
	{
		modelosGRB[ch]->getEnv().set(GRB_DoubleParam_OptimalityTol, precision);
	}
}
void Spcdec3SPHA::CriarVariaveis(int ch)
{
	try 
	{
		// variáveis em cada problema somente para uma cascata
		size_t R = sistema_a->hidreletricasVtr.size();
		// Estagio 1
		//x = [v d]
		int T = sistema_a->GetTt1();
		for (int t = 0; t < T; t++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)	//v
				modelosGRB[ch]->addVar(double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin(t)), double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax(t)), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < cascata[ch].size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
				modelosGRB[ch]->addVar(0, qhmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
		}
		// Estagio 2
		//x = [v d vfol]
		T = sistema_a->GetTt2() - sistema_a->GetTt1();
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (int t = 0; t < T; t++)
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)	//v
					modelosGRB[ch]->addVar(double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin(t+sistema_a->GetTt1()+cen*T)), double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax(t+sistema_a->GetTt1()+cen*T)), 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < cascata[ch].size(); r++)	//d
				{
					double qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
						qhmax = qhmax + double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
					modelosGRB[ch]->addVar(0, qhmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( sistema_a->GetFlagVfol() == true )
				for (size_t r = 0; r < cascata[ch].size(); r++)	//vfol
					//modelosGRB[ch]->addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
					modelosGRB[ch]->addVar(0, double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax() - sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin()), 0.0, GRB_CONTINUOUS, "");
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
void Spcdec3SPHA::CriarRestricoes(int ch)
{
	try
	{
		vetorint LimTipo;
		vetorfloat LimValor;
		LimTipo.resize(0);LimValor.resize(0);
		CMatrizEsparsa MM(0, n[ch]);

		MatrizRestricoesLineares(n[ch], MM, ch);
		MatrizLimitesLineares(n[ch], &LimTipo, &LimValor, ch);

		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		int l;
		for (int i = 0; i < MM.GetNlin(); i++)		// loop no numero de restrições
		{
			l = MM.GetValorLprim(i);
			while ( l != -1 )
			{
				coeficiente = MM.GetValorVal(l);
				variavel = vars[ch][MM.GetValorCol(l)];
				restricao.addTerms( &coeficiente, &variavel, 1);
				l = MM.GetValorLprox(l);
			}
			switch (LimTipo[i]) {
			case 0:
				modelosGRB[ch]->addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modelosGRB[ch]->addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modelosGRB[ch]->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
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
void Spcdec3SPHA::CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int ch)
{
	try 
	{
		//x = [v d vfol]
		int flag2 = int (sistema_a->GetFlagVfol());

		// Termos lineares
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int n_a = n[ch] - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
		int nt;
		int nt_dual = lambda->size() / T;
		int delta_dual = 0;
		double deltaT;
		int delta = 0;
		int cen = 0;
		// vfolga
		int JJ = 0;
		for (size_t i = 0; i < cascata[ch].size(); i++)
			JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
		delta = 2*cascata[ch].size();
		for (int t = 0; t < T; t++)
		{
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*cascata[ch].size();
			else
				nt = (n_a / T);
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{

			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				//	if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				//	{
				deltaT = sistema_a->GetDeltaT2();
				if (sistema_a->GetFlagVfol() == true)
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)
					{
						if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
							vars[ch][delta + r ].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoVfol()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//vfol
					}
				}
			}
			delta = delta + nt;
		}
		
		// Termos da RL
		//x = [v d vfol]
		// lambda [pt ph v d (phmax || res)] -> tamanho: [I R R R (R || 1)] x T
		// Constantes

		// Lineares
		delta = 0;
		delta_dual = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)
			{
				vars[ch][r + delta].set(GRB_DoubleAttr_Obj, lambda->at(cascata[ch].at(r) + sistema_a->termeletricasVtr.size() + 1*sistema_a->hidreletricasVtr.size() + delta_dual));							//v
				vars[ch][r + 1*cascata[ch].size() + delta].set(GRB_DoubleAttr_Obj, lambda->at(cascata[ch].at(r) + sistema_a->termeletricasVtr.size() + 2*sistema_a->hidreletricasVtr.size() + delta_dual));		//d
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*cascata[ch].size();
			else
				nt = (n_a / T);
			delta += nt;
			delta_dual += nt_dual;
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
int Spcdec3SPHA::ResolverProblemaRL(Spcdec3Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int ch)
{
	try
	{
		// Criar função objetivo, a seleçao dos lambdas de interesse é feita dentro da funçao abaixo
		CriarFuncaoObjetivoRL(lambda, ch);
		modelosGRB[ch]->update();	// Atualiza o modelo Gurobi.

		modelosGRB[ch]->reset();

		// Otimizar
		modelosGRB[ch]->optimize();

		// Salvar resultados
		x_ch.clear();
		x_ch.resize(n[ch]);
		int nStatus = modelosGRB[ch]->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The HA subproblem was not solved until optimality! " << nStatus << endl;
		if ((nStatus == 2) || (nStatus == 9))
		{
			fo_ch = double(modelosGRB[ch]->get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n[ch]; i++)
				x_ch[i] = double(vars[ch][i].get(GRB_DoubleAttr_X));
		}
		else
		{
			for (int i = 0; i < n[ch]; i++)
				x_ch[i] = 0;
			fo_ch = 0;
		}
		resultadoGurobi->GravarSolucao(fo_ch, x_ch, nStatus, ch);

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->GravarSolucao(fo_ch, x_ch, e.getErrorCode(), ch);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}

// Caso os subproblemas sejam Easy Components

Spcdec3SPHA::Spcdec3SPHA(CSistema * const sistema_end)
{
	sistema_a = sistema_end;
	int flag2 = int (sistema_a->GetFlagVfol());
	n_cascatas = 0;
	vetorint usinas_fim;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		if ( sistema_a->hidreletricasVtr[r].GetUsinaJusante() == 0)
		{
			usinas_fim.push_back(r);
			n_cascatas++;
		}
	n.resize(n_cascatas);
	cascata.resize(n_cascatas);
	N = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );

	int cc, ult_usina;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
	{
		cc = sistema_a->hidreletricasVtr[r].GetUsinaJusante();
		ult_usina = r;
		while ( cc != 0 )
		{
			ult_usina = cc - 1;
			cc = sistema_a->hidreletricasVtr[ult_usina].GetUsinaJusante();
		}
		for (int ch = 0; ch < n_cascatas; ch++)
		{
			if ( ult_usina == usinas_fim[ch] )
				cascata[ch].push_back(r);			// cada vetor cascata[ch] pode ter tamanho diferente !!
		}	
	}

	for (int ch = 0; ch < n_cascatas; ch++)
		n[ch] = 2 * N * cascata[ch].size() + flag2 * sistema_a->GetNCenarios() * cascata[ch].size();

	Mind.resize(n_cascatas);
	Mval.resize(n_cascatas);
	Mbeg.resize(n_cascatas);
	dtipo.resize(n_cascatas);
	dvalor.resize(n_cascatas);
	Nrestricoes.resize(n_cascatas);
	//var_ub.resize(n_cascatas);
	//var_lb.resize(n_cascatas);

	for (int ch = 0; ch < n_cascatas; ch++)
	{
		//GerarVarBounds(ch);
		GerarRestricoes(ch);
	}
}
void Spcdec3SPHA::GerarRestricoes(int ch)
{
	Mind[ch].resize(0);Mval[ch].resize(0);Mbeg[ch].resize(0);dtipo[ch].resize(0);dvalor[ch].resize(0);
	CMatrizEsparsa MM(0, n[ch]);
	MatrizRestricoesLineares(n[ch], MM, ch);
	MM.SparseMatriz(Mval[ch], Mind[ch], Mbeg[ch]);
	MatrizLimitesLineares(n[ch], &dtipo[ch], &dvalor[ch], ch);
	Nrestricoes[ch] = MM.GetNlin();
}
void Spcdec3SPHA::GerarVarBounds(int ch, double *lbd , double *ubd)
{
	// variáveis em cada problema somente para uma cascata
	size_t R = sistema_a->hidreletricasVtr.size();
	// Estagio 1
	//x = [v d]
	int iub = 0;
	int ilb = 0;
	int T = sistema_a->GetTt1();
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < cascata[ch].size(); r++)	//v
		{
			ubd[iub++] = double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax(t));
			lbd[ilb++] = double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin(t));
		}
		for (size_t r = 0; r < cascata[ch].size(); r++)	//d
		{
			double qhmax = 0;
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
				qhmax = qhmax + double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
			ubd[iub++] = qhmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax();
			lbd[ilb++] = 0;
		}
	}
	// Estagio 2
	//x = [v d vfol]
	T = sistema_a->GetTt2() - sistema_a->GetTt1();
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (int t = 0; t < T; t++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)	//v
			{
				ubd[iub++] = double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax(t+sistema_a->GetTt1()+cen*T));
				lbd[ilb++] = double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin(t+sistema_a->GetTt1()+cen*T));
			}
			for (size_t r = 0; r < cascata[ch].size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
				ubd[iub++] = qhmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax();
				lbd[ilb++] = 0;
			}
		}
		if ( sistema_a->GetFlagVfol() == true )
			for (size_t r = 0; r < cascata[ch].size(); r++)	//vfol
			{
				ubd[iub++] = double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax() - sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin());
				lbd[ilb++] =  0;
			}
	}
}
void Spcdec3SPHA::GerarCoefFuncaoObjetivo(int ch, double *cst)
{
	//x = [v d vfol]
	// Termos lineares
	int delta = 0;
	for (int i = 0; i < n[ch]; i++)		// zerar todos elementos de cst
		cst[i] = 0;

	if (sistema_a->GetFlagVfol() == true)		// atribuir valores somente para vfol
	{
		delta = 2*cascata[ch].size() * sistema_a->GetTt2();
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)
			{
				cst[delta + r ] = - sistema_a->GetCustoVfol()* sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
			}
			delta += 2*cascata[ch].size() * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[ch].size();
		}
	}
}
void Spcdec3SPHA::GetBDesc(int comp , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd)
{
	// fazer uma funçao GetBDesc nessa classe, que passa os valores da matriz jah criada no constructor para os respectivos vetores!!
	// The NDO solver is allowed to query only a subset of the information;
    // that is, by passing 0 to some of the pointers, the oracle will not write
    // the corresponding information
	if ((Bbeg != 0) && (Bind !=0) && (Bval !=0))
	{
		for (size_t i = 0; i < Mbeg[comp].size(); i++)
			Bbeg[i] = Mbeg[comp][i];
		for (size_t i = 0; i < Mind[comp].size(); i++)
		{
			Bind[i] = Mind[comp][i];
			Bval[i] = Mval[comp][i];
		}
	}

	//for (int i = 0; i < Nrestricoes[comp]; i++)			// loop no numero de restriçoes (=num. linhas da matriz)
	//{
	//	switch (dtipo[comp][i]) {
	//	case 0:
	//		rhs[i] = dvalor[comp][i];
	//		lhs[i] = - GRB_INFINITY;
	//		break;
	//	case 1:
	//		rhs[i] = dvalor[comp][i];
	//		lhs[i] = dvalor[comp][i];
	//		break;
	//	case 2:
	//		rhs[i] = GRB_INFINITY;
	//		lhs[i] = dvalor[comp][i];
	//		break;
	//	default:
	//		cout << "Tipo inválido de restrição adicionada" << endl;
	//	}
	//}

	if (rhs != 0)
	{
		for (int i = 0; i < Nrestricoes[comp]; i++)			// loop no numero de restriçoes (=num. linhas da matriz)
		{
			switch (dtipo[comp][i]) {
			case 0:
				rhs[i] = dvalor[comp][i];
				break;
			case 1:
				rhs[i] = dvalor[comp][i];
				break;
			case 2:
				rhs[i] = OPTtypes_di_unipi_it::Inf<double>();
				break;
			default:
				cout << "Tipo inválido de restrição adicionada" << endl;
			}
		}
	}
	if (lhs !=0)
	{
		for (int i = 0; i < Nrestricoes[comp]; i++)			// loop no numero de restriçoes (=num. linhas da matriz)
		{
			switch (dtipo[comp][i]) {
			case 0:
				lhs[i] = - OPTtypes_di_unipi_it::Inf<double>();
				break;
			case 1:
				lhs[i] = dvalor[comp][i];
				break;
			case 2:
				lhs[i] = dvalor[comp][i];
				break;
			default:
				cout << "Tipo inválido de restrição adicionada" << endl;
			}
		}
	}

	if ((lbd != 0) && (ubd != 0))
		GerarVarBounds(comp, lbd, ubd);				// definir os limites das variaveis
	
	if (cst != 0)
		GerarCoefFuncaoObjetivo(comp, cst);			// definir a funçao objetivo (inverter sinal pois o prob. no MP é de max.

	// Escrever em arquivo
	//ofstream inFile( "matrizB.txt", ios::out );   
	//if ( inFile.is_open() )
	//{
	//	for ( size_t l = 0; l < Nrestricoes[comp].size(); l++)
	//	{
	//		inFile << Bcoluna[comp][l] << "\t" << Blinha[comp][l] << "\t" << Bvalor[comp][l] << endl;
	//	}
	//inFile.close();
	//}
	//else
	//	cout << "Unable to open file";

	//ofstream inFile2( "matrizBesp.txt", ios::out );   
	//if ( inFile2.is_open() )
	//{
	//	for ( size_t l = 0; l < n[comp] + 1; l++)
	//		inFile2 << Bbeg[l] << "\t"; 
	//	inFile2 << endl;
	//	for ( size_t l = 0; l < Bcoluna[comp].size(); l++)
	//		inFile2 << Bind[l] << "\t"; 
	//	inFile2 << endl;
	//	for ( size_t l = 0; l < Bcoluna[comp].size(); l++)
	//		inFile2 << Bval[l] << "\t"; 
	//	inFile2 << endl;
	//	inFile2.close();
	//}
	//else
	//	cout << "Unable to open file";

	//ofstream inFile( "matrizB.txt", ios::out );   
	//if ( inFile.is_open() )
	//{
	//	for ( size_t l = 0; l < Nrestricoes[comp]; l++)
	//	{
	//		inFile << lhs[l] << "\t" << rhs[l] << endl;
	//	}
	//inFile.close();
	//}
	//else
	//	cout << "Unable to open file";

	//ofstream inFile2( "matrizBesp.txt", ios::out );   
	//if ( inFile2.is_open() )
	//{
	//	for ( size_t l = 0; l < n[comp]; l++)
	//		inFile2 << cst[l] << "\t" << lbd[l] << "\t" << ubd[l] << endl; 
	//	
	//	inFile2.close();
	//}
	//else
	//	cout << "Unable to open file";

	// Fim (Escreve em arquivo)

}
int Spcdec3SPHA::GetANZ(int comp)
{
	// Retorna o numero de multiplicadores de lagrange considerado em cada subproblema
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	return ( 2 * cascata[comp].size() * T );

	// funçao chamada para cada componente!! qq lambda poderia influenciar na componente
	// deve-se retornar um valor somente para os lambdas q estao em cada subproblema
	// só retornar um valor quando o strt and stp corresponderem ao intervalo do subproblema???

	// com o numero do componente sabe-se quais lambdas sao usados!! mas pq usar strt e stp???

	// as variáveis duais precisam estarem em ordem por cada subproblema???
	//loop de strt até stp conferindo se a var. dual pertence ao subproblema em questao...

}
void Spcdec3SPHA::GetADesc( int comp , int *Abeg , int *Aind , double *Aval )
{
	// Termos da RL
	// lambda [pt ph v d (phmax)]
	// x = [v d vfol]
	int flag3 = int (sistema_a->GetFlagPhmax());
	int nt_dual = sistema_a->termeletricasVtr.size() + (3 + flag3) * sistema_a->hidreletricasVtr.size() + (1 - flag3);
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	int ind = 0;
	int iAbeg = 0;
	Abeg[iAbeg++] = 0;
	int delta_dual = 0;
	for (int t = 0; t < T; t++)		// loop para todas as variáris do subproblema (BNC), deve ser na ordem das colunas, variáveis...
	{
		for (size_t r = 0; r < cascata[comp].size(); r++)		// v
		{
			Aind[ind] = cascata[comp].at(r) + sistema_a->termeletricasVtr.size() + 1*sistema_a->hidreletricasVtr.size() + delta_dual;		// numero da linha correspondente ao lambda
			//Aval[ind++] = 1;
			Aval[ind++] = sistema_a->GetPrecondidioner(cascata[comp].at(r) + sistema_a->termeletricasVtr.size() + 1*sistema_a->hidreletricasVtr.size() + delta_dual);
			Abeg[iAbeg++] = ind;
		}
		for (size_t r = 0; r < cascata[comp].size(); r++)		// d
		{
			Aind[ind] = cascata[comp].at(r) + sistema_a->termeletricasVtr.size() + 2*sistema_a->hidreletricasVtr.size() + delta_dual;
			//Aval[ind++] = 1;
			Aval[ind++] = sistema_a->GetPrecondidioner(cascata[comp].at(r) + sistema_a->termeletricasVtr.size() + 2*sistema_a->hidreletricasVtr.size() + delta_dual);
			Abeg[iAbeg++] = ind;
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		{
			if (sistema_a->GetFlagVfol() == true)		//vfol (n tem lambda na funçao objetivo)
				for (size_t r = 0; r < cascata[comp].size(); r++)
					Abeg[iAbeg++] = ind;		// adiciona-se um elemento em Abeg porem com o valor dos anteriores, indicando que n teve-se adição de elementos n nulos
		}
		delta_dual += nt_dual;
	}

	//// Escrever em arquivo
	//ofstream inFile2( "matrizBesp.txt", ios::out );   
	//if ( inFile2.is_open() )
	//{
	//	for ( size_t l = 0; l < n[comp] + 1; l++)
	//		inFile2 << Abeg[l] << "\t"; 
	//	inFile2 << endl;
	//	for ( size_t l = 0; l < 2 * cascata[comp].size() * T; l++)
	//		inFile2 << Aind[l] << "\t"; 
	//	inFile2 << endl;
	//	for ( size_t l = 0; l < 2 * cascata[comp].size() * T; l++)
	//		inFile2 << Aval[l] << "\t"; 
	//	inFile2 << endl;
	//	inFile2.close();
	//}
	//else
	//	cout << "Unable to open file";
	//// Fim (Escreve em arquivo)
}
int Spcdec3SPHA::GetNNZsRest(int comp)
{
	return ( int(Mind[comp].size()) );
}
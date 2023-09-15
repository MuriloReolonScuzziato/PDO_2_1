#include "SubProblemaCen.h"

// Modelo de otimização de cada cenário
namespace met_DecCen
{
// Funções para criar as matrizes de restrições
// ------------------------------------------------
CMatrizEsparsa MatrizRestDemanda(CSistema * const sistema_end, int n_a)	// Monta matriz da restrição de atendimento a demanda
{
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	CMatrizEsparsa Agh(sistema_end->barrasPtr.size(), sistema_end->hidreletricasPtr.size());
	for ( size_t b = 0; b < sistema_end->barrasPtr.size(); b++)
		for ( size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
			for ( size_t br = 0; br < sistema_end->barrasPtr[b]->hidrosPtr.size(); br++)
				if (sistema_end->barrasPtr[b]->hidrosPtr[br] == sistema_end->hidreletricasPtr[r])
					Agh.InserirElemento(b, r, 1);
	CMatrizEsparsa Agt(sistema_end->barrasPtr.size(), sistema_end->termeletricasPtr.size());
	for ( size_t b = 0; b < sistema_end->barrasPtr.size(); b++)
		for ( size_t i = 0; i < sistema_end->termeletricasPtr.size(); i++)
			for ( size_t bi = 0; bi < sistema_end->barrasPtr[b]->termosPtr.size(); bi++)
				if (sistema_end->barrasPtr[b]->termosPtr[bi] == sistema_end->termeletricasPtr[i])
					Agt.InserirElemento(b, i, 1);
	CMatrizEsparsa B(sistema_end->barrasPtr.size(), sistema_end->barrasPtr.size());
	for ( size_t b = 0; b < sistema_end->barrasPtr.size(); b++)
		for ( size_t bb = 0; bb < sistema_end->barrasPtr.size(); bb++)
			for ( size_t l = 0; l < sistema_end->linhasPtr.size(); l++)
				if (b == bb)
				{
					if ((sistema_end->linhasPtr[l]->de_barra == sistema_end->barrasPtr[b]) || (sistema_end->linhasPtr[l]->para_barra == sistema_end->barrasPtr[b]))
						if ( B.GetElemento(b, bb) != 0 )
							B.SubstituirElemento(b, bb, B.GetElemento(b, bb) + 100/(sistema_end->linhasPtr[l]->GetReatancia()));
						else
							B.InserirElemento(b, bb, 100/(sistema_end->linhasPtr[l]->GetReatancia()));
				}
				else
				{
					if (((sistema_end->linhasPtr[l]->de_barra == sistema_end->barrasPtr[b]) && (sistema_end->linhasPtr[l]->para_barra == sistema_end->barrasPtr[bb])) ||
						((sistema_end->linhasPtr[l]->de_barra == sistema_end->barrasPtr[bb]) && (sistema_end->linhasPtr[l]->para_barra == sistema_end->barrasPtr[b])))
						if ( B.GetElemento(b, bb) != 0 )
							B.SubstituirElemento(b, bb, B.GetElemento(b, bb) - 100/(sistema_end->linhasPtr[l]->GetReatancia()));
						else
							B.InserirElemento(b, bb, - 100/(sistema_end->linhasPtr[l]->GetReatancia()));
				}
	B.RemoverColuna(sistema_end->GetBarraRef() - 1);
	B.MultiplicarPorEscalar( -1);
	CMatrizEsparsa Alin(T * sistema_end->barrasPtr.size(), n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int l = 0;
	int c = 0;
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c, l + sistema_end->barrasPtr.size() - 1, c + sistema_end->termeletricasPtr.size() - 1, &Agt, 0, 0);
		Alin.InserirMatriz(l, c + ((3+flag4)*sistema_end->termeletricasPtr.size()),l + sistema_end->barrasPtr.size() - 1, c + ((3+flag4)*sistema_end->termeletricasPtr.size() + sistema_end->barrasPtr.size()) - 2, &B, 0, 0);
		Alin.InserirMatriz(l, c + (sistema_end->barrasPtr.size() - 1 + (3+flag4)*sistema_end->termeletricasPtr.size()), l + sistema_end->barrasPtr.size() - 1, c + (sistema_end->barrasPtr.size() - 1 + (3+flag4)*sistema_end->termeletricasPtr.size() + sistema_end->hidreletricasPtr.size()) - 1, &Agh, 0, 0);
		l = l + sistema_end->barrasPtr.size();
		c = c + (n_a / T);
	}
	int JJ = 0;
	for (size_t i = 0; i < sistema_end->hidreletricasPtr.size(); i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	c = (3+flag4)*(sistema_end->termeletricasPtr.size()) + sistema_end->barrasPtr.size() - 1 + 5*(sistema_end->hidreletricasPtr.size()) + 3*JJ;
	l = 0;
	CMatrizEsparsa Adef(sistema_end->barrasPtr.size());
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c,l + sistema_end->barrasPtr.size() - 1,c + sistema_end->barrasPtr.size() - 1, &Adef, 0, 0);
		l = l + sistema_end->barrasPtr.size();
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizRestDemandaBarraUnica(CSistema * const sistema_end, int n_a)	// Monta matriz da restrição de atendimento a demanda
{
	int T = sistema_end->GetTt2();
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	size_t I = sistema_end->termeletricasPtr.size();
	size_t R = sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (size_t i = 0; i < sistema_end->hidreletricasPtr.size(); i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	CMatrizEsparsa Alin(0, n_a );
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c = 0;
	for ( int t = 0; t < T; t++)
	{
		a.RemoverTodosElementos();
		for (size_t i = 0; i < I; i++)		// pt
			a.InserirElemento(0, i + c, 1);
		for (size_t r = 0; r < R; r++)		// ph
			a.InserirElemento(0, r + c + (3+flag4)*sistema_end->termeletricasPtr.size(), 1);
		a.InserirElemento(0, (3+flag4)*sistema_end->termeletricasPtr.size() + 5*sistema_end->hidreletricasPtr.size() + 3*JJ + c, 1);		// def
		Alin.JuntarColuna(&a);
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizLimFluxo(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	CMatrizEsparsa Alb(sistema_end->linhasPtr.size(),sistema_end->barrasPtr.size());
	for ( size_t l = 0; l < sistema_end->linhasPtr.size(); l++)
		for ( size_t b = 0; b < sistema_end->barrasPtr.size(); b++)
			if (sistema_end->linhasPtr[l]->de_barra == sistema_end->barrasPtr[b])
				Alb.InserirElemento(l, b, 1);
			else if (sistema_end->linhasPtr[l]->para_barra == sistema_end->barrasPtr[b])
				Alb.InserirElemento(l, b, -1);
	Alb.RemoverColuna(sistema_end->GetBarraRef() - 1);
	CMatrizEsparsa TT(sistema_end->linhasPtr.size(),sistema_end->linhasPtr.size());
	for ( size_t l = 0; l < sistema_end->linhasPtr.size(); l++)
	{
		for ( size_t ll = 0; ll < sistema_end->linhasPtr.size(); ll++)
		{
			if (l == ll)
			{
				TT.InserirElemento(l, ll, 100 / (sistema_end->linhasPtr[l]->GetReatancia()));
			}
		}
	}
	CMatrizEsparsa Alin(T * sistema_end->linhasPtr.size(),n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(ll,c + ((3+flag4)*sistema_end->termeletricasPtr.size()),ll + sistema_end->linhasPtr.size() - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + sistema_end->barrasPtr.size()) - 2, &TT, 0, 0);
		ll = ll + sistema_end->linhasPtr.size();
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizLimPhgMin(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	int JJ = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
		for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetPmin());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ,n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()),l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizLimPhgMax(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	int JJ = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
		for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetPmax() * sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ,n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()),l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizBalHid(CSistema * const sistema_end, int n_a)		// TESTAR ESSA FUNÇÃO PARA DOIS ESTÁGIOS... DEMAIS RESTRIÇÕES PARA DELTAt1 != DELTAt2
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	double deltaT;
	int tempo_viagem;
	int R = sistema_end->hidreletricasPtr.size();
	CMatrizEsparsa Addd(R*T, R);
	CMatrizEsparsa Add(R, R);
	//EyeMatriz(&Addd);
	for (int t = 0; t < T; t++)
	{
		deltaT = sistema_end->GetDeltaT1();
		Add.RemoverTodosElementos();
		for (int r = 0; r < R; r++)
		{
			if (sistema_end->hidreletricasPtr[r]->GetTempoViagem() < deltaT)
				tempo_viagem = 0;
			else
				tempo_viagem = int (sistema_end->hidreletricasPtr[r]->GetTempoViagem() / deltaT);
		if ((tempo_viagem == t) && (sistema_end->hidreletricasPtr[r]->GetUsinaJusante() != 0))
			{
				Add.InserirElemento(sistema_end->hidreletricasPtr[r]->GetUsinaJusante() - 1, r, double (-0.0036 * deltaT));
			}
		}
		Addd.InserirMatriz(t * R,0,(t + 1) * R - 1,R - 1, &Add, 0, 0);
	}
	CMatrizEsparsa Addd2(R * T,R);
	//EyeMatriz(&Addd2);
	for (int t = 0; t < T; t++)
	{
		deltaT = sistema_end->GetDeltaT2();
		Add.RemoverTodosElementos();
		for (int r = 0; r < R; r++)
		{
			if (sistema_end->hidreletricasPtr[r]->GetTempoViagem() < deltaT)
				tempo_viagem = 0;
			else
				tempo_viagem = int (sistema_end->hidreletricasPtr[r]->GetTempoViagem() / deltaT);
			if ((tempo_viagem == t) && (sistema_end->hidreletricasPtr[r]->GetUsinaJusante() != 0))
			{
				Add.InserirElemento(sistema_end->hidreletricasPtr[r]->GetUsinaJusante() - 1, r, double (-0.0036 * deltaT));
			}
		}
		Addd2.InserirMatriz(t * R,0,(t + 1) * R - 1,R - 1, &Add, 0, 0);
	}
	CMatrizEsparsa Ad(R*T);
	CMatrizEsparsa MSoma1(R, R);
	CMatrizEsparsa MSoma2(R, R);
	for (int t = 0; t < T; t++)
	{
		MSoma1.RemoverTodosElementos();
		MSoma2.RemoverTodosElementos();
		if (t < sistema_end->GetTt1())
		{
			deltaT = sistema_end->GetDeltaT1();
			for (int tt = 0; tt <= t; tt++)
			{
				MSoma1.InserirMatriz(0, 0, R - 1, R - 1, &Addd, (t - tt)*R, 0);
				MSoma2.InserirMatriz(0, 0, R - 1, R - 1, &Ad, t*R, tt*R);
				MSoma2.MultiplicarPorEscalar(double (0.0036 * deltaT));
				MSoma1.SomarComMatriz(&MSoma2);
				Ad.InserirMatriz(t*R,tt*R,(t + 1)*R - 1,(tt + 1)*R - 1, &MSoma1, 0, 0);
			}	
		}
		else
		{
			deltaT = sistema_end->GetDeltaT2();
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
	CMatrizEsparsa Av(R);
	CMatrizEsparsa Alin(R*T, n_a);
	CMatrizEsparsa AvNeg(R);
	AvNeg.MultiplicarPorEscalar( -1);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + R),l + R - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 2*R) - 1, &Av, 0, 0);
		for (int tt = 0; tt <= t; tt++)
			Alin.InserirMatriz(l,tt*(n_a / T) + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 2*R),l + R - 1,tt*(n_a / T) + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 3*R) - 1, &Ad, t*R, tt*R);
		if (t > 0)
			Alin.InserirMatriz(l,c - (n_a / T) + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + R),l + R - 1,c - (n_a / T) + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
		l = l + R;
		c = c + (n_a / T);
	}
	//Alin.ImprimirMatrizArquivo();
	return Alin;
}
CMatrizEsparsa MatrizLimqMin(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	int JJ = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
		for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetQmin());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ,n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + JJ,l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizLimqMax(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	int JJ = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_end->hidreletricasPtr.size(); r++)
		for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetQmax() * sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ,n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + JJ,l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + 5*sistema_end->hidreletricasPtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa MatrizTup(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c = 0;
	for (size_t i = 0; i < I; i++)
	{
		c = I + i;
		for (int t = 0; t < T; t++)
		{
			if (sistema_end->termeletricasPtr[i]->GetTUp() > 0)
			{
				for (int tt = 0; tt < sistema_end->termeletricasPtr[i]->GetTUp(); tt++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c, 1);
					if (int(t) - tt - 1 < 0)		// cada if é adiciona um elemento dos dois elementos do lado direito da equação
					{
					}
					else
					{
						a.InserirElemento(0, c - tt*(n_a / T) - 1*(n_a / T), -1);
					}
					if (int(t) - tt - 2 < 0)
					{
					}
					else
					{
						a.InserirElemento(0, c - tt*(n_a / T) - 2*(n_a / T), 1);
					}
					Alin.JuntarColuna(&a);
				}
			}
			else
			{
			}
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizTdown(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c = 0;
	for (size_t i = 0; i < I; i++)
	{
		c = I + i;
		for (int t = 0; t < T; t++)
		{
			if (sistema_end->termeletricasPtr[i]->GetTDown() > 0)
			{
				for (int tt = 0; tt < sistema_end->termeletricasPtr[i]->GetTDown(); tt++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c, 1);
					if (int(t) - tt - 1 < 0)		// cada if é adiciona um elemento dos dois elementos do lado direito da equação
					{
					}
					else
					{
						a.InserirElemento(0, c - tt*(n_a / T) - 1*(n_a / T), -1);
					}
					if (int(t) - tt - 2 < 0)
					{
					}
					else
					{
						a.InserirElemento(0, c - tt*(n_a / T) - 2*(n_a / T), 1);
					}
					Alin.JuntarColuna(&a);
				}
			}
			else
			{
			}
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizRampaUp(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c;
	for (size_t i = 0; i < I; i++)
	{	
		c = i;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			if (t == 0)
			{
				a.InserirElemento(0, c, 1);
			}
			else
			{
				a.InserirElemento(0, c, 1);
				a.InserirElemento(0, c - (n_a / T), -1);
				a.InserirElemento(0, c + I - (n_a / T), sistema_end->termeletricasPtr[i]->GetPmin() - sistema_end->termeletricasPtr[i]->GetRampaUp());
			}
			Alin.JuntarColuna(&a);
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizRampaDown(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c;
	for (size_t i = 0; i < I; i++)
	{	
		c = i;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			if (t == 0)
			{
				a.InserirElemento(0, c, -1);
				a.InserirElemento(0, c + I, sistema_end->termeletricasPtr[i]->GetPmin() - sistema_end->termeletricasPtr[i]->GetRampaDown());
			}
			else
			{
				a.InserirElemento(0, c, -1);
				a.InserirElemento(0, c - (n_a / T), 1);
				a.InserirElemento(0, c + I, sistema_end->termeletricasPtr[i]->GetPmin() - sistema_end->termeletricasPtr[i]->GetRampaDown());
			}
			Alin.JuntarColuna(&a);
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizLimPtMin(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c;
	for (size_t i = 0; i < I; i++)
	{	
		c = i;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + I, - sistema_end->termeletricasPtr[i]->GetPmin());
			Alin.JuntarColuna(&a);
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizLimPtMax(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c;
	for (size_t i = 0; i < I; i++)
	{	
		c = i;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + I, - sistema_end->termeletricasPtr[i]->GetPmax());
			Alin.JuntarColuna(&a);
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizRestCP(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c;
	for (size_t i = 0; i < I; i++)
	{	
		c = i + I;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			if (t == 0)
			{
				a.InserirElemento(0, c + I, 1);
				a.InserirElemento(0, c, - sistema_end->termeletricasPtr[i]->GetCoefCustoPartida());
			}
			else
			{
				a.InserirElemento(0, c + I, 1);
				a.InserirElemento(0, c, - sistema_end->termeletricasPtr[i]->GetCoefCustoPartida());
				a.InserirElemento(0, c - (n_a / T), sistema_end->termeletricasPtr[i]->GetCoefCustoPartida());
			}
			Alin.JuntarColuna(&a);
			c = c + (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa MatrizVmeta(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	int c;
	for (size_t r = 0; r < R; r++)
	{	
		c = r + (3+flag4)*sistema_end->termeletricasPtr.size() + flag1*(sistema_end->barrasPtr.size() - 1) + R;
		c = c + (T - 1)*(n_a / T);
		a.RemoverTodosElementos();
		a.InserirElemento(0, c, 1);
		a.InserirElemento(0, c + 4*R + 3*JJ + flag1*sistema_end->barrasPtr.size() + (1 - flag1), 1);
		Alin.JuntarColuna(&a);
	}
	return Alin;
}
CMatrizEsparsa MatrizBalPotencia(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	size_t I = sistema_end->termeletricasPtr.size();
	size_t B = flag1*(sistema_end->barrasPtr.size() - 1);
	int nt = (n_a - sistema_end->hidreletricasPtr.size())/(T);
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int dd = (3+flag4)*I + B;
	int jj = dd + 5*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, dd + r, 1);
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, -1);
			}
			Alin.JuntarColuna(&a);
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
		dd = dd + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa MatrizBalVazao(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	size_t I = sistema_end->termeletricasPtr.size();
	size_t B = flag1*(sistema_end->barrasPtr.size() - 1);
	int nt = (n_a - sistema_end->hidreletricasPtr.size())/(T);
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int dd = (3+flag4)*I + B + 2*R;
	int jj = dd + JJ + 3*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, dd + r, 1);
			a.InserirElemento(0, dd + R + r, -1);
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, -1);
			}
			Alin.JuntarColuna(&a);
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
		dd = dd + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa MatrizFuncProd(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	size_t B = flag1*(sistema_end->barrasPtr.size() - 1);
	size_t I = sistema_end->termeletricasPtr.size();
	int nt = (n_a - sistema_end->hidreletricasPtr.size())/(T);
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int vv = (3+flag4)*I + B + R;
	int jj = vv + 4*R + JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R ; r++)
		{	
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, jj - JJ + j, 1);															// phg
				a.InserirElemento(0, vv + r, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->CoefPhgV());		// v
				a.InserirElemento(0, jj + j, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->CoefPhgQ());		// q
				a.InserirElemento(0, vv + R + r, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->CoefPhgD());	// d
				if (sistema_end->hidreletricasPtr[r]->GetInflueVert() == 1 )		// Se o s influencia já tenho contabilizado isso no d, caso contrário só subtrair o valor de s de d
				{
				}
				else		// vertimento n influencia no canl de fuga, entao subtraio o valor de s da defluencia
				{
					a.InserirElemento(0, vv + 2*R + r, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->CoefPhgS());	// s
				}
				Alin.JuntarColuna(&a);
			}
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
		vv = vv + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa MatrizPhMax(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	size_t I = sistema_end->termeletricasPtr.size();
	size_t B = flag1*(sistema_end->barrasPtr.size() - 1);
	int nt = (n_a - sistema_end->hidreletricasPtr.size())/(T);
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int dd = (3+flag4)*I + B + 4*R;
	int jj = dd + R + 2*JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, dd + r, 1);
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, - sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetPmax() * sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades());
			}
			Alin.JuntarColuna(&a);
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
		dd = dd + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa MatrizReserva(CSistema * const sistema_end, int n_a)
{
	int flag1 = int (1 - sistema_end->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	size_t I = sistema_end->termeletricasPtr.size();
	size_t B = flag1*(sistema_end->barrasPtr.size() - 1);
	int nt = (n_a - sistema_end->hidreletricasPtr.size())/(T);
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int dd = (3+flag4)*I + B + 4*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		a.RemoverTodosElementos();
		for (size_t r = 0; r < R; r++)
		{
			a.InserirElemento(0, dd + r, 1);
			a.InserirElemento(0, dd - 4*R + r, - 1);
		}
		Alin.JuntarColuna(&a);
		dd = dd + nt;
	}
	return Alin;
}
CMatrizEsparsa MatrizCortesF(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	n_a = n_a - sistema_end->hidreletricasPtr.size();
	int c;
	for (size_t i = 0; i < I; i++)
	{	
		c = i;
		for (int t = 0; t < T; t++)
		{
			for (int n_cort = 0; n_cort < sistema_end->GetFlagAproxCustoT(); n_cort++)		// Adiciona uma restrição por corte
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c, - sistema_end->termeletricasPtr[i]->GetCoefA1(n_cort));
				a.InserirElemento(0, c + I, - sistema_end->termeletricasPtr[i]->GetCoefA0(n_cort));
				a.InserirElemento(0, c + 3*I, 1);
				Alin.JuntarColuna(&a);
			}
			c = c + (n_a / T);
		}
	}
	return Alin;
}
// ------------------------------------------------

// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
// ------------------------------------------------
vetorfloat2 LimRestDemanda(CSistema * const sistema_end, int n_a, const int cen)
{
	int T = sistema_end->GetTt2();
	size_t B = sistema_end->barrasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, B * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
		for (size_t b = 0; b < B; b++)
		{
			double D = 0;
			for (size_t d = 0; d < sistema_end->barrasPtr[b]->demandasPtr.size(); d++)
				D = D + sistema_end->barrasPtr[b]->demandasPtr[d]->GetD(cen, t);
			Lim[t*B + b][0] = 1;
			Lim[t*B + b][1] = D;
		}
	return Lim;
}
vetorfloat2 LimRestDemandaBarraUnica(CSistema * const sistema_end, int n_a, const int cen)
{
	int T = sistema_end->GetTt2();
	size_t B = sistema_end->barrasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		double D = 0;
		for (size_t b = 0; b < B; b++)
		{
			for (size_t d = 0; d < sistema_end->barrasPtr[b]->demandasPtr.size(); d++)
				D = D + sistema_end->barrasPtr[b]->demandasPtr[d]->GetD(cen, t);
		}
		Lim[t][0] = 1;
		Lim[t][1] = D;
	}
	return Lim;
}
vetorfloat2 LimFluxo0(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int L = sistema_end->linhasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, L * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
		for (int l = 0; l < L; l++)
		{
			Lim[t*L + l][0] = 0;
			Lim[t*L + l][1] = sistema_end->linhasPtr[l]->GetCapacidade();
		}
	return Lim;
}
vetorfloat2 LimFluxo2(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int L = sistema_end->linhasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, L * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
		for (int l = 0; l < L; l++)
		{
			Lim[t*L + l][0] = 2;
			Lim[t*L + l][1] = - sistema_end->linhasPtr[l]->GetCapacidade();
		}
	return Lim;
}
vetorfloat2 LimPhgMin(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (size_t r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
				//Lim(jj + j,1) = double (sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetPmax() * sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades());
				Lim[jj + j][0] = 2;
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
	return Lim;
}
vetorfloat2 LimPhgMax(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 LimBalHid(CSistema * const sistema_end, int n_a, const int cen)
{
	int T = sistema_end->GetTt2();
	double deltaT;
	int R = sistema_end->hidreletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		if (t < sistema_end->GetTt1())
			deltaT = sistema_end->GetDeltaT1();
		else
			deltaT = sistema_end->GetDeltaT2();
		for (int r = 0; r < R; r++)
		{
			if (t == 0)
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (sistema_end->hidreletricasPtr[r]->GetV0() + 0.0036*deltaT*sistema_end->hidreletricasPtr[r]->GetAfluencia(cen,t));
			}
			else
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (0.0036*deltaT*sistema_end->hidreletricasPtr[r]->GetAfluencia(cen,t));
			}
		}
	}
	return Lim;
}
vetorfloat2 LimQMin(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (size_t r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
				//Lim(jj + j,1) = double (sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetQmax() * sistema_end->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades());
				Lim[jj + j][0] = 2;
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
	return Lim;
}
vetorfloat2 LimQMax(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t R = sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 LimTup(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	size_t I = sistema_end->termeletricasPtr.size();
	int Tup = 0;
	for (size_t i = 0; i < I; i++)
		Tup = Tup + sistema_end->termeletricasPtr[i]->GetTUp();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, Tup * T, 2);
	IniciaMatriz(&Lim, 0);		// usar ones e n precisar atribuir nos loops abaixo !!!
	int l = 0;
	for (size_t i = 0; i < I; i++)
	{
		for (int t = 0; t < T; t++)
		{
			if (sistema_end->termeletricasPtr[i]->GetTUp() > 0)
			{
				for (int tt = 0; tt < sistema_end->termeletricasPtr[i]->GetTUp(); tt++)
				{
					if (t - tt <= 0)		// if referente ao elemento tt
					{
						if (t - tt > - sistema_end->termeletricasPtr[i]->GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							Lim[l + tt][1] = double (sistema_end->termeletricasPtr[i]->GetU0());
							//Lim(l + tt,1) = 1;
							Lim[l + tt][0] = 2;
						}
						else	// estado da usina antes do perido -x0
						{	
							Lim[l + tt][1] = double ((1 - sistema_end->termeletricasPtr[i]->GetU0()));
							//Lim(l + tt,1) = 1;
							Lim[l + tt][0] = 2;
						}
					}
					else
					{
						Lim[l + tt][0] = 2;
					}
					if (t - tt - 1 <= 0)	// if referente ao elemento tt - 1
					{
						if (t - tt - 1 > - sistema_end->termeletricasPtr[i]->GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							Lim[l + tt][1] = Lim[l + tt][1] - sistema_end->termeletricasPtr[i]->GetU0();
							//Lim(l + tt,1) = 1;
							Lim[l + tt][0] = 2;
						}
						else	// estado da usina antes do perido -x0
						{	
							Lim[l + tt][1] = Lim[l + tt][1] - (1 - sistema_end->termeletricasPtr[i]->GetU0());
							//Lim(l + tt,1) = 1;
							Lim[l + tt][0] = 2;
						}
					}
					else
					{
						//Lim(l + tt,1) = 1;
						Lim[l + tt][0] = 2;
					}
				}
				l = l + sistema_end->termeletricasPtr[i]->GetTUp();
			}
			else
			{
					// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
			}
		}
	}
	return Lim;
}
vetorfloat2 LimTdown(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	int Tdown = 0;
	for (int i = 0; i < I; i++)
		Tdown = Tdown + sistema_end->termeletricasPtr[i]->GetTDown();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, Tdown * T, 2);
	IniciaMatriz(&Lim, 1);		// usar zeros e n precisar atribuir nos loops abaixo !!!
	int l = 0;
	for (int i = 0; i < I; i++)
	{
		for (int t = 0; t < T; t++)
		{
			if (sistema_end->termeletricasPtr[i]->GetTDown() > 0)
			{
				for (int tt = 0; tt < sistema_end->termeletricasPtr[i]->GetTDown(); tt++)
				{
					if (t - tt <= 0)		// if referente ao elemento tt
					{
						if (t - tt > - sistema_end->termeletricasPtr[i]->GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							Lim[l + tt][1] = Lim[l + tt][1] + sistema_end->termeletricasPtr[i]->GetU0();
							//Lim(l + tt,0) = - 1;
							Lim[l + tt][0] = 0;
						}
						else	// estado da usina antes do perido -x0
						{	
							Lim[l + tt][1] = Lim[l + tt][1] + (1 - sistema_end->termeletricasPtr[i]->GetU0());
							//Lim(l + tt,0) = - 1;
							Lim[l + tt][0] = 0;
						}
					}
					else
					{
						//Lim(l + tt,0) = - 1;
						Lim[l + tt][0] = 0;
					}
					if (t - tt - 1 <= 0)	// if referente ao elemento tt - 1
					{
						if (t - tt - 1 > - sistema_end->termeletricasPtr[i]->GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							Lim[l + tt][1] = Lim[l + tt][1] - sistema_end->termeletricasPtr[i]->GetU0();
							//Lim(l + tt,0) = - 1;
							Lim[l + tt][0] = 0;
						}
						else	// estado da usina antes do perido -x0
						{	
							Lim[l + tt][1] = Lim[l + tt][1] - (1 - sistema_end->termeletricasPtr[i]->GetU0());
							//Lim(l + tt,0) = - 1;
							Lim[l + tt][0] = 0;
						}
					}
					else
					{
						//Lim(l + tt,0) = - 1;
						Lim[l + tt][0] = 0;
					}
				}
				l = l + sistema_end->termeletricasPtr[i]->GetTDown();
			}
			else
			{
					// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
			}
		}
	}
	return Lim;
}
vetorfloat2 LimRampaUp(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T, 2);
	IniciaMatriz(&Lim, 1);
	int l = 0;
	for (int i = 0; i < I; i++)
	{
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
			{
				//Lim(l + t,0) = - double (sistema_end->termeletricasPtr[i]->GetPmin() + sistema_end->termeletricasPtr[i]->GetRampaUp() + sistema_end->termeletricasPtr[i]->GetRampaDown());
				Lim[l + t][0] = 0;
				Lim[l + t][1] = sistema_end->termeletricasPtr[i]->GetPmin() + sistema_end->termeletricasPtr[i]->GetPt0() + sistema_end->termeletricasPtr[i]->GetU0() * (sistema_end->termeletricasPtr[i]->GetRampaUp() - sistema_end->termeletricasPtr[i]->GetPmin());
			}
			else
			{
				//Lim(l + t,0) = - double (sistema_end->termeletricasPtr[i]->GetPmin() + sistema_end->termeletricasPtr[i]->GetRampaUp() + sistema_end->termeletricasPtr[i]->GetRampaDown());
				Lim[l + t][0] = 0;
				Lim[l + t][1] = sistema_end->termeletricasPtr[i]->GetPmin();
			}
		}
		l = l + T;
	}
	return Lim;
}
vetorfloat2 LimRampaDown(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T, 2);
	IniciaMatriz(&Lim, 1);
	int l = 0;
	for (int i = 0; i < I; i++)
	{	
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
			{
				//Lim(l + t,0) = - double (sistema_end->termeletricasPtr[i]->GetPmin() + sistema_end->termeletricasPtr[i]->GetRampaUp() + sistema_end->termeletricasPtr[i]->GetRampaDown());
				Lim[l + t][0] = 0;
				Lim[l + t][1] = sistema_end->termeletricasPtr[i]->GetPmin() - sistema_end->termeletricasPtr[i]->GetPt0();
			}
			else
			{
				//Lim(l + t,0) = - double (sistema_end->termeletricasPtr[i]->GetPmin() + sistema_end->termeletricasPtr[i]->GetRampaUp() + sistema_end->termeletricasPtr[i]->GetRampaDown());
				Lim[l + t][0] = 0;
				Lim[l + t][1] = sistema_end->termeletricasPtr[i]->GetPmin();
			}
		}
		l = l + T;
	}
	return Lim;
}
vetorfloat2 LimPtMin(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int i = 0; i < I; i++)
		for (int t = 0; t < T; t++)
		{
			//Lim(i*T + t,1) = double (sistema_end->termeletricasPtr[i]->GetPmax());
			Lim[i*T + t][0] = 2;
		}
	return Lim;
}
vetorfloat2 LimPtMax(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 LimRestCP(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int i = 0; i < I; i++)	
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
			{
				Lim[i*T + t][1] = -sistema_end->termeletricasPtr[i]->GetCoefCustoPartida()*sistema_end->termeletricasPtr[i]->GetU0();
				Lim[i*T + t][0] = 2;
			}
			else
			{
				Lim[i*T + t][1] = 0;
				Lim[i*T + t][0] = 2;
			}
		}
	return Lim;
}
vetorfloat2 LimVmeta(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int R = sistema_end->hidreletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R, 2);
	IniciaMatriz(&Lim, 0);
	for (int r = 0; r < R; r++)
	{
		Lim[r][1] = sistema_end->hidreletricasPtr[r]->GetVMeta();
		//Lim(r,1) = sistema_end->hidreletricasPtr[r]->GetVmax();
		Lim[r][0] = 2;
	}
	return Lim;
}
vetorfloat2 LimBalPotenciaL(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int R = sistema_end->hidreletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 LimBalVazaoL(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int R = sistema_end->hidreletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 LimFuncProdL(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int R = sistema_end->hidreletricasPtr.size();
	int JJ = 0;
	for (int i = 0; i < R; i++)
		JJ = JJ + sistema_end->hidreletricasPtr[i]->GetNGrupos();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, T * JJ, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_end->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				Lim[jj + j][0] = 0;
				Lim[jj + j][1] = sistema_end->hidreletricasPtr[r]->grupoPtr[j]->CoefPhg();
			}
			jj = jj + sistema_end->hidreletricasPtr[r]->GetNGrupos();
		}
	return Lim;
}
vetorfloat2 LimPhMax(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int R = sistema_end->hidreletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 LimReserva(CSistema * const sistema_end, int n_a, const int cen)
{
	int T = sistema_end->GetTt2();
	size_t B = sistema_end->barrasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1 * T, 2);
	IniciaMatriz(&Lim, 0);
	double D;
	for (int t = 0; t < T; t++)
	{
		D = 0;
		for (size_t b = 0; b < B; b++)
		{
			for (size_t d = 0; d < sistema_end->barrasPtr[b]->demandasPtr.size(); d++)
				D = D + sistema_end->barrasPtr[b]->demandasPtr[d]->GetD(cen, t);
		}
		Lim[t][0] = 2;
		Lim[t][1] = D*sistema_end->GetFatorReserva()/100;
	}
	return Lim;
}
vetorfloat2 LimCortesF(CSistema * const sistema_end, int n_a)
{
	int T = sistema_end->GetTt2();
	int I = sistema_end->termeletricasPtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T * sistema_end->GetFlagAproxCustoT(), 2);
	IniciaMatriz(&Lim, 0);
	int c = 0;
	for (int i = 0; i < I; i++)
	{
		for (int t = 0; t < T; t++)
		{
			for (int n_cort = 0; n_cort < sistema_end->GetFlagAproxCustoT(); n_cort++)
			{
				Lim[c][0] = 2;
				c++;
			}
		}
	}
	return Lim;
}
// ------------------------------------------------

// Matriz dos coeficientes e limites das restrições lineares
// ------------------------------------------------
// Caso seja necessário fazer modelagens diferentes para cada estágio é só criar uma condição que depende do T dentro da função
void MatrizRestricoesLineares(CSistema * const sistema_end, int n_a, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes)	
{
	CMatrizEsparsa M(0);
	int count = 0;
	if (sistema_end->GetFlagBarraUnica() == false)
	{
		M = MatrizRestDemanda(sistema_end, n_a);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo(sistema_end, n_a);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo(sistema_end, n_a);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	else
	{
		M = MatrizRestDemandaBarraUnica(sistema_end, n_a);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	M = MatrizLimPhgMin(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPhgMax(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalHid(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimqMin(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimqMax(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizTup(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizTdown(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizRampaUp(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizRampaDown(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPtMin(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPtMax(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizRestCP(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizVmeta(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalPotencia(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalVazao(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizFuncProd(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizPhMax(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizReserva(sistema_end, n_a);
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	if (sistema_end->GetFlagAproxCustoT() > 1)
	{
		M = MatrizCortesF(sistema_end, n_a);
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
}
void MatrizLimitesLineares(CSistema * const sistema_end, int n_a, vetorint * LimTipo, vetorfloat * LimValor, const int cen)
{
	vetorfloat2 L;
	if (sistema_end->GetFlagBarraUnica() == false)
	{
		L = LimRestDemanda(sistema_end, n_a, cen);
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo0(sistema_end, n_a);
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo2(sistema_end, n_a);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	else
	{
		L = LimRestDemandaBarraUnica(sistema_end, n_a, cen);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimPhgMin(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhgMax(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalHid(sistema_end, n_a, cen);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMin(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMax(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimTup(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimTdown(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimRampaUp(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimRampaDown(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPtMin(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPtMax(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimRestCP(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimVmeta(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalPotenciaL(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalVazaoL(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimFuncProdL(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhMax(sistema_end, n_a);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimReserva(sistema_end, n_a, cen);
	AlocarLimites(&L, LimTipo, LimValor);
	if (sistema_end->GetFlagAproxCustoT() > 1)
	{
		L = LimCortesF(sistema_end, n_a);
		AlocarLimites(&L, LimTipo, LimValor);
	}
}
// ------------------------------------------------
};

using namespace met_DecCen;

CSubProblemaCen::CSubProblemaCen(CSistema * const sistema_end, int cenario_a, const GRBEnv &ambiente_gurobi): modeloGRB(ambiente_gurobi)
{
	sistema_a = sistema_end;
	cenario = cenario_a;

	int flag1 = int (1 - sistema_a->GetFlagBarraUnica());
	int flag4;
	if (sistema_end->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasPtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasPtr[i]->GetNGrupos();
	n = sistema_a->GetTt2() * ((3+flag4)*sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 5*(sistema_a->hidreletricasPtr.size()) + 3*JJ + flag1*sistema_a->barrasPtr.size() + (1 - flag1)) + sistema_a->hidreletricasPtr.size();
	x.resize(n);
	n_restricoes = 0;

	// Criar variáveis
	CriarVariaveis();
	modeloGRB.update();	// Atualiza o modelo Gurobi.
	vars = modeloGRB.getVars();
	
	// Adicionar restrições
	CriarRestricoes();
	modeloGRB.update();	// Atualiza o modelo Gurobi.

	// criar vetores e matrizes do modelo somente no constructor, o atributo da clase é o modelo do Gurobi!
}
CSubProblemaCen::CSubProblemaCen(CSistema * const sistema_end, int cenario_a, CSubProblemaCen &cenario0): modeloGRB(cenario0.modeloGRB)
{
	sistema_a = sistema_end;
	cenario = cenario_a;
	vars = modeloGRB.getVars();

	// receber atributos da outra classe copia, por exemplo n_restricoes
	n_restricoes = cenario0.GetN_Restricoes();
	n = cenario0.GetN();
	x.resize(n);

	AtualizarRestricoes();
	modeloGRB.update();
}

CSubProblemaCen::~CSubProblemaCen(void)
{
}

void CSubProblemaCen::CriarVariaveis()
{
	try 
	{
		size_t I = sistema_a->termeletricasPtr.size();
		size_t R = sistema_a->hidreletricasPtr.size();
		size_t B = sistema_a->barrasPtr.size() - 1;
		
		for (int t = 0; t < sistema_a->GetTt2(); t++)
		{
			//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
			for (size_t i = 0; i < I; i++)	//pt
				modeloGRB.addVar(0, double (sistema_a->termeletricasPtr[i]->GetPmax()), 0.0, GRB_CONTINUOUS, "");
			if (sistema_a->GetFlagVarBin() == true)	//u (biário ou continuo)
			{
				for (size_t i = 0; i < I; i++)
					modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)
					modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			//for (size_t i = 0; i < I; i++)	//x
			//	modeloGRB.addVar(0, int (sistema_a->termeletricasPtr[i]->GetX0() + sistema_a->GetTt1()), 0.0, GRB_INTEGER, "");
			for (size_t i = 0; i < I; i++)	//cp
				modeloGRB.addVar(0, sistema_a->termeletricasPtr[i]->GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
			if (sistema_a->GetFlagAproxCustoT() > 1)	//F
			{
				for (size_t i = 0; i < I; i++)
					modeloGRB.addVar(0, double (sistema_a->termeletricasPtr[i]->CustoOperacao(sistema_a->termeletricasPtr[i]->GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagBarraUnica() == false)
				for (size_t b = 0; b < B; b++)	//teta
					modeloGRB.addVar(-6.2832, 6.2832, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetPmax()*sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades();
				modeloGRB.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)	//v
				modeloGRB.addVar(double(sistema_a->hidreletricasPtr[r]->GetVmin()), double(sistema_a->hidreletricasPtr[r]->GetVmax()), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetQmax() * sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades());
				modeloGRB.addVar(0, qhmax + sistema_a->hidreletricasPtr[r]->GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)	//s
			{
				modeloGRB.addVar(0, sistema_a->hidreletricasPtr[r]->GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)	//phmax
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetPmax()*sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades();
				modeloGRB.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < R; r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)	
				{
					modeloGRB.addVar(0, double (sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetPmax() * sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)	
				{
					modeloGRB.addVar(0, double (sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetQmax() * sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if (sistema_a->GetFlagVarBin() == true)	//z (biário ou continuo)
			{
				for (size_t r = 0; r < R; r++)	
					for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)	
						modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)	
						modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagBarraUnica() == false)
			{
				for (size_t b = 0; b < B + 1; b++)	//def
				{
					double cap_d = 0;
					for (size_t db = 0; db < sistema_a->barrasPtr[b]->demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasPtr[b]->demandasPtr[db]->GetD(cenario, t);
					modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				double cap_d = 0;					//def barra única
				for (size_t b = 0; b < B + 1; b++)	
					for (size_t db = 0; db < sistema_a->barrasPtr[b]->demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasPtr[b]->demandasPtr[db]->GetD(cenario, t);
				modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
			}
		}
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)	//vfol
			modeloGRB.addVar(0, double(sistema_a->hidreletricasPtr[r]->GetVmax() - sistema_a->hidreletricasPtr[r]->GetVmin()), 0.0, GRB_CONTINUOUS, "");
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
void CSubProblemaCen::CriarRestricoes()
{
	try
	{
		// Criar modelo do Gurobi
		vetorint indexL, indexC, nnZ, LimTipo;
		vetorfloat indexV, LimValor;
		indexL.resize(0);indexC.resize(0);indexV.resize(0);nnZ.resize(0);LimTipo.resize(0);LimValor.resize(0);
		nnZ.push_back(0);
		MatrizRestricoesLineares(sistema_a, n, &indexL, &indexC, &indexV, &nnZ ,&n_restricoes);
		MatrizLimitesLineares(sistema_a, n, &LimTipo, &LimValor, 0);
		
		GRBLinExpr restricao;
		//double * coeffs;
		//coeffs = new double[n];
		double coeficiente;
		GRBVar variavel;
		//
		for (int i = 0; i < n_restricoes; i++)
		{
			for (int j = nnZ[i]; j < nnZ[i + 1]; j++)
			{
				coeficiente = indexV[j];
				variavel = vars[indexC[j]];
				restricao.addTerms( &coeficiente, &variavel, 1);
			}
			//for (int j = 0; j < n; j++)		demora muito mais tempo
			//	coeffs[j] = 0;
			//for (int j = nnZ[i]; j <  nnZ[i + 1]; j++)
			//	coeffs[indexC[j]] = double (indexV[j]);
			//	//coeffs[indexC[j]] = double (indexV[j]);
			//restricao.addTerms(coeffs, vars, n);
			switch (LimTipo[i])
			{
			case 0:
				modeloGRB.addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modeloGRB.addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modeloGRB.addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
				break;
			default:
				cout << "Tipo inválido de restrição adicionada" << endl;
			}
			restricao.clear();
		}
		//delete coeffs;
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
	//// Escrever em arquivo
	//ofstream inFile( "modelo.txt", ios::out );   
	//if ( inFile.is_open() )
	//{
	//	for (size_t i = 0; i < indexL.size(); i++)
	//	{
	//		inFile << left << setw(5) << i << "	" << setw(5) << indexL[i] << "	" << setw(5) << indexC[i] << "	" << setw(15) << indexV[i];
	//		inFile << endl;
	//	}
	//inFile.close();
	//}
	//else
	//	cout << "Unable to open file";
	//// Fim (Escreve em arquivo)
}
void CSubProblemaCen::CriarFuncaoObjetivoRL(const int cenario, const vetorfloat * const Lambda)
{
	try 
	{
		GRBQuadExpr fo;
		//GRBLinExpr fo;

		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
		double * coeffs;
		coeffs = new double[n];
		int flag1 = int (1 - sistema_a->GetFlagBarraUnica());
		int flag4;
		if (sistema_a->GetFlagAproxCustoT() > 1)
			flag4 = 1;
		else
			flag4 = 0;
		int JJ = 0;
		for (size_t i = 0; i < sistema_a->hidreletricasPtr.size(); i++)
			JJ += sistema_a->hidreletricasPtr[i]->GetNGrupos();
		for (int i = 0; i < n; i++)
			coeffs[i] = 0;
		int nt = (n - sistema_a->hidreletricasPtr.size())/sistema_a->GetTt2();
		double deltaT;
		// Termos lineares
		for (int t = 0; t < sistema_a->GetTt2(); t++)
		{
			if (t < sistema_a->GetTt1())
				deltaT = sistema_a->GetDeltaT1();
			else
				deltaT = sistema_a->GetDeltaT2();
			for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
			{
				if (sistema_a->GetFlagAproxCustoT() > 1)
					coeffs[i + 3*sistema_a->termeletricasPtr.size() + t*nt] = deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//F
				else
				{
					coeffs[i + t*nt] = sistema_a->termeletricasPtr[i]->GetCoefCustoOper(1) * deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//pt
					coeffs[i + sistema_a->termeletricasPtr.size() + t*nt] = sistema_a->termeletricasPtr[i]->GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//u
				}
				coeffs[i + 2*sistema_a->termeletricasPtr.size() + t*nt] = 1 * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//cp
			}
		}
		int delta = nt - flag1*sistema_a->barrasPtr.size() - (1 - flag1);
		for (int t = 0; t < sistema_a->GetTt2(); t++)
		{
			if (t < sistema_a->GetTt1())
				deltaT = sistema_a->GetDeltaT1();
			else
				deltaT = sistema_a->GetDeltaT2();
			if (sistema_a->GetFlagBarraUnica() == false)
			{
				for (size_t b = 0; b < sistema_a->barrasPtr.size(); b++)
					coeffs[delta + b + t*nt] = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//def
			}
			else
			{
				coeffs[delta + t*nt] = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//def
			}
		}
		delta = sistema_a->GetTt2() * nt;
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			coeffs[delta + r] = (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()* sistema_a->GetDeltaT2() * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//vfol
		}
		fo.addTerms(coeffs, vars, n);
		// Termos quadráticos
		if (sistema_a->GetFlagAproxCustoT() > 1)
		{
		}
		else
		{
			for (int i = 0; i < n; i++)
				coeffs[i] = 0;
			for (int t = 0; t < sistema_a->GetTt2(); t++)
			{
				if (t < sistema_a->GetTt1())
					deltaT = sistema_a->GetDeltaT1();
				else
					deltaT = sistema_a->GetDeltaT2();
				for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
					coeffs[i + t*nt] = sistema_a->termeletricasPtr[i]->GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//pt^2
			}
			fo.addTerms(coeffs, vars, vars, n);
		}
		// Termos da RL
		// Constantes

		// Lineares
		int delta_x = 0;
		int delta_subgrad = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)		
		{
			for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
				fo.addTerm(Lambda->at(delta_subgrad + i), vars[delta_x + i]);
			delta_subgrad += sistema_a->termeletricasPtr.size();
			delta_x += (3+flag4) * sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 2*sistema_a->hidreletricasPtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
				fo.addTerm(Lambda->at(delta_subgrad + r), vars[delta_x + r]);
			delta_subgrad += sistema_a->hidreletricasPtr.size();
			delta_x += 3*sistema_a->hidreletricasPtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
			{
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
					fo.addTerm(Lambda->at(delta_subgrad + j), vars[delta_x + j]);
				delta_subgrad += sistema_a->hidreletricasPtr[r]->GetNGrupos();
				delta_x += sistema_a->hidreletricasPtr[r]->GetNGrupos();
			}
			delta_x += 2*JJ + flag1*(sistema_a->barrasPtr.size()) + (1 - flag1);
		}
		delete coeffs;
		modeloGRB.setObjective(fo, GRB_MINIMIZE);
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
void CSubProblemaCen::CriarFuncaoObjetivoRP(const int cenario, const vetorfloat * const Lambda, const vetorfloat * const mi, const vetorfloat * const x_med)
{
	try 
	{
		GRBQuadExpr fo;
		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
		double * coeffs;
		coeffs = new double[n];
		int flag1 = int (1 - sistema_a->GetFlagBarraUnica());
		for (int i = 0; i < n; i++)
			coeffs[i] = 0;
		int nt = (n - sistema_a->hidreletricasPtr.size())/sistema_a->GetTt2();
		double deltaT;
		// Termos lineares
		for (int t = 0; t < sistema_a->GetTt2(); t++)
		{
			if (t < sistema_a->GetTt1())
				deltaT = sistema_a->GetDeltaT1();
			else
				deltaT = sistema_a->GetDeltaT2();
			for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
			{
				if (sistema_a->GetFlagAproxCustoT() > 1)
					coeffs[i + 3*sistema_a->termeletricasPtr.size() + t*nt] = deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//F
				else
				{
					coeffs[i + t*nt] = sistema_a->termeletricasPtr[i]->GetCoefCustoOper(1) * deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//pt
					coeffs[i + sistema_a->termeletricasPtr.size() + t*nt] = sistema_a->termeletricasPtr[i]->GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//u
				}
				coeffs[i + 2*sistema_a->termeletricasPtr.size() + t*nt] = 1 * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//cp
			}
		}
		int delta = nt - flag1*sistema_a->barrasPtr.size() - (1 - flag1);
		for (int t = 0; t < sistema_a->GetTt2(); t++)
		{
			if (t < sistema_a->GetTt1())
				deltaT = sistema_a->GetDeltaT1();
			else
				deltaT = sistema_a->GetDeltaT2();
			if (sistema_a->GetFlagBarraUnica() == false)
			{
				for (size_t b = 0; b < sistema_a->barrasPtr.size(); b++)
					coeffs[delta + b + t*nt] = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//def
			}
			else
			{
				coeffs[delta + t*nt] = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//def
			}
		}
		delta = sistema_a->GetTt2() * nt;
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			coeffs[delta + r] = (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()* sistema_a->GetDeltaT2() * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//vfol
		}
		fo.addTerms(coeffs, vars, n);
		// Termos quadráticos
		if (sistema_a->GetFlagAproxCustoT() > 1)
		{
		}
		else
		{
			for (int i = 0; i < n; i++)
				coeffs[i] = 0;
			for (int t = 0; t < sistema_a->GetTt2(); t++)
			{
				if (t < sistema_a->GetTt1())
					deltaT = sistema_a->GetDeltaT1();
				else
					deltaT = sistema_a->GetDeltaT2();
				for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
					coeffs[i + t*nt] = sistema_a->termeletricasPtr[i]->GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cenario);		//pt^2
			}
			fo.addTerms(coeffs, vars, vars, n);
		}
		// Termos do PH
		// Constantes
		for (int i = 0; i < nt * sistema_a->GetTt1(); i++)
			fo.addConstant(- Lambda->at(i) * x_med->at(i) + pow(x_med->at(i),2) * mi->at(i) / 2 );
		// Lineares
		for (int i = 0; i < nt * sistema_a->GetTt1(); i++)
			fo.addTerm(Lambda->at(i) - x_med->at(i) * mi->at(i), vars[i]);
		// Quadráticos
		for (int i = 0; i < nt * sistema_a->GetTt1(); i++)
			fo.addTerm(mi->at(i) / 2, vars[i], vars[i]);
		delete coeffs;
		modeloGRB.setObjective(fo, GRB_MINIMIZE);
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
void CSubProblemaCen::AtualizarRestricoes()
{
	try
	{
		int flag1 = int (1 - sistema_a->GetFlagBarraUnica());
		// Atualizar restrições de balanço hídrico
		int JJ = 0;
		for (size_t i = 0; i < sistema_a->hidreletricasPtr.size(); i++)
			JJ = JJ + sistema_a->hidreletricasPtr[i]->GetNGrupos();
		int indexBalHid = (flag1*sistema_a->barrasPtr.size() + (1 - flag1) + flag1 * sistema_a->linhasPtr.size() * 2 + JJ * 2) * sistema_a->GetTt2();
		vetorfloat2 atualBalHid = LimBalHid(sistema_a, n, cenario);
		GRBConstr *restricoes;
		restricoes = modeloGRB.getConstrs();
		for (size_t i = indexBalHid; i < indexBalHid + atualBalHid.size(); i++)
			restricoes[i].set(GRB_DoubleAttr_RHS, double (atualBalHid.at(i - indexBalHid).at(1)));
		modeloGRB.update();
		// Atualizar restrições de atendimento a demanda
		int indexDemanda = 0;
		if (sistema_a->GetFlagBarraUnica() == false)
		{
			vetorfloat2 atualDemanda = LimRestDemanda(sistema_a, n, cenario);
			for (size_t i = indexDemanda; i < indexDemanda + atualDemanda.size(); i++)
				restricoes[i].set(GRB_DoubleAttr_RHS, double (atualDemanda.at(i - indexDemanda).at(1)));
		}
		else
		{
			vetorfloat2 atualDemanda = LimRestDemandaBarraUnica(sistema_a, n, cenario);
			for (size_t i = indexDemanda; i < indexDemanda + atualDemanda.size(); i++)
				restricoes[i].set(GRB_DoubleAttr_RHS, double (atualDemanda.at(i - indexDemanda).at(1)));
		}
		modeloGRB.update();
		// Atualizar restrições de reserva
		int indexReserva = n_restricoes - sistema_a->GetTt2();
		vetorfloat2 atualReserva = LimReserva(sistema_a, n, cenario);
		for (size_t i = indexReserva; i < indexReserva + atualReserva.size(); i++)
			restricoes[i].set(GRB_DoubleAttr_RHS, double (atualReserva.at(i - indexReserva).at(1)));
		modeloGRB.update();
		
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
int CSubProblemaCen::ResolverProblemaRL(CResultados * resultadoGurobi, const vetorfloat * const lambda, const int iter)
{
	try
	{
		// Criar função objetivo
		CriarFuncaoObjetivoRL(cenario, lambda);		// para problemas grandes essa atualização da função objetivo demora muito ?
		modeloGRB.update();	// Atualiza o modelo Gurobi.
		
		//// Atualizar valores de afluencias
		//AtualizarRestricoes(cenario);
		//modeloGRB.update();	// Atualiza o modelo Gurobi.

		//modeloGRB.write("prob.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modeloGRB.getEnv().set(GRB_IntParam_OutputFlag, 0);
		
		//modeloGRB.getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modeloGRB.getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		//// fazer valor da iteração como entrada para alterar o mipgap (10%, 5%, 2%, 1%, para it = 0,1,2,3,...)
		//switch (iter)
		//	{
		//	case 0:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.1);	// Define o gap de tolerancia
		//		break;
		//	case 1:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.05);	// Define o gap de tolerancia
		//		break;
		//	case 2:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.02);	// Define o gap de tolerancia
		//		break;
		//	default:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.01);	// Define o gap de tolerancia
		//	}
		
		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.1);	// Define o gap de tolerancia
		
		modeloGRB.getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		modeloGRB.getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		modeloGRB.getEnv().set(GRB_IntParam_Threads, 4);
		//modeloGRB.getEnv().set(GRB_DoubleParam_NodeLimit, 100000);  // Limita o numero de nós a serem explorados
		modeloGRB.getEnv().set(GRB_DoubleParam_TimeLimit, 600);		// Limita o tempo de resolução do problema
		
		// Desativa o presolver
		//modeloGRB.getEnv().set(GRB_IntParam_Presolve,0);
		
		// Reset the model to an unsolved state, discarding any previously computed solution information
		//modeloGRB.reset();

		// Limita o numero de Threads na resolucao do problema
		//modeloGRB.getEnv().set(GRB_IntParam_Threads, 4);
		
		// Otimizar
		modeloGRB.optimize();

		// Arrumar alguma forma de resolver o problema quando dá algum erro!!! (depois do catch ou no Problema_PH.cpp)
		//if (GRB_OPTIMAL != 2)
		//{
		//	modeloGRB.getEnv().set(GRB_IntParam_MIPFocus, 1);
		//	modeloGRB.optimize();
		//}

		// Salvar resultados
		int nStatus = modeloGRB.get(GRB_IntAttr_Status);
		if ((nStatus == 2) || (nStatus == 9))
		{
			fo = double(modeloGRB.get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n; i++)
			{	
				x[i] = double(vars[i].get(GRB_DoubleAttr_X));
			}
		}
		else
		{
			for (int i = 0; i < n; i++)
			{	
				x[i] = 0;
			}
			fo = 0;
		}
		resultadoGurobi->CarregarResultadosDecCen(fo,x,L, nStatus);
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->CarregarResultadosDecCen(fo,x,L, e.getErrorCode());
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}
int CSubProblemaCen::ResolverProblemaRP(CResultados * resultadoGurobi, const vetorfloat * const lambda, const vetorfloat * const mi, const vetorfloat * const x_med, const int iter)
{
	try
	{
		// Criar função objetivo
		CriarFuncaoObjetivoRP(cenario, lambda, mi, x_med);		// para problemas grandes essa atualização da função objetivo demora muito ?
		modeloGRB.update();	// Atualiza o modelo Gurobi.
		
		//// Atualizar valores de afluencias
		//AtualizarRestricoes(cenario);
		//modeloGRB.update();	// Atualiza o modelo Gurobi.

		//modeloGRB.write("prob.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modeloGRB.getEnv().set(GRB_IntParam_OutputFlag, 0);
		
		//modeloGRB.getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modeloGRB.getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		//// fazer valor da iteração como entrada para alterar o mipgap (10%, 5%, 2%, 1%, para it = 0,1,2,3,...)
		//switch (iter)
		//	{
		//	case 0:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.1);	// Define o gap de tolerancia
		//		break;
		//	case 1:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.05);	// Define o gap de tolerancia
		//		break;
		//	case 2:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.02);	// Define o gap de tolerancia
		//		break;
		//	default:
		//		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.01);	// Define o gap de tolerancia
		//	}
		
		modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.1);	// Define o gap de tolerancia
		
		modeloGRB.getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		modeloGRB.getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		modeloGRB.getEnv().set(GRB_IntParam_Threads, 4);
		//modeloGRB.getEnv().set(GRB_DoubleParam_NodeLimit, 100000);  // Limita o numero de nós a serem explorados
		modeloGRB.getEnv().set(GRB_DoubleParam_TimeLimit, 600);		// Limita o tempo de resolução do problema
		
		// Desativa o presolver
		//modeloGRB.getEnv().set(GRB_IntParam_Presolve,0);
		
		// Reset the model to an unsolved state, discarding any previously computed solution information
		modeloGRB.reset();

		// Limita o numero de Threads na resolucao do problema
		//modeloGRB.getEnv().set(GRB_IntParam_Threads, 4);
		
		// Otimizar
		modeloGRB.optimize();

		// Arrumar alguma forma de resolver o problema quando dá algum erro!!! (depois do catch ou no Problema_PH.cpp)
		//if (GRB_OPTIMAL != 2)
		//{
		//	modeloGRB.getEnv().set(GRB_IntParam_MIPFocus, 1);
		//	modeloGRB.optimize();
		//}

		// Salvar resultados
		int nStatus = modeloGRB.get(GRB_IntAttr_Status);
		if ((nStatus == 2) || (nStatus == 9))
		{
			fo = double(modeloGRB.get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n; i++)
			{	
				x[i] = double(vars[i].get(GRB_DoubleAttr_X));
			}
		}
		else
		{
			for (int i = 0; i < n; i++)
			{	
				x[i] = 0;
			}
			fo = 0;
		}
		resultadoGurobi->CarregarResultadosDecCen(fo,x,L, nStatus);
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->CarregarResultadosDecCen(fo,x,L, e.getErrorCode());
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}
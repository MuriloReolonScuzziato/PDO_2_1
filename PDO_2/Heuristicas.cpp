/*--------------------------------------------------------------------------*/
/*---------------------------- File Heuristicas.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para aplicar heuristicas na solução da RL (considerando a decom-
 * posiçao SpcDec_2
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Heuristicas.h"

// Funcoes para criar a matriz de restriçoes
CMatrizEsparsa Heuristicas::MatrizRestDemanda()	// Monta matriz da restrição de atendimento a demanda
{
	int n_a = n;
	CMatrizEsparsa Agh(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
				if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
					Agh.InserirElemento(b, r, 1);
	CMatrizEsparsa Agt(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
				if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
					Agt.InserirElemento(b, i, 1);
	CMatrizEsparsa B(int(sistema_a->barrasVtr.size()), int(sistema_a->barrasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t bb = 0; bb < sistema_a->barrasVtr.size(); bb++)
			for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
				if (b == bb)
				{
					if ((sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[b]) || (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[b]))
						if ( B.GetElemento(b, bb) != 0 )
							B.SubstituirElemento(b, bb, B.GetElemento(b, bb) + 100/(sistema_a->linhasVtr[l].GetReatancia()));
						else
							B.InserirElemento(b, bb, 100/(sistema_a->linhasVtr[l].GetReatancia()));
				}
				else
				{
					if (((sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[b]) && (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[bb])) ||
						((sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[bb]) && (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[b])))
						if ( B.GetElemento(b, bb) != 0 )
							B.SubstituirElemento(b, bb, B.GetElemento(b, bb) - 100/(sistema_a->linhasVtr[l].GetReatancia()));
						else
							B.InserirElemento(b, bb, - 100/(sistema_a->linhasVtr[l].GetReatancia()));
				}
	B.RemoverColuna(sistema_a->GetBarraRef() - 1);
	B.MultiplicarPorEscalar( -1);
	CMatrizEsparsa Alin(T * sistema_a->barrasVtr.size(), n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c, l + sistema_a->barrasVtr.size() - 1, c + sistema_a->termeletricasVtr.size() - 1, &Agt, 0, 0);
		Alin.InserirMatriz(l, c + ((3+flag4)*sistema_a->termeletricasVtr.size()),l + sistema_a->barrasVtr.size() - 1, c + ((3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &B, 0, 0);
		Alin.InserirMatriz(l, c + (sistema_a->barrasVtr.size() - 1 + (3+flag4)*sistema_a->termeletricasVtr.size()), l + sistema_a->barrasVtr.size() - 1, c + (sistema_a->barrasVtr.size() - 1 + (3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size()) - 1, &Agh, 0, 0);
		l = l + sistema_a->barrasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	c = (3+flag4)*(sistema_a->termeletricasVtr.size()) + sistema_a->barrasVtr.size() - 1 + 5*(sistema_a->hidreletricasVtr.size()) + 3*JJ;
	l = 0;
	CMatrizEsparsa Adef(sistema_a->barrasVtr.size());
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c,l + sistema_a->barrasVtr.size() - 1,c + sistema_a->barrasVtr.size() - 1, &Adef, 0, 0);
		l = l + sistema_a->barrasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizRestDemandaBarraUnica()	// Monta matriz da restrição de atendimento a demanda
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = 0;
	for ( int t = 0; t < T; t++)
	{
		a.RemoverTodosElementos();
		for (size_t i = 0; i < I; i++)		// pt
			a.InserirElemento(0, i + c, 1);
		for (size_t r = 0; r < R; r++)		// ph
			a.InserirElemento(0, r + c + (3+flag4)*sistema_a->termeletricasVtr.size(), 1);
		a.InserirElemento(0, (3+flag4)*sistema_a->termeletricasVtr.size() + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + c, 1);		// def
		Alin.JuntarColuna(&a);
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimFluxo()
{
	int n_a = n;
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
	CMatrizEsparsa Alin(T * sistema_a->linhasVtr.size(), n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(ll,c + ((3+flag4)*sistema_a->termeletricasVtr.size()),ll + sistema_a->linhasVtr.size() - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
		ll = ll + sistema_a->linhasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimPhgMin()
{
	int n_a = n;
	int JJ = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ,JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmin());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ, n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()),l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimPhgMax()
{
	int n_a = n;
	int JJ = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ,n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()),l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizBalHid()
{
	int n_a = n;
	double deltaT;
	int tempo_viagem;
	int R = sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Addd(R * T, R);
	CMatrizEsparsa Add(R, R);
	for (int t = 0; t < T; t++)
	{
		deltaT = sistema_a->GetDeltaT1();
		Add.RemoverTodosElementos();
		for (int r = 0; r < R; r++)
		{
			if (sistema_a->hidreletricasVtr[r].GetTempoViagem() < deltaT)
				tempo_viagem = 0;
			else
				tempo_viagem = int (sistema_a->hidreletricasVtr[r].GetTempoViagem() / deltaT);
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[r].GetUsinaJusante() != 0))
			{
				Add.InserirElemento(sistema_a->hidreletricasVtr[r].GetUsinaJusante() - 1, r, double (-0.0036 * deltaT));
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
			if (sistema_a->hidreletricasVtr[r].GetTempoViagem() < deltaT)
				tempo_viagem = 0;
			else
				tempo_viagem = int (sistema_a->hidreletricasVtr[r].GetTempoViagem() / deltaT);
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[r].GetUsinaJusante() != 0))
			{
				Add.InserirElemento(sistema_a->hidreletricasVtr[r].GetUsinaJusante() - 1, r, double (-0.0036 * deltaT));
			}
		}
		Addd2.InserirMatriz(t * R,0,(t + 1) * R - 1,R - 1, &Add, 0, 0);
	}
	CMatrizEsparsa Ad(R * T);
	CMatrizEsparsa MSoma1(R, R);
	CMatrizEsparsa MSoma2(R, R);
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
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				{
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
		}
	}
	CMatrizEsparsa Av(R);
	CMatrizEsparsa Alin(R * T,n_a);
	CMatrizEsparsa AvNeg(R);
	AvNeg.MultiplicarPorEscalar( -1);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &Av, 0, 0);
			for (int tt = 0; tt <= t; tt++)
				Alin.InserirMatriz(l,tt*(n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R),l + R - 1,tt*(n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 3*R) - 1, &Ad, t*R, tt*R);
			if (t > 0)
				Alin.InserirMatriz(l,c - (n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c - (n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
		}
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				{
					Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &Av, 0, 0);
					cc = 0;
					for (int tt = 0; tt <= t; tt++)
					{
						Alin.InserirMatriz(l,cc + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R),l + R - 1,cc + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 3*R) - 1, &Ad, t*R, tt*R);
						if (((tt + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((tt + 1 - sistema_a->GetTt1()) > 0))
							cc = cc + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
						else
							cc = cc + (n_a / T);
					}
					if (t > 0)
						if (t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) == sistema_a->GetTt1())
							Alin.InserirMatriz(l,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
						else
							Alin.InserirMatriz(l,c - (n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c - (n_a / T) + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
				}
		l = l + R;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimqMin()
{
	int n_a = n;
	int JJ = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmin());
			cont++;
		}
	CMatrizEsparsa Alin( T * JJ,n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + JJ,l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimqMax()
{
	int n_a = n;
	int JJ = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ, JJ);
	int cont = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ,n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + JJ,l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 2*JJ),l + JJ - 1,c + ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size()) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizTup()
{
	// T * sum(Tup)
	int n_a = n;
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int Tup1, Tup2, Tup, cen, t_a;
	int c = I;
	int ca = 0;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1()) - flag2*cen*sistema_a->hidreletricasVtr.size();
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}
		for (size_t i = 0; i < I; i++)
		{
			if (sistema_a->termeletricasVtr[i].GetTUp() < sistema_a->GetDeltaT1())
				Tup1 = 0;
			else
				Tup1 = int (sistema_a->termeletricasVtr[i].GetTUp() / sistema_a->GetDeltaT1());
			if (sistema_a->termeletricasVtr[i].GetTUp() < sistema_a->GetDeltaT2())
				Tup2 = 0;
			else
				Tup2 = int (sistema_a->termeletricasVtr[i].GetTUp() / sistema_a->GetDeltaT2());

			if (t < sistema_a->GetTt1())
				Tup = Tup1;
			else
				Tup = Tup2;

			if (Tup > 0)
			{
				for (int tt = 0; tt < Tup; tt++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i, 1);
					if (t_a - tt - 1 >= 0)		// cada if adiciona um elemento dos dois elementos do lado direito da equação
						if (t_a - tt - 1 >= sistema_a->GetTt1())
							a.InserirElemento(0, c + i - tt*(n_a / T) - 1*(n_a / T), -1);
						else
							a.InserirElemento(0, ca + i - tt*(n_a / T) - 1*(n_a / T), -1);
					if (t_a - tt - 2 >= 0)
						if (t_a - tt - 2 >= sistema_a->GetTt1())
							a.InserirElemento(0, c + i - tt*(n_a / T) - 2*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + i - tt*(n_a / T) - 2*(n_a / T), 1);
					
					Alin.JuntarColuna(&a);
				}
			}
			else
			{
			}
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizTdown()
{
	// T * sum(Tdown)
	int n_a = n;
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int Tdown1, Tdown2, Tdown, cen, t_a;
	int c = I;
	int ca = 0;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1()) - flag2*cen*sistema_a->hidreletricasVtr.size();
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}
		for (size_t i = 0; i < I; i++)
		{
			if (sistema_a->termeletricasVtr[i].GetTDown() < sistema_a->GetDeltaT1())
				Tdown1 = 0;
			else
				Tdown1 = int (sistema_a->termeletricasVtr[i].GetTDown() / sistema_a->GetDeltaT1());
			if (sistema_a->termeletricasVtr[i].GetTDown() < sistema_a->GetDeltaT2())
				Tdown2 = 0;
			else
				Tdown2 = int (sistema_a->termeletricasVtr[i].GetTDown() / sistema_a->GetDeltaT2());

			if (t < sistema_a->GetTt1())
				Tdown = Tdown1;
			else
				Tdown = Tdown2;

			if (Tdown > 0)
			{
				for (int tt = 0; tt < Tdown; tt++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i, 1);
					if (t_a - tt - 1 >= 0)		// cada if adiciona um elemento dos dois elementos do lado direito da equação
						if (t_a - tt - 1 >= sistema_a->GetTt1())
							a.InserirElemento(0, c + i - tt*(n_a / T) - 1*(n_a / T), -1);
						else
							a.InserirElemento(0, ca + i - tt*(n_a / T) - 1*(n_a / T), -1);
					if (t_a - tt - 2 >= 0)
						if (t_a - tt - 2 >= sistema_a->GetTt1())
							a.InserirElemento(0, c + i - tt*(n_a / T) - 2*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + i - tt*(n_a / T) - 2*(n_a / T), 1);
					
					Alin.JuntarColuna(&a);
				}
			}
			else
			{
			}
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizRampaUp()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = 0;
	int ca, t_a, cen;
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1()) - flag2*cen*sistema_a->hidreletricasVtr.size();
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}

		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + i, 1);
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
				{
					a.InserirElemento(0, c + i - (n_a / T), -1);
					a.InserirElemento(0, c + i + I - (n_a / T), sistema_a->termeletricasVtr[i].GetPmin() - sistema_a->termeletricasVtr[i].GetRampaUp());
				}
				else
				{
					a.InserirElemento(0, ca + i - (n_a / T), -1);
					a.InserirElemento(0, ca + i + I - (n_a / T), sistema_a->termeletricasVtr[i].GetPmin() - sistema_a->termeletricasVtr[i].GetRampaUp());
				}
			}
			Alin.JuntarColuna(&a);
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizRampaDown()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = 0;
	int ca, t_a, cen;
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1()) - flag2*cen*sistema_a->hidreletricasVtr.size();
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}

		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + i, -1);
			a.InserirElemento(0, c + i + I, sistema_a->termeletricasVtr[i].GetPmin() - sistema_a->termeletricasVtr[i].GetRampaDown());
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
					a.InserirElemento(0, c + i - (n_a / T), 1);
				else
					a.InserirElemento(0, ca + i - (n_a / T), 1);
			}
			Alin.JuntarColuna(&a);
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimPtMin()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + i, 1);
			a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetPmin());
			Alin.JuntarColuna(&a);
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizLimPtMax()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + i, 1);
			a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetPmax());
			Alin.JuntarColuna(&a);
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizRestCP()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = I;
	int ca, t_a, cen;
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1()) - flag2*cen*sistema_a->hidreletricasVtr.size();
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}

		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + i + I, 1);
			a.InserirElemento(0, c + i, - sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
					a.InserirElemento(0, c + i - (n_a / T), sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
				else
					a.InserirElemento(0, ca + i - (n_a / T), sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
			}
			Alin.JuntarColuna(&a);
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizVmeta()
{
	int n_a = n;
	size_t R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c;
	c = (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + R + (sistema_a->GetTt2() - 1)*(n_a / T);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + r, 1);
			if (sistema_a->GetFlagVfol() == true)
				a.InserirElemento(0, c + r + (3+flag3)*R + 3*JJ + flag1*(sistema_a->barrasVtr.size()) + (1 - flag1), 1);
			Alin.JuntarColuna(&a);
		}
		c += (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T) + flag2*R;
	}

	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizBalPotencia()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t B = flag1*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4)*I + B;
	int jj = dd + 5*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, -1);
			}
			Alin.JuntarColuna(&a);
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		dd = dd + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizBalVazao()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t B = flag1*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4)*I + B + 2*R;
	int jj = dd + JJ + 3*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			a.InserirElemento(0, dd + R + r, - 1);
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, -1);
			}
			Alin.JuntarColuna(&a);
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		dd = dd + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizFuncProd()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t B = flag1*(sistema_a->barrasVtr.size() - 1);
	size_t I = sistema_a->termeletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int vv = (3+flag4)*I + B + R;
	int jj = vv + 4*R + JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
		for (size_t r = 0; r < R ; r++)
		{	
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, jj - JJ + j, 1);																		// phg
				a.InserirElemento(0, vv + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgV());					// v
				a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgQ());					// q
				a.InserirElemento(0, vv + R + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgD());				// d
				if (sistema_a->hidreletricasVtr[r].GetInflueVert() == 1 )		// Se o s influencia já tenho contabilizado isso no d, caso contrário só subtrair o valor de s de d
				{
				}
				else		// vertimento n influencia no canal de fuga, entao subtraio o valor de s da defluencia
				{
					a.InserirElemento(0, vv + 2*R + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgS());		// s
				}
				Alin.JuntarColuna(&a);
			}
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		vv = vv + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizPhMax()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t B = flag1*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4)*I + B + 4*R;
	int jj = dd + R + 2*JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
			}
			Alin.JuntarColuna(&a);
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		dd = dd + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizReserva()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t B = flag1*(sistema_a->barrasVtr.size() - 1);
	int nt;
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4)*I + B + 4*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		a.RemoverTodosElementos();
		for (size_t r = 0; r < R; r++)
		{
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			a.InserirElemento(0, dd - 4*R + r, - 1);
		}
		Alin.JuntarColuna(&a);
		dd = dd + nt;
	}
	return Alin;
}
CMatrizEsparsa Heuristicas::MatrizCortesF()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)		// Adiciona uma restrição por corte
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, - sistema_a->termeletricasVtr[i].GetCoefA1(n_cort));
				a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetCoefA0(n_cort));
				a.InserirElemento(0, c + i + 3*I, 1);
				Alin.JuntarColuna(&a);
			}
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Heuristicas::LimRestDemanda()
{
	int n_a = n;
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, B * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
		for (size_t b = 0; b < B; b++)
		{
			double D = 0;
			if ((0 <= t) && (sistema_a->GetTt1() > t))
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);		// cen = 0, ou qq outro cenários, pois os T1 primeiros periodos são iguais (deterministicos)
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
					for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
						D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
		//Lim[t*B + b][0] = 2;
		Lim[t*B + b][0] = 1;
		Lim[t*B + b][1] = D;
		}
	return Lim;
}
vetorfloat2 Heuristicas::LimRestDemandaBarraUnica()
{
	int n_a = n;
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		double D = 0;
		for (size_t b = 0; b < B; b++)
		{
			
			if ((0 <= t) && (sistema_a->GetTt1() > t))
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
					for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
						D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
		}
		Lim[t][0] = 1;
		Lim[t][1] = D;
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimFluxo0()
{
	int n_a = n;
	int L = sistema_a->linhasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, L * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
		for (int l = 0; l < L; l++)
		{
			Lim[t*L + l][0] = 0;
			Lim[t*L + l][1] = sistema_a->linhasVtr[l].GetCapacidade();
		}
	return Lim;
}
vetorfloat2 Heuristicas::LimFluxo2()
{
	int n_a = n;
	int L = sistema_a->linhasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, L * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
		for (int l = 0; l < L; l++)
		{
			Lim[t*L + l][0] = 2;
			Lim[t*L + l][1] = - sistema_a->linhasVtr[l].GetCapacidade();
		}
	return Lim;
}
vetorfloat2 Heuristicas::LimPhgMin()
{
	int n_a = n;
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				//Lim(jj + j,1) = double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Heuristicas::LimPhgMax()
{
	int n_a = n;
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 Heuristicas::LimBalHid()
{
	int n_a = n;
	double deltaT;
	int cenario;
	int periodo;
	int R = sistema_a->hidreletricasVtr.size();
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
				Lim[t*R + r][1] = double (sistema_a->hidreletricasVtr[r].GetV0() + 0.0036*deltaT*sistema_a->hidreletricasVtr[r].GetAfluencia(cenario,periodo));
			}
			else
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (0.0036*deltaT*sistema_a->hidreletricasVtr[r].GetAfluencia(cenario,periodo));
			}
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimQMin()
{
	int n_a = n;
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				//Lim(jj + j,1) = double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Heuristicas::LimQMax()
{
	int n_a = n;
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 Heuristicas::LimTup()
{
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	int Tup1, Tup2, Tup, t_a;
	//
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
	{
		if (t >= sistema_a->GetTt2())
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		else
			t_a = t;
		
		for (int i = 0; i < I; i++)
		{
			if (sistema_a->termeletricasVtr[i].GetTUp() < sistema_a->GetDeltaT1())
				Tup1 = 0;
			else
				Tup1 = int (sistema_a->termeletricasVtr[i].GetTUp() / sistema_a->GetDeltaT1());
			if (sistema_a->termeletricasVtr[i].GetTUp() < sistema_a->GetDeltaT2())
				Tup2 = 0;
			else
				Tup2 = int (sistema_a->termeletricasVtr[i].GetTUp() / sistema_a->GetDeltaT2());

			if (t < sistema_a->GetTt1())
				Tup = Tup1;
			else
				Tup = Tup2;

			if (Tup > 0)		// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
			{
				for (int tt = 0; tt < Tup; tt++)
				{
					if (t_a - tt <= 0)		// if referente ao elemento tt (se ele for constante entra nesse loop)
					{
						if (t_a - tt > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] = double (sistema_a->termeletricasVtr[i].GetU0());
							L[0] = 2;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] = double((1 - sistema_a->termeletricasVtr[i].GetU0()));
							L[0] = 2;
						}
					}
					else
					{
						L[1] = 0;
						L[0] = 2;
					}
					if (t_a - tt - 1 <= 0)	// if referente ao elemento tt - 1 (se ele for constante entra nesse loop)
					{
						if (t_a - tt - 1 > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] += - sistema_a->termeletricasVtr[i].GetU0();
							//L[0] = 2;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] += - (1 - sistema_a->termeletricasVtr[i].GetU0());
							//L[0] = 2;
						}
					}
					else
					{
						//L[1] = 0;
						//L[0] = 2;
					}
					Lim.push_back(L);
					L.clear();L.resize(2);
				}
			}
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimTdown()
{
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	int Tdown1, Tdown2, Tdown, t_a;
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 0;
	L[1] = 1;
	for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
	{
		if (t >= sistema_a->GetTt2())
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		else
			t_a = t;
		
		for (int i = 0; i < I; i++)
		{
			if (sistema_a->termeletricasVtr[i].GetTDown() < sistema_a->GetDeltaT1())
				Tdown1 = 0;
			else
				Tdown1 = int (sistema_a->termeletricasVtr[i].GetTDown() / sistema_a->GetDeltaT1());
			if (sistema_a->termeletricasVtr[i].GetTDown() < sistema_a->GetDeltaT2())
				Tdown2 = 0;
			else
				Tdown2 = int (sistema_a->termeletricasVtr[i].GetTDown() / sistema_a->GetDeltaT2());

			if (t < sistema_a->GetTt1())
				Tdown = Tdown1;
			else
				Tdown = Tdown2;

			if (Tdown > 0)		// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
			{
				for (int tt = 0; tt < Tdown; tt++)
				{
					if (t_a - tt <= 0)		// if referente ao elemento tt (se ele for constante entra nesse loop)
					{
						if (t_a - tt > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] += double (sistema_a->termeletricasVtr[i].GetU0());
							//L[0] = 0;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] += double((1 - sistema_a->termeletricasVtr[i].GetU0()));
							//L[0] = 0;
						}
					}
					else
					{
						//L[1] = 0;
						//L[0] = 0;
					}
					if (t_a - tt - 1 <= 0)	// if referente ao elemento tt - 1 (se ele for constante entra nesse loop)
					{
						if (t_a - tt - 1 > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] += - sistema_a->termeletricasVtr[i].GetU0();
							//L[0] = 0;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] += - (1 - sistema_a->termeletricasVtr[i].GetU0());
							//L[0] = 0;
						}
					}
					else
					{
						//L[1] = 0;
						//L[0] = 0;
					}
					Lim.push_back(L);
					L.clear();L.resize(2);L[0] = 0;L[1] = 1;
				}
			}
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimRampaUp()
{
	// T * I
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 0;
	L[1] = 1;
	for (int t = 0; t < T; t++)
	{
		for (int i = 0; i < I; i++)
		{
			if (t == 0)
				L[1] = sistema_a->termeletricasVtr[i].GetPmin() + sistema_a->termeletricasVtr[i].GetPt0() + sistema_a->termeletricasVtr[i].GetU0() * (sistema_a->termeletricasVtr[i].GetRampaUp() - sistema_a->termeletricasVtr[i].GetPmin());
			else
				L[1] = sistema_a->termeletricasVtr[i].GetPmin();
			Lim.push_back(L);
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimRampaDown()
{
	// T * I
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 0;
	L[1] = 1;
	for (int t = 0; t < T; t++)
	{
		for (int i = 0; i < I; i++)
		{
			if (t == 0)
				L[1] = sistema_a->termeletricasVtr[i].GetPmin() - sistema_a->termeletricasVtr[i].GetPt0();
			else
				L[1] = sistema_a->termeletricasVtr[i].GetPmin();
			Lim.push_back(L);
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimPtMin()
{
	// T * I
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 2;
	L[1] = 0;
	for (int t = 0; t < T; t++)
		for (int i = 0; i < I; i++)
			Lim.push_back(L);
	return Lim;
}
vetorfloat2 Heuristicas::LimPtMax()
{
	// T * I
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, I * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 Heuristicas::LimRestCP()
{
	// T * I
	int n_a = n;
	int I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 2;
	L[1] = 0;
	for (int t = 0; t < T; t++)
	{
		for (int i = 0; i < I; i++)
		{
			if (t == 0)
				L[1] = - sistema_a->termeletricasVtr[i].GetCoefCustoPartida()*sistema_a->termeletricasVtr[i].GetU0();
			else
				L[1] = 0;
			Lim.push_back(L);
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimVmeta()
{
	int n_a = n;
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * sistema_a->GetNCenarios(), 2);
	IniciaMatriz(&Lim, 0);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (int r = 0; r < R; r++)
		{
			Lim[r + R*cen][1] = sistema_a->hidreletricasVtr[r].GetVMeta();
			//Lim(r,1) = sistema_a->hidreletricasVtr[r].GetVmax();
			Lim[r + R*cen][0] = 2;
		}
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimBalPotenciaL()
{
	int n_a = n;
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Heuristicas::LimBalVazaoL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Heuristicas::LimFuncProdL()
{
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				Lim[jj + j][0] = 0;
				Lim[jj + j][1] = sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhg();
				//Lim[jj + j][1] = 0;		// FPH constate
			}
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Heuristicas::LimPhMax()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Heuristicas::LimReserva()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1 * T, 2);
	IniciaMatriz(&Lim, 0);
	double D;
	for (int t = 0; t < T; t++)
	{
		D = 0;
		for (size_t b = 0; b < B; b++)
		{
			if ((0 <= t) && (sistema_a->GetTt1() > t))
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
					for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
						D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
		}
		Lim[t][0] = 2;
		Lim[t][1] = D*sistema_a->GetFatorReserva()/100;
	}
	return Lim;
}
vetorfloat2 Heuristicas::LimCortesF()
{
	// T * I
	int I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 2;
	L[1] = 0;
	for (int t = 0; t < T; t++)
		for (int i = 0; i < I; i++)
			for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)
				Lim.push_back(L);
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void Heuristicas::MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes, int n_modelos)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);		//CMatrizEsparsa M;
	//posicao_rest_din.resize(6);		// seis elementos [ini fin], referentes as restricoes de Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida
	// na matriz esparsa percorrer elemento q n existe, em tese, n adiciona tempo algum, mas percorrer o RHS de tds as restricoes, por isso é bom saber os indices inicial e final das restricoes que sao atualizadas
	
	//posicao_rest_din = new vetorint[n_modelos];	// se esse atributo é inicializado aqui da problema no Release Build, td atributo (ponteiro) da classe deve ser inicializado no constructor
	//??? o problema esta nessa funçao, pois comentando ela n tem-se mais erro
	// colocar optimization maximize speed nas propriedas do projeto depois q encontrar a causa do problema
	// ver funçao criar restriçoes, n é necessario usar vetorfloat * ptr, somente um double * ptr[] qdo sabe-se o tamanho do vetor??? .size()??
		
	A = new CMatrizEsparsa[6];
	int * count = new int[n_modelos];
	//int ind_rest = 0;
	for (int i = 0; i < n_modelos; i++)
	{
		n_restricoes[i] = 0;
		count[i] = 0;
	}
	if ( sistema_a->GetFlagBarraUnica() == false )
	{
		M = MatrizRestDemanda();
		AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
		M = MatrizLimFluxo();
		AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
		M = MatrizLimFluxo();
		AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	}
	else
	{
		M = MatrizRestDemandaBarraUnica();
		AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	}
	M = MatrizLimPhgMin();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizLimPhgMax();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	//posicao_rest_din[0] = n_restricoes[0];
	for (int i = 0; i < n_modelos; i++)			// grava indice_inicio. Todos os subproblemas tem o msm numero de restriçoes antes de adicionar as restricoes de volume meta, portanto para contar o indice da restriçao pode-se usar n_restricoes do modelo 0
		posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizBalHid(); A[0] = M;
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizLimqMin();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizLimqMax();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizTup(); A[1] = M;									
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizTdown(); A[2] = M;
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizRampaUp(); A[3] = M;
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizRampaDown(); A[4] = M;
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizLimPtMin();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizLimPtMax();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizRestCP(); A[5] = M;
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	if ( sistema_a->GetFlagVfol() == true )
	{
		M = MatrizVmeta();
		AlocarRestricoesVmeta(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	}
	M = MatrizBalPotencia();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizBalVazao();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizFuncProd();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizPhMax();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	M = MatrizReserva();
	AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		M = MatrizCortesF();
		AlocarRestricoes(&M, indexL, indexC, indexV, nnZ, n_restricoes, count, n_modelos);
	}
}
void Heuristicas::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor)
{
	vetorfloat2 L;
	lim_iniciais.resize(6);// = new vetorfloat2[6];		// Guardar valores dos limites das restriçoes originais
	if ( sistema_a->GetFlagBarraUnica() == false )
	{
		L = LimRestDemanda();
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimFluxo0();
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimFluxo2();
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
	else
	{
		L = LimRestDemandaBarraUnica();
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
	L = LimPhgMin();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimPhgMax();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimBalHid(); lim_iniciais[0] = L;
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimQMin();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimQMax();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimTup(); lim_iniciais[1] = L;
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimTdown(); lim_iniciais[2] = L;
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimRampaUp(); lim_iniciais[3] = L;
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimRampaDown(); lim_iniciais[4] = L;
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimPtMin();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimPtMax();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimRestCP(); lim_iniciais[5] = L;
	AlocarLimitesSubp(L, LimTipo, LimValor);
	if ( sistema_a->GetFlagVfol() == true )
	{
		L = LimVmeta();
		AlocarLimitesVmetaSubp(L, LimTipo, LimValor);
	}
	L = LimBalPotenciaL();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimBalVazaoL();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimFuncProdL();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimPhMax();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimReserva();
	AlocarLimitesSubp(L, LimTipo, LimValor);
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		L = LimCortesF();
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
}

Heuristicas::Heuristicas(CSistema * const sistema_end, ResultadosConj * const resultadosGurobi_end, double max_time) : ambGRB(GRBEnv())
{
	sistema_a = sistema_end;
	resultadosGurobi = resultadosGurobi_end;
	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	flag1 = int (1 - sistema_a->GetFlagBarraUnica());	
	flag2 = int (sistema_a->GetFlagVfol());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n = T * ((3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	if ( sistema_a->GetFlagBarraUnica() == true )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	x.resize(n);
	x_hat.resize(n);
	x_til.resize(n);
	X = NULL;
	fo = GRB_INFINITY;

	// parametros
	alfa = ALFA;
	beta = BETA;
	gama = GAMA;

	ws = WIN_SIZE;
	wlas = WIN_FUT_SIZE;
	nw = sistema_a->GetTt2() / ws;
	if ( ws > sistema_a->GetTt2() || ws + wlas > sistema_a->GetTt2() )
		cout << "Atençao: numero de janelas incompatível com numero de periodos considerados!!!" << endl;

	flagH1 = BIN_MW;
	flagH2 = BIN_FW;
	flagH3 = NORM_P;

	//if ( (NORM_P != 1) && ((beta == 0 && gama == 1) || (beta == 1 && gama == 0))  )
	//	flagH3 = 1;		// Usar norma 1 sempre que beta=0 e gama=1 ou beta=1 e gama=0
	if (NORM_P == 1)
	{
		if ( ((beta != 0 || gama != 1) && (beta != 1 || gama != 0)) )
		{
			cout << "Norma linear só está implementada para o termo prox. das var. binárias" << endl;
			cout << "Norma alterada para quadratica" << endl;
			flagH3 = 2;
		}
	}

	// Criar modelos do gurobi (variaveis e restrições)
	CriarModelos(max_time);
}
Heuristicas::Heuristicas(CSistema * const sistema_end, ResultadosConj * const resultadosGurobi_end, int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a, double max_time) : ambGRB(GRBEnv())
{
	sistema_a = sistema_end;
	resultadosGurobi = resultadosGurobi_end;
	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	flag1 = int (1 - sistema_a->GetFlagBarraUnica());	
	flag2 = int (sistema_a->GetFlagVfol());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n = T * ((3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	if ( sistema_a->GetFlagBarraUnica() == true )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	x.resize(n);
	x_hat.resize(n);
	x_til.resize(n);
	X = NULL;
	fo = GRB_INFINITY;

	// parametros
	alfa = alfa_a;
	beta = beta_a;
	gama = gama_a;

	ws = ws_a;
	wlas = wlas_a;
	nw = sistema_a->GetTt2() / ws;
	if ( ws > sistema_a->GetTt2() || ws + wlas > sistema_a->GetTt2() )
		cout << "Atençao: numero de janelas incompatível com numero de periodos considerados!!!" << endl;

	flagH1 = bvmw;
	flagH2 = bvfw;
	flagH3 = 2;

	if ( (NORM_P != 1) && ((beta == 0 && gama == 1) || (beta == 1 && gama == 0))  )
		flagH3 = 1;		// Usar norma 1 sempre que beta=0 e gama=1 ou beta=1 e gama=0

	// Criar modelos do gurobi (variaveis e restrições)
	CriarModelos(max_time);
}
Heuristicas::~Heuristicas(void)
{
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
	for (size_t i = 0; i < constr.size(); i++)
		delete constr[i];
	constr.clear();
}
void Heuristicas::CriarModelos(double max_time)
{
	// Os subproblemas sao divididos em 4 categorias:
	int n_subp1 = sistema_a->GetTt1() / ws;		// subproblemas somente de primeiro estágio
	int n_subpA = ( sistema_a->GetTt1() % ws > 0 );		// subproblema de acoplamento (de 1º e 2º est.)
	int n_subp2 = (((sistema_a->GetTt2() - sistema_a->GetTt1()) - (ws - (sistema_a->GetTt1()%ws))*n_subpA) / ws);		// subproblema de segundo estágio somente (um para cada cenário)
	int n_subpFH = (((sistema_a->GetTt2() - sistema_a->GetTt1()) - (ws - (sistema_a->GetTt1()%ws))*n_subpA) % ws > 0 );		// mesmo que o anterior porém sao os problema do fim do horizonte (com tamanhos diferentes)

	// contar qts problemas estocásticos???
	n_modelos = n_subp1 + n_subpA + (n_subp2 + n_subpFH)*sistema_a->GetNCenarios();

	vtr_nos_pri.resize(n_modelos);		// vetor com os numeros dos nos principais de cada subproblema
	vtr_nos_fut.resize(n_modelos);		// vetor com os numeros dos nos futuro de cada subproblema
	int no = 0;
	int cc = 0;
	int cenario = 0;
	if (ws >= sistema_a->GetTt2())
		for (int i = 0; i < T; i++)
			vtr_nos_pri[0].push_back(i);
	else
	{
		while (no < T)
		{
			if (no >= sistema_a->GetTt1())
				cenario = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1()) + 1;
			if ( (no + ws - 1 >= sistema_a->GetTt1() ) && ( no < sistema_a->GetTt1()) )		// condiçao dos subproblemas estocásticos (janelas principais est.)
			{
				for (int i = no; i < sistema_a->GetTt1(); i++)
					vtr_nos_pri[cc].push_back(i);
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				{
					for (int i = sistema_a->GetTt1(); i < no + ws ; i++)
						vtr_nos_pri[cc].push_back(i + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
					for (int i = no + ws; i < min(no + ws + wlas, sistema_a->GetTt2()) ; i++)		// adicionar tb janelas futuro
						vtr_nos_fut[cc].push_back(i + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
				}
				no += ws;
			}
			else if ( (no >= sistema_a->GetTt1()) && ( (no + ws - 1) >= cenario * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1() ) )	// condiçao dos subproblemas de fim do horizonte
			{
				for (int i = no; i < cenario*(sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1(); i++)
					vtr_nos_pri[cc].push_back(i);
				no = cenario*(sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1() + ws - (sistema_a->GetTt1() % ws);
			}
			else		// condição geral
			{
				for (int i = 0; i < ws; i++)
					vtr_nos_pri[cc].push_back(no + i);
				if ( (no + ws + wlas - 1>= sistema_a->GetTt1() ) && ( no + ws - 1 < sistema_a->GetTt1()) )		// condiçao dos subproblemas estocásticos (janelas futuras est.)
				{
					for (int i = no + ws; i < sistema_a->GetTt1(); i++)
						vtr_nos_fut[cc].push_back(i);
					for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
						for (int i = sistema_a->GetTt1(); i < min(no + ws + wlas, sistema_a->GetTt2()) ; i++)
							vtr_nos_fut[cc].push_back(i + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
				}
				else if ( (no >= sistema_a->GetTt1()) && ( (no + ws + wlas - 1) >= cenario * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1() ) )		// condiçao dos subproblemas de fim do horizonte para as janelas futuras
				{
					for (int i = no + ws; i < cenario*(sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->GetTt1(); i++)
						vtr_nos_fut[cc].push_back(i);
				}
				else
				{
					for (int i = 0; i < wlas; i++)
						vtr_nos_fut[cc].push_back(no + ws + i);
				}
				no += ws;
			}
			cc++;
		}
	}
	//// Escrever vetores (debug)
	//for (int i = 0; i < n_modelos; i++)
	//{
	//	for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)
	//		cout << vtr_nos_pri[i][ii] << " ;";
	//	cout << "." << endl;
	//}
	//cout << endl << endl;
	//for (int i = 0; i < n_modelos; i++)
	//{
	//	for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)
	//		cout << vtr_nos_fut[i][ii] << " ;";
	//	cout << "." << endl;
	//}
	modelosGRB.resize(n_modelos);
	vars.resize(n_modelos);
	constr.resize(n_modelos);
	numero_var_subp.resize(n_modelos);
	numero_var_subp_principal.resize(n_modelos);
	posicao_rest_din.resize(n_modelos);


	X = new CMatrizEsparsa(n, 1);
	//A = new CMatrizEsparsa[6];
	//lim_iniciais = new vetorfloat2[6];


	for (int i = 0; i < n_modelos; i++)
	{
		modelosGRB[i] = new GRBModel(ambGRB);		// um modelo de otimização para cada usina
		numero_var_subp[i] = 0;
		numero_var_subp_principal[i] = 0;
	}

	// numero de variaveis dos subproblemas
	int nt;
	for (int i = 0; i < n_modelos; i++)
	{
		for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)
		{
			if ( (flag2 == 1) && (vtr_nos_pri[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
				nt = ( n - sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T + sistema_a->hidreletricasVtr.size();
			else
				nt = n / T;
			numero_var_subp[i] += nt;
			numero_var_subp_principal[i] += nt;
		}
		for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)
		{
			if ( (flag2 == 1) && (vtr_nos_fut[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// tem-se a variável vfol nesse nó!!
				nt = ( n - sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T + sistema_a->hidreletricasVtr.size();
			else
				nt = n / T;
			numero_var_subp[i] += nt;
		}
	}

	// Criar variáveis
	for (int i = 0; i < n_modelos; i++)
	{
		CriarVariaveis(i);
		modelosGRB[i]->update();	// Atualiza o modelo Gurobi.
		vars[i] = modelosGRB[i]->getVars();
	}

	// Adicionar restrições (Monta matrizes do problema inteiro e selecionar conjunto de restriçoes para cada subproblema)
	CriarRestricoes();

	// Definir ajustes do solver
	for (int i = 0; i < n_modelos; i++)
	{
		//modelosGRB[i]->write("Hrst_subprob.lp");	// Escreve modelo em arquivo
		modelosGRB[i]->getEnv().set(GRB_IntParam_OutputFlag, 0);		// Nao escrever detalhes
		//modelosGRB[i].getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[i].getEnv().set(GRB_StringParam_LogFile, "logHrst.txt");		// Escrever log em arquivo
		modelosGRB[i]->getEnv().set(GRB_DoubleParam_MIPGap, 0.01);	// Define o gap de tolerancia
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6); 
		//modelosGRB[i]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm.
		//modelosGRB[i]->getEnv().set(GRB_IntParam_Threads, 1);
		//modelosGRB[i]->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[i].getEnv().set(GRB_IntParam_Method, 0);
		modelosGRB[i]->getEnv().set(GRB_DoubleParam_TimeLimit, max_time*TIME_R/n_modelos);		// Limita o tempo de resolução do problema
	}
}
void Heuristicas::CriarVariaveis(int modelo)
{
	try 
	{
		size_t I = sistema_a->termeletricasVtr.size();
		size_t R = sistema_a->hidreletricasVtr.size();
		size_t B = sistema_a->barrasVtr.size() - 1;
		int JJ = 0;
		for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
			JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
		int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
		int t = 0;
		// Variáveis das janelas principais
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
		{
			int cen = 0;
			if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
				t = vtr_nos_pri[modelo][vtr_i] - cen*(sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t = vtr_nos_pri[modelo][vtr_i];
			//x = [pt u cp F teta ph v d s phmax phg q z def]
			for (size_t i = 0; i < I; i++)	//pt
				modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			if  ( (flagH1 == 0) || (flagH1 == 1) )	//u (fixo e livre continuo ou livre binario)
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			for (size_t i = 0; i < I; i++)	//cp
				modelosGRB[modelo]->addVar(0, sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F		//if (sistema_a->GetFlagInitAproxCT() > 1)
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagBarraUnica() == false)
				for (size_t b = 0; b < B; b++)	//teta
					modelosGRB[modelo]->addVar(-6.2832, 6.2832, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modelosGRB[modelo]->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin()), double(sistema_a->hidreletricasVtr[r].GetVmax()), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				modelosGRB[modelo]->addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
			{
				modelosGRB[modelo]->addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < R; r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modelosGRB[modelo]->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modelosGRB[modelo]->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( (flagH1 == 0) || (flagH1 == 1) )	//z (fixo e livre continuo ou livre binario)
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			if ( sistema_a->GetFlagBarraUnica() == false)
			{
				for (size_t b = 0; b < B + 1; b++)	//def
				{
					double cap_d = 0;
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
					modelosGRB[modelo]->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
					//modelosGRB[modelo]->addVar(-cap_d, 0, 0.0, GRB_CONTINUOUS, "");
					//modelosGRB[modelo]->addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				double cap_d = 0;			//def barra unica
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
				modelosGRB[modelo]->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
			}
			if ( (t + 1 > sistema_a->GetTt1()) && ((t + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
				if ( sistema_a->GetFlagVfol() == true )
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//vfol
						modelosGRB[modelo]->addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
						//modelosGRB[modelo]->addVar(0, 0.0000001, 0.0, GRB_CONTINUOUS, "");
		}

		// Variáveis das janelas futuras
		for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)
		{
			int cen = 0;
			if ( vtr_nos_fut[modelo][vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_nos_fut[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			if ( vtr_nos_fut[modelo][vtr_i] >= sistema_a->GetTt2() )
				t = vtr_nos_fut[modelo][vtr_i] - cen*(sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t = vtr_nos_fut[modelo][vtr_i];
			//x = [pt u cp F teta ph v d s phmax phg q z def]
			for (size_t i = 0; i < I; i++)	//pt
				modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			if ( (flagH2 == 0) || (flagH2 == 1) ) //u (fixo e livre continuo ou livre binario)
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			for (size_t i = 0; i < I; i++)	//cp
				modelosGRB[modelo]->addVar(0, sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F		//if (sistema_a->GetFlagInitAproxCT() > 1)
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagBarraUnica() == false)
				for (size_t b = 0; b < B; b++)	//teta
					modelosGRB[modelo]->addVar(-6.2832, 6.2832, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modelosGRB[modelo]->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin()), double(sistema_a->hidreletricasVtr[r].GetVmax()), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				modelosGRB[modelo]->addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
			{
				modelosGRB[modelo]->addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < R; r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modelosGRB[modelo]->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modelosGRB[modelo]->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( (flagH2 == 0) || (flagH2 == 1) )	//z (fixo e livre continuo ou livre binario)
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			if ( sistema_a->GetFlagBarraUnica() == false)
			{
				for (size_t b = 0; b < B + 1; b++)	//def
				{
					double cap_d = 0;
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
					modelosGRB[modelo]->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
					//modelosGRB[modelo]->addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				double cap_d = 0;			//def barra unica
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
				modelosGRB[modelo]->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
			}
			if ( (t + 1 > sistema_a->GetTt1()) && ((t + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
				if ( sistema_a->GetFlagVfol() == true )
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//vfol
						modelosGRB[modelo]->addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
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
void Heuristicas::CriarRestricoes()
{
	try
	{
		int * n_restricoes = new int[n_modelos];
		vetorint * indexL = new vetorint[n_modelos];
		vetorint * indexC = new vetorint[n_modelos];
		vetorint * nnZ = new vetorint[n_modelos];
		vetorint * LimTipo = new vetorint[n_modelos];
		vetorfloat * indexV = new vetorfloat[n_modelos];
		vetorfloat * LimValor = new vetorfloat[n_modelos];
		//vetorint LimTipo;	// Tds iguais ao primeiro (subp = 0) subproblema (depois, antes de resolver cada subproblema atualizam-se o RHS, right-hand side, das restriçoes para subp > 0)
		//vetorfloat LimValor;
		for (int i = 0; i < n_modelos; i++)
			nnZ[i].push_back(0);

		MatrizRestricoesLineares(indexL, indexC, indexV, nnZ, n_restricoes, n_modelos);
		MatrizLimitesLineares(LimTipo, LimValor);

		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		for (int mod = 0; mod < n_modelos; mod++)
		{
			for (int i = 0; i < n_restricoes[mod]; i++)
			{
				for (int j = nnZ[mod][i]; j < nnZ[mod][i + 1]; j++)
				{
					coeficiente = indexV[mod][j];
					variavel = vars[mod][indexC[mod][j]];
					restricao.addTerms( &coeficiente, &variavel, 1);
				}
				switch (LimTipo[mod][i]) {
				case 0:
					modelosGRB[mod]->addConstr(restricao, GRB_LESS_EQUAL, LimValor[mod][i], "");
					break;
				case 1:
					modelosGRB[mod]->addConstr(restricao, GRB_EQUAL, LimValor[mod][i], "");
					break;
				case 2:
					modelosGRB[mod]->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[mod][i], "");
					break;
				default:
					cout << "Tipo inválido de restrição adicionada" << endl;
				}
				restricao.clear();
			}
			modelosGRB[mod]->update();	// Atualiza o modelo Gurobi
			constr[mod] = modelosGRB[mod]->getConstrs();
		}
		delete n_restricoes;
		delete[] indexL;
		delete[] indexC;
		delete[] nnZ;
		delete[] indexV;
		delete[] LimTipo;
		delete[] LimValor;
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
void Heuristicas::AlocarRestricoes(CMatrizEsparsa * matriz, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes, int * count, int n_modelos)
{
	// Selecionar colunas
	size_t R = sistema_a->hidreletricasVtr.size();
	int nt, col_i, lin_i, cenario;
	int nt_a = ( n - flag2*sistema_a->GetNCenarios()*R ) / T;
	int nt_rest = matriz->GetNlin() / T;
	for (int i = 0; i < n_modelos; i++)		// Uma matriz para cada modelo
	{
		CMatrizEsparsa M_subp( nt_rest*(vtr_nos_pri[i].size() + vtr_nos_fut[i].size()), numero_var_subp[i]);
		col_i = 0;
		for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)	// janelas principais
		{
			if (vtr_nos_pri[i][ii] >= sistema_a->GetTt2())
				cenario = (vtr_nos_pri[i][ii] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			else 
				cenario = 0;
			if ( (flag2 == 1) && (vtr_nos_pri[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
				nt = ( n - sistema_a->GetNCenarios()*R ) / T + R;
			else
				nt = n / T;
			// nt tem o numero de variáveis do nó vtr_nos_pri[i][ii]
			// nt_a tem o numero de variaveis dos nós que não sao do final de horizonte

			// dois loops, o primeiro é para varrer as colunas (ii) e o segundo para percorrer as linhas (iii) da matriz A
			lin_i = ii*nt_rest;
			for (size_t iii = ii; iii < vtr_nos_pri[i].size(); iii++)
			{
				M_subp.InserirMatriz(lin_i, col_i, lin_i + nt_rest - 1, col_i + nt - 1, matriz, nt_rest*vtr_nos_pri[i][iii], nt_a*vtr_nos_pri[i][ii] + flag2*cenario*R);
				lin_i += nt_rest;
			}
			for (size_t iii = 0; iii < vtr_nos_fut[i].size(); iii++)	// segundo loop para as linhas das janelas futuras
			{
				M_subp.InserirMatriz(lin_i, col_i, lin_i + nt_rest - 1, col_i + nt - 1, matriz, nt_rest*vtr_nos_fut[i][iii], nt_a*vtr_nos_pri[i][ii] + flag2*cenario*R);
				lin_i += nt_rest;
			}
			col_i += nt;
		}
		for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)	// janelas futuras
		{
			if (vtr_nos_fut[i][ii] >= sistema_a->GetTt1())
				cenario = (vtr_nos_fut[i][ii] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			else 
				cenario = 0;
			if ( (flag2 == 1) && (vtr_nos_fut[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
				nt = ( n - sistema_a->GetNCenarios()*R ) / T + R;
			else
				nt = n / T;
			
			lin_i = nt_rest*vtr_nos_pri[i].size() + ii*nt_rest;
			for (size_t iii = ii; iii < vtr_nos_fut[i].size(); iii++)
			{
				M_subp.InserirMatriz(lin_i, col_i, lin_i + nt_rest - 1, col_i + nt - 1, matriz, nt_rest*vtr_nos_fut[i][iii], nt_a*vtr_nos_fut[i][ii] + flag2*cenario*R);
				lin_i += nt_rest;
			}
			col_i += nt;
		}
		SparseMatriz(&M_subp, &indexL[i], &indexC[i], &indexV[i], &nnZ[i], &n_restricoes[i], &count[i]);
		n_restricoes[i] += M_subp.GetNlin();
	}
}
void Heuristicas::AlocarRestricoesVmeta(CMatrizEsparsa * matriz, vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes, int * count, int n_modelos)
{
	// Selecionar colunas
	size_t R = sistema_a->hidreletricasVtr.size();
	int nt, col_i, cenario;
	int nt_a = ( n - flag2*sistema_a->GetNCenarios()*R ) / T;
	for (int i = 0; i < n_modelos; i++)		// Uma matriz para cada modelo
	{
		CMatrizEsparsa M_subp(0, numero_var_subp[i]);
		CMatrizEsparsa M_subp_a(R, numero_var_subp[i]);
		col_i = 0;
		for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)	// janelas principais
		{
			if (vtr_nos_pri[i][ii] >= sistema_a->GetTt2())
				cenario = (vtr_nos_pri[i][ii] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			else 
				cenario = 0;
			if ( (vtr_nos_pri[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
			{
				nt = ( n - sistema_a->GetNCenarios()*R ) / T + R;
				M_subp_a.InserirMatriz(0, col_i, R - 1, col_i + nt - 1, matriz, cenario*R, nt_a*vtr_nos_pri[i][ii] + flag2*cenario*R);
				M_subp.JuntarColuna(&M_subp_a);
				M_subp_a.RemoverTodosElementos();
			}
			else
				nt = n / T;
			col_i += nt;
		}
		for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)	// janelas futuras
		{
			if (vtr_nos_fut[i][ii] >= sistema_a->GetTt1())
				cenario = (vtr_nos_fut[i][ii] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			else 
				cenario = 0;
			if ( (vtr_nos_fut[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
			{
				nt = ( n - sistema_a->GetNCenarios()*R ) / T + R;
				M_subp_a.InserirMatriz(0, col_i, R - 1, col_i + nt - 1, matriz, cenario*R, nt_a*vtr_nos_fut[i][ii] + flag2*cenario*R);
				M_subp.JuntarColuna(&M_subp_a);
				M_subp_a.RemoverTodosElementos();
			}
			else
				nt = n / T;
			col_i += nt;
		}
		if (M_subp.GetNnz() != 0)		// nao incluir matriz vazia para nós q n sao do final do horizonte!!!
		{
			SparseMatriz(&M_subp, &indexL[i], &indexC[i], &indexV[i], &nnZ[i], &n_restricoes[i], &count[i]);
			n_restricoes[i] += M_subp.GetNlin();
		}
	}
}
void Heuristicas::AlocarLimitesSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor)
{
	// Selecionar colunas
	int lin_i;
	int nt_rest = limites.size() / T;
	for (int i = 0; i < n_modelos; i++)		// Uma matriz para cada modelo
	{
		vetorfloat2 L;
		L.resize(nt_rest*(vtr_nos_pri[i].size() + vtr_nos_fut[i].size()));
		lin_i = 0;
		for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)	// janelas principais
		{
			for (int iii = 0; iii < nt_rest; iii++)
			{
				L[lin_i + iii].push_back(limites[nt_rest*vtr_nos_pri[i][ii] + iii][0]);
				L[lin_i + iii].push_back(limites[nt_rest*vtr_nos_pri[i][ii] + iii][1]);
			}
			lin_i += nt_rest;
		}
		for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)	// janelas futuras
		{
			for (int iii = 0; iii < nt_rest; iii++)
			{
				L[lin_i + iii].push_back(limites[nt_rest*vtr_nos_fut[i][ii] + iii][0]);
				L[lin_i + iii].push_back(limites[nt_rest*vtr_nos_fut[i][ii] + iii][1]);
			}
			lin_i += nt_rest;
		}
		AlocarLimites(&L, &LimTipo[i], &LimValor[i]);
	}
}
void Heuristicas::AlocarLimitesVmetaSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor)
{
	// Selecionar colunas
	size_t R = sistema_a->hidreletricasVtr.size();
	int cenario, lin_i;
	for (int i = 0; i < n_modelos; i++)		// Uma matriz para cada modelo
	{
		vetorfloat2 L;
		lin_i = 0;
		for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)	// janelas principais
		{
			if (vtr_nos_pri[i][ii] >= sistema_a->GetTt2())
				cenario = (vtr_nos_pri[i][ii] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			else 
				cenario = 0;
			if ( (vtr_nos_pri[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
			{
				L.resize(lin_i + R);
				for (int iii = 0; iii < R; iii++)
				{
					L[lin_i + iii].push_back(limites[cenario*R + iii][0]);
					L[lin_i + iii].push_back(limites[cenario*R + iii][1]);
				}
				lin_i += R;
			}
		}
		for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)	// janelas futuras
		{
			if (vtr_nos_fut[i][ii] >= sistema_a->GetTt2())
				cenario = (vtr_nos_fut[i][ii] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			else 
				cenario = 0;
			if ( (vtr_nos_fut[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
			{
				L.resize(lin_i + R);
				for (int iii = 0; iii < R; iii++)
				{
					L[lin_i + iii].push_back(limites[cenario*R + iii][0]);
					L[lin_i + iii].push_back(limites[cenario*R + iii][1]);
				}
				lin_i += R;
			}
		}
		if (lin_i != 0)		// tem nó de final de horizonte no subproblema!!
			AlocarLimites(&L, &LimTipo[i], &LimValor[i]);
	}
}
void Heuristicas::CriarFuncaoObjetivo(int modelo)
{
	if ( flagH3 == 1 )
		CriarFuncaoObjetivoL(modelo);
	else
		CriarFuncaoObjetivoQ(modelo);
}
void Heuristicas::CriarFuncaoObjetivoQ(int modelo)
{
	// Modela a f.o. do subproblema com termo proximal quadrático
	try 
	{
		int R = sistema_a->hidreletricasVtr.size();
		int I = sistema_a->termeletricasVtr.size();
		int B = sistema_a->barrasVtr.size();
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
				
		// Funçao objetivo normalizada com os valores de cada subproblema
		//CalcularFuncaoObjetivo(modelo) -> só calcula a f.o. referente as janelas principais!!
		//double optimal_LR;
		double optimal_LR = CalcularFuncaoObjetivoRL(modelo);		//calcula a f.o. referente as janelas principais e futuras!!

		int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
		int cen, jj;
		double deltaT;
		int delta = 0;
		int delta_x = 0;
		double constante = 0;
		// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
		vetorint vtr_a;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
			vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
		for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
			vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);
		
		// Variáveis das janelas principais e futuras
		//x = [pt u cp (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		// Termos lineares
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal
		{
			cen = 0;
			if ( vtr_a[vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			delta_x = nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				deltaT = sistema_a->GetDeltaT1();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_Obj, 1*deltaT * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[i + 3*I + delta_x] + (1-beta)*(1-gama)*x_hat[i + 3*I + delta_x])/pow(vars[modelo][i + 3*I + delta].get(GRB_DoubleAttr_UB), 2));	//F
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2));				//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*(1-gama)*x_til[i + I + delta_x] + (1-beta)*gama*x_hat[i + I + delta_x]));		//u
						constante += (1 - alfa)*(beta*gama*pow(x_til[i + 3*I + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + 3*I + delta_x], 2))/pow(vars[modelo][i + 3*I + delta].get(GRB_DoubleAttr_UB), 2);		//F
					}
					else
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2));	//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * alfa/optimal_LR - 2*(1-alfa)*(beta*(1-gama)*x_til[i + I + delta_x] + (1-beta)*gama*x_hat[i + I + delta_x]));		//u
					}
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, 1 * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[i + 2*I + delta_x] + (1-beta)*(1-gama)*x_hat[i + 2*I + delta_x])/pow(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB), 2));		//cp
					constante += (1 - alfa)*(beta*gama*pow(x_til[i + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + delta_x], 2))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2);				//pt
					constante += (1 - alfa)*(beta*(1-gama)*pow(x_til[i + I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + I + delta_x], 2));	//u
					constante += (1 - alfa)*(beta*gama*pow(x_til[i + 2*I + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + 2*I + delta_x], 2))/pow(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB), 2);		//cp
				}
				if (flag1 == 1)		// considerar rede
				{
					for (size_t b = 0; b < B - 1; b++)
					{
						vars[modelo][b + I*(3 + flag4) + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[b + I*(3 + flag4) + delta_x] + (1-beta)*(1-gama)*x_hat[b + I*(3 + flag4) + delta_x])/pow(vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_LB), 2));		//teta
						constante += (1 - alfa)*(beta*gama*pow(x_til[b + I*(3 + flag4) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[b + I*(3 + flag4) + delta_x], 2))/pow(vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_LB), 2);		//teta
					}
					for (size_t b = 0; b < B; b++)
					{
						if ( vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) == 0 )
							vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa/optimal_LR);		//def
						else
						{
							vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x] + (1-beta)*(1-gama)*x_hat[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x])/pow(vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2));		//def
							constante += (1 - alfa)*(beta*gama*pow(x_til[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x], 2))/pow(vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2);		//def
						}
					}
				}
				else
				{
					vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[I*(3 + flag4) + 5*R + 3*JJ + delta_x] + (1-beta)*(1-gama)*x_hat[I*(3 + flag4) + 5*R + 3*JJ + delta_x])/pow(vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2));		//def
					constante += (1 - alfa)*(beta*gama*pow(x_til[I*(3 + flag4) + 5*R + 3*JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[I*(3 + flag4) + 5*R + 3*JJ + delta_x], 2))/pow(vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2);		//def
				}
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].get(GRB_DoubleAttr_UB), 2));				//ph
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_LB), 2));		//v
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB), 2));	//d
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB), 2));	//s
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB), 2));	//phmax
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].get(GRB_DoubleAttr_UB), 2);			//ph
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_LB), 2);		//v
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB), 2);	//d
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB), 2);	//s
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB), 2);	//phmax
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						vars[modelo][j + jj + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[j + jj + delta_x] + (1-beta)*(1-gama)*x_hat[j + jj + delta_x])/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB), 2));							//phg
						vars[modelo][j + jj + JJ + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[j + jj + JJ + delta_x] + (1-beta)*(1-gama)*x_hat[j + jj + JJ + delta_x])/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB), 2));		//q
						vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*(1-gama)*x_til[j + jj + 2*JJ + delta_x] + (1-beta)*gama*x_hat[j + jj + 2*JJ + delta_x]));																//z
						constante += (1 - alfa)*(beta*gama*pow(x_til[j + jj + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[j + jj + delta_x], 2))/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB), 2);					//phg
						constante += (1 - alfa)*(beta*gama*pow(x_til[j + jj + JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[j + jj + JJ + delta_x], 2))/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB), 2);		//q
						constante += (1 - alfa)*(beta*(1-gama)*pow(x_til[j + jj + 2*JJ + delta_x], 2) + (1-beta)*gama*pow(x_hat[j + jj + 2*JJ + delta_x], 2));	//z
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				delta += nt;
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_Obj, 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[i + 3*I + delta_x] + (1-beta)*(1-gama)*x_hat[i + 3*I + delta_x])/pow(vars[modelo][i + 3*I + delta].get(GRB_DoubleAttr_UB), 2));	//F
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2));					//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*(1-gama)*x_til[i + I + delta_x] + (1-beta)*gama*x_hat[i + I + delta_x]));		//u
						constante += (1 - alfa)*(beta*gama*pow(x_til[i + 3*I + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + 3*I + delta_x], 2))/pow(vars[modelo][i + 3*I + delta].get(GRB_DoubleAttr_UB), 2);		//F
					}
					else
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) * alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2));	//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) * alfa/optimal_LR - 2*(1-alfa)*(beta*(1-gama)*x_til[i + I + delta_x] + (1-beta)*gama*x_hat[i + I + delta_x]));		//u
					}
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, 1*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[i + 2*I + delta_x] + (1-beta)*(1-gama)*x_hat[i + 2*I + delta_x])/pow(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB), 2));		//cp
					constante += (1 - alfa)*(beta*gama*pow(x_til[i + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + delta_x], 2))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2);				//pt
					constante += (1 - alfa)*(beta*(1-gama)*pow(x_til[i + I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + I + delta_x], 2));	//u
					constante += (1 - alfa)*(beta*gama*pow(x_til[i + 2*I + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + 2*I + delta_x], 2))/pow(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB), 2);		//cp
				}
				if (flag1 == 1)		// considerar rede
				{
					for (size_t b = 0; b < B - 1; b++)
					{
						vars[modelo][b + I*(3 + flag4) + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[b + I*(3 + flag4) + delta_x] + (1-beta)*(1-gama)*x_hat[b + I*(3 + flag4) + delta_x])/pow(vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_LB), 2));		//teta
						constante += (1 - alfa)*(beta*gama*pow(x_til[b + I*(3 + flag4) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[b + I*(3 + flag4) + delta_x], 2))/pow(vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_LB), 2);		//teta
					}
					for (size_t b = 0; b < B; b++)
					{
						if ( vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) == 0 )
							vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR);		//def
						else
						{
							vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x] + (1-beta)*(1-gama)*x_hat[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x])/pow(vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2));		//def
							constante += (1 - alfa)*(beta*gama*pow(x_til[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta_x], 2))/pow(vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2);		//def
						}
					}
				}
				else
				{
					vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[I*(3 + flag4) + 5*R + 3*JJ + delta_x] + (1-beta)*(1-gama)*x_hat[I*(3 + flag4) + 5*R + 3*JJ + delta_x])/pow(vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2));		//def
					constante += (1 - alfa)*(beta*gama*pow(x_til[I*(3 + flag4) + 5*R + 3*JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[I*(3 + flag4) + 5*R + 3*JJ + delta_x], 2))/pow(vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2);		//def
				}
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].get(GRB_DoubleAttr_UB), 2));							//ph
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_LB), 2));		//v
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB), 2));	//d
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB), 2));	//s
					vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB), 2));	//phmax
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].get(GRB_DoubleAttr_UB), 2);					//ph
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_LB), 2);		//v
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB), 2);	//d
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB), 2);	//s
					constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB), 2);	//phmax
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						vars[modelo][j + jj + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[j + jj + delta_x] + (1-beta)*(1-gama)*x_hat[j + jj + delta_x])/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB), 2));							//phg
						vars[modelo][j + jj + JJ + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*gama*x_til[j + jj + JJ + delta_x] + (1-beta)*(1-gama)*x_hat[j + jj + JJ + delta_x])/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB), 2));		//q
						vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, - 2*(1-alfa)*(beta*(1-gama)*x_til[j + jj + 2*JJ + delta_x] + (1-beta)*gama*x_hat[j + jj + 2*JJ + delta_x]));																//z
						constante += (1 - alfa)*(beta*gama*pow(x_til[j + jj + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[j + jj + delta_x], 2))/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB), 2);					//phg
						constante += (1 - alfa)*(beta*gama*pow(x_til[j + jj + JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[j + jj + JJ + delta_x], 2))/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB), 2);		//q
						constante += (1 - alfa)*(beta*(1-gama)*pow(x_til[j + jj + 2*JJ + delta_x], 2) + (1-beta)*gama*pow(x_hat[j + jj + 2*JJ + delta_x], 2));																		//z
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				if ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)
				{
					if (sistema_a->GetFlagVfol() == true)
						for (int r = 0; r < R; r++)
						{
							vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta].set(GRB_DoubleAttr_Obj, (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR - 2*(1-alfa)*(beta*gama*x_til[r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta_x])/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta].get(GRB_DoubleAttr_UB), 2));		//vfol
							constante += (1 - alfa)*(beta*gama*pow(x_til[r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta].get(GRB_DoubleAttr_UB), 2);		//vfol
						}
					delta += nt + flag2*R;
				}
				else
					delta += nt;
			}
		}
		modelosGRB[modelo]->update();

		// Termos quadráticos
		delta = 0;
		GRBQuadExpr fo;
		//x = [pt u cp (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		double * coeffs;
		coeffs = new double[numero_var_subp[modelo]];
		for (int i = 0; i < numero_var_subp[modelo]; i++)
			coeffs[i] = vars[modelo][i].get(GRB_DoubleAttr_Obj);
		fo.addTerms(coeffs, vars[modelo], numero_var_subp[modelo]);
		//for (int i = 0; i < numero_var_subp[modelo]; i++)		// nao precisaria pois o vetor coeffs é todo completo, a n ser q alguma variavel n seja considerada no termo proximal
		//	coeffs[i] = 0;
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal e futura
		{
			cen = 0;
			if ( vtr_a[vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				deltaT = sistema_a->GetDeltaT1();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
						coeffs[i + 3*I + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + 3*I + delta].get(GRB_DoubleAttr_UB),2);		//F^2
					if (sistema_a->GetFlagInitAproxCT() == 0)
						coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT*alfa/optimal_LR + (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB),2);		//pt^2
					else
						coeffs[i + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB),2);		//pt^2
					
					coeffs[i + I + delta] = (1 - alfa)*(beta*(1-gama) + (1-beta)*gama);		//u^2
					coeffs[i + 2*I + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB),2);		//cp^2
				}
				if (flag1 == 1)		// considerar rede
				{
					for (size_t b = 0; b < B - 1; b++)
						coeffs[b + I*(3 + flag4) + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_LB),2);		//teta
					for (size_t b = 0; b < B; b++)
					{
						if ( vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) != 0 )
							coeffs[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB),2);		//def
						else
							coeffs[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta] = 0;		//def
					}
				}
				else
					coeffs[I*(3 + flag4) + 5*R + 3*JJ + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB),2);		//def

				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].get(GRB_DoubleAttr_UB),2);				//ph
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_LB),2);	//v
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB),2);	//d
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB),2);	//s
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB),2);	//phmax
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						coeffs[j + jj + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB),2);					//phg
						coeffs[j + jj + JJ + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB),2);		//q
						coeffs[j + jj + 2*JJ + delta] = (1 - alfa)*(beta*(1-gama) + (1-beta)*gama);		//z
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				delta += nt;
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
						coeffs[i + 3*I + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + 3*I + delta].get(GRB_DoubleAttr_UB),2);		//F^2
					if (sistema_a->GetFlagInitAproxCT() == 0)
						coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR + (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB),2);		//pt^2
					else
						coeffs[i + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB),2);		//pt^2
					
					coeffs[i + I + delta] = (1 - alfa)*(beta*(1-gama) + (1-beta)*gama);		//u^2
					coeffs[i + 2*I + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB),2);		//cp^2
				}
				if (flag1 == 1)		// considerar rede
				{
					for (size_t b = 0; b < B - 1; b++)
						coeffs[b + I*(3 + flag4) + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4) + delta].get(GRB_DoubleAttr_LB),2);		//teta
					for (size_t b = 0; b < B; b++)
					{
						if ( vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) != 0 )
							coeffs[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB),2);		//def
						else
							coeffs[b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta] = 0;		//def
					}
				}
				else
					coeffs[I*(3 + flag4) + 5*R + 3*JJ + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][I*(3 + flag4) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_UB),2);		//def
				
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + delta].get(GRB_DoubleAttr_UB),2);				//ph
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + R + delta].get(GRB_DoubleAttr_LB),2);	//v
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB),2);	//d
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB),2);	//s
					coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB),2);	//phmax
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						coeffs[j + jj + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB),2);					//phg
						coeffs[j + jj + JJ + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB),2);		//q
						coeffs[j + jj + 2*JJ + delta] = (1 - alfa)*(beta*(1-gama) + (1-beta)*gama);		//z
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}

				if ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)
				{
					if (sistema_a->GetFlagVfol() == true)
						for (int r = 0; r < R; r++)
							coeffs[r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta] = (1 - alfa)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta].get(GRB_DoubleAttr_UB),2);		//vfol
					delta += nt + flag2*R;
				}
				else
					delta += nt;
			}
		}
		fo.addTerms(coeffs, vars[modelo], vars[modelo], numero_var_subp[modelo]);
		modelosGRB[modelo]->setObjective(fo, GRB_MINIMIZE);

		// Termos constantes (devem ser adicionados após o setObjective com GRBQuadExpr, pois setObjective() replaces the entire existing objective)
		//int n_bin_janela_pri = vtr_nos_pri[modelo].size() * (I + JJ);
		//int n_con_janela_pri = numero_var_subp[modelo] - n_bin_janela_pri;
		modelosGRB[modelo]->set(GRB_DoubleAttr_ObjCon, constante);
		modelosGRB[modelo]->update();
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
void Heuristicas::CriarFuncaoObjetivoL(int modelo)
{
	// Modela a f.o. do subproblema com termo proximal linear
	// Só os parametros alfa e beta sao usados.
	try 
	{
		int R = sistema_a->hidreletricasVtr.size();
		int I = sistema_a->termeletricasVtr.size();
		int B = sistema_a->barrasVtr.size();
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
		int cen, jj;
		int delta = 0;
		double deltaT;
		double constante = 0;

		double optimal_LR = CalcularFuncaoObjetivoRL(modelo);
		// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
		vetorint vtr_a;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
			vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
		for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
			vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

		double const * x_aux;
		if ( beta == 1)
			x_aux = &x_til[0];
		else
			x_aux = &x_hat[0];

		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal e futura
		{
			cen = 0;
			if ( vtr_a[vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				deltaT = sistema_a->GetDeltaT1();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_Obj, alfa*1*deltaT/optimal_LR);						//F
						if (x_aux[i + I + delta] == 0)
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, 1*(1-alfa));		//u
						else
						{
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, -1*(1-alfa));		//u
							constante++;
						}
					}
					else
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, alfa*sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT/optimal_LR);			//pt
						if (x_aux[i + I + delta] == 0)
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, 1*(1-alfa) + alfa*sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT/optimal_LR);		//u
						else
						{
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, -1*(1-alfa) + alfa*sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT/optimal_LR);		//u
							constante++;
						}
					}
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, alfa*1/optimal_LR);		//cp
				}
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						if (x_aux[j + jj + 2*JJ + delta] == 0)
							vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, 1*(1-alfa));		//z
						else
						{
							vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, -1*(1-alfa));		//z
							constante++;
						}
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				if (flag1 == 1)
					for (int b = 0; b < B; b++)
						vars[modelo][b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, alfa*sistema_a->GetCustoDeficit()*deltaT/optimal_LR);		//def
				else
					vars[modelo][I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, alfa*sistema_a->GetCustoDeficit()*deltaT/optimal_LR);		//def

				delta += nt;
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_Obj, alfa*1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//F
						if (x_aux[i + I + delta] == 0)
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, 1*(1-alfa));		//u
						else
						{
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, -1*(1-alfa));		//u
							constante++;
						}
					}
					else
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, alfa*sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//pt
						if (x_aux[i + I + delta] == 0)
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, 1*(1 - alfa) + alfa*sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//u
						else
						{
							vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, -1*(1 - alfa) + alfa*sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//u
							constante++;
						}
					}
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, alfa* 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//cp
				}
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						if (x_aux[j + jj + 2*JJ + delta] == 0)
							vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, 1*(1 - alfa));		//z
						else
						{
							vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, -1*(1 - alfa));		//z
							constante++;
						}
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				if (flag1 == 1)
					for (int b = 0; b < B; b++)
						vars[modelo][b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, alfa*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//def
				else
					vars[modelo][I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, alfa*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//def
			
				if ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)		// && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) > 0))
				{
					if (sistema_a->GetFlagVfol() == true)
						for (int r = 0; r < R; r++)
							vars[modelo][r + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta].set(GRB_DoubleAttr_Obj, alfa*(20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)/optimal_LR);		//vfol
					delta += nt + flag2*R;
				}
				else
					delta += nt;
			}
		}
		modelosGRB[modelo]->update();

		// Termos quadráticos (quando existirem na funçao objetivo)
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			GRBQuadExpr fo;
			//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
			double * coeffs;
			coeffs = new double[numero_var_subp[modelo]];
			for (int i = 0; i < numero_var_subp[modelo]; i++)
				coeffs[i] = vars[modelo][i].get(GRB_DoubleAttr_Obj);
			fo.addTerms(coeffs, vars[modelo], numero_var_subp[modelo]);
			delete coeffs;
			double coeficiente;
			GRBVar variavel;

			delta = 0;
			for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal e futura
			{
				cen = 0;
				if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
					cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());

				if ((0 <= vtr_nos_pri[modelo][vtr_i]) && (vtr_nos_pri[modelo][vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
				{
					deltaT = sistema_a->GetDeltaT1();
					for (int i = 0; i < I; i++)
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT*alfa/optimal_LR;		//pt^2
						variavel = vars[modelo][i + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel, 1);
					}
					delta += nt;
				}
				else	// nós do estágio 2
				{
					deltaT = sistema_a->GetDeltaT2();
					for (int i = 0; i < I; i++)
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*alfa/optimal_LR;		//pt^2
						variavel = vars[modelo][i + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel, 1);
					}
			
					if ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)		// && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) > 0))
						delta += nt + flag2*R;
					else
						delta += nt;
				}
			}
			modelosGRB[modelo]->setObjective(fo, GRB_MINIMIZE);
			modelosGRB[modelo]->update();
		}

		// Constantes
		modelosGRB[modelo]->set(GRB_DoubleAttr_ObjCon, constante);
		modelosGRB[modelo]->update();
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
void Heuristicas::AtualizarVariaveis(int janela)
{
	try 
	{
		vetorint * vtr_a;
		vetorint delta_vtr;
		if (janela == 1)
		{
			vtr_a = &vtr_nos_pri[0];
			for (int modelo = 0; modelo < n_modelos; modelo++)
				delta_vtr.push_back(0);
		}
		if (janela == 2)
		{
			vtr_a = &vtr_nos_fut[0];
			for (int modelo = 0; modelo < n_modelos; modelo++)
				delta_vtr.push_back(numero_var_subp_principal[modelo]);
		}

		int R = sistema_a->hidreletricasVtr.size();
		int I = sistema_a->termeletricasVtr.size();
		int B = sistema_a->barrasVtr.size();
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
		int cen, jj;
		int delta;

		for (int modelo = 0; modelo < n_modelos; modelo++)
		{
			delta = delta_vtr[modelo];
			for (size_t vtr_i = 0; vtr_i < vtr_a[modelo].size(); vtr_i++)
			{
				cen = 0;
				if ( vtr_a[modelo][vtr_i] >= sistema_a->GetTt2() )
					cen = (vtr_a[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());

				for (int i = 0; i < I; i++)
				{
					vars[modelo][i + I + delta].set(GRB_DoubleAttr_LB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + I + i]);		//u
					vars[modelo][i + I + delta].set(GRB_DoubleAttr_UB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + I + i]);
				}
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R;
				for (int r = 0; r < R; r++)
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_LB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + j + jj + 2*JJ]);	//z
						vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_UB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + j + jj + 2*JJ]);
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}

				if ( (vtr_a[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_a[modelo][vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) )
					delta += nt + flag2*R;
				else
					delta += nt;
			}

			modelosGRB[modelo]->update();
		}
		//delete[] vtr_a;
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
void Heuristicas::AtualizarRestricoes(int modelo, CMatrizEsparsa * x_spr)
{
	// para t > 0 (modelo > 0) atualizar o RHS dos subproblemas que usam a soluçao anterior como dados de entrada!!
	// matriz A tem seis elementos, em que o índice inicial é dado pelo vetor posicao_rest_din, referente as restricoes de Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida
	// vetor posicao_rest_din é igual para todos os subproblemas, pois esses possuem o msm numero de restriçoes antes de adicionar a restriçao de volume meta!!
	CMatrizEsparsa M;
	int nt_rest;
	for (size_t i = 0; i < posicao_rest_din[modelo].size(); i++)	// loop para as 4 submatrizes q acoplam no tempo, Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida
	{
		M = A[i];
		M.MultiplicarPorMatriz(x_spr);		// o resultado sera os valores de entrada para o subproblema i, visto que a soluçao para os subp > i sao nulas ainda (desconhecidas)
		
		nt_rest = M.GetNlin() / T;		// = R (balanço hidrico), Tup, Tdown, I (rampa up), I (rampa down) e I (custo de partida)
		for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais
			for (int iii = 0; iii < nt_rest; iii++)
				constr[modelo][iii + ii*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_pri[modelo][ii]*nt_rest + iii][1]) - double(M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0)));
				//constr[modelo][iii + ii*nt_rest + posicao_rest_din[i]].set(GRB_DoubleAttr_RHS, constr[modelo][iii + ii*nt_rest + posicao_rest_din[i]].get(GRB_DoubleAttr_RHS) - M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0));
		for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
			for (int iii = 0; iii < nt_rest; iii++)
				constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_fut[modelo][ii]*nt_rest + iii][1]) - double(M.GetElemento(vtr_nos_fut[modelo][ii]*nt_rest + iii, 0)));
				//constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[i]].set(GRB_DoubleAttr_RHS, constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[i]].get(GRB_DoubleAttr_RHS) - M.GetElemento(vtr_nos_fut[modelo][ii]*nt_rest + iii, 0));
	}
	modelosGRB[modelo]->update();	// Atualiza o modelo Gurobi
}
double Heuristicas::CalcularFuncaoObjetivo(int modelo, CMatrizEsparsa * x_spr)
{
	// Calcular funçao objetivo (janelas principais) com a solução x_spr
	double funcao_obj = 0;
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen;
	double deltaT;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int delta_x;

	// Variáveis das janelas principais
	//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
	{
		cen = 0;
		if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
			cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		delta_x = nt*vtr_nos_pri[modelo][vtr_i] + flag2*cen*R;

		if ((0 <= vtr_nos_pri[modelo][vtr_i]) && (vtr_nos_pri[modelo][vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
		{
			deltaT = sistema_a->GetDeltaT1();
			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
					funcao_obj += x_spr->GetElemento(i + 3*I + delta_x, 0) * 1*deltaT;		//F
					//funcao_obj += x[i + 3*I + delta_x] * 1*deltaT;		//F
				else
				{
					funcao_obj += x_spr->GetElemento(i + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
					funcao_obj += x_spr->GetElemento(i + I + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;		//u
				}
				funcao_obj += x_spr->GetElemento(i + 2*I + delta_x, 0) * 1;		//cp
			}
			if (flag1 == 1)
				for (int b = 0; b < B; b++)
					funcao_obj += x_spr->GetElemento(b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT;		//def
			else
				funcao_obj += x_spr->GetElemento(I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT;		//def
		}
		else	// nós do estágio 2
		{
			deltaT = sistema_a->GetDeltaT2();
			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
					funcao_obj += x_spr->GetElemento(i + 3*I + delta_x, 0) * 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//F
				else
				{
					funcao_obj += x_spr->GetElemento(i + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);			//pt
					funcao_obj += x_spr->GetElemento(i + I + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
				}
				funcao_obj += x_spr->GetElemento(i + 2*I + delta_x, 0) * 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//cp
			}
			if (flag1 == 1)
				for (int b = 0; b < B; b++)
					funcao_obj += x_spr->GetElemento(b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
			else
				funcao_obj += x_spr->GetElemento(I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
			
			if ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)		// && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) > 0))
				if (sistema_a->GetFlagVfol() == true)
					for (int r = 0; r < R; r++)
						funcao_obj += x_spr->GetElemento(r + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta_x, 0) * (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
		}
	}

	// Termos quadráticos
	if (sistema_a->GetFlagInitAproxCT() == 0)
	{
		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]

		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
		{
			cen = 0;
			if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			delta_x = nt*vtr_nos_pri[modelo][vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

			if ((0 <= vtr_nos_pri[modelo][vtr_i]) && (vtr_nos_pri[modelo][vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				deltaT = sistema_a->GetDeltaT1();
				for (int i = 0; i < I; i++)
					funcao_obj += pow(x_spr->GetElemento(i + delta_x, 0), 2) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;			//pt^2
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2();
				for (int i = 0; i < I; i++)
					funcao_obj += pow(x_spr->GetElemento(i + delta_x, 0), 2) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);			//pt^2
			}
		}
	}
	return funcao_obj;
}
double Heuristicas::CalcularFuncaoObjetivoRL(int modelo)
{
	// Calcular funçao objetivo (janelas principais e futuras) a partir da soluçao da RL x_hat
	double funcao_obj = 0;
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen;
	double deltaT;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int delta_x;
	// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
	vetorint vtr_a;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
		vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
	for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
		vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

	//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
	for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal e futura
	{
		cen = 0;
		if ( vtr_a[vtr_i] >= sistema_a->GetTt2() )
			cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		delta_x = nt*vtr_a[vtr_i] + flag2*cen*R;

		if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
		{
			deltaT = sistema_a->GetDeltaT1();
			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
					funcao_obj += x_hat[i + 3*I + delta_x] * 1*deltaT;		//F
					//funcao_obj += x[i + 3*I + delta_x] * 1*deltaT;		//F
				else
				{
					funcao_obj += x_hat[i + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
					funcao_obj += x_hat[i + I + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;		//u
				}
				funcao_obj += x_hat[i + 2*I + delta_x] * 1;		//cp
			}
			if (flag1 == 1)
				for (int b = 0; b < B; b++)
					funcao_obj += x_hat[b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT;		//def
			else
				funcao_obj += x_hat[I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT;		//def
		}
		else	// nós do estágio 2
		{
			deltaT = sistema_a->GetDeltaT2();
			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
					funcao_obj += x_hat[i + 3*I + delta_x] * 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//F
				else
				{
					funcao_obj += x_hat[i + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);			//pt
					funcao_obj += x_hat[i + I + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
				}
				funcao_obj += x_hat[i + 2*I + delta_x] * 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//cp
			}
			if (flag1 == 1)
				for (int b = 0; b < B; b++)
					funcao_obj += x_hat[b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
			else
				funcao_obj += x_hat[I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
			
			if ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)		// && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1()) > 0))
				if (sistema_a->GetFlagVfol() == true)
					for (int r = 0; r < R; r++)
						funcao_obj += x_hat[r + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + flag1*B + (1 - flag1) + delta_x] * (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
		}
	}

	// Termos quadráticos
	if (sistema_a->GetFlagInitAproxCT() == 0)
	{
		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]

		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal
		{
			cen = 0;
			if ( vtr_a[vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			delta_x = nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				deltaT = sistema_a->GetDeltaT1();
				for (int i = 0; i < I; i++)
					funcao_obj += pow(x_hat[i + delta_x], 2) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;			//pt^2
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2();
				for (int i = 0; i < I; i++)
					funcao_obj += pow(x_hat[i + delta_x], 2) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);			//pt^2
			}
		}
	}
	return funcao_obj;
}
void Heuristicas::AlocarSolucao(int modelo, CMatrizEsparsa * x_spr)
{
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen;
	int R = sistema_a->hidreletricasVtr.size();
	int delta = 0;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
	{
		if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
			cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		else
			cen = 0;
		if ( ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) > 0) )
		{
			if (sistema_a->GetFlagVfol() == true)
			{
				for (int ii = 0; ii < nt + R; ii++)
					x_spr->InserirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(vars[modelo][ii + delta].get(GRB_DoubleAttr_X)));
				delta += nt + R;
			}
		}
		else
		{
			for (int ii = 0; ii < nt; ii++)
				x_spr->InserirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(vars[modelo][ii + delta].get(GRB_DoubleAttr_X)));
			delta += nt;
		}
	}
}
void Heuristicas::AlterarEstadoVarBin(int modelo, int posicao, int estado)
{
	vars[modelo][posicao].set(GRB_DoubleAttr_LB, estado);
	vars[modelo][posicao].set(GRB_DoubleAttr_UB, estado);
}
void Heuristicas::ResolverInviabilidadeHidro(int modelo)
{
	// Identificar qual hidro esta causando a inviabilidade na restriçao de balanço hídrico e desligar tds suas unidades da usina
	modelosGRB[modelo]->computeIIS();
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int jj;
	int delta = 0;
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais
	{
		for (int iii = 0; iii < R; iii++)
		{
			if ( constr[modelo][iii + ii*R + posicao_rest_din[modelo][0]].get(GRB_IntAttr_IISConstr) == 1 )
			{
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R + 2*JJ;
				for (int rr = 0; rr < iii; rr++)
					jj += sistema_a->hidreletricasVtr[rr].GetNGrupos();	// coloca jj apontando para a variavel z da usina iii

				for (int j = 0; j < sistema_a->hidreletricasVtr[iii].GetNGrupos(); j++)
					AlterarEstadoVarBin(modelo, delta + jj + j, 0);		// posiçao da variavel z (referente a usina iii e ao nó ii do subproblema)
			}
		}
		if ( (vtr_nos_pri[modelo][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][ii] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) )
			delta += nt + flag2*R;
		else
			delta += nt;
	}
	for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
	{
		for (int iii = 0; iii < R; iii++)
		{
			if ( constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*R + posicao_rest_din[modelo][0]].get(GRB_IntAttr_IISConstr) == 1 )
			{
				jj = I*(3 + flag4) + flag1*(B - 1) + 5*R + 2*JJ;
				for (int rr = 0; rr < iii; rr++)
					jj += sistema_a->hidreletricasVtr[rr].GetNGrupos();	// coloca jj apontando para a variavel z da usina iii

				for (int j = 0; j < sistema_a->hidreletricasVtr[iii].GetNGrupos(); j++)
					AlterarEstadoVarBin(modelo, delta + jj + j, 0);
			}
		}
		if ( (vtr_nos_fut[modelo][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[modelo][ii] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) )
			delta += nt + flag2*R;
		else
			delta += nt;
	}
	modelosGRB[modelo]->update();
}

void Heuristicas::ConferirDeficit(int modelo)
{
	vetorfloat def_por_no;
	def_por_no.resize(vtr_nos_pri[modelo].size() + vtr_nos_fut[modelo].size());

	// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
	vetorint vtr_a;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
		vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
	for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
		vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int delta = 0;
	for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)
	{
		def_por_no[vtr_i] = 0;

		if (flag1 == 1)		// considerar rede
		{
			for (size_t b = 0; b < B; b++)
				def_por_no[vtr_i] += vars[modelo][b + I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_X);		//def
		}
		else
			def_por_no[vtr_i] += vars[modelo][I*(3 + flag4) + flag1*(B - 1) + 5*R + 3*JJ + delta].get(GRB_DoubleAttr_X);		//def
		
		// conferir se def != 0, se for ligar termicas nesse nó até o final e dar um break!!
		if ( def_por_no[vtr_i] != 0 )
		{
			vetorint prior_list = ListaPrioridades(modelo, vtr_a[vtr_i], 1);
			for (size_t i = 0; i < prior_list.size(); i++)
			{
				if ( vars[modelo][prior_list[i] + I + delta].get(GRB_DoubleAttr_X) == 0 )		// conferir as usinas na ordem da lista...
				{	
					// ligar usina até o fim do horizonte
					int delta_a = 0;
					for (size_t vtr_ii = vtr_i; vtr_ii < vtr_a.size(); vtr_ii++)
					{
						AlterarEstadoVarBin(modelo, prior_list[i] + I + delta + delta_a, 1);
						if ( ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && (vtr_a[vtr_i] + 1 >= sistema_a->GetTt2()) )
							delta_a += nt + flag2*R;
						else
							delta_a += nt;
					}
					def_por_no[vtr_i] -= sistema_a->termeletricasVtr[prior_list[i]].GetPmin();
					if ( def_por_no[vtr_i] <= 0 )
						break;
				}
			}
			break;
		}


		if ( ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && (vtr_a[vtr_i] + 1 >= sistema_a->GetTt2()) )
			delta += nt + flag2*R;
		else
			delta += nt;
	}


	// ligar térmicas com capacidade sum(Ptmin) suficiente para atender o déficit do primeiro nó
	// manter as térmicas ligadas até o final do horizonte do subproblema
	// liga e tenta resolver, se for inviável deixa quieto...
	//???problema é a q a soluçao do próximo problema com u fixo depende de x_hat...
	//fazer funçao como a do balanço hidraulico para se existir inviabilidade com relaçao ao min. up/down time ou à rampa -> alterar estado da térmica


}
void Heuristicas::ResolverInviabilidadeTermo(int modelo)
{
	// Identificar qual termo esta causando a inviabilidade na restriçao de min. up/down time e rampa e alterar estado das usinas
	
	//modelosGRB[modelo]->computeIIS();
	//int R = sistema_a->hidreletricasVtr.size();
	//int I = sistema_a->termeletricasVtr.size();
	//int B = sistema_a->barrasVtr.size();
	//int JJ = 0;
	//for (int r = 0; r < R; r++)
	//	JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	//int jj;
	//int delta = 0;
	//int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	//for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais
	//{
	//	for (int iii = 0; iii < R; iii++)
	//	{
	//		if ( constr[modelo][iii + ii*R + posicao_rest_din[modelo][0]].get(GRB_IntAttr_IISConstr) == 1 )
	//		{
	//			jj = I*(3 + flag4) + flag1*(B - 1) + 5*R + 2*JJ;
	//			for (int rr = 0; rr < iii; rr++)
	//				jj += sistema_a->hidreletricasVtr[rr]->GetNGrupos();	// coloca jj apontando para a variavel z da usina iii

	//			for (int j = 0; j < sistema_a->hidreletricasVtr[iii]->GetNGrupos(); j++)
	//				AlterarEstadoVarBin(modelo, delta + jj + j, 0);		// posiçao da variavel z (referente a usina iii e ao nó ii do subproblema)
	//		}
	//	}
	//	if ( (vtr_nos_pri[modelo][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][ii] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) )
	//		delta += nt + flag2*R;
	//	else
	//		delta += nt;
	//}
	//for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
	//{
	//	for (int iii = 0; iii < R; iii++)
	//	{
	//		if ( constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*R + posicao_rest_din[modelo][0]].get(GRB_IntAttr_IISConstr) == 1 )
	//		{
	//			jj = I*(3 + flag4) + flag1*(B - 1) + 5*R + 2*JJ;
	//			for (int rr = 0; rr < iii; rr++)
	//				jj += sistema_a->hidreletricasVtr[rr]->GetNGrupos();	// coloca jj apontando para a variavel z da usina iii

	//			for (int j = 0; j < sistema_a->hidreletricasVtr[iii]->GetNGrupos(); j++)
	//				AlterarEstadoVarBin(modelo, delta + jj + j, 0);
	//		}
	//	}
	//	if ( (vtr_nos_fut[modelo][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[modelo][ii] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) )
	//		delta += nt + flag2*R;
	//	else
	//		delta += nt;
	//}
	//modelosGRB[modelo]->update();
}
vetorint Heuristicas::ListaPrioridades(int modelo, int no, int tipo)
{
	// tipo = 1 termos; 2 hidros; - é a lista invertida;
	int cen = 0;
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int jj;
	if ( no >= sistema_a->GetTt2() )
		cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
	vetorint lista;
	vetorfloat2 lista2;
	lista2.resize(I);

	if ( tipo == 1 || tipo == -1 )		// termicas
	{
		for (size_t i = 0; i < I; i++)
		{
			lista2[i].push_back(x_til[i + I + nt*no + flag2*cen*R]);
			lista2[i].push_back(i);
		}
		sort(lista2.begin(), lista2.end());
	}
	else
	{
		for (size_t r = 0; r < R; r++)	// unidades hidrelétricas
		{
			jj = I*(3 + flag4) + flag1*(B - 1) + 5*R + 2*JJ;
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				lista2[j + jj].push_back(x_til[j + jj + nt*no + flag2*cen*R]);
				lista2[j + jj].push_back(j + jj);
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		sort(lista2.begin(), lista2.end());		// ordena em ordem crescente de valores elemento [0] é a unidade que tem o menor z
	}

	for (size_t i = 0; i < lista2.size(); i++)
		lista.push_back(int(lista2[i][1]));
	if ( tipo > 0 )
		reverse(lista.begin(), lista.end());			// retorna o índice das unidades na ordem que devem ser ligadas
		
	return lista;										// retorna o índice das unidades na ordem que devem ser desligadas
}

double Heuristicas::ResolverHeuristica()
{
	//// Receber valores de x_hat e x_til!!!
	//x_hat = resultadosGurobi->GetX_med(0);
	//x_til = resultadosGurobi->GetX_med(1);

	// Imprimir x_hat e x_til
	// Receber valores de x_hat e x_til!!!
	x_hat = resultadosGurobi->GetX_med(0);
	//resultadosGurobi->ExportarXmed("x_hat.txt");
	x_til = resultadosGurobi->GetX_med(1);
	//resultadosGurobi->ExportarXmed("x_til.txt");

	// Atualizar valores limites das variáveis fixas (quando algumas delas forem fixas)
	if ( (flagH1 == 0) || (flagH2 == 0))
	{
		if ( flagH1 == 0 )	// Se as var. binárias forem fixas elas devem ser atualizadas antes da resolução (janela principal)
			AtualizarVariaveis(1);
		if ( flagH2 == 0 )	// Se as var. binárias forem fixas elas devem ser atualizadas antes da resolução (janela futuro)
			AtualizarVariaveis(2);
	}

	CMatrizEsparsa * X_it;
	X_it = new CMatrizEsparsa(n, 1);	// x_spr contem a soluçao de tds os subproblemas anteriores subp <= i
	double fo_it = 0;
	fo_subp.resize(n_modelos);		// somente para conferir resultados!
	nStatus.clear();
	nStatus.resize(n_modelos);
	Status = 0;

	// Resolver subproblemas
	ResolverSubps(X_it, &fo_it);
	// Gravar melhor soluçao durante as iteraçoes da RL
	if ( (Status == 0) && (fo_it <= fo) )	// Só grava soluçao se todos os subproblemas forem resolvidos
	{
		X = X_it;
		fo = fo_it;
	}
	X_it = NULL;
	return fo;

	//if (Status == 0)		// Só grava soluçao se todos os subproblemas forem resolvidos
	//{
	//	if ( gravar_x == true )			// gravar soluçao em x
	//		for (int i = 0; i < n; i++)
	//			x[i] = x_spr->GetElemento(i, 0);

	//	//// imprimir para conferir soluçao da Heuristica
	//	//vetorfloat lixo;
	//	//resultadosGurobi->CarregarResultadosED(fo, x, lixo, Status);
	//	//resultadosGurobi->EscreverArquivo("heuristica.out", 0, 0);
	//}
}
double Heuristicas::ResolverHeuristica(vetorfloat x_a, vetorfloat x_a2)
{
	// Receber valores de x_hat e x_til!!!
	x_hat = x_a;
	x_til = x_a2;

	// Atualizar valores limites das variáveis fixas (quando algumas delas forem fixas)
	if ( (flagH1 == 0) || (flagH2 == 0))
	{
		if ( flagH1 == 0 )	// Se as var. binárias forem fixas elas devem ser atualizadas antes da resolução (janela principal)
			AtualizarVariaveis(1);
		if ( flagH2 == 0 )	// Se as var. binárias forem fixas elas devem ser atualizadas antes da resolução (janela futuro)
			AtualizarVariaveis(2);
	}

	CMatrizEsparsa * X_it;
	X_it = new CMatrizEsparsa(n, 1);	// x_spr contem a soluçao de tds os subproblemas anteriores subp <= i
	double fo_it = 0;
	fo_subp.resize(n_modelos);		// somente para conferir resultados!
	nStatus.clear();
	nStatus.resize(n_modelos);
	Status = 0;

	// Resolver subproblemas
	ResolverSubps(X_it, &fo_it);
	// Gravar melhor soluçao durante as iteraçoes da RL
	if ( (Status == 0) && (fo_it <= fo) )	// Só grava soluçao se todos os subproblemas forem resolvidos
	{
		X = X_it;
		fo = fo_it;
	}
	X_it = NULL;

	// imprimir x
	X->ImprimirMatrizArquivo();
	
	
	return fo;
}
void Heuristicas::ResolverSubps(CMatrizEsparsa * x_spr, double * fo_a)
{
	try
	{
		for (int i = 0; i < n_modelos; i++)
		{
			// Reset the model to an unsolved state, discarding any previously computed solution information.
			modelosGRB[i]->reset();

			// Criar funçao objetivo
			CriarFuncaoObjetivo(i);

			//modelosGRB[i]->write("subp.lp");

			// resolver subproblema i
			modelosGRB[i]->optimize();

			//nStatus[i] = modelosGRB[i]->get(GRB_IntAttr_Status);
			//if ((nStatus[i] == 2) || (nStatus[i] == 9))		// OPTIMAL or TIME_LIMIT
			if (modelosGRB[i]->get(GRB_IntAttr_SolCount))		// Number of solutions found during the most recent optimization, se > 0 fazer...
			{
				//for (int j = 0; j < numero_var_subp[i]; j++)
				//	cout << vars[i][j].get(GRB_DoubleAttr_X) << endl;
				//cout << "f.o.: " << modelosGRB[i]->get(GRB_DoubleAttr_ObjVal) << endl;

				// Alocar solucao do subp em x_spr (com a soluçao somente das janelas principais)
				AlocarSolucao(i, x_spr);
				
				//fo_subp[i] = double(modeloGRB.get(GRB_DoubleAttr_ObjVal));
				//fo_subp[i] = CalcularFuncaoObjetivo(i);
				//cout << double(modelosGRB[i]->get(GRB_DoubleAttr_ObjVal)) << endl;
				*fo_a += CalcularFuncaoObjetivo(i, x_spr);
				fo_subp[i] = *fo_a;
			}
			else
			{
				//modelosGRB[i]->computeIIS();				
				//modelosGRB[i]->write("Hrst_subprob.lp");	
				//modelosGRB[i]->write("Hrst_subprob.ilp");

				// n precisa inserir zeros em x_spr, ele já está iniciado com 0's. Basta n atribuir nenhum valor
				//fo_subp[i] = GRB_INFINITY;
				*fo_a += GRB_INFINITY;
				fo_subp[i] = *fo_a;
				Status++;
				break;
			}
			if ( i < n_modelos - 1)
				AtualizarRestricoes(i + 1, x_spr);		// Para subp > 0 antes de cada subproblema se resolvido o lado direito das restriçoes são atualizados!!
		}
		// parei aqui
		// esse caso esta dando erro!! descobrir o que é!!
		// as variaveis de folga (def e vfol) estao zeradas, acho que é isso que esta atrapalhando as soluçoes da Heuristica
		// testar heuristica com um unico subproblema fixo na soluçao da RL...
		// mesmo com as variáveis de folga liberadas a soluçao é estranha!!


	//cout << "Fim da heuristica" << endl;
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
}
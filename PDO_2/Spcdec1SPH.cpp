#include "Spcdec1SPH.h"

// Criar restrições
// ------------------------------------------------
// Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa Spcdec1SPH::MatrizBalHid()
{
	// T * R
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
		Addd.InserirMatriz(t * R, 0, (t + 1) * R - 1,R - 1, &Add, 0, 0);
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
		Addd2.InserirMatriz(t * R, 0, (t + 1) * R - 1,R - 1, &Add, 0, 0);
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
	CMatrizEsparsa Alin(R * T, n_a);
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
			Alin.InserirMatriz(l, c + R, l + R - 1,c + 2*R - 1, &Av, 0, 0);
			for (int tt = 0; tt <= t; tt++)
				Alin.InserirMatriz(l, tt*(n_a / T) + 2*R, l + R - 1, tt*(n_a / T) + 3*R - 1, &Ad, t*R, tt*R);
			if (t > 0)
				Alin.InserirMatriz(l, c - (n_a / T) + R, l + R - 1, c - (n_a / T) + 2*R - 1, &AvNeg, 0, 0);
		}
		else
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			Alin.InserirMatriz(l, c + R, l + R - 1, c + 2*R - 1, &Av, 0, 0);
			cc = 0;
			for (int tt = 0; tt <= t; tt++)
			{
				Alin.InserirMatriz(l, cc + 2*R, l + R - 1,cc + 3*R - 1, &Ad, t*R, tt*R);
				if (((tt + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((tt + 1 - sistema_a->GetTt1()) > 0))
					cc += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
				else
					cc += (n_a / T);
			}
			if (t > 0)
				if (t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) == sistema_a->GetTt1())
					Alin.InserirMatriz(l, (n_a / T)*sistema_a->GetTt1() - (n_a / T) + R, l + R - 1, (n_a / T)*sistema_a->GetTt1() - (n_a / T) + 2*R - 1, &AvNeg, 0, 0);
				else
					Alin.InserirMatriz(l, c - (n_a / T) + R, l + R - 1, c - (n_a / T) + 2*R - 1, &AvNeg, 0, 0);
		}
		l = l + R;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizVmeta()
{
	int n_a = n;
	size_t R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int c;
	c = R + (sistema_a->GetTt2() - 1)*(n_a / T);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + r, 1);
			if (sistema_a->GetFlagVfol() == true)
				a.InserirElemento(0, c + r + 3*R + 3*JJ, 1);
			Alin.JuntarColuna(&a);
		}
		c += (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T) + flag2*R;
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizFuncProd()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int vv = R;
	int jj = vv + 3*R + JJ;
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
				for (int napr = 0; napr < sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNappFPH(); napr++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, jj - JJ + j, 1);																	// phg
					a.InserirElemento(0, vv + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgV(napr));				// v
					a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgQ(napr));				// q
					a.InserirElemento(0, vv + R + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgD(napr));			// d
					if ((sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgV(napr) * sistema_a->hidreletricasVtr[r].GetVmin() < 0) || (sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgV(napr) * sistema_a->hidreletricasVtr[r].GetVmax() < 0))
						a.InserirElemento(0, jj + j + JJ, sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax());			// z -> adiciona um termo (1 - z)*Pmin na função, para phg <= pmin qdo z = 0; o rhs estava negativo, devido a aproximação linear, fazendo com q a hidro ficasse obrigatoriamente ligada
					if (sistema_a->hidreletricasVtr[r].GetInflueVert() == 1 )		// Se o s influencia já tenho contabilizado isso no d, caso contrário só subtrair o valor de s de d
					{}
					else		// vertimento n influencia no canal de fuga, entao subtraio o valor de s da defluencia
						a.InserirElemento(0, vv + 2*R + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgS(napr));	// s
					Alin.JuntarColuna(&a);
				}
			}
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		vv = vv + nt;
		jj = jj + nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizBalVazao()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = 2*R;
	int jj = dd + JJ + (2)*R;
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
CMatrizEsparsa Spcdec1SPH::MatrizBalPotencia()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = 0;
	int jj = dd + 4*R;
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
CMatrizEsparsa Spcdec1SPH::MatrizLimPhgMin()
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
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size(), l + JJ - 1, c + 4*sistema_a->hidreletricasVtr.size() + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size() + 2*JJ, l + JJ - 1, c + 4*sistema_a->hidreletricasVtr.size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizLimPhgMax()
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
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size(), l + JJ - 1,c + 4*sistema_a->hidreletricasVtr.size() + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size() + 2*JJ, l + JJ - 1,c + 4*sistema_a->hidreletricasVtr.size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizLimqMin()
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
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size() + JJ, l + JJ - 1, c + 4*sistema_a->hidreletricasVtr.size() + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size() + 2*JJ, l + JJ - 1, c + 4*sistema_a->hidreletricasVtr.size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizLimqMax()
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
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size() + JJ, l + JJ - 1, c + 4*sistema_a->hidreletricasVtr.size() + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + 4*sistema_a->hidreletricasVtr.size() + 2*JJ, l + JJ - 1, c + 4*sistema_a->hidreletricasVtr.size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec1SPH::MatrizReserva()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = 0;
	int jj = dd + 4*R + 2*JJ;
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
			a.InserirElemento(0, dd + r, - 1);
			//if (sistema_a->GetFlagPhmax() == 1)
			//	a.InserirElemento(0, dd + 4*R + r, 1);
			//else
			//{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				a.InserirElemento(0, jj + j, sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			//}
		}
		Alin.JuntarColuna(&a);
		dd += nt;
		jj += nt - JJ;
	}
	return Alin;
}

// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Spcdec1SPH::LimBalHid()
{
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
vetorfloat2 Spcdec1SPH::LimVmeta()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * sistema_a->GetNCenarios(), 2);
	IniciaMatriz(&Lim, 0);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (int r = 0; r < R; r++)
		{
			Lim[r + R*cen][1] = sistema_a->hidreletricasVtr[r].GetVMeta();
			Lim[r + R*cen][0] = 2;
		}
	}
	return Lim;
}
vetorfloat2 Spcdec1SPH::LimFuncProdL()
{
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	vetorfloat2 Lim;
	vetorfloat L;
	L.resize(2);
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				for (int napr = 0; napr < sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNappFPH(); napr++)
				{
					L[0] = 0;
					L[1] = sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax();
					Lim.push_back(L);
				}
	return Lim;
}
vetorfloat2 Spcdec1SPH::LimPhgMin()
{
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
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Spcdec1SPH::LimPhgMax()
{
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
vetorfloat2 Spcdec1SPH::LimQMin()
{
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
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Spcdec1SPH::LimQMax()
{
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
vetorfloat2 Spcdec1SPH::LimBalPotenciaL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Spcdec1SPH::LimBalVazaoL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Spcdec1SPH::LimReserva()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1 * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		Lim[t][0] = 2;
		Lim[t][1] = sistema_a->GetReserva(t);
	}
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void Spcdec1SPH::MatrizRestricoesLineares(CMatrizEsparsa &MM)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	M = MatrizLimPhgMin();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPhgMax();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalHid();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimqMin();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimqMax();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizVmeta();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalPotencia();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalVazao();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizFuncProd();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizReserva();
	MM.JuntarColuna(&M); M.ZerarMatriz();
}
void Spcdec1SPH::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor)
{
	vetorfloat2 L;
	L = LimPhgMin();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhgMax();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalHid();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMin();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMax();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimVmeta();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalPotenciaL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalVazaoL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimFuncProdL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimReserva();
	AlocarLimites(&L, LimTipo, LimValor);
}
// ------------------------------------------------

Spcdec1SPH::Spcdec1SPH(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;
	modelosGRB = new GRBModel(ambiente_gurobi);		// um modelo de otimização para todas hidros
	T = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );
	flag2 = int (sistema_a->GetFlagVfol());

	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	n = T * ( 4*sistema_a->hidreletricasVtr.size() + 3*JJ) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();

	// Criar variáveis
	CriarVariaveis();
	modelosGRB->update();	// Atualiza o modelo Gurobi.
	vars = modelosGRB->getVars();

	// Adicionar restrições
	CriarRestricoes();
	modelosGRB->update();	// Atualiza o modelo Gurobi.

	// Definir ajustes do solver
	//modelosGRB->write("probH.lp");	// Escreve modelo em arquivo

	// Nao escrever detalhes
	modelosGRB->getEnv().set(GRB_IntParam_OutputFlag, 0);
		
	//modelosGRB.getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
	//modelosGRB.getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

	modelosGRB->getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);	// Define o gap de tolerancia
	//modelosGRB->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
	//modelosGRB->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
	//modelosGRB->getEnv().set(GRB_IntParam_Threads, 1);
	//modelosGRB->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
	//modelosGRB.getEnv().set(GRB_IntParam_PreSparsify, 1);
	//modelosGRB.getEnv().set(GRB_IntParam_Method, 0);

	//modelosGRB.getEnv().set(GRB_IntParam_MIPFocus, 2);
	//modelosGRB.getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

	//modelosGRB->getEnv().set(GRB_DoubleParam_TimeLimit, 600);		// Limita o tempo de resolução do problema
	modelosGRB->getEnv().set(GRB_DoubleParam_TimeLimit, 5*sistema_a->GetNCenarios());

}
Spcdec1SPH::~Spcdec1SPH(void)
{
	delete modelosGRB;
	delete vars;
}

void Spcdec1SPH::CriarVariaveis()
{
	try 
	{
		// variáveis em cada problema somente para um período
		size_t R = sistema_a->hidreletricasVtr.size();
		// Problema de 1 período para cada usina
		// Estagio 1
		//x = [ph v d s phg q z]
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t r = 0; r < R; r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < R; r++)	//v
				modelosGRB->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t)), double(sistema_a->hidreletricasVtr[r].GetVmax(t)), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < R; r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				modelosGRB->addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < R; r++)	//s
			{
				modelosGRB->addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < R; r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modelosGRB->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modelosGRB->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
		}
		// Estagio 2
		//x = [ph v d s phg q z vfol]
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (int t = 0; t < sistema_a->GetTt2() - sistema_a->GetTt1(); t++)
			{
				for (size_t r = 0; r < R; r++)	//ph
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modelosGRB->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < R; r++)	//v
					modelosGRB->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t+sistema_a->GetTt1()+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), double(sistema_a->hidreletricasVtr[r].GetVmax(t+sistema_a->GetTt1()+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < R; r++)	//d
				{
					double qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
					modelosGRB->addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < R; r++)	//s
				{
					modelosGRB->addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < R; r++)	//phg
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						modelosGRB->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				for (size_t r = 0; r < R; r++)	//q
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						modelosGRB->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
				{
					for (size_t r = 0; r < R; r++)	//z
					{
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							modelosGRB->addVar(0, 1, 0.0, GRB_BINARY, "");
					}
				}
				else
				{
					for (size_t r = 0; r < R; r++)	//z
					{
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							modelosGRB->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					}
				}
			}
			if ( sistema_a->GetFlagVfol() == true )
				for (size_t r = 0; r < R; r++)	//vfol
					modelosGRB->addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
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
void Spcdec1SPH::CriarRestricoes()
{
	try
	{
		vetorint LimTipo;
		vetorfloat LimValor;
		LimTipo.resize(0);LimValor.resize(0);
		CMatrizEsparsa MM(0, n);

		MatrizRestricoesLineares(MM);
		MatrizLimitesLineares(&LimTipo, &LimValor);

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
				variavel = vars[MM.GetValorCol(l)];
				restricao.addTerms( &coeficiente, &variavel, 1);
				l = MM.GetValorLprox(l);
			}
			switch (LimTipo[i]) {
			case 0:
				modelosGRB->addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modelosGRB->addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modelosGRB->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
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
void Spcdec1SPH::CriarFuncaoObjetivoRL(const vetorfloat * const lambda)
{
	try 
	{
		//x = [ph v d s phg q z vfol]
		int flag2 = int (sistema_a->GetFlagVfol());

		// Termos lineares
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		int nt;
		int nt_dual = lambda->size() / T;
		int delta_dual = 0;
		double deltaT;
		int delta = 0;
		// vfolga
		int JJ = 0;
		for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
			JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
		delta = (4*sistema_a->hidreletricasVtr.size() + 3*JJ) * (sistema_a->GetTt2());
		deltaT = sistema_a->GetDeltaT2();
		if (sistema_a->GetFlagVfol() == true)
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					vars[delta + r].set(GRB_DoubleAttr_Obj, (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//vfol
				delta += (4*sistema_a->hidreletricasVtr.size() + 3*JJ) * (sistema_a->GetTt2() - sistema_a->GetTt1()) + flag2*sistema_a->hidreletricasVtr.size();
			}
		}
		//for (int t = 0; t < T; t++)
		//{
		//	if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//		nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		//	else
		//		nt = (n_a / T);
		//	if ((0 <= t) && (sistema_a->GetTt1() > t))
		//	{}
		//	else
		//	{
		//		cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		//		deltaT = sistema_a->GetDeltaT2();
		//		if (sistema_a->GetFlagVfol() == true)
		//		{
		//			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		//					vars[delta + r].set(GRB_DoubleAttr_Obj, (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//vfol
		//		}
		//	}
		//	delta = delta + nt;
		//}
		
		// Termos da RL
		//L = [Lpt Lph]
	
		// Lineares
		delta = 0;
		delta_dual = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				vars[r + delta].set(GRB_DoubleAttr_Obj, - lambda->at(r + sistema_a->termeletricasVtr.size() + delta_dual));			//ph
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
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
int Spcdec1SPH::ResolverProblemaRL(Spcdec1Results * resultadoGurobi, const vetorfloat * const lambda, const int iter)
{
	// Resolve o subproblema
	try
	{
		// Criar função objetivo, a seleçao dos lambdas de interesse é feita dentro da funçao abaixo
		CriarFuncaoObjetivoRL(lambda);
		modelosGRB->update();	// Atualiza o modelo Gurobi.

		// ajustes do solver

		modelosGRB->reset();

		// Otimizar
		modelosGRB->optimize();

		// Salvar resultados
		x.clear();
		x.resize(n);
		int nStatus = modelosGRB->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The H subproblem was not solved until optimality! " << nStatus << endl;
		if (modelosGRB->get(GRB_IntAttr_SolCount))	// se existir ao menos uma solução
		{
			fo = double(modelosGRB->get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n; i++)
				x[i] = double(vars[i].get(GRB_DoubleAttr_X));
		}
		else	// se não existir uma solução
		{
			for (int i = 0; i < n; i++)
				x[i] = 0;
			fo = 0;
		}
		resultadoGurobi->GravarSolucao(fo, x, nStatus, 0);

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->GravarSolucao(fo, x, e.getErrorCode(), 0);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}
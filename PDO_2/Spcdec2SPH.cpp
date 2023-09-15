#include "Spcdec2SPH.h"

// Criar restrições
// ------------------------------------------------
// Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa Spcdec2SPH::MatrizBalHid(int n_a, int ch)
{
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
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetUsinaJusante() != 0))
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
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetUsinaJusante() != 0))
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
					cc += (n_a / T) + flag2*cascata[ch].size();
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
			c += (n_a / T) + flag2*cascata[ch].size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizVmeta(int n_a, int ch)
{
	size_t R = cascata[ch].size();
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
	n_a = n_a - flag2*sistema_a->GetNCenarios()*R;
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*R);
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*R);
	int c;
	c = R + (sistema_a->GetTt2() - 1)*(n_a / T);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + r, 1);
			if (sistema_a->GetFlagVfol() == true)
				a.InserirElemento(0, c + r + (3+flag3)*R + 3*JJ, 1);
			Alin.JuntarColuna(&a);
		}
		c += (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T) + flag2*R;
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizLimPhgMin(int n_a, int ch)
{
	int JJ = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ,JJ);
	int cont = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmin());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ, n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();		// n_a externo n é alterado, pq é passado como cópia para essa função
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size(), l + JJ - 1, c + (4+flag3)*cascata[ch].size() + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size() + 2*JJ, l + JJ - 1, c + (4+flag3)*cascata[ch].size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*cascata[ch].size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizLimPhgMax(int n_a, int ch)
{
	int JJ = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ,JJ);
	int cont = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ, n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();		// n_a externo n é alterado, pq é passado como cópia para essa função
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size(), l + JJ - 1, c + (4+flag3)*cascata[ch].size() + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size() + 2*JJ, l + JJ - 1, c + (4+flag3)*cascata[ch].size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*cascata[ch].size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizLimqMin(int n_a, int ch)
{
	int JJ = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ,JJ);
	int cont = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmin());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ, n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size() + JJ, l + JJ - 1, c + (4+flag3)*cascata[ch].size() + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size() + 2*JJ, l + JJ - 1, c + (4+flag3)*cascata[ch].size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*cascata[ch].size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizLimqMax(int n_a, int ch)
{
	int JJ = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	CMatrizEsparsa A1(JJ);
	CMatrizEsparsa A2(JJ,JJ);
	int cont = 0;
	for (size_t r = 0; r < cascata[ch].size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
		{
			A2.InserirElemento(cont, cont, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
			cont++;
		}
	CMatrizEsparsa Alin(T * JJ, n_a);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
	int l = 0;
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size() + JJ, l + JJ - 1, c + (4+flag3)*cascata[ch].size() + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l, c + (4+flag3)*cascata[ch].size() + 2*JJ, l + JJ - 1, c + (4+flag3)*cascata[ch].size() + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*cascata[ch].size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizBalPotencia(int n_a, int ch)
{
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
	size_t R = cascata[ch].size();
	int nt;
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	int dd = 0;
	int jj = dd + (4+flag3)*R;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*cascata[ch].size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, -1);
			}
			Alin.JuntarColuna(&a);
			jj += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		}
		dd += nt;
		jj += nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizBalVazao(int n_a, int ch)
{
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
	size_t R = cascata[ch].size();
	int nt;
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	int dd = 2*R;
	int jj = dd + JJ + (2+flag3)*R;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*cascata[ch].size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			a.InserirElemento(0, dd + R + r, - 1);
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, -1);
			}
			Alin.JuntarColuna(&a);
			jj += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		}
		dd += nt;
		jj += nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizFuncProd(int n_a, int ch)
{
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
	size_t R = cascata[ch].size();
	int nt;
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	int vv = R;
	int jj = vv + (3+flag3)*R + JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*cascata[ch].size();
			else
				nt = (n_a / T);
		for (size_t r = 0; r < R ; r++)
		{	
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			{
				for (int napr = 0; napr < sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNappFPH(); napr++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, jj - JJ + j, 1);																					// phg
					a.InserirElemento(0, vv + r, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhgV(napr));				// v
					a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhgQ(napr));				// q
					a.InserirElemento(0, vv + R + r, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhgD(napr));			// d
					if ((sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhgV(napr) * sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin() < 0) || (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhgV(napr) * sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax() < 0))
						a.InserirElemento(0, jj + j + JJ, sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax());			// z -> adiciona um termo (1 - z)*Pmin na função, para phg <= pmin qdo z = 0; o rhs estava negativo, devido a aproximação linear, fazendo com q a hidro ficasse obrigatoriamente ligada
					if (sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetInflueVert() == 1 )		// Se o s influencia já tenho contabilizado isso no d, caso contrário só subtrair o valor de s de d
					{}
					else		// vertimento n influencia no canal de fuga, entao subtraio o valor de s da defluencia
						a.InserirElemento(0, vv + 2*R + r, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhgS(napr));	// s
					Alin.JuntarColuna(&a);
				}
			}
			jj += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		}
		vv += nt;
		jj += nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Spcdec2SPH::MatrizPhMax(int n_a, int ch)
{
	n_a = n_a - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
	size_t R = cascata[ch].size();
	int nt;
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*cascata[ch].size());
	int dd = 4*R;
	int jj = dd + R + 2*JJ;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*cascata[ch].size();
			else
				nt = (n_a / T);
			a.InserirElemento(0, dd + r, 1);
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			{
				a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
			}
			Alin.JuntarColuna(&a);
			jj += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		}
		dd += nt;
		jj += nt - JJ;
	}
	return Alin;
}

// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Spcdec2SPH::LimBalHid(int n_a, int ch)
{
	double deltaT;
	int cenario;
	int periodo;
	int R = cascata[ch].size();
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
vetorfloat2 Spcdec2SPH::LimVmeta(int n_a, int ch)
{
	int R = cascata[ch].size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * sistema_a->GetNCenarios(), 2);
	IniciaMatriz(&Lim, 0);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (int r = 0; r < R; r++)
		{
			Lim[r + R*cen][1] = sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVMeta();
			Lim[r + R*cen][0] = 2;
		}
	}
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimPhgMin(int n_a, int ch)
{
	int R = cascata[ch].size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimPhgMax(int n_a, int ch)
{
	int R = cascata[ch].size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimQMin(int n_a, int ch)
{
	int R = cascata[ch].size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	int jj = 0;
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimQMax(int n_a, int ch)
{
	int R = cascata[ch].size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			JJ++;
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, JJ * T, 2);
	IniciaMatriz(&Lim, 0);
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimBalPotenciaL(int n_a, int ch)
{
	int R = cascata[ch].size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimBalVazaoL(int n_a, int ch)
{
	int R = cascata[ch].size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimFuncProdL(int n_a, int ch)
{
	int R = cascata[ch].size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
	vetorfloat2 Lim;
	vetorfloat L;
	L.resize(2);
	for (int t = 0; t < T; t++)
		for (int r = 0; r < R; r++)
			for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
			for (int napr = 0; napr < sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNappFPH(); napr++)
				{
					L[0] = 0;
					L[1] = sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax();
					Lim.push_back(L);
				}
	return Lim;
}
vetorfloat2 Spcdec2SPH::LimPhMax(int n_a, int ch)
{
	int R = cascata[ch].size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void Spcdec2SPH::MatrizRestricoesLineares(int n_a, CMatrizEsparsa &MM, int ch)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	M = MatrizBalHid(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizVmeta(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPhgMin(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPhgMax(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimqMin(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimqMax(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalPotencia(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalVazao(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizFuncProd(n_a, ch);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if (sistema_a->GetFlagPhmax() == 1)
	{
		M = MatrizPhMax(n_a, ch);
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
}
void Spcdec2SPH::MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int ch)
{
	vetorfloat2 L;
	L = LimBalHid(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimVmeta(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhgMin(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhgMax(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMin(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMax(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalPotenciaL(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalVazaoL(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimFuncProdL(n_a, ch);
	AlocarLimites(&L, LimTipo, LimValor);
	if (sistema_a->GetFlagPhmax() == 1)
	{
		L = LimPhMax(n_a, ch);
		AlocarLimites(&L, LimTipo, LimValor);
	}
}
// ------------------------------------------------

Spcdec2SPH::Spcdec2SPH(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;
	flag3 = int (sistema_a->GetFlagPhmax());
	flag2 = int (sistema_a->GetFlagVfol());
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
	T = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );

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
	{
		n[ch] = T * (4+flag3) * cascata[ch].size() + flag2 * sistema_a->GetNCenarios() * cascata[ch].size();
		for (int r = 0; r < cascata[ch].size(); r++)
			n[ch] += T * (3*sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos());
	}

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
		//modelosGRB[ch]->write("probH.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modelosGRB[ch]->getEnv().set(GRB_IntParam_OutputFlag, 0);
		
		//modelosGRB[ch].getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[ch].getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		modelosGRB[ch]->getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);	// Define o gap de tolerancia
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_Threads, 1);
		//modelosGRB[ch]->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[ch].getEnv().set(GRB_IntParam_PreSparsify, 1);
		//modelosGRB[ch].getEnv().set(GRB_IntParam_Method, 0);

		//modelosGRB[ch].getEnv().set(GRB_IntParam_MIPFocus, 2);
		//modelosGRB[ch].getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

		//modelosGRB[ch]->getEnv().set(GRB_DoubleParam_TimeLimit, 600);		// Limita o tempo de resolução do problema
		modelosGRB[ch]->getEnv().set(GRB_DoubleParam_TimeLimit, 10*sistema_a->GetNCenarios());		// Limita o tempo de resolução do problema
	}
}
Spcdec2SPH::~Spcdec2SPH(void)
{
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
}

void Spcdec2SPH::CriarVariaveis(int ch)
{
	try 
	{
		// variáveis em cada problema somente para uma cascata
		// Estagio 1
		//x = [ph v d s (phmax) phg q z]
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades();
				modelosGRB[ch]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < cascata[ch].size(); r++)	//v
				modelosGRB[ch]->addVar(double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin(t)), double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax(t)), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < cascata[ch].size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
				modelosGRB[ch]->addVar(0, qhmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < cascata[ch].size(); r++)	//s
			{
				modelosGRB[ch]->addVar(0, sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagPhmax() == 1)
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)	//phmax
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades();
					modelosGRB[ch]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < cascata[ch].size(); r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
				{
					modelosGRB[ch]->addVar(0, double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < cascata[ch].size(); r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
				{
					modelosGRB[ch]->addVar(0, double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
						modelosGRB[ch]->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
						modelosGRB[ch]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
		}
		// Estagio 2
		//x = [ph v d s (phmax) phg q z vfol]
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (int t = 0; t < sistema_a->GetTt2() - sistema_a->GetTt1(); t++)
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)	//ph
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades();
					modelosGRB[ch]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < cascata[ch].size(); r++)	//v
					modelosGRB[ch]->addVar(double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin(t+sistema_a->GetTt1()+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), double(sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax(t+sistema_a->GetTt1()+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < cascata[ch].size(); r++)	//d
				{
					double qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
						qhmax = qhmax + double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());
					modelosGRB[ch]->addVar(0, qhmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < cascata[ch].size(); r++)	//s
				{
					modelosGRB[ch]->addVar(0, sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagPhmax() == 1)
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)	//phmax
					{
						double phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
							phmax = phmax + sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades();
						modelosGRB[ch]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
					}
				}
				for (size_t r = 0; r < cascata[ch].size(); r++)	//phg
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
					{
						modelosGRB[ch]->addVar(0, double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				for (size_t r = 0; r < cascata[ch].size(); r++)	//q
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
					{
						modelosGRB[ch]->addVar(0, double (sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)
						for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
							modelosGRB[ch]->addVar(0, 1, 0.0, GRB_BINARY, "");
				}
				else
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)
						for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)	
							modelosGRB[ch]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
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
void Spcdec2SPH::CriarRestricoes(int ch)
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
				//4278
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
void Spcdec2SPH::CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int ch)
{
	try 
	{
		//x = [ph v d s (phmax) phg q z vfol]
		// Termos lineares
		int n_a = n[ch] - flag2*sistema_a->GetNCenarios()*cascata[ch].size();
		int nt;
		int nt_dual = lambda->size() / T;
		int delta_dual = 0;
		double deltaT;
		int delta = 0;
		// vfolga
		int JJ = 0;
		for (size_t r = 0; r < cascata[ch].size(); r++)
			JJ += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
		delta = ((4+flag3)*cascata[ch].size() + 3*JJ) * (sistema_a->GetTt2());
		deltaT = sistema_a->GetDeltaT2();
		if (sistema_a->GetFlagVfol() == true)
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)
					vars[ch][delta + r].set(GRB_DoubleAttr_Obj, (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//vfol
				delta += ((4+flag3)*cascata[ch].size() + 3*JJ) * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[ch].size();
			}
		}
	
		// Termos da RL
		//L = [Lpt Lph (Lphmax/Lres)]
		// Lineares
		delta = 0;
		delta_dual = 0;
		if (sistema_a->GetFlagPhmax() == 1)
		{
			for (int t = 0; t < T; t++)
			{
				for (size_t r = 0; r < cascata[ch].size(); r++)
				{
					vars[ch][r + delta].set(GRB_DoubleAttr_Obj, - lambda->at(cascata[ch].at(r) + sistema_a->termeletricasVtr.size() + delta_dual));			//ph
					vars[ch][r + 4*cascata[ch].size() + delta].set(GRB_DoubleAttr_Obj, - lambda->at(cascata[ch].at(r) + sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size() + delta_dual));			//phmax
				}
				if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
					nt = (n_a / T) + flag2*cascata[ch].size();
				else
					nt = (n_a / T);
				delta += nt;
				delta_dual += nt_dual;
			}
		}
		else
		{
			int jj = 0;
			for (int t = 0; t < T; t++)
			{
				jj = delta;
				for (size_t r = 0; r < cascata[ch].size(); r++)
				{
					vars[ch][r + delta].set(GRB_DoubleAttr_Obj, - lambda->at(cascata[ch].at(r) + sistema_a->termeletricasVtr.size() + delta_dual) + lambda->at(sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size() + delta_dual));			//ph
					for (int j = 0; j < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos(); j++)
						vars[ch][4*cascata[ch].size() + 2*JJ + jj + j].set(GRB_DoubleAttr_Obj, - lambda->at(sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size() + delta_dual)*sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[cascata[ch].at(r)].grupoVtr[j].GetNUnidades());		//z
					jj += sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetNGrupos();
				}
				if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
					nt = (n_a / T) + flag2*cascata[ch].size();
				else
					nt = (n_a / T);
				delta += nt;
				delta_dual += nt_dual;
			}
		
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
int Spcdec2SPH::ResolverProblemaRL(Spcdec2Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int ch)
{
	// Resolve o subproblema de uma cascata de usinas
	try
	{
		CriarFuncaoObjetivoRL(lambda, ch);
		modelosGRB[ch]->update();	// Atualiza o modelo Gurobi.

		// ajustes do solver
		modelosGRB[ch]->reset();

		//modelosGRB[ch]->write("SpcDec2_subpH.lp");	// Escreve modelo em arquivo

		// Otimizar
		modelosGRB[ch]->optimize();

		// Salvar resultados
		x_ch.clear();
		x_ch.resize(n[ch]);
		int nStatus = modelosGRB[ch]->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The HE subproblem was not solved until optimality! " << nStatus << endl;
		if (modelosGRB[ch]->get(GRB_IntAttr_SolCount))	// se existir ao menos uma solução
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
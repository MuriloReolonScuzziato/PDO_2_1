#include "Scndec1SPC.h"

// Criar restrições
// ------------------------------------------------
// Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa Scndec1SPC::MatrizRestDemanda()	// Monta matriz da restrição de atendimento a demanda
{
	// T * B
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
	CMatrizEsparsa Agtu(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
	if (sistema_a->GetFlagTbinaryModel() == 1)
	{
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
					if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
						Agtu.InserirElemento(b, i, sistema_a->termeletricasVtr[i].GetPmin());
	}
		
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
	//
	CMatrizEsparsa Alin(T * sistema_a->barrasVtr.size(), n_a);
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Adef(sistema_a->barrasVtr.size());
	int l = 0;
	int c = 0;
	for ( int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l, c, l + sistema_a->barrasVtr.size() - 1, c + sistema_a->termeletricasVtr.size() - 1, &Agt, 0, 0);
		if (sistema_a->GetFlagTbinaryModel() == 1)
			Alin.InserirMatriz(l, c + sistema_a->termeletricasVtr.size(), l + sistema_a->barrasVtr.size() - 1, c + 2*sistema_a->termeletricasVtr.size() - 1, &Agtu, 0, 0);
		Alin.InserirMatriz(l, c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),l + sistema_a->barrasVtr.size() - 1, c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &B, 0, 0);
		Alin.InserirMatriz(l, c + (sistema_a->barrasVtr.size() - 1 + (3+flag4+flag7)*sistema_a->termeletricasVtr.size()), l + sistema_a->barrasVtr.size() - 1, c + (sistema_a->barrasVtr.size() - 1 + (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size()) - 1, &Agh, 0, 0);
		Alin.InserirMatriz(l, c + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + sistema_a->barrasVtr.size() - 1 + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ, l + sistema_a->barrasVtr.size() - 1, c + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + sistema_a->barrasVtr.size() - 1 + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ + sistema_a->barrasVtr.size() - 1, &Adef, 0, 0);
		l = l + sistema_a->barrasVtr.size();
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizRestDemandaBarraUnica()	// Monta matriz da restrição de atendimento a demanda
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
	int c = 0;
	for ( int t = 0; t < T; t++)
	{
		a.RemoverTodosElementos();
		for (size_t i = 0; i < I; i++)		// pt
			a.InserirElemento(0, i + c, 1);
		if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
		{
			for (size_t i = 0; i < I; i++)		// u
				a.InserirElemento(0, i + I + c, sistema_a->termeletricasVtr[i].GetPmin());
		}
		for (size_t r = 0; r < R; r++)		// ph
			a.InserirElemento(0, r + c + (3+flag4+flag7)*sistema_a->termeletricasVtr.size(), 1);
		if ( flag1 == 0)
			a.InserirElemento(0, (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ + c, 1);		// def
		else // flag1 == 2
			for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)		// def
				a.InserirElemento(0, b + c + (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ, 1);		// def
		Alin.JuntarColuna(&a);
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizLimFluxo()
{
	int n_a = n;
	CMatrizEsparsa Alin(T * sistema_a->linhasVtr.size(), n_a);
	if ( flag1 == 1 )
	{
		CMatrizEsparsa Alb(sistema_a->linhasVtr.size(),sistema_a->barrasVtr.size());
		for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				if (sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[b])
					Alb.InserirElemento(l, b, 1);
				else if (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[b])
					Alb.InserirElemento(l, b, -1);			
		Alb.RemoverColuna(sistema_a->GetBarraRef() - 1);
		CMatrizEsparsa TT(sistema_a->linhasVtr.size(), sistema_a->linhasVtr.size());
		for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
			for ( size_t ll = 0; ll < sistema_a->linhasVtr.size(); ll++)
				if (l == ll)
					TT.InserirElemento(l, ll, 100 / (sistema_a->linhasVtr[l].GetReatancia()));
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		int ll = 0;
		int c = 0;
		TT.MultiplicarPorMatriz(&Alb);
		for ( int t = 0; t < T; t++)
		{
			Alin.InserirMatriz(ll,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),ll + sistema_a->linhasVtr.size() - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
			ll = ll + sistema_a->linhasVtr.size();
			c += (n_a / T);
		}
	}
	else // flag1 == 2
	{
		// Montar matriz de incidencia de hidros (ph), térmicas (pt e u) e déficit (def): Agh, Agt, Agtu Adef
		MatrixXd Agh = MatrixXd::Zero(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
					if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
						Agh(b, r) = 1;
		MatrixXd Agt = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
					if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
						Agt(b, i) = 1;
		MatrixXd Agtu = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
		if (sistema_a->GetFlagTbinaryModel() == 1)
		{
			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
						if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
							Agtu(b, i) = sistema_a->termeletricasVtr[i].GetPmin();
		}
		MatrixXd Adef = MatrixXd::Identity(sistema_a->barrasVtr.size(),sistema_a->barrasVtr.size());
		// Multiplicar - Beta pelas matrizes acima para cada conjunto de variáveis e incluir em Alin
		Agh = - sistema_a->Beta * Agh; //Agh = - Beta * Agh;
		Agt = - sistema_a->Beta * Agt; //Agt = - Beta * Agt;
		Agtu = - sistema_a->Beta * Agtu; //Agtu = - Beta * Agtu;
		Adef = - sistema_a->Beta * Adef; //Adef = - Beta * Adef;
		
		//
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		int JJ = 0;
		for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
			JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
		int l = 0;
		int c = 0;
		for ( int t = 0; t < T; t++)
		{
			Alin.InserirMatriz(l, c, l + sistema_a->linhasVtr.size() - 1, c + sistema_a->termeletricasVtr.size() - 1, Agt, 0, 0);
			if (sistema_a->GetFlagTbinaryModel() == 1)
				Alin.InserirMatriz(l, c + sistema_a->termeletricasVtr.size(), l + sistema_a->linhasVtr.size() - 1, c + 2*sistema_a->termeletricasVtr.size() - 1, Agtu, 0, 0);
			Alin.InserirMatriz(l, c + (3+flag4+flag7)*sistema_a->termeletricasVtr.size(), l + sistema_a->linhasVtr.size() - 1, c + (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size() - 1, Agh, 0, 0);
			Alin.InserirMatriz(l, c + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ, l + sistema_a->linhasVtr.size() - 1, c + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ + sistema_a->barrasVtr.size() - 1, Adef, 0, 0);
			l = l + sistema_a->linhasVtr.size();
			c += (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizLimPhgMin()
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
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizLimPhgMax()
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
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizBalHid()
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
				MSoma2.InserirMatriz(0, 0, R - 1, R - 1, &Ad, t*R, tt*R);
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
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int cc = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &Av, 0, 0);
			for (int tt = 0; tt <= t; tt++)
				Alin.InserirMatriz(l,tt*(n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R),l + R - 1,tt*(n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 3*R) - 1, &Ad, t*R, tt*R);
			if (t > 0)
				Alin.InserirMatriz(l,c - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
		}
		else
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &Av, 0, 0);
			cc = 0;
			for (int tt = 0; tt <= t; tt++)
			{
				Alin.InserirMatriz(l,cc + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R),l + R - 1,cc + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 3*R) - 1, &Ad, t*R, tt*R);
				if (((tt + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((tt + 1 - sistema_a->GetTt1()) > 0))
					cc += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
				else
					cc += (n_a / T);
			}
			if (t > 0)
				if (t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) == sistema_a->GetTt1())
					Alin.InserirMatriz(l,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
				else
					Alin.InserirMatriz(l,c - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
		}
		l = l + R;
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizLimqMin()
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
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ,l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizLimqMax()
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
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ,l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizTup()
{
	// Tup é o tempo que a unidade deve ficar ligada após ser ligada, portanto para ficar ligada por 3 periodos, Tup = 2;
	if (flag7 == 0)
	{
		// T * sum(Tup)
		int n_a = n;
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tup1, Tup2, Tup;
		int c = I;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		for (int t = 0; t < T; t++)
		{
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
						if (t - tt - 1 >= 0)		// cada if adiciona um elemento dos dois elementos do lado direito da equação
							a.InserirElemento(0, c + i - tt*(n_a / T) - 1*(n_a / T), -1);
						if (t - tt - 2 >= 0)
							a.InserirElemento(0, c + i - tt*(n_a / T) - 2*(n_a / T), 1);
					
						Alin.JuntarColuna(&a);
					}
				}
			}
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T * I (para usinas com Tup > 0)	nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < Tup
		int n_a = n;
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tup1, Tup2, Tup;
		int c = I;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		for (int t = 0; t < T; t++)
		{
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

				if (Tup > 0)		// só adiciona restrição para Tup > 0
				{
					a.RemoverTodosElementos();
					//if  (t >= Tup + 1)
					if  (t >= Tup)
					{
						a.InserirElemento(0, c + i, -1);		// adiciona termo de u
						for (int tt = 0; tt <= Tup; tt++)
							a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);		// adiciona cada termo do somatório (up_t)
					}
					else if ((Tup >= sistema_a->GetTt2()) && (t == sistema_a->GetTt2() - 1))
					{
						a.InserirElemento(0, c + i, -1);		// adiciona termo de u
						for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
							a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);	// adiciona cada termo do somatório (up_t)
					}
					Alin.JuntarColuna(&a);
				}
				else {}
			}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Scndec1SPC::MatrizTdown()
{
	// O valor de Tdown tem significado similar ao Tup
	if (flag7 == 0)
	{
		// T * sum(Tdown)
		int n_a = n;
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tdown1, Tdown2, Tdown;
		int c = I;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		for (int t = 0; t < T; t++)
		{
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
						if (t - tt - 1 >= 0)		// cada if adiciona um elemento dos dois elementos do lado direito da equação
							a.InserirElemento(0, c + i - tt*(n_a / T) - 1*(n_a / T), -1);
						if (t - tt - 2 >= 0)
							a.InserirElemento(0, c + i - tt*(n_a / T) - 2*(n_a / T), 1);
					
						Alin.JuntarColuna(&a);
					}
				}
			}
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T * I (para usinas com Tdown > 0)	nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < Tdown
		int n_a = n;
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tdown1, Tdown2, Tdown;
		int c = I;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		for (int t = 0; t < T; t++)
		{
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
					a.RemoverTodosElementos();
					if  (t >= Tdown)
					{
						a.InserirElemento(0, c + i, 1);
						for (int tt = 0; tt <= Tdown; tt++)
							a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
					}
					else if ((Tdown >= sistema_a->GetTt2()) && (t == sistema_a->GetTt2() - 1))
					{
						a.InserirElemento(0, c + i, 1);
						for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
							a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
					}
					Alin.JuntarColuna(&a);
				}
			}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Scndec1SPC::MatrizTupDown()
{
	// T * I	(nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < 1
	int n_a = n;
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int c = I;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			if (t >= 1)
			{
				a.InserirElemento(0, c + i, 1);
				a.InserirElemento(0, c + i + I, -1);
				a.InserirElemento(0, c + i + 2*I, 1);
				a.InserirElemento(0, c + i - 1*(n_a / T), -1);
			}
			else
			{
				a.InserirElemento(0, c + i, 1);
				a.InserirElemento(0, c + i + I, -1);
				a.InserirElemento(0, c + i + 2*I, 1);
			}
			Alin.JuntarColuna(&a);
		}
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizRampaUp()
{
	// T * I
	if (flag7 == 0)
	{
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, 1);
				if (t > 0)
				{
					a.InserirElemento(0, c + i - (n_a / T), -1);
					a.InserirElemento(0, c + i + I - (n_a / T), sistema_a->termeletricasVtr[i].GetPmin() - sistema_a->termeletricasVtr[i].GetRampaUp());
				}
				Alin.JuntarColuna(&a);
			}
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, 1);
				if (t > 0)
					a.InserirElemento(0, c + i - (n_a / T), -1);
				Alin.JuntarColuna(&a);
			}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Scndec1SPC::MatrizRampaDown()
{
	// T * I
	if (flag7 == 0)
	{
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, -1);
				a.InserirElemento(0, c + i + I, sistema_a->termeletricasVtr[i].GetPmin() - sistema_a->termeletricasVtr[i].GetRampaDown());
				if (t > 0)
					a.InserirElemento(0, c + i - (n_a / T), 1);
				Alin.JuntarColuna(&a);
			}
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, -1);
				if (t > 0)
					a.InserirElemento(0, c + i - (n_a / T), 1);
				Alin.JuntarColuna(&a);
			}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Scndec1SPC::MatrizLimPtMin()
{
	if (flag7 == 0)
	{
		// T * I
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
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
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T * I
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
		int c = 0;
		c += (n_a / T);						// referente à t = 0
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
			Alin.JuntarColuna(&a);
		}
		for (int t = 1; t < T; t++)			// nó inicial, já adicionadas restrições, por isso t começa em 1
		{
			for (size_t i = 0; i < I; i++)
			{
				if (sistema_a->termeletricasVtr[i].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
					a.InserirElemento(0, c + i - (n_a / T), 1);		// pt-1
					a.InserirElemento(0, c + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u-1
					Alin.JuntarColuna(&a);
				}
				else
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
					a.InserirElemento(0, c + i - (n_a / T), 1);		// pt-1
					a.InserirElemento(0, c + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u-1
					a.InserirElemento(0, c + i - (n_a / T) + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// up-1
					Alin.JuntarColuna(&a);
				}
			}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Scndec1SPC::MatrizLimPtMax()
{
	// T * I
	if (flag7 == 0)
	{
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
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
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T * I
		int n_a = n;
		n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, 1);
				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
				a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
				Alin.JuntarColuna(&a);
			}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Scndec1SPC::MatrizRestCP()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
	int c = I;
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + i + I, 1);
			a.InserirElemento(0, c + i, - sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
			if (t > 0)
				a.InserirElemento(0, c + i - (n_a / T), sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
			Alin.JuntarColuna(&a);
		}
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizVmeta()
{
	int n_a = n;
	size_t R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	for (size_t i = 0; i < R; i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->hidreletricasVtr.size());
	int c;
	c = (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R + (sistema_a->GetTt2() - 1)*(n_a / T);
	for (size_t r = 0; r < R; r++)
	{
		a.RemoverTodosElementos();
		a.InserirElemento(0, c + r, 1);
		if (sistema_a->GetFlagVfol() == true)
			a.InserirElemento(0, c + r + (3+flag3)*R + 3*JJ + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d), 1);
		Alin.JuntarColuna(&a);
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizBalPotencia()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B;
	int jj = dd + (4+flag3)*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		nt = (n_a / T);
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
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
CMatrizEsparsa Scndec1SPC::MatrizBalVazao()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B + 2*R;
	int jj = dd + JJ + (2+flag3)*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
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
CMatrizEsparsa Scndec1SPC::MatrizFuncProd()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	size_t I = sistema_a->termeletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->hidreletricasVtr.size());
	int vv = (3+flag4+flag7)*I + B + R;
	int jj = vv + (3+flag3)*R + JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
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
						a.InserirElemento(0, jj + j + JJ, sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax());				// z -> adiciona um termo (1 - z)*Pmin na função, para phg <= pmin qdo z = 0; o rhs estava negativo, devido a aproximação linear, fazendo com q a hidro ficasse obrigatoriamente ligada
					if (sistema_a->hidreletricasVtr[r].GetInflueVert() == 1 )		// Se o s influencia já tenho contabilizado isso no d, caso contrário só subtrair o valor de s de d
					{}
					else		// vertimento n influencia no canal de fuga, entao subtraio o valor de s da defluencia
						a.InserirElemento(0, vv + 2*R + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgS(napr));		// s
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
CMatrizEsparsa Scndec1SPC::MatrizPhMax()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B + 4*R;
	int jj = dd + R + 2*JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
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
CMatrizEsparsa Scndec1SPC::MatrizReserva()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B;
	int jj = dd + 4*R + 2*JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		nt = (n_a / T);
		a.RemoverTodosElementos();
		for (size_t r = 0; r < R; r++)
		{
			a.InserirElemento(0, dd + r, - 1);
			if (sistema_a->GetFlagPhmax() == 1)
				a.InserirElemento(0, dd + 4*R + r, 1);
			else
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					a.InserirElemento(0, jj + j, sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		Alin.JuntarColuna(&a);
		dd += nt;
		jj += nt - JJ;
	}
	return Alin;
}
CMatrizEsparsa Scndec1SPC::MatrizCortesF()
{
	// T * I
	int n_a = n;
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)		// Adiciona uma restrição por corte
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + i, - sistema_a->termeletricasVtr[i].GetCoefA1(n_cort));
				if ( flag7 == 0 )
					a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetCoefA0(n_cort));
				else
					a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetCoefA0(n_cort) - sistema_a->termeletricasVtr[i].GetCoefA1(n_cort) * sistema_a->termeletricasVtr[i].GetPmin());
				a.InserirElemento(0, c + i + (3+flag7)*I, 1);
				Alin.JuntarColuna(&a);
			}
		}
		c += (n_a / T);
	}
	return Alin;
	// I * T
	//int n_a = n;
	//n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	//size_t I = sistema_a->termeletricasVtr.size();
	//CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->hidreletricasVtr.size());
	//CMatrizEsparsa a(1, n_a + flag2*sistema_a->hidreletricasVtr.size());
	//int c;
	//for (size_t i = 0; i < I; i++)
	//{	
	//	c = i;
	//	for (int t = 0; t < T; t++)
	//	{
	//		for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)		// Adiciona uma restrição por corte
	//		{
	//			a.RemoverTodosElementos();
	//			a.InserirElemento(0, c, - sistema_a->termeletricasVtr[i].GetCoefA1(n_cort));
	//			a.InserirElemento(0, c + I, - sistema_a->termeletricasVtr[i].GetCoefA0(n_cort));
	//			a.InserirElemento(0, c + 3*I, 1);
	//			Alin.JuntarColuna(&a);
	//		}
	//		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
	//			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
	//		else
	//			c = c + (n_a / T);
	//	}
	//}
	//return Alin;
}
// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Scndec1SPC::LimRestDemanda()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, B * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		for (size_t b = 0; b < B; b++)
		{
			double D = 0;
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);		// cen = 0, ou qq outro cenários, pois os T1 primeiros periodos são iguais (deterministicos)
			}
			else
			{
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cenario, t);
			}
		Lim[t*B + b][0] = 1;
		Lim[t*B + b][1] = D;
		}
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimRestDemandaBarraUnica()
{
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
			{
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);
			}
			else
			{
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cenario, t);
			}
		}
		Lim[t][0] = 1;
		Lim[t][1] = D;
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimFluxo0()
{
	vetorfloat2 Lim;
	size_t L = sistema_a->linhasVtr.size();
	size_t B = sistema_a->barrasVtr.size();
	DimensionarMatriz(&Lim, L * T, 2);
	if ( flag1 == 1 )
	{
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
			for (int l = 0; l < L; l++)
			{
				Lim[t*L + l][0] = 0;
				Lim[t*L + l][1] = sistema_a->linhasVtr[l].GetCapacidade();
			}
	}
	else	// flag1 == 2
	{
		VectorXd Ad(sistema_a->barrasVtr.size());
		for (int t = 0; t < T; t++)
		{
			Ad.resize(sistema_a->barrasVtr.size());
			for (size_t b = 0; b < B; b++)
			{
				double D = 0;
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cenario, t);
				Ad(b) = D;
			}
			Ad = - sistema_a->Beta * Ad; //Ad = - Beta * Ad;		// aqui o vetor passa a ser de tamanho L

			for (int l = 0; l < L; l++)
			{
				Lim[t*L + l][0] = 0;
				Lim[t*L + l][1] = sistema_a->linhasVtr[l].GetCapacidade() + Ad(l);
			}
		}
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimFluxo2()
{
	vetorfloat2 Lim;
	size_t L = sistema_a->linhasVtr.size();
	size_t B = sistema_a->barrasVtr.size();
	DimensionarMatriz(&Lim, L * T, 2);
	if ( flag1 == 1 )
	{
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
			for (int l = 0; l < L; l++)
			{
				Lim[t*L + l][0] = 2;
				Lim[t*L + l][1] = - sistema_a->linhasVtr[l].GetCapacidade();
			}
	}
	else	// flag1 == 2
	{
		VectorXd Ad(sistema_a->barrasVtr.size());
		for (int t = 0; t < T; t++)
		{
			Ad.resize(sistema_a->barrasVtr.size());
			for (size_t b = 0; b < B; b++)
			{
				double D = 0;
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cenario, t);
				Ad(b) = D;
			}
			Ad = - sistema_a->Beta * Ad; //Ad = - Beta * Ad;		// aqui o vetor passa a ser de tamanho L

			for (int l = 0; l < L; l++)
			{
				Lim[t*L + l][0] = 2;
				Lim[t*L + l][1] = - sistema_a->linhasVtr[l].GetCapacidade() + Ad(l);
			}
		}
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimPhgMin()
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
				//Lim(jj + j,1) = double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimPhgMax()
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
vetorfloat2 Scndec1SPC::LimBalHid()
{
	double deltaT;
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0); 
	for (int t = 0; t < T; t++)
	{
		if ((0 <= t) && (sistema_a->GetTt1() > t))
			deltaT = sistema_a->GetDeltaT1();
		else
			deltaT = sistema_a->GetDeltaT2();

		for (int r = 0; r < R; r++)
		{
			if (t == 0)
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (sistema_a->hidreletricasVtr[r].GetV0() + 0.0036*deltaT*sistema_a->hidreletricasVtr[r].GetAfluencia(cenario, t));
			}
			else
			{
				Lim[t*R + r][0] = 1;
				Lim[t*R + r][1] = double (0.0036*deltaT*sistema_a->hidreletricasVtr[r].GetAfluencia(cenario, t));
			}
		}
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimQMin()
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
vetorfloat2 Scndec1SPC::LimQMax()
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
vetorfloat2 Scndec1SPC::LimTup()
{
	if (flag7 == 0)
	{
		int n_a = n;
		int I = sistema_a->termeletricasVtr.size();
		int Tup1, Tup2, Tup;
		//
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
		{
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
						if (t - tt <= 0)		// if referente ao elemento tt (se ele for constante entra nesse loop)
						{
							if (t - tt > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
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
						if (t - tt - 1 <= 0)	// if referente ao elemento tt - 1 (se ele for constante entra nesse loop)
						{
							if (t - tt - 1 > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
							{
								L[1] += - sistema_a->termeletricasVtr[i].GetU0();
							}
							else	// estado da usina antes do perido -x0
							{	
								L[1] += - (1 - sistema_a->termeletricasVtr[i].GetU0());
							}
						}
						Lim.push_back(L);
						L.clear();L.resize(2);
					}
				}
			}
		}
		return Lim;
	}
	else
	{
		// T * I
		int n_a = n;
		int I = sistema_a->termeletricasVtr.size();
		int Tup1, Tup2, Tup;
		//
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0; L[1] = 0;
		for (int t = 0; t < T; t++)
		{
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
					Lim.push_back(L);
			}
		}
		return Lim;
	}
}
vetorfloat2 Scndec1SPC::LimTdown()
{
	if (flag7 == 0)
	{
		int n_a = n;
		int I = sistema_a->termeletricasVtr.size();
		int Tdown1, Tdown2, Tdown;
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0;
		L[1] = 1;
		for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
		{
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
						if (t - tt <= 0)		// if referente ao elemento tt (se ele for constante entra nesse loop)
						{
							if (t - tt > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
								L[1] += double (sistema_a->termeletricasVtr[i].GetU0());
							else	// estado da usina antes do perido -x0
								L[1] += double((1 - sistema_a->termeletricasVtr[i].GetU0()));
						}
						if (t - tt - 1 <= 0)	// if referente ao elemento tt - 1 (se ele for constante entra nesse loop)
						{
							if (t - tt - 1 > - sistema_a->termeletricasVtr[i].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
								L[1] += - sistema_a->termeletricasVtr[i].GetU0();
							else	// estado da usina antes do perido -x0
								L[1] += - (1 - sistema_a->termeletricasVtr[i].GetU0());
						}
						Lim.push_back(L);
						L.clear();L.resize(2);L[0] = 0;L[1] = 1;
					}
				}
			}
		}
		return Lim;
	}
	else
	{
		int I = sistema_a->termeletricasVtr.size();
		int Tdown1, Tdown2, Tdown;
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0;L[1] = 1;
		for (int t = 0; t < T; t++)
		{
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
					Lim.push_back(L);
			}
		}
		return Lim;
	}
}
vetorfloat2 Scndec1SPC::LimTupDown()
{
	// T * I	(nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < 2 ou t < Tup + 2 ou t < Tdown + 2
	size_t I = sistema_a->termeletricasVtr.size();
	vetorfloat2 Lim;
	vetorfloat L, L0;
	L.resize(2); L0.resize(2);
	L0[0] = 1;
	L[0] = 1;L[1] = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t i = 0; i < I; i++)
		{
			if (t >= 1)
				Lim.push_back(L);
			else
			{
				L0[1] = sistema_a->termeletricasVtr[i].GetU0();
				Lim.push_back(L0);
			}
		}
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimRampaUp()
{
	if (flag7 == 0)
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
	else
	{
		// T * I
		int n_a = n;
		int I = sistema_a->termeletricasVtr.size();
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0;
		for (int t = 0; t < T; t++)
		{
			for (int i = 0; i < I; i++)
			{
				if (t == 0)
					L[1] = sistema_a->termeletricasVtr[i].GetRampaUp() + (sistema_a->termeletricasVtr[i].GetPt0() - sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetU0());
				else
					L[1] = sistema_a->termeletricasVtr[i].GetRampaUp();
				Lim.push_back(L);
			}
		}
		return Lim;
	}
}
vetorfloat2 Scndec1SPC::LimRampaDown()
{
	if (flag7 == 0)
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
	else
	{
		// T * I
		int n_a = n;
		int I = sistema_a->termeletricasVtr.size();
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0;
		for (int t = 0; t < T; t++)
		{
			for (int i = 0; i < I; i++)
			{
				if (t == 0)
					L[1] = sistema_a->termeletricasVtr[i].GetRampaDown() - (sistema_a->termeletricasVtr[i].GetPt0() - sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetU0());
				else
					L[1] = sistema_a->termeletricasVtr[i].GetRampaDown();
				Lim.push_back(L);
			}
		}
		return Lim;
	}
}
vetorfloat2 Scndec1SPC::LimPtMin()
{
	if (flag7 == 0)
	{
		// T * I
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
	else
	{
		// T * I
		int I = sistema_a->termeletricasVtr.size();
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0; L[1] = 0;
		for (size_t i = 0; i < I; i++)
		{
			if (sistema_a->termeletricasVtr[i].GetU0() == 1)
				if (sistema_a->termeletricasVtr[i].GetX0() > 1)  // indica que a usina foi ligada antes do periodo t = -1, portanto up-1 = 0;
					L[1] = sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPt0();		// valor de pt0 está em valor de ptmin a ptmax e n de 0 a ptmax-ptmin
			Lim.push_back(L);
			L[0] = 0; L[1] = 0;
		}
		for (int t = 1; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
				Lim.push_back(L);
		}
		return Lim;
	}
}
vetorfloat2 Scndec1SPC::LimPtMax()
{
	if (flag7 == 0)
	{
		// T * I
		int I = sistema_a->termeletricasVtr.size();
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, I * T, 2);
		IniciaMatriz(&Lim, 0);
		return Lim;
	}
	else
	{
		// T * I
		int I = sistema_a->termeletricasVtr.size();
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0; L[1] = 0;
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				Lim.push_back(L);
			}
		}
		return Lim;
	}
}
vetorfloat2 Scndec1SPC::LimRestCP()
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
vetorfloat2 Scndec1SPC::LimVmeta()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R, 2);
	IniciaMatriz(&Lim, 0);
	for (int r = 0; r < R; r++)
	{
		Lim[r][1] = sistema_a->hidreletricasVtr[r].GetVMeta();
		Lim[r][0] = 2;
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimBalPotenciaL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Scndec1SPC::LimBalVazaoL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Scndec1SPC::LimFuncProdL()
{
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (int i = 0; i < R; i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
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
vetorfloat2 Scndec1SPC::LimPhMax()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Scndec1SPC::LimReserva()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1 * T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		Lim[t][0] = 2;
		// função de reserva é calculada para todos os cenários
		if (t < sistema_a->GetTt1())
			Lim[t][1] = sistema_a->GetReserva(t);
		else
			Lim[t][1] = sistema_a->GetReserva(t + cenario*(sistema_a->GetTt2() - sistema_a->GetTt1()));
	}
	return Lim;
}
vetorfloat2 Scndec1SPC::LimCortesF()
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
void Scndec1SPC::MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	int count = 0;
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		M = MatrizRestDemanda(); i_rest_cen[0] = *n_restricoes;
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	else
	{
		M = MatrizRestDemandaBarraUnica(); i_rest_cen[0] = *n_restricoes;
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		if ( sistema_a->GetFlagModeloRede() == 2 )
		{
			M = MatrizLimFluxo();
			SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
			*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
			M = MatrizLimFluxo();
			SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
			*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		}
	}
	M = MatrizLimPhgMin();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPhgMax();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimqMin();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimqMax();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalHid(); i_rest_cen[1] = *n_restricoes;
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizTup();									
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizTdown();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	if ( sistema_a->GetFlagTbinaryModel() == 1 )
	{
		M = MatrizTupDown();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	M = MatrizRampaUp();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizRampaDown();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPtMax();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizLimPtMin();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		M = MatrizRestCP();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	M = MatrizVmeta();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalPotencia();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizBalVazao();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	M = MatrizFuncProd();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	if (sistema_a->GetFlagPhmax() == 1)
	{
		M = MatrizPhMax();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	M = MatrizReserva(); i_rest_cen[2] = *n_restricoes;
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		M = MatrizCortesF();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
}
void Scndec1SPC::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor)
{
	vetorfloat2 L;
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		L = LimRestDemanda();
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo0();
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo2();
		AlocarLimites(&L, LimTipo, LimValor);
	}
	else
	{
		L = LimRestDemandaBarraUnica();
		AlocarLimites(&L, LimTipo, LimValor);
		if ( sistema_a->GetFlagModeloRede() == 2 )
		{
			L = LimFluxo0();
			AlocarLimites(&L, LimTipo, LimValor);
			L = LimFluxo2();
			AlocarLimites(&L, LimTipo, LimValor);
		}
	}
	L = LimPhgMin();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPhgMax();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMin();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimQMax();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalHid();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimTup();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimTdown();
	AlocarLimites(&L, LimTipo, LimValor);
	if ( sistema_a->GetFlagTbinaryModel() == 1 )
	{
		L = LimTupDown();
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimRampaUp();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimRampaDown();
	AlocarLimites(&L, LimTipo, LimValor);
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		L = LimRestCP();
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimPtMax();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPtMin();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimVmeta();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalPotenciaL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalVazaoL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimFuncProdL();
	AlocarLimites(&L, LimTipo, LimValor);
	if (sistema_a->GetFlagPhmax() == 1)
	{
		L = LimPhMax();
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimReserva();
	AlocarLimites(&L, LimTipo, LimValor);
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		L = LimCortesF();
		AlocarLimites(&L, LimTipo, LimValor);
	}
}
// ------------------------------------------------

Scndec1SPC::Scndec1SPC(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;
	cenario = 0;
	modeloGRB = new GRBModel(ambiente_gurobi);		// um modelo de otimização para o cenário
	i_rest_cen.resize(3);
	T = sistema_a->GetTt2();
	flag1 = int (sistema_a->GetFlagModeloRede());	
	flag2 = int (sistema_a->GetFlagVfol());
	flag3 = int (sistema_a->GetFlagPhmax());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());
	int JJ = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();

	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size();
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= T * (sistema_a->barrasVtr.size() - 1);

	// Criar variáveis
	CriarVariaveis();
	modeloGRB->update();	// Atualiza o modelo Gurobi.
	vars = modeloGRB->getVars();

	// Adicionar restrições
	CriarRestricoes();
	modeloGRB->update();	// Atualiza o modelo Gurobi.

	// Fixar condições iniciais
	FixarCondIniciais();
	modeloGRB->update();	// Atualiza o modelo Gurobi.

	nStatus = 2;
	sol_mipgap = 1e-4;
	sol_mipgap_curr = sol_mipgap;
	modeloGRB->getEnv().set(GRB_DoubleParam_IntFeasTol, sol_mipgap);

	// Definir ajustes do solver
	//modelosGRB->write("probH.lp");									// Escreve modelo em arquivo
	modeloGRB->getEnv().set(GRB_IntParam_OutputFlag, 0);				// Nao escrever detalhes
	//modelosGRB.getEnv().set(GRB_IntParam_LogToConsole, 0);			// Nao escrever detalhes na tela
	//modelosGRB.getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo
	modeloGRB->getEnv().set(GRB_DoubleParam_MIPGap, sol_mipgap);		// Define o gap de tolerancia
	//modeloGRB->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9);
	//modeloGRB->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);
	//modeloGRB->getEnv().set(GRB_DoubleParam_IntFeasTol, 1e-9);
	//modeloGRB->getEnv().set(GRB_IntParam_ScaleFlag, 0);				// Desabilita o model scaling
	//modelosGRB->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
	modeloGRB->getEnv().set(GRB_IntParam_Threads, 4);
	//modelosGRB->getEnv().set(GRB_IntParam_Presolve, 0);				//Desabilitar o presolve
	//modelosGRB.getEnv().set(GRB_IntParam_PreSparsify, 1);
	modeloGRB->getEnv().set(GRB_IntParam_Method, 0);
	//modelosGRB.getEnv().set(GRB_IntParam_MIPFocus, 2);
	//modelosGRB.getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);
	modeloGRB->getEnv().set(GRB_DoubleParam_TimeLimit, 60);			// Limita o tempo de resolução do problema
}
Scndec1SPC::Scndec1SPC(CSistema * const sistema_end, int cenario_a, Scndec1SPC &cenario0)
{
	sistema_a = sistema_end;
	cenario = cenario_a;
	// receber atributos da outra classe copia
	n = cenario0.n;
	T = cenario0.T;
	flag1 = cenario0.flag1;
	flag2 = cenario0.flag2;
	flag3 = cenario0.flag3;
	flag4 = cenario0.flag4;
	flag7 = cenario0.flag7;
	nStatus = cenario0.nStatus;
	sol_mipgap = cenario0.sol_mipgap;
	i_rest_cen = cenario0.i_rest_cen;
	modeloGRB = new GRBModel(*cenario0.modeloGRB);
	vars = modeloGRB->getVars();
	x.resize(n);

	//atualizar restrições que dependam dos cenários
	AtualizarRestricoes();
}
Scndec1SPC::~Scndec1SPC(void)
{
	delete modeloGRB;
	delete vars;
}

void Scndec1SPC::CriarVariaveis()
{
	try 
	{
		size_t I = sistema_a->termeletricasVtr.size();
		size_t R = sistema_a->hidreletricasVtr.size();
		size_t B = sistema_a->barrasVtr.size() - 1;
		//x = [pt u cp F teta ph v d s phmax phg q z def]
		//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model
		// Estagio 1 e 2
		for (int t = 0; t < T; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				for (size_t i = 0; i < I; i++)	//pt
					modeloGRB->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//pt
					modeloGRB->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagVarBin() == true)	//u (binario ou continuo)
			{
				for (size_t i = 0; i < I; i++)	//u
					modeloGRB->addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//u
					modeloGRB->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				if (sistema_a->GetFlagVarBin() == true)
				{
					for (size_t i = 0; i < I; i++)		//up
						modeloGRB->addVar(0, 1, 0.0, GRB_BINARY, "");
					for (size_t i = 0; i < I; i++)		//ud
						modeloGRB->addVar(0, 1, 0.0, GRB_BINARY, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)		//up
						modeloGRB->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					for (size_t i = 0; i < I; i++)		//ud
						modeloGRB->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//cp
					modeloGRB->addVar(0, sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
			{
				for (size_t i = 0; i < I; i++)
					modeloGRB->addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagModeloRede() == 1)
				for (size_t b = 0; b < B; b++)	//teta
					modeloGRB->addVar(-3.14159, 3.14159, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modeloGRB->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modeloGRB->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin()), double(sistema_a->hidreletricasVtr[r].GetVmax()), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				modeloGRB->addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
			{
				modeloGRB->addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagPhmax() == 1)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modeloGRB->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modeloGRB->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modeloGRB->addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
			{
				for (size_t r = 0; r < R; r++)	//z
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modeloGRB->addVar(0, 1, 0.0, GRB_BINARY, "");
				}
			}
			else
			{
				for (size_t r = 0; r < R; r++)	//z
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modeloGRB->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( sistema_a->GetFlagModeloRede() > 0)
			{
				for (size_t b = 0; b < B + 1; b++)	//def
				{
					double cap_d = 0;
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cenario, t);
					modeloGRB->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				double cap_d = 0;			//def barra unica
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cenario, t);
				modeloGRB->addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
			}
		}
		if ( sistema_a->GetFlagVfol() == true )
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//vfol
				modeloGRB->addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
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
void Scndec1SPC::CriarRestricoes()
{
	try
	{
		vetorint indexL, indexC, nnZ, LimTipo;
		vetorfloat indexV, LimValor;
		int n_restricoes = 0;
		indexL.resize(0);indexC.resize(0);indexV.resize(0);nnZ.resize(0);LimTipo.resize(0);LimValor.resize(0);
		nnZ.push_back(0);
		MatrizRestricoesLineares(&indexL, &indexC, &indexV, &nnZ ,&n_restricoes);
		MatrizLimitesLineares(&LimTipo, &LimValor);

		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		for (int i = 0; i < n_restricoes; i++)
		{
			for (int j = nnZ[i]; j < nnZ[i + 1]; j++)
			{
				coeficiente = indexV[j];
				variavel = vars[indexC[j]];
				restricao.addTerms( &coeficiente, &variavel, 1);
			}
			switch (LimTipo[i]) {
			case 0:
				modeloGRB->addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modeloGRB->addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modeloGRB->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
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
void Scndec1SPC::AtualizarRestricoes()
{
	try
	{
		GRBConstr * restricoes;
		restricoes = modeloGRB->getConstrs();
		// Atualizar restrições de atendimento a demanda
		int indexDemanda = i_rest_cen[0];
		if (sistema_a->GetFlagModeloRede() == 1)
		{
			vetorfloat2 atualDemanda = LimRestDemanda();
			for (size_t i = indexDemanda; i < indexDemanda + atualDemanda.size(); i++)
				restricoes[i].set(GRB_DoubleAttr_RHS, double (atualDemanda[i - indexDemanda][1]));
		}
		else
		{
			vetorfloat2 atualDemanda = LimRestDemandaBarraUnica();
			for (size_t i = indexDemanda; i < indexDemanda + atualDemanda.size(); i++)
				restricoes[i].set(GRB_DoubleAttr_RHS, double (atualDemanda[i - indexDemanda][1]));
			// Se o modelo de rede compacto é usado os limites dos fluxos também devem ser atualizados!
			if ( sistema_a->GetFlagModeloRede() == 2 )
			{
				indexDemanda += atualDemanda.size();
				vetorfloat2 atualfluxo0 = LimFluxo0();
				for (size_t i = indexDemanda; i < indexDemanda + atualfluxo0.size(); i++)
					restricoes[i].set(GRB_DoubleAttr_RHS, double (atualfluxo0[i - indexDemanda][1]));
				indexDemanda += atualfluxo0.size();
				vetorfloat2 atualfluxo2 = LimFluxo2();
				for (size_t i = indexDemanda; i < indexDemanda + atualfluxo2.size(); i++)
					restricoes[i].set(GRB_DoubleAttr_RHS, double (atualfluxo2[i - indexDemanda][1]));
			}
		}
		// Atualizar restrições de balanço hídrico
		int JJ = 0;
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int indexBalHid = i_rest_cen[1];
		vetorfloat2 atualBalHid = LimBalHid();
		for (size_t i = indexBalHid; i < indexBalHid + atualBalHid.size(); i++)
			restricoes[i].set(GRB_DoubleAttr_RHS, double (atualBalHid[i - indexBalHid][1]));
		// Atualizar restrições de reserva
		int indexReserva = i_rest_cen[2];
		vetorfloat2 atualReserva = LimReserva();
		for (size_t i = indexReserva; i < indexReserva + atualReserva.size(); i++)
			restricoes[i].set(GRB_DoubleAttr_RHS, double (atualReserva[i - indexReserva][1]));
		modeloGRB->update();
		delete restricoes;
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
void Scndec1SPC::CriarFuncaoObjetivoRL(const vetorfloat &coefRL)
{
	// Funçao Lagrangiana
	try 
	{
		int R = sistema_a->hidreletricasVtr.size();
		int I = sistema_a->termeletricasVtr.size();
		int B = sistema_a->barrasVtr.size();
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int nt = ( n - flag2*sistema_a->hidreletricasVtr.size() ) / T;
		int jj;
		double deltaT;
		int delta = 0;
		int flag1a = 0;		// referente à var. teta
		if (flag1 == 1)
			flag1a = 1;
		int flag1d = 1;		// referente à var. def
		if (flag1 == 0)
			flag1d = 0;
		//x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		//x = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		// Lambda:
		//L = [Lpt Lu Lcp (LF) (Lteta) Lph Lv Ld Ls (Lphmax) Lphg Lq Lz (Ldef)] -> (.) tamanho dos vetores que dependem de flags
		//L = [Lpt Lu Lup Lud (LF) (Lteta) Lph Lv Ld Ls (Lphmax) Lphg Lq Lz (Ldef)] -> (.) tamanho dos vetores que dependem de flags

		// Termos lineares
		for (int t = 0; t < T; t++)
		{
			if ((0 <= t) && (t < sistema_a->GetTt1()))		// nós do estágio 1
			{
				// as restrições de não antecipatividade só são adicionadas para o primeiro estágio
				deltaT = sistema_a->GetDeltaT1() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cenario);
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						vars[i + (3+flag7)*I + delta].set(GRB_DoubleAttr_Obj, deltaT + coefRL[i + (3+flag7)*I + delta]);	//F
						if (sistema_a->GetFlagTbinaryModel() == 0)
							vars[i + delta].set(GRB_DoubleAttr_Obj, coefRL[i + delta]);				//pt
						else
							vars[i + delta].set(GRB_DoubleAttr_Obj, coefRL[i + delta]);				//pt
						vars[i + I + delta].set(GRB_DoubleAttr_Obj, coefRL[i + I + delta]);			//u
					}
					else
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
						{
							vars[i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT + coefRL[i + delta]);	//pt
							vars[i + I + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT + coefRL[i + I + delta]);		//u
						}
						else
						{
							vars[i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT + coefRL[i + delta]);	//pt
							vars[i + I + delta].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT + coefRL[i + I + delta]);		//u
						}
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						vars[i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT + coefRL[i + 2*I + delta]);		//cp
					else
					{
						vars[i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + coefRL[i + 2*I + delta]);		//up
						vars[i + 3*I + delta].set(GRB_DoubleAttr_Obj, coefRL[i + 3*I + delta]);		//ud
					}
				}
				if (flag1 > 0)		// considerar rede
				{
					if (flag1 == 1)
					{
						for (size_t b = 0; b < B - 1; b++)
							vars[b + I*(3 + flag4 + flag7) + delta].set(GRB_DoubleAttr_Obj, coefRL[b + I*(3 + flag4 + flag7) + delta]);		//teta
					}
					for (size_t b = 0; b < B; b++)
						vars[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT + coefRL[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta]);		//def
				}
				else
					vars[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT + coefRL[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta]);		//def
				jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
				for (int r = 0; r < R; r++)
				{
					vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].set(GRB_DoubleAttr_Obj, coefRL[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta]);				//ph
					vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].set(GRB_DoubleAttr_Obj, coefRL[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta]);		//v
					vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].set(GRB_DoubleAttr_Obj, coefRL[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta]);	//d
					vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].set(GRB_DoubleAttr_Obj, coefRL[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta]);	//s
					if (sistema_a->GetFlagPhmax() == 1)
						vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].set(GRB_DoubleAttr_Obj, coefRL[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta]);	//phmax
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						vars[j + jj + delta].set(GRB_DoubleAttr_Obj, coefRL[j + jj + delta]);							//phg
						vars[j + jj + JJ + delta].set(GRB_DoubleAttr_Obj, coefRL[j + jj + JJ + delta]);				//q
						vars[j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, coefRL[j + jj + 2*JJ + delta]);			//z
					}
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
			}
			else		// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cenario);
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
						vars[i + (3+flag7)*I + delta].set(GRB_DoubleAttr_Obj, deltaT);	//F
					else
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
						{
							vars[i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT);	//pt
							vars[i + I + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT);		//u
						}
						else
						{
							vars[i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT);	//pt
							vars[i + I + delta].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT);		//u
						}
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						vars[i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT);		//cp
					else
						vars[i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida());		//up
				}
				if (flag1 > 0)		// considerar rede
				{
					for (size_t b = 0; b < B; b++)
						vars[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT);		//def
				}
				else
					vars[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT);		//def
				if ((t >= sistema_a->GetTt1()) && ((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
					if (sistema_a->GetFlagVfol() == true)
						for (int r = 0; r < R; r++)
							vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].set(GRB_DoubleAttr_Obj, (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT);		//vfol
			}
			delta += nt;
		}
		
		modeloGRB->update();

		// Termos quadráticos
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			delta = 0;
			GRBQuadExpr fo;
			//x = [pt u cp (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
			//x = [pt u up ud (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
			double * coeffs;
			coeffs = new double[n];
			for (int i = 0; i < n; i++)
				coeffs[i] = vars[i].get(GRB_DoubleAttr_Obj);
			fo.addTerms(coeffs, vars, n);
			for (int i = 0; i < n; i++)
				coeffs[i] = 0;
			for (int t = 0; t < T; t++)
			{
				if ((0 <= t) && (t) < sistema_a->GetTt1())		// nós do estágio 1
					deltaT = sistema_a->GetDeltaT1();
				else		// nós do estágio 2
					deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cenario);

				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
						coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;		//pt^2
					else
					{
						cout << "CriarFuncaoObjetivo2: Termo quadrático não implementado para essa modelagem!!" << endl;
						break;
						// a2*(pt^2 + 2*u*pt*pt_min + u^2*pt_min^2) = a2*(pt^2 + 2*u*pt*pt_min + u*pt_min^2)
					}
				}
				delta += nt;
			}
			fo.addTerms(coeffs, vars, vars, n);
			modeloGRB->setObjective(fo, GRB_MINIMIZE);
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during Scndec1SPC::CriarFuncaoObjetivoRL()" << endl;	
	}
}
void Scndec1SPC::FixarCondIniciais()
{
	int n_a = n - flag2*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int c = 0;
	int TUr = 0;
	if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
	{
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				// Verifica condições iniciais e fixa valores das variáveis caso necessário
				if (sistema_a->termeletricasVtr[i].GetU0() == 0)
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTDown() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				else
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTUp() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				if ((TUr >= 1) && (t + 1 <= TUr))
				{
					//u
					vars[i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					vars[i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//up
					vars[i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
					//ud
					vars[i + 3*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[i + 3*I + c].set(GRB_DoubleAttr_UB, 0);
				}
			}
		c += (n_a / T);
		}
	}
	else
	{
		for (int t = 0; t < T; t++)
		{
			for (size_t i = 0; i < I; i++)
			{
				// Verifica condições iniciais e fixa valores das variáveis caso necessário
				if (sistema_a->termeletricasVtr[i].GetU0() == 0)
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTDown() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				else
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTUp() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				if ((TUr >= 1) && (t + 1 <= TUr))
				{
					//u
					vars[i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					vars[i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//cp
					vars[i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
				}
			}
		c += (n_a / T);
		}
	}
}
int Scndec1SPC::ResolverProblemaRL(Scndec1Results * resultadoGurobi, const vetorfloat2 &lambda, const int iter)
{
	// Resolve o subproblema
	try
	{
		// vetor lambda é uma matriz, em q cada linha são os lambdas de cada cenário
		// o lambda de entrada aqui  é um vetor somente, referente ao cenário

		// coeficientes da RL
		vetorfloat coefRL;
		if (cenario == sistema_a->GetNCenarios() - 1)
		{
			coefRL.resize(lambda[cenario - 1].size());
			double sumL;
			for (size_t i = 0; i < coefRL.size(); i++)
			{
				sumL = 0;
				for (size_t j = 0; j < lambda.size(); j++)		// loop no numero de cenários
					sumL += lambda[j][i];
				coefRL[i] = sistema_a->hidreletricasVtr[0].GetProbAfluencia(cenario)*sumL;
			}
		}
		else
		{
			coefRL.resize(lambda[cenario].size());
			double sumL;
			for (size_t i = 0; i < coefRL.size(); i++)
			{
				sumL = 0;
				for (size_t j = 0; j < lambda.size(); j++)		// loop no numero de cenários
					sumL += lambda[j][i];
				coefRL[i] = - lambda[cenario][i] + sistema_a->hidreletricasVtr[0].GetProbAfluencia(cenario)*sumL;
			}
		}
		// Criar função objetivo, a seleçao dos coeficientes que dependem dos lambda de interesse é feita acima
		CriarFuncaoObjetivoRL(coefRL);
		
		modeloGRB->reset();

		//modeloGRB->write("SPC_subproblem.lp");	

		// Otimizar
		modeloGRB->optimize();

		// Salvar resultados
		x.clear();
		x.resize(n);
		nStatus = modeloGRB->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The subproblem " << cenario << " was not solved until optimality! " << nStatus << " Gap: " << modeloGRB->get(GRB_DoubleAttr_MIPGap) << endl;
		if ((nStatus == 2) || (nStatus == 9))
		{
			fo = double(modeloGRB->get(GRB_DoubleAttr_ObjVal));
			if (sistema_a->GetFlagVarBin())
				lb = double(modeloGRB->get(GRB_DoubleAttr_ObjBound));
			else
				lb = fo;
			if (sistema_a->GetFlagVarBin() == 1)
				sol_mipgap_curr = modeloGRB->get(GRB_DoubleAttr_MIPGap);
			else
				sol_mipgap_curr = 0;
			///cout << "Mipgap = " << sol_mipgap << endl;
			for (int i = 0; i < n; i++)
				x[i] = double(vars[i].get(GRB_DoubleAttr_X));
		}
		else
		{
			cout << "Inviabilidade no subproblema " << cenario << endl;
			modeloGRB->computeIIS();				
			modeloGRB->write("SPC_subproblem.lp");	
			modeloGRB->write("SPC_subproblem.ilp");
			for (int i = 0; i < n; i++)
			{	
				x[i] = 0;
			}
			fo = 0;
			lb = 0;
		}
		resultadoGurobi->GravarSolucao(fo, x, nStatus, cenario, lb);

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->GravarSolucao(fo, x, e.getErrorCode(), cenario, lb);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}
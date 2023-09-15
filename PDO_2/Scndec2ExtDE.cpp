#include "Scndec2ExtDE.h"

// ARRUMAR FUNÇÕES PARA DELTAt1 != DELTAt2 (todas que usam GetDeltaT)
CMatrizEsparsa Scndec2ExtDE::MatrizRestDemanda()	// Monta matriz da restrição de atendimento a demanda
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
CMatrizEsparsa Scndec2ExtDE::MatrizRestDemandaBarraUnica()	// Monta matriz da restrição de atendimento a demanda
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
CMatrizEsparsa Scndec2ExtDE::MatrizLimFluxo()
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
CMatrizEsparsa Scndec2ExtDE::MatrizLimPhgMin()
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
CMatrizEsparsa Scndec2ExtDE::MatrizLimPhgMax()
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
CMatrizEsparsa Scndec2ExtDE::MatrizBalHid()
{
	// T * R
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
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
CMatrizEsparsa Scndec2ExtDE::MatrizLimqMin()
{
	int n_a = n;
	int JJ = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
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
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ,l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec2ExtDE::MatrizLimqMax()
{
	int n_a = n;
	int JJ = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
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
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ,l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Scndec2ExtDE::MatrizTup()
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
CMatrizEsparsa Scndec2ExtDE::MatrizTdown()
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
CMatrizEsparsa Scndec2ExtDE::MatrizTupDown()
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
CMatrizEsparsa Scndec2ExtDE::MatrizRampaUp()
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
CMatrizEsparsa Scndec2ExtDE::MatrizRampaDown()
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
CMatrizEsparsa Scndec2ExtDE::MatrizLimPtMin()
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
		// referente à t = 0
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
			Alin.JuntarColuna(&a);
		}
		c += (n_a / T);
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
CMatrizEsparsa Scndec2ExtDE::MatrizLimPtMax()
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
CMatrizEsparsa Scndec2ExtDE::MatrizRestCP()
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
CMatrizEsparsa Scndec2ExtDE::MatrizVmeta()
{
	int n_a = n;
	size_t R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	n_a = n_a - flag2*sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->hidreletricasVtr.size());
	int c;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
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
CMatrizEsparsa Scndec2ExtDE::MatrizBalPotencia()
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
CMatrizEsparsa Scndec2ExtDE::MatrizBalVazao()
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
CMatrizEsparsa Scndec2ExtDE::MatrizFuncProd()
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
CMatrizEsparsa Scndec2ExtDE::MatrizPhMax()
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
CMatrizEsparsa Scndec2ExtDE::MatrizReserva()
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
CMatrizEsparsa Scndec2ExtDE::MatrizCortesF()
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
CMatrizEsparsa Scndec2ExtDE::MatrizRestNAnt()
{
	// x_w = sum_w(p_w * x_w)
	// para algumas variaveis de primeiro estagio
	// tamanho do vetor: no de n em n para cada cenario
	CMatrizEsparsa Alin(0, no);
	CMatrizEsparsa a(1, no);
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
	{
		//x = [pt u d phg] das restrições de não antecipatividade
		int deltax = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, i + deltax + n_cen*n, 1);	// x_w
				for (int j = 0; j < sistema_a->GetNCenarios(); j++)	// loop para todos elementos da restricao
					a.InserirElemento(0, i + deltax + j*n, - sistema_a->hidreletricasVtr[0].GetProbAfluencia(j));	// sum_w(p_w * x_w)
				Alin.JuntarColuna(&a);
			}
			deltax += sistema_a->termeletricasVtr.size();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, i + deltax + n_cen*n, 1);	// x_w
				for (int j = 0; j < sistema_a->GetNCenarios(); j++)	// loop para todos elementos da restricao
					a.InserirElemento(0, i + deltax + j*n, - sistema_a->hidreletricasVtr[0].GetProbAfluencia(j));	// sum_w(p_w * x_w)
				Alin.JuntarColuna(&a);
			}
			deltax += 2*sistema_a->termeletricasVtr.size();

			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				deltax += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				deltax += sistema_a->termeletricasVtr.size();
		
			if ( sistema_a->GetFlagModeloRede() == 1)
				deltax += sistema_a->barrasVtr.size() - 1;
			deltax += 2*sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, r + deltax + n_cen*n, 1);	// x_w
				for (int j = 0; j < sistema_a->GetNCenarios(); j++)	// loop para todos elementos da restricao
					a.InserirElemento(0, r + deltax + j*n, - sistema_a->hidreletricasVtr[0].GetProbAfluencia(j));	// sum_w(p_w * x_w)
				Alin.JuntarColuna(&a);
			}
			deltax += 2*sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, j + deltax + n_cen*n, 1);	// x_w
					for (int jj = 0; jj < sistema_a->GetNCenarios(); jj++)	// loop para todos elementos da restricao
						a.InserirElemento(0, j + deltax + jj*n, - sistema_a->hidreletricasVtr[0].GetProbAfluencia(jj));	// sum_w(p_w * x_w)
					Alin.JuntarColuna(&a);
				}
				deltax += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
				deltax += 2*sistema_a->hidreletricasVtr[r].GetNGrupos();

			if ( sistema_a->GetFlagModeloRede() > 0)
				deltax += sistema_a->barrasVtr.size();
			else
				deltax += 1;
		}
	}
	return Alin;
}
// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Scndec2ExtDE::LimRestDemanda(int &cenario)
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
			//if ((0 <= t) && (sistema_a->GetTt1() > t))
			//{
			//	for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
			//		D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);		// cen = 0, ou qq outro cenários, pois os T1 primeiros periodos são iguais (deterministicos)
			//}
			//else
			//{
			for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
				D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cenario, t);
			//}
		Lim[t*B + b][0] = 1;
		Lim[t*B + b][1] = D;
		}
	}
	return Lim;
}
vetorfloat2 Scndec2ExtDE::LimRestDemandaBarraUnica(int &cenario)
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
			for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
				D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cenario, t);
		}
		Lim[t][0] = 1;
		Lim[t][1] = D;
	}
	return Lim;
}
vetorfloat2 Scndec2ExtDE::LimFluxo0(int &cenario)
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
vetorfloat2 Scndec2ExtDE::LimFluxo2(int &cenario)
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
vetorfloat2 Scndec2ExtDE::LimPhgMin()
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
vetorfloat2 Scndec2ExtDE::LimPhgMax()
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
vetorfloat2 Scndec2ExtDE::LimBalHid(int &cenario)
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
vetorfloat2 Scndec2ExtDE::LimQMin()
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
vetorfloat2 Scndec2ExtDE::LimQMax()
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
vetorfloat2 Scndec2ExtDE::LimTup()
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
vetorfloat2 Scndec2ExtDE::LimTdown()
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
vetorfloat2 Scndec2ExtDE::LimTupDown()
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
vetorfloat2 Scndec2ExtDE::LimRampaUp()
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
vetorfloat2 Scndec2ExtDE::LimRampaDown()
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
vetorfloat2 Scndec2ExtDE::LimPtMin()
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
vetorfloat2 Scndec2ExtDE::LimPtMax()
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
vetorfloat2 Scndec2ExtDE::LimRestCP()
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
vetorfloat2 Scndec2ExtDE::LimVmeta()
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
vetorfloat2 Scndec2ExtDE::LimBalPotenciaL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Scndec2ExtDE::LimBalVazaoL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Scndec2ExtDE::LimFuncProdL()
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
vetorfloat2 Scndec2ExtDE::LimPhMax()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Scndec2ExtDE::LimReserva(int &cenario)
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
vetorfloat2 Scndec2ExtDE::LimCortesF()
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
vetorfloat2 Scndec2ExtDE::LimRestNAnt()
{
	vetorfloat2 Lim;
	vetorfloat L; L.resize(2);
	L[0] = 1;L[1] = 0;
	
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	//x = [pt u d phg] das restrições de não antecipatividade
	int n1 = sistema_a->GetTt1() * ( 2*sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size() + JJ);
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		for (int i = 0; i < n1; i++)		// loop para tds variaveis de primeiro estagio
			Lim.push_back(L);
	return Lim;
}
// Matriz dos coeficientes e limites das restrições lineares
// Caso seja necessário fazer modelagens diferentes para cada estágio é só criar uma condição que depende do T dentro da função
void Scndec2ExtDE::MatrizRestricoesLineares(vetorint * indexL, vetorint * indexC, vetorfloat * indexV, vetorint * nnZ, int * n_restricoes, int &cenario)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	//int count = 0;
	int count = nnZ->at(*n_restricoes);
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		M = MatrizRestDemanda();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		M = MatrizLimFluxo();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	else	// == 0 e == 3(nesse caso os limites de fluxo somente são adicionados se violados)
	{
		M = MatrizRestDemandaBarraUnica();
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
	M = MatrizBalHid();
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
		M = MatrizTupDown();		// inviabilidade aqui no cenario 0!!!
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
	M = MatrizReserva();
	SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
	*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		M = MatrizCortesF();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
	}
	rest_o = *n_restricoes;
	// Matrizes acima são geradas para cada cenário, duplicadas na funcao CriarRestricoes()
	// e ao final incluir restrição de não antecipatividade!!
	if (cenario == sistema_a->GetNCenarios() - 1)
	{
		M = MatrizRestNAnt();
		SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		rest_tot = *n_restricoes;
	}
}
void Scndec2ExtDE::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor, int &cenario)
{
	vetorfloat2 L;
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		L = LimRestDemanda(cenario);
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo0(cenario);
		AlocarLimites(&L, LimTipo, LimValor);
		L = LimFluxo2(cenario);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	else	// == 0 e == 3(nesse caso os limites de fluxo somente são adicionados se violados)
	{
		L = LimRestDemandaBarraUnica(cenario);
		AlocarLimites(&L, LimTipo, LimValor);
		if ( sistema_a->GetFlagModeloRede() == 2 )
		{
			L = LimFluxo0(cenario);
			AlocarLimites(&L, LimTipo, LimValor);
			L = LimFluxo2(cenario);
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
	L = LimBalHid(cenario);
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
	L = LimPtMax();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPtMin();
	AlocarLimites(&L, LimTipo, LimValor);
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		L = LimRestCP();
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimVmeta();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalPotenciaL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimBalVazaoL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimFuncProdL();
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimReserva(cenario);
	AlocarLimites(&L, LimTipo, LimValor);
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		L = LimCortesF();
		AlocarLimites(&L, LimTipo, LimValor);
	}
	if (cenario == sistema_a->GetNCenarios() - 1)
	{
		L = LimRestNAnt();
		AlocarLimites(&L, LimTipo, LimValor);
	}
}

Scndec2ExtDE::Scndec2ExtDE(CSistema * const sistema_end) :ambGRB(GRBEnv()),modeloGRB(GRBModel(ambGRB))
{
	sistema_a = sistema_end;
	
	T = sistema_a->GetTt2();
	flag1 = int (sistema_a->GetFlagModeloRede());	
	flag2 = int (sistema_a->GetFlagVfol());
	flag3 = int (sistema_a->GetFlagPhmax());
	if (flag3 == 1)
		cout << "Atencao: phmax = 1 não implementado para esta decomposição!!!" << endl;
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size();
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= T * (sistema_a->barrasVtr.size() - 1);
	
	no = n * sistema_a->GetNCenarios();
	//x.resize(no);

	// Criar variáveis
	CriarVariaveis();
	modeloGRB.update();	// Atualiza o modelo Gurobi.
	vars = modeloGRB.getVars();
	
	// Adicionar restrições
	CriarRestricoes();
	modeloGRB.update();	// Atualiza o modelo Gurobi.
	constr = modeloGRB.getConstrs();
	
	// Fixar condições iniciais
	FixarCondIniciais();
	modeloGRB.update();
}
Scndec2ExtDE::~Scndec2ExtDE(void)
{
	delete vars;
	delete constr;
}

void Scndec2ExtDE::CriarVariaveis()
{
	try 
	{
		size_t I = sistema_a->termeletricasVtr.size();
		size_t R = sistema_a->hidreletricasVtr.size();
		size_t B = sistema_a->barrasVtr.size() - 1;
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			// Estagio 1 e 2
			for (int t = 0; t < sistema_a->GetTt2(); t++)
			{
				//x = [pt u cp F teta ph v d s (phmax) phg q z def vfol]
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				{
					for (size_t i = 0; i < I; i++)	//pt
						modeloGRB.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)	//pt
						modeloGRB.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagVarBin() == true)	//u (binario ou continuo)
				{
					for (size_t i = 0; i < I; i++)	//u
						modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)	//u
						modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				{
					if (sistema_a->GetFlagVarBin() == true)
					{
						for (size_t i = 0; i < I; i++)		//up
							modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
						for (size_t i = 0; i < I; i++)		//ud
							modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
					}
					else
					{
						for (size_t i = 0; i < I; i++)		//up
							modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
						for (size_t i = 0; i < I; i++)		//ud
							modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					}
				}
				else
				{
					for (size_t i = 0; i < I; i++)	//cp
						modeloGRB.addVar(0, sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				{
					for (size_t i = 0; i < I; i++)
						modeloGRB.addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
				}
				if ( sistema_a->GetFlagModeloRede() == 1)
					for (size_t b = 0; b < B; b++)	//teta
						modeloGRB.addVar(-3.14159, 3.14159, 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modeloGRB.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					modeloGRB.addVar(double(sistema_a->hidreletricasVtr[r].GetVmin()), double(sistema_a->hidreletricasVtr[r].GetVmax()), 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
				{
					double qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
					modeloGRB.addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
				{
					modeloGRB.addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < R; r++)	//phg
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						modeloGRB.addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				for (size_t r = 0; r < R; r++)	//q
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						modeloGRB.addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
				{
					for (size_t r = 0; r < R; r++)	//z
					{
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
					}
				}
				else
				{
					for (size_t r = 0; r < R; r++)	//z
					{
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					}
				}
				if ( sistema_a->GetFlagModeloRede() > 0)
				{
					for (size_t b = 0; b < B + 1; b++)	//def
					{
						double cap_d = 0;
						for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
							cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
						modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
					}
				}
				else
				{
					double cap_d = 0;			//def barra unica
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
							cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t);
					modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( sistema_a->GetFlagVfol() == true )
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//vfol
					modeloGRB.addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
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
void Scndec2ExtDE::CriarRestricoes()
{
	try
	{
		int n_restricoes = 0;
		vetorint indexL, indexC, nnZ, LimTipo;
		vetorfloat indexV, LimValor;
		indexL.resize(0);indexC.resize(0);indexV.resize(0);nnZ.resize(0);LimTipo.resize(0);LimValor.resize(0);
		nnZ.push_back(0);
		int tam_ant = 0;
		int delta_tam = 0;
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			MatrizRestricoesLineares(&indexL, &indexC, &indexV, &nnZ ,&n_restricoes, cen);
			// ajustar colunas
			if (cen == 0)
				delta_tam = indexC.size() - tam_ant;
			else
				for (size_t i = tam_ant; i < tam_ant + delta_tam; i++)		// poderia ser ateh indexC.size(), mas no ultimo cenario tem-se as colunas das restrições de n antecipatividade q n devem ser deslocadas
					indexC[i] += n * cen;
			tam_ant = indexC.size();
			MatrizLimitesLineares(&LimTipo, &LimValor, cen);
		}
		GRBLinExpr restricao;
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
			switch (LimTipo[i]) {
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
void Scndec2ExtDE::CriarFuncaoObjetivo()
{
	try 
	{
		int R = sistema_a->hidreletricasVtr.size();
		int I = sistema_a->termeletricasVtr.size();
		int B = sistema_a->barrasVtr.size();
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int nt = ( n - flag2*sistema_a->hidreletricasVtr.size() ) / T;
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

		// todos os cenarios... com probabilidades para os termos do primeiro estagio tb
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		{
			delta = 0;
			// Termos lineares
			for (int t = 0; t < T; t++)
			{
				if ((0 <= t) && (t < sistema_a->GetTt1()))		// nós do estágio 1
					deltaT = sistema_a->GetDeltaT1() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen);
				else		// nós do estágio 2
					deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen);
			
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
						vars[i + (3+flag7)*I + delta + n_cen*n].set(GRB_DoubleAttr_Obj, deltaT);	//F
					else
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
						{
							vars[i + delta + n_cen*n].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT);	//pt
							vars[i + I + delta + n_cen*n].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT);		//u
						}
						else
						{
							vars[i + delta + n_cen*n].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT);	//pt
							vars[i + I + delta + n_cen*n].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT);		//u
						}
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						vars[i + 2*I + delta + n_cen*n].set(GRB_DoubleAttr_Obj, deltaT);		//cp
					else
						vars[i + 2*I + delta + n_cen*n].set(GRB_DoubleAttr_Obj, deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida());		//up
				}
				if (sistema_a->GetFlagModeloRede() > 0)		// considerar rede
				{
					for (size_t b = 0; b < B; b++)
						vars[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta + n_cen*n].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT);		//def
				}
				else
					vars[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta + n_cen*n].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT);		//def
				if ((t >= sistema_a->GetTt1()) && ((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
					if (sistema_a->GetFlagVfol() == true)
						for (int r = 0; r < R; r++)
							vars[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta + n_cen*n].set(GRB_DoubleAttr_Obj, (20/(0.0036*sistema_a->GetDeltaT2()))*sistema_a->GetCustoDeficit()*deltaT);		//vfol
				delta += nt;
			}
		}
		modeloGRB.update();
		
		// Termos quadráticos
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			GRBQuadExpr fo;
			double * coeffs;
			coeffs = new double[n];
			for (int i = 0; i < n; i++)
				coeffs[i] = vars[i].get(GRB_DoubleAttr_Obj);
			fo.addTerms(coeffs, vars, n);
			for (int i = 0; i < n; i++)
				coeffs[i] = 0;
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			{
				delta = 0;
				for (int t = 0; t < T; t++)
				{
					if ((0 <= t) && (t) < sistema_a->GetTt1())		// nós do estágio 1
						deltaT = sistema_a->GetDeltaT1() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen);
					else		// nós do estágio 2
						deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen);

					for (int i = 0; i < I; i++)
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
							coeffs[i + delta + n_cen*n] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;		//pt^2
						else
						{
							cout << "CriarFuncaoObjetivo2: Termo quadrático não implementado para essa modelagem!!" << endl;
							break;
							// a2*(pt^2 + 2*u*pt*pt_min + u^2*pt_min^2) = a2*(pt^2 + 2*u*pt*pt_min + u*pt_min^2)
						}
					}
					delta += nt;
				}
			}
			fo.addTerms(coeffs, vars, vars, n);
			modeloGRB.setObjective(fo, GRB_MINIMIZE);
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
void Scndec2ExtDE::FixarCondIniciais()
{
	int n_a = n - flag2*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	{
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
						vars[i + I + c + n_cen*n].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
						vars[i + I + c + n_cen*n].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
						//up
						vars[i + 2*I + c + n_cen*n].set(GRB_DoubleAttr_LB, 0);
						vars[i + 2*I + c + n_cen*n].set(GRB_DoubleAttr_UB, 0);
						//ud
						vars[i + 3*I + c + n_cen*n].set(GRB_DoubleAttr_LB, 0);
						vars[i + 3*I + c + n_cen*n].set(GRB_DoubleAttr_UB, 0);
					}
				}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
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
						vars[i + I + c + n_cen*n].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
						vars[i + I + c + n_cen*n].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
						//cp
						vars[i + 2*I + c + n_cen*n].set(GRB_DoubleAttr_LB, 0);
						vars[i + 2*I + c + n_cen*n].set(GRB_DoubleAttr_UB, 0);
					}
				}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
			}
		}
	}
}
double Scndec2ExtDE::ResolverProblema(LMRow lambda)
{
	try
	{
		// Criar função objetivo
		CriarFuncaoObjetivo();
		modeloGRB.update();	// Atualiza o modelo Gurobi.

		//modeloGRB.write("probED_ext.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modeloGRB.getEnv().set(GRB_IntParam_OutputFlag, 0);

		//modeloGRB.getEnv().set(GRB_IntParam_LogToConsole, 0);				// Nao escrever detalhes na tela
		//modeloGRB.getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		//modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.0001);			// Define o gap de tolerancia
		//modeloGRB.getEnv().set(GRB_IntParam_ScaleFlag, 0);				// Desabilita o model scaling
		//modeloGRB.getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modeloGRB.getEnv().set(GRB_IntParam_Threads, 4);
		//modeloGRB.getEnv().set(GRB_IntParam_Presolve, 0);					//Desabilitar o presolve
		//modeloGRB.getEnv().set(GRB_IntParam_PreSparsify, 1);
		
		
		modeloGRB.getEnv().set(GRB_IntParam_Method, 2);
		// (-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent)
		// o método de barreira para este problema resulta em multiplicadores muito ruins! (errado para a ultima restrição?)
		
		modeloGRB.getEnv().set(GRB_IntParam_Crossover, 0);
		//modeloGRB.getEnv().set(GRB_IntParam_MIPFocus, 2);
		//modeloGRB.getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

		modeloGRB.getEnv().set(GRB_DoubleParam_TimeLimit, 1800);		// Limita o tempo de resolução do problema
		// fazer esse tempo como uma razão do tempo máximo do Bundle (ou melhor do tamanho do problema)!!
		// ou melhor, como essa etapa é essencial para o bom desempenho do bundle, deixar tempo ilimitado (= tempo máximo do bundle)
		
		//modeloGRB.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9);
		//modeloGRB.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);	// Define tolerancia de otimalidade
		//modeloGRB.reset();

		// Otimizar
		modeloGRB.optimize();

		//// ??? imprimir a solução do modelo estendido para comparar com a solução do ED continuo sem as restrições de n antecipatividade
		////
		//vetorfloat x;
		//for (int i = 0; i < no; i++)
		//	x.push_back(vars[i].get(GRB_DoubleAttr_X));
		//ofstream * inFile;
		//inFile = new ofstream( "x.txt", ios::out );
		//if ( inFile->is_open() )
		//{	
		//	for (int i = 0; i < no; i++)
		//		*inFile << std::scientific << setprecision(10) << x[i] << endl;
		//}
		//else
		//	cout << "Unable to open file";
		//inFile->close();
		////

		// Salvar resultados
		int nStatus = modeloGRB.get(GRB_IntAttr_Status);
		if ((nStatus == 2) || (nStatus == 9) || (nStatus == 13))
		{
			cout << setprecision(15) << modeloGRB.get(GRB_IntAttr_Status) << " : " << modeloGRB.get(GRB_DoubleAttr_ObjVal) << " \n";
			fo = double(modeloGRB.get(GRB_DoubleAttr_ObjVal));
			for (int i = rest_o; i < rest_tot; i++)		// Somente valores dos L das restriçoes auxiliares
				lambda[i - rest_o] = double(constr[i].get(GRB_DoubleAttr_Pi));

			//ofstream * inFile;
			//inFile = new ofstream( "LM.txt", ios::out );
			//if ( inFile->is_open() )
			//{
			//	for (int i = rest_o; i < rest_tot; i++)		// Somente valores dos L das restriçoes auxiliares
			//	{	
			//		lambda[i - rest_o] = double(constr[i].get(GRB_DoubleAttr_Pi));
			//		*inFile << lambda[i - rest_o] << endl;
			//	}
			//	inFile->close();
			//}
			//else
			//	cout << "Unable to open file";

			//ifstream inFile( "LM_LR.txt", ios::in );
			//if ( !inFile )                                                            
			//{                                                                               
			//	cout << "File "<< "LM_LR.txt" << " could not be opened" << endl;
			//	exit( 1 );
			//}
			//int cc = 0;
			//while ( ! inFile.eof() )
			//{
			//	inFile >> lambda[cc++];
			//}
			//inFile.close();

		}
		else
		{
			cout << "Inviabilidade no problema extendido " << endl;
			modeloGRB.computeIIS();				
			modeloGRB.write("probED_ext.lp");	
			modeloGRB.write("probED_ext.ilp");
			for (int i = rest_o; i < rest_tot; i++)
			{	
				lambda[i - rest_o] = 0;
				//lambda->push_back(0);
			}
			fo = 0;
		}
		//cout << setprecision(12) << endl;

		//ofstream * inFile;
		//inFile = new ofstream( "Lambda.txt", ios::out );
		//if ( inFile->is_open() )
		//{	
		//	*inFile << fo << endl;
		//	for (int i = rest_o; i < rest_tot; i++)
		//		*inFile << lambda[i - rest_o] << endl;
		//}
		//else
		//	cout << "Unable to open file";
		//inFile->close();
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}

return ( fo );
}
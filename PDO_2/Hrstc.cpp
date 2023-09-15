/*--------------------------------------------------------------------------*/
/*---------------------------- File Hrstc.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para aplicar heuristicas na solução da RL
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Hrstc.h"

// Funcoes para criar a matriz de restriçoes

// Algumas dessas funções são diferentes do ED, pois são "adaptadas" para a divisão por períodos de tempo

CMatrizEsparsa Hrstc::MatrizRestDemanda()	// Monta matriz da restrição de atendimento a demanda
{
	// T * B
	int n_a = n;
	CMatrizEsparsa Agh(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
			Agh.InserirElemento(b, sistema_a->barrasVtr[b].hidrosPtr[br]->GetIdentUsina(), 1);
	//	for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
	//		for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
	//			if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
	//				Agh.InserirElemento(b, r, 1);
	CMatrizEsparsa Agt(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
			Agt.InserirElemento(b, sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina(), 1);
		//for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
		//	for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
		//		if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
		//			Agt.InserirElemento(b, i, 1);
	CMatrizEsparsa Agtu(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
	if (sistema_a->GetFlagTbinaryModel() == 1)
	{
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
				Agtu.InserirElemento(b, sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina(), sistema_a->termeletricasVtr[sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina()].GetPmin());
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
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizRestDemandaBarraUnica()	// Monta matriz da restrição de atendimento a demanda
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
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizLimFluxo()
{
	int n_a = n;
	CMatrizEsparsa Alin(T * sistema_a->linhasVtr.size(), n_a);
	if ( flag1 == 1 )
	{
		CMatrizEsparsa Alb(sistema_a->linhasVtr.size(),sistema_a->barrasVtr.size());
		for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		{
			Alb.InserirElemento(l, sistema_a->linhasVtr[l].de_barra->GetIdentBarra(), 1);
			Alb.InserirElemento(l, sistema_a->linhasVtr[l].para_barra->GetIdentBarra(), -1);
		}
		Alb.RemoverColuna(sistema_a->GetBarraRef() - 1);
		CMatrizEsparsa TT(sistema_a->linhasVtr.size(), sistema_a->linhasVtr.size());
		for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
			TT.InserirElemento(l, l, 100 / (sistema_a->linhasVtr[l].GetReatancia()));
		n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		int ll = 0;
		int c = 0;
		TT.MultiplicarPorMatriz(&Alb);
		for ( int t = 0; t < T; t++)
		{
			Alin.InserirMatriz(ll, c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),ll + sistema_a->linhasVtr.size() - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
			ll = ll + sistema_a->linhasVtr.size();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
	}
	else // flag1 == 2
	{
		MatrixXd Agh = MatrixXd::Zero(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
				Agh(b, sistema_a->barrasVtr[b].hidrosPtr[br]->GetIdentUsina()) = 1;
		MatrixXd Agt = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
				Agt(b, sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina()) = 1;
		MatrixXd Agtu = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
		if (sistema_a->GetFlagTbinaryModel() == 1)
		{
			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
					Agtu(b, sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina()) = sistema_a->termeletricasVtr[sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina()].GetPmin();
		}
		MatrixXd Adef = MatrixXd::Identity(sistema_a->barrasVtr.size(),sistema_a->barrasVtr.size());
		// Multiplicar - Beta pelas matrizes acima para cada conjunto de variáveis e incluir em Alin
		Agh = - sistema_a->Beta * Agh; //Agh = - Beta * Agh;
		Agt = - sistema_a->Beta * Agt; //Agt = - Beta * Agt;
		Agtu = - sistema_a->Beta * Agtu; //Agtu = - Beta * Agtu;
		Adef = - sistema_a->Beta * Adef; //Adef = - Beta * Adef;

		//
		n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizLimPhgMin()
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizLimPhgMax()
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizBalHid()
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
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[r].GetUsinaJusante() != 0) && (tempo_viagem < sistema_a->GetTt2()))
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
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[r].GetUsinaJusante() != 0) && (tempo_viagem < sistema_a->GetTt2()))
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
	CMatrizEsparsa Alin(R * T,n_a);
	CMatrizEsparsa AvNeg(R);
	AvNeg.MultiplicarPorEscalar( -1);
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	int l = 0;
	int c = 0;
	int cc = 0;
	int flag1a = 0;
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
		//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		//	if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
		//	{
		//		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &Av, 0, 0);
		//		cc = 0;
		//		for (int tt = 0; tt <= t; tt++)
		//		{
		//			Alin.InserirMatriz(l,cc + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R),l + R - 1,cc + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 3*R) - 1, &Ad, t*R, tt*R);
		//			if (((tt + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((tt + 1 - sistema_a->GetTt1()) > 0))
		//				cc += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		//			else
		//				cc += (n_a / T);
		//		}
		//		if (t > 0)
		//			if (t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) == sistema_a->GetTt1())
		//				Alin.InserirMatriz(l,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,(n_a / T)*sistema_a->GetTt1() - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
		//			else
		//				Alin.InserirMatriz(l,c - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R),l + R - 1,c - (n_a / T) + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R) - 1, &AvNeg, 0, 0);
		//	}
		l = l + R;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizLimqMin()
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ,l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizLimqMax()
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < T; t++)
	{
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + JJ,l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 2*JJ - 1, &A1, 0, 0);
		Alin.InserirMatriz(l,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 2*JJ),l + JJ - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size())) + 3*JJ - 1, &A2, 0, 0);
		l = l + JJ;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c = c + (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizTup()
{
	// Tup é o tempo que a unidade deve ficar ligada após ser ligada, portanto para ficar ligada por 3 periodos, Tup = 2;
	if (flag7 == 0)
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
	else
	{
		// T * I (para usinas com Tup > 0)	nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < Tup + 2
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

				if (Tup > 0)		// só adiciona restrição para Tup > 0
				{
					a.RemoverTodosElementos();
					// Na Heuristica deve-se adicionar essa restrição desde o inicio para não causar inviabilidades no problema!! Uma alternatica seria usar cortes de viabilidade...
					// These being constraints looking back in time, they have to be adapted when you are "at the beginning of time". But "adapted" does not mean "removed altogether".
					a.InserirElemento(0, c + i, -1);		// adiciona termo de u
					for (int tt = 0; tt <= min(t_a,Tup); tt++)
					{
						if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
							a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + i + I - tt*(n_a / T), 1);
					}
					Alin.JuntarColuna(&a);
				}
				//// Código antigo:
				//{
				//	a.RemoverTodosElementos();
				//	//if  (t >= Tup + 1)
				//	if  (t_a >= Tup)
				//	{
				//		a.InserirElemento(0, c + i, -1);		// adiciona termo de u
				//		for (int tt = 0; tt <= Tup; tt++)
				//		{
				//			if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
				//				a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);
				//			else
				//				a.InserirElemento(0, ca + i + I - tt*(n_a / T), 1);
				//		}
				//	}
				//	else if ((Tup >= sistema_a->GetTt2()) && (t_a == sistema_a->GetTt2() - 1))
				//	{
				//		a.InserirElemento(0, c + i, -1);		// adiciona termo de u
				//		for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
				//		{
				//			if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
				//				a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);
				//			else
				//				a.InserirElemento(0, ca + i + I - tt*(n_a / T), 1);
				//		}
				//	}
				//	Alin.JuntarColuna(&a);
				//}
				else {}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Hrstc::MatrizTdown()
{
	// O valor de Tdown tem significado similar ao Tup
	if (flag7 == 0)
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
	else
	{
		// T * I (para usinas com Tdown > 0)	nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < Tdown + 2
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
					a.RemoverTodosElementos();
					// Na Heuristica deve-se adicionar essa restrição desde o inicio para não causar inviabilidades no problema!! Uma alternatica seria usar cortes de viabilidade...
					// These being constraints looking back in time, they have to be adapted when you are "at the beginning of time". But "adapted" does not mean "removed altogether".
					a.InserirElemento(0, c + i, 1);
					for (int tt = 0; tt <= min(t_a,Tdown); tt++)
					{
						if (t_a - tt >= sistema_a->GetTt1())
							a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + i + 2*I - tt*(n_a / T), 1);
					}
					// Código antigo:
					//if  (t_a >= Tdown)
					//{
					//	a.InserirElemento(0, c + i, 1);
					//	for (int tt = 0; tt <= Tdown; tt++)
					//	{
					//		if (t_a - tt >= sistema_a->GetTt1())
					//			a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
					//		else
					//			a.InserirElemento(0, ca + i + 2*I - tt*(n_a / T), 1);
					//	}
					//}
					//else if ((Tdown >= sistema_a->GetTt2()) && (t_a == sistema_a->GetTt2() - 1))
					//{
					//	a.InserirElemento(0, c + i, 1);
					//	for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
					//	{
					//		if (t_a - tt >= sistema_a->GetTt1())
					//			a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
					//		else
					//			a.InserirElemento(0, ca + i + 2*I - tt*(n_a / T), 1);
					//	}
					//}
					Alin.JuntarColuna(&a);
				}
				else {}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Hrstc::MatrizTupDown()
{
	// T * I	(nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < 2 ou t < Tup + 2 ou t < Tdown + 2
	int n_a = n;
	size_t I = sistema_a->termeletricasVtr.size();
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int cen, t_a;
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
			a.RemoverTodosElementos();
			if (t >= 1)
			{
				a.InserirElemento(0, c + i, 1);
				a.InserirElemento(0, c + i + I, -1);
				a.InserirElemento(0, c + i + 2*I, 1);
				if (t_a - 1 >= sistema_a->GetTt1())
					a.InserirElemento(0, c + i - 1*(n_a / T), -1);
				else
					a.InserirElemento(0, ca + i - 1*(n_a / T), -1);
			}
			else
			{
				a.InserirElemento(0, c + i, 1);
				a.InserirElemento(0, c + i + I, -1);
				a.InserirElemento(0, c + i + 2*I, 1);
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
CMatrizEsparsa Hrstc::MatrizRampaUp()
{
	// T * I
	if (flag7 == 0)
	{
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
	else
	{
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
						a.InserirElemento(0, c + i - (n_a / T), -1);
					else
						a.InserirElemento(0, ca + i - (n_a / T), -1);
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
}
CMatrizEsparsa Hrstc::MatrizRampaDown()
{
	// T * I
	if (flag7 == 0)
	{
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
	else
	{
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
}
CMatrizEsparsa Hrstc::MatrizLimPtMin()
{
	if (flag7 == 0)
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
	else
	{
		// T * I
		int n_a = n;
		n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		size_t I = sistema_a->termeletricasVtr.size();
		CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
		CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
		int ca, t_a, cen;
		int c = 0;
		c += (n_a / T);						// referente à t = 0
		// t = 0
		for (size_t i = 0; i < I; i++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
			Alin.JuntarColuna(&a);
		}
		for (int t = 1; t < T; t++)			// nó inicial, já adicionadas restrições, por isso t começa em 1
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
				if (sistema_a->termeletricasVtr[i].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
					if (t_a - 1 >= sistema_a->GetTt1())		// todos periodos menos do inicio do segundo estágio (incluindo o de final de horizonte)
					{
						a.InserirElemento(0, c + i - (n_a / T), 1);		// pt-1
						a.InserirElemento(0, c + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u-1
					}
					else									// periodos no inicio do segundo estágio
					{
						a.InserirElemento(0, ca + i - (n_a / T), 1);
						a.InserirElemento(0, ca + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					}
					Alin.JuntarColuna(&a);
				}
				else
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
					if (t_a - 1 >= sistema_a->GetTt1())		// todos periodos menos do inicio do segundo estágio (incluindo o de final de horizonte)
					{
						a.InserirElemento(0, c + i - (n_a / T), 1);		// pt-1
						a.InserirElemento(0, c + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u-1
						a.InserirElemento(0, c + i - (n_a / T) + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// up-1
					}
					else									// periodos no inicio do segundo estágio
					{
						a.InserirElemento(0, ca + i - (n_a / T), 1);
						a.InserirElemento(0, ca + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
						a.InserirElemento(0, ca + i - (n_a / T) + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					}
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
}
CMatrizEsparsa Hrstc::MatrizLimPtMax()
{
	// T * I
	if (flag7 == 0)
	{
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
	else
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
					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
					Alin.JuntarColuna(&a);
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Hrstc::MatrizRestCP()
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
	// I * T
	//int n_a = n;
	//n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	//size_t I = sistema_a->termeletricasVtr.size();
	//CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	//CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	//CMatrizEsparsa Alin3(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	//int c;
	//for (size_t i = 0; i < I; i++)
	//{	
	//	Alin.ZerarMatriz(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	//	c = i + I;
	//	for (int t = 0; t < T; t++)
	//	{
	//		a.RemoverTodosElementos();
	//		if (t == 0)
	//		{
	//			a.InserirElemento(0, c + I, 1);
	//			a.InserirElemento(0, c, - sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
	//			//a.InserirElemento(0, c, - 1);
	//		}
	//		else
	//		{
	//			a.InserirElemento(0, c + I, 1);
	//			a.InserirElemento(0, c, - sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
	//			a.InserirElemento(0, c - (n_a / T), sistema_a->termeletricasVtr[i].GetCoefCustoPartida());
	//			//a.InserirElemento(0, c, - 1);
	//			//a.InserirElemento(0, c - (n_a / T), 1);
	//		}
	//		Alin.JuntarColuna(&a);
	//		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
	//			c = c + (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
	//		else
	//			c = c + (n_a / T);
	//	}
	//	CMatrizEsparsa Alin2((sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())), n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	//	CMatrizEsparsa A1(sistema_a->GetTt2() - sistema_a->GetTt1(), (n_a / T)*sistema_a->GetTt1());
	//	CMatrizEsparsa A2((sistema_a->GetTt2() - sistema_a->GetTt1()), (n_a / T)*sistema_a->GetTt2() - (n_a / T)*sistema_a->GetTt1());
	//	Alin2.InserirMatriz(0, 0, sistema_a->GetTt1() - 1, (n_a / T)*sistema_a->GetTt1() - 1, &Alin, 0, 0);
	//	A1.InserirMatriz(0, 0, sistema_a->GetTt2() - sistema_a->GetTt1() - 1, (n_a / T)*sistema_a->GetTt1() - 1, &Alin, sistema_a->GetTt1(), 0);
	//	A2.InserirMatriz(0, 0, sistema_a->GetTt2() - sistema_a->GetTt1() - 1, (n_a / T)*sistema_a->GetTt2() - (n_a / T)*sistema_a->GetTt1() - 1, &Alin, sistema_a->GetTt1(), (n_a / T)*sistema_a->GetTt1());
	//	int linha = 0;
	//	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	//	{
	//		linha = sistema_a->GetTt1() + (cen)*(sistema_a->GetTt2() - sistema_a->GetTt1());
	//		Alin2.InserirMatriz(linha, 0, linha + (sistema_a->GetTt2() - sistema_a->GetTt1()) - 1, (n_a / T)*sistema_a->GetTt1() - 1, &A1, 0, 0);
	//		Alin2.InserirMatriz(linha, flag2*cen*sistema_a->hidreletricasVtr.size() + (n_a / T)*(sistema_a->GetTt1()+cen*(sistema_a->GetTt2() - sistema_a->GetTt1())), linha + (sistema_a->GetTt2() - sistema_a->GetTt1()) - 1, flag2*cen*sistema_a->hidreletricasVtr.size() + (n_a / T)*(sistema_a->GetTt2() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1())) - 1, &A2, 0, 0);
	//	}
	//	Alin3.JuntarColuna(&Alin2);
	//}
	//return Alin3;
}
CMatrizEsparsa Hrstc::MatrizVmeta()
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	c = (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R + (sistema_a->GetTt2() - 1)*(n_a / T);
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		for (size_t r = 0; r < R; r++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c + r, 1);
			if (sistema_a->GetFlagVfol() == true)
				a.InserirElemento(0, c + r + (3+flag3)*R + 3*JJ + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d), 1);
			Alin.JuntarColuna(&a);
		}
		c += (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T) + flag2*R;
	}
	return Alin;
}
CMatrizEsparsa Hrstc::MatrizBalPotencia()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B;
	int jj = dd + (4+flag3)*R;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
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
CMatrizEsparsa Hrstc::MatrizBalVazao()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B + 2*R;
	int jj = dd + JJ + (2+flag3)*R;
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
CMatrizEsparsa Hrstc::MatrizFuncProd()
{
	// R * Tapp * T
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	size_t I = sistema_a->termeletricasVtr.size();
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int vv = (3+flag4+flag7)*I + B + R;
	int jj = vv + (3+flag3)*R + JJ;
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
CMatrizEsparsa Hrstc::MatrizPhMax()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B + 4*R;
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
CMatrizEsparsa Hrstc::MatrizReserva()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t R = sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	size_t B = flag1a*(sistema_a->barrasVtr.size() - 1);
	int nt;
	int JJ = 0;
	for (size_t i = 0; i < R; i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	CMatrizEsparsa Alin(0,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	CMatrizEsparsa a(1,n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int dd = (3+flag4+flag7)*I + B;
	int jj = dd + 4*R + 2*JJ;
	int cc = 0;
	for (int t = 0; t < T; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
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
CMatrizEsparsa Hrstc::MatrizCortesF()
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
				if ( flag7 == 0 )
					a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetCoefA0(n_cort));
				else
					a.InserirElemento(0, c + i + I, - sistema_a->termeletricasVtr[i].GetCoefA0(n_cort) - sistema_a->termeletricasVtr[i].GetCoefA1(n_cort) * sistema_a->termeletricasVtr[i].GetPmin());
				a.InserirElemento(0, c + i + (3+flag7)*I, 1);
				Alin.JuntarColuna(&a);
			}
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
	// I * T
	//int n_a = n;
	//n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	//size_t I = sistema_a->termeletricasVtr.size();
	//CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	//CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
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
vetorfloat2 Hrstc::LimRestDemanda()
{
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
vetorfloat2 Hrstc::LimRestDemandaBarraUnica()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, T, 2);
	IniciaMatriz(&Lim, 0);
	int cen = 0;
	for (int t = 0; t < T; t++)
	{
		double D = 0;
		for (size_t b = 0; b < B; b++)
		{
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
			}
		}
		Lim[t][0] = 1;
		Lim[t][1] = D;
	}
	return Lim;
}
vetorfloat2 Hrstc::LimFluxo0()
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
		int cen = 0;
		VectorXd Ad(sistema_a->barrasVtr.size());
		for (int t = 0; t < T; t++)
		{
			Ad.resize(sistema_a->barrasVtr.size());
			if (sistema_a->GetTt1() <= t)
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			for (size_t b = 0; b < B; b++)
			{
				double D = 0;
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
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
vetorfloat2 Hrstc::LimFluxo2()
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
		int cen = 0;
		VectorXd Ad(sistema_a->barrasVtr.size());
		for (int t = 0; t < T; t++)
		{
			Ad.resize(sistema_a->barrasVtr.size());
			if (sistema_a->GetTt1() <= t)
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			for (size_t b = 0; b < B; b++)
			{
				double D = 0;
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
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
vetorfloat2 Hrstc::LimPhgMin()
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
vetorfloat2 Hrstc::LimPhgMax()
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
vetorfloat2 Hrstc::LimBalHid()
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
vetorfloat2 Hrstc::LimQMin()
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
				//Lim(jj + j,1) = double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				Lim[jj + j][0] = 2;
			jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	return Lim;
}
vetorfloat2 Hrstc::LimQMax()
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
vetorfloat2 Hrstc::LimTup()
{
	if (flag7 == 0)
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
vetorfloat2 Hrstc::LimTdown()
{
	if (flag7 == 0)
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
vetorfloat2 Hrstc::LimTupDown()
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
vetorfloat2 Hrstc::LimRampaUp()
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
vetorfloat2 Hrstc::LimRampaDown()
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
vetorfloat2 Hrstc::LimPtMin()
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
		// (T - 1 + n_cen) * I (porém em t = 0 não são adicionadas restrições e em t = T2 são adicionadas 2 por usina)
		int I = sistema_a->termeletricasVtr.size();
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0; L[1] = 0;
		// t = 0
		for (size_t i = 0; i < I; i++)
		{
			if (sistema_a->termeletricasVtr[i].GetU0() == 1)
				if (sistema_a->termeletricasVtr[i].GetX0() > 1)  // indica que a usina foi ligada antes do periodo t = -1, portanto up-1 = 0;
					L[1] = sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPt0();
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
vetorfloat2 Hrstc::LimPtMax()
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
vetorfloat2 Hrstc::LimRestCP()
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
	// I * T
	//int n_a = n;
	//int I = sistema_a->termeletricasVtr.size();
	//vetorfloat2 Lim;
	//DimensionarMatriz(&Lim, I * T, 2);
	//IniciaMatriz(&Lim, 0);
	//for (int i = 0; i < I; i++)	
	//	for (int t = 0; t < T; t++)
	//	{
	//		if (t == 0)
	//		{
	//			Lim[i*T + t][1] = - sistema_a->termeletricasVtr[i].GetCoefCustoPartida()*sistema_a->termeletricasVtr[i].GetU0();
	//			//Lim[i*T + t][1] = - double (sistema_a->termeletricasVtr[i].GetU0());
	//			Lim[i*T + t][0] = 2;
	//		}
	//		else
	//		{
	//			Lim[i*T + t][1] = 0;
	//			Lim[i*T + t][0] = 2;
	//		}
	//	}
	//return Lim;
}
vetorfloat2 Hrstc::LimVmeta()
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
			//Lim(r,1) = sistema_a->hidreletricasVtr[r].GetVmax();
			Lim[r + R*cen][0] = 2;
		}
	}
	return Lim;
}
vetorfloat2 Hrstc::LimBalPotenciaL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Hrstc::LimBalVazaoL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Hrstc::LimFuncProdL()
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
vetorfloat2 Hrstc::LimPhMax()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 Hrstc::LimReserva()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, 1 * T, 2);
	IniciaMatriz(&Lim, 0);
	//double D;
	//int cen;
	for (int t = 0; t < T; t++)
	{
		//D = 0;
		//for (size_t b = 0; b < B; b++)
		//{
		//	if ((0 <= t) && (sistema_a->GetTt1() > t))
		//		for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
		//			D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);
		//	else
		//	{
		//		cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		//		for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
		//			D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
		//	}
		//}
		Lim[t][0] = 2;
		Lim[t][1] = sistema_a->GetReserva(t);
	}
	return Lim;
}
vetorfloat2 Hrstc::LimCortesF()
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
	// I * T
	//int I = sistema_a->termeletricasVtr.size();
	//vetorfloat2 Lim;
	//DimensionarMatriz(&Lim, I * T * sistema_a->GetFlagInitAproxCT(), 2);
	//IniciaMatriz(&Lim, 0);
	//int c = 0;
	//for (int i = 0; i < I; i++)
	//{
	//	for (int t = 0; t < T; t++)
	//	{
	//		for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)
	//		{
	//			Lim[c][0] = 2;
	//			c++;
	//		}
	//	}
	//}
	//return Lim;
}
// Matriz dos coeficientes e limites das restrições lineares
void Hrstc::MatrizRestricoesLineares(CMatrizEsparsa  * MM, CMatrizEsparsa &MMpl)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	// Para cada tipo de restrição cria-se uma matriz com todas as restrições, depois essas são alocadas para cada subproblema AlocarRestricoes()
	// matriz A contem todas as linhas das restrições com acoplamento temporal, para ser multiplicada pela solução do subp. anterior e gerar as condições iniciais do subp. atual
	// posicao_rest_din[modelo][i] contem (para cada modelo) a posição inicial das restrições (na matriz do subproblema) referente ao elemento i da matriz A

	// as restrições com acoplamento temporal (< t) sao registradas para serem alteradas em cada subproblema

	CMatrizEsparsa M(0);		//CMatrizEsparsa M;
	//posicao_rest_din.resize(6);		// seis elementos [ini fin], referentes as restricoes de Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida
	// na matriz esparsa percorrer elemento q n existe, em tese, n adiciona tempo algum, mas percorrer o RHS de tds as restricoes, por isso é bom saber os indices inicial e final das restricoes que sao atualizadas
	A = new CMatrizEsparsa[6 + flag7];
	// modelo PL
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		M = MatrizRestDemanda();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		M = MatrizLimFluxo();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		//M = MatrizLimFluxo();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
	}
	else	// == 0 e == 3(nesse caso os limites de fluxo somente são adicionados se violados)
	{
		M = MatrizRestDemandaBarraUnica();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		if ( sistema_a->GetFlagModeloRede() == 2 )
		{
			M = MatrizLimFluxo();
			MMpl.JuntarColuna(&M);
			AlocarRestricoes(&M, MM);
			//cout << "02" << endl;
			//M = MatrizLimFluxo();		// nao precisar calcular denovo a mesma matriz!!
			MMpl.JuntarColuna(&M);
			AlocarRestricoes(&M, MM);
		}
	}
	M = MatrizLimPhgMin();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	M = MatrizLimPhgMax();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	M = MatrizLimqMin();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	M = MatrizLimqMax();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	for (int i = 0; i < n_modelos; i++)			// grava indice_inicio. Todos os subproblemas tem o msm numero de restriçoes antes de adicionar as restricoes de volume meta e ptmin
		posicao_rest_din[i].push_back(MM[i].GetNlin());
		//posicao_rest_din[i].push_back(n_restricoes[i]);
	M = MatrizBalHid(); A[0] = M; posicao_rest_din_Ext[0] = MMpl.GetNlin();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(MM[i].GetNlin());
	M = MatrizTup(); A[1] = M; posicao_rest_din_Ext[1] = MMpl.GetNlin();
	MMpl.JuntarColuna(&M); 
	AlocarRestricoes(&M, MM);
	for (int i = 0; i < n_modelos; i++)
		posicao_rest_din[i].push_back(MM[i].GetNlin());
	M = MatrizTdown(); A[2] = M; posicao_rest_din_Ext[2] = MMpl.GetNlin();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizRampaUp(); A[3] = M; posicao_rest_din_Ext[3] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizRampaDown(); A[4] = M; posicao_rest_din_Ext[4] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizRestCP(); A[5] = M; posicao_rest_din_Ext[5] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		M = MatrizLimPtMax();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		M = MatrizLimPtMin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
	}
	else
	{
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizTupDown(); A[3] = M; posicao_rest_din_Ext[3] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizRampaUp(); A[4] = M; posicao_rest_din_Ext[4] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizRampaDown(); A[5] = M; posicao_rest_din_Ext[5] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		M = MatrizLimPtMax();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
		for (int i = 0; i < n_modelos; i++)
			posicao_rest_din[i].push_back(MM[i].GetNlin());
		M = MatrizLimPtMin(); A[6] = M; posicao_rest_din_Ext[6] = MMpl.GetNlin();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
	}
	M = MatrizVmeta();
	MMpl.JuntarColuna(&M);
	AlocarRestricoesVmeta(&M, MM);
	M = MatrizBalPotencia();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	M = MatrizBalVazao();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	M = MatrizFuncProd();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	if (sistema_a->GetFlagPhmax() == 1)
	{
		M = MatrizPhMax();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
	}
	M = MatrizReserva();
	MMpl.JuntarColuna(&M);
	AlocarRestricoes(&M, MM);
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		M = MatrizCortesF();
		MMpl.JuntarColuna(&M);
		AlocarRestricoes(&M, MM);
	}
}
void Hrstc::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor, vetorint * LimTipo_pl, vetorfloat * LimValor_pl)
{
	vetorfloat2 L;
	lim_iniciais.resize(6 + flag7);// = new vetorfloat2[6];		// Guardar valores dos limites das restriçoes originais
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		L = LimRestDemanda();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimFluxo0();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimFluxo2();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
	else
	{
		L = LimRestDemandaBarraUnica();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		if ( sistema_a->GetFlagModeloRede() == 2 )
		{
			L = LimFluxo0();
			AlocarLimites(&L, LimTipo_pl, LimValor_pl);
			AlocarLimitesSubp(L, LimTipo, LimValor);
			L = LimFluxo2();
			AlocarLimites(&L, LimTipo_pl, LimValor_pl);
			AlocarLimitesSubp(L, LimTipo, LimValor);
		}
	}
	L = LimPhgMin();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimPhgMax();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimQMin();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimQMax();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimBalHid(); lim_iniciais[0] = L;
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimTup(); lim_iniciais[1] = L;
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimTdown(); lim_iniciais[2] = L;
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		L = LimRampaUp(); lim_iniciais[3] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimRampaDown(); lim_iniciais[4] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimRestCP(); lim_iniciais[5] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimPtMax();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimPtMin();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
	else
	{
		L = LimTupDown(); lim_iniciais[3] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimRampaUp(); lim_iniciais[4] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimRampaDown(); lim_iniciais[5] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimPtMax();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
		L = LimPtMin(); lim_iniciais[6] = L;
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
	L = LimVmeta();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesVmetaSubp(L, LimTipo, LimValor);
	L = LimBalPotenciaL();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimBalVazaoL();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	L = LimFuncProdL();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	if (sistema_a->GetFlagPhmax() == 1)
	{
		L = LimPhMax();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
	L = LimReserva();
	AlocarLimites(&L, LimTipo_pl, LimValor_pl);
	AlocarLimitesSubp(L, LimTipo, LimValor);
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		L = LimCortesF();
		AlocarLimites(&L, LimTipo_pl, LimValor_pl);
		AlocarLimitesSubp(L, LimTipo, LimValor);
	}
}

Hrstc::Hrstc(CSistema * const sistema_end, Results * const resultadosGurobi_end, double max_time, istream *iStrm) : ambGRB(GRBEnv()),modeloPL(GRBModel(ambGRB))
{
	// ler parametros do arquivo
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , ALFA , double( 1.0 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , BETA , double( 0.5 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , GAMA , double( 1.0 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , NORM_P , int( 1 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , APPROACH , int( 1 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , INITCUT , int( 0 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , itBenders , int( 500 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , itBendersF , int( 0 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , WIN_SIZE , int( 1 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , TIME_R , int( 5 ) );
	OPTtypes_di_unipi_it::DfltdSfInpt( iStrm , coeftol , double( 1e-12 ) );
	WIN_FUT_SIZE = 0;
	BIN_MW = 2;
	BIN_FW = 1;

	// coeftol é usado para ponderar os coeficientes do corte, para evitar inviabilidades (erro numerico) devido a grande diferença na ordem dos coeficientes do corte
	//coeftol = 1e-8;		// 1e-6
	sistema_a = sistema_end;
	resultadosGurobi = resultadosGurobi_end;
	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
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
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= T * (sistema_a->barrasVtr.size() - 1);

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

	if (sistema_a->GetFlagVarBin() == 0)
	{
		flagH1 = 1;
		flagH2 = 1;
	}
	else
	{
		flagH1 = BIN_MW;
		flagH2 = BIN_FW;
	}
	flagH3 = NORM_P;

	// Criar modelos do gurobi (variaveis e restrições)
	CriarModelos(max_time);
}
Hrstc::Hrstc(CSistema * const sistema_end, Results * const resultadosGurobi_end, int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a, double max_time) : ambGRB(GRBEnv()),modeloPL(GRBModel(ambGRB))
{
	sistema_a = sistema_end;
	resultadosGurobi = resultadosGurobi_end;
	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	flag1 = int (sistema_a->GetFlagModeloRede());	
	flag2 = int (sistema_a->GetFlagVfol());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= T * (sistema_a->barrasVtr.size() - 1);

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

	if (sistema_a->GetFlagVarBin() == 0)
	{
		flagH1 = 1;
		flagH2 = 1;
	}
	else
	{
		flagH1 = bvmw;
		flagH2 = bvfw;
	}
	flagH3 = 2;

	if ( (NORM_P != 1) && ((beta == 0 && gama == 1) || (beta == 1 && gama == 0))  )
		flagH3 = 1;		// Usar norma 1 sempre que beta=0 e gama=1 ou beta=1 e gama=0

	// Criar modelos do gurobi (variaveis e restrições)
	CriarModelos(max_time);
}
Hrstc::~Hrstc(void)
{
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
	for (size_t i = 0; i < constr.size(); i++)
		delete constr[i];
	constr.clear();
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	delete varPL;
	delete X;
	delete[] A;
}
void Hrstc::CriarModelos(double max_time)
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
	nrest_orig.resize(n_modelos);		// vetor com o numero de restrições de cada subproblema
	
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
	fo_subp.resize(n_modelos);
	ff_subp.resize(n_modelos);
	modelosGRB.resize(n_modelos);
	vars.resize(n_modelos);
	constr.resize(n_modelos);
	numero_var_subp.resize(n_modelos);
	numero_var_subp_principal.resize(n_modelos);
	numero_var_ad.resize(n_modelos);
	numero_var_ad2.resize(n_modelos);
	posicao_rest_din.resize(n_modelos);
	debug_alfa_B2.resize(sistema_a->GetNCenarios());
	posicao_rest_din_Ext.resize(6 + flag7);
	restr_mod.resize(n_modelos);
	acomplamentos.resize(n_modelos - 1);		// não precisa-se desse vetor para o subproblema 0
	//PI = CMatrizEsparsa(1, n);
	L.resize(n_modelos - 1);					// a cada iteração é adicionada uma linha
	//L_bin.resize(n_modelos - 1);
	rhs_orig.resize(n_modelos - 1);
	//rhs_orig_bin.resize(n_modelos - 1);
	for (int i = 0; i < n_modelos - 1; i++)
	{
		L[i] = CMatrizEsparsa(0, n);
		//L_bin[i] = CMatrizEsparsa(0, n);
		rhs_orig[i].resize(0);
		//rhs_orig_bin[i].resize(0);
	}
	itB = 0;
	///DeterminarIndiceVarBin();
	//RHS = 0;
	foRL.resize(n_modelos);
	fo_ub_a = GRB_INFINITY;
	ind_var_fol.resize(n_modelos);
	//Cvfol = 20/(0.0036*sistema_a->GetDeltaT2())*sistema_a->GetCustoDeficit();
	//Cvfol = 20/0.0036*sistema_a->GetCustoDeficit();

	Cvfol = sistema_a->GetCustoVfol();
	// o problema do UB < LB estava aqui, com o Cvfol diferente na RL e na RP
	///Cvfol = (20/0.0036)*sistema_a->GetCustoDeficit();
	cout << "Valor da penalidade do vfol na Heurística é = " << Cvfol << endl;
	PenVF = 20*sistema_a->GetCustoDeficit();		// penalidade para as variáveis de folga no backward (adicionadas caso alguma restrição seja inviável)
	//PenVF = sistema_a->GetCustoVfol();
	//PenVF = (0.2/0.0036)*sistema_a->GetCustoDeficit();

	SelecionarVarProx();

	X = new CMatrizEsparsa(n, 1);
	//A = new CMatrizEsparsa[6];
	//lim_iniciais = new vetorfloat2[6];

	for (int i = 0; i < n_modelos; i++)
	{
		foRL[i] = 0;
		modelosGRB[i] = new GRBModel(ambGRB);		// um modelo de otimização para cada janela principal
		numero_var_subp[i] = 0;
		numero_var_subp_principal[i] = 0;
	}

	// numero de variaveis e var. acumuladas dos subproblemas
	// para as janelas principais identificar os indices dos subproblemas para a etapa backward
	int nt;
	nvar_acum.push_back(0);
	for (int i = 0; i < n_modelos; i++)
	{
		for (size_t ii = 0; ii < vtr_nos_pri[i].size(); ii++)
		{
			if ( (flag2 == 1) && (vtr_nos_pri[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
				nt = ( n - sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T + sistema_a->hidreletricasVtr.size();
			else
				nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
			numero_var_subp[i] += nt;
			numero_var_subp_principal[i] += nt;
		}
		for (size_t ii = 0; ii < vtr_nos_fut[i].size(); ii++)
		{
			if ( (flag2 == 1) && (vtr_nos_fut[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// tem-se a variável vfol nesse nó!!
				nt = ( n - sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T + sistema_a->hidreletricasVtr.size();
			else
				nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
			numero_var_subp[i] += nt;
		}
		nvar_acum.push_back(numero_var_subp_principal[i] + nvar_acum[i]);

		if (vtr_nos_pri[i][vtr_nos_pri[i].size() - 1] < sistema_a->GetTt1() - 1)	// subproblema é da 3 etapa
			subpB3etapa.push_back(i);
		else if (vtr_nos_pri[i][0] <= sistema_a->GetTt1() - 1)	// subproblema é da 2 etapa
		{
			subpB2etapa.push_back(i);
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				subpB2etapa.push_back(i + 1 + cen*(n_subp2+n_subpFH));
		}
		else // subproblema é de 1 etapa
			subpB1etapa.push_back(i);
	}
	//subpB3etapa: vetor com subproblemas do 1º estágio
	//subpB2etapa: elemento 0 é o subp de acoplam. e os n_cenarios elementos são os subproblemas futuros à ele
	//subpB1etapa: vetor com subproblemas do 2º estágio, em que a cada this.size / n_cenarios tem-se um subp de final de horizonte

	// Criar variáveis
	for (int i = 0; i < n_modelos; i++)
	{
		CriarVariaveis(i);
		modelosGRB[i]->update();	// Atualiza o modelo Gurobi.
		vars[i] = modelosGRB[i]->getVars();
	}
	
	// Criar variáveis modelo PL de todo o horizonte
	CriarVariaveisPL();
	modeloPL.update();
	varPL = modeloPL.getVars();

	// Adicionar restrições (Monta matrizes do problema inteiro e selecionar conjunto de restriçoes para cada subproblema)
	CriarRestricoes();
	
	// Criar restrições proximais
	if ((flagH3 == 1) && (alfa != 1))
		for (int i = 0; i < n_modelos; i++)
			CriarRestriçõesProximais(i);
	else
		for (int i = 0; i < n_modelos; i++)
		{
			numero_var_ad[i] = 0;
			numero_var_ad2[i] = 0;
		}

	// Fixar condições iniciais 
	for (int i = 0; i < n_modelos; i++)
	{
		FixarCondIniciais(i);
		modelosGRB[i]->update();	// Atualiza o modelo Gurobi.
		nrest_orig[i] = modelosGRB[i]->get(GRB_IntAttr_NumConstrs);	// numero de restrições de cada subproblema
	}

	// Fixar condições iniciais modelo PL
	FixarCondIniciaisPL();
	modeloPL.update();
	
	// Definir ajustes do solver
	for (int i = 0; i < n_modelos; i++)
	{
		//modelosGRB[i]->write("Hrst_subprob.lp");	// Escreve modelo em arquivo
		modelosGRB[i]->getEnv().set(GRB_IntParam_OutputFlag, 0);		// Nao escrever detalhes
		//modelosGRB[i].getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[i].getEnv().set(GRB_StringParam_LogFile, "logHrst.txt");		// Escrever log em arquivo
		modelosGRB[i]->getEnv().set(GRB_DoubleParam_MIPGap, 0.0001);	// Define o gap de tolerancia
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9); 
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_IntFeasTol, 1e-9); 
		//modelosGRB[i]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm.
		//modelosGRB[i]->getEnv().set(GRB_IntParam_Threads, 1);
		//modelosGRB[i]->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[i].getEnv().set(GRB_IntParam_Method, 0);
		//modelosGRB[i]->getEnv().set(GRB_DoubleParam_TimeLimit, max_time*TIME_R/n_modelos);		// Limita o tempo de resolução do problema
		modelosGRB[i]->getEnv().set(GRB_DoubleParam_TimeLimit, 10);		// Limita o tempo de resolução dos subproblemas
	}

	// Função objetivo e ajustes do modelo PL
	CriarFuncaoObjetivoPL();
	modeloPL.update();	// Atualiza o modelo Gurobi.
	
	modeloPL.getEnv().set(GRB_IntParam_OutputFlag, 0);		// Nao escrever detalhes
	//modeloPL.getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
	//modeloPL.getEnv().set(GRB_StringParam_LogFile, "log_ED.txt");		// Escrever log em arquivo
	//modeloPL.getEnv().set(GRB_DoubleParam_MIPGap, 0.000001);	// Define o gap de tolerancia
	//modeloPL.getEnv().set(GRB_IntParam_Threads, 4);
	modeloPL.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6); 
	modeloPL.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-6);
	//modeloPL.getEnv().set(GRB_DoubleParam_IntFeasTol, 1e-9);
	//modeloPL.getEnv().set(GRB_IntParam_Method, 2);
	//modeloPL.getEnv().set(GRB_IntParam_Crossover, 0);
	modeloPL.getEnv().set(GRB_DoubleParam_TimeLimit, 60);

	// Identificação acoplamentos
	for (int i = 1; i < n_modelos; i++)
		IdentAcoplamentos(i);

	// Adicinoar variáveis de folga ao subproblemas
	n_var_folga.resize(n_modelos);
	for (int i = 0; i < n_modelos; i++)
		n_var_folga[i] = 0;

	// Adição da variável que representa o custo futuro
	for (int i = 0; i < n_modelos; i++)
	{
		// Adicionar variável que representa a aproximação dos subproblemas futuros
		if ( i == subpB2etapa[0] )		// subproblema de acoplamento
		{
			for ( int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				modelosGRB[i]->addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "");
		}
		else
			modelosGRB[i]->addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "");
		//nvar_acum[mod + 1] += 1;	// não incluir alfa pois o vetor acomplamentos n considera o alfa, ele que endereça as var. de cada subp.
		//numero_var_subp_principal[mod] += 1;
		//numero_var_subp[mod] += 1;
		// variavel alfa não é considerada no x_spr nem no numero de váriaveis!! Como é só uma variável n precisa adicionar esse numero nos vetores!
		modelosGRB[i]->update();
		delete vars[i];
		vars[i] = modelosGRB[i]->getVars();
	}
}

void Hrstc::CriarVariaveis(int modelo)
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
			//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model

			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				for (size_t i = 0; i < I; i++)	//pt
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//pt
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				if  ( (flagH1 == 0) || (flagH1 == 1) )	// (fixo e livre continuo ou livre binario)
				{
					for (size_t i = 0; i < I; i++)		//u
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					for (size_t i = 0; i < I; i++)		//up
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					for (size_t i = 0; i < I; i++)		//ud
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)		//u
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
					for (size_t i = 0; i < I; i++)		//up
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
					for (size_t i = 0; i < I; i++)		//ud
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
				}
			}
			else
			{
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
			}
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F		//if (sistema_a->GetFlagInitAproxCT() > 1)
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagModeloRede() == 1)
				for (size_t b = 0; b < B; b++)	//teta
					modelosGRB[modelo]->addVar(-3.14159, 3.14159, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modelosGRB[modelo]->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), double(sistema_a->hidreletricasVtr[r].GetVmax(t+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), 0.0, GRB_CONTINUOUS, "");
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
			if (sistema_a->GetFlagPhmax() == 1)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
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
			if ( sistema_a->GetFlagModeloRede() > 0)
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
			if ( (t + 1 > sistema_a->GetTt1()) && ((t + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
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
			//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model

			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				for (size_t i = 0; i < I; i++)	//pt
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//pt
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				if  ( (flagH2 == 0) || (flagH2 == 1) )	// (fixo e livre continuo ou livre binario)
				{
					for (size_t i = 0; i < I; i++)		//u
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					for (size_t i = 0; i < I; i++)		//up
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					for (size_t i = 0; i < I; i++)		//ud
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)		//u
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
					for (size_t i = 0; i < I; i++)		//up
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
					for (size_t i = 0; i < I; i++)		//ud
						modelosGRB[modelo]->addVar(0, 1, 0.0, GRB_BINARY, "");
				}
			}
			else
			{
				if  ( (flagH2 == 0) || (flagH2 == 1) )	//u (fixo e livre continuo ou livre binario)
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
			}
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F		//if (sistema_a->GetFlagInitAproxCT() > 1)
			{
				for (size_t i = 0; i < I; i++)
					modelosGRB[modelo]->addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagModeloRede() == 1)
				for (size_t b = 0; b < B; b++)	//teta
					modelosGRB[modelo]->addVar(-3.14159, 3.14159, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modelosGRB[modelo]->addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), double(sistema_a->hidreletricasVtr[r].GetVmax(t+cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))), 0.0, GRB_CONTINUOUS, "");
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
			if (sistema_a->GetFlagPhmax() == 1)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modelosGRB[modelo]->addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
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
			if ( sistema_a->GetFlagModeloRede() > 0)
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
			if ( (t + 1 > sistema_a->GetTt1()) && ((t + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
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
			cout << "Exception during CriarVariaveis()" << endl;	
	}
}
void Hrstc::CriarVariaveisPL()
{
	try 
	{
		size_t I = sistema_a->termeletricasVtr.size();
		size_t R = sistema_a->hidreletricasVtr.size();
		size_t B = sistema_a->barrasVtr.size() - 1;
		// Estagio 1
		int T1 = sistema_a->GetTt1();
		for (int t = 0; t < T1; t++)
		{
			//x = [pt u cp F teta ph v d s phmax phg q z def]
			//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				for (size_t i = 0; i < I; i++)	//pt
					modeloPL.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//pt
					modeloPL.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagVarBin() == true)	//u (binario ou continuo)
			{
				for (size_t i = 0; i < I; i++)
					modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)
					modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				if (sistema_a->GetFlagVarBin() == true)
				{
					for (size_t i = 0; i < I; i++)		//up
						modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
					for (size_t i = 0; i < I; i++)		//ud
						modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)		//up
						modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					for (size_t i = 0; i < I; i++)		//ud
						modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//cp
					modeloPL.addVar(0, sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F		//if (sistema_a->GetFlagInitAproxCT() > 1)
			{
				for (size_t i = 0; i < I; i++)
					modeloPL.addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagModeloRede() == 1)
				for (size_t b = 0; b < B; b++)	//teta
					modeloPL.addVar(-3.14159, 3.14159, 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
			{
				double phmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
				modeloPL.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modeloPL.addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t)), double(sistema_a->hidreletricasVtr[r].GetVmax(t)), 0.0, GRB_CONTINUOUS, "");
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
			{
				double qhmax = 0;
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
				modeloPL.addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
			{
				modeloPL.addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagPhmax() == 1)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modeloPL.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//phg
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modeloPL.addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			for (size_t r = 0; r < R; r++)	//q
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
				{
					modeloPL.addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
				}
			}
			if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagModeloRede() > 0)
			{
				for (size_t b = 0; b < B + 1; b++)	//def
				{
					double cap_d = 0;
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(0, t);
					modeloPL.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				double cap_d = 0;			//def barra unica
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(0, t);
				modeloPL.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
			}
		}
		// Estagio 2
		int T2 = sistema_a->GetTt2() - sistema_a->GetTt1();
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (int t = 0; t < T2; t++)
			{
				//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				{
					for (size_t i = 0; i < I; i++)	//pt
						modeloPL.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)	//pt
						modeloPL.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagVarBin() == true)	//u (binario ou continuo)
				{
					for (size_t i = 0; i < I; i++)	//u
						modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
				}
				else
				{
					for (size_t i = 0; i < I; i++)	//u
						modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				{
					if (sistema_a->GetFlagVarBin() == true)
					{
						for (size_t i = 0; i < I; i++)		//up
							modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
						for (size_t i = 0; i < I; i++)		//ud
							modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
					}
					else
					{
						for (size_t i = 0; i < I; i++)		//up
							modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
						for (size_t i = 0; i < I; i++)		//ud
							modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					}
				}
				else
				{
					for (size_t i = 0; i < I; i++)	//cp
						modeloPL.addVar(0, sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				{
					for (size_t i = 0; i < I; i++)
						modeloPL.addVar(0, double (sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
				}
				if ( sistema_a->GetFlagModeloRede() == 1)
					for (size_t b = 0; b < B; b++)	//teta
						modeloPL.addVar(-3.14159, 3.14159, 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modeloPL.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					modeloPL.addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t+sistema_a->GetTt1()+cen*T2)), double(sistema_a->hidreletricasVtr[r].GetVmax(t+sistema_a->GetTt1()+cen*T2)), 0.0, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
				{
					double qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						qhmax = qhmax + double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
					modeloPL.addVar(0, qhmax + sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
				{
					modeloPL.addVar(0, sistema_a->hidreletricasVtr[r].GetSmax(), 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagPhmax() == 1)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
					{
						double phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						modeloPL.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
					}
				}
				for (size_t r = 0; r < R; r++)	//phg
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						modeloPL.addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				for (size_t r = 0; r < R; r++)	//q
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						modeloPL.addVar(0, double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades()), 0.0, GRB_CONTINUOUS, "");
					}
				}
				if (sistema_a->GetFlagVarBin() == true)	//z (binario ou continuo)
				{
					for (size_t r = 0; r < R; r++)	//z
					{
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							modeloPL.addVar(0, 1, 0.0, GRB_BINARY, "");
					}
				}
				else
				{
					for (size_t r = 0; r < R; r++)	//z
					{
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							modeloPL.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
					}
				}
				if ( sistema_a->GetFlagModeloRede() > 0)
				{
					for (size_t b = 0; b < B + 1; b++)	//def
					{
						double cap_d = 0;
						for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
							cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t + sistema_a->GetTt1());
						modeloPL.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
					}
				}
				else
				{
					double cap_d = 0;			//def barra unica
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
							cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t + sistema_a->GetTt1());
					modeloPL.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( sistema_a->GetFlagVfol() == true )
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//vfol
					modeloPL.addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during CriarVariaveisPL()" << endl;	
	}
}

void Hrstc::CriarRestricoes()
{
	try
	{
		///cout << "CriarRestricoes() " << n_modelos << endl;
		///std::cin.ignore();
		CMatrizEsparsa * MM = new CMatrizEsparsa[n_modelos];
		for (int i = 0; i < n_modelos; i++)
			MM[i].ZerarMatriz(0, numero_var_subp[i]);
		vetorint * LimTipo = new vetorint[n_modelos];
		vetorfloat * LimValor = new vetorfloat[n_modelos];
		//vetorint LimTipo;	// Tds iguais ao primeiro (subp = 0) subproblema (depois, antes de resolver cada subproblema atualizam-se o RHS, right-hand side, das restriçoes para subp > 0)
		// modelo PL
		CMatrizEsparsa MMpl(0, n);
		vetorint LimTipo_PL;
		vetorfloat LimValor_PL;
		///cout << "Matriz" << endl;
		MatrizRestricoesLineares(MM, MMpl);
		///cout << "Limites" << endl;
		MatrizLimitesLineares(LimTipo, LimValor, &LimTipo_PL, &LimValor_PL);
		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		int l;
		for (int mod = 0; mod < n_modelos; mod++)
		{
			///cout << "modelo " << mod << endl;
			for (int i = 0; i < MM[mod].GetNlin(); i++)		// loop no numero de restrições
			{
				l = MM[mod].GetValorLprim(i);
				while ( l != -1 )
				{
					coeficiente = MM[mod].GetValorVal(l);
					variavel = vars[mod][MM[mod].GetValorCol(l)];
					restricao.addTerms( &coeficiente, &variavel, 1);
					l = MM[mod].GetValorLprox(l);
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
			//restr_mod[mod] = n_restricoes[mod];
		}
		delete[] MM;
		delete[] LimTipo;
		delete[] LimValor;

		// modelo PL
		restricao.clear();
		for (int i = 0; i < MMpl.GetNlin(); i++)
		{
			l = MMpl.GetValorLprim(i);
			while ( l != -1 )
			{
				coeficiente = MMpl.GetValorVal(l);
				variavel = varPL[MMpl.GetValorCol(l)];
				restricao.addTerms( &coeficiente, &variavel, 1);
				l = MMpl.GetValorLprox(l);
			}
			switch (LimTipo_PL[i]) {
			case 0:
				modeloPL.addConstr(restricao, GRB_LESS_EQUAL, LimValor_PL[i], "");
				break;
			case 1:
				modeloPL.addConstr(restricao, GRB_EQUAL, LimValor_PL[i], "");
				break;
			case 2:
				modeloPL.addConstr(restricao, GRB_GREATER_EQUAL, LimValor_PL[i], "");
				break;
			default:
				cout << "Tipo inválido de restrição adicionada" << endl;
			}
			restricao.clear();
		}
		modeloPL.update();
		MMpl.ZerarMatriz();
		LimTipo_PL.clear();
		LimValor_PL.clear();
	}
	catch(GRBException e) 
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
		cout << "Exception during CriarRestricoes()" << endl;	
	}
}
void Hrstc::AlocarRestricoes(CMatrizEsparsa * matriz, CMatrizEsparsa * MM)
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
			if ( (vtr_nos_pri[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T + flag2*R;
			else
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T;
			// nt tem o numero de variáveis do nó vtr_nos_pri[i][ii]
			// nt_a tem o numero de variaveis dos nós que não sao do final de horizonte

			// dois loops, o primeiro é para varrer as colunas (ii) e o segundo para percorrer as linhas (iii) da matriz A
			lin_i = ii*nt_rest;
			for (size_t iii = ii; iii < vtr_nos_pri[i].size(); iii++)		// começa em ii, pq n tem-se restrições de acoplamento com o futuro (matriz triangular inferior), todas restrições são escritas como dependendo do passado
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
			if ( (vtr_nos_fut[i][ii] + 1 > sistema_a->GetTt1()) && ((vtr_nos_fut[i][ii] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T + flag2*R;
			else
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T;
			
			lin_i = nt_rest*vtr_nos_pri[i].size() + ii*nt_rest;
			for (size_t iii = ii; iii < vtr_nos_fut[i].size(); iii++)
			{
				M_subp.InserirMatriz(lin_i, col_i, lin_i + nt_rest - 1, col_i + nt - 1, matriz, nt_rest*vtr_nos_fut[i][iii], nt_a*vtr_nos_fut[i][ii] + flag2*cenario*R);
				lin_i += nt_rest;
			}
			col_i += nt;
		}
		//SparseMatriz(&M_subp, &indexL[i], &indexC[i], &indexV[i], &nnZ[i], &n_restricoes[i], &count[i]);
		//n_restricoes[i] += M_subp.GetNlin();
		MM[i].JuntarColuna(&M_subp);
	}
}
void Hrstc::AlocarRestricoesVmeta(CMatrizEsparsa * matriz, CMatrizEsparsa * MM)
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
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T + flag2*R;
				M_subp_a.InserirMatriz(0, col_i, R - 1, col_i + nt - 1, matriz, cenario*R, nt_a*vtr_nos_pri[i][ii] + flag2*cenario*R);
				M_subp.JuntarColuna(&M_subp_a);
				M_subp_a.RemoverTodosElementos();
			}
			else
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T;
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
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T + flag2*R;
				M_subp_a.InserirMatriz(0, col_i, R - 1, col_i + nt - 1, matriz, cenario*R, nt_a*vtr_nos_fut[i][ii] + flag2*cenario*R);
				M_subp.JuntarColuna(&M_subp_a);
				M_subp_a.RemoverTodosElementos();
			}
			else
				nt = ( n - flag2*sistema_a->GetNCenarios()*R ) / T;
			col_i += nt;
		}
		if (M_subp.GetNnz() != 0)		// nao incluir matriz vazia para nós q n sao do final do horizonte!!!
		{
			//SparseMatriz(&M_subp, &indexL[i], &indexC[i], &indexV[i], &nnZ[i], &n_restricoes[i], &count[i]);
			//n_restricoes[i] += M_subp.GetNlin();
			MM[i].JuntarColuna(&M_subp);
		}
	}
}
void Hrstc::AlocarLimitesSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor)
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
void Hrstc::AlocarLimitesVmetaSubp(vetorfloat2 &limites, vetorint * LimTipo, vetorfloat * LimValor)
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
void Hrstc::CriarRestriçõesProximais(int modelo)
{
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	vetorint vtr_a;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
		vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
	for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
		vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int nC = ((2 + flag4 - flag7)*I + B - 1 + (4+flag3)*R + 2*JJ + B);		// numero de variáveis continuas por periodo
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		nC -= (B - 1 + B - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		nC -= (B - 1);
	int n_mod = numero_var_subp[modelo];
	int delta = 0;
	int delta_p = 0;
	int jj, jj_p, cen;
	double custo;
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	try
	{
		// criar variáveis adicionais
		// termos proximais em torno de x_til
		numero_var_ad[modelo] = 0;
		numero_var_ad2[modelo] = 0;
		custo = double((1-alfa)*beta*gama < 1e-6 ? 0 : (1-alfa)*beta*gama);
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal/futuro
		{
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
			{}		// não adiciona-se termos proximais ao ultimo periodo do horizonte
			else
			{
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					for (size_t i = 0; i < I; i++)	//pt
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				else
					for (size_t i = 0; i < I; i++)	//pt
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if (sistema_a->GetFlagTbinaryModel() == 0)		// "old" model
					for (size_t i = 0; i < I; i++)	//cp
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					for (size_t i = 0; i < I; i++)
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if ( sistema_a->GetFlagModeloRede() == 1)
					for (size_t b = 0; b < B - 1; b++)	//teta
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if (sistema_a->GetFlagPhmax() == 1)
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < R; r++)	//phg
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < R; r++)	//q
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if ( sistema_a->GetFlagModeloRede() > 0)
					for (size_t b = 0; b < B; b++)	//def
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				else
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				//for (int i = 0; i < nC; i++)
				//	modelosGRB[modelo]->addVar(0, 9999, custo, GRB_CONTINUOUS, "");		// já adiciona custo de 1 na f.o., limite máximo de 1 pois todas variáveis são escalonadas em seus valores máximos
				numero_var_ad[modelo] += nC;
			}
		}
		// termos proximais em torno de x_hat
		custo = double((1-alfa)*(1-beta)*(1-gama) < 1e-6 ? 0 : (1-alfa)*(1-beta)*(1-gama));
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal/futuro
		{
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
			{}		// não adiciona-se termos proximais ao ultimo periodo do horizonte
			else
			{
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					for (size_t i = 0; i < I; i++)	//pt
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				else
					for (size_t i = 0; i < I; i++)	//pt
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if (sistema_a->GetFlagTbinaryModel() == 0)		// "modern" model
					for (size_t i = 0; i < I; i++)	//cp
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					for (size_t i = 0; i < I; i++)
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if ( sistema_a->GetFlagModeloRede() == 1)
					for (size_t b = 0; b < B - 1; b++)	//teta
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if (sistema_a->GetFlagPhmax() == 1)
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < R; r++)	//phg
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				for (size_t r = 0; r < R; r++)	//q
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				if ( sistema_a->GetFlagModeloRede() > 0)
					for (size_t b = 0; b < B; b++)	//def
						modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				else
					modelosGRB[modelo]->addVar(0, 1, custo, GRB_CONTINUOUS, "");
				//for (int i = 0; i < nC; i++)
				//	modelosGRB[modelo]->addVar(0, 9999, custo, GRB_CONTINUOUS, "");		// já adiciona custo de 1 na f.o., limite máximo de 1 pois todas variáveis são escalonadas em seus valores máximos
				numero_var_ad2[modelo] += nC;
			}
		}
		modelosGRB[modelo]->update();	// Atualiza o modelo Gurobi.
		delete vars[modelo];
		vars[modelo] = modelosGRB[modelo]->getVars();

		// Gravar numero de restrições originais
		restr_mod[modelo] = modelosGRB[modelo]->get(GRB_IntAttr_NumConstrs);

		// criar restrições adicionais
		GRBLinExpr restricao;
		for (int prox_term = 0; prox_term < 2; prox_term++)
		{	
			delta = 0;
			for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal/futuro
			{
				if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				{	// não adiciona-se termos proximais ao ultimo periodo do horizonte
					delta += nt + flag2*R;
					//delta_p += nC + flag2*R;
				}		
				else
				{
					// loop para cada tipo de variável, tem q pular as binarias! 2 restrições por variável
					if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
						cen = 0;
					else		// nós do estágio 2
						cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					//x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
					//x = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
					for (int i = 0; i < I; i++)
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
						{
							// pt
							restricao = vars[modelo][n_mod + i + delta_p] + vars[modelo][i + delta]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + i + delta_p] - vars[modelo][i + delta]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							// cp
							restricao = vars[modelo][n_mod + i  + I + delta_p] + vars[modelo][i + 2*I + delta]/abs(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + i  + I + delta_p] - vars[modelo][i + 2*I + delta]/abs(vars[modelo][i + 2*I + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						}
						else
						{
							// pt
							restricao = vars[modelo][n_mod + i + delta_p] + vars[modelo][i + delta]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + i + delta_p] - vars[modelo][i + delta]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						}
						if (sistema_a->GetFlagInitAproxCT() > 1)
						{
							// F
							restricao = vars[modelo][n_mod + (2-flag7)*I + i + delta_p] + vars[modelo][i + (3+flag7)*I + delta]/abs(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + (2-flag7)*I + i + delta_p] - vars[modelo][i + (3+flag7)*I + delta]/abs(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						}
					}
					if (sistema_a->GetFlagModeloRede() > 0)		// considerar rede
					{
						if ( sistema_a->GetFlagModeloRede() == 1)
						{
							for (size_t b = 0; b < B - 1; b++)
							{
								// teta
								restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + b + delta_p] + vars[modelo][b + I*(3 + flag4 + flag7) + delta]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB));
								modelosGRB[modelo]->addConstr(restricao, GRB_LESS_EQUAL, 0, "");
								restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + b + delta_p] - vars[modelo][b + I*(3 + flag4 + flag7) + delta]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB));
								modelosGRB[modelo]->addConstr(restricao, GRB_LESS_EQUAL, 0, "");
							}
						}
						for (size_t b = 0; b < B; b++)
						{
							//def
							if ( vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) > 0 )
							{
								restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + (4+flag3)*R + 2*JJ + b + delta_p] + vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB));
								modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
								restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + (4+flag3)*R + 2*JJ + b + delta_p] - vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB));
								modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							}
						}
					}
					else
					{
						//def
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + (4+flag3)*R + 2*JJ + delta_p] + vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta]/abs(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + (4+flag3)*R + 2*JJ + delta_p] - vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta]/abs(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
					}
					jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
					jj_p = I*(2 + flag4 - flag7) + flag1a*(B - 1) + (4+flag3)*R;
					for (int r = 0; r < R; r++)
					{
						// ph
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + r + delta_p] + vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + r + delta_p] - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						// v
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + R + r + delta_p] + vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + R + r + delta_p] - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						// d
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + 2*R + r + delta_p] + vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + 2*R + r + delta_p] - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						// s
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + 3*R + r + delta_p] + vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + 3*R + r + delta_p] - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB));
						modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						if (sistema_a->GetFlagPhmax() == 1)
						{
							// phmax
							restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + 4*R + r + delta_p] + vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + (2+flag4-flag7)*I + flag1a*(B - 1) + 4*R + r + delta_p] - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						}
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						{
							// phg
							restricao = vars[modelo][n_mod + j + jj_p + delta_p] + vars[modelo][j + jj + delta]/abs(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + j + jj_p + delta_p] - vars[modelo][j + jj + delta]/abs(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							// q
							restricao = vars[modelo][n_mod + j + jj_p + JJ + delta_p] + vars[modelo][j + jj + JJ + delta]/abs(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
							restricao = vars[modelo][n_mod + j + jj_p + JJ + delta_p] - vars[modelo][j + jj + JJ + delta]/abs(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB));
							modelosGRB[modelo]->addConstr(restricao, GRB_GREATER_EQUAL, 0, "");
						}
						jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
						jj_p += sistema_a->hidreletricasVtr[r].GetNGrupos();
					}
					delta += nt;
					delta_p += nC;
				}
			}
		}
		modelosGRB[modelo]->update();	// Atualiza o modelo Gurobi
		delete constr[modelo];
		constr[modelo] = modelosGRB[modelo]->getConstrs();
		restricao.clear();
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during CriarRestriçõesProximais" << endl;	
	}
}
void Hrstc::AtualizarRHSproximais(int modelo)
{
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();

	vetorint vtr_a;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
		vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
	for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
		vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int nC = ((2 + flag4 - flag7)*I + B - 1 + (4+flag3)*R + 2*JJ + B);		// numero de variável continua por periodo
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		nC -= (B - 1 + B - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		nC -= (B - 1);

	int n_mod = numero_var_subp[modelo];
	int delta = 0;
	int delta_x = 0;
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	int jj, cen;
	int contad = restr_mod[modelo];		// contador incremental das restrições
	double * ptr_x;
	try
	{
		for (int prox_term = 0; prox_term < 2; prox_term++)
		{	
			if (prox_term == 0)		// primeiro loop é para os termos proximais em x_til
				ptr_x = &x_til[0];
			else					// o segundo para os termos em x_hat
				ptr_x = &x_hat[0];
			delta = 0;
			for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal/futuro
			{
				if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				{	// não adiciona-se termos proximais ao ultimo periodo do horizonte
					delta += nt + flag2*R;
				}		
				else
				{
					// loop para cada tipo de variável, tem q pular as binarias! 2 restrições por variável
					if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
						cen = 0;
					else		// nós do estágio 2
						cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					delta_x = nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

					//x = [pt u cp (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
					//x = [pt u up ud (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
					for (int i = 0; i < I; i++)
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
						{
							// pt
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[i + delta_x]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB)));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[i + delta_x]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB)));
							// cp
							// podem ocorrer situações que o upper bound de cp seja 0, devido à fixação das soluções iniciais!
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[i + 2*I + delta_x]/abs(sistema_a->termeletricasVtr[i].GetCoefCustoPartida()));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[i + 2*I + delta_x]/abs(sistema_a->termeletricasVtr[i].GetCoefCustoPartida()));
						}
						else
						{
							// pt
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[i + delta_x]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB)));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[i + delta_x]/abs(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB)));
						}
						if (sistema_a->GetFlagInitAproxCT() > 1)
						{
							// F
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[i + (3+flag7)*I + delta_x]/abs(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB)));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[i + (3+flag7)*I + delta_x]/abs(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB)));
						}
					}
					if (sistema_a->GetFlagModeloRede() > 0)		// considerar rede
					{
						if ( sistema_a->GetFlagModeloRede() == 1)
						{
							for (size_t b = 0; b < B - 1; b++)
							{
								// teta
								constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[b + I*(3 + flag4 + flag7) + delta_x]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB)));
								constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[b + I*(3 + flag4 + flag7) + delta_x]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB)));
							}
						}
						for (size_t b = 0; b < B; b++)
						{
							//def
							if ( vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) > 0 )
							{
								constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB)));
								constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x]/abs(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB)));
							}
						}
					}
					else
					{
						//def
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x]/abs(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB)));
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x]/abs(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB)));
					}
					jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
					for (int r = 0; r < R; r++)
					{
						// ph
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB)));
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB)));
						// v
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB)));
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB)));
						//if (ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x] < vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB))
						//{
						//	cout << "x_til incorreto!!" << endl;
						//	std::cin.ignore();
						//	cout << r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x << endl;
						//	cout << ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x - 2] << " ";
						//	cout << ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x - 1] << " ";
						//	cout << ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x - 0] << " ";
						//	cout << ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x + 1] << " ";
						//	cout << ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x + 2] << endl;
						//}
						
						//if (r == 5)
						//	cout << ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x] << " : " << - ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x] << endl;
						// d
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB)));
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB)));
						// s
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB)));
						constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB)));
						if (sistema_a->GetFlagPhmax() == 1)
						{
							// phmax
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB)));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta_x]/abs(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB)));
						}
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						{
							// phg
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[j + jj + delta_x]/abs(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB)));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[j + jj + delta_x]/abs(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB)));
							// q
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, ptr_x[j + jj + JJ + delta_x]/abs(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB)));
							constr[modelo][contad++].set(GRB_DoubleAttr_RHS, - ptr_x[j + jj + JJ + delta_x]/abs(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB)));
						}
						jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
					}
					delta += nt;
				}
			}
		}
		modelosGRB[modelo]->update();	// Atualiza o modelo Gurobi
		delete constr[modelo];
		constr[modelo] = modelosGRB[modelo]->getConstrs();

		
		
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during AtualizarRHSproximais" << endl;	
	}
}

void Hrstc::CriarFuncaoObjetivo(int modelo, bool backward)
{
	// Aqui tem-se a opção de usar as normas 2, 1 e infinita para os termos proximais
	// Norma infinita ainda não implementada!

	if (flagH3 == 1)
	{
		// Zerar função objetivo
		GRBLinExpr obj;
		modelosGRB[modelo]->setObjective(obj);
		modelosGRB[modelo]->update();

		CriarFuncaoObjetivoL(modelo, backward);
		if (ALFA != 1)
			AtualizarRHSproximais(modelo);
	}
	else
	{
		// Zerar função objetivo
		GRBLinExpr obj;
		modelosGRB[modelo]->setObjective(obj);
		modelosGRB[modelo]->update();

		CriarFuncaoObjetivoQ(modelo, backward);
	}
}
void Hrstc::CriarFuncaoObjetivoQ(int modelo, bool backward)
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
		int flag1a = 0;
		if (flag1 == 1)
			flag1a = 1;
		int flag1d = 1;		// referente à var. def
		if (flag1 == 0)
			flag1d = 0;
		// Funçao objetivo normalizada com os valores de cada subproblema
		//CalcularFuncaoObjetivo(modelo) -> só calcula a f.o. referente as janelas principais!!
		//double optimal_LR = 1; cout << "optimal_LR = 1!" << endl;
		double optimal_LR;
		if (backward == true || alfa == 1)
			optimal_LR = 1;
		else
			optimal_LR = foRL[modelo];		//calcula a f.o. referente as janelas principais e futuras!!

		int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
		int cen, jj;
		double deltaT;
		int delta = 0;
		int delta_x = 0;
		double constante = 0;
		double alfa_a;
		// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
		vetorint vtr_a;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
			vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
		for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
			vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);
		
		// Variáveis das janelas principais e futuras
		//x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		//x = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		// ident_var[3+flag7+flag4+flag1a+4+flag3+3+flag1d]

		// Termos lineares
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal e futura
		{
			// Se for um periodo de final de horizonte considerar somente a função objetivo original!
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				alfa_a = 1.0;
			else
				alfa_a = alfa;

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				cen = 0;
				deltaT = sistema_a->GetDeltaT1();
			}
			else		// nós do estágio 2
			{
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
			}
			delta_x = nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					vars[modelo][i + (3+flag7)*I + delta].set(GRB_DoubleAttr_Obj, 1*deltaT * alfa_a/optimal_LR - ident_var[3+flag7]*2*(1-alfa_a)*(beta*gama*x_til[i + (3+flag7)*I + delta_x] + (1-beta)*(1-gama)*x_hat[i + (3+flag7)*I + delta_x])/pow(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB), 2));	//F
					if (sistema_a->GetFlagTbinaryModel() == 0)
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, - ident_var[0]*2*(1-alfa_a)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2));				//pt
					else
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, - ident_var[0]*2*(1-alfa_a)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB), 2));				//pt
					vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, ident_var[1]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + I + delta_x])));		//u
					constante += ident_var[3+flag7]*(1 - alfa_a)*(beta*gama*pow(x_til[i + (3+flag7)*I + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + (3+flag7)*I + delta_x], 2))/pow(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB), 2);		//F
				}
				else
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * alfa_a/optimal_LR - ident_var[0]*2*(1-alfa_a)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2));	//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * alfa_a/optimal_LR + ident_var[1]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + I + delta_x])));		//u
					}
					else
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * alfa_a/optimal_LR - ident_var[0]*2*(1-alfa_a)*(beta*gama*x_til[i + delta_x] + (1-beta)*(1-gama)*x_hat[i + delta_x])/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB), 2));	//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT*alfa_a/optimal_LR + ident_var[1]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + I + delta_x])));		//u
					}
				}
				if (sistema_a->GetFlagTbinaryModel() == 0)
				{
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT * alfa_a/optimal_LR - ident_var[2]*2*(1-alfa_a)*(beta*gama*x_til[i + 2*I + delta_x] + (1-beta)*(1-gama)*x_hat[i + 2*I + delta_x])/pow(sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 2));		//cp
					constante += ident_var[0]*(1 - alfa_a)*(beta*gama*pow(x_til[i + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + delta_x], 2))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB), 2);				//pt
					constante += ident_var[2]*(1 - alfa_a)*(beta*gama*pow(x_til[i + 2*I + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + 2*I + delta_x], 2))/pow(sistema_a->termeletricasVtr[i].GetCoefCustoPartida(), 2);		//cp
				}
				else
				{
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() * alfa_a/optimal_LR + ident_var[2]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + 2*I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + 2*I + delta_x])));		//up
					vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_Obj, ident_var[3]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + 3*I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + 3*I + delta_x])));		//ud
					constante += ident_var[0]*(1 - alfa_a)*(beta*gama*pow(x_til[i + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[i + delta_x], 2))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB), 2);				//pt
					constante += ident_var[2]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[i + 2*I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + 2*I + delta_x], 2));		//up
					constante += ident_var[3]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[i + 3*I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + 3*I + delta_x], 2));		//ud
				}
				constante += ident_var[1]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[i + I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + I + delta_x], 2));	//u
			}
			if (sistema_a->GetFlagModeloRede() > 0)		// considerar rede
			{
				if ( sistema_a->GetFlagModeloRede() == 1)
				{
					for (size_t b = 0; b < B - 1; b++)
					{
						vars[modelo][b + I*(3 + flag4 + flag7) + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag4+flag7]*2*(1-alfa_a)*(beta*gama*x_til[b + I*(3 + flag4 + flag7) + delta_x] + (1-beta)*(1-gama)*x_hat[b + I*(3 + flag4 + flag7) + delta_x])/pow(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB), 2));		//teta
						constante += ident_var[3+flag4+flag7]*(1 - alfa_a)*(beta*gama*pow(x_til[b + I*(3 + flag4 + flag7) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[b + I*(3 + flag4 + flag7) + delta_x], 2))/pow(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB), 2);		//teta
					}
				}
				for (size_t b = 0; b < B; b++)
				{
					if ( vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) == 0 )
						vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa_a/optimal_LR);		//def
					else
					{
						vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa_a/optimal_LR - ident_var[3+flag7+flag4+flag1a+4+flag3+3]*2*(1-alfa_a)*(beta*gama*x_til[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x] + (1-beta)*(1-gama)*x_hat[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x])/pow(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2));		//def
						constante += ident_var[3+flag7+flag4+flag1a+4+flag3+3]*(1 - alfa_a)*(beta*gama*pow(x_til[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x], 2))/pow(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2);		//def
					}
				}
			}
			else
			{
				vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa_a/optimal_LR - ident_var[3+flag7+flag4+flag1a+4+flag3+3]*2*(1-alfa_a)*(beta*gama*x_til[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x] + (1-beta)*(1-gama)*x_hat[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x])/pow(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2));		//def
				constante += ident_var[3+flag7+flag4+flag1a+4+flag3+3]*(1 - alfa_a)*(beta*gama*pow(x_til[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x], 2))/pow(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB), 2);		//def
			}
			jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
			for (int r = 0; r < R; r++)
			{
				vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a]*2*(1-alfa_a)*(beta*gama*x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta_x])/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB), 2));				//ph
				vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a+1]*2*(1-alfa_a)*(beta*gama*x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x])/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB), 2));		//v
				vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a+2]*2*(1-alfa_a)*(beta*gama*x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB), 2));	//d
				vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a+3]*2*(1-alfa_a)*(beta*gama*x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB), 2));	//s
				constante += ident_var[3+flag7+flag4+flag1a]*(1 - alfa_a)*(beta*gama*pow(x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB), 2);			//ph
				constante += ident_var[3+flag7+flag4+flag1a+1]*(1 - alfa_a)*(beta*gama*pow(x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB), 2);		//v
				constante += ident_var[3+flag7+flag4+flag1a+2]*(1 - alfa_a)*(beta*gama*pow(x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB), 2);	//d
				constante += ident_var[3+flag7+flag4+flag1a+3]*(1 - alfa_a)*(beta*gama*pow(x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB), 2);	//s
				if (sistema_a->GetFlagPhmax() == 1)
				{
					vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a+4]*2*(1-alfa_a)*(beta*gama*x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta_x])/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB), 2));	//phmax
					constante += ident_var[3+flag7+flag4+flag1a+4]*(1 - alfa_a)*(beta*gama*pow(x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB), 2);	//phmax
				}
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					vars[modelo][j + jj + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a+4+flag3]*2*(1-alfa_a)*(beta*gama*x_til[j + jj + delta_x] + (1-beta)*(1-gama)*x_hat[j + jj + delta_x])/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB), 2));							//phg
					vars[modelo][j + jj + JJ + delta].set(GRB_DoubleAttr_Obj, - ident_var[3+flag7+flag4+flag1a+4+flag3+1]*2*(1-alfa_a)*(beta*gama*x_til[j + jj + JJ + delta_x] + (1-beta)*(1-gama)*x_hat[j + jj + JJ + delta_x])/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB), 2));		//q
					vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+2]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[j + jj + 2*JJ + delta_x]) + (1-beta)*gama*(1-2*x_hat[j + jj + 2*JJ + delta_x])));														//z
					constante += ident_var[3+flag7+flag4+flag1a+4+flag3]*(1 - alfa_a)*(beta*gama*pow(x_til[j + jj + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[j + jj + delta_x], 2))/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB), 2);					//phg
					constante += ident_var[3+flag7+flag4+flag1a+4+flag3+1]*(1 - alfa_a)*(beta*gama*pow(x_til[j + jj + JJ + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[j + jj + JJ + delta_x], 2))/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB), 2);		//q
					constante += ident_var[3+flag7+flag4+flag1a+4+flag3+2]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[j + jj + 2*JJ + delta_x], 2) + (1-beta)*gama*pow(x_hat[j + jj + 2*JJ + delta_x], 2));	//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
			{
				if (sistema_a->GetFlagVfol() == true)
					for (int r = 0; r < R; r++)
						vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].set(GRB_DoubleAttr_Obj, Cvfol*deltaT*alfa_a/optimal_LR);		//vfol
						//vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].set(GRB_DoubleAttr_Obj, Cvfol*deltaT*alfa_a/optimal_LR - 2*(1-alfa_a)*(beta*gama*x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta_x] + (1-beta)*(1-gama)*x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta_x])/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].get(GRB_DoubleAttr_UB), 2));		//vfol
						//constante += (1 - alfa_a)*(beta*gama*pow(x_til[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta_x], 2) + (1-beta)*(1-gama)*pow(x_hat[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta_x], 2))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].get(GRB_DoubleAttr_UB), 2);		//vfol
				delta += nt + flag2*R;
			}
			else
				delta += nt;
		}
		// adicionar alfa da FCF
		if ( modelo == subpB2etapa[0] )		// subproblema de acoplamento
		{
			for ( int cenn = 0; cenn < sistema_a->GetNCenarios(); cenn++)
				vars[modelo][numero_var_subp[modelo] + n_var_folga[modelo] + cenn].set(GRB_DoubleAttr_Obj, alfa_a/optimal_LR);		// variavel que representa a aproximação do custo dos suproblemas futuros
		}
		else
			vars[modelo][numero_var_subp[modelo] + n_var_folga[modelo]].set(GRB_DoubleAttr_Obj, alfa_a/optimal_LR);		// variavel que representa a aproximação do custo dos suproblemas futuros
		modelosGRB[modelo]->update();

		// Termos quadráticos
		delta = 0;
		GRBQuadExpr fo;
		//x = [pt u cp (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		//x = [pt u up ud (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		double * coeffs;
		//coeffs = new double[numero_var_subp[modelo] + 1];		// +1 por causa da FCF
		//for (int i = 0; i < numero_var_subp[modelo] + 1; i++)
		//	coeffs[i] = vars[modelo][i].get(GRB_DoubleAttr_Obj);
		//fo.addTerms(coeffs, vars[modelo], numero_var_subp[modelo] + 1);
		//for (int i = 0; i < numero_var_subp[modelo] + 1; i++)		// nao precisaria pois o vetor coeffs é todo completo, a n ser q alguma variavel n seja considerada no termo proximal
		//	coeffs[i] = 0;

		coeffs = new double[modelosGRB[modelo]->get(GRB_IntAttr_NumVars)];		// +1 por causa da FCF
		for (int i = 0; i < modelosGRB[modelo]->get(GRB_IntAttr_NumVars); i++)
			coeffs[i] = vars[modelo][i].get(GRB_DoubleAttr_Obj);
		fo.addTerms(coeffs, vars[modelo], modelosGRB[modelo]->get(GRB_IntAttr_NumVars));
		for (int i = 0; i < modelosGRB[modelo]->get(GRB_IntAttr_NumVars); i++)		// nao precisaria pois o vetor coeffs é todo completo, a n ser q alguma variavel n seja considerada no termo proximal
			coeffs[i] = 0;

		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal e futura
		{
			// Se for um periodo de final de horizonte considerar somente a função objetivo original!
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				alfa_a = 1.0;
			else
				alfa_a = alfa;

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				cen = 0;
				deltaT = sistema_a->GetDeltaT1();
			}
			else		// nós do estágio 2
			{
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
			}

			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					coeffs[i + (3 + flag7)*I + delta] = ident_var[3+flag7]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + (3+flag7)*I + delta].get(GRB_DoubleAttr_UB),2);		//F^2
				}
				if (sistema_a->GetFlagInitAproxCT() == 0)
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
						coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT*alfa_a/optimal_LR + ident_var[0]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB),2);		//pt^2
					else
					{
						cout << "CriarFuncaoObjetivo2: Termo quadrático não implementado para essa modelagem!!" << endl;
						break;
						// a2*(pt^2 + 2*u*pt*pt_min + u^2*pt_min^2) = a2*(pt^2 + 2*u*pt*pt_min + u*pt_min^2)
					}
				}
				else
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
						coeffs[i + delta] = ident_var[0]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB),2);		//pt^2
					else
						coeffs[i + delta] = ident_var[0]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][i + delta].get(GRB_DoubleAttr_UB) - vars[modelo][i + delta].get(GRB_DoubleAttr_LB),2);		//pt^2
				}
				if (sistema_a->GetFlagTbinaryModel() == 0)
				{
					coeffs[i + 2*I + delta] = ident_var[2]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(sistema_a->termeletricasVtr[i].GetCoefCustoPartida(),2);		//cp^2
				}
				else
				{
					// já foram considerados como linear anteriormente, pois u^2 = u
					//coeffs[i + 2*I + delta] = (1 - alfa_a)*(beta*(1-gama) + (1-beta)*gama);		//up^2
					//coeffs[i + 3*I + delta] = (1 - alfa_a)*(beta*(1-gama) + (1-beta)*gama);		//ud^2
				}
				//coeffs[i + I + delta] = (1 - alfa_a)*(beta*(1-gama) + (1-beta)*gama);		//u^2
			}
			if (sistema_a->GetFlagModeloRede() > 0)		// considerar rede
			{
				if ( sistema_a->GetFlagModeloRede() == 1)
				{
					for (size_t b = 0; b < B - 1; b++)
					{
						coeffs[b + I*(3 + flag4 + flag7) + delta] = ident_var[3+flag7+flag4]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_UB) - vars[modelo][b + I*(3 + flag4 + flag7) + delta].get(GRB_DoubleAttr_LB),2);		//teta
					}
				}
				for (size_t b = 0; b < B; b++)
				{
					if ( vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB) != 0 )
					{
						coeffs[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta] = ident_var[3+flag7+flag4+flag1a+4+flag3+3]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB),2);		//def
					}
					else
					{
						coeffs[b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta] = 0;		//def
					}
				}
			}
			else
			{
				coeffs[I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta] = ident_var[3+flag7+flag4+flag1a+4+flag3+3]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].get(GRB_DoubleAttr_UB),2);		//def
			}

			jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
			for (int r = 0; r < R; r++)
			{
				coeffs[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta] = ident_var[3+flag7+flag4+flag1a]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + delta].get(GRB_DoubleAttr_UB),2);				//ph
				coeffs[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta] = ident_var[3+flag7+flag4+flag1a+1]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_UB) - vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + R + delta].get(GRB_DoubleAttr_LB),2);	//v
				coeffs[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta] = ident_var[3+flag7+flag4+flag1a+2]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 2*R + delta].get(GRB_DoubleAttr_UB),2);	//d
				coeffs[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta] = ident_var[3+flag7+flag4+flag1a+3]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 3*R + delta].get(GRB_DoubleAttr_UB),2);	//s
				if (sistema_a->GetFlagPhmax() == 1)
				{
					coeffs[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta] = ident_var[3+flag7+flag4+flag1a+4]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + 4*R + delta].get(GRB_DoubleAttr_UB),2);	//phmax
				}
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					coeffs[j + jj + delta] = ident_var[3+flag7+flag4+flag1a+4+flag3]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][j + jj + delta].get(GRB_DoubleAttr_UB),2);					//phg
					coeffs[j + jj + JJ + delta] = ident_var[3+flag7+flag4+flag1a+4+flag3+1]*(1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][j + jj + JJ + delta].get(GRB_DoubleAttr_UB),2);		//q
					//coeffs[j + jj + 2*JJ + delta] = (1 - alfa_a)*(beta*(1-gama) + (1-beta)*gama);		//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
			{
				if (sistema_a->GetFlagVfol() == true)
					for (int r = 0; r < R; r++)
					{
						//coeffs[r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta] = (1 - alfa_a)*(beta*gama + (1-beta)*(1-gama))/pow(vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].get(GRB_DoubleAttr_UB),2);		//vfol
					}
				delta += nt + flag2*R;
			}
			else
				delta += nt;
		}
		fo.addTerms(coeffs, vars[modelo], vars[modelo], modelosGRB[modelo]->get(GRB_IntAttr_NumVars));
		modelosGRB[modelo]->setObjective(fo, GRB_MINIMIZE);

		// Termos constantes (devem ser adicionados após o setObjective com GRBQuadExpr, pois setObjective() replaces the entire existing objective)
		//int n_bin_janela_pri = vtr_nos_pri[modelo].size() * (I + JJ);
		//int n_con_janela_pri = numero_var_subp[modelo] - n_bin_janela_pri;
		modelosGRB[modelo]->set(GRB_DoubleAttr_ObjCon, constante);
		modelosGRB[modelo]->update();
		delete coeffs;
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during CriarFuncaoObjetivoQ()" << endl;	
	}
}
void Hrstc::CriarFuncaoObjetivoL(int modelo, bool backward)
{
	// Norma inifinita não implemntada! Somente a norma 1.

	// Modela a f.o. do subproblema com termo proximal linear (norma 1 ou infinita)
	// Para as variáveis continuas:
	// Na norma 1 cria-se uma variável e uma restrição a mais para cada variável continua, portanto nC var. e nC rest. a mais
	// Na norma infinita cria-se somente uma variável a mais e nC restrições por variável continua
	// Para as variáveis binárias:
	// Em ambas as normas utiliza-se a norma quadrática, visto que por ser uma variável binária essa função se torna linear

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
		int delta_x = 0;
		double deltaT;
		double constante = 0;
		double alfa_a;
		int flag1a = 0;
		if (flag1 == 1)
			flag1a = 1;
		int flag1d = 1;		// referente à var. def
		if (flag1 == 0)
			flag1d = 0;
		double optimal_LR;
		if  (backward == true || alfa == 1)
			optimal_LR = 1;
		else
			optimal_LR = foRL[modelo];
		// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
		vetorint vtr_a;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
			vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
		for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
			vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

		//x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		//x = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		// ident_var[3+flag7+flag4+flag1a+4+flag3+3+flag1d]*

		// Termos lineares
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)
		{
			// Se for um periodo de final de horizonte considerar somente a função objetivo original!
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				alfa_a = 1.0;
			else
				alfa_a = alfa;

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				cen = 0;
				deltaT = sistema_a->GetDeltaT1();
			}
			else		// nós do estágio 2
			{
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
			}
			delta_x = nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					vars[modelo][i + (3+flag7)*I + delta].set(GRB_DoubleAttr_Obj, 1*deltaT * alfa_a/optimal_LR);	//F
					vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, ident_var[1]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + I + delta_x])));		//u
				}
				else
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * alfa_a/optimal_LR);	//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * alfa_a/optimal_LR + ident_var[1]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + I + delta_x])));		//u
					}
					else
					{
						vars[modelo][i + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * alfa_a/optimal_LR);	//pt
						vars[modelo][i + I + delta].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT*alfa_a/optimal_LR + ident_var[1]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + I + delta_x])));		//u
					}
				}
				if (sistema_a->GetFlagTbinaryModel() == 0)
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT * alfa_a/optimal_LR);		//cp
				else
				{
					vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_Obj, deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() * alfa_a/optimal_LR + ident_var[2]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + 2*I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + 2*I + delta_x])));		//up
					vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_Obj, ident_var[3]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[i + 3*I + delta_x]) + (1-beta)*gama*(1-2*x_hat[i + 3*I + delta_x])));		//ud
					constante += ident_var[2]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[i + 2*I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + 2*I + delta_x], 2));		//up
					constante += ident_var[3]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[i + 3*I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + 3*I + delta_x], 2));		//ud
				}
				constante += ident_var[1]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[i + I + delta_x], 2) + (1-beta)*gama*pow(x_hat[i + I + delta_x], 2));	//u
			}
			if (sistema_a->GetFlagModeloRede() > 0)		// considerar rede
			{
				for (size_t b = 0; b < B; b++)
					vars[modelo][b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa_a/optimal_LR);		//def
			}
			else
				vars[modelo][I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta].set(GRB_DoubleAttr_Obj, sistema_a->GetCustoDeficit()*deltaT * alfa_a/optimal_LR);		//def
			jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
			for (int r = 0; r < R; r++)
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					vars[modelo][j + jj + 2*JJ + delta].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+2]*(1-alfa_a)*(beta*(1-gama)*(1-2*x_til[j + jj + 2*JJ + delta_x]) + (1-beta)*gama*(1-2*x_hat[j + jj + 2*JJ + delta_x])));	//z
					constante += ident_var[3+flag7+flag4+flag1a+4+flag3+2]*(1 - alfa_a)*(beta*(1-gama)*pow(x_til[j + jj + 2*JJ + delta_x], 2) + (1-beta)*gama*pow(x_hat[j + jj + 2*JJ + delta_x], 2));	//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
			{
				if (sistema_a->GetFlagVfol() == true)
					for (int r = 0; r < R; r++)
						vars[modelo][r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta].set(GRB_DoubleAttr_Obj, Cvfol*deltaT*alfa_a/optimal_LR);		//vfol
				delta += nt + flag2*R;
			}
			else
				delta += nt;
		}
		if ( modelo == subpB2etapa[0] )		// subproblema de acoplamento
		{
			for ( int cenn = 0; cenn < sistema_a->GetNCenarios(); cenn++)
				vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + numero_var_ad2[modelo] + n_var_folga[modelo] + cenn].set(GRB_DoubleAttr_Obj, alfa_a/optimal_LR);		// variavel que representa a aproximação do custo dos suproblemas futuros
		}
		else
			vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + numero_var_ad2[modelo] + n_var_folga[modelo]].set(GRB_DoubleAttr_Obj, alfa_a/optimal_LR);		// variavel que representa a aproximação do custo dos suproblemas futuros
		modelosGRB[modelo]->update();

		// Termos quadráticos
		delta = 0;
		GRBQuadExpr fo;
		//x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		//x = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		double * coeffs;
		coeffs = new double[modelosGRB[modelo]->get(GRB_IntAttr_NumVars)];
		for (int i = 0; i < modelosGRB[modelo]->get(GRB_IntAttr_NumVars); i++)
			coeffs[i] = vars[modelo][i].get(GRB_DoubleAttr_Obj);
		fo.addTerms(coeffs, vars[modelo], modelosGRB[modelo]->get(GRB_IntAttr_NumVars));
		for (int i = 0; i < modelosGRB[modelo]->get(GRB_IntAttr_NumVars); i++)
			coeffs[i] = 0;
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)
		{
			// Se for um periodo de final de horizonte considerar somente a função objetivo original!
			if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				alfa_a = 1.0;
			else
				alfa_a = alfa;

			if ((0 <= vtr_a[vtr_i]) && (vtr_a[vtr_i]) < sistema_a->GetTt1())		// nós do estágio 1
			{
				cen = 0;
				deltaT = sistema_a->GetDeltaT1();
			}
			else		// nós do estágio 2
			{
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
			}

			for (int i = 0; i < I; i++)
			{
				if (sistema_a->GetFlagInitAproxCT() == 0)
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
						coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT*alfa_a/optimal_LR;		//pt^2
					else
					{
						cout << "CriarFuncaoObjetivo2: Termo quadrático não implementado para essa modelagem!!" << endl;
						break;
						// a2*(pt^2 + 2*u*pt*pt_min + u^2*pt_min^2) = a2*(pt^2 + 2*u*pt*pt_min + u*pt_min^2)
					}
				}
			}
		}
		fo.addTerms(coeffs, vars[modelo], vars[modelo], modelosGRB[modelo]->get(GRB_IntAttr_NumVars));		// só inclui até o numero de variaveis originais, tem problema?
		modelosGRB[modelo]->setObjective(fo, GRB_MINIMIZE);
		delete[] coeffs;

		// Termos constantes (devem ser adicionados após o setObjective com GRBQuadExpr, pois setObjective() replaces the entire existing objective)
		//int n_bin_janela_pri = vtr_nos_pri[modelo].size() * (I + JJ);
		//int n_con_janela_pri = numero_var_subp[modelo] - n_bin_janela_pri;
		modelosGRB[modelo]->set(GRB_DoubleAttr_ObjCon, constante);
		modelosGRB[modelo]->update();


		// Definir custo das variáveis adicionais (usadas para representar o termo proximal linear, módulo)
		if (ALFA != 1)
		{
			double CTPL = 0;		// custo do termo proximal linear
			if  (!backward)
				CTPL = (1-alfa)*beta*gama;
			// termos proximais em torno de x_til
			for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal/futuro
			{
				if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				{}		// não adiciona-se termos proximais ao ultimo periodo do horizonte
				else
				{
					// terá q ser feito um loop para cada tipo de variável (para selecioná-las nos termos proximais)
					int ii = 0;
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
						for (size_t i = 0; i < I; i++)	//pt
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[0]*CTPL);
					else
						for (size_t i = 0; i < I; i++)	//pt
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[0]*CTPL);
					if (sistema_a->GetFlagTbinaryModel() == 0)		// "old" model
						for (size_t i = 0; i < I; i++)	//cp
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[2]*CTPL);
					if (sistema_a->GetFlagInitAproxCT() > 1)	//F
						for (size_t i = 0; i < I; i++)
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7]*CTPL);
					if ( sistema_a->GetFlagModeloRede() == 1)
						for (size_t b = 0; b < B - 1; b++)	//teta
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
						vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+1]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
						vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+2]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
						vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+3]*CTPL);
					if (sistema_a->GetFlagPhmax() == 1)
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4]*CTPL);
					for (size_t r = 0; r < R; r++)	//phg
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3]*CTPL);
					for (size_t r = 0; r < R; r++)	//q
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+1]*CTPL);
					if ( sistema_a->GetFlagModeloRede() > 0)
						for (size_t b = 0; b < B; b++)	//def
							vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+3]*CTPL);
					else
						vars[modelo][numero_var_subp[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+3]*CTPL);
					//for (int i = 0; i < numero_var_ad[modelo]; i++)
					//	vars[modelo][numero_var_subp[modelo] + i].set(GRB_DoubleAttr_Obj, CTPL);
				}
			}
			// termos proximais em torno de x_hat
			CTPL = 0;
			if  (!backward)
				CTPL = (1-alfa)*(1-beta)*(1-gama);
			for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)		// referente as var. da janela principal/futuro
			{
				if ((vtr_a[vtr_i] >= sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0))
				{}		// não adiciona-se termos proximais ao ultimo periodo do horizonte
				else
				{
					int ii = 0;
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
						for (size_t i = 0; i < I; i++)	//pt
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[0]*CTPL);
					else
						for (size_t i = 0; i < I; i++)	//pt
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[0]*CTPL);
					if (sistema_a->GetFlagTbinaryModel() == 0)		// "modern" model
						for (size_t i = 0; i < I; i++)	//cp
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[2]*CTPL);
					if (sistema_a->GetFlagInitAproxCT() > 1)	//F
						for (size_t i = 0; i < I; i++)
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7]*CTPL);
					if ( sistema_a->GetFlagModeloRede() == 1)
						for (size_t b = 0; b < B - 1; b++)	//teta
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
						vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+1]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
						vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+2]*CTPL);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
						vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+3]*CTPL);
					if (sistema_a->GetFlagPhmax() == 1)
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4]*CTPL);
					for (size_t r = 0; r < R; r++)	//phg
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3]*CTPL);
					for (size_t r = 0; r < R; r++)	//q
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+1]*CTPL);
					if ( sistema_a->GetFlagModeloRede() > 0)
						for (size_t b = 0; b < B; b++)	//def
							vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+3]*CTPL);
					else
						vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + ii++].set(GRB_DoubleAttr_Obj, ident_var[3+flag7+flag4+flag1a+4+flag3+3]*CTPL);
					//for (int i = 0; i < numero_var_ad2[modelo]; i++)
					//	vars[modelo][numero_var_subp[modelo] + numero_var_ad[modelo] + i].set(GRB_DoubleAttr_Obj, CTPL);
				}
			}
			modelosGRB[modelo]->update();	// Atualiza o modelo Gurobi.

		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during CriarFuncaoObjetivoL()" << endl;	
	}
}
void Hrstc::CriarFuncaoObjetivoPL()
{
	try 
	{
		GRBLinExpr fo;
		//x = [pt u cp F teta ph v d s phmax phg q z def vfol] ou
		//x = [pt u up ud F teta ph v d s phmax phg q z def vfol]
		double coeficiente;
		GRBVar variavel;
		// Termos lineares
		int nt;
		int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		double deltaT;
		int delta = 0;
		int cen = 0;
		int flag1a = 0;
		if (flag1 == 1)
			flag1a = 1;
		int flag1d = 1;		// referente à var. def
		if (flag1 == 0)
			flag1d = 0;
		for (int t = 0; t < T; t++)
		{
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				deltaT = sistema_a->GetDeltaT1();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						coeficiente = 1*deltaT;		//F
						variavel = varPL[i + (3+flag7)*sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					else
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;		//pt
						variavel = varPL[i + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
						if (sistema_a->GetFlagTbinaryModel() == 0)
							coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;	//u
						else
							coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT + sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;		//u
						variavel = varPL[i + sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						coeficiente = 1;														//cp
					else
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		//up
					variavel = varPL[i + 2*sistema_a->termeletricasVtr.size() + delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				deltaT = sistema_a->GetDeltaT2();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						coeficiente = 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//F
						variavel = varPL[i + (3+flag7)*sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					else
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//pt
						variavel = varPL[i + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
						if (sistema_a->GetFlagTbinaryModel() == 0)
							coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
						else
							coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) + sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
						variavel = varPL[i + sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						coeficiente = 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//cp
					else
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoPartida() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//up
					variavel = varPL[i + 2*sistema_a->termeletricasVtr.size() + delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
			}
			delta = delta + nt;
		}
		
		// Deficit e vfolga
		int JJ = 0;
			for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
				JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
		delta = (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ;
		for (int t = 0; t < T; t++)
		{
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				deltaT = sistema_a->GetDeltaT1();
				if (sistema_a->GetFlagModeloRede() > 0)
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					{
						coeficiente = sistema_a->GetCustoDeficit()*deltaT;		//def
						variavel = varPL[delta + b];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
				}
				else
				{
					coeficiente = sistema_a->GetCustoDeficit()*deltaT;		//def
					variavel = varPL[delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				deltaT = sistema_a->GetDeltaT2();
				if (sistema_a->GetFlagModeloRede() > 0)
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					{
						coeficiente = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
						variavel = varPL[delta + b];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
				}
				else
				{
					coeficiente = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
					variavel = varPL[delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
				if (sistema_a->GetFlagVfol() == true)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					{
						if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
						{
							coeficiente = Cvfol*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
							variavel = varPL[delta + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d) + r];
							fo.addTerms(&coeficiente, &variavel, 1);
						}
					}
				}
			}
			delta = delta + nt;
		}

		// Termos quadráticos
		if (sistema_a->GetFlagInitAproxCT() == 0)
			cout << "Heuristica - CriarFuncaoObjetivoPL: Termo quadrático não suportado para essa modelagem!!" << endl;

		modeloPL.setObjective(fo, GRB_MINIMIZE);
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during CriarFuncaoObjetivoPL()" << endl;	
	}
}

void Hrstc::FixarCondIniciais(int modelo)
{
	// Fixa condições iniciais de cada subproblema, quando for o caso
	vetorint vtr_a;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
		vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
	for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
		vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

	int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int c = 0;
	int TUr = 0;
	int t_a = 0;
	if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
	{
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)
		{
			if (vtr_a[vtr_i] >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (vtr_a[vtr_i] - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = vtr_a[vtr_i];
			for (size_t i = 0; i < I; i++)
			{
				// Verifica condições iniciais e fixa valores das variáveis caso necessário
				if (sistema_a->termeletricasVtr[i].GetU0() == 0)
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTDown() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				else
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTUp() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				if ((TUr >= 1) && (t_a + 1 <= TUr))
				{
					//u
					vars[modelo][i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					vars[modelo][i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//up
					vars[modelo][i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[modelo][i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
					//ud
					vars[modelo][i + 3*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[modelo][i + 3*I + c].set(GRB_DoubleAttr_UB, 0);
				}
			}

			if (((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);		// quantidade de variáveis por periodo do problema, portanto é o delta para cada periodo do subp. também
		}
	}
	else
	{
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)
		{
			if (vtr_a[vtr_i] >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (vtr_a[vtr_i] - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = vtr_a[vtr_i];
			for (size_t i = 0; i < I; i++)
			{
				// Verifica condições iniciais e fixa valores das variáveis caso necessário
				if (sistema_a->termeletricasVtr[i].GetU0() == 0)
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTDown() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				else
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTUp() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				if ((TUr >= 1) && (t_a + 1 <= TUr))
				{
					//u
					vars[modelo][i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					vars[modelo][i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//cp
					vars[modelo][i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[modelo][i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
				}
			}
		if (((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
		}
	}
}
void Hrstc::FixarCondIniciaisPL()
{
	int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	size_t I = sistema_a->termeletricasVtr.size();
	int c = 0;
	int TUr = 0;
	int t_a = 0;
	if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
	{
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = t;
			for (size_t i = 0; i < I; i++)
			{
				// Verifica condições iniciais e fixa valores das variáveis caso necessário
				if (sistema_a->termeletricasVtr[i].GetU0() == 0)
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTDown() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				else
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTUp() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				if ((TUr >= 1) && (t_a + 1 <= TUr))
				{
					//u
					varPL[i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					varPL[i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//up
					varPL[i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					varPL[i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
					//ud
					varPL[i + 3*I + c].set(GRB_DoubleAttr_LB, 0);
					varPL[i + 3*I + c].set(GRB_DoubleAttr_UB, 0);
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
			if (t >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = t;
			for (size_t i = 0; i < I; i++)
			{
				// Verifica condições iniciais e fixa valores das variáveis caso necessário
				if (sistema_a->termeletricasVtr[i].GetU0() == 0)
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTDown() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				else
					TUr = max(0, sistema_a->termeletricasVtr[i].GetTUp() + 1 - sistema_a->termeletricasVtr[i].GetX0());
				if ((TUr >= 1) && (t_a + 1 <= TUr))
				{
					//u
					varPL[i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					varPL[i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//cp
					varPL[i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					varPL[i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
				}
			}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
		}
	}
}
void Hrstc::AtualizarVariaveis(int janela)
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
		int flag1a = 0;		// referente à var. teta
		if (flag1 == 1)
			flag1a = 1;
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
				if ( sistema_a->GetFlagTbinaryModel() == 0 )
				{
					for (int i = 0; i < I; i++)
					{
						vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_LB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + 2*I + i]);		//up
						vars[modelo][i + 2*I + delta].set(GRB_DoubleAttr_UB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + 2*I + i]);
						vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_LB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + 3*I + i]);		//ud
						vars[modelo][i + 3*I + delta].set(GRB_DoubleAttr_UB, x_hat[vtr_a[modelo][vtr_i]*nt + flag2*cen*R + 3*I + i]);
					}
				}
				jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R;
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
	}
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during AtualizarVariaveis()" << endl;	
	}
}
void Hrstc::AtualizarRestricoes(int modelo, CMatrizEsparsa * x_spr)
{
	// para t > 0 (modelo > 0) atualizar o RHS dos subproblemas que usam a soluçao anterior como dados de entrada!!
	// matriz A contem as linhas das restrições com acoplamento temporal:
	// -> no caso da modelagem antiga são 6 elementos, em que o índice inicial é dado pelo vetor posicao_rest_din
	// referente as restricoes de Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida (restrições com acoplamento temporal)
	// -> no caso da modelagem nova são 7 elementos
	// referente as restricoes de Balanço hídrico, min. up/down time, extra up/down time, up/down ramp rate e limites min geração (restrições com acoplamento temporal)

	CMatrizEsparsa M;
	int nt_rest;
	double rhs_ac = 0;
	if (sistema_a->GetFlagTbinaryModel() == 0)
	{
		for (size_t i = 0; i < posicao_rest_din[modelo].size(); i++)	// loop para as 6 submatrizes q acoplam no tempo, Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida
		{
			M = A[i];
			M.MultiplicarPorMatriz(x_spr);		// o resultado sera os valores de entrada para o subproblema i, visto que a soluçao para os subp > i sao nulas ainda (desconhecidas)
			nt_rest = M.GetNlin() / T;		// = R (balanço hidrico), Tup, Tdown, I (rampa up), I (rampa down) e I (custo de partida)
			for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais	
				for (int iii = 0; iii < nt_rest; iii++)	
					constr[modelo][iii + ii*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_pri[modelo][ii]*nt_rest + iii][1]) - double(M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0)));
			for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
				for (int iii = 0; iii < nt_rest; iii++)
					constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_fut[modelo][ii]*nt_rest + iii][1]) - double(M.GetElemento(vtr_nos_fut[modelo][ii]*nt_rest + iii, 0)));
		}
	}
	else
	{
		for (size_t i = 0; i < posicao_rest_din[modelo].size(); i++)	// loop para as 7 submatrizes q acoplam no tempo
		{
			M = A[i];
			M.MultiplicarPorMatriz(x_spr);		// o resultado sera os valores de entrada para o subproblema i, visto que a soluçao para os subp > i sao nulas ainda (desconhecidas)
			nt_rest = M.GetNlin() / T;		// = R (balanço hidrico), sum(i) para Tup > 0, sum(i) para Tdown > 0, I Tupdown, I (rampa up), I (rampa down), I (limite de ptmin)
			
			if ( ((i == 1) || (i == 2) || (i == 3)) && (sistema_a->GetFlagVarBin() == 1) )		//restrições em que o rhs é inteiro
			{
				for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais	
					for (int iii = 0; iii < nt_rest; iii++)	
						constr[modelo][iii + ii*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, floor(lim_iniciais[i][vtr_nos_pri[modelo][ii]*nt_rest + iii][1] + 0.5) - floor(M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0) + 0.5));
				for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
					for (int iii = 0; iii < nt_rest; iii++)
						constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, floor(lim_iniciais[i][vtr_nos_fut[modelo][ii]*nt_rest + iii][1] + 0.5) - floor(M.GetElemento(vtr_nos_fut[modelo][ii]*nt_rest + iii, 0) + 0.5));
			}
			else if ((i == 6) && (sistema_a->GetFlagVarBin() == 1))		// rhs dessa restrição não pode dar < 0
			{
				for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais	
					for (int iii = 0; iii < nt_rest; iii++)
					{
						rhs_ac = M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0);
						//cout << modelo << " | " << rhs_ac << endl;
						if ((0.001 > rhs_ac) && (rhs_ac > 0))
							rhs_ac = 0.0;
						else if (rhs_ac >= 0.001)
							cout << "Valor de rhs_ac em AtualizarRestricoes() esta errado!" << " Modelo " << modelo << endl;
						constr[modelo][iii + ii*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_pri[modelo][ii]*nt_rest + iii][1]) - rhs_ac);
					}
				for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
					for (int iii = 0; iii < nt_rest; iii++)
					{
						rhs_ac = M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0);
						if ((0.001 > rhs_ac) && (rhs_ac > 0))
							rhs_ac = 0.0;
						else if (rhs_ac >= 0.001)
							cout << "Valor de rhs_ac em AtualizarRestricoes() esta errado!" << " Modelo " << modelo << endl;
						constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_fut[modelo][ii]*nt_rest + iii][1]) - rhs_ac);
					}
			}
			else
			{
				for (int ii = 0; ii < vtr_nos_pri[modelo].size(); ii++)		// janelas principais	
					for (int iii = 0; iii < nt_rest; iii++)	
						constr[modelo][iii + ii*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_pri[modelo][ii]*nt_rest + iii][1]) - double(M.GetElemento(vtr_nos_pri[modelo][ii]*nt_rest + iii, 0)));
				for (int ii = 0; ii < vtr_nos_fut[modelo].size(); ii++)		// janelas futuro
					for (int iii = 0; iii < nt_rest; iii++)
						constr[modelo][iii + (ii + vtr_nos_pri[modelo].size())*nt_rest + posicao_rest_din[modelo][i]].set(GRB_DoubleAttr_RHS, double (lim_iniciais[i][vtr_nos_fut[modelo][ii]*nt_rest + iii][1]) - double(M.GetElemento(vtr_nos_fut[modelo][ii]*nt_rest + iii, 0)));
			}
		}

	}
	modelosGRB[modelo]->update();	// Atualiza o modelo Gurobi
}

double Hrstc::ResolverSubpLinear(int modelo, CMatrizEsparsa * x_spr)
{
	// GRBmodel *fixed = GRBfixedmodel(modelosGRB[modelo]);
	//Create the fixed model associated with a MIP model. The MIP model must have a solution loaded (e.g., after a call to GRBoptimize). In the fixed model, each integer variable is fixed to the value that variable takes in the MIP solution.
	//This routine returns the computed model. If there is a problem, the routine returns NULL.
	//GRBmodel *fixed = modelosGRB[modelo]->fixedModel();

	// calcular função objetivo de cada subproblema resolvendo um PL com as var. binárias da solução fixas!
	// usar o mesmo modelo do gurobi: fixar var. binarias e deixá-las contínuas e depois reverter alterações
	// comparar com o custo obtido com a resolução do PL completo, deve ser o mesmo!
	
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen, cen_subp, delta_x, delta_subp, jj;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	vetorfloat limites;
	cen_subp = 0;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
	{
		cen = 0;
		if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
			cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		delta_x = nt*vtr_nos_pri[modelo][vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();
		delta_subp = nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size();

		// fixar somente as variáveis binárias
		for (int it = 0; it < I; it++)
		{
			//gravar limites
			limites.push_back(vars[modelo][it + I + delta_subp].get(GRB_DoubleAttr_LB));
			limites.push_back(vars[modelo][it + I + delta_subp].get(GRB_DoubleAttr_UB));
			//u
			vars[modelo][it + I + delta_subp].set(GRB_CharAttr_VType, 'C');
			vars[modelo][it + I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(it + I + delta_x, 0));
			vars[modelo][it + I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(it + I + delta_x, 0));
		}
		jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				//gravar limites
				limites.push_back(vars[modelo][j + jj + delta_subp].get(GRB_DoubleAttr_LB));
				limites.push_back(vars[modelo][j + jj + delta_subp].get(GRB_DoubleAttr_UB));
				//z
				vars[modelo][j + jj + delta_subp].set(GRB_CharAttr_VType, 'C');
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(j + jj + delta_x, 0));
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(j + jj + delta_x, 0));
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		if ( (vtr_nos_pri[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
			cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
	}
	modelosGRB[modelo]->update();

	double funcao_obj = 0;
	try
	{
		modelosGRB[modelo]->reset();
		// Criar funçao objetivo
		///CriarFuncaoObjetivoPLs(modelo);
		// Otimizar
		modelosGRB[modelo]->optimize();

		// Salvar resultados
		if ((modelosGRB[modelo]->get(GRB_IntAttr_Status) == 2) || (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 9))
		{
			AlocarSolucao(modelo, x_spr);			
			funcao_obj = double(modelosGRB[modelo]->get(GRB_DoubleAttr_ObjVal));
		}
		else
		{
			modelosGRB[modelo]->computeIIS();				
			modelosGRB[modelo]->write("Hrst_subprobLP.lp");	
			modelosGRB[modelo]->write("Hrst_subprobLP.ilp");
			// n precisa inserir zeros em x_spr, ele já está iniciado com 0's. Basta n atribuir nenhum valor
			funcao_obj = GRB_INFINITY;
		}

	} catch(GRBException e) {
	cout << "Error code = " << e.getErrorCode() << endl;
	cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during ResolverSubpLinear()" << endl;
	}

	// Reverter variáveis binárias e seus limites
	cen_subp = 0;
	delta_x = 0;		// delta_x usado para percorrer vetor limites, que já está na ordem apropriada
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
	{
		delta_subp = nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size();

		for (int it = 0; it < I; it++)
		{
			//u
			vars[modelo][it + I + delta_subp].set(GRB_CharAttr_VType, 'B');
			vars[modelo][it + I + delta_subp].set(GRB_DoubleAttr_LB, limites[delta_x++]);
			vars[modelo][it + I + delta_subp].set(GRB_DoubleAttr_UB, limites[delta_x++]);
		}
		jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				//z
				vars[modelo][j + jj + delta_subp].set(GRB_CharAttr_VType, 'B');
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_LB, limites[delta_x++]);
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_UB, limites[delta_x++]);
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		if ( (vtr_nos_pri[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
			cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
	}
	modelosGRB[modelo]->update();




	return ( funcao_obj );
}
double Hrstc::ResolverPL(CMatrizEsparsa * x_spr)
{

	////// arquivo para debug;
	////ofstream debugF( "debug.txt", ios::out );
	////if ( debugF.is_open() )
	////{
	////	debugF << std::scientific << setprecision(10);
	////	debugF << "var  /" << "w    /" << "t    /" << "valor " << endl;
	////}
	////else
	////	cout << "Unable to open file";

	// fixar soluções necessarias ou n, resolver para calcular f.o. apenas ou para rodar modelo completo!
	// loop paras todos subproblemas e fixar
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen, cen_subp, delta_x, jj;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int i = 0; i < n_modelos; i++)
	{
		vetorint vtr_a;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[i].size(); vtr_i++)		// referente as var. da janela principal
			vtr_a.push_back(vtr_nos_pri[i][vtr_i]);
		//for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[i].size(); vtr_i++)		// referente as var. da janela futura
		//	vtr_a.push_back(vtr_nos_fut[i][vtr_i]);

		cen_subp = 0;
		for (size_t vtr_i = 0; vtr_i < vtr_a.size(); vtr_i++)
		{
			cen = 0;
			if ( vtr_a[vtr_i] >= sistema_a->GetTt2() )
				cen = (vtr_a[vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			delta_x = nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();

			//// fixar somente as variáveis de decisão (pt, ph)	-> não obedece a restrição de atend. à demanda
			//for (int it = 0; it < I; it++)
			//{
			//	varPL[it + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(it + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//	varPL[it + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(it + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//}
			//delta_x += I*(3 + flag4 + flag7) + flag1a*(B - 1);
			//for (int r = 0; r < R; r++)
			//{
			//	varPL[r + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(r + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//	varPL[r + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(r + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//}
			//delta_x -= I*(3 + flag4 + flag7) + flag1a*(B - 1);


			//// fixar somente as variáveis binárias
			for (int it = 0; it < I; it++)
			{
				//u
				varPL[it + I + delta_x].set(GRB_CharAttr_VType, 'C');
				varPL[it + I + delta_x].set(GRB_DoubleAttr_LB, floor (x_spr->GetElemento(it + I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
				varPL[it + I + delta_x].set(GRB_DoubleAttr_UB, floor (x_spr->GetElemento(it + I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
				//varPL[it + I + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(it + I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
				//varPL[it + I + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(it + I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
				////if ( debugF.is_open() )
				////	debugF << "u:\t" << cen << "\t" << vtr_a[vtr_i] << "\t" << x_spr->GetElemento(it + I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) << "\t" << floor (x_spr->GetElemento(it + I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5) << endl;

			}
			if (sistema_a->GetFlagTbinaryModel() == 1)
			{
				for (int it = 0; it < I; it++)
				{
					//up
					varPL[it + 2*I + delta_x].set(GRB_CharAttr_VType, 'C');
					varPL[it + 2*I + delta_x].set(GRB_DoubleAttr_LB, floor (x_spr->GetElemento(it + 2*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
					varPL[it + 2*I + delta_x].set(GRB_DoubleAttr_UB, floor (x_spr->GetElemento(it + 2*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
					//varPL[it + 2*I + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(it + 2*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
					//varPL[it + 2*I + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(it + 2*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
					////if ( debugF.is_open() )
					////	debugF << "up:\t" << cen << "\t" << vtr_a[vtr_i] << "\t" << x_spr->GetElemento(it + 2*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) << "\t" << floor (x_spr->GetElemento(it + 2*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5) << endl;
					
					//ud
					varPL[it + 3*I + delta_x].set(GRB_CharAttr_VType, 'C');
					varPL[it + 3*I + delta_x].set(GRB_DoubleAttr_LB, floor (x_spr->GetElemento(it + 3*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
					varPL[it + 3*I + delta_x].set(GRB_DoubleAttr_UB, floor (x_spr->GetElemento(it + 3*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
					//varPL[it + 3*I + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(it + 3*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
					//varPL[it + 3*I + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(it + 3*I + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
				}
			}
			jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
			for (int r = 0; r < R; r++)
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					//z
					varPL[j + jj + delta_x].set(GRB_CharAttr_VType, 'C');
					varPL[j + jj + delta_x].set(GRB_DoubleAttr_LB, floor (x_spr->GetElemento(j + jj + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
					varPL[j + jj + delta_x].set(GRB_DoubleAttr_UB, floor (x_spr->GetElemento(j + jj + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0) + 0.5));
					//varPL[j + jj + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(j + jj + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
					//varPL[j + jj + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(j + jj + nt*vtr_a[vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size(), 0));
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}

			// fixar tds variaveis
			//for (int ii = 0; ii < nt; ii++)
			//{
			//	if (varPL[ii + delta_x].get(GRB_CharAttr_VType) == 'B')		// se for variável binária
			//	{
			//		// alterar para se tornar um modelo PL, o problema é como identificar as variáveis...
			//		//varPL[ii + delta_x].set(GRB_CharAttr_VType, 'C');
			//		varPL[ii + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(ii + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//		varPL[ii + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(ii + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//	}
			//	else
			//	{
			//		//varPL[ii + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(ii + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//		//varPL[ii + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(ii + nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//	}
			//}

			//if ( (vtr_a[vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
			//{
			//	for (int r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	// vfol
			//	{
			//		//varPL[r + nt + delta_x].set(GRB_DoubleAttr_LB, x_spr->GetElemento(r + nt*(vtr_i+1) + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//		//varPL[r + nt + delta_x].set(GRB_DoubleAttr_UB, x_spr->GetElemento(r + nt*(vtr_i+1) + flag2*cen_subp*sistema_a->hidreletricasVtr.size(), 0));
			//	}
			//	cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
			//}
		}
	}
	modeloPL.update();


	// arquivo para debug;
	////debugF.close();

	double f_obj;
	try
	{
		// Otimizar
		modeloPL.optimize();

		//modeloPL.write("prob_hrst_DEL.lp");	

		// Salvar resultados
		if (modeloPL.get(GRB_IntAttr_Status) == 2)
		{
			//for (int i = 0; i < n; i++)
			//	x_spr->SubstituirElemento(i, 0, varPL[i].get(GRB_DoubleAttr_X));
			f_obj = double(modeloPL.get(GRB_DoubleAttr_ObjVal));
		}
		else if (modeloPL.get(GRB_IntAttr_Status) == 9)
		{
			modeloPL.getEnv().set(GRB_IntParam_ScaleFlag, 0);
			modeloPL.reset();
			modeloPL.optimize();

			if (modeloPL.get(GRB_IntAttr_Status) == 2)
				f_obj = double(modeloPL.get(GRB_DoubleAttr_ObjVal));
			else
				f_obj = GRB_INFINITY;

			modeloPL.getEnv().set(GRB_IntParam_ScaleFlag, 1);
		}
		else if (modeloPL.get(GRB_IntAttr_Status) == 3)
		{
			//modeloPL.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-4); 
			//modeloPL.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-4);
			//modeloPL.reset();
			//modeloPL.optimize();

			//if (modeloPL.get(GRB_IntAttr_Status) == 2)
			//	f_obj = double(modeloPL.get(GRB_DoubleAttr_ObjVal));
			//else
			//{
				// status = 3, inviável mesmo n é erro numérico...
				f_obj = GRB_INFINITY;
			//}

			//modeloPL.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6); 
			//modeloPL.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-6);
		}
		else
		{
			cout << "Modelo PL da Heuristica status = " << modeloPL.get(GRB_IntAttr_Status) << endl;
			modeloPL.computeIIS();				
			modeloPL.write("probDEL.lp");	
			modeloPL.write("probDEL.ilp");
			for (int i = 0; i < n; i++)
				x_spr->SubstituirElemento(i, 0, 0);
			f_obj = GRB_INFINITY;
		}
		

	} catch(GRBException e) {
	cout << "Error code = " << e.getErrorCode() << endl;
	cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during ResolverPL()" << endl;
	}

	return ( f_obj );
}
void Hrstc::CalcularFuncaoObjetivoRL(vetorfloat & fo)
{
	// Calcular funçao objetivo (janelas principais e futuras) a partir da soluçao da RL x_hat, sem os termos da RL (lambda*g(x))
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	for (int modelo = 0; modelo < n_modelos; modelo++)
	{
		funcao_obj = 0;
		// Criar um vetor com as janelas principais e futuras juntas, para n precisa repetir os loops abaixo
		vetorint vtr_a;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)		// referente as var. da janela principal
			vtr_a.push_back(vtr_nos_pri[modelo][vtr_i]);
		for (size_t vtr_i = 0; vtr_i < vtr_nos_fut[modelo].size(); vtr_i++)		// referente as var. da janela futura
			vtr_a.push_back(vtr_nos_fut[modelo][vtr_i]);

		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
		//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model
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
						funcao_obj += x_hat[i + (3+flag7)*I + delta_x] * 1*deltaT;		//F
					else
					{
						funcao_obj += x_hat[i + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
						if (sistema_a->GetFlagTbinaryModel() == 0)
							funcao_obj += x_hat[i + I + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;		//u
						else
							funcao_obj += x_hat[i + I + delta_x] * (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin());		//u
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						funcao_obj += x_hat[i + 2*I + delta_x] * 1;		//cp
					else
						funcao_obj += x_hat[i + 2*I + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		//up
				}
				if (sistema_a->GetFlagModeloRede() > 0)
					for (int b = 0; b < B; b++)
						funcao_obj += x_hat[b + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT;		//def
				else
					funcao_obj += x_hat[I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT;		//def
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2();
				for (int i = 0; i < I; i++)
				{
					if (sistema_a->GetFlagInitAproxCT() > 1)
						funcao_obj += x_hat[i + (3 + flag7)*I + delta_x] * 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//F
					else
					{
						funcao_obj += x_hat[i + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);			//pt
						if (sistema_a->GetFlagTbinaryModel() == 0)
							funcao_obj += x_hat[i + I + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
						else
							funcao_obj += x_hat[i + I + delta_x] * (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT*sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
					}
					if (sistema_a->GetFlagTbinaryModel() == 0)
						funcao_obj += x_hat[i + 2*I + delta_x] * 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//cp
					else
						funcao_obj += x_hat[i + 2*I + delta_x] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//up
				}
				if (sistema_a->GetFlagModeloRede() > 0)
					for (int b = 0; b < B; b++)
						funcao_obj += x_hat[b + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
				else
					funcao_obj += x_hat[I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta_x] * sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
			
				if ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)		// && ((vtr_a[vtr_i] + 1 - sistema_a->GetTt1()) > 0))
					if (sistema_a->GetFlagVfol() == true)
						for (int r = 0; r < R; r++)
							funcao_obj += x_hat[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta_x] * Cvfol*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
			}
		}

		// Termos quadráticos
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
			if (sistema_a->GetFlagTbinaryModel() == 1)
				cout << "CalcularFuncaoObjetivoRL: Termo quadrático não implementado para essa modelagem!!" << endl;

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
		//if ( abs(funcao_obj) == 0 )
		if ( abs(funcao_obj) < 1e-4 )
			fo[modelo] = 1;
		else
			fo[modelo] = funcao_obj;
	}
}
void Hrstc::AlocarSolucao(int modelo, CMatrizEsparsa * x_spr)
{
	// ao alocar solução observar IntFeasTol para as var. binárias e FeasibilityTol para as continuas (arredondar!)
	// arredondar para duas casa após a vírgula! ?? Ver com o Antonio como deixar dependendo da precisão do solver (truncar valor até a precisão)

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen;
	int R = sistema_a->hidreletricasVtr.size();
	int delta = 0;
	double valor;
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
				{
					//cout << ii << " : " << vars[modelo][ii + delta].get(GRB_DoubleAttr_X) << endl;
					valor = min(max(vars[modelo][ii + delta].get(GRB_DoubleAttr_X),vars[modelo][ii + delta].get(GRB_DoubleAttr_LB)),vars[modelo][ii + delta].get(GRB_DoubleAttr_UB));
					x_spr->SubstituirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(valor > 0 ? double( valor <= 1e-10 ? 0 : valor) : valor));

					//x_spr->SubstituirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(vars[modelo][ii + delta].get(GRB_DoubleAttr_X)));
					//x_spr->SubstituirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(vars[modelo][ii + delta].get(GRB_DoubleAttr_X)*1 + 0.0)/1);
				}
				delta += nt + R;
			}
		}
		else
		{
			for (int ii = 0; ii < nt; ii++)
			{
				//cout << ii << " : " << vars[modelo][ii + delta].get(GRB_DoubleAttr_X) << endl;
				valor = min(max(vars[modelo][ii + delta].get(GRB_DoubleAttr_X),vars[modelo][ii + delta].get(GRB_DoubleAttr_LB)),vars[modelo][ii + delta].get(GRB_DoubleAttr_UB));
				x_spr->SubstituirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(valor > 0 ? double( valor <= 1e-10 ? 0 : valor) : valor));
				
				//x_spr->SubstituirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(vars[modelo][ii + delta].get(GRB_DoubleAttr_X)));
				//x_spr->SubstituirElemento(vtr_nos_pri[modelo][vtr_i]*nt + flag2*cen*R + ii, 0, double(vars[modelo][ii + delta].get(GRB_DoubleAttr_X)*1 + 0.0)/1);
			}
			delta += nt;
		}
	}
	
}

double Hrstc::ResolverHeuristica()
{
	///cout << "ResolverHeuristica() : ";
	
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

	CMatrizEsparsa * X_it_b;
	X_it_b = new CMatrizEsparsa(n, 1);	// x_spr contem a soluçao de tds os subproblemas anteriores subp <= i

	double fo_it = 0;
	double fo_ub = 0;
	double fo_ub_med = 0;
	double fo_lb = 0;
	
	nStatus.clear();
	nStatus.resize(n_modelos);
	Status = 1;
	int iter = 0;
	CalcularFuncaoObjetivoRL(foRL);

	if (APPROACH)		// com benders
	{
		// Primeira iteração, calcula-se cortes para a solução do problema Ext
		if ((itB == 0) && (INITCUT == 1))
		{
			X_it->RemoverTodosElementos();
			//fo_ub = 0;
			//fo_lb = 0;
			ResolverBackward0(X_it);
			
			iter++;
			itB++;
			Status = 2;
		}
		else
		{
			// Receber valores de x_hat e x_til!!!
			if (alfa != 1.0)
			{
				x_hat = resultadosGurobi->GetX_med(0);
				x_til = resultadosGurobi->GetX_med(1);
			}
			else
			{
				x_hat = vetorfloat(n, 0);
				x_til = vetorfloat(n, 0);
			}

			// Demais iterações
			while (iter < itBenders)
			{
				X_it->RemoverTodosElementos();
				fo_ub = 0;

				// Forward
				if (ResolverForward(X_it, &fo_ub))
				{
					cout << "Inviabilidade no Forward!" << endl;
					Status++;
					break;
				}
				fo_lb = fo_subp[0];

				#if (DEBUG)
				//--------------------------------
				// Conferir LB e UB
				if (fo_ub < fo_lb)
					cout << "LB > UB!" << endl;
				// comparar custos individuais!
				// etapa 3, primeiro estágio
				for (size_t ccc = 0; ccc < subpB3etapa.size(); ccc++)
					if (fo_subp[subpB3etapa[ccc]+1] - ff_subp[subpB3etapa[ccc]] < - 0.01)
						cout << "approx. " << subpB3etapa[ccc] << ": " << fo_subp[subpB3etapa[ccc]+1] - ff_subp[subpB3etapa[ccc]] << endl;
				// etapa 2, subp de acoplamento
				double vlr_esp = 0;
				for (size_t ccc = 1; ccc < subpB2etapa.size(); ccc++)
				{
					vlr_esp += fo_subp[subpB2etapa[ccc]];
					if (fo_subp[subpB2etapa[ccc]] - debug_alfa_B2[ccc - 1] < - 0.01)
					{
						cout << "modelo " << subpB2etapa[0] << " aprox. = " << ccc-1 << ": " << fo_subp[subpB2etapa[ccc]] - debug_alfa_B2[ccc - 1] << endl;
					}
				}
				if (vlr_esp - ff_subp[subpB2etapa[0]] < - 0.01 )
					cout << "approx. " << subpB2etapa[0] << vlr_esp - ff_subp[subpB2etapa[0]] << endl;
				// etapa 1, segundo estágio
				int n_supbB1_por_cen = subpB1etapa.size() / sistema_a->GetNCenarios();
				for (size_t ccc = 0; ccc < subpB1etapa.size(); ccc++)
					if ((subpB1etapa[ccc] >= subpB1etapa[0]) && ((subpB1etapa[ccc] - subpB2etapa[0]) % n_supbB1_por_cen == 0))	
					{}// se n for de final de horizonte
					else
						if (fo_subp[subpB1etapa[ccc]+1] - ff_subp[subpB1etapa[ccc]] < - 0.01)
							cout << "approx. " << subpB1etapa[ccc] << ": " << fo_subp[subpB1etapa[ccc]+1] - ff_subp[subpB1etapa[ccc]] << endl;
				//--------------------------------
				#endif

				// Backward
				if (ResolverBackward(X_it))
				{
					cout << "Inviabilidade no Backward!" << endl;
					Status++;
					break;
				}

				if ( abs ((fo_ub - fo_lb) / fo_ub) <= 1e-7)
				{
					cout << "Convergiu! Iter = " << iter << endl;
					break;
				}
				// Se fo_ub não se altear para o processo
				if ( abs ((fo_ub_med - fo_ub) / fo_ub) <= 1e-7)
				{
					cout << "Benders saturado" << endl;
					break;
				}
				// Gravar melhor soluçao durante as iteraçoes da RL
				if ( fo_ub <= fo_ub_a )	// Só grava soluçao se a fo for menor que a atual (menor de todas as iterações de Benders)
				{
					///*X = *X_it;
					*X_it_b = *X_it;
					fo_ub_a = fo_ub;
					Status = 0;
				}
				fo_ub_med = (fo_ub+fo_ub_med)/2;
				iter++;
				itB++;
			}
		}
		if (Status == 0)	// nova solução em X
		{
			if (sistema_a->GetFlagVarBin() == 1)
			{
				///double fo_PL = ResolverPL(X);
				double fo_PL = ResolverPL(X_it_b);
				//double fo_PL = ResolverPL(X_it);		// poderia resolver para todo X_it encontrado após o Benders (pois não quer dizer que fo_ub <= fo_ub_a seja melhor, devido às penalidades...)

				//cout << "fo = " << fo << " ; fo_PL = " << fo_PL << endl;
				if (fo_PL <= fo)
				{
					// solução do forward pode ser inviável, pois adicionam-se variáveis de folga para restrições impossiveis de serem atendidas (de acoplamento)
					// assim pode-se gerar uma solução binária que é inviável no forward (para os subproblemas a solução é viável, mas com um custo elavado devido a penalidade, para criar os cortes)
					// se ao aplicar essa sol. binária no PL ainda for inviável, então ainda n tem-se uma solução primal da Heurística
					// uma inviabilidade no forward deveria ser tratada diferente, parando o processo e incluindo um corte de viabilidade naquele subproblema

					// como a solução obinária n é fixada no backward eu posso n ver essa penalidade por inviabilidade (ou vejo sim?!? pois resolvo no backward o mesmo problema, mesmo v0, só que com CF)
					if (fo_PL <= fo_ub)		// Somente se o PL gerar uma solução melhor que o Forward (pois se PL for inviável a sol é INF) e a solução final é a do forward com var. de folga
					{
						fo = fo_PL;
						// usado o X para passar o melhor X_it dentro de todas iterações do benders para o PL!!!
						X->RemoverTodosElementos();
						for (int i = 0; i < n; i++)
							X->SubstituirElemento(i, 0, varPL[i].get(GRB_DoubleAttr_X));
					}
				}
			}
			else
			{
				fo = fo_ub_a;
				///*X = *X_it;
				*X = *X_it_b;
			}
		}
		// Imprimir log
		if ( log_auxiliar->is_open() )
		{
			// Ajustes no .h
			*log_auxiliar << char(9) << "Heuristic:" << char(9) << itB << char(9) << fo_lb << char(9) << fo_ub << char(9) << fo_ub_a << char(9) << fo;
			///*log_auxiliar << endl;
		}
		else
			cout << "Unable to open file";
	}
	else		// sem benders
	{
		// Receber valores de x_hat e x_til!!!
		if (alfa != 1.0)
		{
			x_hat = resultadosGurobi->GetX_med(0);
			x_til = resultadosGurobi->GetX_med(1);
		}
		else
		{
			x_hat = vetorfloat(n, 0);
			x_til = vetorfloat(n, 0);
		}

		X_it->RemoverTodosElementos();

		if (ResolverForward(X_it, &fo_ub))
		{
			cout << "Inviabilidade no Forward!" << endl;
			Status++;
		}

		//cout << fo_ub << " | " << fo_ub_a << " | " << fo << endl;
		// Gravar melhor soluçao durante as iteraçoes da RL
		if ( fo_ub <= fo_ub_a )	// Só grava soluçao se todos os subproblemas forem resolvidos
		{
			fo_ub_a = fo_ub;
			Status = 0;
		}
		if (Status == 0)	// nova solução em X
		{
			double fo_PL = ResolverPL(X_it);
			if (fo_PL <= fo)
			{
				fo = fo_PL;
				X->RemoverTodosElementos();
				for (int i = 0; i < n; i++)
					X->SubstituirElemento(i, 0, varPL[i].get(GRB_DoubleAttr_X));
			}
			else
				cout << "situação n esperada, X foi substituido!" << endl;
		}
		// Imprimir log
		if ( log_auxiliar->is_open() )
		{
			if ((itB == 0) && (INITCUT == 1))	// Primeira iteração, calcula-se cortes para a solução do problema Ext
			{
				*log_auxiliar << char(9) << "Heuristic:" << char(9) << itB << char(9) << "Cortes do modelo Ext.";
				///*log_auxiliar << endl;
			}
			else
			{
				// Ajustes no .h
				*log_auxiliar << char(9) << "Heuristic:" << char(9) << fo_ub << char(9) << fo_ub_a << char(9) << fo;
				///*log_auxiliar << endl;
			}
		}
		else
			cout << "Unable to open file";
	}

	delete X_it_b;
	delete X_it;
	return fo;
}
double Hrstc::ResolverHeuristica(bool so_heuristica)
{
	// Termo proximal é a melhor solução anterior!
	//alfa = 1.0;

	x_hat = vetorfloat(n, 0);
	x_til = vetorfloat(n, 0);
	for (int i = 0; i < n; i++)
	{
		x_hat[i] = X->GetElemento(i, 0);
		x_til[i] = X->GetElemento(i, 0);
	}

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
	//CMatrizEsparsa X_it(n, 1);
	double fo_it = 0;
	double fo_ub = 0;
	double fo_lb = 0;
	
	nStatus.clear();
	nStatus.resize(n_modelos);
	Status = 1;
	int iter = 0;
	CalcularFuncaoObjetivoRL(foRL);

	while (iter < itBendersF)
	{
		alfa *= 1.01;
		if (alfa >= 1)
			alfa = 1;
		cout << alfa << endl;

		X_it->RemoverTodosElementos();
		fo_ub = 0;

		// Forward
		if (ResolverForward(X_it, &fo_ub))
		{
			cout << "Inviabilidade no Forward!" << endl;
			Status++;
			break;
		}
		//cout << endl;
		fo_lb = fo_subp[0];

		// Backward
		if (ResolverBackward(X_it))
		{
			cout << "Inviabilidade no Backward!" << endl;
			Status++;
			break;
		}

		// Imprimir log
		if ( log_auxiliar->is_open() )
		{
			// Ajustes no .h
			*log_auxiliar << char(9) << "Heuristic:" << char(9) << itB << char(9) << fo_lb << char(9) << fo_ub << char(9) << fo_ub_a << char(9) << fo;
			///*log_auxiliar << endl;
		}
		else
			cout << "Unable to open file";
		
		cout << itB << ":: " << setprecision(15) << fo_lb << " | " << fo_ub << " | " << fo_ub_a << " | " << fo << endl;
		if ( abs ((fo_ub - fo_lb) / fo_ub) <= 1e-7)
		{
			cout << "Convergiu! Iter = " << iter << endl;
			//cout << "ub = " << setprecision(15) << fo_ub << endl;
			//cout << "lb = " << setprecision(15) << fo_lb << endl;
			break;
		}
		// Se fo_ub não se altear para o processo
		if ( abs ((fo_ub_a - fo_ub) / fo_ub) <= 1e-7)
		{
			cout << "Benders saturado" << endl;
			//cout << "ub = " << setprecision(15) << fo_ub << endl;
			//cout << "lb = " << setprecision(15) << fo_lb << endl;
			break;
		}

		// testar primeiro o estocastico 46b continuo (ok, demora muito nas últimas iterações (cutting planes) porém tds valores coerentes)
		// testar as var. binarias como complicadas: fixas e atualizando o rhs_bin!!


		//fo = fo_ub;
		// Gravar melhor soluçao durante as iteraçoes da RL
		if ( fo_ub <= fo_ub_a )	// Só grava soluçao se todos os subproblemas forem resolvidos
		{
			//*X = *X_it;
			fo_ub_a = fo_ub;
			Status = 0;
			for (int i = 0; i < n; i++)
			{
				x_hat[i] = X->GetElemento(i, 0);
				x_til[i] = X->GetElemento(i, 0);
			}
		}
		//*X = *X_it;
		iter++;
		itB++;
	}
	if (Status == 0)	// nova solução em X
	{
		double fo_PL = ResolverPL(X_it);
		if (fo_PL <= fo)
		{
			fo = fo_PL;
			for (int i = 0; i < n; i++)
				X->SubstituirElemento(i, 0, varPL[i].get(GRB_DoubleAttr_X));
		}
	}
	cout << fo << endl;
	delete X_it;

	return fo;
}

void Hrstc::EscreverX()
{
	// imprimir x_hat e x_til para to iteração (duas colunas por iteração)
	// comparar com solução impressa x* (do ED)
	// tb testar rodar heurisitca com x*

	x_hat = resultadosGurobi->GetX_med(0);
	//resultadosGurobi->ExportarXmed("x_hat.txt");
	x_til = resultadosGurobi->GetX_med(1);

	ofstream * inFile;
	//inFile = new ofstream( "x_hrstc.txt", ios::app );
	inFile = new ofstream( "x_hrstc.txt", ios::out );
	if ( inFile->is_open() )                                                            
	{
		*inFile << std::scientific << setprecision(10);
		for (size_t i = 0; i < x_hat.size(); i++)
		{
			*inFile << x_hat[i] << char(9) << x_til[i];
			*inFile << endl;
		}
		*inFile << endl;
	}
	else
		cout << "Unable to open file";
	inFile->close();
	delete inFile;

	//// Escrever x_hat e x_til
	//ofstream * inFile;
	//inFile = new ofstream( "x_hat_til.txt", ios::out );
	//if ( inFile->is_open() )                                                            
	//{
	//	*inFile << std::scientific << setprecision(10);
	//	for (size_t i = 0; i < n; i++)
	//		*inFile << x_hat[i] << char(9) << x_til[i] << endl;
	//}
	//else
	//	cout << "Unable to open file";
	//inFile->close();

	//// Ler solução x
	//ifstream inFile( "x.txt", ios::in );   
	//if ( !inFile )                                                            
	//	cout << "File x.txt could not be opened" << endl;
	//int j = 0;
	//while ( ! inFile.eof() )
	//	inFile >> x_hat[j++];
	//inFile.close();
	//x_til = x_hat;

	//// calcular fo de x_til e x_hat
	//for (int i = 0; i < n; i++)
	//	X_it->SubstituirElemento(i, 0, x_til[i]);
	//double fo_xtil = 0;
	//for (int i = 0; i < n_modelos; i++)
	//	fo_xtil+= CalcularFuncaoObjetivo(i, X_it);
	//cout << fo_xtil << " | ";
	//
	//for (int i = 0; i < n; i++)
	//	X_it->SubstituirElemento(i, 0, x_hat[i]);
	//fo_xtil = 0;
	//for (int i = 0; i < n_modelos; i++)
	//	fo_xtil+= CalcularFuncaoObjetivo(i, X_it);
	//cout << fo_xtil << endl;

	//// Escrever x_spr
	//ofstream * inFile;
	//inFile = new ofstream( "x_spr.txt", ios::out );
	//if ( inFile->is_open() )                                                            
	//{
	//	*inFile << std::scientific << setprecision(10);
	//	for (size_t i = 0; i < n; i++)
	//	{
	//		*inFile << X_it->GetElemento(i, 0);
	//		*inFile << endl;
	//	}
	//	*inFile << endl;
	//}
	//else
	//	cout << "Unable to open file";
	//inFile->close();
}
// Benders //
int Hrstc::ResolverSubp(int modelo, bool backward)
{
	try
	{
		// Reset the model to an unsolved state, discarding any previously computed solution information.
		modelosGRB[modelo]->reset();
		
		// Criar funçao objetivo
		CriarFuncaoObjetivo(modelo, backward);

		//if (modelo == 0)
		//	modelosGRB[modelo]->write("Hrst_subprob2.lp");
			
		// resolver subproblema i
		modelosGRB[modelo]->optimize();

		//nStatus[i] = modelosGRB[i]->get(GRB_IntAttr_Status);
		//if ((nStatus[i] == 2) || (nStatus[i] == 9))		// OPTIMAL or TIME_LIMIT
		if (modelosGRB[modelo]->get(GRB_IntAttr_SolCount))		// Number of solutions found during the most recent optimization, se > 0 fazer...
			return (0);
		else
		{
			if (APPROACH)
			{
				//if ((itB == 3) && (modelo == 46) && (backward == true))
				//{
				//	cout << modelosGRB[modelo]->get(GRB_IntAttr_NumVars) << endl;
				//	modelosGRB[modelo]->write("Hrst_subprob_antes.lp");
				//}


				if ((modelosGRB[modelo]->get(GRB_IntAttr_Status) == 9) || (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 4) || (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 12))
				{
					#if (DEBUG)
						cout << "Status code = " << modelosGRB[modelo]->get(GRB_IntAttr_Status) << " modelo = "<< modelo << " backward = " << backward << "itB = " << itB << endl;
						modelosGRB[modelo]->write("Hrst_subprob_N.lp");	
					#endif
					modelosGRB[modelo]->reset();
					modelosGRB[modelo]->getEnv().set(GRB_IntParam_ScaleFlag, 0);
					modelosGRB[modelo]->update();
					modelosGRB[modelo]->optimize();
					modelosGRB[modelo]->getEnv().set(GRB_IntParam_ScaleFlag, 1);
					if (modelosGRB[modelo]->get(GRB_IntAttr_SolCount) > 0)		// Number of solutions found during the most recent optimization, se > 0 fazer...
						return (0);
				}
				////modelosGRB[modelo]->reset();
				////modelosGRB[modelo]->getEnv().set(GRB_IntParam_ScaleFlag, 0);
				////modelosGRB[modelo]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9);
				////modelosGRB[modelo]->update();
				////modelosGRB[modelo]->optimize();
				////if (modelosGRB[modelo]->get(GRB_IntAttr_SolCount) > 0)		// Number of solutions found during the most recent optimization, se > 0 fazer...
				////{
				////	modelosGRB[modelo]->getEnv().set(GRB_IntParam_ScaleFlag, 1);
				////	modelosGRB[modelo]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6);
				////	return (0);
				////}

				//This can happen when a model is right on the edge of feasibility.  One simplex basis might give a dual unbounded ray that barely satisfies tolerances, while another might give an optimal solution that barely satisfies tolerance.  IIS solves a model that isn't identical to your original model, so it can get a slightly different result.
				//I think you'll find that loosening the FeasibilityTol will gives a feasible solution to your original model, and tightening it will give an IIS.
				// https://groups.google.com/forum/#!topic/gurobi/Vg3i8vLqLW8
				int cont_inv = 0;
				while (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 3)		// enquanto o modelo for inviavel
				{
					////if (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 12)
					////	break;
					AdicionarVarFol(modelo);
					modelosGRB[modelo]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6);		// adicionado pela questão acima
					modelosGRB[modelo]->reset();
					modelosGRB[modelo]->optimize();
					cont_inv++;
					if (cont_inv > 100)		// limite para adicionar variaveis de folga
						break;
				}
				// salvar a solução primeiro depois remover as variáveis de folga
				if (modelosGRB[modelo]->get(GRB_IntAttr_SolCount) > 0)
				{
					#if (DEBUG)
						cout << "|=Resolvida com variavel de folga." << backward << endl;
					#endif
					////modelosGRB[modelo]->getEnv().set(GRB_IntParam_ScaleFlag, 1);
					////modelosGRB[modelo]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6);
					return (0);
				}
				if (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 5)	// UNBOUNDED
				{
					modelosGRB[modelo]->reset();
					////modelosGRB[modelo]->getEnv().set(GRB_IntParam_ScaleFlag, 1);
					////modelosGRB[modelo]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-6);
					modelosGRB[modelo]->update();
					modelosGRB[modelo]->optimize();
					if (modelosGRB[modelo]->get(GRB_IntAttr_SolCount) > 0)
						return (0);
				}
				else
				{
					// Se mesmo com variável de folga e sem scale flag não for resolvido
					cout << "Inviabilidade " << modelosGRB[modelo]->get(GRB_IntAttr_Status) << " no subproblema " << modelo << " nao resolvida. " << backward << " itB = " << itB << endl;
					//cout << "enter para gravar arquivos... " << endl;
					//std::cin.ignore();
					//modelosGRB[modelo]->computeIIS();				
					//modelosGRB[modelo]->write("Hrst_subprob.lp");	
					//modelosGRB[modelo]->write("Hrst_subprob.ilp");
					//cout << "enter para continuar... " << endl;
					//std::cin.ignore();

					// anular ultimo corte adicionado, para evitar outros problemas
					int last_const = modelosGRB[modelo]->get(GRB_IntAttr_NumConstrs) - 1;
					constr[modelo][last_const].set(GRB_DoubleAttr_RHS, 0);
					modelosGRB[modelo]->update();
					modelosGRB[modelo]->optimize();
					if (modelosGRB[modelo]->get(GRB_IntAttr_SolCount) > 0)
						return (0);

					//modelosGRB[modelo]->write("Hrst_subprob.lp");
					//modelosGRB[modelo]->computeIIS();
					//modelosGRB[modelo]->write("Hrst_subprob.ilp");
					//if (modelosGRB[modelo]->get(GRB_IntAttr_Status) == 3)
					//{
					//	modelosGRB[modelo]->computeIIS();				
					//	modelosGRB[modelo]->write("Hrst_subprob.lp");	
					//	modelosGRB[modelo]->write("Hrst_subprob.ilp");
					//	std::cin.ignore();

					//	// Nesses casos que da inviabilidade o coeficiente das var. binárias no corte está muito alto, na ordem de 10^12
					//	// Acredito que pq ele venha se acumulando devido a restrição de min. up/down time
					//	// Uma alternativa poderia ser construir o corte diferente, sem alfa (como um corte de viabilidade)
					//	// ou
					//	// colocar um limite máximo para o coeficiente do corte, Se > (1e6 * RHS) , = 9.99e+10
					//	// ou o máximo ser 9.99e+10, pois o alfa é de ordem 10^0.
					//}
					return (1);
				}
			}
			else
			{
				cout << "Inviabilidade " << modelosGRB[modelo]->get(GRB_IntAttr_Status) << " no subproblema " << modelo;
				modelosGRB[modelo]->computeIIS();				
				modelosGRB[modelo]->write("Hrst_subprob.lp");	
				modelosGRB[modelo]->write("Hrst_subprob.ilp");
				return (1);
			}
			return (1);
		}
	//cout << "Fim da heuristica" << endl;
	} catch(GRBException e) {
		modelosGRB[modelo]->write("Hrst_subprob.lp");
		cout << "Excep.: modelo : " << modelo << " Status : " << modelosGRB[modelo]->get(GRB_IntAttr_Status) << endl;
		cout << "ResolverSubp() Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		return (1);
	} catch(...) {
		cout << "Exception during ResolverSubp()" << endl;
		return (1);
	}
}
int Hrstc::ResolverForward(CMatrizEsparsa * x_spr, double * fo_a)
{
	CMatrizEsparsa x_op(n, 1);
	int subp_status = 0;
	for (int i = 0; i < n_modelos; i++)
	{
		if (APPROACH)
		{
			//  Atualizar rhs do corte com as variáveis que não são do subproblema atual!!
			// no subp0 as unicas var. q n são dele são as var. binárias complicadas dos subp. futuros
			int n_supbB1_por_cen = subpB1etapa.size() / sistema_a->GetNCenarios();		// numero de subproblemas por cenário da 1ª etapa
			if ((i >= subpB1etapa[0]) && ((i - subpB2etapa[0]) % n_supbB1_por_cen == 0))
			{}	// subproblema de final de horizonte, nao tem cortes
			else if (i > 0)
				AtualizarCortes(x_op, i);
		}

		// subp. de final de horizonte n sao atualizados!! n tem cortes!!
		//cout << i;
		
		////int nnvar = modelosGRB[i]->get(GRB_IntAttr_NumVars);

		subp_status = ResolverSubp(i);
		//if ((i == 29 || i == 30) && (itB == 3))
		//	ImprimirSolSub(i);

		//cout << "..";
		if ( subp_status == 1)		// Modelo inviavel
		{
			// n precisa inserir zeros em x_spr, ele já está iniciado com 0's. Basta n atribuir nenhum valor
			fo_subp[i] = GRB_INFINITY;
			Status++;
			return ( 1 );
		}
		else
		{
			// Alocar solucao do subp em x_spr (com a soluçao somente das janelas principais)
			AlocarSolucao(i, x_spr);		// aprox. da FCF n é alocada em x_spr

			// alocar solução (deste modelo) em x_op
			for (int ii = nvar_acum[i]; ii < nvar_acum[i + 1]; ii++)
				x_op.SubstituirElemento(ii, 0, x_spr->GetElemento(ii, 0));
			//x_op->InserirMatriz(nvar_acum[i], 0, nvar_acum[i + 1] - 1, 0, x_spr, nvar_acum[i], 0);

			if (APPROACH)
			{
				// Alocar aproximação do custo futuro
				if (i == subpB2etapa[0])
				{
					double valor = 0;
					for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					{
						valor += double(vars[i][numero_var_subp[i] + numero_var_ad[i] + numero_var_ad2[i] + n_var_folga[i] + (sistema_a->GetNCenarios() - 1 - cen)].get(GRB_DoubleAttr_X));
						debug_alfa_B2[cen] = double(vars[i][numero_var_subp[i] + numero_var_ad[i] + numero_var_ad2[i] + n_var_folga[i] + (sistema_a->GetNCenarios() - 1 - cen)].get(GRB_DoubleAttr_X));
						// os cortes são adicionados na ordem inversa, portanto a var. referente ao cen 0 é a ultima adicionada no subproblema (assim como os cortes)
					}
					ff_subp[i] = valor;		// Zsup_i = f.obj_i - alfa (n precisa multiplicar pela probabilidade, pois ela já foi considerada em cada subproblema)
				}
				else
					ff_subp[i] = double(vars[i][numero_var_subp[i] + numero_var_ad[i] + numero_var_ad2[i] + n_var_folga[i]].get(GRB_DoubleAttr_X));
			
				// Alocar função objetivo, essa função objetivo tem os termos proximais e as janelas futuro
				if (alfa != 1)
					fo_subp[i] = CalcularFuncaoObjetivo(i, x_spr) + ff_subp[i];
				else
					fo_subp[i] = double(modelosGRB[i]->get(GRB_DoubleAttr_ObjVal));
			
				// remover variaveis de folga quando existirem
				if (ind_var_fol[i].size() != 0)
					RemoverVarFol(i);
			}
			else
			{
				fo_subp[i] = double(modelosGRB[i]->get(GRB_DoubleAttr_ObjVal));
			}
			
			//cout << setprecision(15) << i << " | " << CalcularFuncaoObjetivo(i, x_spr) + ff_subp[i] << " | " << double(modelosGRB[i]->get(GRB_DoubleAttr_ObjVal)) << endl;
			//if ( CalcularFuncaoObjetivo(i, x_spr) + ff_subp[i] > 1.0000001*double(modelosGRB[i]->get(GRB_DoubleAttr_ObjVal)) || CalcularFuncaoObjetivo(i, x_spr) + ff_subp[i] < 0.9999999*double(modelosGRB[i]->get(GRB_DoubleAttr_ObjVal)))
			//	cout << "diferença subp." << i << endl;
					
			*fo_a += fo_subp[i] - ff_subp[i];
			//cout << "Subproblema " << i << ": " << fo_subp[i] << endl;
			//cout << "Subp. " << i << " [u10, up10, ud10] = [" << vars[i][23].get(GRB_DoubleAttr_X) << ", " << vars[i][37].get(GRB_DoubleAttr_X) << ", " << vars[i][51].get(GRB_DoubleAttr_X) << "]" << endl;
			//*fo_a += fo_subp[i] - ff_subp[i];
		}
		if ( i < n_modelos - 1)
			AtualizarRestricoes(i + 1, x_spr);		// Para subp > 0 antes de cada subproblema se resolvido o lado direito das restriçoes são atualizados!!
	}

	return ( 0 );
	//// Escrever x_spr
	//ofstream * inFile;
	//inFile = new ofstream( "x_spr_PL.txt", ios::out );
	//if ( inFile->is_open() )                                                            
	//{
	//	*inFile << std::scientific << setprecision(10);
	//	for (size_t i = 0; i < n; i++)
	//	{
	//		*inFile << x_spr->GetElemento(i, 0);
	//		*inFile << endl;
	//	}
	//	*inFile << endl;
	//}
	//else
	//	cout << "Unable to open file";
	//inFile->close();
}
int Hrstc::ResolverBackward(CMatrizEsparsa * x_spr)
{
	// Quando consideram-se os cortes de Benders as janelas futuras perdem o sentido, portanto não são consideradas!
	double alfa_0 = alfa;
	alfa = 1.0;		// no backward não são considerados os termos proximais
	int subp_status = 0;
	int mod_ant;
	int n_supbB1_por_cen = subpB1etapa.size() / sistema_a->GetNCenarios();		// numero de subproblemas por cenário da 1ª etapa

	GRBLinExpr restricao;
	double coeficiente;
	GRBVar variavel;
	double rhs_cte, rhs_var, rhs_bin;
	
	vetorfloat RHS_acop;

	// Etapa 1,2 e 3 (problemas de segundo estágio, problema de acoplamento e problemas de primeiro estágio)
	for (int mod = n_modelos - 1; mod > 0; mod--)
	{
		vetorint lim_var_bin;

		// Tornar subpr MILP em LP
		if (sistema_a->GetFlagVarBin() == 1)
		{
			//if (itB > 70)		// talvez colocar para ativar se o LB ficar 0 por muito tempo!!!
			//	SubpMILP2LPFixing(x_spr, mod, &lim_var_bin);
			//else
				SubpMILP2LP(x_spr, mod, &lim_var_bin);
		}

		// resolver subproblema
		subp_status = ResolverSubp(mod, true);
		//if ((mod == 29 || mod == 30) && (itB == 3))
		//	ImprimirSolSub(mod);

		if ( subp_status == 1)
		{
			cout << "Inviabilidade no backward, n devia ocorrer!!!" << endl;
			// Tornar subpr LP em MILP
			SubpLP2MILP(x_spr, mod, &lim_var_bin);
			alfa = alfa_0;

			return ( 1 );
		}
		else
		{
			// Alocar função objetivo, essa função objetivo NÃO tem os termos proximais e as janelas futuro
			fo_subp[mod] = double(modelosGRB[mod]->get(GRB_DoubleAttr_ObjVal));
			// Alocar aproximação do custo futuro
			if (mod == subpB2etapa[0])
			{
				double valor = 0;
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					valor += double(vars[mod][cen + numero_var_subp[mod] + numero_var_ad[mod] + numero_var_ad2[mod] + n_var_folga[mod]].get(GRB_DoubleAttr_X));
				ff_subp[mod] = valor;
			}
			else
				ff_subp[mod] = double(vars[mod][numero_var_subp[mod] + numero_var_ad[mod] + numero_var_ad2[mod] + n_var_folga[mod]].get(GRB_DoubleAttr_X));

			// Identificar subproblema anterior
			if ((mod >= subpB1etapa[0]) && ((mod - (subpB2etapa[0]+1)) % n_supbB1_por_cen == 0))		// se este modelo for de inicio do 2 estágio, acoplamento com o subproblema que é da segunda etapa
				mod_ant = subpB2etapa[0];
			else	// acoplamento com o subproblema anterior dentro do segundo ou do primeiro estágio
				mod_ant = mod - 1;

			// armazenar coeficientes de cada variável para cada subproblema (para todo acoplamento do subproblema mod)
			CMatrizEsparsa L_it(1, n);
			for (size_t ii = 0; ii < acomplamentos[mod - 1].size(); ii++)		// loop no numero de acoplamentos (L_it recebe, para cada variável, o mult. da restrição de acoplamento e o coeficiente da variável)
				L_it.SubstituirElemento(0, int(acomplamentos[mod - 1][ii][1]), L_it.GetElemento(0, int(acomplamentos[mod - 1][ii][1])) + constr[mod][int(acomplamentos[mod - 1][ii][0])].get(GRB_DoubleAttr_Pi) * acomplamentos[mod - 1][ii][3]);
			
			//// -------------------------------------------------------------------
			//// Contabilizar var. bin. do subproblema no corte (só são somadas no rhs variável e na atualização dos cortes)
			//// No caso das var. binárias do subproblema atual o multiplicador é dado pelo - (negativo do) reduced cost delas!!
			//for (size_t ii = 0; ii < index_var_bin[mod].size(); ii++)
			//{
			//	if ( x_spr->GetElemento(index_var_bin[mod][ii] + nvar_acum[mod], 0) >= 0.5)		// se estiver no max = - RC
			//		L_it.SubstituirElemento(0, index_var_bin[mod][ii] + nvar_acum[mod_ant], L_it.GetElemento(0, index_var_bin[mod][ii] + nvar_acum[mod]) - vars[mod][index_var_bin[mod][ii]].get(GRB_DoubleAttr_RC));
			//	else		// se estiver no min = RC
			//		L_it.SubstituirElemento(0, index_var_bin[mod][ii] + nvar_acum[mod_ant], L_it.GetElemento(0, index_var_bin[mod][ii] + nvar_acum[mod]) + vars[mod][index_var_bin[mod][ii]].get(GRB_DoubleAttr_RC));
			//}
			//// -------------------------------------------------------------------

			// atribuir valores para o vetor de coeficientes acumulados
			//double mult_corte = 0;
			if ( (mod >= subpB1etapa[0]) && ((mod - subpB2etapa[0]) % n_supbB1_por_cen == 0) )		// subproblema de final de horizonte
			{}
			else if (mod == subpB2etapa[0])		// subproblema de acoplamento (com várias realizações)
			{
				// multiplicador de todos os cortes do subp. mod (de cada cenário)
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				{
					CMatrizEsparsa m_corte(1, itB + 1);
					for (int iter = 0; iter < itB + 1; iter++)		
						m_corte.SubstituirElemento(0, iter, constr[mod][cen + iter*sistema_a->GetNCenarios() + nrest_orig[mod]].get(GRB_DoubleAttr_Pi)); // * sistema_a->hidreletricasVtr[0].GetProbAfluencia(sistema_a->GetNCenarios() - 1 - cen));
					if (L[subpB2etapa[(sistema_a->GetNCenarios() - 1 - cen) + 1] - 1].GetNnz() != 0)
					{
						m_corte.MultiplicarPorMatriz(&L[subpB2etapa[(sistema_a->GetNCenarios() - 1 - cen) + 1] - 1]);
						L_it.SomarComMatriz(&m_corte);
					}
				}
			}
			else
			{
				// multiplicador de todos os cortes do subp. mod
				CMatrizEsparsa m_corte(1, itB + 1);
				for (int iter = 0; iter < itB + 1; iter++)		
					m_corte.SubstituirElemento(0, iter, constr[mod][iter + nrest_orig[mod]].get(GRB_DoubleAttr_Pi));		// receber o multiplicador do corte
				if (L[mod].GetNnz() != 0)
				{
					// L é de tamanho (iter x n) e m_corte é (1 x iter); L armazena os coeficientes das variáveis anteriores os subproblema t-1 (que vai receber o corte corrente)
					m_corte.MultiplicarPorMatriz(&L[mod]);		// L[mod_apr - 1], que nesses casos = L[mod]
					L_it.SomarComMatriz(&m_corte);
				}
			}
			// Construir corte para o subproblema anterior (mod - 1)
			rhs_cte = 0;
			rhs_var = 0;
			rhs_bin = 0;
			// incluir variável de aproximação do CF
			coeficiente = 1;
			if ((mod >= subpB1etapa[0]) && ((mod - (subpB2etapa[0]+1)) % n_supbB1_por_cen == 0))
			{
				int cen = (mod - (subpB2etapa[0]+1)) / n_supbB1_por_cen;
				variavel = vars[mod_ant][numero_var_subp[mod_ant] + numero_var_ad[mod_ant] + numero_var_ad2[mod_ant] + n_var_folga[mod_ant] + (sistema_a->GetNCenarios() - 1 - cen)];		// ultimas variáveis
			}
			else
				variavel = vars[mod_ant][numero_var_subp[mod_ant] + numero_var_ad[mod_ant] + numero_var_ad2[mod_ant] + n_var_folga[mod_ant]];		// ultima variável
			restricao.addTerms( &coeficiente, &variavel, 1);
			rhs_cte += fo_subp[mod];

			// Como o pi acumulado já foi guardado em L, pode-se somar com o PI do subproblema para montar o corte!
			//PI.SomarComMatriz(&L_it_subp);

			// Percorrer a matriz PI para adicionar cortes (ou L da iteração corrente)!!!
			CMatrizEsparsa L_ac(1, n);
			int l;
			// Calcular o RHS primeiro para avaliar a ordem dos coeficientes
			l = L_it.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
			while ( l != -1 )	// percorre todos elementos nnz da linha
			{
				if ( (L_it.GetValorCol(l) >= nvar_acum[mod_ant]) && (L_it.GetValorCol(l) < nvar_acum[mod_ant + 1]) )		// se existir acoplamento com o subproblema imediatamente anterior
					rhs_cte += L_it.GetValorVal(l) * x_spr->GetElemento(L_it.GetValorCol(l), 0);							// rhs_cte só leva em conta acoplamento com subp. imediatamente anterior
				l = L_it.GetValorLprox(l);		
			}
			l = L_it.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
			while ( l != -1 )	// percorre todos elementos nnz da linha
			{
				if ( (L_it.GetValorCol(l) >= nvar_acum[mod_ant]) && (L_it.GetValorCol(l) < nvar_acum[mod_ant + 1]) )		// se existir acoplamento com o subproblema imediatamente anterior
				{
					// Conferir ordem dos coeficientes, se for menor que coeftol não adiciona-se o termo
					if (abs(L_it.GetValorVal(l) * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB)) >= coeftol * rhs_cte)
					///if ( (abs(L_it.GetValorVal(l) * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB)) >= coeftol * rhs_cte) && (abs(L_it.GetValorVal(l) * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB)) <= rhs_cte / coeftol) )
					{
						coeficiente = L_it.GetValorVal(l);
						variavel = vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]];		// endereço deve ser referente à variável do modelo
						restricao.addTerms( &coeficiente, &variavel, 1);
						//rhs_cte += coeficiente * x_spr->GetElemento(L_it.GetValorCol(l), 0);
					}
				}
				else
				{
					rhs_var += L_it.GetValorVal(l) * x_spr->GetElemento(L_it.GetValorCol(l), 0);
					L_ac.InserirElemento(0, L_it.GetValorCol(l), L_it.GetValorVal(l));
				}
				l = L_it.GetValorLprox(l);		
			}

			// vetor de coeficientes para cada subproblema é o valor do PI para aquela iteração!
			L[mod - 1].JuntarColuna(&L_ac);			// matriz do subproblema corrente, com os coeficientes da aproximação do subproblema anterior a este
			// Armazenar todos os coeficientes menos do subp. anterior
			rhs_cte += rhs_var;

			rhs_orig[mod - 1].push_back(rhs_cte);		// para atualizar o corte necessita-se do rhs_original
			
			// remover variaveis de folga quando existirem
			if (ind_var_fol[mod].size() != 0)
				RemoverVarFol(mod);

			modelosGRB[mod_ant]->addConstr(restricao, GRB_GREATER_EQUAL, rhs_cte - rhs_var, "");			// subtrair rhs_var, pois tem-se o lhs_constante que seria igual ao rhs_var.
			modelosGRB[mod_ant]->update();
			delete constr[mod_ant];
			constr[mod_ant] = modelosGRB[mod_ant]->getConstrs();
			// ao adicionar restrições ao subproblema da 2 etapa, a ordem é de n_cen até 0!!

			// Tornar subpr LP em MILP 
			if (sistema_a->GetFlagVarBin() == 1)
				SubpLP2MILP(x_spr, mod, &lim_var_bin);

			restricao.clear();
		}
		// as var. primais de cada problema só são utilizadas para montar o corte do subproblema anterior
	}

	alfa = alfa_0;
	return ( 0 );
}
int Hrstc::ResolverBackward0(CMatrizEsparsa * x_spr)
{
	// Backward0 para construir cortes com a solução do ED Ext.

	// Quando consideram-se os cortes de Benders as janelas futuras perdem o sentido, portanto não são consideradas!
	//double alfa_0 = alfa;
	//alfa = 1.0;		// no backward não são considerados os termos proximais
	int subp_status = 0;
	int mod_ant;
	int n_supbB1_por_cen = subpB1etapa.size() / sistema_a->GetNCenarios();		// numero de subproblemas por cenário da 1ª etapa

	GRBLinExpr restricao;
	double coeficiente;
	GRBVar variavel;
	double rhs_cte, rhs_var, rhs_bin;
	//
	//vetorfloat RHS_acop;
	
	// Etapa 1,2 e 3 (problemas de segundo estágio, problema de acoplamento e problemas de primeiro estágio)
	for (int mod = n_modelos - 1; mod > 0; mod--)
	{
		// problema já resolvido
		resultadosGurobi->AlocarXmed(x_spr);		// aloca solução do problema estendido em x_spr

		if ( (mod >= subpB1etapa[0]) && ((mod - subpB2etapa[0]) % n_supbB1_por_cen == 0) )			// subproblema de final de horizonte
			fo_subp[mod] = CalcularFuncaoObjetivo(mod, x_spr);
		else if (mod == subpB2etapa[0])		// subproblemas de acoplamento e da etapa 3 dependem dos cenários
		{
			fo_subp[mod] = CalcularFuncaoObjetivo(mod, x_spr);
			for (size_t i = 1; i < subpB2etapa.size(); i++)
				fo_subp[mod] += fo_subp[subpB2etapa[i]];
		}
		else
			fo_subp[mod] = fo_subp[mod + 1] + CalcularFuncaoObjetivo(mod, x_spr);
		//if (mod == 0)
		//{
		//	fo_subp[mod] = fo_subp[mod + 1] + CalcularFuncaoObjetivo(mod, x_spr);
		//	cout << fo_subp[mod] << endl;
		//}

		// Identificar subproblema anterior

		if ((mod >= subpB1etapa[0]) && ((mod - (subpB2etapa[0]+1)) % n_supbB1_por_cen == 0))		// se este modelo for de inicio do 2 estágio, acoplamento com o subproblema que é da segunda etapa
			mod_ant = subpB2etapa[0];
		else	// acoplamento com o subproblema anterior dentro do segundo ou do primeiro estágio
			mod_ant = mod - 1;

		// armazenar coeficientes de cada variável para cada subproblema (para o subproblema anterior)
		CMatrizEsparsa L_it(1, n);
		for (size_t ii = 0; ii < acomplamentos[mod - 1].size(); ii++)		// loop no numero de acoplamentos (se algum for com o subp anterior adicioná-lo no corte)
		{
			//if ( (int(acomplamentos[mod - 1][ii][1]) >= nvar_acum[mod_ant]) && (int(acomplamentos[mod - 1][ii][1]) < nvar_acum[mod_ant + 1]) )
			//	L_it.SubstituirElemento(0, int(acomplamentos[mod - 1][ii][1]), L_it.GetElemento(0, int(acomplamentos[mod - 1][ii][1])) + resultadosGurobi->GetLambda(int(acomplamentos[mod - 1][ii][4])) * acomplamentos[mod - 1][ii][3]);
			//else
				L_it.SubstituirElemento(0, int(acomplamentos[mod - 1][ii][1]), L_it.GetElemento(0, int(acomplamentos[mod - 1][ii][1])) + resultadosGurobi->GetLambda(int(acomplamentos[mod - 1][ii][4])) * acomplamentos[mod - 1][ii][3]);
		}

		// atribuir valores para o vetor de coeficientes acumulados
		// (considerar todos cortes ativos, m_corte = 1)
		double mult_corte = 0;
		if ( (mod >= subpB1etapa[0]) && ((mod - subpB2etapa[0]) % n_supbB1_por_cen == 0) )		// subproblema de final de horizonte
		{}
		else if (mod == subpB2etapa[0])		// subproblemas de acoplamento e da etapa 3 dependem dos cenários
		{
			// multiplicador de todos os cortes do subp. mod e para todos os cenários
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				CMatrizEsparsa m_corte(1, itB + 1);
				for (int iter = 0; iter < itB + 1; iter++)		
				{
					m_corte.SubstituirElemento(0, iter, 1);
				}
				if (L[subpB2etapa[(sistema_a->GetNCenarios() - 1 - cen) + 1]].GetNnz() != 0)
				{
					m_corte.MultiplicarPorMatriz(&L[subpB2etapa[(sistema_a->GetNCenarios() - 1 - cen) + 1] - 1]);
					L_it.SomarComMatriz(&m_corte);
				}
			}
		}
		else
		{
			// multiplicador de todos os cortes do subp. mod
			CMatrizEsparsa m_corte(1, itB + 1);
			for (int iter = 0; iter < itB + 1; iter++)		
			{
				m_corte.SubstituirElemento(0, iter, 1);
			}
			if (L[mod].GetNnz() != 0)
			{
				m_corte.MultiplicarPorMatriz(&L[mod]);		// L[mod_apr - 1], que nesses casos = L[mod]
				L_it.SomarComMatriz(&m_corte);
			}
		}

		// Construir corte para o subproblema anterior (mod - 1)
		restricao.clear();
		rhs_cte = 0;
		rhs_var = 0;
		rhs_bin = 0;
		// incluir variável de aproximação do CF
		coeficiente = 1;
		if ((mod >= subpB1etapa[0]) && ((mod - (subpB2etapa[0]+1)) % n_supbB1_por_cen == 0))
		{
			int cen = (mod - subpB1etapa[0]) / n_supbB1_por_cen;
			variavel = vars[mod_ant][numero_var_subp[mod_ant] + numero_var_ad[mod_ant] + numero_var_ad2[mod_ant] + n_var_folga[mod_ant] + (sistema_a->GetNCenarios() - 1 - cen)];		// ultima variável
		}
		else
			variavel = vars[mod_ant][numero_var_subp[mod_ant] + numero_var_ad[mod_ant] + numero_var_ad2[mod_ant] + n_var_folga[mod_ant]];		// ultima variável
		restricao.addTerms( &coeficiente, &variavel, 1);
		rhs_cte += fo_subp[mod];

		// Percorrer a matriz PI para adicionar cortes (ou L da iteração corrente)!!!
		CMatrizEsparsa L_ac(1, n);
		int l;

		// Calcular o RHS primeiro para avaliar a ordem dos coeficientes
		l = L_it.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
		while ( l != -1 )	// percorre todos elementos nnz da linha
		{
			if ( (L_it.GetValorCol(l) >= nvar_acum[mod_ant]) && (L_it.GetValorCol(l) < nvar_acum[mod_ant + 1]) )		// se existir acoplamento com o subproblema anterior
				rhs_cte += L_it.GetValorVal(l) * x_spr->GetElemento(L_it.GetValorCol(l), 0);
			l = L_it.GetValorLprox(l);		
		}
		l = L_it.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
		while ( l != -1 )	// percorre todos elementos nnz da linha
		{
			if ( (L_it.GetValorCol(l) >= nvar_acum[mod_ant]) && (L_it.GetValorCol(l) < nvar_acum[mod_ant + 1]) )		// se existir acoplamento com o subproblema anterior
			{
				// Conferir ordem dos coeficientes, se for menor que coeftol não adiciona-se o termo
				if (abs(L_it.GetValorVal(l) * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB)) >= coeftol * rhs_cte)
				///if ( (abs(L_it.GetValorVal(l) * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB)) >= coeftol * rhs_cte) && (abs(L_it.GetValorVal(l) * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB)) <= (rhs_cte/coeftol)) )
				{
					coeficiente = L_it.GetValorVal(l);
					//if (coeficiente * vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]].get(GRB_DoubleAttr_UB) < coeftol)
					//	coeficiente = 0;
					variavel = vars[mod_ant][L_it.GetValorCol(l) - nvar_acum[mod_ant]];		// endereço deve ser referente à variável do modelo
					restricao.addTerms( &coeficiente, &variavel, 1);
					//rhs_cte += coeficiente * x_spr->GetElemento(L_it.GetValorCol(l), 0);
				}
			}
			else
			{
				rhs_var += L_it.GetValorVal(l) * x_spr->GetElemento(L_it.GetValorCol(l), 0);
				L_ac.InserirElemento(0, L_it.GetValorCol(l), L_it.GetValorVal(l));
			}
			l = L_it.GetValorLprox(l);		
		}


		// vetor de coeficientes para cada subproblema é o valor do PI para aquela iteração!
		L[mod - 1].JuntarColuna(&L_ac);			// matriz do subproblema corrente, com os coeficientes da aproximação do subproblema anterior a este
		// Armazenar todos os coeficientes menos do subp. anterior
		rhs_cte += rhs_var;

		rhs_orig[mod - 1].push_back(rhs_cte);		// para atualizar o corte necessita-se do rhs_original
		modelosGRB[mod_ant]->addConstr(restricao, GRB_GREATER_EQUAL, rhs_cte - rhs_var, "");			// não precisa subtrair rhs_var, pois tem-se o lhs_constante que seria igual ao rhs_var.
		modelosGRB[mod_ant]->update();
		delete constr[mod_ant];
		constr[mod_ant] = modelosGRB[mod_ant]->getConstrs();
			
	}

	// Etapa 2 (problema de acoplamento)
	// Etapa 3 (problemas de primeiro estágio)
	// nao tem como deixar essas 3 etapas em um loop só?!? com ifs,...
	return ( 0 );
}

void Hrstc::IdentAcoplamentos(int modelo)
{
	// Cria vetor com a posição da linha (referente à restrição do subproblema) e da variável (referente ao vetor x_spr) que acopla o subprolema com os subproblemas < t.
	// assim como o numero do subproblema referente à variável adicionada
	// e o coeficiente que multiplica a variável!
	vetorfloat2 acplmnt;
	vetorfloat lc;
	lc.resize(5);

	// coluna limite (soma das variáveis dos subproblemas anteriores)
	int n_col  = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen = 0;
	int t_init = 0;
	int index_subp = 0;
	double coef = 0;

	CMatrizEsparsa M;
	int nt_rest;
	for (size_t i = 0; i < posicao_rest_din[modelo].size(); i++)	// loop para as 6 ou 7 submatrizes q acoplam no tempo, Balanço hídrico, min. up/down time, up/down ramp rate e custo de partida
	{
		M = A[i];
		nt_rest = M.GetNlin() / T;		// = R (balanço hidrico), Tup, Tdown, I (rampa up), I (rampa down) e I (custo de partida)
										// = R (balanço hidrico), sum(i) para Tup > 0, sum(i) para Tdown > 0, I Tupdown, I (rampa up), I (rampa down), I (limite de ptmin)
		// varrer linhas e colunas
		for (int t = vtr_nos_pri[modelo][0]; t < (vtr_nos_pri[modelo].size() + vtr_nos_fut[modelo].size() + vtr_nos_pri[modelo][0]); t++)
		{
			for (int linha = t*nt_rest; linha < (t + 1)*nt_rest; linha++)
			{
				if ( t < sistema_a->GetTt2() )
				{
					for (int coluna = 0; coluna < n_col*t; coluna++)
					{
						coef = M.GetElemento(linha, coluna);
						if ( coef != 0 )
						{
							// identificar a qual problema pertence a var. referente à coluna
							index_subp = modelo - 1;		// referente ao subproblema corrente
							while (coluna < nvar_acum[index_subp])
								index_subp--;
							lc[0] = posicao_rest_din[modelo][i] + linha - t*nt_rest;
							lc[1] = coluna;
							lc[2] = index_subp;		// subproblema referente à variavel coluna
							// aqui para atribuir o coeficiente que vai para o corte depende da convenção, para restrições de >= ou == o coef. tem o mesmo sinal da restrição
							//if (constr[modelo][lc[0]].get(GRB_CharAttr_Sense) == GRB_LESS_EQUAL)	nao pois o multiplicador tambem tem sinal trocado!!!
							//	lc[3] = - coef;
							//else
								lc[3] = coef;
							lc[4] =	posicao_rest_din_Ext[i] + linha;		// posição da restrição da matriz original de restrições (modelo Ext)
							acplmnt.push_back(lc);
						}
					}
				}
				else
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					t_init = sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1());
					for (int coluna = 0; coluna < n_col*sistema_a->GetTt1(); coluna++)
					{
						coef = M.GetElemento(linha, coluna);
						if ( coef != 0 )
						{
							// identificar a qual problema pertence a var. referente à coluna
							index_subp = modelo - 1;		// referente ao subproblema corrente
							while (coluna < nvar_acum[index_subp])
								index_subp--;
							lc[0] = posicao_rest_din[modelo][i] + linha - t*nt_rest;
							lc[1] = coluna;
							lc[2] = index_subp;		// subproblema referente à variavel coluna
							//if (constr[modelo][lc[0]].get(GRB_CharAttr_Sense) == GRB_LESS_EQUAL)
							//	lc[3] = - coef;
							//else
								lc[3] = coef;
							lc[4] =	posicao_rest_din_Ext[i] + linha;		// posição da restrição da matriz original de restrições (modelo Ext)
							acplmnt.push_back(lc);
						}
					}
					for (int coluna = n_col*t_init + cen*flag2*sistema_a->hidreletricasVtr.size(); coluna < n_col*t + cen*flag2*sistema_a->hidreletricasVtr.size(); coluna++)
					{
						coef = M.GetElemento(linha, coluna);
						if ( coef != 0 )
						{
							// identificar a qual problema pertence a var. referente à coluna
							index_subp = modelo - 1;		// referente ao subproblema corrente
							while (coluna < nvar_acum[index_subp])
								index_subp--;
							lc[0] = posicao_rest_din[modelo][i] + linha - t*nt_rest;
							lc[1] = coluna;
							lc[2] = index_subp;		// subproblema referente à variavel coluna
							//if (constr[modelo][lc[0]].get(GRB_CharAttr_Sense) == GRB_LESS_EQUAL)
							//	lc[3] = - coef;
							//else
								lc[3] = coef;
							lc[4] =	posicao_rest_din_Ext[i] + linha;		// posição da restrição da matriz original de restrições (modelo Ext)
							acplmnt.push_back(lc);
						}
					}
				}
			}
		}
	}
	// Vetores L e acoplamentos tem uma posição a menos do que o numero de subproblemas, pois não existem para o subproblema 0
	acomplamentos[modelo - 1] = acplmnt;
	//L[modelo - 1] = CMatrizEsparsa(0, acplmnt.size());		antigo
}
void Hrstc::SubpMILP2LP(CMatrizEsparsa * x_spr, int modelo, vetorint * limites)
{
	// Usar fixedmodel do gurobi não é o suficiente aqui, pois precisamos adicionar uma restrição para fixar as variáveis com um epsilon,
	// sem isso as restrições com somente var. binárias não terão valores de multiplicadores...
	// os multiplicadores das var. binárias (fixadas) do subproblema atual são dadas pelo reduced cost!

	//double eps = 1 / (20*sistema_a->GetCustoDeficit()/0.0036);
	//double eps = 1 / sistema_a->GetCustoDeficit();
	double eps = 0;

	// Somente fixando as var. bin. (sem considerá-las nos cortes) o lower bound fica maior que a solução ótima do problema
	// Considerando nos cortes, o problema é qual solução usar quando se esta resolvendo os subproblemas no forward

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen, cen_subp, delta_x, delta_subp, jj;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	cen_subp = 0;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
	{
		cen = 0;
		if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
			cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		delta_x = nt*vtr_nos_pri[modelo][vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();
		delta_subp = nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size();

		// fixar somente as variáveis binárias
		if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
		{
			for (int i = 0; i < I; i++)
			{
				//gravar limites	(gravar limites pois as variáveis já podem ter sido fixas pelas condições iniciais!!)
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 2*I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 2*I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 3*I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 3*I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				//u
				vars[modelo][i + I + delta_subp].set(GRB_CharAttr_VType, 'C');
				//vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + I + delta_x, 0) - eps);
				//vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + I + delta_x, 0) + eps);
				//up
				vars[modelo][i + 2*I + delta_subp].set(GRB_CharAttr_VType, 'C');
				//vars[modelo][i + 2*I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + 2*I + delta_x, 0) - eps);
				//vars[modelo][i + 2*I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + 2*I + delta_x, 0) + eps);
				//ud
				vars[modelo][i + 3*I + delta_subp].set(GRB_CharAttr_VType, 'C');
				//vars[modelo][i + 3*I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + 3*I + delta_x, 0) - eps);
				//vars[modelo][i + 3*I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + 3*I + delta_x, 0) + eps);
			}
		}
		else
		{
			for (int i = 0; i < I; i++)
			{
				//gravar limites	(gravar limites pois as variáveis já podem ter sido fixas pelas condições iniciais!!)
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				//u
				vars[modelo][i + I + delta_subp].set(GRB_CharAttr_VType, 'C');
				//vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + I + delta_x, 0) - eps);
				//vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + I + delta_x, 0) + eps);
			}
		}
		jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				//gravar limites
				limites->push_back(int(floor(vars[modelo][j + jj + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][j + jj + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				//z
				vars[modelo][j + jj + delta_subp].set(GRB_CharAttr_VType, 'C');
				//vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(j + jj + delta_x, 0) - eps);
				//vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(j + jj + delta_x, 0) + eps);
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		if ( (vtr_nos_pri[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
			cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
	}
	modelosGRB[modelo]->update();
}
void Hrstc::SubpMILP2LPFixing(CMatrizEsparsa * x_spr, int modelo, vetorint * limites)
{
	// Usar fixedmodel do gurobi não é o suficiente aqui, pois precisamos adicionar uma restrição para fixar as variáveis com um epsilon,
	// sem isso as restrições com somente var. binárias não terão valores de multiplicadores...
	// os multiplicadores das var. binárias (fixadas) do subproblema atual são dadas pelo reduced cost!

	double eps = 1 / PenVF;
	//double eps = 1 / (20*sistema_a->GetCustoDeficit()/0.0036);
	//double eps = 1 / sistema_a->GetCustoDeficit();
	//double eps = 0;
	//double eps = 0.01;

	// Somente fixando as var. bin. (sem considerá-las nos cortes) o lower bound fica maior que a solução ótima do problema
	// Considerando nos cortes, o problema é qual solução usar quando se esta resolvendo os subproblemas no forward

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen, cen_subp, delta_x, delta_subp, jj;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	cen_subp = 0;
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
	{
		cen = 0;
		if ( vtr_nos_pri[modelo][vtr_i] >= sistema_a->GetTt2() )
			cen = (vtr_nos_pri[modelo][vtr_i] - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		delta_x = nt*vtr_nos_pri[modelo][vtr_i] + flag2*cen*sistema_a->hidreletricasVtr.size();
		delta_subp = nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size();

		// fixar somente as variáveis binárias
		if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
		{
			for (int i = 0; i < I; i++)
			{
				//gravar limites	(gravar limites pois as variáveis já podem ter sido fixas pelas condições iniciais!!)
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 2*I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 2*I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 3*I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + 3*I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				//u
				vars[modelo][i + I + delta_subp].set(GRB_CharAttr_VType, 'C');
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + I + delta_x, 0) - eps);
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + I + delta_x, 0) + eps);
				//up
				vars[modelo][i + 2*I + delta_subp].set(GRB_CharAttr_VType, 'C');
				vars[modelo][i + 2*I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + 2*I + delta_x, 0) - eps);
				vars[modelo][i + 2*I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + 2*I + delta_x, 0) + eps);
				//ud
				vars[modelo][i + 3*I + delta_subp].set(GRB_CharAttr_VType, 'C');
				vars[modelo][i + 3*I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + 3*I + delta_x, 0) - eps);
				vars[modelo][i + 3*I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + 3*I + delta_x, 0) + eps);
			}
		}
		else
		{
			for (int i = 0; i < I; i++)
			{
				//gravar limites	(gravar limites pois as variáveis já podem ter sido fixas pelas condições iniciais!!)
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][i + I + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				//u
				vars[modelo][i + I + delta_subp].set(GRB_CharAttr_VType, 'C');
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(i + I + delta_x, 0) - eps);
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(i + I + delta_x, 0) + eps);
			}
		}
		jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				//gravar limites
				limites->push_back(int(floor(vars[modelo][j + jj + delta_subp].get(GRB_DoubleAttr_LB) + 0.5)));
				limites->push_back(int(floor(vars[modelo][j + jj + delta_subp].get(GRB_DoubleAttr_UB) + 0.5)));
				//z
				vars[modelo][j + jj + delta_subp].set(GRB_CharAttr_VType, 'C');
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_LB, x_spr->GetElemento(j + jj + delta_x, 0) - eps);
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_UB, x_spr->GetElemento(j + jj + delta_x, 0) + eps);
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		if ( (vtr_nos_pri[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
			cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
	}
	modelosGRB[modelo]->update();
}
void Hrstc::SubpLP2MILP(CMatrizEsparsa * x_spr, int modelo, vetorint * limites)
{
	// Deixar var. binárias em seus limites originas novamente

	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen_subp, delta_x, delta_subp, jj;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	cen_subp = 0;
	char tipo_var;
	if (sistema_a->GetFlagVarBin() == 0)
		tipo_var = 'C';
	else
		tipo_var = 'B';

	// Reverter variáveis binárias e seus limites
	cen_subp = 0;
	delta_x = 0;		// delta_x usado para percorrer vetor limites, que já está na ordem apropriada
	for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
	{
		delta_subp = nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size();
		if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
		{
			for (int i = 0; i < I; i++)
			{
				//u
				vars[modelo][i + I + delta_subp].set(GRB_CharAttr_VType, tipo_var);
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_LB, limites->at(delta_x++));
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_UB, limites->at(delta_x++));
				//up
				vars[modelo][i + 2*I + delta_subp].set(GRB_CharAttr_VType, tipo_var);
				vars[modelo][i + 2*I + delta_subp].set(GRB_DoubleAttr_LB, limites->at(delta_x++));
				vars[modelo][i + 2*I + delta_subp].set(GRB_DoubleAttr_UB, limites->at(delta_x++));
				//ud
				vars[modelo][i + 3*I + delta_subp].set(GRB_CharAttr_VType, tipo_var);
				vars[modelo][i + 3*I + delta_subp].set(GRB_DoubleAttr_LB, limites->at(delta_x++));
				vars[modelo][i + 3*I + delta_subp].set(GRB_DoubleAttr_UB, limites->at(delta_x++));
			}
		}
		else
		{
			for (int i = 0; i < I; i++)
			{
				//u
				vars[modelo][i + I + delta_subp].set(GRB_CharAttr_VType, tipo_var);
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_LB, limites->at(delta_x++));
				vars[modelo][i + I + delta_subp].set(GRB_DoubleAttr_UB, limites->at(delta_x++));
			}
		}
		jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				//z
				vars[modelo][j + jj + delta_subp].set(GRB_CharAttr_VType, tipo_var);
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_LB, limites->at(delta_x++));
				vars[modelo][j + jj + delta_subp].set(GRB_DoubleAttr_UB, limites->at(delta_x++));
			}
			jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
		if ( (vtr_nos_pri[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
			cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
	}
	modelosGRB[modelo]->update();
}
void Hrstc::DeterminarIndiceVarBin()
{
	index_var_bin.resize(n_modelos);
	double eps = 0;
	int nt = ( n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size() ) / T;
	int cen_subp, delta_subp, jj;
	int R = sistema_a->hidreletricasVtr.size();
	int I = sistema_a->termeletricasVtr.size();
	int B = sistema_a->barrasVtr.size();
	int JJ = 0;
	for (int r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	for (int modelo = 0; modelo < n_modelos; modelo++)
	{
		delta_subp = 0;
		cen_subp = 0;
		for (size_t vtr_i = 0; vtr_i < vtr_nos_pri[modelo].size(); vtr_i++)
		{
			delta_subp = nt*vtr_i + flag2*cen_subp*sistema_a->hidreletricasVtr.size();

			// fixar somente as variáveis binárias
			if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
			{
				for (int i = 0; i < I; i++)
				{
					// indice das var. binarias
					index_var_bin[modelo].push_back(i + I + delta_subp);
					index_var_bin[modelo].push_back(i + 2*I + delta_subp);
					index_var_bin[modelo].push_back(i + 3*I + delta_subp);
				}
			}
			else
			{
				for (int i = 0; i < I; i++)
				{
					// indice das var. binarias
					index_var_bin[modelo].push_back(i + I + delta_subp);
				}
			}
			jj = I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 2*JJ;
			for (int r = 0; r < R; r++)
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					// indice das var. binarias
					index_var_bin[modelo].push_back(j + jj + delta_subp);
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if ( (vtr_nos_pri[modelo][vtr_i] + 1 > sistema_a->GetTt1()) && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1()) == 0 ) )	// nó de final do horizonte
				cen_subp ++;		// só adiciona valor quando tem um nó que é de final de horizonte no subproblema!
		}
	}
}
void Hrstc::AtualizarCortes(CMatrizEsparsa & x_spr, int modelo)
{

	// n seria mais fácil multiplicar L_it por x_spr? para calcular o rhs + fo_subp[mod]?!?
	// o rhs de todos os cortes deve ser atualizado pelo lhs variável, que são todas variáveis do corte que n sõa decisões no subproblema atual
	// calcular lhs_var
	// como identificar quais variáveis estão no corte que devem ser atualizadas? pelo L, tds var. q n forem do subproblema em questão
	// a matriz L de cada iteração já possui os coeficientes atualizados (lambda + mi*lambda_-1,...) para cada subproblema!

	// atualizar rhs_var para a solução atual e para o subproblema modelo
	double lhs_var = 0;
	int modelo_apr;		// modelo aproximado pelos cortes (subp. do futuro!)
	CMatrizEsparsa PI_temp(1, n);
	int l;

	if (modelo != subpB2etapa[0])		// acoplamento com o subproblema futuro dentro do segundo ou do primeiro estágio, exceto subp de final de horizonte
	{
		modelo_apr = modelo + 1;
		// atualizar cada corte (restrição) adicionado = itB
		for (int iter = 0; iter < itB; iter++)		// n atualizar o corte da ultima iteração, jah está correto! (após o backward)
		{
			lhs_var = 0;
			PI_temp.RemoverTodosElementos();
			PI_temp.InserirMatriz(0, 0, 0, n - 1, &L[modelo_apr - 1], iter, 0);
			// multiplicar PI_temp pelos elementos que n são decisão do subproblema corrente
			l = PI_temp.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
			while ( l != -1 )	// percorre todos elementos nnz da linha
			{
				lhs_var += PI_temp.GetValorVal(l) * x_spr.GetElemento(PI_temp.GetValorCol(l), 0);
				l = PI_temp.GetValorLprox(l);		
			}
			constr[modelo][iter + nrest_orig[modelo]].set(GRB_DoubleAttr_RHS, rhs_orig[modelo_apr - 1][iter] - lhs_var);
		}

	}
	else	// subproblema de acoplamento (tem um corte para cada cenário)
	{
		for ( int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			modelo_apr = subpB2etapa[(sistema_a->GetNCenarios() - 1 - cen) + 1];
			// atualizar cada corte (restrição) adicionado = itB * n_cenarios para o subp. de acoplamento
			for (int iter = 0; iter < itB; iter++)
			{
				lhs_var = 0;
				PI_temp.RemoverTodosElementos();
				PI_temp.InserirMatriz(0, 0, 0, n - 1, &L[modelo_apr - 1], iter, 0);
				// multiplicar PI_temp pelos elementos que n são decisão do subproblema corrente
				l = PI_temp.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
				while ( l != -1 )	// percorre todos elementos nnz da linha
				{
					lhs_var += PI_temp.GetValorVal(l) * x_spr.GetElemento(PI_temp.GetValorCol(l), 0);
					l = PI_temp.GetValorLprox(l);		
				}
				constr[modelo][cen + iter*sistema_a->GetNCenarios() + nrest_orig[modelo]].set(GRB_DoubleAttr_RHS, rhs_orig[modelo_apr - 1][iter] - lhs_var);
			}
		}
	}
	modelosGRB[modelo]->update();

}
void Hrstc::AtualizarCortesBin(CMatrizEsparsa * x_spr, int modelo)
{
	//// atualizar rhs_var para a solução atual e para o subproblema modelo
	//double lhs_bin = 0;
	//int modelo_apr;		// modelo aproximado pelos cortes (subp. do futuro!)
	//CMatrizEsparsa PI_temp(1, n);
	//int l;

	//if (modelo != subpB2etapa[0])		// acoplamento com o subproblema futuro dentro do segundo ou do primeiro estágio, exceto subp de final de horizonte
	//{
	//	modelo_apr = modelo + 1;
	//	// atualizar cada corte (restrição) adicionado = itB
	//	for (int iter = 0; iter < itB; iter++)		// n atualizar o corte da ultima iteração, jah está correto! (após o backward)
	//	{
	//		lhs_bin = 0;
	//		PI_temp.RemoverTodosElementos();
	//		PI_temp.InserirMatriz(0, 0, 0, n - 1, &L_bin[modelo_apr - 1], iter, 0);
	//		// multiplicar PI_temp pelos elementos que n são decisão do subproblema corrente
	//		l = PI_temp.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
	//		while ( l != -1 )	// percorre todos elementos nnz da linha
	//		{
	//			lhs_bin += PI_temp.GetValorVal(l) * x_spr->GetElemento(PI_temp.GetValorCol(l), 0);
	//			l = PI_temp.GetValorLprox(l);		
	//		}
	//		constr[modelo][iter + nrest_orig[modelo]].set(GRB_DoubleAttr_RHS, rhs_orig_bin[modelo_apr - 1][iter] - lhs_bin);
	//	}

	//}
	//else	// subproblema de acoplamento (tem um corte para cada cenário)
	//{
	//	for ( int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	//	{
	//		modelo_apr = subpB2etapa[(sistema_a->GetNCenarios() - 1 - cen) + 1];
	//		// atualizar cada corte (restrição) adicionado = itB * n_cenarios para o subp. de acoplamento
	//		for (int iter = 0; iter < itB; iter++)
	//		{
	//			lhs_bin = 0;
	//			PI_temp.RemoverTodosElementos();
	//			PI_temp.InserirMatriz(0, 0, 0, n - 1, &L_bin[modelo_apr - 1], iter, 0);
	//			// multiplicar PI_temp pelos elementos que n são decisão do subproblema corrente
	//			l = PI_temp.GetValorLprim(0);		// eu sei q a matriz PI tem somente uma linha
	//			while ( l != -1 )	// percorre todos elementos nnz da linha
	//			{
	//				lhs_bin += PI_temp.GetValorVal(l) * x_spr->GetElemento(PI_temp.GetValorCol(l), 0);
	//				l = PI_temp.GetValorLprox(l);		
	//			}
	//			constr[modelo][cen + iter*sistema_a->GetNCenarios() + nrest_orig[modelo]].set(GRB_DoubleAttr_RHS, rhs_orig_bin[modelo_apr - 1][iter] - lhs_bin);
	//		}
	//	}
	//}
	//modelosGRB[modelo]->update();

}
double Hrstc::CalcularFuncaoObjetivo(int modelo, CMatrizEsparsa * x_spr)
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
	int flag1a = 0;
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	// Variáveis das janelas principais
	//x = [pt u cp (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
	//x = [pt u up ud (F) (teta) ph v d s phmax phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
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
					funcao_obj += x_spr->GetElemento(i + (3+flag7)*I + delta_x, 0) * 1*deltaT;		//F
				else
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
					{
						funcao_obj += x_spr->GetElemento(i + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
						funcao_obj += x_spr->GetElemento(i + I + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;		//u
					}
					else
					{
						funcao_obj += x_spr->GetElemento(i + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
						funcao_obj += x_spr->GetElemento(i + I + delta_x, 0) * (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT;		//u
					}
				}
				if (sistema_a->GetFlagTbinaryModel() == 0)
					funcao_obj += x_spr->GetElemento(i + 2*I + delta_x, 0) * deltaT;		//cp
				else
					funcao_obj += x_spr->GetElemento(i + 2*I + delta_x, 0) * deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		//up
			}
			if (sistema_a->GetFlagModeloRede() > 0)
				for (int b = 0; b < B; b++)
					funcao_obj += x_spr->GetElemento(b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT;		//def
			else
				funcao_obj += x_spr->GetElemento(I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT;		//def
		}
		else	// nós do estágio 2
		{
			deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
			for (int i = 0; i < I; i++)
			{

				if (sistema_a->GetFlagInitAproxCT() > 1)
					funcao_obj += x_spr->GetElemento(i + (3+flag7)*I + delta_x, 0) * 1*deltaT;		//F
				else
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
					{
						funcao_obj += x_spr->GetElemento(i + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
						funcao_obj += x_spr->GetElemento(i + I + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;		//u
					}
					else
					{
						funcao_obj += x_spr->GetElemento(i + delta_x, 0) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;			//pt
						funcao_obj += x_spr->GetElemento(i + I + delta_x, 0) * (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*sistema_a->termeletricasVtr[i].GetPmin())*deltaT;		//u
					}
				}
				if (sistema_a->GetFlagTbinaryModel() == 0)
					funcao_obj += x_spr->GetElemento(i + 2*I + delta_x, 0) * deltaT;		//cp
				else
					funcao_obj += x_spr->GetElemento(i + 2*I + delta_x, 0) * deltaT * sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		//up
			}
			if (sistema_a->GetFlagModeloRede() > 0)
				for (int b = 0; b < B; b++)
					funcao_obj += x_spr->GetElemento(b + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT;		//def
			else
				funcao_obj += x_spr->GetElemento(I*(3 + flag4 + flag7) + (4+flag3)*R + 3*JJ + delta_x, 0) * sistema_a->GetCustoDeficit()*deltaT;		//def
			if ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0)		// && ((vtr_nos_pri[modelo][vtr_i] + 1 - sistema_a->GetTt1()) > 0))
				if (sistema_a->GetFlagVfol() == true)
					for (int r = 0; r < R; r++)
						funcao_obj += x_spr->GetElemento(r + I*(3 + flag4 + flag7) + flag1a*(B - 1) + (4+flag3)*R + 3*JJ + flag1d*B + (1 - flag1d) + delta_x, 0) * Cvfol*deltaT;		//vfol
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
				{
					if (sistema_a->GetFlagTbinaryModel() == 0)
						funcao_obj += pow(x_spr->GetElemento(i + delta_x, 0), 2) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;			//pt^2
					else
						cout << "CriarFuncaoObjetivo2: Termo quadrático não implementado para essa modelagem!!" << endl;
				}
			}
			else	// nós do estágio 2
			{
				deltaT = sistema_a->GetDeltaT2()* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
				for (int i = 0; i < I; i++)
					funcao_obj += pow(x_spr->GetElemento(i + delta_x, 0), 2) * sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;			//pt^2
			}
		}
	}
	return funcao_obj;
}
void Hrstc::RemoverVarFol(int modelo)
{
	// remover variaveis de folga quando existirem
	GRBVar var_fol;
	for (size_t cc = 0; cc < ind_var_fol[modelo].size(); cc++)
	{
		var_fol = modelosGRB[modelo]->getVar(ind_var_fol[modelo][cc]);
		modelosGRB[modelo]->remove(var_fol);
	}
	modelosGRB[modelo]->update();
	delete vars[modelo];
	vars[modelo] = modelosGRB[modelo]->getVars();
	ind_var_fol[modelo].clear();
}
void Hrstc::AdicionarVarFol(int modelo)
{
	// adicionar variaveis de folga
	// duas variávais por restrição: varfol+ varfol-
	modelosGRB[modelo]->getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9);
	modelosGRB[modelo]->computeIIS();	
	//if ((itB == 3) && (modelo == 46))
	//{
	//	cout << modelosGRB[modelo]->get(GRB_IntAttr_NumVars) << endl;
	//	modelosGRB[modelo]->write("Hrst_subprob.ilp");
	//}
	double * coef = new double;
	*coef = 1;
	double * coef2 = new double;
	*coef2 = -1;
	int cont = 0;
	// somente para as restrições de acoplamento temporal
	int nt_rest;
	for ( size_t c = 0; c < posicao_rest_din[modelo].size(); c++)
	{
		nt_rest = A[c].GetNlin() / T;
		for (int cc = posicao_rest_din[modelo][c]; cc < posicao_rest_din[modelo][c]+nt_rest; cc++)
		{
			// identificar quais restrições
			if (constr[modelo][cc].get(GRB_IntAttr_IISConstr))	
			{
				modelosGRB[modelo]->addVar(0, GRB_INFINITY, PenVF, GRB_CONTINUOUS, 1, &constr[modelo][cc], coef);
				ind_var_fol[modelo].push_back(modelosGRB[modelo]->get(GRB_IntAttr_NumVars) + cont);
				cont++;
				modelosGRB[modelo]->addVar(0, GRB_INFINITY, PenVF, GRB_CONTINUOUS, 1, &constr[modelo][cc], coef2);
				ind_var_fol[modelo].push_back(modelosGRB[modelo]->get(GRB_IntAttr_NumVars) + cont);
				cont++;
			}
		}
	}
	// somente para a restrição de balanço hídrico
	//for (int cc = posicao_rest_din[modelo][0]; cc < posicao_rest_din[modelo][1]; cc++)
	//{
	//	// identificar quais restrições
	//	if (constr[modelo][cc].get(GRB_IntAttr_IISConstr))	
	//	{
	//		modelosGRB[modelo]->addVar(0, GRB_INFINITY, PenVF, GRB_CONTINUOUS, 1, &constr[modelo][cc], coef);
	//		ind_var_fol[modelo].push_back(modelosGRB[modelo]->get(GRB_IntAttr_NumVars) + cont);
	//		cont++;
	//		modelosGRB[modelo]->addVar(0, GRB_INFINITY, PenVF, GRB_CONTINUOUS, 1, &constr[modelo][cc], coef2);
	//		ind_var_fol[modelo].push_back(modelosGRB[modelo]->get(GRB_IntAttr_NumVars) + cont);
	//		cont++;
	//	}
	//}
	//
	if ( WIN_SIZE > 1)
		cout << "Função AdicionarVarFol() deve ser modificada para janelas maiores!" << endl;		// nt_rest para cada periodo...
	//for (int cc = 0; cc < nrest_orig[modelo]; cc++)
	//{
	//	// identificar quais restrições
	//	if (constr[modelo][cc].get(GRB_IntAttr_IISConstr))
	//	{
	//		modelosGRB[modelo]->addVar(0, GRB_INFINITY, PenVF, GRB_CONTINUOUS, 1, &constr[modelo][cc], coef);
	//		ind_var_fol[modelo].push_back(modelosGRB[modelo]->get(GRB_IntAttr_NumVars) + cont);
	//		cont++;
	//		modelosGRB[modelo]->addVar(0, GRB_INFINITY, PenVF, GRB_CONTINUOUS, 1, &constr[modelo][cc], coef2);
	//		ind_var_fol[modelo].push_back(modelosGRB[modelo]->get(GRB_IntAttr_NumVars) + cont);
	//		cont++;
	//	}
	//}
	modelosGRB[modelo]->update();
	delete coef;
	delete coef2;
}
void Hrstc::SelecionarVarProx()
{
	//x = [pt u cp F teta ph v d s (phmax) phg q z def]
	//x = [pt u up ud F teta ph v d s (phmax) phg q z def]	"modern" model
	// ident_var = [ 0/1 ... ]
		
	// pt
	ident_var.push_back(1);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("pt");
	// u
	ident_var.push_back(1);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("u");
	if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
	{
		// up
		ident_var.push_back(0);
		if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("up");
		// ud
		ident_var.push_back(0);
		if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("ud");
	}
	else
	{
		// cp
		ident_var.push_back(0);
		if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("cp");
	}
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		// F
		ident_var.push_back(0);
		if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("F");
	}
	if ( sistema_a->GetFlagModeloRede() == 1)
	{
		// teta
		ident_var.push_back(0);
		if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("teta");
	}
	// ph
	ident_var.push_back(0);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("ph");
	// v
	ident_var.push_back(1);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("v");
	// d
	ident_var.push_back(0);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("d");
	// s
	ident_var.push_back(0);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("s");
	if (sistema_a->GetFlagPhmax() == 1)
	{
		// phmax
		ident_var.push_back(0);
		if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("phmax");
	}
	// phg
	ident_var.push_back(1);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("phg");
	// q
	ident_var.push_back(0);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("q");
	// z
	ident_var.push_back(0);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("z");
	// def
	ident_var.push_back(1);
	if (ident_var[ident_var.size() - 1]) stringVariaveis.push_back("def");
}
void Hrstc::SetLog(ofstream * log)
{
	log_auxiliar = log;
	//Escrever Ajustes da Heuristica
	if ( log_auxiliar->is_open() )
	{
		// Ajustes no .h
		*log_auxiliar << "Heuristic settings: " << char(9) << APPROACH << char(9) << WIN_SIZE << char(9) << WIN_FUT_SIZE;
		*log_auxiliar << char(9) << ALFA << char(9) << BETA << char(9) << GAMA << char(9) << INITCUT << char(9) << BIN_MW;
		*log_auxiliar << char(9) << BIN_FW << char(9) << NORM_P << char(9) << TIME_R << endl;

		// Varíáveis do termo proximal
		*log_auxiliar << "Proximal terms: [ ";
		for (size_t i = 0; i < stringVariaveis.size(); i++)
			*log_auxiliar << stringVariaveis[i] << " ";
		*log_auxiliar << "]" << endl;
		*log_auxiliar << endl;
	}
	else
		cout << "Unable to open file";
}
void Hrstc::ImprimirSolSub(int modelo)
{
	ostringstream n_sub;
	n_sub << modelo;
	modelosGRB[modelo]->write("Hrst_subp" + n_sub.str() + ".lp");
	// Escrever solulção do subproblema
	ofstream * inFile;
	inFile = new ofstream( "subp" + n_sub.str() + "sol.txt", ios::out );
	if ( inFile->is_open() )                                                            
	{
		*inFile << std::scientific << setprecision(10);
		*inFile << "obj = " << char(9) << modelosGRB[modelo]->get(GRB_DoubleAttr_ObjVal) << endl;
		*inFile << endl;
		for (size_t i = 0; i < modelosGRB[modelo]->get(GRB_IntAttr_NumVars); i++)
			*inFile << "C" << i << char(9) << vars[modelo][i].get(GRB_DoubleAttr_X) << endl;
		*inFile << endl;
		for (size_t i = 0; i < modelosGRB[modelo]->get(GRB_IntAttr_NumConstrs); i++)
			*inFile << "R" << i << char(9) << constr[modelo][i].get(GRB_DoubleAttr_Pi) << endl;
	}
	else
		cout << "Unable to open file";
	inFile->close();

}
void Hrstc::ImprimirSolSub(int modelo, int nvar)
{
	ostringstream n_sub;
	n_sub << modelo;
	modelosGRB[modelo]->write("Hrst_subp" + n_sub.str() + ".lp");
	// Escrever solulção do subproblema
	ofstream * inFile;
	inFile = new ofstream( "subp" + n_sub.str() + "sol.txt", ios::out );
	if ( inFile->is_open() )                                                            
	{
		*inFile << std::scientific << setprecision(10);
		*inFile << "obj = " << char(9) << modelosGRB[modelo]->get(GRB_DoubleAttr_ObjVal) << endl;
		*inFile << endl;
		for (size_t i = 0; i < nvar; i++)
			*inFile << "C" << i << char(9) << vars[modelo][i].get(GRB_DoubleAttr_X) << endl;
	}
	else
		cout << "Unable to open file";
	inFile->close();

}
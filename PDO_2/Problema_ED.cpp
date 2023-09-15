#include "Problema_ED.h"

// ARRUMAR FUNÇÕES PARA DELTAt1 != DELTAt2 (todas que usam GetDeltaT)
CMatrizEsparsa CProblema_ED::MatrizRestDemanda()	// Monta matriz da restrição de atendimento a demanda
{
	// T * B
	int n_a = n;
	CMatrizEsparsa Agh(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
			Agh.InserirElemento(b, sistema_a->barrasVtr[b].hidrosPtr[br]->GetIdentUsina(), 1);
	CMatrizEsparsa Agt(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
	for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
			Agt.InserirElemento(b, sistema_a->barrasVtr[b].termosPtr[bi]->GetIdentUsina(), 1);
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
CMatrizEsparsa CProblema_ED::MatrizRestDemandaBarraUnica()	// Monta matriz da restrição de atendimento a demanda
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
CMatrizEsparsa CProblema_ED::MatrizLimFluxo()
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
			Alin.InserirMatriz(ll, c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()), ll + sistema_a->linhasVtr.size() - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
			ll = ll + sistema_a->linhasVtr.size();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
	}
	else // flag1 == 2
	{
		////
		//// essas matrizes não são tão grandes para serem tratas como sparsas, porém Alin é
		//MatrixXd B = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->barrasVtr.size()));
		//MatrixXd St = MatrixXd::Zero(int(sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));
		//MatrixXd TT = MatrixXd::Zero(int(sistema_a->linhasVtr.size()), int(sistema_a->linhasVtr.size()));		
		//// algoritmo para formar a matriz de incidencia linhaxbarra e Ybarra matriz B
		//int ifr, ito;
		//for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		//{
		//	ifr = sistema_a->linhasVtr[l].GetIndDeBarra() - 1;
		//	ito = sistema_a->linhasVtr[l].GetIndParaBarra() - 1;
		//	St(l, ifr) = 1;
		//	St(l, ito) = -1;
		//	B(ifr, ifr) = B(ifr, ifr) + 100/(sistema_a->linhasVtr[l].GetReatancia());
		//	B(ifr, ito) = B(ifr, ito) - 100/(sistema_a->linhasVtr[l].GetReatancia());
		//	B(ito, ifr) = B(ito, ifr) - 100/(sistema_a->linhasVtr[l].GetReatancia());
		//	B(ito, ito) = B(ito, ito) + 100/(sistema_a->linhasVtr[l].GetReatancia());
		//	for ( size_t ll = 0; ll < sistema_a->linhasVtr.size(); ll++)
		//		if (l == ll)
		//			TT(l, ll) = 100 / (sistema_a->linhasVtr[l].GetReatancia());
		//}
		////
		////MatrixXd Beta(int(sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));
		//Beta = TT*St;
		//// Eliminar linha e coluna de referencia de B antes de invertê-la (em seguida incluí-se uma linha e uma coluna com 0's)
		//int refbus = sistema_a->GetBarraRef() - 1;
		//// Remover linha
		//int numRows = B.rows() - 1;
		//int numCols = B.cols();
		//B.block(refbus, 0, numRows - refbus, numCols) = B.block(refbus + 1, 0, numRows - refbus, numCols);
		//B.conservativeResize(numRows,numCols);
		//// Remover coluna
		//numRows = B.rows();
		//numCols = B.cols() - 1;
		//B.block(0, refbus, numRows, numCols - refbus) = B.block(0, refbus + 1, numRows, numCols - refbus);
		//B.conservativeResize(numRows,numCols);
		//// Inverter matriz
		////MatrixXd Ieye = MatrixXd::Identity(int(sistema_a->barrasVtr.size()), int(sistema_a->barrasVtr.size()));
		////B = B.householderQr().solve(Ieye);
		//B = B.inverse();
		//// Incluir linha de 0's na posição da barra de referencia
		//numRows = B.rows() + 1;
		//numCols = B.cols();
		//B.conservativeResize(numRows,numCols);
		//MatrixXd temp = B.block(refbus, 0, numRows - refbus - 1, numCols);
		//B.block(refbus + 1, 0, numRows - refbus - 1, numCols) = temp;
		//for (int i = 0; i < B.cols(); i++)
		//	B(refbus, i) = 0;
		//// Incluir coluna de 0's na posição da barra de referencia
		//numRows = B.rows();
		//numCols = B.cols() + 1;
		//B.conservativeResize(numRows,numCols);
		//temp = B.block(0, refbus, numRows, numCols - refbus - 1);
		//B.block(0, refbus + 1, numRows, numCols - refbus - 1) = temp;
		//for (int i = 0; i < B.rows(); i++)
		//	B(i, refbus) = 0;
		//Beta = Beta*B;		// Beta é uma matriz densa!!
		////
		
		// Montar matriz de incidencia de hidros (ph), térmicas (pt e u) e déficit (def): Agh, Agt, Agtu Adef
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
CMatrizEsparsa CProblema_ED::MatrizLimPhgMin()
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
CMatrizEsparsa CProblema_ED::MatrizLimPhgMax()
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
CMatrizEsparsa CProblema_ED::MatrizBalHid()
{
	// T * R
	int n_a = n;
	double deltaT;
	int tempo_viagem;
	int R = sistema_a->hidreletricasVtr.size();
	CMatrizEsparsa Addd(R * T, R);
	CMatrizEsparsa Add(R, R);

	// codigo abaixo tem problema quando o tempo de viagem é maior que T_2: corrigido, linhas 329 e 347 (corrigido no ED, na ScnDec3, SpcDec3 e Hrst)

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
			if ((tempo_viagem == t) && (sistema_a->hidreletricasVtr[r].GetUsinaJusante() != 0) && (tempo_viagem < sistema_a->GetTt2()) )
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
			//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			//{
			//	if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
			//	{
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
CMatrizEsparsa CProblema_ED::MatrizLimqMin()
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
CMatrizEsparsa CProblema_ED::MatrizLimqMax()
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
CMatrizEsparsa CProblema_ED::MatrizTup()
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
		// T * I (para usinas com Tup > 0)	nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < Tup
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
					//if  (t >= Tup + 1)
					if  (t_a >= Tup)
					{
						a.InserirElemento(0, c + i, -1);		// adiciona termo de u
						for (int tt = 0; tt <= Tup; tt++)
						{
							if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
								a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);
							else
								a.InserirElemento(0, ca + i + I - tt*(n_a / T), 1);
						}
					}
					else if ((Tup >= sistema_a->GetTt2()) && (t_a == sistema_a->GetTt2() - 1))
					{
						a.InserirElemento(0, c + i, -1);		// adiciona termo de u
						for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
						{
							if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
								a.InserirElemento(0, c + i + I - tt*(n_a / T), 1);
							else
								a.InserirElemento(0, ca + i + I - tt*(n_a / T), 1);
						}
					}
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
CMatrizEsparsa CProblema_ED::MatrizTdown()
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
		// T * I (para usinas com Tdown > 0)	nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < Tdown
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
					//if  (t >= Tdown + 1)
					if  (t_a >= Tdown)
					{
						a.InserirElemento(0, c + i, 1);
						for (int tt = 0; tt <= Tdown; tt++)
						{
							if (t_a - tt >= sistema_a->GetTt1())
								a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
							else
								a.InserirElemento(0, ca + i + 2*I - tt*(n_a / T), 1);
						}
					}
					else if ((Tdown >= sistema_a->GetTt2()) && (t_a == sistema_a->GetTt2() - 1))
					{
						a.InserirElemento(0, c + i, 1);
						for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
						{
							if (t_a - tt >= sistema_a->GetTt1())
								a.InserirElemento(0, c + i + 2*I - tt*(n_a / T), 1);
							else
								a.InserirElemento(0, ca + i + 2*I - tt*(n_a / T), 1);
						}
					}
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
CMatrizEsparsa CProblema_ED::MatrizTupDown()
{
	// T * I	(nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < 1
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
CMatrizEsparsa CProblema_ED::MatrizRampaUp()
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
CMatrizEsparsa CProblema_ED::MatrizRampaDown()
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
CMatrizEsparsa CProblema_ED::MatrizLimPtMin()
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
					//if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))		// final de horizonte, adicionam-se 2 restrições
					//{
					//	a.RemoverTodosElementos();
					//	a.InserirElemento(0, c + i, 1);		// pt
					//	a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u
					//	Alin.JuntarColuna(&a);
					//}
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
					//if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))		// final de horizonte, adicionam-se 2 restrições
					//{
					//	a.RemoverTodosElementos();
					//	a.InserirElemento(0, c + i, 1);		// pt
					//	a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u
					//	a.InserirElemento(0, c + i + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// up
					//	Alin.JuntarColuna(&a);
					//}
				}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
		return Alin;

		// metodo antigo com t+1

		//// (T - 1 + n_cen) * I (porém em t = 0 não são adicionadas restrições e em t = T2 são adicionadas 2 por usina)
		//int n_a = n;
		//n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		//size_t I = sistema_a->termeletricasVtr.size();
		//CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
		//CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
		//int ca, t_a, cen;
		//int c = 0;
		//c += (n_a / T);						// referente à t = 0
		//for (int t = 1; t < T; t++)			// nó inicial, não são adicionadas restrições, por isso t começa em 1
		//{
		//	if (t >= sistema_a->GetTt2())
		//	{
		//		cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		//		ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1()) - flag2*cen*sistema_a->hidreletricasVtr.size();
		//		t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		//	}
		//	else
		//	{
		//		ca = c;
		//		t_a = t;
		//	}

		//	for (size_t i = 0; i < I; i++)
		//	{
		//		if (sistema_a->termeletricasVtr[i].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
		//		{
		//			a.RemoverTodosElementos();
		//			a.InserirElemento(0, c + i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
		//			if (t_a - 1 >= sistema_a->GetTt1())		// todos periodos menos do inicio do segundo estágio (incluindo o de final de horizonte)
		//			{
		//				a.InserirElemento(0, c + i - (n_a / T), 1);		// pt-1
		//				a.InserirElemento(0, c + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u-1
		//			}
		//			else									// periodos no inicio do segundo estágio
		//			{
		//				a.InserirElemento(0, ca + i - (n_a / T), 1);
		//				a.InserirElemento(0, ca + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//			}
		//			Alin.JuntarColuna(&a);
		//			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))		// final de horizonte, adicionam-se 2 restrições
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c + i, 1);		// pt
		//				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u
		//				Alin.JuntarColuna(&a);
		//			}
		//		}
		//		else
		//		{
		//			a.RemoverTodosElementos();
		//			a.InserirElemento(0, c + i + 3*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// ud
		//			if (t_a - 1 >= sistema_a->GetTt1())		// todos periodos menos do inicio do segundo estágio (incluindo o de final de horizonte)
		//			{
		//				a.InserirElemento(0, c + i - (n_a / T), 1);		// pt-1
		//				a.InserirElemento(0, c + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u-1
		//				a.InserirElemento(0, c + i - (n_a / T) + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// up-1
		//			}
		//			else									// periodos no inicio do segundo estágio
		//			{
		//				a.InserirElemento(0, ca + i - (n_a / T), 1);
		//				a.InserirElemento(0, ca + i - (n_a / T) + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//				a.InserirElemento(0, ca + i - (n_a / T) + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//			}
		//			Alin.JuntarColuna(&a);
		//			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))		// final de horizonte, adicionam-se 2 restrições
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c + i, 1);		// pt
		//				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// u
		//				a.InserirElemento(0, c + i + 2*I, (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));		// up
		//				Alin.JuntarColuna(&a);
		//			}
		//		}
		//	}
		//	if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//		c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		//	else
		//		c += (n_a / T);
		//}
		//return Alin;

		// metodo mais antigo com t+1

		//// sum(n_aberturas) * I
		//int n_a = n;
		//n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		//size_t I = sistema_a->termeletricasVtr.size();
		//CMatrizEsparsa Alin(0, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
		//CMatrizEsparsa a(1, n_a + flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
		//int c = 0;
		//for (int t = 0; t < T; t++)
		//{
		//	for (size_t i = 0; i < I; i++)
		//	{
		//		if (sistema_a->termeletricasVtr[i].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
		//		{
		//			if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//			{
		//				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//				{
		//					a.RemoverTodosElementos();
		//					a.InserirElemento(0, c + i, 1);
		//					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//					a.InserirElemento(0, c + i + 3*I + (n_a / T)*(1 + n_cen*(sistema_a->GetTt2() - sistema_a->GetTt1())) + flag2*n_cen*sistema_a->hidreletricasVtr.size(), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//					Alin.JuntarColuna(&a);
		//				}
		//			}
		//			else if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))		// final de horizonte
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c + i, 1);
		//				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//				Alin.JuntarColuna(&a);
		//			}
		//			else		// periodos restantes
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c + i, 1);
		//				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//				a.InserirElemento(0, c + i + 3*I + (n_a / T), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//				Alin.JuntarColuna(&a);
		//			}
		//		}
		//		else
		//		{
		//			if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//			{
		//				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//				{
		//					a.RemoverTodosElementos();
		//					a.InserirElemento(0, c + i, 1);
		//					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//					a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//					a.InserirElemento(0, c + i + 3*I + (n_a / T)*(1 + n_cen*(sistema_a->GetTt2() - sistema_a->GetTt1())) + flag2*n_cen*sistema_a->hidreletricasVtr.size(), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//					Alin.JuntarColuna(&a);
		//				}
		//			}
		//			else if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c + i, 1);
		//				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//				a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//				Alin.JuntarColuna(&a);
		//			}
		//			else
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c + i, 1);
		//				a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
		//				a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//				a.InserirElemento(0, c + i + 3*I + (n_a / T), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
		//				Alin.JuntarColuna(&a);
		//			}
		//		}
		//	}
		//	if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//		c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		//	else
		//		c += (n_a / T);
		//}
		//return Alin;
	}
}
CMatrizEsparsa CProblema_ED::MatrizLimPtMax()
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
				//if (sistema_a->termeletricasVtr[i].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
				//{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i, 1);
					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
					Alin.JuntarColuna(&a);
				//}
				//else
				//{
				//	a.RemoverTodosElementos();		// restrição já adicionada na função ptmin
				//	Alin.JuntarColuna(&a);
				//}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa CProblema_ED::MatrizLimPtMinMax()	// nao esta sendo usado
{
	// (T*n_aberturas) * I ????
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
			if (sistema_a->termeletricasVtr[i].GetTUp() == 0)		// units that can be online for just one period
			{
				if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
				{
					for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					{
						a.RemoverTodosElementos();
						a.InserirElemento(0, c + i, 1);
						a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
						a.InserirElemento(0, c + i + 3*I + (n_a / T)*(1 + n_cen*(sistema_a->GetTt2() - sistema_a->GetTt1())) + flag2*n_cen*sistema_a->hidreletricasVtr.size(), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
						Alin.JuntarColuna(&a);
					}
				}
				else if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				{
					a.RemoverTodosElementos();
					//a.InserirElemento(0, c + i, 1);
					//a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					Alin.JuntarColuna(&a);
				}
				else
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i, 1);
					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					a.InserirElemento(0, c + i + 3*I + (n_a / T), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
					Alin.JuntarColuna(&a);
				}
			}
			else
			{
				if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
				{
					for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					{
						a.RemoverTodosElementos();
						a.InserirElemento(0, c + i, 1);
						a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
						a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
						a.InserirElemento(0, c + i + 3*I + (n_a / T)*(1 + n_cen*(sistema_a->GetTt2() - sistema_a->GetTt1())) + flag2*n_cen*sistema_a->hidreletricasVtr.size(), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
						Alin.JuntarColuna(&a);
					}
				}
				else if (t == 0)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i, 1);
					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					a.InserirElemento(0, c + i + 3*I + (n_a / T), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
					Alin.JuntarColuna(&a);
				}
				else if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				{
					a.RemoverTodosElementos();
					Alin.JuntarColuna(&a);
				}
				else
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c + i, 1);
					a.InserirElemento(0, c + i + I, - (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					a.InserirElemento(0, c + i + 2*I, sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
					a.InserirElemento(0, c + i + 3*I + (n_a / T), sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin());
					Alin.JuntarColuna(&a);
				}
			}
		}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa CProblema_ED::MatrizRestCP()
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
CMatrizEsparsa CProblema_ED::MatrizVmeta()
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
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	//for (int cen = sistema_a->GetNCenarios() - 1 ; cen >= 0; cen--)
	//{
	//	for (size_t r = 0; r < R; r++)
	//	{	
	//		c = r + (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + R + (sistema_a->GetNCenarios() - 1 - cen)*R;
	//		c = c + (T - 1)*(n_a / T) - cen * (sistema_a->GetTt2() - sistema_a->GetTt1())*(n_a / T);
	//		a.RemoverTodosElementos();
	//		a.InserirElemento(0, c, 1);
	//		if (sistema_a->GetFlagVfol() == true)
	//			a.InserirElemento(0, c + 4*R + 3*JJ + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d), 1);
	//		Alin.JuntarColuna(&a);
	//	}
	//}
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
CMatrizEsparsa CProblema_ED::MatrizBalPotencia()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
CMatrizEsparsa CProblema_ED::MatrizBalVazao()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
CMatrizEsparsa CProblema_ED::MatrizFuncProd()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
					a.InserirElemento(0, jj - JJ + j, 1);																		// phg
					a.InserirElemento(0, vv + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgV(napr));					// v
					a.InserirElemento(0, jj + j, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgQ(napr));					// q
					a.InserirElemento(0, vv + R + r, - sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhgD(napr));				// d
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
CMatrizEsparsa CProblema_ED::MatrizPhMax()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
CMatrizEsparsa CProblema_ED::MatrizReserva()
{
	int n_a = n;
	n_a = n_a - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
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
CMatrizEsparsa CProblema_ED::MatrizCortesF()
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
vetorfloat2 CProblema_ED::LimRestDemanda()
{
	size_t B = sistema_a->barrasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, B * T, 2);
	IniciaMatriz(&Lim, 0);
	int cen = 0;
	for (int t = 0; t < T; t++)
	{
		for (size_t b = 0; b < B; b++)
		{
			double D = 0;
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(0, t);		// cen = 0, ou qq outro cenários, pois os T1 primeiros periodos são iguais (deterministicos)
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
					D = D + sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
			}
		//Lim[t*B + b][0] = 2;
		Lim[t*B + b][0] = 1;
		Lim[t*B + b][1] = D;
		}
	}
	return Lim;
}
vetorfloat2 CProblema_ED::LimRestDemandaBarraUnica()
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
vetorfloat2 CProblema_ED::LimFluxo0()
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
vetorfloat2 CProblema_ED::LimFluxo2()
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
vetorfloat2 CProblema_ED::LimPhgMin()
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
vetorfloat2 CProblema_ED::LimPhgMax()
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
vetorfloat2 CProblema_ED::LimBalHid()
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
vetorfloat2 CProblema_ED::LimQMin()
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
vetorfloat2 CProblema_ED::LimQMax()
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
vetorfloat2 CProblema_ED::LimTup()
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
vetorfloat2 CProblema_ED::LimTdown()
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
vetorfloat2 CProblema_ED::LimTupDown()
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
vetorfloat2 CProblema_ED::LimRampaUp()
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
vetorfloat2 CProblema_ED::LimRampaDown()
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
vetorfloat2 CProblema_ED::LimPtMin()
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
		// t = 0
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

		//// (T - 1 + n_cen) * I (porém em t = 0 não são adicionadas restrições e em t = T2 são adicionadas 2 por usina)
		//int I = sistema_a->termeletricasVtr.size();
		//vetorfloat2 Lim;
		//vetorfloat L; L.resize(2);
		//L[0] = 0; L[1] = 0;
		//for (int t = 1; t < T; t++)
		//{
		//	for (size_t i = 0; i < I; i++)
		//	{
		//		Lim.push_back(L);
		//		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))		// final de horizonte, adicionam-se 2 restrições
		//			Lim.push_back(L);
		//	}
		//}
		//return Lim;

		//// sum(t*n_aberturas) * I
		//int I = sistema_a->termeletricasVtr.size();
		//vetorfloat2 Lim;
		//vetorfloat L; L.resize(2);
		//L[0] = 0; L[1] = 0;
		//for (int t = 0; t < T; t++)
		//{
		//	for (size_t i = 0; i < I; i++)
		//	{
		//		if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//		{
		//			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//				Lim.push_back(L);
		//		}
		//		else
		//			Lim.push_back(L);
		//	}
		//}
		//return Lim;
	}
}
vetorfloat2 CProblema_ED::LimPtMax()
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
vetorfloat2 CProblema_ED::LimRestCP()
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
vetorfloat2 CProblema_ED::LimVmeta()
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
vetorfloat2 CProblema_ED::LimBalPotenciaL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 CProblema_ED::LimBalVazaoL()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 CProblema_ED::LimFuncProdL()
{
	int R = sistema_a->hidreletricasVtr.size();
	int JJ = 0;
	//int total_app = 0;
	for (int i = 0; i < R; i++)
	{
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
		//for (int j = 0; j < sistema_a->hidreletricasVtr[i].GetNGrupos(); j++)
		//	total_app += sistema_a->hidreletricasVtr[i].grupoVtr[j].GetNappFPH();
	}
	vetorfloat2 Lim;
	//DimensionarMatriz(&Lim, JJ * T * total_app, 2);
	vetorfloat L;
	L.resize(2);
	//int jj = 0;
	for (int t = 0; t < T; t++)
	{
		for (int r = 0; r < R; r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			{
				for (int napr = 0; napr < sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNappFPH(); napr++)
				{
					L[0] = 0;
					L[1] = sistema_a->hidreletricasVtr[r].grupoVtr[j].CoefPhg(napr) + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax();
					Lim.push_back(L);
				}
			}
			//jj = jj + sistema_a->hidreletricasVtr[r].GetNGrupos();
		}
	}
	return Lim;
}
vetorfloat2 CProblema_ED::LimPhMax()
{
	int R = sistema_a->hidreletricasVtr.size();
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, R * T, 2);
	IniciaMatriz(&Lim, 0);
	for (size_t i = 0; i < Lim.size(); i++)
		Lim[i][0] = 1;
	return Lim;
}
vetorfloat2 CProblema_ED::LimReserva()
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
vetorfloat2 CProblema_ED::LimCortesF()
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
// Caso seja necessário fazer modelagens diferentes para cada estágio é só criar uma condição que depende do T dentro da função
void CProblema_ED::MatrizRestricoesLineares(CMatrizEsparsa &MM)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	if ( sistema_a->GetFlagModeloRede() == 1 )
	{
		M = MatrizRestDemanda();
		//SparseMatriz(&M, indexL, indexC, indexV, nnZ, n_restricoes, &count);
		//*n_restricoes = *n_restricoes + M.GetNlin(); M.ZerarMatriz();
		MM.JuntarColuna(&M); M.ZerarMatriz();
		M = MatrizLimFluxo();
		MM.JuntarColuna(&M); M.ZerarMatriz();
		M = MatrizLimFluxo();
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
	else	// == 0 e == 3(nesse caso os limites de fluxo somente são adicionados se violados)
	{
		M = MatrizRestDemandaBarraUnica();
		MM.JuntarColuna(&M); M.ZerarMatriz();
		if ( sistema_a->GetFlagModeloRede() == 2 )
		{
			M = MatrizLimFluxo();
			MM.JuntarColuna(&M); M.ZerarMatriz();
			M = MatrizLimFluxo();
			MM.JuntarColuna(&M); M.ZerarMatriz();
		}
	}
	M = MatrizLimPhgMin();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPhgMax();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimqMin();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimqMax();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalHid();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizTup();									
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizTdown();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if ( sistema_a->GetFlagTbinaryModel() == 1 )
	{
		M = MatrizTupDown();
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
	M = MatrizRampaUp();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizRampaDown();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPtMax();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPtMin();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		M = MatrizRestCP();
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
	M = MatrizVmeta();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalPotencia();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizBalVazao();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizFuncProd();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if (sistema_a->GetFlagPhmax() == 1)
	{
		M = MatrizPhMax();
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
	M = MatrizReserva();
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		M = MatrizCortesF();
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
}
void CProblema_ED::MatrizLimitesLineares(vetorint * LimTipo, vetorfloat * LimValor)
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

CProblema_ED::CProblema_ED(CSistema * const sistema_ptr) :ambGRB(GRBEnv()),modeloGRB(GRBModel(ambGRB))
{
	sistema_a = sistema_ptr;
	resultados = new Resultados_ED(sistema_a);
	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	flag1 = int (sistema_a->GetFlagModeloRede());
	flag2 = int (sistema_a->GetFlagVfol());
	flag3 = int (sistema_a->GetFlagPhmax());
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
	Cvfol = sistema_a->GetCustoVfol();
	///Cvfol = 20*sistema_a->GetCustoDeficit();
	// os problemas estendidos (ext) usam Cvfol = (20/0.0036)*sistema_a->GetCustoDeficit()

	// Criar variáveis
	CriarVariaveis();
	modeloGRB.update();	// Atualiza o modelo Gurobi.
	vars = modeloGRB.getVars();
	
	//if (flag1 >= 2)
	//	Beta.resize(int(sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));		// matriz usada no modelo compacto
	//else
	//	Beta.resize(0,0);

	// Adicionar restrições
	CriarRestricoes();
	modeloGRB.update();

	// Fixar condições iniciais
	FixarCondIniciais();
	modeloGRB.update();
}
CProblema_ED::~CProblema_ED(void)
{
	delete vars;
	delete resultados;
}

void CProblema_ED::CriarVariaveis()
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
					modeloGRB.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()), 0.0, GRB_CONTINUOUS, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)	//pt
					modeloGRB.addVar(0, double (sistema_a->termeletricasVtr[i].GetPmax()), 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagVarBin() == true)	//u (binario ou continuo)
			{
				for (size_t i = 0; i < I; i++)
					modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t i = 0; i < I; i++)
					modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
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
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F		//if (sistema_a->GetFlagInitAproxCT() > 1)
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
				//modeloGRB.addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
			}
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				modeloGRB.addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t)), double(sistema_a->hidreletricasVtr[r].GetVmax(t)), 0.0, GRB_CONTINUOUS, "");
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
				//modeloGRB.addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
			}
			if (sistema_a->GetFlagPhmax() == 1)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
				{
					double phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					modeloGRB.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
				}
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
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modeloGRB.addVar(0, 1, 0.0, GRB_BINARY, "");
			}
			else
			{
				for (size_t r = 0; r < R; r++)
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
						modeloGRB.addVar(0, 1, 0.0, GRB_CONTINUOUS, "");
			}
			if ( sistema_a->GetFlagModeloRede() > 0)
			{
				for (size_t b = 0; b < B + 1; b++)	//def
				{
					double cap_d = 0;
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(0, t);
					modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			else
			{
				double cap_d = 0;			//def barra unica
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
						cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(0, t);
				modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
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
					//modeloGRB.addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					modeloGRB.addVar(double(sistema_a->hidreletricasVtr[r].GetVmin(t+sistema_a->GetTt1()+cen*T2)), double(sistema_a->hidreletricasVtr[r].GetVmax(t+sistema_a->GetTt1()+cen*T2)), 0.0, GRB_CONTINUOUS, "");
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
					//modeloGRB.addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
				}
				if (sistema_a->GetFlagPhmax() == 1)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
					{
						double phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax = phmax + sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						modeloGRB.addVar(0, double (phmax), 0.0, GRB_CONTINUOUS, "");
					}
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
							cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t + sistema_a->GetTt1());
						modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
					}
				}
				else
				{
					double cap_d = 0;			//def barra unica
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						for (size_t db = 0; db < sistema_a->barrasVtr[b].demandasPtr.size(); db++)
							cap_d = cap_d + sistema_a->barrasVtr[b].demandasPtr[db]->GetD(cen, t + sistema_a->GetTt1());
					modeloGRB.addVar(0, cap_d, 0.0, GRB_CONTINUOUS, "");
				}
			}
			if ( sistema_a->GetFlagVfol() == true )
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//vfol
					modeloGRB.addVar(0, double(sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()), 0.0, GRB_CONTINUOUS, "");
					//modeloGRB.addVar(0, 0, 0.0, GRB_CONTINUOUS, "");
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
void CProblema_ED::CriarRestricoes()
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
		MM.ZerarMatriz();
		LimTipo.clear();
		LimValor.clear();
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
void CProblema_ED::CriarFuncaoObjetivo()
{
	try 
	{
		GRBQuadExpr fo;
		//x = [pt u cp F teta ph v d s phmax phg q z def vfol]
		double coeficiente;
		GRBVar variavel;
		// Termos lineares
		int nt;
		int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		double deltaT;
		int delta = 0;
		int cen = 0;
		int flag1a = 0;		// referente à var. teta
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
						variavel = vars[i + 3*sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					else
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;		//pt
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT;		//u
						variavel = vars[i + sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					coeficiente = 1;		//cp
					variavel = vars[i + 2*sistema_a->termeletricasVtr.size() + delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
				//for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				//{
				//	coeffs[r + (3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 3*sistema_a->hidreletricasVtr.size() + delta ] = 1;		//s: 1e-6
				//}
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
						variavel = vars[i + 3*sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					else
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//pt
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
						variavel = vars[i + sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					coeficiente = 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//cp
					variavel = vars[i + 2*sistema_a->termeletricasVtr.size() + delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
				//for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				//{
				//	coeffs[r + (3+flag4)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 3*sistema_a->hidreletricasVtr.size() + delta ] = 1e-6;		//s
				//}
			}
			delta = delta + nt;
		}
		
		// Deficit e vfolga
		int JJ = 0;
		for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
			JJ = JJ + sistema_a->hidreletricasVtr[i].GetNGrupos();
		delta = (3+flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ;
		for (int t = 0; t < T; t++)
		{
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
			else
				nt = (n_a / T);
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				deltaT = sistema_a->GetDeltaT1();
				if (sistema_a->GetFlagModeloRede() > 0 )
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					{
						coeficiente = sistema_a->GetCustoDeficit()*deltaT;		//def
						variavel = vars[delta + b];
						fo.addTerms(&coeficiente, &variavel, 1);
						///coeffs[delta + b] = sistema_a->GetCustoDeficit()*deltaT;		//def
					}
				}
				else
				{
					coeficiente = sistema_a->GetCustoDeficit()*deltaT;		//def
					variavel = vars[delta];
					fo.addTerms(&coeficiente, &variavel, 1);
					///coeffs[delta] = sistema_a->GetCustoDeficit()*deltaT;		//def
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
						variavel = vars[delta + b];
						fo.addTerms(&coeficiente, &variavel, 1);
						///coeffs[delta + b] = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
					}
				}
				else
				{
					coeficiente = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
					variavel = vars[delta];
					fo.addTerms(&coeficiente, &variavel, 1);
					///coeffs[delta] = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
				}
				if (sistema_a->GetFlagVfol() == true)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					{
						if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
						{
							coeficiente = Cvfol*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
							variavel = vars[delta + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d) + r];
							fo.addTerms(&coeficiente, &variavel, 1);
							///coeffs[delta + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d) + r ] = Cvfol*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
						}
					}
				}
			}
			delta = delta + nt;
		}

		// Termos quadráticos
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			delta = 0;
			///for (int i = 0; i < n; i++)
			///	coeffs[i] = 0;
			for (int t = 0; t < T; t++)
			{
				if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
					nt = (n_a / T) + int(flag2*sistema_a->hidreletricasVtr.size());
				else
					nt = (n_a / T);
				if ((0 <= t) && (sistema_a->GetTt1() > t))
				{
					deltaT = sistema_a->GetDeltaT1();
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;		//pt^2
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel, 1);
						///coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;		//pt^2
					}
				}
				else
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					deltaT = sistema_a->GetDeltaT2();
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//pt^2
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel, 1);
						///coeffs[i + delta] = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//pt^2
					}
				}
				delta = delta + nt;
			}
			///fo.addTerms(coeffs, vars, vars, n);
		}
		///delete coeffs;
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
void CProblema_ED::CriarFuncaoObjetivo2()
{
	try 
	{
		//GRBExpr * fo;
		//if (sistema_a->GetFlagInitAproxCT() == 0)
		//	fo = new GRBQuadExpr;
		//else
		//	fo = new GRBLinExpr;

		GRBQuadExpr fo;
		//x = [pt u up ud F teta ph v d s phmax phg q z def vfol]
		double coeficiente;
		GRBVar variavel, variavel2;
		// Termos lineares
		int nt;
		int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
		double deltaT;
		int delta = 0;
		int cen = 0;
		int flag1a = 0;		// referente à var. teta
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
						variavel = vars[i + 4*sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					else
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;		//pt
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
						if (sistema_a->GetFlagInitAproxCT() == 0)	// termo quadrático de u (u^2 = u), portanto ele é adicionado aqui
						{
							coeficiente = deltaT*(sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetPmin()*(sistema_a->termeletricasVtr[i].GetCoefCustoOper(1) + sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)));		//u
							variavel = vars[i + sistema_a->termeletricasVtr.size() + delta];
							fo.addTerms(&coeficiente, &variavel, 1);
						}
						else
						{
							coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT + sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT;		//u
							variavel = vars[i + sistema_a->termeletricasVtr.size() + delta];
							fo.addTerms(&coeficiente, &variavel, 1);
						}
					}
					coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		//up
					variavel = vars[i + 2*sistema_a->termeletricasVtr.size() + delta];
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
						variavel = vars[i + 4*sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
					else
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//pt
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, 1);
						if (sistema_a->GetFlagInitAproxCT() == 0)	// termo quadrático de u (u^2 = u), portanto ele é adicionado aqui
						{
							coeficiente = deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)*(sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[i].GetPmin()*(sistema_a->termeletricasVtr[i].GetCoefCustoOper(1) + sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)));		//u
							variavel = vars[i + sistema_a->termeletricasVtr.size() + delta];
							fo.addTerms(&coeficiente, &variavel, 1);
						}
						else
						{
							coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) + sistema_a->termeletricasVtr[i].GetPmin()*sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//u
							variavel = vars[i + sistema_a->termeletricasVtr.size() + delta];
							fo.addTerms(&coeficiente, &variavel, 1);
						}
					}
					coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoPartida() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//up
					variavel = vars[i + 2*sistema_a->termeletricasVtr.size() + delta];
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
						variavel = vars[delta + b];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
				}
				else
				{
					coeficiente = sistema_a->GetCustoDeficit()*deltaT;		//def
					variavel = vars[delta];
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
						variavel = vars[delta + b];
						fo.addTerms(&coeficiente, &variavel, 1);
					}
				}
				else
				{
					coeficiente = sistema_a->GetCustoDeficit()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//def
					variavel = vars[delta];
					fo.addTerms(&coeficiente, &variavel, 1);
				}
				if (sistema_a->GetFlagVfol() == true)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					{
						if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
						{
							coeficiente = Cvfol*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//vfol
							variavel = vars[delta + flag1d*(sistema_a->barrasVtr.size()) + (1 - flag1d) + r];
							fo.addTerms(&coeficiente, &variavel, 1);
						}
					}
				}
			}
			delta = delta + nt;
		}

		// Termos quadráticos
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			delta = 0;
			for (int t = 0; t < T; t++)
			{
				if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
					nt = (n_a / T) + int(flag2*sistema_a->hidreletricasVtr.size());
				else
					nt = (n_a / T);
				if ((0 <= t) && (sistema_a->GetTt1() > t))
				{
					deltaT = sistema_a->GetDeltaT1();
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;		//pt^2
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel, 1);
					}
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						coeficiente = 2*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT*sistema_a->termeletricasVtr[i].GetPmin();		//pt*u
						variavel = vars[i + delta];
						variavel2 = vars[i + sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel2, 1);
					}
				}
				else
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					deltaT = sistema_a->GetDeltaT2() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						coeficiente = sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT;		//pt^2
						variavel = vars[i + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel, 1);
					}
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						coeficiente = 2*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*deltaT*sistema_a->termeletricasVtr[i].GetPmin();		//pt*u
						variavel = vars[i + delta];
						variavel2 = vars[i + sistema_a->termeletricasVtr.size() + delta];
						fo.addTerms(&coeficiente, &variavel, &variavel2, 1);
					}
				}
				delta = delta + nt;
			}
		}

		//// Termos quadráticos
		//if (sistema_a->GetFlagInitAproxCT() == 0)
		//	cout << "CriarFuncaoObjetivo2: Termo quadrático não implementado para essa modelagem!!" << endl;

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
void CProblema_ED::FixarCondIniciais()
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
					vars[i + I + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[i].GetU0());
					vars[i + I + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[i].GetU0());
					//cp
					vars[i + 2*I + c].set(GRB_DoubleAttr_LB, 0);
					vars[i + 2*I + c].set(GRB_DoubleAttr_UB, 0);
				}
			}
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			c += (n_a / T);
		}
	}
}
int CProblema_ED::ResolverProblema()
{
	try
	{
		// Criar função objetivo
		if (flag7 == 0)
			CriarFuncaoObjetivo();
		else
			CriarFuncaoObjetivo2();
		modeloGRB.update();	// Atualiza o modelo Gurobi.

		//modeloGRB.write("probDE.lp");	// Escreve modelo em arquivo
		//modeloGRB.write("probDE.mps");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		//modeloGRB.getEnv().set(GRB_IntParam_OutputFlag, 0);

		modeloGRB.getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		modeloGRB.getEnv().set(GRB_StringParam_LogFile, "log_ED.txt");		// Escrever log em arquivo

		//modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.0001);	// Define o gap de tolerancia
		//modeloGRB.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9); 
		//modeloGRB.getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modeloGRB.getEnv().set(GRB_IntParam_Threads, 4);	//12 para o servidor
		modeloGRB.getEnv().set(GRB_IntParam_Method, 2);
		modeloGRB.getEnv().set(GRB_IntParam_Crossover, 0);

		//modeloGRB.getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modeloGRB.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9); 
		//modeloGRB.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);
		//modeloGRB.getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modeloGRB.getEnv().set(GRB_IntParam_PreSparsify, 1);

		//modeloGRB.getEnv().set(GRB_DoubleParam_IterationLimit, 100000);	// Limita o numero de iterações do método simplex
		//modeloGRB.getEnv().set(GRB_DoubleParam_NodefileStart, 0.5);
		//modeloGRB.getEnv().set(GRB_IntParam_MIPFocus, 2);
		//modeloGRB.getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

		//modeloGRB.getEnv().set(GRB_DoubleParam_TimeLimit, 600*pow(sistema_a->GetNCenarios(), 0.5));		// Limita o tempo de resolução do problema 36000
		if (sistema_a->GetNCenarios() == 1)
		{
			modeloGRB.getEnv().set(GRB_DoubleParam_TimeLimit, 1200);
			modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.001);
		}
		else
		{
			modeloGRB.getEnv().set(GRB_DoubleParam_TimeLimit, 7200);
			modeloGRB.getEnv().set(GRB_DoubleParam_MIPGap, 0.01);
		}

		// modeloGRB.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-2);	// Define tolerancia de otimalidade
		//modeloGRB.reset();
		
		// Create a callback object and associate it with the model
		CallbackED * cb;
		if ((sistema_a->GetFlagMaxAproxCT() >= 2) || (sistema_a->GetFlagModeloRede() == 3))
		{
			//modeloGRB.getEnv().set(GRB_IntParam_PreCrush, 1);
			modeloGRB.getEnv().set(GRB_IntParam_LazyConstraints, 1);
			cb = new CallbackED(vars, n, sistema_a);
			modeloGRB.setCallback(cb);
		}

		// Otimizar
		modeloGRB.optimize();
		
		//modeloGRB.write("probDE.lp");	

		// Salvar resultados
		int nStatus = modeloGRB.get(GRB_IntAttr_Status);
		double gap;
		if (sistema_a->GetFlagVarBin() == 1)
			gap = modeloGRB.get(GRB_DoubleAttr_MIPGap);
		else
			gap = 0;
		cout << modeloGRB.get(GRB_IntAttr_Status) << " \n";
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
			modeloGRB.computeIIS();				
			modeloGRB.write("probDE.lp");	
			modeloGRB.write("probDE.ilp");

			for (int i = 0; i < n; i++)
				x[i] = 0;
			fo = 0;
		}
		resultados->CarregarResultadosED(fo, x, L, nStatus, gap);

	//// Escrever x
	//ofstream * inFile;
	//inFile = new ofstream( "x.txt", ios::out );
	//if ( inFile->is_open() )                                                            
	//{
	//	*inFile << std::scientific << setprecision(10);
	//	for (size_t i = 0; i < n; i++)
	//		*inFile << x[i] << endl;
	//}
	//else
	//	cout << "Unable to open file";
	//inFile->close();

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultados->CarregarResultadosED(fo, x, L, e.getErrorCode(), 1);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
	return 0;
}
void CProblema_ED::ConferirCoerenciaDados()
{
	// Demanda e reserva girante por periodo
	double dem_b, dem_max, dem_min;
	dem_max = 0;
	dem_min = 999999;
	int cen;
	vetorfloat dem_cen;
	dem_cen.resize(sistema_a->GetNCenarios() + 1);
	cout << std::scientific << setprecision(9) << "Demanda e reserva girante por periodo" << endl;
	for (int t = 0; t < T; t++)
	{
		dem_b = 0;
		if ((0 <= t) && (sistema_a->GetTt1() > t))
			cen = 0;
		else
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
				dem_b += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
		if (dem_b >= dem_max)
			dem_max = dem_b;
		if (dem_b <= dem_min)
			dem_min = dem_b; 

		if ((0 <= t) && (sistema_a->GetTt1() > t))
			dem_cen[cen] += dem_b;
		else
			dem_cen[cen + 1] += dem_b;
		cout << "Periodo: " << t << " : " << dem_b << " Reserva girante: " << sistema_a->GetReserva(t) << endl;
	}
	cout << endl;
	cout << "Demanda maxima: " << dem_max << " Demanda minima: " << dem_min << endl;

	// Capacidade instalada
	double cap_h = 0;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			cap_h += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
	double cap_t = 0;
	for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
		cap_t += sistema_a->termeletricasVtr[i].GetPmax();
	cout << "Capacidade instalada" << endl;
	cout << "Hidro: " << cap_h << " Termo: " << cap_t << " Total: " << cap_h + cap_t << endl;

	// Energia capaz de ser produzida com as hidros, sem otimização (considerando que turbina-se o máximo até que o reservatório fique no nível da meta)
	// Nos primeiros periodos turbina-se o máximo até que v = vmeta, a partir dai turbina-se somente a afluencia  (respeitando-se o limite mínimo)
	double EH_total = 0;
	double ger_min;
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
	{
		vetorfloat vaz_disp;		// vazao disponivel total para geração (por cenário)
		vaz_disp.resize(sistema_a->GetNCenarios());
		// total de afluencia por cenário
		vetorfloat t_aflu;		// vetor com soma das afluencias no primeiro estágio (elemento 0) e no segundo estágio (elementos 1,...,N, em que cada um representa cada cenário)
		for (int i = 0; i < sistema_a->GetNCenarios() + 1; i++)
			t_aflu.push_back(0);
		for (int t = 0; t < T; t++)
		{
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				cen = 0;
				t_aflu[0] += sistema_a->hidreletricasVtr[r].GetAfluencia(cen, t);
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				t_aflu[cen + 1] += sistema_a->hidreletricasVtr[r].GetAfluencia(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
			}
		}
		// Verificar se a afluencia será o suficiente: calcular vazão disponível
		// quando meta > v0
		if (sistema_a->hidreletricasVtr[r].GetVMeta() > sistema_a->hidreletricasVtr[r].GetV0())
		{
			// afluencia total necessária para atingir a meta: (considerando periodos de 1 hora)
			double aflu_nec = (sistema_a->hidreletricasVtr[r].GetVMeta() - sistema_a->hidreletricasVtr[r].GetV0())/0.0036*sistema_a->GetTt2();
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			{
				if ( t_aflu[0] + t_aflu[n_cen+1] <= aflu_nec )
				{
					cout << "Usina " << r << " , afluencia insuficiente para atingir a meta no cenario " << n_cen << endl;
					vaz_disp[n_cen] = 0;
				}
				else	// calcular vazão disponivel (por cenário) que poderia ser usada para gerar energia
					vaz_disp[n_cen] = t_aflu[0] + t_aflu[n_cen+1] - aflu_nec;		// água disponivel
			}
		}
		else	// quando meta <= v0
		{
			// afluencia restante descontando a meta: (considerando periodos de 1 hora)
			double aflu_res = (sistema_a->hidreletricasVtr[r].GetV0() - sistema_a->hidreletricasVtr[r].GetVMeta())/(0.0036);
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
				vaz_disp[n_cen] = t_aflu[0] + t_aflu[n_cen+1] + aflu_res;
		}

		// escrever q a energia e capac. das termicas depende das restrições, portanto é menor!!
		// o calculo da energia das hidros tb é aprox. usando á agua disponivel e considerando que qdo pode-se gerar na vazao maxima gera-se a potencia maxima

		// alocar a vazão disponível nos periodos (gerando sempre com a vazão máxima) até o final do horizonte ou acabar a água
		// usando vetor vaz_disp por cenário (essa podendo ser = 0)
		// Simplificação: fphg (v,qmax,d,s) = phgmax!!

		double qhmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			qhmax += double (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
		double phmax = 0;
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
			phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
		ger_min = 999999;		// geração minima entre os cenários
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		{
			double ger = 0;
			if (vaz_disp[n_cen] / qhmax >= sistema_a->GetTt2())
				ger = sistema_a->GetTt2() * phmax;
			else
				ger = vaz_disp[n_cen] / qhmax * phmax;		// qts periodos se usar qhmax, a geração máximo é portanto porcentagem dessas
			if (ger <= ger_min)		// considera-se a menor geração entre os cenários possiveis
				ger_min = ger;
		}
		EH_total += ger_min;
	}
	cout << endl;
	// Energias (demanda, hidro e termo)
	double ener_dem_max = 0;
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	{
		if (dem_cen[0] + dem_cen[n_cen] >= ener_dem_max)
			ener_dem_max = dem_cen[0] + dem_cen[n_cen];
		cout << "Cenario: " << n_cen << " Energia: " << dem_cen[0] + dem_cen[n_cen] << endl;
	}
	cout << "Energia"  << endl;
	cout << "Hidro: " << EH_total << " Termo: " << cap_t*sistema_a->GetTt2() << " Demanda: " << ener_dem_max << endl;
	cout << endl;
}
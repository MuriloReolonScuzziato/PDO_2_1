#include "Spcdec3Results.h"

/*CMatrizEsparsa Spcdec3Results::MatrizCalcFluxo(int n_a)
{
	CMatrizEsparsa Alb(int (sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			if (sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[b])
				Alb.InserirElemento(int (l), int(b), 1);
			else if (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[b])
				Alb.InserirElemento(int (l), int(b), -1);			
	Alb.RemoverColuna(sistema_a->GetBarraRef() - 1);
	CMatrizEsparsa TT(int (sistema_a->linhasVtr.size()), int (sistema_a->linhasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
	{
		for ( size_t ll = 0; ll < sistema_a->linhasVtr.size(); ll++)
		{
			if (l == ll)
			{
				TT.InserirElemento(int(l), int(ll), 100 / (sistema_a->linhasVtr[l].GetReatancia()));
			}
		}
	}
	CMatrizEsparsa Alin(int (N * sistema_a->linhasVtr.size()), n_a);
	n_a = n_a - int(flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < N; t++)
	{
		Alin.InserirMatriz(ll,c + int ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),ll + int (sistema_a->linhasVtr.size()) - 1,c + int((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
		ll = ll + int(sistema_a->linhasVtr.size());
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / N) + int(flag2*sistema_a->hidreletricasVtr.size());
		else
			c += (n_a / N);
	}
	return Alin;
}
CMatrizEsparsa Spcdec3Results::MatrizCalcFluxo_cen(int n_a)
{
	CMatrizEsparsa Alb(int(sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			if (sistema_a->linhasVtr[l].de_barra == &sistema_a->barrasVtr[b])
				Alb.InserirElemento(int(l), int(b), 1);
			else if (sistema_a->linhasVtr[l].para_barra == &sistema_a->barrasVtr[b])
				Alb.InserirElemento(int(l), int(b), -1);
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
	CMatrizEsparsa Alin(N * sistema_a->linhasVtr.size(),n_a);
	n_a = n_a - sistema_a->hidreletricasVtr.size();
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < N; t++)
	{
		Alin.InserirMatriz(ll,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),ll + sistema_a->linhasVtr.size() - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
		ll = ll + sistema_a->linhasVtr.size();
		c += (n_a / N);
	}
	return Alin;
}*/

Spcdec3Results::Spcdec3Results(CSistema * const sistema_end)
{
	sistema_a = sistema_end;
	inFile = NULL;
	flag1 = int (sistema_a->GetFlagModeloRede());
	flag2 = int (sistema_a->GetFlagVfol());
	flag3 = int (sistema_a->GetFlagPhmax());
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = sistema_a->GetFlagTbinaryModel();
	
	R = sistema_a->hidreletricasVtr.size();
	I = sistema_a->termeletricasVtr.size();
	CH = 0;
	vetorint usinas_fim;
	for (size_t r = 0; r < R; r++)
		if ( sistema_a->hidreletricasVtr[r].GetUsinaJusante() == 0)
		{
			usinas_fim.push_back(r);
			CH++;
		}
	cascata.resize(CH);
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
		for (int ch = 0; ch < CH; ch++)
		{
			if ( ult_usina == usinas_fim[ch] )
				cascata[ch].push_back(r);
		}	
	}
	N = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );

	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	n = N * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	na = N * (sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size());
	nd = N * (sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size() + (1 - flag3));		// numero de variáveis duais (e n de variaveis auxiliares/duplicadas)
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= N * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= N * (sistema_a->barrasVtr.size() - 1);
	if ( sistema_a->GetFlagVfol() == false )	// Remover vfol
		n -= sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	x.resize(n);
	x_med.resize(n);
	xa.resize(na);
	lambda.resize(0);
	obj = 0;
	obj_subp.resize(CH + N*R + I + N);
	for (int i = 0; i <  CH + R*N + I + N; i++)
		obj_subp[i] = 0;
	ptr_x_subp.resize(CH + R*N + I + N);
	ptr_x_subp_agg.resize(1);
	ptr_x_subp_agg[0].resize(n + na);

	//AjustarPreconditioners();

	//ptr_x_hat = new vetorfloat[CH + R*N + I + N];		// iniciar com esse tamanho e depois fazer um resize em SetComponents()
	//ptr_x_til = new vetorfloat[CH + R*N + I + N];

	//decomp_str = 0;

	// Criar vetor coluna com o nome das variáveis e tamanho do numero de variáveis em cada período
	ostringstream aux;
	for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
	{
		aux << i + 1;
		stringVariaveis.push_back("pt" + aux.str());
		aux.str("");
	}
	for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
	{
		aux << i + 1;
		stringVariaveis.push_back("u" + aux.str());
		aux.str("");
	}
	if (sistema_a->GetFlagTbinaryModel() == 0)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
		{
			aux << i + 1;
			stringVariaveis.push_back("cp" + aux.str());
			aux.str("");
		}
	}
	else
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
		{
			aux << i + 1;
			stringVariaveis.push_back("up" + aux.str());
			aux.str("");
		}
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
		{
			aux << i + 1;
			stringVariaveis.push_back("ud" + aux.str());
			aux.str("");
		}
	}
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
		{
			aux << i + 1;
			stringVariaveis.push_back("F" + aux.str());
			aux.str("");
		}
	}
	if (sistema_a->GetFlagModeloRede() == 1)
	{
		for (size_t b = 1; b < sistema_a->barrasVtr.size() + 1; b++)
		{
			aux << b + 1;
			stringVariaveis.push_back("teta" + aux.str());
			aux.str("");
		}
		//remover barra de referencia
		stringVariaveis.erase(stringVariaveis.end() - (sistema_a->barrasVtr.size() - sistema_a->GetBarraRef()) - 1);
	}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
	{
		aux << r + 1;
		stringVariaveis.push_back("ph" + aux.str());
		aux.str("");
	}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
	{
		aux << r + 1;
		stringVariaveis.push_back("v" + aux.str());
		aux.str("");
	}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
	{
		aux << r + 1;
		stringVariaveis.push_back("d" + aux.str());
		aux.str("");
	}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
	{
		aux << r + 1;
		stringVariaveis.push_back("s" + aux.str());
		aux.str("");
	}
	if (sistema_a->GetFlagPhmax() == 1)
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		{
			aux << r + 1;
			stringVariaveis.push_back("phmax" + aux.str());
			aux.str("");
		}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			aux << j + 1 << "," << r + 1;
			stringVariaveis.push_back("phg" + aux.str());
			aux.str("");
		}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			aux << j + 1 << "," << r + 1;
			stringVariaveis.push_back("q" + aux.str());
			aux.str("");
		}
	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
		{
			aux << j + 1 << "," << r + 1;
			stringVariaveis.push_back("z" + aux.str());
			aux.str("");
		}
	if (sistema_a->GetFlagModeloRede() > 0)	
	{
		for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
		{
			aux << b + 1;
			stringVariaveis.push_back("def" + aux.str());
			aux.str("");
		}
	}
	else
	{
		stringVariaveis.push_back("def");
	}
	if (sistema_a->GetFlagVfol() == true)
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		{
			aux << r + 1;
			stringVariaveis.push_back("vfol" + aux.str());
			aux.str("");
		}	
}
Spcdec3Results::~Spcdec3Results(void)
{
	delete inFile;
}

//void Spcdec3Results::EscreverArquivo(string nome_arquivo, double tempo_modelo, double tempo_resol, string saida)
//{
//	int flag2 = int (sistema_a->GetFlagVfol());
//	size_t n = x.size();
//	size_t nt = stringVariaveis.size() - flag2*sistema_a->hidreletricasVtr.size();
//	int num_cen = sistema_a->GetNCenarios();
//	//num_cen = 1;		// para escrever resultado de um cenario
//	int T1 = sistema_a->GetTt1();
//	int T2 = sistema_a->GetTt2();
//	//CreateDirectoryA(varString.c_str(), NULL);		// criar pasta
//	inFile = new ofstream( saida + nome_arquivo, ios::out );
//	if ( inFile->is_open() )
//	{
//		if (status != 2)
//		{
//			*inFile << "Solução não ótima! " << "Status final = " << status << char(9) << "Gap = " << mipgap << endl;
//			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
//			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
//			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
//			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagModeloRede() << char(9) << sistema_a->GetFlagVfol() << char(9) << sistema_a->GetFlagPhmax() << char(9) << sistema_a->GetFlagInitAproxCT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
//		}
//		else
//		{
//			*inFile << "Solução ótima! " << "Status final = " << status << char(9) << "Gap = " << mipgap << endl;
//			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
//			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
//			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
//			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagModeloRede() << char(9) << sistema_a->GetFlagVfol() << char(9) << sistema_a->GetFlagPhmax() << char(9) << sistema_a->GetFlagInitAproxCT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
//		}
//		*inFile << endl;
//		// Estágio 1
//		*inFile << setw(10) << left << setprecision(6) << "Estágio 1" << endl;
//		*inFile << setw(10) << left << setprecision(6) << "Periodo";
//		for (int t = 0; t < T1; t++)
//			*inFile << setw(14) << left << char(9) << t + 1;
//		for (size_t i = 0; i < nt; i++)
//		{
//			*inFile << endl;
//			*inFile << setw(10) << left << stringVariaveis[i];
//			for (int t = 0; t < T1; t++)
//				*inFile << char(9) << setw(15) << right << double(x[i + t*nt] > 0 ? double (x[i + t*nt] <= 1e-10 ? 0 : x[i + t*nt]) : x[i + t*nt]);
//		}
//		// Estágio 2
//		for (int n_c = 0; n_c < num_cen; n_c++)
//		{
//			*inFile << endl << endl;
//			*inFile << setw(10) << left << setprecision(6) << "Estágio 2 - Cenário "<< n_c + 1 << endl;
//			*inFile << setw(10) << left << setprecision(6) << "Periodo";
//			for (int t = T1; t < T2; t++)
//				*inFile << setw(14) << left << char(9) << t + 1;
//			for (size_t i = 0; i < nt; i++)
//			{
//				*inFile << endl;
//				*inFile << setw(10) << left << stringVariaveis[i];
//				for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1); t++)
//				{
//					*inFile << char(9) << setw(15) << right << double(x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] > 0 ? double (x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] <= 1e-10 ? 0 : x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]) : x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]);
//				}
//			}
//			if (sistema_a->GetFlagVfol() == 1)
//				for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
//				{
//					*inFile << endl;
//					*inFile << setw(10) << left << stringVariaveis[i + nt];
//					for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1) - 1; t++)
//						*inFile << char(9) << setw(15) << right << " - ";
//					*inFile << char(9) << setw(15) << right << double (x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] > 0 ? double (x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] <= 1e-10 ? 0 : x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]) : x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]);
//				}
//		}
//		// Fluxo nas linhas
//		vetorfloat fluxo;
//		fluxo.resize(sistema_a->linhasVtr.size()*(T1 + num_cen * (T2 - T1)));
//		if (sistema_a->GetFlagModeloRede() == 1)
//		{
//			vetorfloat teta;
//			CMatrizEsparsa M(0);
//			M = MatrizCalcFluxo(int (x.size()) );		//Multiplicar matriz MatrizLimFluxo pelo vetor x
//			fluxo = M.MultiplicarPorVetorDenso(&x);
//		}
//		else	// sistema_a->GetFlagModeloRede() != 1
//		{
//			int delta = 0;
//			int cen = 0;
//			double D = 0;
//			int JJ = 0;
//			for (size_t i = 0; i < R; i++)
//				JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
//			MatrixXd Agh = MatrixXd::Zero(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
//			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
//				for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
//					for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
//						if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
//							Agh(b, r) = 1;
//			MatrixXd Agt = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
//			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
//				for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
//					for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
//						if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
//							Agt(b, i) = 1;
//			MatrixXd Agtu = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
//			if (sistema_a->GetFlagTbinaryModel() == 1)
//			{
//				for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
//					for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
//						for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
//							if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
//								Agtu(b, i) = sistema_a->termeletricasVtr[i].GetPmin();
//			}
//			MatrixXd Adef = MatrixXd::Identity(sistema_a->barrasVtr.size(),sistema_a->barrasVtr.size());
//			Agh = - sistema_a->Beta * Agh;
//			Agt = - sistema_a->Beta * Agt;
//			Agtu = - sistema_a->Beta * Agtu;
//			Adef = - sistema_a->Beta * Adef;
//			VectorXd Ad(sistema_a->barrasVtr.size());
//			VectorXd gg(sistema_a->linhasVtr.size());
//			VectorXd g, Rhs;
//			VectorXd flow(sistema_a->linhasVtr.size());
//			for (int t = 0; t < N; t++)
//			{
//				Ad.resize(sistema_a->barrasVtr.size());
//				if (sistema_a->GetTt1() <= t)
//					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
//				// Determinar o vetor de demandas
//				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
//				{
//					D = 0;
//					for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
//						D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
//					Ad(b) = D;
//				}
//				Rhs = - sistema_a->Beta * ( Ad );
//				// Determinar o vetor de gerações
//				gg = VectorXd::Zero(sistema_a->linhasVtr.size());
//				g.resize(R);
//				for (size_t r = 0; r < R; r++)
//					g(r) = x[r + delta + (3+flag4+flag7)*I];
//				gg = gg + Agh * g;
//				g.resize(I);
//				for (size_t i = 0; i < I; i++)
//					g(i) = x[i + delta];
//				gg = gg + Agt * g;
//				g.resize(I);
//				if (sistema_a->GetFlagTbinaryModel() == 1)
//					for (size_t i = 0; i < I; i++)
//						g(i) = x[i + delta + I];
//				gg = gg + Agtu * g;
//				g.resize(sistema_a->barrasVtr.size());
//				if ( sistema_a->GetFlagModeloRede() == 0 )		// modelo de barra unica, deficit alocado igualmente entre todas barras, o que não é verdade!
//				{
//					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
//						g(b) = x[delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ] / sistema_a->barrasVtr.size();
//				}
//				else
//				{
//					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
//						g(b) = x[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
//				}
//				gg = gg + Adef * g;
//				// Conferir limites das linhas extrapolados para cada periodo
//				flow = - gg + Rhs;		// Beta ( d - g): fluxo nas linhas para o periodo t
//				for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
//					fluxo[l + t*sistema_a->linhasVtr.size()] = flow(l);
//				if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
//					delta += nt + flag2*R;
//				else
//					delta += nt;
//			}
//		}
//		// Escrever fluxo
//		*inFile << endl << endl;
//		if ( sistema_a->GetFlagModeloRede() == 0 )	
//		{
//			*inFile << "Atenção: Modelo de barra única, caso exista déficit ele será alocado igualmente entre as barras! Portanto o resultado pode ser diferente do modelo de rede.";
//			*inFile << endl << endl;
//		}
//		*inFile << setw(10) << left << setprecision(5) << "Fluxo nas linhas (MW)" << endl;
//		*inFile << setw(10) << left << setprecision(5) << "Estágio 1" << endl;
//		for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
//		{
//			*inFile << setprecision(5) << "L" ;
//			*inFile << setprecision(5) << l + 1;
//			*inFile << setprecision(1) << std::fixed << " (" << sistema_a->linhasVtr[l].GetCapacidade() << ")" << setw(8);
//			*inFile << setprecision(5) << std::scientific;
//			for (int t = 0; t < T1; t++)
//			{
//				if ((fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade()) ||
//					(fluxo[l + t*sistema_a->linhasVtr.size()] <= - 0.9999*sistema_a->linhasVtr[l].GetCapacidade()))
//					*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()] << "*";
//				else
//					*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()];
//
//			}
//			*inFile << endl;
//		}
//		for (int n_c = 0; n_c < num_cen; n_c++)
//		{
//			*inFile << setw(10) << left << setprecision(6) << "Estágio 2 - Cenário "<< n_c + 1 << endl;
//			for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
//			{
//				*inFile << setprecision(5) << "L" ;
//				*inFile << setprecision(5) << l + 1;
//				*inFile << setprecision(1) << std::fixed << " (" << sistema_a->linhasVtr[l].GetCapacidade() << ")" << setw(8);
//				*inFile << setprecision(5) << std::scientific;
//				for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1); t++)
//				{
//					if ((fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade()) ||
//						(fluxo[l + t*sistema_a->linhasVtr.size()] <= - 0.9999*sistema_a->linhasVtr[l].GetCapacidade()))
//						*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()] << "*";
//					else
//						*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()];
//				}
//				*inFile << endl;
//			}
//		}
//		inFile->close();
//	}
//	else
//		cout << "Unable to open file";
//}
//void Spcdec3Results::GravarSolHeuristica(double obj_a, CMatrizEsparsa * x_spr, double LB)
//{
//	obj = obj_a;
//	mipgap = (obj_a - LB) / (LB + 1);
//
//	for (size_t i = 0; i < n; i++)
//		x[i] = x_spr->GetElemento(i, 0);
//
//	CIO = 0;		// Custo imediato de operação (custo do primeiro estágio), considerando o custo de operação e partida das termelétricas
//	CTL = 0;		// Custo de operação das termelétricas linear
//	CTQ = 0;		// Custo de operação das termelétricas quadrático
//	CFO = 0;		// Custo futuro de operação ( valor incremental da água, da função de custo futuro)
//	DEF = 0;		// Somatorio dos déficits
//	VFOL = 0;		// Somatorio dos vfolga
//
//	int nt;
//	//n = x.size();
//	int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
//	double deltaT;
//	int delta = 0;
//	int flag1a = 0;		// referente à var. teta
//	if (flag1 == 1)
//		flag1a = 1;
//	int flag1d = 1;		// referente à var. def
//	if (flag1 == 0)
//		flag1d = 0;
//	int cen;
//	for (int t = 0; t < N; t++)
//	{
//		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
//			nt = (n_a / N) + flag2*sistema_a->hidreletricasVtr.size();
//		else
//			nt = (n_a / N);
//		if ((0 <= t) && (sistema_a->GetTt1() > t))
//		{
//			deltaT = sistema_a->GetDeltaT1();
//			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
//			{
//				if (sistema_a->GetFlagTbinaryModel() == 0)
//				{
//					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT;		// Custo de operaçao quadrático
//					if (sistema_a->GetFlagInitAproxCT() > 1)
//					{
//						CTL += x[i + 3*sistema_a->termeletricasVtr.size() + delta];		// Custo de operaçao linear
//						CIO += x[i + 2*sistema_a->termeletricasVtr.size() + delta] + x[i + 3*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
//					}
//					else
//						CIO += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
//				}
//				else
//				{
//					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT;		// Custo de operaçao quadrático
//					if (sistema_a->GetFlagInitAproxCT() > 1)
//					{
//						CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta];		// Custo de operaçao linear
//						CIO += x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
//					}
//					else
//						CIO += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		// Custo de 1 estagio
//				}
//			}
//		}
//		else
//		{
//			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
//			deltaT = sistema_a->GetDeltaT2();
//			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
//			{
//				if (sistema_a->GetFlagTbinaryModel() == 0)
//				{
//					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta])*deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
//					if (sistema_a->GetFlagInitAproxCT() > 1)
//					{
//						CTL += x[i + 3*sistema_a->termeletricasVtr.size() + delta] * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
//						CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] + x[i + 3*sistema_a->termeletricasVtr.size() + delta]) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
//					}
//					else
//						CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) + x[i + 2*sistema_a->termeletricasVtr.size() + delta] ) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
//				}
//				else
//				{
//					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
//					if (sistema_a->GetFlagInitAproxCT() > 1)
//					{
//						CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta]* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
//						CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta])* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
//					}
//					else
//						CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida())* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
//				}
//			}
//		}
//		delta = delta + nt;
//	}
//
//	delta = (n_a / N) - (1 - flag1d) - flag1d*sistema_a->barrasVtr.size(); //(n_a / N) = 34 52
//	for (int t = 0; t < N; t++)
//	{
//		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
//			nt = (n_a / N) + flag2*sistema_a->hidreletricasVtr.size();
//		else
//			nt = (n_a / N);
//		if ((0 <= t) && (sistema_a->GetTt1() > t))
//		{
//			deltaT = sistema_a->GetDeltaT1();
//			for (size_t i = 0; i < (1 - flag1d) + flag1d*sistema_a->barrasVtr.size(); i++)
//				DEF += x[i + delta];		// Def
//		}
//		else
//		{
//			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
//			deltaT = sistema_a->GetDeltaT2();
//			for (size_t i = 0; i < (1 - flag1d) + flag1d*sistema_a->barrasVtr.size(); i++)
//				DEF += x[i + delta];					// Def
//		}
//		delta = delta + nt;
//	}
//
//	if (sistema_a->GetFlagVfol() == 1)
//	{
//		delta = (n_a / N) * sistema_a->GetTt2();
//		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
//		{
//			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
//				VFOL += x[delta++];					// vfol
//			delta += (n_a / N) * (sistema_a->GetTt2() - sistema_a->GetTt1());
//		}
//	}
//	else
//		VFOL = 0;
//}
//void Spcdec3Results::GravarSolExtDE(GRBVar * vars, GRBConstr * constrs, int nvar, int nvar_a, int nconstr)
//{
//	// gravar em x, xa e L os resultados do smart start
//	for  (int i = 0; i < nvar; i++)
//		x[i] = vars[i].get(GRB_DoubleAttr_X);
//	for  (int i = 0; i < nvar_a; i++)
//		xa[i] = vars[i + nvar].get(GRB_DoubleAttr_X);
//	lambda.resize(nconstr);
//	for  (int i = 0; i < nconstr; i++)
//		lambda[i] = constrs[i].get(GRB_DoubleAttr_Pi);
//	// lambda: restrições devem estar na mesma ordem no Ext DE e na heuristica
//	//cout << "Spcdec3Results: conferir se restrições do modelo Ext estão na mesma ordem que na Heuristica!" << endl;
//}
//void Spcdec3Results::GravarSolExtDE(int nvar, int nvar_a, int nconstr)
//{
//	// gravar em x, xa e L os resultados do smart start
//	for  (int i = 0; i < nvar; i++)
//		x[i] = 0;
//	for  (int i = 0; i < nvar_a; i++)
//		xa[i] = 0;
//	lambda.resize(nconstr);
//	for  (int i = 0; i < nconstr; i++)
//		lambda[i] = 0;
//}

void Spcdec3Results::AlocarXmed(CMatrizEsparsa * x_spr)
{
	// Alocar x_med
	for  (int i = 0; i < n; i++)
		x_spr->SubstituirElemento(i, 0, x[i]);
	
	int delta, deltaa;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	delta = 0;
	deltaa = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	//
	for (int t = 0; t < N; t++)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)		//pt
			x_spr->SubstituirElemento(i + delta, 0, (x[i + delta] + xa[i + deltaa]) / 2);
		delta += (3+flag7+flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			x_spr->SubstituirElemento(r + delta, 0, (x[r + delta] + xa[r + deltaa]) / 2);		//ph
			x_spr->SubstituirElemento(r + sistema_a->hidreletricasVtr.size() + delta, 0, (x[r + sistema_a->hidreletricasVtr.size() + delta] + xa[r + sistema_a->hidreletricasVtr.size() + deltaa]) / 2);		//v
			x_spr->SubstituirElemento(r + 2*sistema_a->hidreletricasVtr.size() + delta, 0, (x[r + 2*sistema_a->hidreletricasVtr.size() + delta] + xa[r + 2*sistema_a->hidreletricasVtr.size() + deltaa]) / 2);	//d
			if (sistema_a->GetFlagPhmax() == 1)
				x_spr->SubstituirElemento(r + 4*sistema_a->hidreletricasVtr.size() + delta, 0, (x[r + 4*sistema_a->hidreletricasVtr.size() + delta] + xa[r + 3*sistema_a->hidreletricasVtr.size() + deltaa]) / 2);	//phmax
		}
		delta -= (3+flag7+flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa -= sistema_a->termeletricasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
		else
			delta += nt;
		deltaa += nat;
	}
}
// SpcDec3

// Fazer outra funçao para guardar resultados dos subproblemas, por componente (guardar x no formato de cada subproblema)
// fazer tb outra funçao para calcular o subgradiente de acordo com o formato do x acima
// só colocar os dados nos formatos do x e xa (e x_med) originais no resultado final!!!

// subproblema HA comp = id1;					// id1 é o índice da cascata x = [v d vfol]
// subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó	//x = [ph v d s (phmax) phg q z]
// subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica	//x = [pt u up ud F]
// subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó	//x = [pt teta ph (phmax) def]

void Spcdec3Results::SetComponents( int * comp_inf, bool aggr )
{
	Aggrgtd = aggr;
	// tamanho de ptr_x_til dependo do tipo de modelo
	// tamanho de cada componente de ptr_x_subp poderia ser definido antes, mas o numero de variáveis de cada subproblema só é obtido após criar os objetos subp.
	if ( Aggrgtd )
	{
		ptr_x_til.resize(1);
		ptr_x_til[0].resize(n + na);		// vetor [X Xa]
	}
	else
	{
		vetorfloat2().swap(ptr_x_subp_agg);		// se o modelo for desagregado zerar vetor ptr_x_subp_agg, vetor n é usado
		ptr_x_til.resize(CH + N*R + I + N);
		for (int i = 0; i < CH + N*R + I + N; i++)
		{
			ptr_x_til[i].resize(comp_inf[i]);
			ptr_x_subp[i].resize(comp_inf[i]);		// deve ser definido caso usa-se easy components, pois a soluçao é alocada e n copiada do subp.
		}
	}
}
void Spcdec3Results::GravarSolucao(double obj_a, vetorfloat &x_a, int status_a, int comp)
{
	// Grava o resultado de cada subproblema
	// para o subproblema HA comp = ch;
	// para o subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó
	// para o subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica
	// para o subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó
	status = status_a;
	// wFi = 0 é o termo constante entao wFi = 1 corresponde à comp = 0;
	// conferir se ptr_x_subp e x_a tem o mesmo tamanho
	if (ptr_x_subp[comp].size() == x_a.size())
	{
		for (size_t i = 0; i < ptr_x_subp[comp].size(); i++)
			ptr_x_subp[comp][i] = x_a[i];
	}
	else
		cout << "Error GravarSolucao(): vetores de tamanho diferente!" << endl;
	//ptr_x_subp[comp] = x_a;		// talvez fique mais rápido receber em cada subproblema o endereço do vetorfloat de ptr_x_subp e preenche-lo (clear, resize e atribuir valor)

	obj_subp[comp] = obj_a;
}
void Spcdec3Results::GetSubGradiente(int wFi_a, double * SubG)
{
	for (int i = 0; i < nd; i++)
		SubG[i] = 0;
	
	// O vetor de subgradientes tem a dimensao do vetor lambda, portanto a posicao do subg q sera preenchido depende do tipo de problema
	// Selecionar tipo do subproblema pelo valor componente: (Hidraulico, Hidreletrico, Termeletrico, Demanda)
	// 1 <= wFi (comp) <= GetNrFi();
	// lambda [pt ph v d (phmax || res)] -> tamanho: [I R R R (R || 1)] x T
	int delta = 0;
	int deltaa = 0;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (nd / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	if ( (wFi_a > 0) && (wFi_a <= CH) )	//Hidraulico ( 1 <= comp <= CH )	//x = [v d (vfol)]
	{
		int ch = wFi_a - 1;
		nt = (ptr_x_subp[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
		deltaa = sistema_a->termeletricasVtr.size();
		for (int t = 0; t < N; t++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)
			{
				SubG[cascata[ch].at(r) + sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[ch][r + delta];						//v
				SubG[cascata[ch].at(r) + 2*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[ch][r + cascata[ch].size() + delta];	//d
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*cascata[ch].size();
			else
				delta += nt;
			deltaa += nat;
		}
	}
	else if ( (wFi_a > CH) && (wFi_a <= CH + N * R) )		//Hidreletrico	//x = [ph v d s (phmax) phg q z]
	{
		int r = (wFi_a - CH - 1) % int (sistema_a->hidreletricasVtr.size());
	    int no = (wFi_a - CH - 1) / int (sistema_a->hidreletricasVtr.size());
		
		deltaa = sistema_a->termeletricasVtr.size();
		deltaa += no * nat;

		SubG[r + deltaa] =  ptr_x_subp[wFi_a - 1][0];											//ph
		SubG[r + sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[wFi_a - 1][1];		//v
		SubG[r + 2*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[wFi_a - 1][2];		//d
		if (sistema_a->GetFlagPhmax() == 1)
			SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[wFi_a - 1][4];		//phmax
		else
		{
			double som = 0;
			for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				som += ptr_x_subp[wFi_a - 1][4 + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + j]*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();		// z
			SubG[3*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[wFi_a - 1][0] + som;		//reserva
		}
	}
	else if ( (wFi_a > CH + N * R) && (wFi_a <= CH + N * R + I) )		//Termeletrico	//x = [pt u cp F] ou [pt u up ud F]
	{
		int i = wFi_a - CH - N * R - 1;
		deltaa = 0;
		delta = 0;
		nt = ( ptr_x_subp[wFi_a - 1].size() ) / N;
		if (sistema_a->GetFlagTbinaryModel() == 0)
		{
			for (int t = 0; t < N; t++)
			{
				SubG[i + deltaa] = ptr_x_subp[wFi_a - 1][delta];		//pt
				delta += nt;
				deltaa += nat;
			}
		}
		else
		{
			for (int t = 0; t < N; t++)
			{
				SubG[i + deltaa] = ptr_x_subp[wFi_a - 1][delta] + ptr_x_subp[wFi_a - 1][delta + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//pt + u*pt_min
				delta += nt;
				deltaa += nat;
			}
		}
	}
	else if ( (wFi_a > CH + N * R + I) && (wFi_a <= CH + N * R + I + N) )		//Demanda	//x = [pt teta ph phmax def]
	{
		int no = wFi_a - CH - N * R - I - 1;
		int cen = 0;
		if (no >= sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
			cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		deltaa = no * nat;

		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			SubG[i + deltaa] = - ptr_x_subp[wFi_a - 1][i];		//pt
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - ptr_x_subp[wFi_a - 1][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)];			//ph
			if (sistema_a->GetFlagPhmax() == 1)
				SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[wFi_a - 1][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];	//phmax
		}
	}

	for (int i = 0; i < nd; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);
}
void Spcdec3Results::GetSubGradiente(double * SubG)
{
	int deltaa = 0;
	int delta = 0;
	int cen = 0;
	int nt;
	int nat = (nd / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int t = 0; t < N; t++)
	{
		if (sistema_a->GetFlagTbinaryModel() == 0)
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				SubG[i + deltaa] = - ptr_x_subp[t + CH + N*R + I][i] + ptr_x_subp[i + CH + N*R][t*(3 + flag4)];		// - pta + pt
		else
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				SubG[i + deltaa] = - ptr_x_subp[t + CH + N*R + I][i] + ptr_x_subp[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_subp[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
		deltaa += sistema_a->termeletricasVtr.size();
		double som = 0;
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - ptr_x_subp[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_subp[CH + t*R + r][0];		//ph
			if (sistema_a->GetFlagPhmax() == 1)
				SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()] + ptr_x_subp[CH + t*R + r][4];		//phmax
			else
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					som += ptr_x_subp[CH + t*R + r][4 + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + j]*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();		// z
				som += - ptr_x_subp[CH + t*R + r][0];
			}
		}
		if (sistema_a->GetFlagPhmax() == 0)
			SubG[3*sistema_a->hidreletricasVtr.size() + deltaa] = som;		//reserva: sum_r(ph_r - sum_j(phg_max*z_j))
		for (size_t ch = 0; ch < cascata.size(); ch++)
		{
			nt = (ptr_x_subp[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
			if ( t >= sistema_a->GetTt2() )
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				delta = t*nt + flag2*cen*cascata[ch].size();
			}	
			else
				delta = t*nt;

			for (size_t r = 0; r < cascata[ch].size(); r++)
			{
				SubG[cascata[ch].at(r) + sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[ch][r + delta] + ptr_x_subp[CH + t*R + cascata[ch].at(r)][1];		//v
				SubG[cascata[ch].at(r) + 2*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[ch][r + cascata[ch].size() + delta] + ptr_x_subp[CH + t*R + cascata[ch].at(r)][2];		//d
			}
		}
		deltaa -= sistema_a->termeletricasVtr.size();
		deltaa += nat;
	}

	for (int i = 0; i < nd; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);
}
double Spcdec3Results::GetFobj()
{
	//obj = 0;
	//inFile = new ofstream( "obj_subprb.txt", ios::out );
	//if ( inFile->is_open() )
	//{
	//	for (int sp = 0; sp < CH + R*N + I + N; sp++)
	//	{
	//		obj += obj_subp[sp];
	//		*inFile << obj_subp[sp] << endl;
	//	}
	//	inFile->close();
	//}
	//else
	//	cout << "Unable to open file";

	obj = 0;
	for (int sp = 0; sp < CH + R*N + I + N; sp++)
		obj += obj_subp[sp];
	return obj;
}
vetorfloat Spcdec3Results::GetCompX(int comp)
{
	if ( Aggrgtd )
	{
		AlocarPtrXHat();
		return ptr_x_subp_agg[0];
		
		//return ptr_x_hat[0];		// vetor [X Xa], assim o ptr_x_til vai ter a mesma estrutura do ptr_x_hat
	}
	else
	{
		return ptr_x_subp[comp];
		//ptr_x_hat[comp] = ptr_x_subp[comp];
		//return ptr_x_hat[comp];
	}
}	
vetorfloat Spcdec3Results::GetX_med()
{
	// Alocar x_med a partir de x e xa
	AlocarX_med();

	return x_med;
}
vetorfloat Spcdec3Results::GetX_med(bool til)
{
	// Alocar x_med a partir de ptr_x_hat ou ptr_x_til (estrutura depende se o modelo é agregado ou desagregado)
	// ptr_x_hat ou ptr_x_til tem a mesma estrutura!!
	// A cada chamada do SetGiName o ptr_x_hat aqui nessa classe é atualizado, ou seja, ptr_x_hat é atualizado right after Fi() and all the GetGi()

	// hat_or_til indica a partir de qual solução é calculado o x_med
	// 0 a partir da soluçao da RL, ptr_x_hat
	// 1 a partir da soluçao convexifica, ptr_x_til
	vetorfloat * ptr_x_a;
	if ( til )
		ptr_x_a = &ptr_x_til[0];
	else
		if ( Aggrgtd )
			ptr_x_a = &ptr_x_subp_agg[0];
		else
			ptr_x_a = &ptr_x_subp[0];

	// Para as variáveis duplicadas x_med = (x + xa) / 2
	//x = [pt u cp F teta ph v d s (phmax) phg q z def] ou
	//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model
	//xa = [pta pha va da (phmaxa)]
	int delta, deltaa;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	if ( Aggrgtd )
	{
		for (int i = 0; i < n; i++)		// for para x_med[i] = ptr_x_hat[0][i];
			x_med[i] = ptr_x_a[0][i];
		// Atualizar váriáveis duplicadas, x_med[i] = (ptr_x_hat[0][i] + ptr_x_hat[0][n + i]) / 2;
		delta = 0;
		deltaa = 0;
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)		//pt
					x_med[i + delta] = (ptr_x_a[0][i + delta] + ptr_x_a[0][n + i + deltaa]) / 2;
			}
			else
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)		//pt
					x_med[i + delta] = (max(ptr_x_a[0][n + i + deltaa] - sistema_a->termeletricasVtr[i].GetPmin(), 0.0) + ptr_x_a[0][i + delta]) / 2;
					//x_med[i + delta] = max(double(((ptr_x_a[0][i + delta] + ptr_x_a[0][i + 1 + delta]*sistema_a->termeletricasVtr[i].GetPmin()) + ptr_x_a[0][n + i + deltaa]) / 2 - ptr_x_a[0][i + 1 + delta]*sistema_a->termeletricasVtr[i].GetPmin()), double(0));
				// pt é a variável original e nao pta, então o valor final é: (pt + pt_min*u + pta) / 2 - pt_min*u, porém pode dar negativo, por isso o max().
			}
			delta += (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
			deltaa += sistema_a->termeletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				x_med[r + delta] = (ptr_x_a[0][r + delta] + ptr_x_a[0][n + r + deltaa]) / 2;		//ph
				x_med[r + sistema_a->hidreletricasVtr.size() + delta] = (ptr_x_a[0][r + sistema_a->hidreletricasVtr.size() + delta] + ptr_x_a[0][n + r + sistema_a->hidreletricasVtr.size() + deltaa]) / 2;				//v
				x_med[r + 2*sistema_a->hidreletricasVtr.size() + delta] = (ptr_x_a[0][r + 2*sistema_a->hidreletricasVtr.size() + delta] + ptr_x_a[0][n + r + 2*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//d
				if (sistema_a->GetFlagPhmax() == 1)
					x_med[r + 4*sistema_a->hidreletricasVtr.size() + delta] = (ptr_x_a[0][r + 4*sistema_a->hidreletricasVtr.size() + delta] + ptr_x_a[0][n + r + 3*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//phmax
			}
			delta -= (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
			deltaa -= sistema_a->termeletricasVtr.size();
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
			deltaa += nat;
		}
	}
	else
	{
		// subproblema HA comp = ch;					//x = [v d vfol]
		// subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó	//x = [ph v d s (phmax) phg q z]
		// subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica	//x = [pt u cp F] ou x = [pt u up ud F]
		// subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó	//x = [pt teta ph (phmax) def]
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int jj, cen;
		delta = 0;
		deltaa = 0;
		//double valor = 0;
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)		// flag7 = 0
			{
				for (int i = 0; i < I; i++)
				{
					x_med[i + delta] = (ptr_x_a[t + CH + N*R + I][i] + ptr_x_a[i + CH + N*R][t*(3 + flag4)]) / 2;		//pt
					x_med[i + I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4) + 1];									//u		
					x_med[i + 2*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4) + 2];									//cp
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)
					for (int i = 0; i < I; i++)
						x_med[i + 3*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4) + 3];								//F
			}
			else
			{
				for (int i = 0; i < I; i++)
				{
					// o valor de pt aqui deve estar entre 0 e ptmax-ptmin! Pois na heuristica a variavel é modelada assim
					//x_med[i + delta] = max(double((ptr_x_a[t + CH + N*R + I][i] + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin()) / 2 - ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin()), double(0));		//pt
					x_med[i + delta] = (max(ptr_x_a[t + CH + N*R + I][i] - sistema_a->termeletricasVtr[i].GetPmin(), 0.0) + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7)]) / 2;		//pt
					x_med[i + I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1];									//u		
					x_med[i + 2*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 2];									//up
					x_med[i + 3*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 3];									//ud
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)
					for (int i = 0; i < I; i++)
						x_med[i + (3 + flag7)*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 3 + flag7];								//F
			}
			jj = I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R;
			for (int r = 0; r < R; r++)
			{
				x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + delta] = (ptr_x_a[t + CH + N*R + I][r + I + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_a[CH + t*R + r][0]) / 2;				//ph
				x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 3*R + delta] = ptr_x_a[CH + t*R + r][3];		//s
				if (sistema_a->GetFlagPhmax() == 1)
					x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + delta] = (ptr_x_a[t + CH + N*R + I][r + I + flag1a*(sistema_a->barrasVtr.size() - 1) + R] + ptr_x_a[CH + t*R + r][4]) / 2;	//phmax
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					x_med[j + jj + delta] = ptr_x_a[CH + t*R + r][j + (4 + flag3)];															//phg
					x_med[j + jj + JJ + delta] = ptr_x_a[CH + t*R + r][j + sistema_a->hidreletricasVtr[r].GetNGrupos() + (4 + flag3)];		//q
					x_med[j + jj + 2*JJ + delta] = ptr_x_a[CH + t*R + r][j + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + (4 + flag3)];	//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if (sistema_a->GetFlagModeloRede() > 0)
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					x_med[b + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta] = ptr_x_a[t + CH + N*R + I][b + I + flag1a*(sistema_a->barrasVtr.size() - 1) + (1+flag3)*R];	//def
			}
			else
				x_med[I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta] = ptr_x_a[t + CH + N*R + I][I + (1+flag3)*R];	//def
			if (sistema_a->GetFlagModeloRede() == 1)
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
					x_med[b + I*(3 + flag4 + flag7) + delta] = ptr_x_a[t + CH + N*R + I][b + I];	//teta
			}
			for (size_t ch = 0; ch < cascata.size(); ch++)
			{
				nat = (ptr_x_a[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
				if ( t >= sistema_a->GetTt2() )
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					deltaa = t*nat + flag2*cen*cascata[ch].size();
				}	
				else
					deltaa = t*nat;

				for (size_t r = 0; r < cascata[ch].size(); r++)
				{
					//valor = (ptr_x_a[ch][r + deltaa] + ptr_x_a[CH + t*R + cascata[ch].at(r)][1]) / 2;
					//x_med[cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + R + delta] = double ( valor < sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin() ? sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmin() : double (valor > sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax() ? sistema_a->hidreletricasVtr[cascata[ch].at(r)].GetVmax() : valor) );							//v
					//x_med[cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + R + delta] = double ( valor );							//v
					x_med[cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + R + delta] = (ptr_x_a[ch][r + deltaa] + ptr_x_a[CH + t*R + cascata[ch].at(r)][1]) / 2;							//v

					//if (cascata[ch].at(r) == 5)
					//	cout << t << "= " << ptr_x_a[ch][r + deltaa] << " : " << ptr_x_a[CH + t*R + cascata[ch].at(r)][1] << " : " << (ptr_x_a[ch][r + deltaa] + ptr_x_a[CH + t*R + cascata[ch].at(r)][1]) / 2 << endl;
					x_med[cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R + delta] = (ptr_x_a[ch][r + cascata[ch].size() + deltaa] + ptr_x_a[CH + t*R + cascata[ch].at(r)][2]) / 2;	//d
				}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
		}
		if (sistema_a->GetFlagVfol() == true)
		{
			for (size_t ch = 0; ch < cascata.size(); ch++)
			{
				delta = nt * (sistema_a->GetTt2() - 1);
				nat = (ptr_x_a[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
				deltaa = nat * (sistema_a->GetTt2() - 1);
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)
						x_med[cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + flag1d*sistema_a->barrasVtr.size() + (1 - flag1d) + delta] = ptr_x_a[ch][r + 2*cascata[ch].size() + deltaa];		//vfol
					delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
					deltaa += nat * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[ch].size();
				}
			}
		}
	}
	//delete[] ptr_x_a;
	return x_med;
}
void Spcdec3Results::AdicionarLBXtil(CMatrizEsparsa * x_spr, double mult)
{
	// Adicionar solução setada como LB no Bundle na combinação da solução convexificada

	// Para as variáveis duplicadas x = xa += x_spr * mult
	//x = [pt u cp F teta ph v d s (phmax) phg q z def] ou
	//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model
	//xa = [pta pha va da (phmaxa)]
	int delta, deltaa;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	if ( Aggrgtd )
	{
		if (mult >= 0.99999999)
			for (int i = 0; i < n; i++)
				ptr_x_til[0][i] = x_spr->GetElemento(i, 0) * mult;
		else
			for (int i = 0; i < n; i++)
				ptr_x_til[0][i] += x_spr->GetElemento(i, 0) * mult;
	}
	else
	{
		// subproblema HA comp = ch;					//x = [v d vfol]
		// subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó	//x = [ph v d s (phmax) phg q z]
		// subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica	//x = [pt u cp F] ou x = [pt u up ud F]
		// subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó	//x = [pt teta ph (phmax) def]
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int jj, cen;
		delta = 0;
		deltaa = 0;
		//double valor = 0;
		if (mult >= 0.99999999)
		{
			for (size_t i = 0; i < ptr_x_til.size(); i++)
				for (size_t ii = 0; ii < ptr_x_til[i].size(); ii++)
					ptr_x_til[i][ii] = 0;
		}
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)		// flag7 = 0
			{
				for (int i = 0; i < I; i++)
				{
					ptr_x_til[i + CH + N*R][t*(3 + flag4)] += x_spr->GetElemento(i + delta, 0) * mult;				//pt
					ptr_x_til[t + CH + N*R + I][i] += x_spr->GetElemento(i + delta, 0) * mult;						//pta
					ptr_x_til[i + CH + N*R][t*(3 + flag4) + 1] += x_spr->GetElemento(i + I + delta, 0) * mult;		//u	
					ptr_x_til[i + CH + N*R][t*(3 + flag4) + 2] += x_spr->GetElemento(i + 2*I + delta, 0) * mult;	//cp	
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)
					for (int i = 0; i < I; i++)
						 ptr_x_til[i + CH + N*R][t*(3 + flag4) + 3] += x_spr->GetElemento(i + 3*I + delta, 0) * mult;	//F
			}
			else
			{
				for (int i = 0; i < I; i++)
				{
					// o valor de pt aqui deve estar entre 0 e ptmax-ptmin! Pois na heuristica a variavel é modelada assim
					//x_med[i + delta] = max(double((ptr_x_a[t + CH + N*R + I][i] + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin()) / 2 - ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin()), double(0));		//pt
					ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7)] += x_spr->GetElemento(i + delta, 0) * mult;				//pt
					ptr_x_til[t + CH + N*R + I][i] += (x_spr->GetElemento(i + delta, 0) + x_spr->GetElemento(i + I + delta, 0)*sistema_a->termeletricasVtr[i].GetPmin())* mult;		//pta
					ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7) + 1] += x_spr->GetElemento(i + I + delta, 0) * mult;		//u	
					ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7) + 2] += x_spr->GetElemento(i + 2*I + delta, 0) * mult;	//up	
					ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7) + 3] += x_spr->GetElemento(i + 3*I + delta, 0) * mult;	//ud	
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)
					for (int i = 0; i < I; i++)
						ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7) + 3 + flag7] += x_spr->GetElemento(i + (3 + flag7)*I + delta, 0) * mult;	//F
			}
			jj = I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R;
			for (int r = 0; r < R; r++)
			{
				ptr_x_til[t + CH + N*R + I][r + I + flag1a*(sistema_a->barrasVtr.size() - 1)] += x_spr->GetElemento(r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + delta, 0) * mult;				//pha
				ptr_x_til[CH + t*R + r][0] += x_spr->GetElemento(r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + delta, 0) * mult;					//ph
				ptr_x_til[CH + t*R + r][3] += x_spr->GetElemento(r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 3*R + delta, 0) * mult;				//s
				if (sistema_a->GetFlagPhmax() == 1)
				{
					ptr_x_til[CH + t*R + r][4] += x_spr->GetElemento(r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + delta, 0) * mult;				//phmax
					ptr_x_til[t + CH + N*R + I][r + I + flag1a*(sistema_a->barrasVtr.size() - 1) + R] += x_spr->GetElemento(r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + delta, 0) * mult;				//phmaxa
				}
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					ptr_x_til[CH + t*R + r][j + (4 + flag3)] += x_spr->GetElemento(j + jj + delta, 0) * mult;				//phg
					ptr_x_til[CH + t*R + r][j + sistema_a->hidreletricasVtr[r].GetNGrupos() + (4 + flag3)] += x_spr->GetElemento(j + jj + JJ + delta, 0) * mult;				//q
					ptr_x_til[CH + t*R + r][j + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + (4 + flag3)] += x_spr->GetElemento(j + jj + 2*JJ + delta, 0) * mult;			//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if (sistema_a->GetFlagModeloRede() > 0)
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					ptr_x_til[t + CH + N*R + I][b + I + flag1a*(sistema_a->barrasVtr.size() - 1) + (1+flag3)*R] += x_spr->GetElemento(b + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta, 0) * mult;			//def
			}
			else
				ptr_x_til[t + CH + N*R + I][I + (1+flag3)*R] += x_spr->GetElemento(I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + delta, 0) * mult;			//def
			if (sistema_a->GetFlagModeloRede() == 1)
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
					ptr_x_til[t + CH + N*R + I][b + I] += x_spr->GetElemento(b + I*(3 + flag4 + flag7) + delta, 0) * mult;			//teta
			}
			for (size_t ch = 0; ch < cascata.size(); ch++)
			{
				nat = (ptr_x_til[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
				if ( t >= sistema_a->GetTt2() )
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					deltaa = t*nat + flag2*cen*cascata[ch].size();
				}	
				else
					deltaa = t*nat;

				for (size_t r = 0; r < cascata[ch].size(); r++)
				{
					ptr_x_til[ch][r + deltaa] += x_spr->GetElemento(cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + R + delta, 0) * mult;							//va
					ptr_x_til[CH + t*R + cascata[ch].at(r)][1] += x_spr->GetElemento(cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + R + delta, 0) * mult;			//v
					ptr_x_til[ch][r + cascata[ch].size() + deltaa] += x_spr->GetElemento(cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R + delta, 0) * mult;		//da
					ptr_x_til[CH + t*R + cascata[ch].at(r)][2] += x_spr->GetElemento(cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R + delta, 0) * mult;			//d
				}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
		}
		if (sistema_a->GetFlagVfol() == true)
		{
			for (size_t ch = 0; ch < cascata.size(); ch++)
			{
				delta = nt * (sistema_a->GetTt2() - 1);
				nat = (ptr_x_til[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
				deltaa = nat * (sistema_a->GetTt2() - 1);
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)
						ptr_x_til[ch][r + 2*cascata[ch].size() + deltaa] += x_spr->GetElemento(cascata[ch].at(r) + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*R + 3*JJ + flag1d*sistema_a->barrasVtr.size() + (1 - flag1d) + delta, 0) * mult;			//vfol
					delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
					deltaa += nat * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[ch].size();
				}
			}
		}
	}
}

void Spcdec3Results::AlocarX_med()
{
	// Alocar x_med a partir de x e xa
	AlocarXeXa();
	// Para as variáveis duplicadas x_med = (x + xa) / 2
	//x = [pt u cp F teta ph v d s phmax phg q z def]
	//xa = [pta pha va da phmaxa]
	
	// Variáveis duplicadas dependem do tipo de decomposiçao!!!!! (implementar isso)

	x_med = x;
	int delta, deltaa;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	delta = 0;
	deltaa = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	//
	for (int t = 0; t < N; t++)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)		//pt
			x_med[i + delta] = (x[i + delta] + xa[i + deltaa]) / 2;
		delta += (3+flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			x_med[r + delta] = (x[r + delta] + xa[r + deltaa]) / 2;		//ph
			x_med[r + sistema_a->hidreletricasVtr.size() + delta] = (x[r + sistema_a->hidreletricasVtr.size() + delta] + xa[r + sistema_a->hidreletricasVtr.size() + deltaa]) / 2;				//v
			x_med[r + 2*sistema_a->hidreletricasVtr.size() + delta] = (x[r + 2*sistema_a->hidreletricasVtr.size() + delta] + xa[r + 2*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//d
			if (sistema_a->GetFlagPhmax() == 1)
				x_med[r + 4*sistema_a->hidreletricasVtr.size() + delta] = (x[r + 4*sistema_a->hidreletricasVtr.size() + delta] + xa[r + 3*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//phmax
		}
		delta -= (3+flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa -= sistema_a->termeletricasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
		else
			delta += nt;
		deltaa += nat;
	}
}
void Spcdec3Results::AlocarPtrXHat()
{
	// Alocar ptr_x_hat a partir de ptr_x_subp para o caso agregado
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();

	int delta, deltaa, delta_subP;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int nt_subP;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int comp = 0; comp < CH; comp++)	//Hidraulico
	{
		// comp é o índice da cascata
		int id1 = comp;
		deltaa = sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size();
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() - flag2*sistema_a->GetNCenarios()*cascata[id1].size()) / N;
		// Alocar resultados do subproblema
		// xa e x
		for (int t = 0; t < N; t++)
		{
			for (size_t r = 0; r < cascata[id1].size(); r++)
			{
				ptr_x_subp_agg[0][n + cascata[id1].at(r) + deltaa] = ptr_x_subp[comp][r + delta_subP];		//v
				ptr_x_subp_agg[0][n + cascata[id1].at(r) + sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + cascata[id1].size() + delta_subP];		//d
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta_subP += nt_subP + flag2*cascata[id1].size() ;
			else
				delta_subP += nt_subP;
			deltaa += nat;
		}
		if (flag2 == 1)
		{
			delta = nt * sistema_a->GetTt2();
			delta_subP = nt_subP * sistema_a->GetTt2();
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < cascata[id1].size(); r++)
					ptr_x_subp_agg[0][delta + cascata[id1].at(r)] = ptr_x_subp[comp][delta_subP + r];		//vfol
				delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
				delta_subP += nt_subP * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[id1].size();
			}
		}
	}
	for (int comp = CH; comp < CH + N * R; comp++)		//Hidreletrico
	{
		// id1 é o índice da usina e id2 é o índice do nó
		int id1 = (comp - CH) % int (sistema_a->hidreletricasVtr.size());
		int id2 = (comp - CH) / int (sistema_a->hidreletricasVtr.size());

		delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa = 0;
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() / N);
		int Jr = 0;
		for (size_t r = 0; r < id1; r++)
			Jr += sistema_a->hidreletricasVtr[r].GetNGrupos();

		delta += id2 * nt;
		if (id2 >= sistema_a->GetTt1())
			delta += ((id2 - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1())) * flag2*sistema_a->hidreletricasVtr.size();
		// Alocar resultados do subproblema
		// xa e x
		ptr_x_subp_agg[0][id1 + delta] = ptr_x_subp[comp][0];		//ph
		ptr_x_subp_agg[0][id1 + 1*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][1];		//v
		ptr_x_subp_agg[0][id1 + 2*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][2];		//d
		ptr_x_subp_agg[0][id1 + 3*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][3];		//s
		if (sistema_a->GetFlagPhmax() == 1)
			ptr_x_subp_agg[0][id1 + 4*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][4];		//phmax
		for (int j = 0; j < sistema_a->hidreletricasVtr[id1].GetNGrupos(); j++)
		{
			ptr_x_subp_agg[0][j + (4+flag3)*sistema_a->hidreletricasVtr.size() + Jr + delta] = ptr_x_subp[comp][j + (4+flag3)];					//phg
			ptr_x_subp_agg[0][j + (4+flag3)*sistema_a->hidreletricasVtr.size() + JJ + Jr + delta] = ptr_x_subp[comp][j + sistema_a->hidreletricasVtr[id1].GetNGrupos() + (4+flag3)];		//q
			ptr_x_subp_agg[0][j + (4+flag3)*sistema_a->hidreletricasVtr.size() + 2*JJ + Jr + delta] = ptr_x_subp[comp][j + 2*sistema_a->hidreletricasVtr[id1].GetNGrupos() + (4+flag3)];	//z
		}
	}
	for (int comp = CH + N * R; comp < CH + N * R + I; comp++)		//Termeletrico
	{
		// id1 é o índice da termeletrica e id2 é nada
		int id1 = comp - CH - N * R;
		delta = 0;
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() / N);
		// Alocar resultados do subproblema
		// xa e x
		for (int t = 0; t < N; t++)
		{
			ptr_x_subp_agg[0][id1 + delta] = ptr_x_subp[comp][delta_subP];														//pt
			ptr_x_subp_agg[0][id1 + 1*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][1 + delta_subP];			//u
			ptr_x_subp_agg[0][id1 + 2*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][2 + delta_subP];			//cp ou up
			if (flag7 == 1)
				ptr_x_subp_agg[0][id1 + 3*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][3 + delta_subP];			//ud
			if (sistema_a->GetFlagInitAproxCT() > 1)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					ptr_x_subp_agg[0][id1 + (3 + flag7)*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][3 + flag7 + delta_subP];	//F
			}

			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
			delta_subP += nt_subP;
		}
	}
	for (int comp = CH + N * R + I; comp < CH + N * R + I + N; comp++)		//Demanda
	{
		// id1 é o índice do nó
		int id1 = comp - CH - N * R - I;
		//x_a = [pt teta ph phmax def]
		delta = (3 + flag4)*sistema_a->termeletricasVtr.size();
		delta += id1 * nt;
		nt_subP = (ptr_x_subp[comp].size() / N);
		if (id1 >= sistema_a->GetTt1())
			delta += ((id1 - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1())) * flag2*sistema_a->hidreletricasVtr.size();
		deltaa = id1 * nat;
		// Alocar resultados do subproblema
		// xa e x
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			ptr_x_subp_agg[0][n + i + deltaa] = ptr_x_subp[comp][i];		//pt
		if (flag1 == 1)		// considerar rede
		{
			for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
				ptr_x_subp_agg[0][b + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size()];	//teta
			for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				ptr_x_subp_agg[0][b + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (1+flag3)*sistema_a->hidreletricasVtr.size()];	//def
		}
		else
			ptr_x_subp_agg[0][flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + (1+flag3)*sistema_a->hidreletricasVtr.size()];	//def
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			ptr_x_subp_agg[0][n + r + sistema_a->termeletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)];		//ph
			if (sistema_a->GetFlagPhmax() == 1)
				ptr_x_subp_agg[0][n + r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];		//phmax
		}
	}
}
void Spcdec3Results::AlocarXeXa()
{
	// Alocar x e xa a partir de ptr_x_subp
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();

	int delta, deltaa, delta_subP;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int nt_subP;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int comp = 0; comp < CH; comp++)	//Hidraulico
	{
		// comp é o índice da cascata
		int id1 = comp;
		deltaa = sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size();
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() - flag2*sistema_a->GetNCenarios()*cascata[id1].size()) / N;
		// Alocar resultados do subproblema
		// xa e x
		for (int t = 0; t < N; t++)
		{
			for (size_t r = 0; r < cascata[id1].size(); r++)
			{
				xa[cascata[id1].at(r) + deltaa] = ptr_x_subp[comp][r + delta_subP];		//v
				xa[cascata[id1].at(r) + sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + cascata[id1].size() + delta_subP];		//d
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta_subP += nt_subP + flag2*cascata[id1].size() ;
			else
				delta_subP += nt_subP;
			deltaa += nat;
		}
		if (flag2 == 1)
		{
			delta = nt * sistema_a->GetTt2();
			delta_subP = nt_subP * sistema_a->GetTt2();
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < cascata[id1].size(); r++)
					x[delta + cascata[id1].at(r)] = ptr_x_subp[comp][delta_subP + r];		//vfol
				delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
				delta_subP += nt_subP * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[id1].size();
			}
		}
	}
	for (int comp = CH; comp < CH + N * R; comp++)		//Hidreletrico
	{
		// id1 é o índice da usina e id2 é o índice do nó
		int id1 = (comp - CH) % int (sistema_a->hidreletricasVtr.size());
		int id2 = (comp - CH) / int (sistema_a->hidreletricasVtr.size());

		delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa = 0;
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() / N);
		int Jr = 0;
		for (size_t r = 0; r < id1; r++)
			Jr += sistema_a->hidreletricasVtr[r].GetNGrupos();

		delta += id2 * nt;
		if (id2 >= sistema_a->GetTt1())
			delta += ((id2 - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1())) * flag2*sistema_a->hidreletricasVtr.size();
		// Alocar resultados do subproblema
		// xa e x
		x[id1 + delta] = ptr_x_subp[comp][0];		//ph
		x[id1 + 1*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][1];		//v
		x[id1 + 2*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][2];		//d
		x[id1 + 3*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][3];		//s
		if (sistema_a->GetFlagPhmax() == 1)
			x[id1 + 4*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][4];		//phmax
		for (int j = 0; j < sistema_a->hidreletricasVtr[id1].GetNGrupos(); j++)
		{
			x[j + (4+flag3)*sistema_a->hidreletricasVtr.size() + Jr + delta] = ptr_x_subp[comp][j + 5];					//phg
			x[j + (4+flag3)*sistema_a->hidreletricasVtr.size() + JJ + Jr + delta] = ptr_x_subp[comp][j + sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];		//q
			x[j + (4+flag3)*sistema_a->hidreletricasVtr.size() + 2*JJ + Jr + delta] = ptr_x_subp[comp][j + 2*sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];	//z
		}
	}
	for (int comp = CH + N * R; comp < CH + N * R + I; comp++)		//Termeletrico
	{
		// id1 é o índice da termeletrica e id2 é nada
		int id1 = comp - CH - N * R;
		delta = 0;
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() / N);
		// Alocar resultados do subproblema
		// xa e x
		for (int t = 0; t < N; t++)
		{
			x[id1 + delta] = ptr_x_subp[comp][delta_subP];														//pt
			x[id1 + 1*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][1 + delta_subP];			//u
			x[id1 + 2*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][2 + delta_subP];			//cp ou up
			if (flag7 == 1)
				x[id1 + 3*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][3 + delta_subP];			//ud
			if (sistema_a->GetFlagInitAproxCT() > 1)
				x[id1 + (3 + flag7)*sistema_a->termeletricasVtr.size() + delta] = ptr_x_subp[comp][3 + flag7 + delta_subP];		//F

			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
			delta_subP += nt_subP;
		}
	}
	for (int comp = CH + N * R + I; comp < CH + N * R + I + N; comp++)		//Demanda
	{
		// id1 é o índice do nó
		int id1 = comp - CH - N * R - I;
		//x_a = [pt teta ph phmax def]
		delta = (3 + flag4)*sistema_a->termeletricasVtr.size();
		delta += id1 * nt;
		nt_subP = (ptr_x_subp[comp].size() / N);
		if (id1 >= sistema_a->GetTt1())
			delta += ((id1 - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1())) * flag2*sistema_a->hidreletricasVtr.size();
		deltaa = id1 * nat;
		// Alocar resultados do subproblema
		// xa e x
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			xa[i + deltaa] = ptr_x_subp[comp][i];		//pt
		if (flag1 == 1)		// considerar rede
		{
			for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
				x[b + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size()];	//teta
			for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				x[b + flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		}
		else
			x[flag1a*(sistema_a->barrasVtr.size() - 1) + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			xa[r + sistema_a->termeletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)];		//ph
			if (sistema_a->GetFlagPhmax() == 1)
				xa[r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];		//phmax
		}
	}
}

void Spcdec3Results::ExportarX(string nome_arquivo)
{
	inFile = new ofstream( nome_arquivo, ios::out );
	if ( inFile->is_open() )
	{
		for (size_t i = 0; i < x.size(); i++)
			*inFile << std::scientific << setprecision(9) << char(9) << setw(15) << right << x[i] << endl;
		inFile->close();
	}
	else
		cout << "Unable to open file";
	//delete inFile;
	inFile = NULL;
}
void Spcdec3Results::ExportarXA(string nome_arquivo)
{
	inFile = new ofstream( nome_arquivo, ios::out );
	if ( inFile->is_open() )
	{
		for (size_t i = 0; i < xa.size(); i++)
			*inFile << std::scientific << setprecision(6) << char(9) << setw(15) << right << xa[i] << endl;
		inFile->close();
	}
	else
		cout << "Unable to open file";
	//delete inFile;
	inFile = NULL;
}
void Spcdec3Results::ExportarXmed(string nome_arquivo)
{
	inFile = new ofstream( nome_arquivo, ios::out );
	if ( inFile->is_open() )
	{
		for (size_t i = 0; i < x_med.size(); i++)
			*inFile << std::scientific << setprecision(9) << char(9) << setw(15) << right << x_med[i] << endl;
		inFile->close();
	}
	else
		cout << "Unable to open file";
	//delete inFile;
	inFile = NULL;
}

void Spcdec3Results::ZerarSolucoes()
{
	x.clear();
	x_med.clear();
	xa.clear();
	obj_subp.clear();
	ptr_x_subp.clear();
	ptr_x_subp_agg.clear();
	ptr_x_til.clear();

	x.resize(n);
	x_med.resize(n);
	xa.resize(na);
	obj = 0;
	obj_subp.resize(CH + N*R + I + N);
	for (int i = 0; i <  CH + R*N + I + N; i++)
		obj_subp[i] = 0;
	ptr_x_subp.resize(CH + R*N + I + N);
	ptr_x_subp_agg.resize(1);
	ptr_x_subp_agg[0].resize(n + na);
}
void Spcdec3Results::CalcularNormaSubgXtil(vetorfloat &norma1, vetorfloat &norma2, vetorfloat &normaInf)
{
	// Um valor de norma para cada tipo de restrição relaxada, norma 1, 2 e Inf;
	int deltaa = 0;
	int delta = 0;
	int cen = 0;
	int nt;
	int nat = (nd / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	double dif = 0;
	// zerar elementos do vetor norma Inf
	for (int i = 0; i < 5; i++)
		normaInf[i] = 0;
	for (int t = 0; t < N; t++)
	{
		if (sistema_a->GetFlagTbinaryModel() == 0)
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				dif = - ptr_x_til[t + CH + N*R + I][i] + ptr_x_til[i + CH + N*R][t*(3 + flag4)];		// - pta + pt
				norma1[0] += abs(dif);
				norma2[0] += pow(dif, 2);
				if (normaInf[0] < abs(dif))
					normaInf[0] = abs(dif);
			}
		else
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				dif = - ptr_x_til[t + CH + N*R + I][i] + ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
				norma1[0] += abs(dif);
				norma2[0] += pow(dif, 2);
				if (normaInf[0] < abs(dif))
					normaInf[0] = abs(dif);
			}
		deltaa += sistema_a->termeletricasVtr.size();
		double som = 0;
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			dif = - ptr_x_til[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_til[CH + t*R + r][0];		//ph
			norma1[1] += abs(dif);
			norma2[1] += pow(dif, 2);
			if (normaInf[1] < abs(dif))
				normaInf[1] = abs(dif);
			if (sistema_a->GetFlagPhmax() == 1)
			{
				dif = - ptr_x_til[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()] + ptr_x_til[CH + t*R + r][4];		//phmax
				norma1[4] += abs(dif);
				norma2[4] += pow(dif, 2);
				if (normaInf[4] < abs(dif))
					normaInf[4] = abs(dif);
			}
			else
			{
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					som += ptr_x_til[CH + t*R + r][4 + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + j]*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();		// z
				som += - ptr_x_til[CH + t*R + r][0];
			}
		}
		if (sistema_a->GetFlagPhmax() == 0)
		{
			dif = min(0.0, som - sistema_a->GetReserva(t));		//reserva: sum_r(sum_j(phg_max*z_j) - ph_r) - Res
			// each positive entry of the vector Ax - b is set to zero
			norma1[4] += abs(dif);
			norma2[4] += pow(dif, 2);
			if (normaInf[4] < abs(dif))
				normaInf[4] = abs(dif);
		}
		for (size_t ch = 0; ch < cascata.size(); ch++)
		{
			nt = (ptr_x_til[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
			if ( t >= sistema_a->GetTt2() )
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				delta = t*nt + flag2*cen*cascata[ch].size();
			}	
			else
				delta = t*nt;

			for (size_t r = 0; r < cascata[ch].size(); r++)
			{
				dif = - ptr_x_til[ch][r + delta] + ptr_x_til[CH + t*R + cascata[ch].at(r)][1];		//v
				norma1[2] += abs(dif);
				norma2[2] += pow(dif, 2);
				if (normaInf[2] < abs(dif))
					normaInf[2] = abs(dif);
				dif = - ptr_x_til[ch][r + cascata[ch].size() + delta] + ptr_x_til[CH + t*R + cascata[ch].at(r)][2];		//d
				norma1[3] += abs(dif);
				norma2[3] += pow(dif, 2);
				if (normaInf[3] < abs(dif))
					normaInf[3] = abs(dif);
			}
		}
		deltaa -= sistema_a->termeletricasVtr.size();
		deltaa += nat;
	}
	for (size_t i = 0; i < norma2.size(); i++)
	{
		norma2[i] = sqrt(norma2[i]);
	}
}
void Spcdec3Results::CalcularNormaSubgXtil(vetorfloat &normas)
{
	// Um valor de norma para cada tipo de restrição relaxada, norma 1, 2 e Inf;
	int delta = 0;
	int nt;
	int cen;
	int nat = (nd / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	double dif = 0;
	// zerar elementos do vetor
	for (int i = 0; i < normas.size(); i++)
		normas[i] = 0;

	// lambda [pt ph v d (phmax || res)]
	// x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)]

	if ( Aggrgtd )
	{
		// aqui ptr_x_til tem a mesma estrutura do ptr_x_subp_agg
		int delta_sg = 0;
		int JJ = 0;
		for (size_t r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		nt = (n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;
		int jj = 0;
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					dif = - ptr_x_til[0][n + i + delta_sg] + ptr_x_til[0][i + delta];		// - pta + pt
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			else
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					dif = - ptr_x_til[0][n + i + delta_sg] + ptr_x_til[0][i + delta] + ptr_x_til[0][i + I + delta]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			double som = 0;
			delta_sg += I;
			delta += (3+flag4+flag7)*I + flag1a*(sistema_a->barrasVtr.size() - 1);
			jj = delta;
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				dif = - ptr_x_til[0][n + r + delta_sg] + ptr_x_til[0][r + delta];		//ph
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
				dif = - ptr_x_til[0][n + r + R + delta_sg] + ptr_x_til[0][r + R + delta];		//v
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
				dif = - ptr_x_til[0][n + r + 2*R + delta_sg] + ptr_x_til[0][r + 2*R + delta];		//d
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
				if (sistema_a->GetFlagPhmax() == 1)
				{
					dif = - ptr_x_til[0][n + r + 4*R + delta_sg] + ptr_x_til[0][r + 4*R + delta];		//phmax
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
					delta += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				else
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						som += ptr_x_til[0][r + 4*R + 2*JJ + jj + j]*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();		// z
					som += - ptr_x_til[0][r + delta];
					jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
			}
			delta_sg -= I;
			delta -= (3+flag4+flag7)*I + flag1a*(sistema_a->barrasVtr.size() - 1);
			delta -= JJ;
			if (sistema_a->GetFlagPhmax() == 0)
			{
				dif = min(0.0, som - sistema_a->GetReserva(t));		//reserva: sum_r(sum_j(phg_max*z_j) - ph_r) - Res
				//dif = som - sistema_a->GetReserva(t);		//reserva: sum_r(sum_j(phg_max*z_j) - ph_r) - Res
				// each positive entry of the vector Ax - b is set to zero
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
			}
			else
				delta_sg += R;
			
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size();
			else
				delta += nt;
			delta_sg += I + 3*R;		// delta para o numero de variaveis duplicadas xa
		}
		normas[1] = sqrt(normas[1]);
	}
	else
	{
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					dif = - ptr_x_til[t + CH + N*R + I][i] + ptr_x_til[i + CH + N*R][t*(3 + flag4)];		// - pta + pt
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			else
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					dif = - ptr_x_til[t + CH + N*R + I][i] + ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_til[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			double som = 0;
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				dif = - ptr_x_til[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_til[CH + t*R + r][0];		//ph
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
				if (sistema_a->GetFlagPhmax() == 1)
				{
					dif = - ptr_x_til[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()] + ptr_x_til[CH + t*R + r][4];		//phmax
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
				else
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						som += ptr_x_til[CH + t*R + r][4 + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + j]*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();		// z
					som += - ptr_x_til[CH + t*R + r][0];
				}
			}
			if (sistema_a->GetFlagPhmax() == 0)
			{
				dif = min(0.0, som - sistema_a->GetReserva(t));		//reserva: sum_r(sum_j(phg_max*z_j) - ph_r) - Res
				//dif = som - sistema_a->GetReserva(t);		//reserva: sum_r(sum_j(phg_max*z_j) - ph_r) - Res
				// each positive entry of the vector Ax - b is set to zero
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
			}
			for (size_t ch = 0; ch < cascata.size(); ch++)
			{
				nt = (ptr_x_til[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
				if ( t >= sistema_a->GetTt2() )
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					delta = t*nt + flag2*cen*cascata[ch].size();
				}	
				else
					delta = t*nt;

				for (size_t r = 0; r < cascata[ch].size(); r++)
				{
					dif = - ptr_x_til[ch][r + delta] + ptr_x_til[CH + t*R + cascata[ch].at(r)][1];		//v
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
					dif = - ptr_x_til[ch][r + cascata[ch].size() + delta] + ptr_x_til[CH + t*R + cascata[ch].at(r)][2];		//d
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			}
		}
		normas[1] = sqrt(normas[1]);
	}
}
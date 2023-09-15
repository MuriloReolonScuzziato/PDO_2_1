#include "ResultadosConj.h"

CMatrizEsparsa ResultadosConj::MatrizCalcFluxo(int n_a)
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
CMatrizEsparsa ResultadosConj::MatrizCalcFluxo_cen(int n_a)
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
}

ResultadosConj::ResultadosConj(CSistema * const sistema_end)
{
	sistema_a = sistema_end;
	inFile = NULL;
	flag1 = int (1 - sistema_a->GetFlagBarraUnica());
	flag2 = int (sistema_a->GetFlagVolumeMeta());
	flag3 = int (sistema_a->GetFlagValorAgua());
	if (sistema_a->GetFlagAproxCustoT() > 1)
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
	n = N * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	na = N * (sistema_a->termeletricasVtr.size() + 4*sistema_a->hidreletricasVtr.size());
	if ( sistema_a->GetFlagBarraUnica() == true )		// Remover teta, def e deixar uma variável de deficit por período
		n -= N * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	if ( sistema_a->GetFlagVolumeMeta() == false )	// Remover vfol
		n -= sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
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
	if (sistema_a->GetFlagAproxCustoT() > 1)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size() ; i++)
		{
			aux << i + 1;
			stringVariaveis.push_back("F" + aux.str());
			aux.str("");
		}
	}
	int bb, bbb;
	if (sistema_a->GetFlagBarraUnica() == false)
	{
		for (size_t b = 1; b < sistema_a->barrasVtr.size() + 1; b++)
		{
			//aux = "teta";
			//bb = b/10;
			//bbb = b%10;
			//if ( b < 10)
			//	aux.push_back(char (48 + b));
			//else
			//{
			//	aux.push_back(char (48 + bb));
			//	aux.push_back(char (48 + bbb));
			//}
			//stringVariaveis.push_back(aux);

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
	if (sistema_a->GetFlagBarraUnica() == false)	
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
	if (sistema_a->GetFlagVolumeMeta() == true)
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		{
			aux << r + 1;
			stringVariaveis.push_back("vfol" + aux.str());
			aux.str("");
		}	
}

ResultadosConj::~ResultadosConj(void)
{
	delete inFile;
}

// ED
void ResultadosConj::CarregarResultadosED(double obj_a, vetorfloat x_a, vetorfloat lambda_a, int status_a, double mipgap_a)
{
	obj = obj_a;
	x = x_a;
	lambda = lambda_a;
	status = status_a;
	mipgap = mipgap_a;

	CIO = 0;		// Custo imediato de operação (custo do primeiro estágio), considerando o custo de operação e partida das termelétricas
	CTL = 0;		// Custo de operação das termelétricas linear
	CTQ = 0;		// Custo de operação das termelétricas quadrático
	CFO = 0;		// Custo futuro de operação ( valor incremental da água, da função de custo futuro)
	DEF = 0;		// Somatorio dos déficits
	VFOL = 0;		// Somatorio dos vfolga

	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	int nt;
	n = x.size();
	int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	double deltaT;
	int delta = 0;
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
				if (sistema_a->GetFlagTbinaryModel() == 0)
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT;		// Custo de operaçao quadrático
					if (sistema_a->GetFlagAproxCustoT() > 1)
					{
						CTL += x[i + 3*sistema_a->termeletricasVtr.size() + delta];		// Custo de operaçao linear
						CIO += x[i + 2*sistema_a->termeletricasVtr.size() + delta] + x[i + 3*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
					}
					else
						CIO += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
				}
				else
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT;		// Custo de operaçao quadrático
					if (sistema_a->GetFlagAproxCustoT() > 1)
					{
						CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta];		// Custo de operaçao linear
						CIO += x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
					}
					else
						CIO += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		// Custo de 1 estagio
				}
			}
		}
		else
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				{
					deltaT = sistema_a->GetDeltaT2();
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					{
						if (sistema_a->GetFlagTbinaryModel() == 0)
						{
							CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta])*deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
							if (sistema_a->GetFlagAproxCustoT() > 1)
							{
								CTL += x[i + 3*sistema_a->termeletricasVtr.size() + delta] * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
								CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] + x[i + 3*sistema_a->termeletricasVtr.size() + delta]) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
							}
							else
								CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) + x[i + 2*sistema_a->termeletricasVtr.size() + delta] ) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
						}
						else
						{
							CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
							if (sistema_a->GetFlagAproxCustoT() > 1)
							{
								CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta]* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
								CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta])* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
							}
							else
								CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida())* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
						}
					}
				}
			}
		}
		delta = delta + nt;
	}
	
	//if (sistema_a->GetFlagValorAgua() == 1)
	//	{
	//		nt = (n_a / T);
	//		delta = nt*sistema_a->GetTt2() - (nt - ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()));
	//		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	//		{
	//			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
	//			{
	//				CFO += - sistema_a->hidreletricasVtr[r].GetCustoAgua() * x[r + delta] * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
	//			}
	//			delta = delta + nt*(sistema_a->GetTt2() - sistema_a->GetTt1()) + flag2*sistema_a->hidreletricasVtr.size();
	//		}
	//	}

	delta = (n_a / T) - (1 - flag1) - flag1*sistema_a->barrasVtr.size(); //(n_a / T) = 34 52
	for (int t = 0; t < T; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			nt = (n_a / T);
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			deltaT = sistema_a->GetDeltaT1();
			for (size_t i = 0; i < (1 - flag1) + flag1*sistema_a->barrasVtr.size(); i++)
				DEF += x[i + delta];		// Def
		}
		else
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				{
					deltaT = sistema_a->GetDeltaT2();
					for (size_t i = 0; i < (1 - flag1) + flag1*sistema_a->barrasVtr.size(); i++)
						DEF += x[i + delta];					// Def
				}
			}
		}
		delta = delta + nt;
	}

	delta = (n_a / T) * sistema_a->GetTt2();
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			VFOL += x[delta++];					// vfol
		delta += (n_a / T) * (sistema_a->GetTt2() - sistema_a->GetTt1());
	}
}
void ResultadosConj::EscreverArquivo(string nome_arquivo, double tempo_modelo, double tempo_resol)
{
	int flag2 = int (sistema_a->GetFlagVolumeMeta());
	size_t n = x.size();
	size_t nt = stringVariaveis.size() - flag2*sistema_a->hidreletricasVtr.size();
	int num_cen = sistema_a->GetNCenarios();
	//num_cen = 1;		// para escrever resultado de um cenario
	int T1 = sistema_a->GetTt1();
	int T2 = sistema_a->GetTt2();
	inFile = new ofstream( nome_arquivo, ios::out );
	if ( inFile->is_open() )
	{
		if (status != 2)
		{
			*inFile << "Solução não ótima! " << "Status final = " << status << char(9) << "MIPgap = " << mipgap << endl;
			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagBarraUnica() << char(9) << sistema_a->GetFlagVolumeMeta() << char(9) << sistema_a->GetFlagValorAgua() << char(9) << sistema_a->GetFlagAproxCustoT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
		}
		else
		{
			*inFile << "Solução ótima! " << "Status final = " << status << char(9) << "MIPgap = " << mipgap << endl;
			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagBarraUnica() << char(9) << sistema_a->GetFlagVolumeMeta() << char(9) << sistema_a->GetFlagValorAgua() << char(9) << sistema_a->GetFlagAproxCustoT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
		}
		*inFile << endl;
		// Estágio 1
		*inFile << setw(10) << left << setprecision(6) << "Estágio 1" << endl;
		*inFile << setw(10) << left << setprecision(6) << "Periodo";
		for (int t = 0; t < T1; t++)
			*inFile << setw(14) << left << char(9) << t + 1;
		for (size_t i = 0; i < nt; i++)
		{
			*inFile << endl;
			*inFile << setw(10) << left << stringVariaveis[i];
			for (int t = 0; t < T1; t++)
				*inFile << char(9) << setw(15) << right << x[i + t*nt];
		}
		// Estágio 2
		for (int n_c = 0; n_c < num_cen; n_c++)
		{
			*inFile << endl << endl;
			*inFile << setw(10) << left << setprecision(6) << "Estágio 2 - Cenário "<< n_c + 1 << endl;
			*inFile << setw(10) << left << setprecision(6) << "Periodo";
			for (int t = T1; t < T2; t++)
				*inFile << setw(14) << left << char(9) << t + 1;
			for (size_t i = 0; i < nt; i++)
			{
				*inFile << endl;
				*inFile << setw(10) << left << stringVariaveis[i];
				for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1); t++)
				{
					*inFile << char(9) << setw(15) << right << x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()];
				}
			}
			if (sistema_a->GetFlagVolumeMeta() == 1)
				for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
				{
					*inFile << endl;
					*inFile << setw(10) << left << stringVariaveis[i + nt];
					for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1) - 1; t++)
						*inFile << char(9) << setw(15) << right << " - ";
					*inFile << char(9) << setw(15) << right << x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()];
				}
		}
		// Fluxo nas linhas
		if (sistema_a->GetFlagBarraUnica() == false)
		{
			vetorfloat teta, fluxo;
			fluxo.resize(sistema_a->linhasVtr.size()*(T1 + num_cen * (T2 - T1)) );
			CMatrizEsparsa M(0);
			M = MatrizCalcFluxo(int (x.size()) );		//Multiplicar matriz MatrizLimFluxo pelo vetor x
			fluxo = M.MultiplicarPorVetorDenso(&x);
		
			*inFile << endl << endl;
			*inFile << setw(10) << left << setprecision(5) << "Fluxo nas linhas (MW)" << endl;
			*inFile << setw(10) << left << setprecision(5) << "Estágio 1" << endl;
			for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
			{
				*inFile << setprecision(5) << "L" ;
				*inFile << setprecision(5) << l + 1;
				*inFile << setprecision(1) << std::fixed << " (" << sistema_a->linhasVtr[l].GetCapacidade() << ")" << setw(8);
				*inFile << setprecision(5) << std::scientific;
				for (int t = 0; t < T1; t++)
				{
					if (fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade())
						*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()] << "*";
					else
						*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()];

				}
				*inFile << endl;
			}
			for (int n_c = 0; n_c < num_cen; n_c++)
			{
				*inFile << setw(10) << left << setprecision(6) << "Estágio 2 - Cenário "<< n_c + 1 << endl;
				for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
				{
					*inFile << setprecision(5) << "L" ;
					*inFile << setprecision(5) << l + 1;
					*inFile << setprecision(1) << std::fixed << " (" << sistema_a->linhasVtr[l].GetCapacidade() << ")" << setw(8);
					*inFile << setprecision(5) << std::scientific;
					for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1); t++)
					{
						if ((fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade()) ||
							(fluxo[l + t*sistema_a->linhasVtr.size()] <= - 0.9999*sistema_a->linhasVtr[l].GetCapacidade()))
							*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()] << "*";
						else
							*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()];
					}
					*inFile << endl;
				}
			}
			//*inFile << endl << endl;
			//for (int l = 0; l < M.GetNlin(); l++)
			//	*inFile << setprecision(6) << "L" << setw(3) << l + 1 << ": " << setw(10) << left << fluxo[l] << endl;
		}
		inFile->close();
	}
	else
		cout << "Unable to open file";
	
	//inFile = NULL;
	//delete inFile;
}

// DecEsp
void ResultadosConj::CarregarResultadosDecEsp(double obj_a, vetorfloat x_a, int status_a, int tipoSubproblema, int id1, int id2)
{
	status = status_a;

	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();

	int delta, deltaa, delta_subP;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int nt_subP = (x_a.size() / N);
	// Tipo do subproblema: (0 = Hidraulico, 1 = Hidreletrico, 2 = Termeletrico, 3 = Demanda)
	switch (tipoSubproblema)
	{
	case 0:
		{
			// id1 é o índice da cascata e id2 é nada
			//delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + flag1*sistema_a->barrasVtr.size() + (1 - flag1);
			deltaa = sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size();
			delta_subP = 0;
			nt_subP = (x_a.size() - flag2*sistema_a->GetNCenarios()*cascata[id1].size()) / N;
			// Alocar resultados do subproblema
			// xa e x
			for (int t = 0; t < N; t++)
			{
				for (size_t r = 0; r < cascata[id1].size(); r++)
				{
					xa[cascata[id1].at(r) + deltaa] = x_a[r + delta_subP];		//v
					xa[cascata[id1].at(r) + sistema_a->hidreletricasVtr.size() + deltaa] = x_a[r + cascata[id1].size() + delta_subP];		//d
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
						x[delta + cascata[id1].at(r)] = x_a[delta_subP + r];		//vfol
					delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
					delta_subP += nt_subP * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[id1].size();
				}
			}
			// fo
			obj_subp[id1] = obj_a;
			break;
		}
	case 1:
		{
			// id1 é o índice da usina e id2 é o índice do nó
			delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
			deltaa = 0;
			delta_subP = 0;
			int Jr = 0;
			for (size_t r = 0; r < id1; r++)
				Jr += sistema_a->hidreletricasVtr[r].GetNGrupos();

			delta += id2 * nt;
			if (id2 >= sistema_a->GetTt1())
				delta += ((id2 - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1())) * flag2*sistema_a->hidreletricasVtr.size();
			// Alocar resultados do subproblema
			// xa e x
			x[id1 + delta] = x_a[0];		//ph
			x[id1 + 1*sistema_a->hidreletricasVtr.size() + delta] = x_a[1];		//v
			x[id1 + 2*sistema_a->hidreletricasVtr.size() + delta] = x_a[2];		//d
			x[id1 + 3*sistema_a->hidreletricasVtr.size() + delta] = x_a[3];		//s
			x[id1 + 4*sistema_a->hidreletricasVtr.size() + delta] = x_a[4];		//phmax
			for (int j = 0; j < sistema_a->hidreletricasVtr[id1].GetNGrupos(); j++)
			{
				x[j + 5*sistema_a->hidreletricasVtr.size() + Jr + delta] = x_a[j + 5];					//phg
				x[j + 5*sistema_a->hidreletricasVtr.size() + JJ + Jr + delta] = x_a[j + sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];		//q
				x[j + 5*sistema_a->hidreletricasVtr.size() + 2*JJ + Jr + delta] = x_a[j + 2*sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];	//z
			}
			// fo
			obj_subp[CH + id2*R + id1] = obj_a;
			break;
		}
	case 2:
		{
			// id1 é o índice da termeletrica e id2 é nada
			delta = 0;
			delta_subP = 0;
			// Alocar resultados do subproblema
			// xa e x
			for (int t = 0; t < N; t++)
			{
				x[id1 + delta] = x_a[delta_subP];														//pt
				x[id1 + 1*sistema_a->termeletricasVtr.size() + delta] = x_a[1 + delta_subP];			//u
				x[id1 + 2*sistema_a->termeletricasVtr.size() + delta] = x_a[2 + delta_subP];			//cp
				if (sistema_a->GetFlagAproxCustoT() > 1)
				{
					x[id1 + 3*sistema_a->termeletricasVtr.size() + delta] = x_a[3 + delta_subP];	//F
				}

				if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
					delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
				else
					delta += nt;
				delta_subP += nt_subP;
			}
			// fo
			obj_subp[CH + N * R + id1] = obj_a;
			break;
		}
	case 3:
		{
			// id1 é o índice do nó e id2 é nada
			//x_a = [pt teta ph phmax def]
			delta = (3 + flag4)*sistema_a->termeletricasVtr.size();
			delta += id1 * nt;
			if (id1 >= sistema_a->GetTt1())
				delta += ((id1 - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1())) * flag2*sistema_a->hidreletricasVtr.size();
			deltaa = id1 * nat;
			// Alocar resultados do subproblema
			// xa e x
			
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				xa[i + deltaa] = x_a[i];		//pt
			if (flag1 == 1)		// considerar rede
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
					x[b + delta] = x_a[b + sistema_a->termeletricasVtr.size()];	//teta
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					x[b + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = x_a[b + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
			}
			else
				x[flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = x_a[sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				xa[r + sistema_a->termeletricasVtr.size() + deltaa] = x_a[r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1)];		//ph
				xa[r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + deltaa] = x_a[r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];		//phmax
			}
			// fo
			obj_subp[CH + N * R + I + id1] = obj_a;
			break;
		}
	}
	
}
double ResultadosConj::GetFobj()
{
	obj = 0;
	for (int sp = 0; sp < CH + R*N + I + N; sp++)
		obj += obj_subp[sp];
	return obj;
}

void ResultadosConj::GetSubGradiente(int comp, double * SubG)
{
	for (int i = 0; i < na; i++)
		SubG[i] = 0;
	// Selecionar tipo do subproblema pelo valor componente: (Hidraulico, Hidreletrico, Termeletrico, Demanda)
	// 1 <= wFi (comp) <= GetNrFi();
	int delta = 0;
	int deltaa = 0;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);

	if ( (comp > 0) && (comp <= CH) )	//Hidraulico ( 1 <= comp <= CH )
	{
		int ch = comp - 1;
		deltaa = sistema_a->termeletricasVtr.size();
		for (int t = 0; t < N; t++)
		{
			for (size_t r = 0; r < cascata[ch].size(); r++)
			{
				SubG[cascata[ch].at(r) + sistema_a->hidreletricasVtr.size() + deltaa] = - xa[cascata[ch].at(r) + sistema_a->hidreletricasVtr.size() + deltaa];			//v
				SubG[cascata[ch].at(r) + 2*sistema_a->hidreletricasVtr.size() + deltaa] = - xa[cascata[ch].at(r) + 2*sistema_a->hidreletricasVtr.size() + deltaa];		//d
			}
			deltaa += nat;
		}
	}
	else if ( (comp > CH) && (comp <= CH + N * R) )		//Hidreletrico
	{
		int r = (comp - CH - 1) % int (sistema_a->hidreletricasVtr.size());
	    int no = (comp - CH - 1) / int (sistema_a->hidreletricasVtr.size());
		int cen = 0;
		if (no >= sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
			cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		
		delta = (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
		deltaa = sistema_a->termeletricasVtr.size();
		delta += no * nt + cen * flag2 * sistema_a->hidreletricasVtr.size();
		deltaa += no * nat;
		
		SubG[r + deltaa] =  x[r + delta];		//ph
		SubG[r + sistema_a->hidreletricasVtr.size() + deltaa] = x[r + sistema_a->hidreletricasVtr.size() + delta];			//v
		SubG[r + 2*sistema_a->hidreletricasVtr.size() + deltaa] = x[r + 2*sistema_a->hidreletricasVtr.size() + delta];		//d
		SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = x[r + 4*sistema_a->hidreletricasVtr.size() + delta];		//phmax
	}
	else if ( (comp > CH + N * R) && (comp <= CH + N * R + I) )		//Termeletrico
	{
		int i = comp - CH - N * R - 1;
		deltaa = 0;
		delta = 0;
		for (int t = 0; t < N; t++)
		{
			SubG[i + deltaa] = x[i + delta];		//pt

			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
			deltaa += nat;
		}
	}
	else if ( (comp > CH + N * R + I) && (comp <= CH + N * R + I + N) )		//Demanda
	{
		int no = comp - CH - N * R - I - 1;
		int cen = 0;
		if (no >= sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
			cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		//delta = no * nt + cen * flag2 * sistema_a->hidreletricasVtr.size();
		deltaa = no * nat;

		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			SubG[i + deltaa] = - xa[i + deltaa];		//pt
		//delta += (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - xa[r + deltaa];		//ph
			SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = - xa[r + 3*sistema_a->hidreletricasVtr.size() + deltaa];		//phmax
		}
	}
}
void ResultadosConj::GetSubGradiente(double * SubG)
{
	AlocarXeXa();
	int delta = 0;
	int deltaa = 0;
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	for (int t = 0; t < N; t++)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			SubG[i + deltaa] = - xa[i + deltaa] + x[i + delta];		//pt
		delta += (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - xa[r + deltaa] + x[r + delta];		//ph
			SubG[r + sistema_a->hidreletricasVtr.size() + deltaa] = - xa[r + sistema_a->hidreletricasVtr.size() + deltaa] + x[r + sistema_a->hidreletricasVtr.size() + delta];		//v
			SubG[r + 2*sistema_a->hidreletricasVtr.size() + deltaa] = - xa[r + 2*sistema_a->hidreletricasVtr.size() + deltaa] + x[r + 2*sistema_a->hidreletricasVtr.size() + delta];		//d
			SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = - xa[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] + x[r + 4*sistema_a->hidreletricasVtr.size() + delta];		//phmax
		}
		delta -= ((3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1));
		deltaa -= sistema_a->termeletricasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
		else
			delta += nt;
		deltaa += nat;
	}
}

void ResultadosConj::ExportarX(string nome_arquivo)
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
void ResultadosConj::ExportarXA(string nome_arquivo)
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
void ResultadosConj::ExportarXmed(string nome_arquivo)
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
vetorfloat ResultadosConj::GetX_med()
{
	// Alocar x_med a partir de x e xa
	AlocarX_med();

	return x_med;
}
void ResultadosConj::AlocarX_med()
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
	//
	for (int t = 0; t < N; t++)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)		//pt
			x_med[i + delta] = (x[i + delta] + xa[i + deltaa]) / 2;
		delta += (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			x_med[r + delta] = (x[r + delta] + xa[r + deltaa]) / 2;		//ph
			x_med[r + sistema_a->hidreletricasVtr.size() + delta] = (x[r + sistema_a->hidreletricasVtr.size() + delta] + xa[r + sistema_a->hidreletricasVtr.size() + deltaa]) / 2;				//v
			x_med[r + 2*sistema_a->hidreletricasVtr.size() + delta] = (x[r + 2*sistema_a->hidreletricasVtr.size() + delta] + xa[r + 2*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//d
			x_med[r + 4*sistema_a->hidreletricasVtr.size() + delta] = (x[r + 4*sistema_a->hidreletricasVtr.size() + delta] + xa[r + 3*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//phmax
		}
		delta -= (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
		deltaa -= sistema_a->termeletricasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
		else
			delta += nt;
		deltaa += nat;
	}
}

// DecEsp_2.0

// Fazer outra funçao para guardar resultados dos subproblemas, por componente (guardar x no formato de cada subproblema)
// fazer tb outra funçao para calcular o subgradiente de acordo com o formato do x acima
// só colocar os dados nos formatos do x e xa (e x_med) originais no resultado final!!!

// subproblema HA comp = ch;					//x = [v d vfol]
// subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó	//x = [ph v d s phmax phg q z]
// subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica	//x = [pt u cp F]
// subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó	//x = [pt teta ph phmax def]

void ResultadosConj::SetComponents( int * comp_inf, bool aggr )
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
void ResultadosConj::GravarSolucao(double obj_a, vetorfloat x_a, int status_a, int comp)
{
	// Grava o resultado de cada subproblema
	// para o subproblema HA comp = ch;
	// para o subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó
	// para o subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica
	// para o subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó
	status = status_a;
	// wFi = 0 é o termo constante entao wFi = 1 corresponde à comp = 0;
	ptr_x_subp[comp] = x_a;		// talvez fique mais rápido receber em cada subproblema o endereço do vetorfloat de ptr_x_subp e preenche-lo (clear, resize e atribuir valor)
	obj_subp[comp] = obj_a;
}
void ResultadosConj::GetSubGradiente2(int wFi_a, double * SubG)
{
	for (int i = 0; i < na; i++)
		SubG[i] = 0;
	
	// O vetor de subgradientes tem a dimensao do vetor lambda, portanto a posicao do subg q sera preenchido depende do tipo de problema
	// Selecionar tipo do subproblema pelo valor componente: (Hidraulico, Hidreletrico, Termeletrico, Demanda)
	// 1 <= wFi (comp) <= GetNrFi();
	// lambda [pt ph v d phmax]
	int delta = 0;
	int deltaa = 0;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);

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
	else if ( (wFi_a > CH) && (wFi_a <= CH + N * R) )		//Hidreletrico	//x = [ph v d s phmax phg q z]
	{
		int r = (wFi_a - CH - 1) % int (sistema_a->hidreletricasVtr.size());
	    int no = (wFi_a - CH - 1) / int (sistema_a->hidreletricasVtr.size());
		
		deltaa = sistema_a->termeletricasVtr.size();
		deltaa += no * nat;

		SubG[r + deltaa] =  ptr_x_subp[wFi_a - 1][0];											//ph
		SubG[r + sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[wFi_a - 1][1];		//v
		SubG[r + 2*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[wFi_a - 1][2];		//d
		SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[wFi_a - 1][4];		//phmax
	}
	else if ( (wFi_a > CH + N * R) && (wFi_a <= CH + N * R + I) )		//Termeletrico	//x = [pt u cp F]
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
			SubG[r + deltaa] = - ptr_x_subp[wFi_a - 1][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1)];			//ph
			SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[wFi_a - 1][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];	//phmax
		}
	}
}
void ResultadosConj::GetSubGradiente2(double * SubG)
{
	int deltaa = 0;
	int delta = 0;
	int cen = 0;
	int nt;
	int nat = (na / N);

	for (int t = 0; t < N; t++)
	{
		if (sistema_a->GetFlagTbinaryModel() == 0)
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				SubG[i + deltaa] = - ptr_x_subp[t + CH + N*R + I][i] + ptr_x_subp[i + CH + N*R][t*(3 + flag4)];		//pt
		else
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				SubG[i + deltaa] = - ptr_x_subp[t + CH + N*R + I][i] + ptr_x_subp[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_subp[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - ptr_x_subp[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1)] + ptr_x_subp[CH + t*R + r][0];		//ph
			SubG[r + 3*sistema_a->hidreletricasVtr.size() + deltaa] = - ptr_x_subp[t + CH + N*R + I][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()] + ptr_x_subp[CH + t*R + r][4];		//phmax
		}
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
}

vetorfloat ResultadosConj::GetCompX(int comp)
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
void ResultadosConj::AlocarPtrXHat()
{
	// Alocar ptr_x_hat a partir de ptr_x_subp para o caso agregado
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();

	int delta, deltaa, delta_subP;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int nt_subP;
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

		delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
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
		ptr_x_subp_agg[0][id1 + 4*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][4];		//phmax
		for (int j = 0; j < sistema_a->hidreletricasVtr[id1].GetNGrupos(); j++)
		{
			ptr_x_subp_agg[0][j + 5*sistema_a->hidreletricasVtr.size() + Jr + delta] = ptr_x_subp[comp][j + 5];					//phg
			ptr_x_subp_agg[0][j + 5*sistema_a->hidreletricasVtr.size() + JJ + Jr + delta] = ptr_x_subp[comp][j + sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];		//q
			ptr_x_subp_agg[0][j + 5*sistema_a->hidreletricasVtr.size() + 2*JJ + Jr + delta] = ptr_x_subp[comp][j + 2*sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];	//z
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
			if (sistema_a->GetFlagAproxCustoT() > 1)
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
				ptr_x_subp_agg[0][b + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		}
		else
			ptr_x_subp_agg[0][flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			ptr_x_subp_agg[0][n + r + sistema_a->termeletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1)];		//ph
			ptr_x_subp_agg[0][n + r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];		//phmax
		}
	}
}
void ResultadosConj::AlocarXeXa()
{
	// Alocar x e xa a partir de ptr_x_subp
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();

	int delta, deltaa, delta_subP;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int nt_subP;
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

		delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
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
		x[id1 + 4*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][4];		//phmax
		for (int j = 0; j < sistema_a->hidreletricasVtr[id1].GetNGrupos(); j++)
		{
			x[j + 5*sistema_a->hidreletricasVtr.size() + Jr + delta] = ptr_x_subp[comp][j + 5];					//phg
			x[j + 5*sistema_a->hidreletricasVtr.size() + JJ + Jr + delta] = ptr_x_subp[comp][j + sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];		//q
			x[j + 5*sistema_a->hidreletricasVtr.size() + 2*JJ + Jr + delta] = ptr_x_subp[comp][j + 2*sistema_a->hidreletricasVtr[id1].GetNGrupos() + 5];	//z
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
			if (sistema_a->GetFlagAproxCustoT() > 1)
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
				x[b + flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		}
		else
			x[flag1*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			xa[r + sistema_a->termeletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1)];		//ph
			xa[r + sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];		//phmax
		}
	}
}
vetorfloat ResultadosConj::GetX_med(bool hat_or_til)
{
	// Alocar x_med a partir de ptr_x_hat ou ptr_x_til (estrutura depende se o modelo é agregado ou desagregado)
	// ptr_x_hat ou ptr_x_til tem a mesma estrutura!!
	// A cada chamada do SetGiName o ptr_x_hat aqui nessa classe é atualizado, ou seja, ptr_x_hat é atualizado right after Fi() and all the GetGi()

	// hat_or_til indica a partir de qual solução é calculado o x_med
	// 0 a partir da soluçao da RL, ptr_x_hat
	// 1 a partir da soluçao convexifica, ptr_x_til
	vetorfloat * ptr_x_a;
	if ( hat_or_til )
		ptr_x_a = &ptr_x_til[0];
	else
		if ( Aggrgtd )
			ptr_x_a = &ptr_x_subp_agg[0];
		else
			ptr_x_a = &ptr_x_subp[0];

	// Para as variáveis duplicadas x_med = (x + xa) / 2
	//x = [pt u cp F teta ph v d s phmax phg q z def]
	//xa = [pta pha va da phmaxa]
	int delta, deltaa;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
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
					x_med[i + delta] = ((ptr_x_a[0][i + delta] + ptr_x_a[0][i + 1 + delta]*sistema_a->termeletricasVtr[i].GetPmin()) + ptr_x_a[0][n + i + deltaa]) / 2;
			}
			delta += (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
			deltaa += sistema_a->termeletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				x_med[r + delta] = (ptr_x_a[0][r + delta] + ptr_x_a[0][n + r + deltaa]) / 2;		//ph
				x_med[r + sistema_a->hidreletricasVtr.size() + delta] = (ptr_x_a[0][r + sistema_a->hidreletricasVtr.size() + delta] + ptr_x_a[0][n + r + sistema_a->hidreletricasVtr.size() + deltaa]) / 2;				//v
				x_med[r + 2*sistema_a->hidreletricasVtr.size() + delta] = (ptr_x_a[0][r + 2*sistema_a->hidreletricasVtr.size() + delta] + ptr_x_a[0][n + r + 2*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//d
				x_med[r + 4*sistema_a->hidreletricasVtr.size() + delta] = (ptr_x_a[0][r + 4*sistema_a->hidreletricasVtr.size() + delta] + ptr_x_a[0][n + r + 3*sistema_a->hidreletricasVtr.size() + deltaa]) / 2;		//phmax
			}
			delta -= (3+flag4)*sistema_a->termeletricasVtr.size() + flag1*(sistema_a->barrasVtr.size() - 1);
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
		// subproblema HE comp = CH + id2*R + id1;		// id1 é o índice da usina e id2 é o índice do nó	//x = [ph v d s phmax phg q z]
		// subproblema T comp = CH + N * R + id1;		// id1 é o índice da termeletrica	//x = [pt u cp F]
		// subproblema D comp = CH + N * R + I + id1;	// id1 é o índice do nó	//x = [pt teta ph phmax def]
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int jj, cen;
		delta = 0;
		deltaa = 0;
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)
			{
				for (int i = 0; i < I; i++)
				{
					x_med[i + delta] = (ptr_x_a[t + CH + N*R + I][i] + ptr_x_a[i + CH + N*R][t*(3 + flag4)]) / 2;		//pt
					x_med[i + I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4) + 1];									//u		
					x_med[i + 2*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4) + 2];									//cp
				}
				if (sistema_a->GetFlagAproxCustoT() > 1)
					for (int i = 0; i < I; i++)
						x_med[i + 3*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4) + 3];								//F
			}
			else
			{
				for (int i = 0; i < I; i++)
				{
					x_med[i + delta] = (ptr_x_a[t + CH + N*R + I][i] + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7)] + ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin()) / 2;		//pt
					x_med[i + I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 1];									//u		
					x_med[i + 2*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 2];									//up
					x_med[i + 3*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 3];									//ud
				}
				if (sistema_a->GetFlagAproxCustoT() > 1)
					for (int i = 0; i < I; i++)
						x_med[i + (3 + flag7)*I + delta] = ptr_x_a[i + CH + N*R][t*(3 + flag4 + flag7) + 3 + flag7];								//F
			}

			jj = I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R;
			for (int r = 0; r < R; r++)
			{
				x_med[r + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + delta] = (ptr_x_a[t + CH + N*R + I][r + I + flag1*(sistema_a->barrasVtr.size() - 1)] + ptr_x_a[CH + t*R + r][0]) / 2;				//ph
				x_med[r + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 3*R + delta] = ptr_x_a[CH + t*R + r][3];		//s
				x_med[r + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 4*R + delta] = (ptr_x_a[t + CH + N*R + I][r + I + flag1*(sistema_a->barrasVtr.size() - 1) + R] + ptr_x_a[CH + t*R + r][4]) / 2;	//phmax
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					x_med[j + jj + delta] = ptr_x_a[CH + t*R + r][j + 5];															//phg
					x_med[j + jj + JJ + delta] = ptr_x_a[CH + t*R + r][j + sistema_a->hidreletricasVtr[r].GetNGrupos() + 5];		//q
					x_med[j + jj + 2*JJ + delta] = ptr_x_a[CH + t*R + r][j + 2*sistema_a->hidreletricasVtr[r].GetNGrupos() + 5];	//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if (flag1 == 1)		// considerar rede
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
					x_med[b + I*(3 + flag4) + delta] = ptr_x_a[t + CH + N*R + I][b + I];	//teta
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					x_med[b + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta] = ptr_x_a[t + CH + N*R + I][b + I + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R];	//def
			}
			else
				x_med[I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + delta] = ptr_x_a[t + CH + N*R + I][I + 2*R];	//def

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
					x_med[cascata[ch].at(r) + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + R + delta] = (ptr_x_a[ch][r + deltaa] + ptr_x_a[CH + t*R + cascata[ch].at(r)][1]) / 2;							//v
					x_med[cascata[ch].at(r) + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 2*R + delta] = (ptr_x_a[ch][r + cascata[ch].size() + deltaa] + ptr_x_a[CH + t*R + cascata[ch].at(r)][2]) / 2;	//d
				}
			}

			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			else
				delta += nt;
		}
		if (flag2 == 1)
		{
			delta = nt * (sistema_a->GetTt2() - 1);
			for (size_t ch = 0; ch < cascata.size(); ch++)
			{
				nat = (ptr_x_a[ch].size() - flag2*sistema_a->GetNCenarios()*cascata[ch].size()) / N;
				deltaa = nat * (sistema_a->GetTt2() - 1);
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				{
					for (size_t r = 0; r < cascata[ch].size(); r++)
						x_med[cascata[ch].at(r) + I*(3 + flag4) + flag1*(sistema_a->barrasVtr.size() - 1) + 5*R + 3*JJ + flag1*sistema_a->barrasVtr.size() + (1 - flag1) + delta] = ptr_x_a[ch][r + 2*cascata[ch].size() + deltaa];		//vfol
					delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
					deltaa += nat * (sistema_a->GetTt2() - sistema_a->GetTt1()) + cascata[ch].size();
				}
			}
		}
	}
	//delete[] ptr_x_a;
	return x_med;
}
void ResultadosConj::ZerarSolucoes()
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
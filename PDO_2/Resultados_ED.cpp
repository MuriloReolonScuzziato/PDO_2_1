#include "Resultados_ED.h"

CMatrizEsparsa Resultados_ED::MatrizCalcFluxo(int n_a)
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
		for ( size_t ll = 0; ll < sistema_a->linhasVtr.size(); ll++)
			if (l == ll)
				TT.InserirElemento(int(l), int(ll), 100 / (sistema_a->linhasVtr[l].GetReatancia()));
	CMatrizEsparsa Alin(int (N * sistema_a->linhasVtr.size()), n_a);
	n_a = n_a - int(flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < N; t++)
	{
		Alin.InserirMatriz(ll, c + int ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()), ll + int (sistema_a->linhasVtr.size()) - 1, c + int((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
		ll = ll + int(sistema_a->linhasVtr.size());
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / N) + int(flag2*sistema_a->hidreletricasVtr.size());
		else
			c += (n_a / N);
	}
	return Alin;
}
CMatrizEsparsa Resultados_ED::MatrizCalcFluxo_cen(int n_a)
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
void Resultados_ED::CriarStringVars()
{
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
		for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
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
	{
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		{
			aux << r + 1;
			stringVariaveis.push_back("phmax" + aux.str());
			aux.str("");
		}
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
	{
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size() ; r++)
		{
			aux << r + 1;
			stringVariaveis.push_back("vfol" + aux.str());
			aux.str("");
		}	
	}
}
Resultados_ED::Resultados_ED(CSistema * const sistema_end)
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
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= N * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= N * (sistema_a->barrasVtr.size() - 1);
	if ( sistema_a->GetFlagVfol() == false )		// Remover vfol
		n -= sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	x.resize(n);
	obj = 0;
	
	// Criar vetor coluna com o nome das variáveis e tamanho do numero de variáveis em cada período
	CriarStringVars();
}
Resultados_ED::~Resultados_ED(void)
{
	delete inFile;
}

// ED
void Resultados_ED::CarregarResultadosED(double obj_a, vetorfloat x_a, vetorfloat lambda_a, int status_a, double mipgap_a)
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
	int nt, cen;
	//n = x.size();
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
					if (sistema_a->GetFlagInitAproxCT() > 1)
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
					if (sistema_a->GetFlagInitAproxCT() > 1)
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
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			deltaT = sistema_a->GetDeltaT2();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				if (sistema_a->GetFlagTbinaryModel() == 0)
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta])*deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
					if (sistema_a->GetFlagInitAproxCT() > 1)
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
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta]* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
						CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta])* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
					}
					else
						CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida())* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
				}
			}
		}
		delta = delta + nt;
	}
	
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	delta = (n_a / T) - (1 - flag1d) - flag1d*sistema_a->barrasVtr.size(); //(n_a / T) = 34 52
	for (int t = 0; t < T; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			nt = (n_a / T) + flag2*sistema_a->hidreletricasVtr.size();
		else
			nt = (n_a / T);
		for (size_t i = 0; i < (1 - flag1d) + flag1d*sistema_a->barrasVtr.size(); i++)
			DEF += x[i + delta];		// Def
		delta = delta + nt;
	}

	if (sistema_a->GetFlagVfol() == 1)
	{
		delta = (n_a / T) * sistema_a->GetTt2();
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
		
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				VFOL += x[delta++];					// vfol
			delta += (n_a / T) * (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
	}
	else
		VFOL = 0;
}
void Resultados_ED::EscreverArquivo(string nome_arquivo, double tempo_modelo, double tempo_resol, string saida)
{
	int flag2 = int (sistema_a->GetFlagVfol());
	size_t n = x.size();
	size_t nt = stringVariaveis.size() - flag2*sistema_a->hidreletricasVtr.size();
	int num_cen = sistema_a->GetNCenarios();
	//num_cen = 1;		// para escrever resultado de um cenario
	int T1 = sistema_a->GetTt1();
	int T2 = sistema_a->GetTt2();
	//CreateDirectoryA(varString.c_str(), NULL);		// criar pasta
	inFile = new ofstream( saida + nome_arquivo, ios::out );
	if ( inFile->is_open() )
	{
		if (status != 2)
		{
			*inFile << "Solução não ótima! " << "Status final = " << status << char(9) << "MIPgap = " << mipgap << endl;
			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagModeloRede() << char(9) << sistema_a->GetFlagVfol() << char(9) << sistema_a->GetFlagPhmax() << char(9) << sistema_a->GetFlagInitAproxCT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
		}
		else
		{
			*inFile << "Solução ótima! " << "Status final = " << status << char(9) << "MIPgap = " << mipgap << endl;
			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagModeloRede() << char(9) << sistema_a->GetFlagVfol() << char(9) << sistema_a->GetFlagPhmax() << char(9) << sistema_a->GetFlagInitAproxCT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
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
			if (sistema_a->GetFlagVfol() == 1)
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
		vetorfloat fluxo;
		fluxo.resize(sistema_a->linhasVtr.size()*(T1 + num_cen * (T2 - T1)));
		if (sistema_a->GetFlagModeloRede() == 1)
		{
			vetorfloat teta;
			CMatrizEsparsa M(0);
			M = MatrizCalcFluxo(int (x.size()) );		//Multiplicar matriz MatrizLimFluxo pelo vetor x
			fluxo = M.MultiplicarPorVetorDenso(&x);
		}
		else	// sistema_a->GetFlagModeloRede() != 1
		{
			int delta = 0;
			int cen = 0;
			double D = 0;
			int JJ = 0;
			for (size_t i = 0; i < R; i++)
				JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
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
			Agh = - sistema_a->Beta * Agh;
			Agt = - sistema_a->Beta * Agt;
			Agtu = - sistema_a->Beta * Agtu;
			Adef = - sistema_a->Beta * Adef;
			VectorXd Ad(sistema_a->barrasVtr.size());
			VectorXd gg(sistema_a->linhasVtr.size());
			VectorXd g, Rhs;
			VectorXd flow(sistema_a->linhasVtr.size());
			for (int t = 0; t < N; t++)
			{
				Ad.resize(sistema_a->barrasVtr.size());
				if (sistema_a->GetTt1() <= t)
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				// Determinar o vetor de demandas
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				{
					D = 0;
					for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
						D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
					Ad(b) = D;
				}
				Rhs = - sistema_a->Beta * ( Ad );
				// Determinar o vetor de gerações
				gg = VectorXd::Zero(sistema_a->linhasVtr.size());
				g.resize(R);
				for (size_t r = 0; r < R; r++)
					g(r) = x[r + delta + (3+flag4+flag7)*I];
				gg = gg + Agh * g;
				g.resize(I);
				for (size_t i = 0; i < I; i++)
					g(i) = x[i + delta];
				gg = gg + Agt * g;
				g.resize(I);
				if (sistema_a->GetFlagTbinaryModel() == 1)
					for (size_t i = 0; i < I; i++)
						g(i) = x[i + delta + I];
				gg = gg + Agtu * g;
				g.resize(sistema_a->barrasVtr.size());
				if ( sistema_a->GetFlagModeloRede() == 0 )		// modelo de barra unica, deficit alocado igualmente entre todas barras, o que não é verdade!
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						g(b) = x[delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ] / sistema_a->barrasVtr.size();
				}
				else
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						g(b) = x[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
				}
				gg = gg + Adef * g;
				// Conferir limites das linhas extrapolados para cada periodo
				flow = - gg + Rhs;		// Beta ( d - g): fluxo nas linhas para o periodo t
				for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
					fluxo[l + t*sistema_a->linhasVtr.size()] = flow(l);
				if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
					delta += nt + flag2*R;
				else
					delta += nt;
			}
		}
		// Escrever fluxo
		*inFile << endl << endl;
		if ( sistema_a->GetFlagModeloRede() == 0 )	
		{
			*inFile << "Atenção: Modelo de barra única, caso exista déficit ele será alocado igualmente entre as barras! Portanto o resultado pode ser diferente do modelo de rede.";
			*inFile << endl << endl;
		}
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
				if ((fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade()) ||
					(fluxo[l + t*sistema_a->linhasVtr.size()] <= - 0.9999*sistema_a->linhasVtr[l].GetCapacidade()))
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
		inFile->close();
	}
	else
		cout << "Unable to open file";
	
	//inFile = NULL;
	//delete inFile;
}
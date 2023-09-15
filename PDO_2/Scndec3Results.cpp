#include "Scndec3Results.h"

Scndec3Results::Scndec3Results(CSistema * const sistema_end)
{
	sistema_a = sistema_end;
	inFile = NULL;
	flag1 = int (sistema_a->GetFlagModeloRede());
	flag2 = int (sistema_a->GetFlagVfol());
	flag3 = int (sistema_a->GetFlagPhmax());
	if (flag3 == 1)
		cout << "Scndec3Results: classe não preparada para modelo com variavel phmax! Desnecessario." << endl;

	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = sistema_a->GetFlagTbinaryModel();
	R = sistema_a->hidreletricasVtr.size();
	I = sistema_a->termeletricasVtr.size();
	N = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	n = N * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	nd = sistema_a->GetTt1() * (2*sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size());
	// numero de var. por cenário ou por subproblema!
	n_sub = sistema_a->GetTt2() * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size();
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
	{
		n -= N * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
		n_sub -= sistema_a->GetTt2() * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	}
	else if ( flag1 >= 2)	// Remover somente teta
	{
		n -= N * (sistema_a->barrasVtr.size() - 1);
		n_sub -= sistema_a->GetTt2() * (sistema_a->barrasVtr.size() - 1);
	}

	nd *= (sistema_a->GetNCenarios() - 1);
	x.resize(n_sub * sistema_a->GetNCenarios());
	x_med.resize(n);
	obj = 0;
	obj_subp.resize(sistema_a->GetNCenarios());
	epsilon_subp.resize(sistema_a->GetNCenarios());
	for (int i = 0; i <  sistema_a->GetNCenarios(); i++)
		obj_subp[i] = 0;
	ptr_x_subp.resize(sistema_a->GetNCenarios());
	//ptr_x_subp_agg.resize(1);
	//ptr_x_subp_agg[0].resize(n_cen * sistema_a->GetNCenarios());
	
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
Scndec3Results::~Scndec3Results(void)
{
	delete inFile;
}

void Scndec3Results::AlocarXmed(CMatrizEsparsa * x_spr)
{
	// Alocar x_med a partir de x
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nt1 = nt * sistema_a->GetTt1();		// numero de var. de primeiro estágio

	// variáveis de primeiro estágio (valor médio/esperado)
	// até ndt o valor das variáveis é dada pela média dos subproblemas
	for (int i = 0; i < nt1; i++)
	{
		x_spr->SubstituirElemento(i, 0, 0);
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			x_spr->SubstituirElemento(i, 0, x_spr->GetElemento(i, 0) + sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * x[n_cen*n_sub + i]);
	}
	// variáveis de segundo estágio
	// cada cenário tem sua própria variável, organizadas por ordem crescente no vetor x_med
	int n_var_2_est = (sistema_a->GetTt2() - sistema_a->GetTt1()) * nt + flag2*sistema_a->hidreletricasVtr.size();
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		for (int i = 0; i < n_var_2_est; i++)
			x_spr->SubstituirElemento(nt1 + i + n_var_2_est*n_cen, 0, x[n_cen*n_sub + nt1 + i]);
}

// ScnDec3

// subproblemas cen = 0,...,Ncen:
//X_cen = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
//X_cen = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
// Lambda:
//L_cen = [Lpt Lu Lv] -> (.) tamanho dos vetores que dependem de flags

void Scndec3Results::SetComponents( int * comp_inf, bool aggr )
{
	//Aggrgtd = aggr;
	// tamanho de ptr_x_til dependo do tipo de modelo
	// tamanho de cada componente de ptr_x_subp poderia ser definido antes, mas o numero de variáveis de cada subproblema só é obtido após criar os objetos subp.
	if ( Aggrgtd )
	{
		ptr_x_til.resize(1);
		ptr_x_til[0].resize(comp_inf[0] * sistema_a->GetNCenarios());		// vetor [X_0 X_1 ... X_Ncen]
		// comp_inf[i] = numero de variáveis de cada subproblema (cenário), igual para tds subp.
		// comp_inf[0] = n_sub
		//int n_cen = 0;
		//for (int i = 0; i < sistema_a->GetNCenarios(); i++)
		//	n_cen += comp_inf[i];
		//ptr_x_til[0].resize(n_cen);		// vetor [X_0 X_1 ... X_Ncen]
	}
	else
	{
		//vetorfloat2().swap(ptr_x_subp_agg);		// se o modelo for desagregado zerar vetor ptr_x_subp_agg, vetor n é usado
		ptr_x_til.resize(sistema_a->GetNCenarios());
		for (int i = 0; i < sistema_a->GetNCenarios(); i++)
		{
			ptr_x_til[i].resize(comp_inf[i]);
			ptr_x_subp[i].resize(comp_inf[i]);		// deve ser definido caso usa-se easy components, pois a soluçao é alocada e n copiada do subp.
		}
	}
}
void Scndec3Results::GravarSolucao(double obj_a, vetorfloat x_a, int status_a, int comp, double lb_a)
{
	// Grava o resultado de cada subproblema
	status = status_a;
	// wFi = 0 é o termo constante entao wFi = 1 corresponde à comp = 0;
	ptr_x_subp[comp] = x_a;		// talvez fique mais rápido receber em cada subproblema o endereço do vetorfloat de ptr_x_subp e preenche-lo (clear, resize e atribuir valor)
	
	// FiOracle exato
	//obj_subp[comp] = obj_a;
	
	// FiOracle inexato
	obj_subp[comp] = lb_a;
	epsilon_subp[comp] = obj_a - lb_a;
}
void Scndec3Results::GetSubGradiente(int wFi_a, double * SubG)
{
	//for (int i = 0; i < nd; i++)	// não precisa, pois é um vetor denso
	//	SubG[i] = 0;
	// O vetor de subgradientes tem a dimensao do vetor lambda, portanto a posicao do subg q sera preenchido depende do tipo de problema
	// Selecionar tipo do subproblema pelo valor componente
	int ndt = (nd / (sistema_a->GetNCenarios() - 1));		// numero de var. duais por cenário, diferente do numero de variáveis por subproblema!
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
	{
		// aqui o loop é somente para algumas variaveis de primeiro estágio
		//x = [pt u v] das restrições de não antecipatividade
		int deltax = 0;
		int deltasg = 0;
		if ((wFi_a - 1) == n_cen)		// subgradiente referente ao multiplicador da componente (wFi_a - 1)
		{
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
					SubG[i + deltasg + ndt*n_cen] = (1 - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) ) * ptr_x_subp[wFi_a - 1][i + deltax];
				deltax += sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
					SubG[i + deltasg + ndt*n_cen] = (1 - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) ) * ptr_x_subp[wFi_a - 1][i + deltax];
				deltax += 2*sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					deltax += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					deltax += sistema_a->termeletricasVtr.size();
				if ( sistema_a->GetFlagModeloRede() == 1)
					deltax += sistema_a->barrasVtr.size() - 1;
				deltax += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					SubG[r + deltasg + ndt*n_cen] = (1 - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) ) * ptr_x_subp[wFi_a - 1][r + deltax];
				deltax += 3*sistema_a->hidreletricasVtr.size();
				deltasg += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
					deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
				if ( sistema_a->GetFlagModeloRede() > 0)
					deltax += sistema_a->barrasVtr.size();
				else
					deltax += 1;
			}
		}
		else							// subgradiente referente aos multiplicadores das demais componentes
		{
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
					SubG[i + deltasg + ndt*n_cen] = - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) * ptr_x_subp[wFi_a - 1][i + deltax];
				deltax += sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
					SubG[i + deltasg + ndt*n_cen] = - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) * ptr_x_subp[wFi_a - 1][i + deltax];
				deltax += 2*sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					deltax += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					deltax += sistema_a->termeletricasVtr.size();
				if ( sistema_a->GetFlagModeloRede() == 1)
					deltax += sistema_a->barrasVtr.size() - 1;
				deltax += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					SubG[r + deltasg + ndt*n_cen] = - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) * ptr_x_subp[wFi_a - 1][r + deltax];
				deltax += 3*sistema_a->hidreletricasVtr.size();
				deltasg += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
					deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
				if ( sistema_a->GetFlagModeloRede() > 0)
					deltax += sistema_a->barrasVtr.size();
				else
					deltax += 1;
			}
		}
	}
	for (int i = 0; i < nd; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);
}
void Scndec3Results::GetSubGradiente(double * SubG)
{
	int ndt = (nd / (sistema_a->GetNCenarios() - 1));		// numero de var. duais por cenário, diferente do numero de variáveis por subproblema!
	vetorfloat valor_esperado;
	int deltax = 0;
	int deltasg = 0;
	for (int t = 0; t < sistema_a->GetTt1(); t++)
	{
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
		{
			valor_esperado.push_back(0);
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
				valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_subp[n_cen][i + deltax];
		}
		deltax += sistema_a->termeletricasVtr.size();
		deltasg += sistema_a->termeletricasVtr.size();
		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
		{
			valor_esperado.push_back(0);
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
				valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_subp[n_cen][i + deltax];
		}
		deltax += 2*sistema_a->termeletricasVtr.size();
		deltasg += sistema_a->termeletricasVtr.size();
		if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			deltax += sistema_a->termeletricasVtr.size();
		if (sistema_a->GetFlagInitAproxCT() > 1)	//F
			deltax += sistema_a->termeletricasVtr.size();
		if ( sistema_a->GetFlagModeloRede() == 1)
			deltax += sistema_a->barrasVtr.size() - 1;
		deltax += sistema_a->hidreletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
		{
			valor_esperado.push_back(0);
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
				valor_esperado[r + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_subp[n_cen][r + deltax];
		}
		deltax += 3*sistema_a->hidreletricasVtr.size();
		deltasg += sistema_a->hidreletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
			deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
		if ( sistema_a->GetFlagModeloRede() > 0)
			deltax += sistema_a->barrasVtr.size();
		else
			deltax += 1;
	}
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
	{
		deltax = 0;
		deltasg = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				SubG[i + deltasg + ndt*n_cen] = ptr_x_subp[n_cen][i + deltax] - valor_esperado[i + deltasg];
			deltax += sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
				SubG[i + deltasg + ndt*n_cen] = ptr_x_subp[n_cen][i + deltax] - valor_esperado[i + deltasg];
			deltax += 2*sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				deltax += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				deltax += sistema_a->termeletricasVtr.size();
			if ( sistema_a->GetFlagModeloRede() == 1)
				deltax += sistema_a->barrasVtr.size() - 1;
			deltax += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				SubG[r + deltasg + ndt*n_cen] = ptr_x_subp[n_cen][r + deltax] - valor_esperado[r + deltasg];
			deltax += 3*sistema_a->hidreletricasVtr.size();
			deltasg += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
				deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
			if ( sistema_a->GetFlagModeloRede() > 0)
				deltax += sistema_a->barrasVtr.size();
			else
				deltax += 1;
		}
	}
	for (int i = 0; i < nd; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);
}
void Scndec3Results::GetSubGradiente(int wFi_a, double * SubG, vetorfloat x_name)
{
	// Calcular subgradientes a partir de x_name!!
	int ndt = (nd / (sistema_a->GetNCenarios() - 1));		// numero de var. duais por cenário, diferente do numero de variáveis por subproblema!
	if (wFi_a > sistema_a->GetNCenarios())		// agregado
	{
		// transformar o vetor x_name agregado em um vetor com x desagregado por cenário, com a mesma estrutura que ptr_x_subp
		vetorfloat2 ptr_x_name_subp;
		ptr_x_name_subp.resize(sistema_a->GetNCenarios());
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			for (int i = 0; i < n_sub; i++)
				ptr_x_name_subp[n_cen].push_back(x_name[i + n_sub*n_cen]);
		//
		vetorfloat valor_esperado;
		int deltax = 0;
		int deltasg = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_name_subp[n_cen][i + deltax];
			}
			deltax += sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_name_subp[n_cen][i + deltax];
			}
			deltax += 2*sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				deltax += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				deltax += sistema_a->termeletricasVtr.size();
			if ( sistema_a->GetFlagModeloRede() == 1)
				deltax += sistema_a->barrasVtr.size() - 1;
			deltax += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[r + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_name_subp[n_cen][r + deltax];
			}
			deltax += 3*sistema_a->hidreletricasVtr.size();
			deltasg += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
				deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
			if ( sistema_a->GetFlagModeloRede() > 0)
				deltax += sistema_a->barrasVtr.size();
			else
				deltax += 1;
		}
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		{
			deltax = 0;
			deltasg = 0;
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
					SubG[i + deltasg + ndt*n_cen] = ptr_x_name_subp[n_cen][i + deltax] - valor_esperado[i + deltasg];
				deltax += sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
					SubG[i + deltasg + ndt*n_cen] = ptr_x_name_subp[n_cen][i + deltax] - valor_esperado[i + deltasg];
				deltax += 2*sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					deltax += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					deltax += sistema_a->termeletricasVtr.size();
				if ( sistema_a->GetFlagModeloRede() == 1)
					deltax += sistema_a->barrasVtr.size() - 1;
				deltax += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					SubG[r + deltasg + ndt*n_cen] = ptr_x_name_subp[n_cen][r + deltax] - valor_esperado[r + deltasg];
				deltax += 3*sistema_a->hidreletricasVtr.size();
				deltasg += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
					deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
				if ( sistema_a->GetFlagModeloRede() > 0)
					deltax += sistema_a->barrasVtr.size();
				else
					deltax += 1;
			}
		}
	}
	else		// desagregado
	{
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		{
			// aqui o loop é somente para algumas variaveis de primeiro estágio
			//x = [pt u d phg] das restrições de não antecipatividade
			int deltax = 0;
			int deltasg = 0;
			if ((wFi_a - 1) == n_cen)		// subgradiente referente ao multiplicador da componente (wFi_a - 1)
			{
				for (int t = 0; t < sistema_a->GetTt1(); t++)
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
						SubG[i + deltasg + ndt*n_cen] = (1 - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) ) * x_name[i + deltax];
					deltax += sistema_a->termeletricasVtr.size();
					deltasg += sistema_a->termeletricasVtr.size();
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
						SubG[i + deltasg + ndt*n_cen] = (1 - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) ) * x_name[i + deltax];
					deltax += 2*sistema_a->termeletricasVtr.size();
					deltasg += sistema_a->termeletricasVtr.size();
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
						deltax += sistema_a->termeletricasVtr.size();
					if (sistema_a->GetFlagInitAproxCT() > 1)	//F
						deltax += sistema_a->termeletricasVtr.size();
					if ( sistema_a->GetFlagModeloRede() == 1)
						deltax += sistema_a->barrasVtr.size() - 1;
					deltax += sistema_a->hidreletricasVtr.size();
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						SubG[r + deltasg + ndt*n_cen] = (1 - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) ) * x_name[r + deltax];
					deltax += 3*sistema_a->hidreletricasVtr.size();
					deltasg += sistema_a->hidreletricasVtr.size();
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
						deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
					if ( sistema_a->GetFlagModeloRede() > 0)
						deltax += sistema_a->barrasVtr.size();
					else
						deltax += 1;
				}
			}
			else							// subgradiente referente aos multiplicadores das demais componentes
			{
				for (int t = 0; t < sistema_a->GetTt1(); t++)
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
						SubG[i + deltasg + ndt*n_cen] = - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) * x_name[i + deltax];
					deltax += sistema_a->termeletricasVtr.size();
					deltasg += sistema_a->termeletricasVtr.size();
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
						SubG[i + deltasg + ndt*n_cen] = - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) * x_name[i + deltax];
					deltax += 2*sistema_a->termeletricasVtr.size();
					deltasg += sistema_a->termeletricasVtr.size();
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
						deltax += sistema_a->termeletricasVtr.size();
					if (sistema_a->GetFlagInitAproxCT() > 1)	//F
						deltax += sistema_a->termeletricasVtr.size();
					if ( sistema_a->GetFlagModeloRede() == 1)
						deltax += sistema_a->barrasVtr.size() - 1;
					deltax += sistema_a->hidreletricasVtr.size();
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						SubG[r + deltasg + ndt*n_cen] = - sistema_a->hidreletricasVtr[0].GetProbAfluencia(wFi_a - 1) * x_name[r + deltax];
					deltax += 3*sistema_a->hidreletricasVtr.size();
					deltasg += sistema_a->hidreletricasVtr.size();
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
						deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
					if ( sistema_a->GetFlagModeloRede() > 0)
						deltax += sistema_a->barrasVtr.size();
					else
						deltax += 1;
				}
			}
		}
	}
	for (int i = 0; i < nd; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);

}


double Scndec3Results::GetFobj()
{
	obj = 0;
	for (int sp = 0; sp < sistema_a->GetNCenarios(); sp++)
		obj += obj_subp[sp];
	return obj;
}
//
vetorfloat Scndec3Results::GetCompX(int comp)
{
	if ( Aggrgtd )
	{
		vetorfloat ptr_x_agg;
		ptr_x_agg.resize(n_sub * sistema_a->GetNCenarios());
		// alocar ptr_x_agg, junta o resultado de todos os subproblemas em um vetor só
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			for (int i = 0; i < n_sub; i++)
				ptr_x_agg[i + n_sub*n_cen] = ptr_x_subp[n_cen][i];

		return ptr_x_agg;
	}
	else
	{
		return ptr_x_subp[comp];
	}
}	
void Scndec3Results::AlocarX_med()
{
	// Alocar x_med a partir de ptr_x_subp
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nt1 = nt * sistema_a->GetTt1();		// numero de var. de primeiro estágio

	// variáveis de primeiro estágio (valor médio/esperado)
	// até ndt o valor das variáveis é dada pela média dos subproblemas
	for (int i = 0; i < nt1; i++)
	{
		x_med[i] = 0;
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			x_med[i] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_subp[n_cen][i];
	}
	// variáveis de segundo estágio
	// cada cenário tem sua própria variável, organizadas por ordem crescente no vetor x_med
	int n_var_2_est = (sistema_a->GetTt2() - sistema_a->GetTt1()) * nt + flag2*sistema_a->hidreletricasVtr.size();
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		for (int i = 0; i < n_var_2_est; i++)
			x_med[nt1 + i + n_var_2_est*n_cen] = ptr_x_subp[n_cen][nt1 + i];
}
vetorfloat Scndec3Results::GetX_med(bool til)
{
	// Alocar x_med a partir de ptr_x_subp ou ptr_x_til
	// A cada chamada do SetGiName o ptr_x_hat aqui nessa classe é atualizado, ou seja, ptr_x_hat é atualizado right after Fi() and all the GetGi()

	// hat_or_til indica a partir de qual solução é calculado o x_med
	// 0 a partir da soluçao da RL, ptr_x_hat (ptr_x_subp)
	// 1 a partir da soluçao convexifica, ptr_x_til
	vetorfloat * ptr_x_a;
	if ( til )
		ptr_x_a = &ptr_x_til[0];
	else
		ptr_x_a = &ptr_x_subp[0];

	// Alocar x_med a partir de ptr_x_subp
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nt1 = nt * sistema_a->GetTt1();		// numero de var. de primeiro estágio
	// variáveis de primeiro estágio (valor médio/esperado)
	// até ndt o valor das variáveis é dada pela média dos subproblemas
	// ndt = nt * sistema_a->GetTt1()
	for (int i = 0; i < nt1; i++)
	{
		x_med[i] = 0;
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			x_med[i] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_a[n_cen][i];
	}
	// variáveis de segundo estágio
	// cada cenário tem sua própria variável, organizadas por ordem crescente no vetor x_med
	int n_var_2_est = (sistema_a->GetTt2() - sistema_a->GetTt1()) * nt + flag2*sistema_a->hidreletricasVtr.size();
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		for (int i = 0; i < n_var_2_est; i++)
			x_med[nt1 + i + n_var_2_est*n_cen] = ptr_x_a[n_cen][nt1 + i];

	return x_med;
}
void Scndec3Results::AdicionarLBXtil(CMatrizEsparsa * x_spr, double mult)
{
	// Adicionar solução setada como LB no Bundle na combinação da solução convexificada

	// Para as variáveis replicadas (de primeiro estágio) x11 = x12 += x_spr * mult
	//x = [pt u cp F teta ph v d s (phmax) phg q z def] ou
	//x = [pt u up ud F teta ph v d s phmax phg q z def]	"modern" model
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;

	int nT1 = nt * sistema_a->GetTt1();		// numero de var. de primeiro estágio
	int nT2 = (sistema_a->GetTt2() - sistema_a->GetTt1()) * nt + flag2*sistema_a->hidreletricasVtr.size(); // numero de var. de segundo estágio

	if ( Aggrgtd )
	{
		if (mult >= 0.99999999)
		{
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			{
				for (int i = 0; i < nT1; i++)
					ptr_x_til[0][i + (nT1+nT2)*n_cen] = x_spr->GetElemento(i, 0) * mult;
				for (int i = nT1; i < nT2 + nT1; i++)
					ptr_x_til[0][i + (nT1+nT2)*n_cen] = x_spr->GetElemento(i + nT2*n_cen, 0) * mult;
			}
		}
		else
		{
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			{
				for (int i = 0; i < nT1; i++)
					ptr_x_til[0][i + (nT1+nT2)*n_cen] += x_spr->GetElemento(i, 0) * mult;
				for (int i = nT1; i < nT2 + nT1; i++)
					ptr_x_til[0][i + (nT1+nT2)*n_cen] += x_spr->GetElemento(i + nT2*n_cen, 0) * mult;
			}
		}
	}
	else
	{
		// subproblemas cen = 0,...,Ncen:
		//X_cen = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		//X_cen = [pt u up ud (F) (teta) ph v d s (phmax) phg q z (def) (vfol)] -> (.) tamanho dos vetores que dependem de flags
		
		if (mult >= 0.99999999)
		{
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			{
				for (int i = 0; i < nT1; i++)
					ptr_x_til[n_cen][i] = 0;
				for (int i = nT1; i < nT2 + nT1; i++)
					ptr_x_til[n_cen][i] = 0;
			}
		}
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		{
			for (int i = 0; i < nT1; i++)
				ptr_x_til[n_cen][i] += x_spr->GetElemento(i, 0) * mult;
			for (int i = nT1; i < nT2 + nT1; i++)
				ptr_x_til[n_cen][i] += x_spr->GetElemento(i + nT2*n_cen, 0) * mult;
		}
	}
}
void Scndec3Results::ExportarXmed(string nome_arquivo)
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
void Scndec3Results::ZerarSolucoes()
{
	x_med.clear();
	obj_subp.clear();
	ptr_x_subp.clear();
	ptr_x_til.clear();

	x_med.resize(n);
	obj = 0;
	obj_subp.resize(sistema_a->GetNCenarios());
	for (int i = 0; i <  sistema_a->GetNCenarios(); i++)
		obj_subp[i] = 0;
	ptr_x_subp.resize(sistema_a->GetNCenarios());
}
void Scndec3Results::CalcularNormaSubgXtil(vetorfloat &norma)
{
	// Um valor de norma para cada tipo de restrição relaxada, norma 1 e 2;
	// numero de tipos de restrições = nd / (sistema_a->GetNCenarios() - 1) * T1)
	// x = [pt u v]
	double dif = 0;
	vetorfloat valor_esperado;
	int ndt = (nd / (sistema_a->GetNCenarios() - 1));		// numero de var. duais por cenário, diferente do numero de variáveis por subproblema!
	int deltax = 0;
	int deltasg = 0;
	if ( Aggrgtd )
	{
		// calcula o valor esperado
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[0][i + deltax + n_sub*n_cen];
			}
			deltax += sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[0][i + deltax + n_sub*n_cen];
			}
			deltax += 2*sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				deltax += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				deltax += sistema_a->termeletricasVtr.size();
			if ( sistema_a->GetFlagModeloRede() == 1)
				deltax += sistema_a->barrasVtr.size() - 1;
			deltax += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[r + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[0][r + deltax + n_sub*n_cen];
			}
			deltax += 3*sistema_a->hidreletricasVtr.size();
			deltasg += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
				deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
			if ( sistema_a->GetFlagModeloRede() > 0)
				deltax += sistema_a->barrasVtr.size();
			else
				deltax += 1;
		}
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		{
			deltax = 0;
			deltasg = 0;
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				{
					dif = ptr_x_til[0][i + deltax + n_sub*n_cen] - valor_esperado[i + deltasg];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				deltax += sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
				{
					dif = ptr_x_til[0][i + deltax + n_sub*n_cen] - valor_esperado[i + deltasg];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				deltax += 2*sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					deltax += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					deltax += sistema_a->termeletricasVtr.size();
				if ( sistema_a->GetFlagModeloRede() == 1)
					deltax += sistema_a->barrasVtr.size() - 1;
				deltax += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				{
					dif = ptr_x_til[0][r + deltax + n_sub*n_cen] - valor_esperado[r + deltasg];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				deltax += 3*sistema_a->hidreletricasVtr.size();
				deltasg += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
					deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
				if ( sistema_a->GetFlagModeloRede() > 0)
					deltax += sistema_a->barrasVtr.size();
				else
					deltax += 1;
			}
		}
		norma[1] = sqrt(norma[1]);
	}
	else
	{
		// calcula o valor esperado
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[n_cen][i + deltax];
			}
			deltax += sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[i + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[n_cen][i + deltax];
			}
			deltax += 2*sistema_a->termeletricasVtr.size();
			deltasg += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				deltax += sistema_a->termeletricasVtr.size();
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				deltax += sistema_a->termeletricasVtr.size();
			if ( sistema_a->GetFlagModeloRede() == 1)
				deltax += sistema_a->barrasVtr.size() - 1;
			deltax += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
			{
				valor_esperado.push_back(0);
				for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
					valor_esperado[r + deltasg] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[n_cen][r + deltax];
			}
			deltax += 3*sistema_a->hidreletricasVtr.size();
			deltasg += sistema_a->hidreletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
				deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
			if ( sistema_a->GetFlagModeloRede() > 0)
				deltax += sistema_a->barrasVtr.size();
			else
				deltax += 1;
		}
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		{
			deltax = 0;
			deltasg = 0;
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				{
					dif = ptr_x_til[n_cen][i + deltax] - valor_esperado[i + deltasg];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				deltax += sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
				{
					dif = ptr_x_til[n_cen][i + deltax] - valor_esperado[i + deltasg];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				deltax += 2*sistema_a->termeletricasVtr.size();
				deltasg += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					deltax += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					deltax += sistema_a->termeletricasVtr.size();
				if ( sistema_a->GetFlagModeloRede() == 1)
					deltax += sistema_a->barrasVtr.size() - 1;
				deltax += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				{
					dif = ptr_x_til[n_cen][r + deltax] - valor_esperado[r + deltasg];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				deltax += 3*sistema_a->hidreletricasVtr.size();
				deltasg += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q, z
					deltax += 3*sistema_a->hidreletricasVtr[r].GetNGrupos();
				if ( sistema_a->GetFlagModeloRede() > 0)
					deltax += sistema_a->barrasVtr.size();
				else
					deltax += 1;
			}
		}
		norma[1] = sqrt(norma[1]);
	}
}
void Scndec3Results::CalcularNormaSubgXtilTV(vetorfloat &norma)
{
	// zerar elementos do vetor
	for (int i = 0; i < norma.size(); i++)
		norma[i] = 0;
	double dif = 0;
	vetorfloat valor_esperado;
	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	int nd_total = sistema_a->GetTt1() * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size());
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		nd_total -= sistema_a->GetTt1() * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		nd_total -= sistema_a->GetTt1() * (sistema_a->barrasVtr.size() - 1);
	int ndt = nd_total;		// numero de var. duais por cenário, diferente do numero de variáveis por subproblema!
	int delta = 0;
	if ( Aggrgtd )
	{
		// Nesse caso o vetor ptr_x_til tem somente um componente!!
		// As variáveis de cada cenário estão agrupadas em ordem nesse vetor
		for (int i = 0; i < ndt; i++)		// calcula o valor esperado para todas as variáveis (a posição das variáveis duplicadas é de 0 a ndt)
		{
			valor_esperado.push_back(0);
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
				valor_esperado[i] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[0][i + n_sub*n_cen];
		}
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		{
			for (int i = 0; i < ndt; i++)
			{
				dif = ptr_x_til[0][i + n_sub*n_cen] - valor_esperado[i];
				norma[0] += abs(dif);
				norma[1] += pow(dif, 2);
				if (norma[2] < abs(dif))
					norma[2] = abs(dif);
			}
		}
		norma[1] = sqrt(norma[1]);
	}
	else
	{
		// Um valor de norma para cada tipo de restrição relaxada, norma 1 e 2;
		// numero de tipos de restrições = nd / (sistema_a->GetNCenarios() - 1) * T1)
		// x = [pt u (up ud / cp) (F) (teta) ph v d s () phg q z (def)]
		for (int i = 0; i < ndt; i++)		// calcula o valor esperado para todas as variáveis
		{
			valor_esperado.push_back(0);
			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
				valor_esperado[i] += sistema_a->hidreletricasVtr[0].GetProbAfluencia(n_cen) * ptr_x_til[n_cen][i];
		}
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		{
			delta = 0;
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				{
					dif = ptr_x_til[n_cen][i + delta] - valor_esperado[i + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//u
				{
					dif = ptr_x_til[n_cen][i + delta] - valor_esperado[i + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->termeletricasVtr.size();
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	// up / cp
				{
					dif = ptr_x_til[n_cen][i + delta] - valor_esperado[i + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->termeletricasVtr.size();
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	// ud
					{
						dif = ptr_x_til[n_cen][i + delta] - valor_esperado[i + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->termeletricasVtr.size();
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	// ud
					{
						dif = ptr_x_til[n_cen][i + delta] - valor_esperado[i + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->termeletricasVtr.size();
				}
				if ( sistema_a->GetFlagModeloRede() == 1)
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)	// teta
					{
						dif = ptr_x_til[n_cen][b + delta] - valor_esperado[b + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->barrasVtr.size() - 1;
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
				{
					dif = ptr_x_til[n_cen][r + delta] - valor_esperado[r + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
				{
					dif = ptr_x_til[n_cen][r + delta] - valor_esperado[r + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
				{
					dif = ptr_x_til[n_cen][r + delta] - valor_esperado[r + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//s
				{
					dif = ptr_x_til[n_cen][r + delta] - valor_esperado[r + delta];
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
				}
				delta += sistema_a->hidreletricasVtr.size();
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phg
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
					{
						dif = ptr_x_til[n_cen][j + delta] - valor_esperado[j + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//q
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						dif = ptr_x_til[n_cen][j + delta] - valor_esperado[j + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//z
				{
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)	
					{
						dif = ptr_x_til[n_cen][j + delta] - valor_esperado[j + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->hidreletricasVtr[r].GetNGrupos();
				}
				if ( sistema_a->GetFlagModeloRede() > 0)
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)	// def
					{
						dif = ptr_x_til[n_cen][b + delta] - valor_esperado[b + delta];
						norma[0] += abs(dif);
						norma[1] += pow(dif, 2);
						if (norma[2] < abs(dif))
							norma[2] = abs(dif);
					}
					delta += sistema_a->barrasVtr.size();
				}
				else
				{
					dif = ptr_x_til[n_cen][delta] - valor_esperado[delta];	//def
					norma[0] += abs(dif);
					norma[1] += pow(dif, 2);
					if (norma[2] < abs(dif))
						norma[2] = abs(dif);
					delta += 1;
				}
			}
		}
		norma[1] = sqrt(norma[1]);
	}
}
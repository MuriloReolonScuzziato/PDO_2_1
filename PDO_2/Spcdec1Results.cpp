#include "Spcdec1Results.h"


Spcdec1Results::Spcdec1Results(CSistema * const sistema_end)
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
	N = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() );

	int JJ = 0;
	for (size_t r = 0; r < R; r++)
		JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
	n = N * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	na = N * (sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size());
	// na = nd nesta decomposição
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= N * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= N * (sistema_a->barrasVtr.size() - 1);
	x.resize(n);
	x_med.resize(n);
	xa.resize(na);
	obj = 0;
	obj_subp.resize(1 + I + N);
	for (int i = 0; i <  1 + I + N; i++)
		obj_subp[i] = 0;
	ptr_x_subp.resize(1 + I + N);
	ptr_x_subp_agg.resize(1);
	ptr_x_subp_agg[0].resize(n + na);
	
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
Spcdec1Results::~Spcdec1Results(void)
{
	delete inFile;
}

// SpcDec1

// Fazer outra funçao para guardar resultados dos subproblemas, por componente (guardar x no formato de cada subproblema)
// fazer tb outra funçao para calcular o subgradiente de acordo com o formato do x acima
// só colocar os dados nos formatos do x e xa (e x_med) originais no resultado final!!!

// subproblema H comp = 1;						// x = [ph v d s phg q z vfol] poderia ser um por cascata (ch) se fosse criada a var. phmax para dividir a restrição de reserva
// subproblema T comp = 1 + id1;				// id1 é o índice da termeletrica	//x = [pt u up ud F]
// subproblema D comp = 1 + I + id1;			// id1 é o índice do nó	//x = [pt teta ph def]

void Spcdec1Results::SetComponents( int * comp_inf, bool aggr )
{
	//Aggrgtd = aggr;
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
		ptr_x_til.resize(1 + I + N);
		for (int i = 0; i < 1 + I + N; i++)
		{
			ptr_x_til[i].resize(comp_inf[i]);
			ptr_x_subp[i].resize(comp_inf[i]);		// deve ser definido caso usa-se easy components, pois a soluçao é alocada e n copiada do subp.
		}
	}
}
void Spcdec1Results::GravarSolucao(double obj_a, vetorfloat x_a, int status_a, int comp)
{
	// Grava o resultado de cada subproblema
	// subproblema H comp = 1;
	// subproblema T comp = 1 + id1;				// id1 é o índice da termeletrica	//x = [pt u up ud F]
	// subproblema D comp = 1 + I + id1;			// id1 é o índice do nó	//x = [pt teta ph def]
	status = status_a;
	// wFi = 0 é o termo constante entao wFi = 1 corresponde à comp = 0;
	ptr_x_subp[comp] = x_a;		// talvez fique mais rápido receber em cada subproblema o endereço do vetorfloat de ptr_x_subp e preenche-lo (clear, resize e atribuir valor)
	obj_subp[comp] = obj_a;
}
void Spcdec1Results::GetSubGradiente(int wFi_a, double * SubG)
{
	for (int i = 0; i < na; i++)
		SubG[i] = 0;
	
	// O vetor de subgradientes tem a dimensao do vetor lambda, portanto a posicao do subg q sera preenchido depende do tipo de problema
	// Selecionar tipo do subproblema pelo valor componente: (Hidreletrico, Termeletrico, Demanda)
	// 1 <= wFi (comp) <= GetNrFi();
	// lambda [pt ph]
	int delta = 0;
	int deltaa = 0;
	int nt;
	int nat = (na / N);
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	if ( (wFi_a > 0) && (wFi_a <= 1) )			//Hidreletrico ( 1 <= comp <= 1 )	//x = [ph v d s phg q z (vfol)]
	{
		nt = (ptr_x_subp[0].size() - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;
		deltaa = sistema_a->termeletricasVtr.size();
		for (int t = 0; t < N; t++)
		{
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				SubG[r + deltaa] = ptr_x_subp[0][r + delta];						//ph
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size();
			else
				delta += nt;
			deltaa += nat;
		}
	}
	else if ( (wFi_a > 1) && (wFi_a <= 1 + I) )		//Termeletrico	//x = [pt u cp F] ou [pt u up ud F]
	{
		int i = wFi_a - 1 - 1;
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
	else if ( (wFi_a > 1 + I) && (wFi_a <= 1 + I + N) )		//Demanda	//x = [pta teta pha def]
	{
		int no = wFi_a - 1 - I - 1;
		int cen = 0;
		if (no >= sistema_a->GetTt1())		// identifica o cenário e o período de tempo do nó
			cen = (no - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
		deltaa = no * nat;

		for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			SubG[i + deltaa] = - ptr_x_subp[wFi_a - 1][i];		//pta
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - ptr_x_subp[wFi_a - 1][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)];			//pha
		}
	}

	for (int i = 0; i < na; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);
}
void Spcdec1Results::GetSubGradiente(double * SubG)
{
	int deltaa = 0;
	int delta = 0;
	int cen = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int nat = (na / N);
	int nt = (ptr_x_subp[0].size() - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;

	for (int t = 0; t < N; t++)
	{
		if (sistema_a->GetFlagTbinaryModel() == 0)
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				SubG[i + deltaa] = - ptr_x_subp[t + 1 + I][i] + ptr_x_subp[i + 1][t*(3 + flag4)];		// - pta + pt
		else
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				SubG[i + deltaa] = - ptr_x_subp[t + 1 + I][i] + ptr_x_subp[i + 1][t*(3 + flag4 + flag7)] + ptr_x_subp[i + 1][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
		deltaa += sistema_a->termeletricasVtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			SubG[r + deltaa] = - ptr_x_subp[t + 1 + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_subp[0][r + delta];		// - pha + ph
		}
		deltaa -= sistema_a->termeletricasVtr.size();
		deltaa += nat;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size();
		else
			delta += nt;
	}

	for (int i = 0; i < na; i++)
		SubG[i] *= sistema_a->GetPrecondidioner(i);
}
double Spcdec1Results::GetFobj()
{
	obj = 0;
	for (int sp = 0; sp < 1 + I + N; sp++)
		obj += obj_subp[sp];
	return obj;
}
vetorfloat Spcdec1Results::GetCompX(int comp)
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
vetorfloat Spcdec1Results::GetX_med()
{
	// Alocar x_med a partir de x e xa
	AlocarX_med();

	return x_med;
}
vetorfloat Spcdec1Results::GetX_med(bool til)
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
	//x = [pt u cp F teta ph v d s phg q z def]
	//xa = [pta pha]
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
					x_med[i + delta] = ((ptr_x_a[0][i + delta] + ptr_x_a[0][i + 1 + delta]*sistema_a->termeletricasVtr[i].GetPmin()) + ptr_x_a[0][n + i + deltaa]) / 2;
			}
			delta += (3+flag4+flag7)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
			deltaa += sistema_a->termeletricasVtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				x_med[r + delta] = (ptr_x_a[0][r + delta] + ptr_x_a[0][n + r + deltaa]) / 2;		//ph
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
		// subproblema H comp = 1;						// x = [ph v d s phg q z vfol] poderia ser um por cascata (ch) se fosse criada a var. phmax para dividir a restrição de reserva
		// subproblema T comp = 1 + id1;				// id1 é o índice da termeletrica	//x = [pt u up ud F]
		// subproblema D comp = 1 + I + id1;			// id1 é o índice do nó	//x = [pt teta ph def]
		int JJ = 0;
		for (int r = 0; r < R; r++)
			JJ += sistema_a->hidreletricasVtr[r].GetNGrupos();
		int jj;
		delta = 0;
		deltaa = 0;
		nat = (ptr_x_a[0].size() - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)
			{
				for (int i = 0; i < I; i++)
				{
					x_med[i + delta] = (ptr_x_a[t + 1 + I][i] + ptr_x_a[i + 1][t*(3 + flag4)]) / 2;				//pt
					x_med[i + I + delta] = ptr_x_a[i + 1][t*(3 + flag4) + 1];									//u		
					x_med[i + 2*I + delta] = ptr_x_a[i + 1][t*(3 + flag4) + 2];									//cp
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)
					for (int i = 0; i < I; i++)
						x_med[i + 3*I + delta] = ptr_x_a[i + 1][t*(3 + flag4) + 3];								//F
			}
			else
			{
				for (int i = 0; i < I; i++)
				{
					x_med[i + delta] = (max(ptr_x_a[t + 1 + I][i] - sistema_a->termeletricasVtr[i].GetPmin(), 0.0) + ptr_x_a[i + 1][t*(3 + flag4 + flag7)]) / 2;		//pt
					x_med[i + I + delta] = ptr_x_a[i + 1][t*(3 + flag4 + flag7) + 1];									//u		
					x_med[i + 2*I + delta] = ptr_x_a[i + 1][t*(3 + flag4 + flag7) + 2];									//up
					x_med[i + 3*I + delta] = ptr_x_a[i + 1][t*(3 + flag4 + flag7) + 3];									//ud
				}
				if (sistema_a->GetFlagInitAproxCT() > 1)
					for (int i = 0; i < I; i++)
						x_med[i + (3 + flag7)*I + delta] = ptr_x_a[i + 1][t*(3 + flag4 + flag7) + 3 + flag7];								//F
			}

			jj = 0;
			for (int r = 0; r < R; r++)
			{
				x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + delta] = (ptr_x_a[t + 1 + I][r + I + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_a[0][r + deltaa]) / 2;				//ph
				x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 1*R + delta] = ptr_x_a[0][r + R + deltaa];						//v
				x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*R + delta] = ptr_x_a[0][r + 2*R + deltaa];					//d
				x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 3*R + delta] = ptr_x_a[0][r + 3*R + deltaa];					//s
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					x_med[j + jj + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + delta] = ptr_x_a[0][j + jj + 4*R + deltaa];					//phg
					x_med[j + jj + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + JJ + delta] = ptr_x_a[0][j + jj + 4*R + JJ + deltaa];		//q
					x_med[j + jj + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + 2*JJ + delta] = ptr_x_a[0][j + jj + 4*R + 2*JJ + deltaa];	//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if (sistema_a->GetFlagModeloRede() == 1)
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size() - 1; b++)
					x_med[b + I*(3 + flag4 + flag7) + delta] = ptr_x_a[t + 1 + I][b + I];	//teta
			}

			if (sistema_a->GetFlagModeloRede() > 0)
			{
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					x_med[b + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + 3*JJ + delta] = ptr_x_a[t + 1 + I][b + I + flag1a*(sistema_a->barrasVtr.size() - 1) + R];	//def
			}
			else
				x_med[I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + 3*JJ + delta] = ptr_x_a[t + 1 + I][I + R];	//def
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			{
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
				deltaa += nat + flag2*sistema_a->hidreletricasVtr.size();
			}
			else
			{
				delta += nt;
				deltaa += nat;
			}
		}
		if (flag2 == 1)
		{
			delta = nt * (sistema_a->GetTt2() - 1);
			deltaa = nat * (sistema_a->GetTt2() - 1);
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					x_med[r + I*(3 + flag4 + flag7) + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*R + 3*JJ + flag1d*sistema_a->barrasVtr.size() + (1 - flag1d) + delta] = ptr_x_a[0][r + 4*R + 3*JJ + deltaa];		//vfol
				delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
				deltaa += nat * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
			}
		}
	}
	return x_med;
}

void Spcdec1Results::AlocarX_med()
{
	// Alocar x_med a partir de x e xa
	AlocarXeXa();
	// Para as variáveis duplicadas x_med = (x + xa) / 2
	//x = [pt u cp F teta ph v d s phg q z def]
	//xa = [pta pha]
	
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
void Spcdec1Results::AlocarPtrXHat()
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
	for (int comp = 0; comp < 1; comp++)	//Hidreletrico
	{
		delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;
		for (int t = 0; t < N; t++)
		{
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				ptr_x_subp_agg[0][r + delta] = ptr_x_subp[comp][r + delta_subP];																					//ph
				ptr_x_subp_agg[0][r + 1*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][r + 1*sistema_a->hidreletricasVtr.size() + delta_subP];		//v
				ptr_x_subp_agg[0][r + 2*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][r + 2*sistema_a->hidreletricasVtr.size() + delta_subP];		//d
				ptr_x_subp_agg[0][r + 3*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][r + 3*sistema_a->hidreletricasVtr.size() + delta_subP];		//s
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					ptr_x_subp_agg[0][j + 5*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][j + 5*sistema_a->hidreletricasVtr.size() + delta_subP];					//phg
					ptr_x_subp_agg[0][j + 5*sistema_a->hidreletricasVtr.size() + JJ + delta] = ptr_x_subp[comp][j + JJ + 5*sistema_a->hidreletricasVtr.size() + delta_subP];		//q
					ptr_x_subp_agg[0][j + 5*sistema_a->hidreletricasVtr.size() + 2*JJ + delta] = ptr_x_subp[comp][j + 2*JJ + 5*sistema_a->hidreletricasVtr.size() + delta_subP];	//z
				}
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			{
				delta_subP += nt_subP + flag2*sistema_a->hidreletricasVtr.size() ;
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			}
			else
			{
				delta_subP += nt_subP;
				delta += nt;
			}
		}
		if (flag2 == 1)
		{
			delta = nt * sistema_a->GetTt2();
			delta_subP = nt_subP * sistema_a->GetTt2();
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					x[delta + r] = ptr_x_subp[comp][delta_subP + r];		//vfol
				delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
				delta_subP += nt_subP * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
			}
		}
	}
	for (int comp = 1; comp < 1 + I; comp++)		//Termeletrico
	{
		// id1 é o índice da termeletrica e id2 é nada
		int id1 = comp - 1;
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
	for (int comp = 1 + I; comp < 1 + I + N; comp++)		//Demanda
	{
		// id1 é o índice do nó
		int id1 = comp - 1 - I;
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
				ptr_x_subp_agg[0][b + flag1a*(sistema_a->barrasVtr.size() - 1) + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];	//def
		}
		else
			ptr_x_subp_agg[0][flag1a*(sistema_a->barrasVtr.size() - 1) + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + sistema_a->hidreletricasVtr.size()];	//def
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			ptr_x_subp_agg[0][n + r + sistema_a->termeletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)];		//ph
		}
	}
}
void Spcdec1Results::AlocarXeXa()
{
	// Alocar x e xa a partir de ptr_x_subp
	int JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();

	int delta, deltaa, delta_subP;
	int nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N);
	int nat = (na / N);
	int nt_subP;
	int jj;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	for (int comp = 0; comp < 1; comp++)	//Hidreletrico
	{
		delta = (3 + flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		delta_subP = 0;
		nt_subP = (ptr_x_subp[comp].size() - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;
		// Alocar resultados do subproblema
		// xa e x
		for (int t = 0; t < N; t++)
		{
			jj = 0;
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
			{
				x[r + delta] = ptr_x_subp[comp][r + delta_subP];																					//ph
				x[r + 1*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][r + 1*sistema_a->hidreletricasVtr.size() + delta_subP];		//v
				x[r + 2*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][r + 2*sistema_a->hidreletricasVtr.size() + delta_subP];		//d
				x[r + 3*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][r + 3*sistema_a->hidreletricasVtr.size() + delta_subP];		//s
				for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
				{
					x[j + jj + 4*sistema_a->hidreletricasVtr.size() + delta] = ptr_x_subp[comp][j + jj + 4*sistema_a->hidreletricasVtr.size() + delta_subP];					//phg
					x[j + jj + 4*sistema_a->hidreletricasVtr.size() + JJ + delta] = ptr_x_subp[comp][j + jj + JJ + 4*sistema_a->hidreletricasVtr.size() + delta_subP];		//q
					x[j + jj + 4*sistema_a->hidreletricasVtr.size() + 2*JJ + delta] = ptr_x_subp[comp][j + jj + 2*JJ + 4*sistema_a->hidreletricasVtr.size() + delta_subP];	//z
				}
				jj += sistema_a->hidreletricasVtr[r].GetNGrupos();
			}
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			{
				delta_subP += nt_subP + flag2*sistema_a->hidreletricasVtr.size() ;
				delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
			}
			else
			{
				delta_subP += nt_subP;
				delta += nt;
			}
		}
		if (flag2 == 1)
		{
			delta = nt * sistema_a->GetTt2();
			delta_subP = nt_subP * sistema_a->GetTt2();
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					x[delta + r] = ptr_x_subp[comp][delta_subP + r];		//vfol
				delta += nt * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
				delta_subP += nt_subP * (sistema_a->GetTt2() - sistema_a->GetTt1()) + sistema_a->hidreletricasVtr.size();
			}
		}
	}
	for (int comp = 1; comp < 1 + I; comp++)		//Termeletrico
	{
		// id1 é o índice da termeletrica e id2 é nada
		int id1 = comp - 1;
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
	for (int comp = 1 + I; comp < 1 + I + N; comp++)		//Demanda
	{
		// id1 é o índice do nó
		int id1 = comp - 1 - I;
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
				x[b + flag1a*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][b + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		}
		else
			x[flag1a*(sistema_a->barrasVtr.size() - 1) + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + delta] = ptr_x_subp[comp][sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1) + 2*sistema_a->hidreletricasVtr.size()];	//def
		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			xa[r + sistema_a->termeletricasVtr.size() + deltaa] = ptr_x_subp[comp][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)];		//ph
		}
	}
}

void Spcdec1Results::ExportarX(string nome_arquivo)
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
void Spcdec1Results::ExportarXA(string nome_arquivo)
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
void Spcdec1Results::ExportarXmed(string nome_arquivo)
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

void Spcdec1Results::ZerarSolucoes()
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
	obj_subp.resize(1 + I + N);
	for (int i = 0; i <  1 + I + N; i++)
		obj_subp[i] = 0;
	ptr_x_subp.resize(1 + I + N);
	ptr_x_subp_agg.resize(1);
	ptr_x_subp_agg[0].resize(n + na);
}
void Spcdec1Results::CalcularNormaSubgXtil(vetorfloat &norma1, vetorfloat &norma2)
{
	// Um valor de norma para cada tipo de restrição relaxada, norma 1 e 2;
	int deltaa = 0;
	int delta = 0;
	int cen = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	double dif = 0;
	int nat = (na / N);
	int nt = (ptr_x_subp[0].size() - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / N;

	for (int t = 0; t < N; t++)
	{
		if (sistema_a->GetFlagTbinaryModel() == 0)
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				dif = - ptr_x_til[t + 1 + I][i] + ptr_x_til[i + 1][t*(3 + flag4)];		// - pta + pt
				norma1[0] += abs(dif);
				norma2[0] += pow(dif, 2);
			}
		else
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				dif = - ptr_x_til[t + 1 + I][i] + ptr_x_til[i + 1][t*(3 + flag4 + flag7)] + ptr_x_til[i + 1][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
				norma1[0] += abs(dif);
				norma2[0] += pow(dif, 2);
			}
		deltaa += sistema_a->termeletricasVtr.size();

		for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
		{
			dif = - ptr_x_til[t + 1 + I][r + sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_til[0][r + delta];		// - pha + ph
			norma1[1] += abs(dif);
			norma2[1] += pow(dif, 2);
		}
		deltaa -= sistema_a->termeletricasVtr.size();
		deltaa += nat;
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size();
		else
			delta += nt;
	}
	for (size_t i = 0; i < norma2.size(); i++)
	{
		norma2[i] = sqrt(norma2[i]);
	}
}
void Spcdec1Results::CalcularNormaSubgXtil(vetorfloat &normas)
{
	// Um valor de norma para cada tipo de restrição relaxada, norma 1, 2 e Inf;
	int delta = 0;
	int cen = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int nt;
	int nat = (na / N);		// nd = na nesta decomposição
	double dif = 0;
	// zerar elementos do vetor
	for (int i = 0; i < normas.size(); i++)
		normas[i] = 0;

	// lambda [pt ph (phmax || res)]
	// x = [pt u cp (F) (teta) ph v d s (phmax) phg q z (def) (vfol)]

	if ( Aggrgtd )
	{
		// aqui ptr_x_til tem a mesma estrutura do ptr_x_subp_agg
		cout << " Spcdec2Results::CalcularNormaSubgXtil : calculo da norma para o modelo agregado não verificada!" << endl;
		
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
			}
			delta_sg -= I;
			delta -= (3+flag4+flag7)*I + flag1a*(sistema_a->barrasVtr.size() - 1);
			
			if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
				delta += nt + flag2*sistema_a->hidreletricasVtr.size();
			else
				delta += nt;
			delta_sg += I + R;		// delta para o numero de variaveis duplicadas xa
		}
		normas[1] = sqrt(normas[1]);
	}
	else
	{
		for (int t = 0; t < N; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 0)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					dif = - ptr_x_til[t + 1 + I][i] + ptr_x_til[i + 1][t*(3 + flag4)];		// - pta + pt
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			}
			else
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				{
					dif = - ptr_x_til[t + 1 + I][i] + ptr_x_til[i + 1][t*(3 + flag4 + flag7)] + ptr_x_til[i + 1][t*(3 + flag4 + flag7) + 1]*sistema_a->termeletricasVtr[i].GetPmin();		//- pta + pt + u*pt_min
					normas[0] += abs(dif);
					normas[1] += pow(dif, 2);
					if (normas[2] < abs(dif))
						normas[2] = abs(dif);
				}
			}
			nt = (ptr_x_subp[0].size() - flag2*sistema_a->GetNCenarios()*R) / N;
			if ( t >= sistema_a->GetTt2() )
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				delta = t*nt + flag2*cen*R;
			}	
			else
				delta = t*nt;
			for (size_t r = 0; r < R; r++)
			{
				dif = - ptr_x_til[t + 1 + I][r + I + flag1a*(sistema_a->barrasVtr.size() - 1)] + ptr_x_til[0][r + delta];		// - pha + ph
				normas[0] += abs(dif);
				normas[1] += pow(dif, 2);
				if (normas[2] < abs(dif))
					normas[2] = abs(dif);
			}
		}
		normas[1] = sqrt(normas[1]);
	}
}
void Spcdec1Results::AlocarXmed(CMatrizEsparsa * x_spr)
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
			x_spr->SubstituirElemento(r + delta, 0, (x[r + delta] + xa[r + deltaa]) / 2);		//ph
		delta -= (3+flag7+flag4)*sistema_a->termeletricasVtr.size() + flag1a*(sistema_a->barrasVtr.size() - 1);
		deltaa -= sistema_a->termeletricasVtr.size();
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			delta += nt + flag2*sistema_a->hidreletricasVtr.size() ;
		else
			delta += nt;
		deltaa += nat;
	}
}
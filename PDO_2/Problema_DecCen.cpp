#include "Problema_DecCen.h"

using namespace met_DecCen;

CProblema_DecCen::CProblema_DecCen(CSistema * const sistema_end, CResultados * const resultadosGurobi_end) :ambGRB(GRBEnv())
{
	sistema_a = sistema_end;
	resultadosGurobi = resultadosGurobi_end;
	fdual_RL = 0;
	fdual_RP = 0;
	norma_sg = 0;
	norma_g	= 0;

	inicio_cenarios.resize(sistema_a->GetNCenarios()); fim_cenarios.resize(sistema_a->GetNCenarios());
	tempo_cenarios.resize(sistema_a->GetNCenarios());
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		tempo_cenarios[cen] = 0;
	tempoRL = 0; tempoModelo = 0; tempoResolucao = 0; tempoGravar_Escrever = 0;

	JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasPtr.size(); i++)
		JJ += sistema_a->hidreletricasPtr[i]->GetNGrupos();
	int n_var_cen;
	flag1 = int (1 - sistema_a->GetFlagBarraUnica());
	flag4;
	if (sistema_a->GetFlagAproxCustoT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	n_var_cen = sistema_a->GetTt2() * ((3+flag4)*sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 5*(sistema_a->hidreletricasPtr.size()) + 3*JJ + flag1*(sistema_a->barrasPtr.size()) + (1 - flag1)) + sistema_a->hidreletricasPtr.size();
	n_var_est1 = (n_var_cen - sistema_a->hidreletricasPtr.size())/sistema_a->GetTt2()*sistema_a->GetTt1();
	x.resize(sistema_a->GetNCenarios());
	for (int i = 0; i < sistema_a->GetNCenarios(); i++)
		x[i].resize(n_var_cen);		// o tamanho de x é : numero de cenário x (numero variaveis em um cenario (ate T2))
	// variáveis com restrição de não antecipatividade pt, d, phg. (x = [pt u cp F teta ph v d s phmax phg q z def vfol])
	n_var_duplicadas = int(sistema_a->GetTt1() * (sistema_a->termeletricasPtr.size() + sistema_a->hidreletricasPtr.size() + JJ));	// numero de variáveis de não antecipatividade para um cenário
	Lambda.resize(sistema_a->GetNCenarios());	// o tamanho de Lambda é : numero de cenário x (numero variaveis de não antecipatividade em um cenario (ate T1), ou seja, do estagio 1)
	//mi.resize(sistema_a->GetNCenarios());
	mi.resize(n_var_duplicadas);	// o tamanho de mi é o numero variaveis de não antecipatividade em um cenario (ate T1)
	x_med.resize(n_var_duplicadas);	// o tamanho de x_med é o numero variaveis de não antecipatividade em um cenario (ate T1)
	Lambda_med.resize(n_var_duplicadas);
	subgrad.resize(sistema_a->GetNCenarios());
	grad.resize(sistema_a->GetNCenarios());
	for (int i = 0; i < sistema_a->GetNCenarios(); i++)
	{
		Lambda[i].resize(n_var_duplicadas);
		subgrad[i].resize(n_var_duplicadas);
		grad[i].resize(n_var_duplicadas);
		//mi[i].resize(n_var_duplicadas);
	}

	//bundle = new CBundle(ambGRB, sistema_a->GetNCenarios() * n_var_duplicadas);
	//bundle->oracle.bind(this, &CProblema_DecCen::oracle);
}

CProblema_DecCen::~CProblema_DecCen(void)
{
}

void CProblema_DecCen::oracle(const vetorfloat2 &lambda, double &fdual, vetorfloat2 &subgradiente)
{
	// Resoluçao dos problemas de cada cenario
	// ------------------------------------------------
	fdual = 0;
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		QueryPerformanceCounter(&inicio_cenarios[cen]);
		//
		problemaCenario[cen]->ResolverProblemaRL(resultadosGurobi, &lambda[cen], it_RL);
		//
		QueryPerformanceCounter(&fim_cenarios[cen]);
		tempo_cenarios[cen] += double(fim_cenarios[cen].QuadPart - inicio_cenarios[cen].QuadPart)/double(clockPerSec.QuadPart);
		x[cen] = resultadosGurobi->GetX();
		fdual += resultadosGurobi->GetFobj();
	}
	// ------------------------------------------------

	// Resolução do subproblema de coordenação
	// ------------------------------------------------
	double Ld_med, Lphg_med, Lpt_med;
	int delta_x = 0;
	int delta_subgrad = 0;
	double fdual_C = 0;
	for (int t = 0; t < sistema_a->GetTt1(); t++)		
	{
		for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
		{
			Lpt_med = 0;
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				Lpt_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(x[cen][delta_x + i]);
			//Lpt_med = Lpt_med/sistema_a->GetNCenarios();
			if ( Lpt_med >= 0 )
				x_med[delta_subgrad + i] = sistema_a->termeletricasPtr[i]->GetPmax();		// pt_med
			else
				x_med[delta_subgrad + i] = 0;		// pt_med
			fdual_C += - Lpt_med * x_med[delta_subgrad + i];
		}
		delta_subgrad += sistema_a->termeletricasPtr.size();
		delta_x += (3+flag4) * sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 2*sistema_a->hidreletricasPtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			Ld_med = 0;
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				Ld_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(x[cen][delta_x + r]);
			if ( Ld_med >= 0 )
				x_med[delta_subgrad + r] = sistema_a->hidreletricasPtr[r]->GetDmax();		// d_med
			else
				x_med[delta_subgrad + r] = 0;		// d_med
			fdual_C += - Ld_med * x_med[delta_subgrad + r];
		}
		delta_subgrad += sistema_a->hidreletricasPtr.size();
		delta_x += 3*sistema_a->hidreletricasPtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				Lphg_med = 0;
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					Lphg_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(x[cen][delta_x + j]);
				if ( Lphg_med >= 0 )
					x_med[delta_subgrad + j] = sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetPmax() * sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades();		// phg_med
				else
					x_med[delta_subgrad + j] = 0;		// phg_med
				fdual_C += - Lphg_med * x_med[delta_subgrad + j];
			}
			delta_subgrad += sistema_a->hidreletricasPtr[r]->GetNGrupos();
			delta_x += sistema_a->hidreletricasPtr[r]->GetNGrupos();
		}
		delta_x += 2*JJ + flag1*(sistema_a->barrasPtr.size()) + (1 - flag1);
	}
	fdual += fdual_C;
	// ------------------------------------------------

	// Calculo do Subgradiente
	// ------------------------------------------------
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		delta_subgrad = 0;
		delta_x = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
				subgradiente[cen][delta_subgrad + i] = x[cen][delta_x + i] - x_med[delta_subgrad + i];
			delta_subgrad += sistema_a->termeletricasPtr.size();
			delta_x += (3+flag4) * sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 2*sistema_a->hidreletricasPtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
				subgradiente[cen][delta_subgrad + r] = x[cen][delta_x + r] - x_med[delta_subgrad + r];
			delta_subgrad += sistema_a->hidreletricasPtr.size();
			delta_x += 3*sistema_a->hidreletricasPtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
			{
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
					subgradiente[cen][delta_subgrad + j] = x[cen][delta_x + j] - x_med[delta_subgrad + j];
				delta_subgrad += sistema_a->hidreletricasPtr[r]->GetNGrupos();
				delta_x += sistema_a->hidreletricasPtr[r]->GetNGrupos();
			}
			delta_x += 2*JJ + flag1*(sistema_a->barrasPtr.size()) + (1 - flag1);
		}
	}
	// ------------------------------------------------
}
void CProblema_DecCen::oracle(const VectorXd &lambda_e, double &fdual, VectorXd &subgradiente_e)
{
	// Converter tipo Lambda
	// ------------------------------------------------
	// Usar Lambda e subgrad (atribuitos da classe)
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		for (int i = 0; i < n_var_duplicadas; i++)
			Lambda[cen][i] = lambda_e(cen*n_var_duplicadas + i);
	// ------------------------------------------------

	// Calculo do valor esperado de lambda (Lambda_med)
	// ------------------------------------------------
	double Ld_med, Lphg_med, Lpt_med;
	int delta_subgrad = 0;
	double fdual_C = 0;
	for (int t = 0; t < sistema_a->GetTt1(); t++)		
	{
		for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
		{
			Lpt_med = 0;
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				Lpt_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(Lambda[cen][delta_subgrad + i]);
			Lambda_med[delta_subgrad + i] = Lpt_med;		// Lpt_med
		}
		delta_subgrad += sistema_a->termeletricasPtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			Ld_med = 0;
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				Ld_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(Lambda[cen][delta_subgrad + r]);
			Lambda_med[delta_subgrad + r] = Ld_med;		// Ld_med
		}
		delta_subgrad += sistema_a->hidreletricasPtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				Lphg_med = 0;
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					Lphg_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(Lambda[cen][delta_subgrad + j]);
				Lambda_med[delta_subgrad + j] = Lphg_med;		// Lphg_med
			}
			delta_subgrad += sistema_a->hidreletricasPtr[r]->GetNGrupos();
		}
	}
	// Atualizar: Lambda = Lambda - Lambda_med
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		for (int i = 0; i < n_var_duplicadas; i++)
			Lambda[cen][i] = Lambda[cen][i] - Lambda_med[i];
	// ------------------------------------------------
	
	// Resoluçao dos problemas de cada cenario
	// ------------------------------------------------
	fdual_RL = 0;
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		QueryPerformanceCounter(&inicio_cenarios[cen]);
		//
		problemaCenario[cen]->ResolverProblemaRL(resultadosGurobi, &Lambda[cen], it_RL);
		//
		QueryPerformanceCounter(&fim_cenarios[cen]);
		tempo_cenarios[cen] += double(fim_cenarios[cen].QuadPart - inicio_cenarios[cen].QuadPart)/double(clockPerSec.QuadPart);
		x[cen] = resultadosGurobi->GetX();
		fdual_RL += resultadosGurobi->GetFobj();
	}
	// ------------------------------------------------

	// Calculo do valor esperado de x (x_med)
	// ------------------------------------------------
	double d_med, phg_med, pt_med;
	int delta_x = 0;
	delta_subgrad = 0;
	for (int t = 0; t < sistema_a->GetTt1(); t++)		
	{
		for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
		{
			pt_med = 0;
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				pt_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(x[cen][delta_x + i]);
			x_med[delta_subgrad + i] = pt_med;		// pt_med
		}
		delta_subgrad += sistema_a->termeletricasPtr.size();
		delta_x += (3+flag4) * sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 2*sistema_a->hidreletricasPtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			d_med = 0;
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				d_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(x[cen][delta_x + r]);
			x_med[delta_subgrad + r] = d_med;		// d_med
		}
		delta_subgrad += sistema_a->hidreletricasPtr.size();
		delta_x += 3*sistema_a->hidreletricasPtr.size();
		for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
		{
			for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
			{
				phg_med = 0;
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					phg_med += sistema_a->hidreletricasPtr[0]->GetProbAfluencia(cen)*(x[cen][delta_x + j]);
				x_med[delta_subgrad + j] = phg_med;		// phg_med
					
			}
			delta_subgrad += sistema_a->hidreletricasPtr[r]->GetNGrupos();
			delta_x += sistema_a->hidreletricasPtr[r]->GetNGrupos();
		}
		delta_x += 2*JJ + flag1*(sistema_a->barrasPtr.size()) + (1 - flag1);
	}
	// ------------------------------------------------

	// Calculo do Subgradiente
	// ------------------------------------------------
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	{
		delta_subgrad = 0;
		delta_x = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
				subgrad[cen][delta_subgrad + i] = x[cen][delta_x + i] - x_med[delta_subgrad + i];
			delta_subgrad += sistema_a->termeletricasPtr.size();
			delta_x += (3+flag4) * sistema_a->termeletricasPtr.size() + flag1*(sistema_a->barrasPtr.size() - 1) + 2*sistema_a->hidreletricasPtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
				subgrad[cen][delta_subgrad + r] = x[cen][delta_x + r] - x_med[delta_subgrad + r];
			delta_subgrad += sistema_a->hidreletricasPtr.size();
			delta_x += 3*sistema_a->hidreletricasPtr.size();
			for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
			{
				for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
					subgrad[cen][delta_subgrad + j] = x[cen][delta_x + j] - x_med[delta_subgrad + j];
				delta_subgrad += sistema_a->hidreletricasPtr[r]->GetNGrupos();
				delta_x += sistema_a->hidreletricasPtr[r]->GetNGrupos();
			}
			delta_x += 2*JJ + flag1*(sistema_a->barrasPtr.size()) + (1 - flag1);
		}
	}
	// ------------------------------------------------

	// Converter tipo subgradiente
	// ------------------------------------------------
	for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		for (int i = 0; i < n_var_duplicadas; i++)
			subgradiente_e(cen*n_var_duplicadas + i) = - subgrad[cen][i];
	// ------------------------------------------------
	fdual = - fdual_RL;
}
void CProblema_DecCen::subgradiente(int &itmax)
{
	int n_var_dup_pper = n_var_duplicadas/sistema_a->GetTt1();
	int n_var_pper = n_var_est1/sistema_a->GetTt1();
	//int delta_x = 0;
	//int delta_subgrad = 0;

	//// Calculo da tolerancia da norma do gradiente, depende do numero e do valor maximo das variaveis relaxadas
	//// ------------------------------------------------
	//double tol_sg;
	double media_movel;
	//double dhmax, phgmax;
	vetorfloat media_fo;
	media_fo.resize(5);		// media movel dos ultimos 5 valores da função objetivo
	for (size_t i = 0; i < media_fo.size(); i++)
		media_fo[i] = 0;
	//tol_sg = 0;	// Relaxando as variaveis pt, d e phg
	//double erro_rel = 0.001;	// valor de diferença relativa máxima permitida por variável (0.001 = 0.1% do máximo valor)
	//for (size_t i = 0; i < sistema_a->termeletricasPtr.size(); i++)
	//	tol_sg += double (pow(erro_rel*sistema_a->termeletricasPtr[i]->GetPmax(),2));
	//for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
	//{
	//	dhmax = sistema_a->hidreletricasPtr[r]->GetSmax();
	//	for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
	//		dhmax = dhmax + sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetQmax() * sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades();
	//tol_sg += double (pow(erro_rel*dhmax,2));
	//}
	//for (size_t r = 0; r < sistema_a->hidreletricasPtr.size(); r++)
	//	for (int j = 0; j < sistema_a->hidreletricasPtr[r]->GetNGrupos(); j++)
	//		tol_sg += double (pow(erro_rel*sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetPmax() * sistema_a->hidreletricasPtr[r]->grupoPtr[j]->GetNUnidades(),2));
	//tol_sg = sqrt(tol_sg * sistema_a->GetTt1() * sistema_a->GetNCenarios());
	//// ------------------------------------------------

	// Processo iterativo do método do subgradiente
	// ------------------------------------------------
	for (it_RL = 0; it_RL < itmax; it_RL++)
	{
		// Resoluçao dos problemas de cada cenario
		// ------------------------------------------------
		//vetorfloat lambda, subgradiente;		// converter tamanho dos vetores Lambda e subgradiente
		//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		//{
		//	for (int i = 0; i < n_var_duplicadas; i++)
		//	{
		//		lambda.push_back(Lambda[cen][i]);
		//		subgradiente.push_back(subgrad[cen][i]);
		//	}
		//}
		QueryPerformanceCounter(&inicio_Resol);
		oracle(Lambda, fdual_RL, subgrad);
		QueryPerformanceCounter(&fim_Resol);
		//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		//{
		//	for (int i = 0; i < n_var_duplicadas; i++)
		//	{
		//		Lambda[cen][i] = lambda[i + cen*n_var_duplicadas];
		//		subgrad[cen][i] = subgradiente[i + cen*n_var_duplicadas];
		//	}
		//}
		tempoResolucao += double(fim_Resol.QuadPart - inicio_Resol.QuadPart)/double(clockPerSec.QuadPart);
		// ------------------------------------------------
		
		// Calcular critério de parada (norma do subgradiente e media movel)
		// ------------------------------------------------
		norma_sg = 0;
		for (size_t i = 0; i < subgrad.size(); i++)
			for (size_t ii = 0; ii < subgrad[i].size(); ii++)
				norma_sg += pow(subgrad[i][ii],2);
		norma_sg = sqrt(norma_sg);
		//
		media_movel = 0;
		for (size_t i = 0; i < media_fo.size(); i++)
			media_movel += media_fo[i];
		media_movel = media_movel/media_fo.size();
		media_fo[it_RL%5] = fdual_RL;		// Inserir valor da f.o. no vetor de f.o's
		// ------------------------------------------------

		resultadosGurobi->CarregarResIterativoDecCen(fdual_RL, norma_sg, &x, &Lambda, &subgrad, &x_med, &mi, n_var_pper);
		cout << it_RL << " : " << norma_sg << endl;

		// Verificar convergencia
		// ------------------------------------------------
		//if ((norma_sg <= tol_sg) && (abs(media_movel - fdual_RL) <= fdual_RL*0.001))		// Ajuste do 0.0001
		if ((it_RL == itmax - 1) || (abs(media_movel - fdual_RL) <= fdual_RL*0.0001))		// Ajuste do 0.0001
		{	
			//resultadosGurobi->CarregarResultadosProcessoIterativo(fdual_RL, norma_sg, &x, &Lambda, &subgrad, &x_med, &mi, sistema_a, n_var_pper);
			//resultadosGurobi->EscreverTelaProcessoIterativo(1);
			cout << "Convergiu!" << endl;
			break;
		}
		// ------------------------------------------------

		// Atualizar Lambda
		// ------------------------------------------------
		double alfa;
		if (it_RL == 0)
			alfa = 100/norma_sg;
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			for (int i = 0; i < n_var_duplicadas; i++)
			{
				Lambda[cen][i] += alfa*(double(1)/(it_RL+1))*subgrad[cen][i];
				//Lambda[cen][i] += subgrad[cen][i];					// fica oscilando e não converge
				//Lambda[cen][i] += 0.1*subgrad[cen][i]/norma_sg;		// Aplicar armijo para encontrar um tamanho de passo!!! neste caso alfa = 0.1
				//Lambda[cen][i] += subgrad[cen][i]/norma_sg;		// o valor da norma é muito grande, assim o algoritmo n converge
			}
		// ------------------------------------------------
	}
	// ------------------------------------------------
}

void CProblema_DecCen::ResolverRL()
{
	//// Parametro e valores iniciais
	//// ------------------------------------------------
	//int itmax = 100;			// maximo de iteraçoes
	////
	//bool importar_lambda;
	//double lambda_inicial;
	//if (sistema_a->GetFlagVarBin() == false)
	//{
	//	importar_lambda = false;		// importar lambda inical (true/false)
	//	lambda_inicial = 0.0;		// lambda inical
	//}
	//else
	//{
	//	importar_lambda = true;
	//	lambda_inicial = 0.0;
	//}
	//// ------------------------------------------------
	//bundle->npmax = itmax;
	//bundle->itmax = itmax;
	//bundle->setTol(10);
	//bundle->npmax = 100;
	//bundle->tmin = 1e-10;
	////bundle->modPD = 0;
	//bundle->mf = 0.1;

	//// Criação dos modelos de otimização do gurobi para cada cenário
	//// ------------------------------------------------
	//QueryPerformanceFrequency(&clockPerSec);
	//QueryPerformanceCounter(&inicio_Modelo);
	//problemaCenario.push_back(new CSubProblemaCen(sistema_a, 0, ambGRB));
	//for (int i = 1; i < sistema_a->GetNCenarios(); i++)
	//	problemaCenario.push_back(new CSubProblemaCen(sistema_a, i, *problemaCenario[0]));
	//QueryPerformanceCounter(&fim_Modelo);
	//tempoModelo = double(fim_Modelo.QuadPart - inicio_Modelo.QuadPart)/double(clockPerSec.QuadPart);
	//// ------------------------------------------------

	//// Resolver RL
	//// ------------------------------------------------
	//// Lambda inicial
	//if (importar_lambda == true)
	//{
	//	ifstream inFile( "lambda.txt", ios::in );   
	//	if ( !inFile )                                                            
	//	{                                                                               
	//		cout << "File "<< "lambda.txt" << " could not be opened" << endl;
	//		exit( 1 );
	//	}
	//	while ( ! inFile.eof() )
	//	{
	//		for (int i = 0 ; i < sistema_a->GetNCenarios() * n_var_duplicadas; i++)
	//			inFile >> bundle->xEig(i);
	//	}
	//	inFile.close();
	//}
	//else
	//{
	//	//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	//	//	for (int i = 0; i < n_var_duplicadas; i++)
	//	//		Lambda[cen][i] = lambda_inicial;
	//	for (int i = 0; i < sistema_a->GetNCenarios() * n_var_duplicadas; i++)
	//		bundle->xEig(i) = lambda_inicial;
	//		//bundle->xEig(i) = (double(rand())/RAND_MAX)*lambda_inicial;
	//}
	////
	//QueryPerformanceCounter(&inicio_RL);
	/////subgradiente(itmax);
	//bundle->optimize();
	//for (int i = 0; i < bundle->fdual.size(); i++)
	//	bundle->fdual(i) = - bundle->fdual(i);
	//QueryPerformanceCounter(&fim_RL);
	////
	//tempoRL = double(fim_RL.QuadPart - inicio_RL.QuadPart)/double(clockPerSec.QuadPart);
	////
	////QueryPerformanceCounter(&inicio_Gravar_Escrever);
	////QueryPerformanceCounter(&fim_Gravar_Escrever);
	////tempoGravar_Escrever = double(fim_Gravar_Escrever.QuadPart - inicio_Gravar_Escrever.QuadPart)/double(clockPerSec.QuadPart);
	//tempoGravar_Escrever = 0;
	////
	//resultadosGurobi->EscreverArquivoDecCen("out_temp.txt", tempoRL, &fdual_RL, &x, &Lambda, &norma_sg, &x_med);		// Escreve resultados da última iteração do RL
	//resultadosGurobi->EscreverArquivoIterativoDecCen("out_temp2.txt", bundle->normS, bundle->fdual, tempoRL, tempoModelo, tempoResolucao, tempoGravar_Escrever, tempo_cenarios);
	//if (sistema_a->GetFlagVarBin() == false)
	//	resultadosGurobi->ExportarSol("lambda.txt", &Lambda);
	////
	//// ------------------------------------------------
}


// na RL criar subproblemas e depois deletá-los para criá-los novamente na RP ??? ou só alterar ?
// entrada de dados na RP (executar só RP ou com valores da RL)
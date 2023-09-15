/*--------------------------------------------------------------------------*/
/*-------------------------- File Scndec3Fi.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Scndec3Fi is a "concrete" class which implements the
 * interface defined by the abstract base class FiOracle. Used to solve the
 * Scenario Decomposition scheme 3
 *
 * \version 0.50
 *
 * \date 16 - 05 - 2012
 *
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- INCLUDES -------------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Scndec3Fi.h"

#include "OPTvect.h"

#include <math.h>

#include "NDOSlver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
	using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTANTS -----------------------------------*/
/*--------------------------------------------------------------------------*/

static cIndex InINF = Inf<Index>();

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

Scndec3Fi::Scndec3Fi( CSistema * const sistema_end , Scndec3Results * const resultadosGurobi_end, int tipo_dual)
	:
        FiOracle(), ambGRB(GRBEnv())
{
	sistema_a = sistema_end;
	resultadosGurobi = resultadosGurobi_end;
	
	fdual_RL = 0;
	norma_sg = 0;
	JJ = 0;
	for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
		JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
	flag1 = int (sistema_a->GetFlagModeloRede());
	flag2 = int (sistema_a->GetFlagVfol());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());

	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 4*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	na = sistema_a->GetTt1() * (2*sistema_a->termeletricasVtr.size() + sistema_a->hidreletricasVtr.size());
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)
		n -= T * (sistema_a->barrasVtr.size() - 1);

	sistema_a->SetPreconditioners(GenPreconditioners());		// Define os condicionadores
	x.resize(n);
	//Lambda.resize(na*(sistema_a->GetNCenarios() - 1));
	//subgrad.resize(na*(sistema_a->GetNCenarios() - 1));
	// ordem dos vetores: L[0] L[1],..., L[n_cen], em que L é o vetor nas linhas abaixo, de tamanho na

	nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / T);
	//nat = (na / T);
	
	// Inicialização para o Bundle:
	NumVar = na*(sistema_a->GetNCenarios() - 1);		// Numero de multiplicadores de Lagrange
	NumQVar = 0;
	NumQCmp = 0;
	if ( tipo_dual == 0)	// aggregated
	{
		NumPCmp = 1;
		NumECmp = 0;
		EasyComp = 0;
		Aggrgtd = true;
	}
	else		// disaggregated
	{
		NumPCmp = sistema_a->GetNCenarios();		// Dividir subproblemas em varios componentes
		NumECmp = 0;
		EasyComp = 0;
		Aggrgtd = false;
	}
	// Nesse tipo de decomposição não existem componentes "easy"
	
	NumHCmp = NumQCmp + NumPCmp - NumECmp;

	resultadosGurobi->SetAgrModel(Aggrgtd);

	// lambda:
	//L = [Lpt Lu Ld Lphg] -> (.) tamanho dos vetores que dependem de flags
	// Subproblemas:
	//1 para cada cenário

	last_wFi = 0;
	New_wFi = new bool[ NumHCmp ];

	SetLambda( 0 );
	LamB = 0;
	LamBd = NumVar;

	SGBse1 = new Index[ NumVar + 1 ];
	LowerB = - Inf<HpNum>();

	heuristica = NULL;
	Hrstt = 0;
	//last_setprcs.resize(sistema_a->GetNCenarios());
	//for (size_t i = 0; i < last_setprcs.size(); i++)
	//	last_setprcs[i] = 2;

	iter_count = 0;
	t_init = 0;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void Scndec3Fi::SetData( cLMRow UB  )
{
	UpBnd = UB;    
 }

void Scndec3Fi::CriarModelos( void )
{
	// cenário 0
	vtr_subproblemaC.push_back(new Scndec3SPC(sistema_a, ambGRB));

	// demais cenários
	for (int n_cen = 1; n_cen < sistema_a->GetNCenarios(); n_cen++)
		vtr_subproblemaC.push_back(new Scndec3SPC(sistema_a, n_cen, *vtr_subproblemaC[0]));
	
	int * comp_inf;		// vetor com informaçoes do n. variaveis de cada componente
	comp_inf = new int[sistema_a->GetNCenarios()];
	// mesmo no modelo agregado deve-se saber essa informaçao acima, pois os subproblemas ainda sao resolvidos separadamente!!

	for (int i = 0; i < sistema_a->GetNCenarios(); i++)
		comp_inf[i] = vtr_subproblemaC[i]->GetNVarPorComp();		// nessa decomposição todos os componentes tem o mesmo numero de variáveis

	resultadosGurobi->SetComponents(comp_inf, Aggrgtd);
	delete comp_inf;
}

void Scndec3Fi::SetHeuristic ( Hrstc *heuristic )
{
	heuristica = heuristic;
}

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

bool Scndec3Fi::GetUC( cIndex i )
{
	return( true );
}

LMNum Scndec3Fi::GetUB( cIndex i )
{
	return( Inf<LMNum>() );
}

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/

void Scndec3Fi::SetLambda( cLMRow Lmbd )
{
	Lam1 = Lmbd;
	for( Index i = 0 ; i < NumHCmp ; i++ )
		New_wFi[ i ] = true;
 }

/*--------------------------------------------------------------------------*/

void Scndec3Fi::SetLamBase( cIndex_Set LmbdB , cIndex LmbdBD )
{
	LamB = LmbdB;
	LamBd = LmbdBD;
}

bool Scndec3Fi::SetPrecision( HpNum Eps )
{
	//// manipulando reset()
	//vetorint setprcs;
	//setprcs.resize(sistema_a->GetNCenarios());
	//int count_true = 0;
	//for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	//{
	//	if (vtr_subproblemaC[n_cen]->GetMipgap() <= Eps)
	//	{
	//		setprcs[n_cen] = 0;
	//		//if (vtr_subproblemaC[n_cen]->GetResetStt() == 0)		// reativar reset()
	//		//	vtr_subproblemaC[n_cen]->SetResetStt(1);
	//	}
	//	else	// > Eps
	//	{
	//		if (vtr_subproblemaC[n_cen]->GetStatus() == 2)
	//		{
	//			setprcs[n_cen] = 1;
	//			vtr_subproblemaC[n_cen]->SetMaxMipgap( Eps );
	//		}
	//		else	// != 2
	//		{
	//			if (vtr_subproblemaC[n_cen]->GetResetStt() == 1)
	//			{
	//				setprcs[n_cen] = 1;
	//				vtr_subproblemaC[n_cen]->SetResetStt(0);
	//				vtr_subproblemaC[n_cen]->SetMaxMipgap( Eps );
	//			}
	//			else
	//				setprcs[n_cen] = -1;
	//		}
	//	}
	//	if (setprcs[n_cen] > 0)
	//		count_true++;
	//}
	//
	//cout << "SetPrecision called: Eps = " << Eps  << ". Situação dos subp = [";
	//for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	//	cout << setprcs[n_cen] << " ";
	//cout << "]" << endl;
	//
	//cout << "Reset() = [";
	//for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	//	cout << vtr_subproblemaC[n_cen]->GetResetStt() << " ";
	//cout << "]" << endl;

	//cout << "Mipgap = [";
	//for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	//	cout << vtr_subproblemaC[n_cen]->GetMipgap() << " ";
	//cout << "]" << endl;

	//cout << "Status = [";
	//for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	//	cout << vtr_subproblemaC[n_cen]->GetStatus() << " ";
	//cout << "]" << endl;
	//
	//if (count_true > 0)
	//	return (true);
	//else
	//	return (false);

	// -------------------------------------------------------------------------------------------
	// sem manipular reset()
	double total_mipgap = 0;
	int total_status = 0;
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
	{
		total_mipgap += vtr_subproblemaC[n_cen]->GetMipgap();
		if ( vtr_subproblemaC[n_cen]->GetStatus() != 2)
			total_status++;
	}

	//cout << total_mipgap << " : " << total_status << endl;
	// media do mipgap (soma é muito alta para problemas com muitos cenários), ou o valor esperado
	// valor esperado pois a imprecisão em cada subproblema é multiplicada pela sua probabilidade
	total_mipgap = total_mipgap / sistema_a->GetNCenarios();

	cout << "SetPrecision called: Eps = " << Eps  << ". Prec. atual = " << total_mipgap << endl;

	if ( ( total_mipgap > Eps ) && (total_status == 0) )
	{
		//cout << "Precision changed to: " << Eps << endl;
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			vtr_subproblemaC[n_cen]->SetMaxMipgap( Eps );
		return ( true );
	}
	else
		return ( false );

}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

HpNum Scndec3Fi::Fi( cIndex wFi )
{
	// Resolver subproblemas aqui!!!
	HpNum temp = 0;
	
	// timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Start();

	if( wFi == 0 )
	{
		if( Fit )
			Fit->Stop();
		return( 0 );
	}

	//imprimir o lambda
	//ofstream log_lambda("lambda.txt", ofstream::app );
	//log_lambda << wFi << char(9);

	// Converter lambda
	vetorfloat2 lambda;
	lambda.resize(sistema_a->GetNCenarios() - 1);
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios() - 1; n_cen++)
		for (int i = 0; i < na; i++)
		{
			///lambda[n_cen].push_back(Lam1[i + n_cen*na]*sistema_a->GetPrecondidioner(i));
			lambda[n_cen].push_back(Lam1[i + n_cen*na]);
			//log_lambda << Lam1[i + n_cen*na] << char(9);
		}

	//log_lambda << endl;

	// lambda:
	//L = [Lpt Lu Lv] -> (.) tamanho dos vetores que dependem de flags
	// Subproblemas:
	// 1 para cada cenário

	// Resoluçao dos subproblemas
	// ------------------------------------------------
	// Disaggregated: the full function Fi except that of the "easy" components or
	// Aggregated model
	if ((wFi > NumPCmp) || (Aggrgtd == true))
	{
		for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
			vtr_subproblemaC[n_cen]->ResolverProblemaRL(resultadosGurobi, lambda, 0);
		temp = - resultadosGurobi->GetFobj();
	}
	else if ( (wFi > 0) && (wFi <= NumPCmp ) )
	{
		vtr_subproblemaC[wFi - 1]->ResolverProblemaRL(resultadosGurobi, lambda, 0);
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
		//cout << "Subp" << wFi << " = " << temp << endl;
	}
	else
	{
		throw NDOException( "Scndec3Fi: Fi called fon an easy component." );	// Essa funçao nao pode ser chamada para easy components
		return( - Inf<HpNum>() );
	}
	//x = resultadosGurobi->GetX();
	//fdual_RL = - resultadosGurobi->GetFobj();
	// ------------------------------------------------

	// timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Stop();

	return temp;
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool Scndec3Fi::NewGi( cIndex wFi )
{
	last_wFi = wFi;
	if( wFi == 0 )
		throw NDOException( "NewGi: wFi = 0" );
	bool x;
	int iNew_wFi;
	if( wFi <= NumPCmp )
	{
		iNew_wFi = wFi - 1;
		x  = New_wFi[ iNew_wFi ];		// New_wFi tem o tamanho das hard componentes
		New_wFi[ iNew_wFi ] = false;
		return( x );
	}
	if( wFi > NumPCmp )
		return (true);

	return( false );
 }

/*--------------------------------------------------------------------------*/

Index Scndec3Fi::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name , cIndex strt , Index stp )
{
	// timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Start();

	if( stp > NumVar )
		stp = NumVar;

	//if( Name < MBN )
	//	throw NDOException( "GetGi: past information is not recorded" );
	
	Index temp;
	
	if( Name < MBN )
	{
		cout << "NoK1" << endl;
		SGBse = 0;		// Vetor de subgradientes denso
		resultadosGurobi->GetSubGradiente(comp_names[Name], SubG, ptr_x_names[Name]);
		temp = NumVar;
	}

	if( Name == MBN )
	{
		SGBse = &InINF;
		if( Fit )
			Fit->Stop();
		return( 0 );
	}

	if( Name > MBN )
	{
		// Calculo do Subgradiente
		// ------------------------------------------------
		// Disaggregated com Inf<Index>() >= wFi > NumHCmp or
		// Aggregated
		if ((last_wFi > NumPCmp) || (Aggrgtd == true))		
		{
			SGBse = 0;		// Vetor de subgradientes denso
			//resultadosGurobi->GetSubGradiente(SubG);
			resultadosGurobi->GetSubGradiente(SubG);
			temp = NumVar;
		}
		else	// Disaggregated (last_wFi <= NumPCmp) - vetor de subgradientes esparso
		{
			// Essa funçao n é chamada para as easy components, pois é chamado NewGi() antes
			resultadosGurobi->GetSubGradiente(last_wFi, SubG);
			// turn SubG from a "dense" NumVar-vector to a "sparse" one - - - - - - - - -
			Index_Set SGBse2 = Sparsify( SubG , SGBse1 , NumVar );
			*SGBse2 = Inf<Index>();
			temp = SGBse2 - SGBse1;

			SGBse = SGBse1;
		}
	}

	// timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Stop();

	return( temp );
 }

/*--------------------------------------------------------------------------*/

HpNum Scndec3Fi::GetVal( cIndex Name )
{
	if ( Name < MBN )
		return ( epsilon_names[Name] );
	if ( Name == MBN )
		return( 0 );
	if ( Name > MBN )
		return ( resultadosGurobi->GetEpsilon(last_wFi - 1) );

	//if( Name < MBN )
	//	throw( NDOException( "GetVal: past information is not recorded" ) );

	return( 0 );

	// Let me mention one thing, though, that may not be obvious: you have two
	// different choice about "who the function value is".
	// The first choice is that the function value is the (Lagrangian) value of
	// the feasible integer solution you got. In this case the corresponding 
	// subgradient is (treated as) a 0-subgradient, which is not because that
	// value is not the true function value, which is larger (assuming you 
	// minimize Fi, which means your optimization problem is a max). This is 
	// when you get negative linearization errors.
	// You can avoid them, though. You simply have to report as the function 
	// value the *best upper bound* that Gurobi gives you when you terminate
	// the B&C beforehand. This is larger than the other value, and durely a
	// correct upper approximation of the true Fi-value; hence, you never have
	// negative linearization errors. In this case, however, GetVal() has to 
	// return the difference between the two values, because the subgradient 
	// you give is instead attached to the integer solution and therefore it's
	// not a 0-subgradient but an eps-subgradient, with eps precisely that difference.
 }

/*--------------------------------------------------------------------------*/

void Scndec3Fi::SetGiName( cIndex Name )
{
	if ( !ptr_x_names[Name].empty() )
		ptr_x_names[Name].clear();
	ptr_x_names[Name] = resultadosGurobi->GetCompX(last_wFi - 1);

	// também gravar um vetor com o numero da componente que o item name pertence!!
	comp_names[Name] = last_wFi;
	epsilon_names[Name] = resultadosGurobi->GetEpsilon(last_wFi - 1);
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

HpNum Scndec3Fi::GetLowerBound( cIndex wFi )
{
	// Chamada para wFi = Inf e depois para cada componente wFi = 1,...,GetNrFi()
	// This is called at every iteration, right after Fi() and all the GetGi()
	// A feasible solution gives an upper bound for the minimization problem, and since we invert the sign to minimize, a lower bound on our Fi()
	if (wFi == 0)
		return( - Inf<HpNum>() );
	if( wFi > GetNrFi() )
	{
		// Para aplicar a heuristica em somente algumas iteraçoes é só colocar um contador aqui e aplicar somente qdo contador == NumIteHeur else só retornar GetFO();
		if ( heuristica != NULL )		
		{
			// timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			if( Hrstt )
				Hrstt->Start();
			// Não é necessário chamar aqui, já foi chamado em GetFiStatus()!!!!
			//ConvSol();		// precisa ser chamado caso a soluçao convexificada seja usada na heuristica
			//cout << "GetLowerBound iter= " << Slvr->NrIter() << endl;
			if ((Slvr->NrIter() == 1) && ((heuristica->GetInitCut() == 0) || (heuristica->GetBenders() == 0)))	// se não forem adicionados cortes do modelo estendido na iteração 0 do Bundle, não se faz nada (x_til e x_hat não são definidos ainda)
			{
				*log_auxiliar << endl;
				LowerB = - Inf<HpNum>();
			}
			else
			{
				if ( log_auxiliar->is_open() )
				{
					*log_auxiliar << char(9) << "GetLB()" << char(9) << Slvr->NDOTime();
				}
				else
					cout << "Unable to open file";

				LowerB = - heuristica->ResolverHeuristica();
				// adicionar tempo da heuristica (somente quando ResolverHeuristica() é chamado)
				if ( log_auxiliar->is_open() )
				{
					*log_auxiliar << char(9) << Slvr->NDOTime();
					*log_auxiliar << endl;
				}
				else
					cout << "Unable to open file";
			}
			
			// qdo a soluçao é infinito passar o Inf<HpNum>()
			if ( LowerB <= - GRB_INFINITY)
				LowerB = - Inf<HpNum>();

			// timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			if( Hrstt )
				Hrstt->Stop();
			// Na iterção 0 e primeira 1 (usadas para calcular parametros do Bundle) adiciona-se somente os cortes??? alterar parametros??
			//if (Slvr->NrIter() > 10)
				return( LowerB );
			//else
			//	return( - Inf<HpNum>() );
		}
		else
			return( - Inf<HpNum>() );
	}
	else
		return( - Inf<HpNum>() );
}

FiOracle::FiStatus Scndec3Fi::GetFiStatus( Index wFi )
{
	// Fazer calculo da medida de inviabilidade e imprimir em arquivo auxiliar...
	// <iter> <Fi> <Sigma> <t> <n1> <n2> <nInf>
	if (wFi > GetNrFi())
	{
		// Calcula a solução convexificada
		// Não é necessário chamar aqui se GetLowerBound() (heuristica) for chamado!!!!
		if (Slvr->NrIter() != 0)
			ConvSol();

		vetorfloat norma(3, 0);	// n1, n2 e nINF para todas as restrições (independente do tipo)
		resultadosGurobi->CalcularNormaSubgXtil(norma);

		//// Calcular a norma considerando todas variáveis de primeiro estágio
		//vetorfloat norma_compl(3, 0);
		//resultadosGurobi->CalcularNormaSubgXtilTV(norma_compl);

		// Dependendo do modelo dual a norma (e o desvio) utilizada como criterio de para é diferente
		double * norm;
		double desvio;
		if (Aggrgtd)
		{
			norm = &norma[1];
			desvio = 1*pow(NumVar, 0.5);		// esse desvio fica muito baixo!!
			//desvio = 2*NumVar;
		}
		else
		{
			norm = &norma[0];
			HpNum epsB;
			Slvr->GetPar(Bundle::kEpsLin, epsB);
			desvio = epsB*2e4*NumVar;	// com Epslin = 1e-4 tem-se 2*NumVar
			//desvio = 2*NumVar;
		}

		if ( log_auxiliar->is_open() )
		{
			// Fi
			*log_auxiliar << "LR:" << char(9) << Slvr->NrIter() << char(9) << - Slvr->ReadFiVal();
			// Escrever tb valor de sigma e t
			if( BSlvr )
				*log_auxiliar << char(9) << BSlvr->ReadSigma() << char(9) << BSlvr->Readt();
			for (size_t i = 0; i < norma.size(); i++)
				*log_auxiliar << char(9) << norma[i];
			// escrever também o tempo por iteração... gráfico com a evolução do algoritmo
			*log_auxiliar << char(9) << Slvr->NDOTime();
			
			if ( iter_count == 1)
				*log_auxiliar << endl;
			else if (heuristica == NULL)
				*log_auxiliar << endl;
		}
		else
			cout << "Unable to open file";

		// Aplica t inicial na segunda iteração (usar o iter_cout pois a iteração se repete quando kFiChgd)
		//HpNum tinicial;
		if (iter_count == 1)
		{
			// na interação 0 os suproblemas não são resolvido, portanto x = 0. Assim t_init deve ser calculado na iteração 1
			//t_init = (1 + abs(fdual_cont))/(5*norma[0]);		// fdual_cont é o valor do prob. extendido continuo
			t_init = (1 + abs(-Slvr->ReadFiVal()))/(500*(*norm));	// aqui então utilização o valor da função dual avaliada para lambda da relax. continua
			//cout << "norma :" << norma[0] << endl << "t0 :" << t_init << endl;
			// Arredondar t_init
			// Contar ordem de t_init para deixá-lo com 2 algarismos significativos
			int ordem = 0;
			double result = 0;
			if (t_init >= 1)	// maior que 1, entao a ordem é definida dividindo
			{
				result = t_init;
				ordem--;
				while ( result >= 1 )
				{
					ordem++;
					result = result / 10;
				}
			}
			else	// t_init entre 0 e 1, entao 
			{
				result = t_init;
				while ( result < 1 )
				{
					ordem--;
					result = result * 10;
				}
			}
			ordem = -ordem + 2;			// 2 é o numero de algarismos significativos desejados para t_init
			HpNum tmax, tmin;
			Slvr->GetPar(Bundle::ktMinor, tmin);
			Slvr->GetPar(Bundle::ktMaior, tmax);
			t_init = min( tmax , max( tmin , ceil(t_init*pow(10.0, ordem) - 0.5)/pow(10.0, ordem) ) );

			// Tolerancia para o criterio de parada
			// 0.1 desvio medio requerido para cada restrição relaxada x numero de restrições(ou var. duais)
			//tol = t_init*desvio/(1 + fdual_cont);
			//tol = desvio/(5*(*norm));
			
			//Slvr->GetPar(Bundle::ktInit, tinicial);
			//tinicial = t_init;
			//tinicial = 0.1;
			if ( t_init < 1)
				tol = t_init*desvio/(1 + abs(fdual_cont));
			else
				tol = desvio/(1 + abs(fdual_cont));

			//cout << tol << endl;
			//tol = 1e-5;		// justificada pelo erro aceitável em sigma e tb em t*g, na casa das dezenas

			Slvr->SetPar(Bundle::ktStar, 50*t_init);		// altera tStar mas no log fica o valor inicial
			Slvr->SetPar(Bundle::ktCurr, t_init);	// esse t_init não é bom!!!!!

			// parei aqui
			// a solução dos subproblemas na it. 2 é diferentes da sol. na iteração 1 (é o mesmo ponto de lambda, mas o oracle entrega uma sol. diferente?!?)
			///Slvr->SetPar(Bundle::ktCurr, 1e-8);
			
			//HpNum test1;
			//Slvr->GetPar(Bundle::ktStar, test1);
			//cout << "tStar = " << test1 << " ; tCurr = " << BSlvr->Readt() << " ; " << endl;


			iter_count++;
			return ( kFiChgd );
		}

		// Calcula o t inicial da primeira iteração
		if (Slvr->NrIter() == 0)
			iter_count++;
		
		// Criterio de parada pelo valor relativo do sigma-subgradente e do sigma
		else if ( (abs(BSlvr->ReadSigma()/(1 - Slvr->ReadFiVal())) <= tol) && (abs((BSlvr->Readt())*(*norm)/(1 - Slvr->ReadFiVal())) <= tol) )
		//else if ( (abs(BSlvr->ReadSigma()/(1 - Slvr->ReadFiVal())) <= tol) && (abs((BSlvr->Readt()/tinicial)*(*norm)/(1 - Slvr->ReadFiVal())) <= tol) )
		//else if ( (abs(BSlvr->ReadSigma()/(1 - Slvr->ReadFiVal())) <= tol) && (abs((*norm)/(1 - Slvr->ReadFiVal())) <= tol) )
		{
			///if ( abs((*norm)/(1 - Slvr->ReadFiVal())) <= desvio/(1 + abs(fdual_cont)) )		// casos em que t pode ser muito baixo
				return( kFiStop );
			///else
			///	return( kFiCont );
		}
		// Criterio de parada pelo gap entre Fi e o valor ótimo
		else if (((-LowerB + Slvr->ReadFiVal()) / (-Slvr->ReadFiVal()) <= 1) && (heuristica != NULL))// se a diferença for menor que 10%
			return(kFiStop);
		else
		{
			//// ajuste da tolerancia de cada subproblema pelos valores de sigma e norma
			//if (Slvr->NrIter() > 1)
			//{
			//	HpNum new_mipgap = max( abs(BSlvr->ReadSigma()/(1 - Slvr->ReadFiVal())), abs((*norm)/(1 - Slvr->ReadFiVal())) ) / 10;
			//	bool aa = SetPrecision(new_mipgap);
			//}
			return( kFiCont );
		}
		
		//return( kFiNorm );
	}

	// Chamada para wFi = Inf e depois para cada componente wFi = 1,...,GetNrFi()
	// called at every iteration after that the master problem has been solved but before SetLambda() and Fi() are
	return( kFiNorm );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
	
void Scndec3Fi::Deleted( cIndex Name )
{
	ptr_x_names[Name].clear();
	//vetorfloat().swap(ptr_x_names[Name]);
	epsilon_names[Name] = NULL;
}

void Scndec3Fi::Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm )
{
	throw NDOException( "Aggregate: Solucao convexificada nao agregada!!" );
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

Scndec3Fi::~Scndec3Fi()
{
	for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		delete vtr_subproblemaC[n_cen];

	BSlvr = NULL;
	heuristica = NULL;
	log_auxiliar = NULL;
	resultadosGurobi = NULL;
	sistema_a = NULL;

	delete[] New_wFi;
	delete[] SGBse1;

	vtr_subproblemaC.clear();

	ptr_x_names.clear();
	epsilon_names.clear();
	comp_names.clear();

	delete Hrstt;
	ambGRB.~GRBEnv();
 }
/*--------------------------------------------------------------------------*/
/*--------------------------------- OTHER ----------------------------------*/
/*--------------------------------------------------------------------------*/

void Scndec3Fi::ConvSol(  )
{
	// Calcula a soluçao convexificada e armazena na classe Scndec1Results (ptr_x_til)
	vetorfloat * ptr_x_a = resultadosGurobi->GetPtrXtil();		//[ptr_x_til] => numero de componentes x numero de var. em cada componente 
	cHpRow mult;
	Index D;
	cIndex_Set I;
	// compute ptr_x_til = Sum{i \in D[]} Theta[ i ] * ptr_x_names[ i ] - - -
	if (Aggrgtd == true)		// modelo agregado
	{
		mult = Slvr->ReadMult( I, D, 1 );
		if (D != 0)
		{
			if( I )
			{
				for (size_t j = 0; j < ptr_x_names[ I[0] ].size(); j++)		// loop nas variáveis de cada componente
				{
					ptr_x_a[0][j] = 0;
					for( Index i = 0; i < D; i++ )		// loop nos multiplicadores (names)
						ptr_x_a[0][j] += mult[i] * ptr_x_names[ I[ i ] ][j];
						//ptr_x_a[comp][j] += mult[i];
				}
			}
			else
			{
				for (size_t j = 0; j < ptr_x_names[ 0 ].size(); j++)
				{
					ptr_x_a[0][j] = 0;
					for( Index i = 0; i < D; i++ )
						ptr_x_a[0][j] += mult[i] * ptr_x_names[ i ][j];
				}
			}
		}
		else
			if (Slvr->ReadLBMult() != 1)
				throw( NDOException( "Scndec1Fi::ConvSol: D == 0" ) );
	}
	else	// modelo desagregado: calculo para cada componente separado
	{
		for (int comp = 0; comp < NumPCmp; comp++)	
		{
			mult = Slvr->ReadMult( I, D, comp + 1 );
			if (D != 0)
			{
				if( I )
					for (size_t j = 0; j < ptr_x_names[ I[0] ].size(); j++)		// loop nas variáveis de cada componente
					{
						ptr_x_a[comp][j] = 0;
						for( Index i = 0; i < D; i++ )		// loop nos multiplicadores (names)
							ptr_x_a[comp][j] += mult[i] * ptr_x_names[ I[ i ] ][j];
					}
				else
					for (size_t j = 0; j < ptr_x_names[ 0 ].size(); j++)
					{
						ptr_x_a[comp][j] = 0;
						for( Index i = 0; i < D; i++ )
							ptr_x_a[comp][j] += mult[i] * ptr_x_names[ i ][j];
					}
			}
			else
				if (Slvr->ReadLBMult() != 1)
					throw( NDOException( "Scndec1Fi::ConvSol: D == 0" ) );
		}
	}

	// Adicionar componentes do Lower Bound!!
	if (Slvr->ReadLBMult() != 0)
	{
		HpNum multiLB = Slvr->ReadLBMult();
		resultadosGurobi->AdicionarLBXtil(heuristica->GetX(), multiLB);
	}

	ptr_x_a = NULL;
	delete ptr_x_a;
 }

HpNum Scndec3Fi::GetHeuristicSolution(  )
{
	return heuristica->GetFO();
}

vetorfloat Scndec3Fi::GenPreconditioners()
{
	// aqui é um pouco diferente, pois a ordem dos multiplicadores é diferente, porém o tamanho do vetor de precond. é de na*num_cenarios

	vetorfloat precond;
	// lambda:
	//L = [Lpt Lu Lcp (LF) (Lteta) Lph Lv Ld Ls (Lphmax) Lphg Lq Lz (Ldef)] -> (.) tamanho dos vetores que dependem de flags
	//L = [Lpt Lu Lup Lud (LF) (Lteta) Lph Lv Ld Ls (Lphmax) Lphg Lq Lz (Ldef)] -> (.) tamanho dos vetores que dependem de flags
	int flag3 = int (sistema_a->GetFlagPhmax());
	// Kind of Preconditioner
	switch (sistema_a->GetFlagPrecond())
	{
	case 1:		//probability of the node
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios() - 1; cen++)
			{
				for (size_t i = 0; i < na; i++)
					precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));
			}
			break;
		}
	case 2:		//square root of the probability of the node
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios() - 1; cen++)
			{
				for (size_t i = 0; i < na; i++)
					precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));
			}
			break;
		}
	case 3:		//maximum value of the relaxed constraint
		{
			size_t I = sistema_a->termeletricasVtr.size();
			size_t R = sistema_a->hidreletricasVtr.size();
			size_t B = sistema_a->barrasVtr.size() - 1;
			for (int cen = 0; cen < sistema_a->GetNCenarios() - 1; cen++)
			{
				//x = [pt u v]
				//x = [pt u v]	"modern" model
				for (int t = 0; t < sistema_a->GetTt1(); t++)
				{
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					{
						for (size_t i = 0; i < I; i++)	//pt
							precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					}
					else
					{
						for (size_t i = 0; i < I; i++)	//pt
							precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
					}
					for (size_t i = 0; i < I; i++)	//u
						precond.push_back(1 / 1);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						precond.push_back(1 / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
				}
			}
			break;
		}
	case 4:		//combination between 1 and 3
		{
			size_t I = sistema_a->termeletricasVtr.size();
			size_t R = sistema_a->hidreletricasVtr.size();
			size_t B = sistema_a->barrasVtr.size() - 1;
			for (int cen = 0; cen < sistema_a->GetNCenarios() - 1; cen++)
			{
				//x = [pt u d phg]
				for (int t = 0; t < sistema_a->GetTt1(); t++)
				{
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					{
						for (size_t i = 0; i < I; i++)	//pt
							precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					}
					else
					{
						for (size_t i = 0; i < I; i++)	//pt
							precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
					}
					for (size_t i = 0; i < I; i++)	//u
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / 1);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					{
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
					}

				}
			}
			break;
		}
	case 5:		//combination between 2 and 3
		{
			size_t I = sistema_a->termeletricasVtr.size();
			size_t R = sistema_a->hidreletricasVtr.size();
			size_t B = sistema_a->barrasVtr.size() - 1;
			for (int cen = 0; cen < sistema_a->GetNCenarios() - 1; cen++)
			{
				//x = [pt u d phg]
				for (int t = 0; t < sistema_a->GetTt1(); t++)
				{
					if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					{
						for (size_t i = 0; i < I; i++)	//pt
							precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					}
					else
					{
						for (size_t i = 0; i < I; i++)	//pt
							precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
					}
					for (size_t i = 0; i < I; i++)	//u
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / 1);
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					{
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
					}
				}
			}
			break;
		}
	default:	//none of them
		{
			for (int cen = 0; cen < sistema_a->GetNCenarios() - 1; cen++)
			{
				for (size_t i = 0; i < na; i++)
					precond.push_back(1);
			}
			break;
		}
	}
	return precond;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------ End File Scndec3Fi.cpp ------------------------*/
/*--------------------------------------------------------------------------*/



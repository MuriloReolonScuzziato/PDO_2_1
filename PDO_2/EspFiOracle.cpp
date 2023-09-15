/*--------------------------------------------------------------------------*/
/*-------------------------- File EspFiOracle.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * TestFi is a simple example of a "concrete" class which implements the
 * interface defined by the abstract base class FiOracle.
 *
 * \version 0.50
 *
 * \date 16 - 05 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n

 * \author Andrea Nerli \n
 *
 * Copyright &copy 2001 - 2012 by Antonio Frangioni.
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- INCLUDES -------------------------------------*/
/*--------------------------------------------------------------------------*/

#include "EspFiOracle.h"

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

EspFiOracle::EspFiOracle( CSistema * const sistema_end , ResultadosConj * const resultadosGurobi_end, int tipo_dual)
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
	int n_var_cen;
	flag1 = int (1 - sistema_a->GetFlagBarraUnica());
	flag2 = int (sistema_a->GetFlagVfol());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());

	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + 5*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	na = T * (sistema_a->termeletricasVtr.size() + 4*sistema_a->hidreletricasVtr.size());
	if ( sistema_a->GetFlagBarraUnica() == true )		// Remover teta, def e deixar uma variável de deficit por período
		n -= (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() ) * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	
	x.resize(n);
	xa.resize(na);
	Lambda.resize(na);
	subgrad.resize(na);
	
	nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / T);
	nat = (na / T);
	
	// Inicialização para o Bundle:
	NumVar = na;		// Numero de multiplicadores de Lagrange
	NumQVar = 0;
	NumQCmp = 0;
	if ( tipo_dual == 0)
	{
		NumPCmp = 1;
		NumECmp = 0;
		EasyComp = 0;
		Aggrgtd = true;
	}
	else if ( tipo_dual == 1 )
	{
		NumPCmp = resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size() + resultadosGurobi->GetN();		// Dividir subproblemas em varios componentes ou 1 (desagregado ou agregado)
		NumECmp = 0;
		EasyComp = 0;
		Aggrgtd = false;
	}
	else
	{
		NumPCmp = resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size() + resultadosGurobi->GetN();
		Aggrgtd = false;
		if ( tipo_dual == 2 )
		{
			NumECmp = resultadosGurobi->GetCH();								// Subproblemas HA são easy ones
			EasyComp = 1;
		}
		else if ( tipo_dual == 3 )
		{
			NumECmp = resultadosGurobi->GetN();									// Subproblemas D são easy ones
			EasyComp = 2;
		}
		else
		{
			NumECmp = resultadosGurobi->GetCH() + resultadosGurobi->GetN();		// Subproblemas HA e D são easy ones
			EasyComp = 3;
		}
	}
	
	NumHCmp = NumQCmp + NumPCmp - NumECmp;

	resultadosGurobi->SetAgrModel(Aggrgtd);

	// lambda [pt ph v d phmax]
	// Subproblemas HA, HE, T, D
	last_wFi = 0;
	New_wFi = new bool[ NumHCmp ];
	//New_wFi = new bool[ NumPCmp ];		// tem o tamanho de todas as componentes, mas só para ter a posiçoes correspondentes aos wFi

	SetLambda( 0 );
	LamB = 0;
	LamBd = NumVar;

	SGBse1 = new Index[ NumVar + 1 ];
	LowerB = - Inf<HpNum>();

	heuristica = NULL;
	//ptr_x_til = new vetorfloat[NumHCmp];
	Hrstt = 0;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void EspFiOracle::SetData( cLMRow UB  )
{
	UpBnd = UB;    
	//UnCnstr = UC;
 }

void EspFiOracle::CriarModelos( void )
{
	// criar modelos que sao easy components de forma difente, sem o modelo do gurobi, somente com as informaçoes que serao passadas para o MP!!!
	if (EasyComp == 1 || EasyComp == 3)		// Criado como easy component
		subproblemaHA = new ConjSPHA(sistema_a);
	else
		subproblemaHA = new ConjSPHA(sistema_a, ambGRB);
	subproblemaHE = new ConjSPHE(sistema_a, ambGRB);
	subproblemaT = new ConjSPT(sistema_a, ambGRB);
	if (EasyComp == 2 || EasyComp == 3)		// Criado como easy component
		subproblemaD = new ConjSPD(sistema_a);
	else
		subproblemaD = new ConjSPD(sistema_a, ambGRB);

	int * comp_inf;		// vetor com informaçoes do n. variaveis de cada componente
	comp_inf = new int[resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size() + resultadosGurobi->GetN()];
	// mesmo no modelo agregado deve-se saber essa informaçao acima, pois os subproblemas ainda sao resolvidos separadamente!!

	int i_ant = 0;
	for (int i = i_ant; i < resultadosGurobi->GetCH(); i++)
		comp_inf[i] = subproblemaHA->GetNVarPorComp(i);
	i_ant += resultadosGurobi->GetCH();
	for (int i = i_ant; i < i_ant + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size(); i++)
		comp_inf[i] = subproblemaHE->GetNVarPorComp( (i - i_ant) % sistema_a->hidreletricasVtr.size() );
	i_ant += resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size();
	for (int i = i_ant; i < i_ant + sistema_a->termeletricasVtr.size(); i++)
		comp_inf[i] = subproblemaT->GetNVarPorComp(i - i_ant);
	i_ant += sistema_a->termeletricasVtr.size();
	for (int i = i_ant; i < i_ant + resultadosGurobi->GetN(); i++)
		comp_inf[i] = subproblemaD->GetNVarPorComp(i - i_ant);

	resultadosGurobi->SetComponents(comp_inf, Aggrgtd);
	delete comp_inf;
}

void EspFiOracle::SetHeuristic ( Heuristicas *heuristic )
{
	heuristica = heuristic;
	//heuristica = new Heuristicas(sistema_a, resultadosGurobi, tempo_max);
}

//void EspFiOracle::DefinirAjustesHeuristica ( int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a, double tempo_max )
//{
//	heuristica = new Heuristicas(sistema_a, resultadosGurobi, ws_a, wlas_a, bvmw, bvfw, alfa_a, beta_a, gama_a, tempo_max);
//}

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

// somente quando tem-se easy components

Index EspFiOracle::GetBNC( cIndex wFi )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetNVarPorComp(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetNVarPorComp(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	return( 0 );
}

Index EspFiOracle::GetBNR( cIndex wFi )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetNRestricoes(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetNRestricoes(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	throw( NDOException( "EspFiOracle::GetBNR: this component is not easy" ) );
	//return( 0 );
}

Index EspFiOracle::GetBNZ( cIndex wFi )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetNNZsRest(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetNNZsRest(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	throw( NDOException( "EspFiOracle::GetBNZ: this component is not easy" ) );
	//return( 0 );
}

void EspFiOracle::GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			subproblemaHA->GetBDesc( wFi - 1, Bbeg, Bind, Bval, lhs, rhs, cst, lbd, ubd);
		else if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			subproblemaD->GetBDesc(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()), Bbeg, Bind, Bval, lhs, rhs, cst, lbd, ubd);
		else	// difficult component
			throw( NDOException( "EspFiOracle::GetBDesc: this component is not easy" ) );
	}
}

Index EspFiOracle::GetANZ( cIndex wFi, cIndex strt, Index stp )
{
	// Considerando que strt = 0 e stp = Inf<Index>() sempre.
	if ( strt != 0 || stp != Inf<Index>() )
		throw( NDOException( "EspFiOracle::GetANZ: dynamic generation of multipliers not implemented" ) );

	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetANZ(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetANZ(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	throw( NDOException( "EspFiOracle::GetANZ: this component is not easy" ) );
	//return( 0 );
}

void EspFiOracle::GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval , cIndex strt, Index stp )
{
	// Considerando que strt = 0 e stp = Inf<Index>() sempre.
	if ( strt != 0 || stp != Inf<Index>() )
		throw( NDOException( "EspFiOracle::GetANZ: dynamic generation of multipliers not implemented" ) );

	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			subproblemaHA->GetADesc(wFi - 1, Abeg, Aind, Aval);
		else if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			subproblemaD->GetADesc(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()), Abeg, Aind, Aval);
		else	// difficult component
			throw( NDOException( "EspFiOracle::GetANZ: this component is not easy" ) );
	}
}

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/

void EspFiOracle::SetLambda( cLMRow Lmbd )
{
	Lam1 = Lmbd;
	for( Index i = 0 ; i < NumHCmp ; i++ )
		New_wFi[ i ] = true;
 }

/*--------------------------------------------------------------------------*/

void EspFiOracle::SetLamBase( cIndex_Set LmbdB , cIndex LmbdBD )
{
	LamB = LmbdB;
	LamBd = LmbdBD;
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

HpNum EspFiOracle::Fi( cIndex wFi )
{
	if( wFi == 0 )
		return( 0 );

	// Resolver subproblemas aqui!!!

	HpNum temp = 0;
	
	// timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Start();

	// Converter lambda
	vetorfloat lambda;
	lambda.resize(NumVar);
	for( Index i = 0 ; i < NumVar ; i++ )
		lambda[i] = Lam1[ i ]*sistema_a->GetPrecondidioner(i);

	// lambda [pt ph v d phmax]
	// Subproblemas HA, HE, T, D

	// Disaggregated: the full function Fi except that of the "easy" components or
	// Aggregated model
	if ((wFi > NumPCmp) || (Aggrgtd == true))
	{
		// Resoluçao dos subproblemas
		// ------------------------------------------------
		if (EasyComp == 0 || EasyComp == 2)		// Subp. HA só é resolvido se não for easy component
			for (int n = 0; n < subproblemaHA->GetNCascatas(); n++)
				subproblemaHA->ResolverProblemaRL(resultadosGurobi, &lambda, 0, n);		// Subproblema Hidraulico
		for (int n = 0; n < (sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())); n++)
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				subproblemaHE->ResolverProblemaRL(resultadosGurobi, &lambda, 0, r, n);	// Subproblema Hidreletrico
		for (int n = 0; n < sistema_a->termeletricasVtr.size(); n++)
			subproblemaT->ResolverProblemaRL(resultadosGurobi, &lambda, 0, n);		// Subproblema Termeletrico
		if (EasyComp == 0 || EasyComp == 1)
			for (int n = 0; n < (sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())); n++)
				subproblemaD->ResolverProblemaRL(resultadosGurobi, &lambda, 0, n);		// Subproblema Demanda
		// ------------------------------------------------
		temp = - resultadosGurobi->GetFobj();
	}		
	else if ( (EasyComp == 0 || EasyComp == 2) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico ( 1 <= comp <= CH )
	{
		subproblemaHA->ResolverProblemaRL(resultadosGurobi, &lambda, 0, wFi - 1);
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else if ( (wFi > resultadosGurobi->GetCH()) && (wFi <= resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size()) )		//Hidreletrico
	{
		int r = (wFi - resultadosGurobi->GetCH() - 1) % int (sistema_a->hidreletricasVtr.size());
		int no = (wFi - resultadosGurobi->GetCH() - 1) / int (sistema_a->hidreletricasVtr.size());
		subproblemaHE->ResolverProblemaRL(resultadosGurobi, &lambda, 0, r, no);
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else if ( (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size()) && (wFi <= resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) )		//Termeletrico
	{
		int i = wFi - resultadosGurobi->GetCH() - resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() - 1;
		subproblemaT->ResolverProblemaRL(resultadosGurobi, &lambda, 0, i);
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else if ( (EasyComp == 0 || EasyComp == 1) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp ) )		//Demanda
	{
		int no = wFi - resultadosGurobi->GetCH() - resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() - sistema_a->termeletricasVtr.size() - 1;
		subproblemaD->ResolverProblemaRL(resultadosGurobi, &lambda, 0, no);
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else
	{
		throw NDOException( "EspFiOracle: Fi called fon an easy component." );	// Essa funçao nao pode ser chamada para easy components
		return( - Inf<HpNum>() );
	}
	//x = resultadosGurobi->GetX();
	//xa = resultadosGurobi->GetXa();
	//fdual_RL = - resultadosGurobi->GetFobj();
	// ------------------------------------------------

	// timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Stop();

	return temp;

	//throw NDOException( "TestFi: Fi called fon an easy component." );
	//return( - Inf<HpNum>() );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool EspFiOracle::NewGi( cIndex wFi )
{
	if( wFi == 0 )
		throw NDOException( "NewGi: wFi = 0" );
	if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico ( 1 <= comp <= CH )
		throw NDOException( "EspFiOracle: NewGi called fon an easy component." );
	if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp ) )		//Demanda
		throw NDOException( "EspFiOracle: NewGi called fon an easy component." );
	last_wFi = wFi;
	bool x;
	int iNew_wFi;
	if( wFi <= NumPCmp )		// wFi nao é easy, pois passou pelas condiçoes acima
	{
		if (EasyComp == 1 || EasyComp == 3)		// se o subp. HA é easy o indice do New_wFi correspondente à wFi é diferente
			iNew_wFi = wFi - 1 - resultadosGurobi->GetCH();		// subtrair de wFi o numero easy decomponentes 
		else
			iNew_wFi = wFi - 1;

		x  = New_wFi[ iNew_wFi ];		// New_wFi tem o tamanho das hard componentes
		New_wFi[ iNew_wFi ] = false;
		return( x );
		//x  = New_wFi[ wFi - 1 ];
		//New_wFi[ wFi - 1 ] = false;
		//return( x );
	}
	if( wFi > NumPCmp )
		return (true);

	return( false );
 }

/*--------------------------------------------------------------------------*/

Index EspFiOracle::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name , cIndex strt , Index stp )
{
	if( Name < MBN )
	throw NDOException( "GetGi: past information is not recorded" );

	if( stp > NumVar )
		stp = NumVar;

	if( Name == MBN )
	{
		SGBse = &InINF;
		return( 0 );
	}

	Index temp;
	
	// timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Start();

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

	// timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Stop();

	return( temp );
 }

/*--------------------------------------------------------------------------*/

HpNum EspFiOracle::GetVal( cIndex Name )
{
 if( Name < MBN )
  throw( NDOException( "GetVal: past information is not recorded" ) );

 return( 0 );
 }

/*--------------------------------------------------------------------------*/

void EspFiOracle::SetGiName( cIndex Name )
{
	if ( !ptr_x_names[Name].empty() )
		ptr_x_names[Name].clear();
		//vetorfloat().swap(ptr_x_names[Name]);

	ptr_x_names[Name] = resultadosGurobi->GetCompX(last_wFi - 1);
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

HpNum EspFiOracle::GetLowerBound( cIndex wFi )
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
			ConvSol();		// precisa ser chamado caso a soluçao convexificada seja usada na heuristica
			LowerB = - heuristica->ResolverHeuristica();
			//LowerB = - 963576.12369;

			// qdo a soluçao é infinito passar o Inf<HpNum>()!!!
			// se o problema persistir n passar o lower bound (em nenhuma iteraçao) para o bundle!!
			if ( LowerB <= - GRB_INFINITY)
				LowerB = - Inf<HpNum>();

			// timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			if( Hrstt )
				Hrstt->Stop();
			return( LowerB );
		}
		else
			return( - Inf<HpNum>() );
	}
	else
		return( - Inf<HpNum>() );
}

FiOracle::FiStatus EspFiOracle::GetFiStatus( Index wFi )
{
	// Chamada para wFi = Inf e depois para cada componente wFi = 1,...,GetNrFi()
	// called at every iteration after that the master problem has been solved but before SetLambda() and Fi() are
	return( kFiNorm );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
	
void EspFiOracle::Deleted( cIndex Name )
{
	ptr_x_names[Name].clear();
	//vetorfloat().swap(ptr_x_names[Name]);
}

void EspFiOracle::Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm )
{
	throw NDOException( "Aggregate: Solucao convexificada nao agregada!!" );
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

EspFiOracle::~EspFiOracle()
{
 delete[] New_wFi;
 delete[] SGBse1;

 delete subproblemaHA;
 delete subproblemaHE;
 delete subproblemaT;
 delete subproblemaD;
 //delete heuristica;

 delete Hrstt;
 }
/*--------------------------------------------------------------------------*/
/*--------------------------------- OTHER ----------------------------------*/
/*--------------------------------------------------------------------------*/

void EspFiOracle::ConvSol(  )
{
	// Calcula a soluçao convexificada e armazena na classe ResultadosConj (ptr_x_til)
	vetorfloat * ptr_x_a = resultadosGurobi->GetPtrXtil();		//[ptr_x_til] => numero de componentes x numero de var. em cada componente 
	vetorfloat * ptr_x_a_hat = resultadosGurobi->GetPtrXhat();		//[ptr_x_hat] => deve ser preenchido no caso de easy component
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
		//else
		//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
	}
	else	// modelo desagregado: calculo para cada componente separado, pois algumas podem ser easy components
	{
		if (EasyComp == 0 || EasyComp == 2)		
		{//Hidraulico not easy
			for (int comp = 0; comp < resultadosGurobi->GetCH(); comp++)
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
				//else
				//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
			}
		}
		else
		{//Hidraulico easy
			for (int comp = 0; comp < resultadosGurobi->GetCH(); comp++)
			{
				mult = Slvr->ReadMult( I, D, comp + 1 );
				if (D != 0)
				{
					if( I )
					{
						ptr_x_a_hat[comp] = vetorfloat(ptr_x_a_hat[comp].size(), 0);	// preencher tds elementos de ptr_x_hat com zero, pois o mult é esparso!!
						for( Index i = 0; i < D; i++ )		// read  x[ wFi ]*
						{
							ptr_x_a[comp][I[ i ]] = mult[i];
							ptr_x_a_hat[comp][I[ i ]] = mult[i];
						}
					}
					else
						for( Index i = 0; i < D; i++ )
						{
							ptr_x_a[comp][ i ] = mult[i];
							ptr_x_a_hat[comp][ i ] = mult[i];
						}
				}
				//else
				//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
			}
		}
		for (int comp = resultadosGurobi->GetCH(); comp < resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size(); comp++)	
		{//Hidreletrico
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
			//else
			//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
		}
		for (int comp = resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size(); comp < resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size(); comp++)
		{//Termeletrico
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
			//else
			//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
		}
		if (EasyComp == 0 || EasyComp == 1)		
		{//Demanda not easy
			for (int comp = resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size(); comp < NumPCmp; comp++)
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
				//else
				//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
			}
		}
		else
		{//Demanda easy
			for (int comp = resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size(); comp < NumPCmp; comp++)
			{
				mult = Slvr->ReadMult( I, D, comp + 1 );
				if (D != 0)
				{
					if( I )
					{
						ptr_x_a_hat[comp] = vetorfloat(ptr_x_a_hat[comp].size(), 0);
						for( Index i = 0; i < D; i++ )		// read  x[ wFi ]*
						{
							ptr_x_a[comp][I[ i ]] = mult[i];
							ptr_x_a_hat[comp][I[ i ]] = mult[i];
						}
					}
					else
						for( Index i = 0; i < D; i++ )
						{
							ptr_x_a[comp][ i ] = mult[i];
							ptr_x_a_hat[comp][ i ] = mult[i];
						}
				}
				//else
				//	throw( NDOException( "EspFiOracle::ConvSol: D == 0" ) );
			}
		}
	}

	ptr_x_a = NULL;
	ptr_x_a_hat = NULL;
	delete ptr_x_a;
 }

HpNum EspFiOracle::GetHeuristicSolution(  )
{
	return heuristica->GetFO();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------ End File EspFiOracle.cpp ------------------------*/
/*--------------------------------------------------------------------------*/



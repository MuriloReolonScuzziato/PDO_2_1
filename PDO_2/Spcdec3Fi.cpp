/*--------------------------------------------------------------------------*/
/*-------------------------- File Spcdec3Fi.cpp --------------------------*/
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

#include "Spcdec3Fi.h"

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

Spcdec3Fi::Spcdec3Fi( CSistema * const sistema_end , Spcdec3Results * const resultadosGurobi_end, int tipo_dual)
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
	flag3 = int (sistema_a->GetFlagPhmax());
	flag4;
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());

	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	n = T * ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size() - 1 + (4+flag3)*sistema_a->hidreletricasVtr.size() + 3*JJ + sistema_a->barrasVtr.size()) + flag2 * sistema_a->hidreletricasVtr.size() * sistema_a->GetNCenarios();
	na = T * (sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size());
	nd = T * (sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size() + (1 - flag3));		// numero de variáveis duais (e n de variaveis auxiliares/duplicadas)
	// [Lpt Lph Lv Ld (Lphmax || Lres)] x t
	if ( flag1 == 0 )		// Remover teta, def e deixar uma variável de deficit por período
		n -= T * (sistema_a->barrasVtr.size() - 1 + sistema_a->barrasVtr.size() - 1);
	else if ( flag1 >= 2)	// Remover somente teta
		n -= T * (sistema_a->barrasVtr.size() - 1);
	
	sistema_a->SetPreconditioners(GenPreconditioners());		// Define os condicionadores
	x.resize(n);
	xa.resize(na);
	Lambda.resize(nd);
	subgrad.resize(nd);
	
	nt = ((n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size()) / T);
	nat = (na / T);
	
	// Inicialização para o Bundle:
	NumVar = nd;		// Numero de multiplicadores de Lagrange
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

	// lambda [pt ph v d (phmax)]
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

	//Fi5.resize(5);		// Media movel, valores zerados
	//Fi10.resize(10);
	//Fi20.resize(20);

	iter_count = 0;
	t_init = 0;

	SHAt = 0;
	SHEt = 0;
	STt = 0;
	SDt = 0;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void Spcdec3Fi::SetData( cLMRow UB  )
{
	UpBnd = UB;    
	//UnCnstr = UC;
 }

void Spcdec3Fi::CriarModelos( void )
{
	// criar modelos que sao easy components de forma difente, sem o modelo do gurobi, somente com as informaçoes que serao passadas para o MP!!!
	if (EasyComp == 1 || EasyComp == 3)		// Criado como easy component
		subproblemaHA = new Spcdec3SPHA(sistema_a);
	else
		subproblemaHA = new Spcdec3SPHA(sistema_a, ambGRB);
	subproblemaHE = new Spcdec3SPHE(sistema_a, ambGRB);
	subproblemaT = new Spcdec3SPT(sistema_a, ambGRB);
	if (EasyComp == 2 || EasyComp == 3)		// Criado como easy component
		subproblemaD = new Spcdec3SPD(sistema_a);
	else
		subproblemaD = new Spcdec3SPD(sistema_a, ambGRB);

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

void Spcdec3Fi::SetHeuristic ( Hrstc *heuristic )
{
	heuristica = heuristic;
	//heuristica = new Hrstc(sistema_a, resultadosGurobi, tempo_max);
}

void Spcdec3Fi::SetSubpPrecision( double precision )
{
	if (EasyComp == 1 || EasyComp == 3)		// Criado como easy component
	{} // precisao definida no solver do problema dual
	else
		subproblemaHA->SetPrecision( precision );
	subproblemaHE->SetPrecision( precision );
	subproblemaT->SetPrecision( precision );
	if (EasyComp == 2 || EasyComp == 3)		// Criado como easy component
	{} // precisao definida no solver do problema dual
	else
		subproblemaD->SetPrecision( precision );
}

void Spcdec3Fi::SetSubpTimes()
{
	tempoS.resize(4);
	SHAt = new OPTtimers();
	SHEt = new OPTtimers();
	STt = new OPTtimers();
	SDt = new OPTtimers();
}

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

bool Spcdec3Fi::GetUC( cIndex i )
{
	// retorna true para os multiplicadores que forem irrestritos em sinal (restrição de igualdade relaxada) e false para multiplicadores positivos (restrição de desigualdade relaxada)
	if ( (flag3 == 0) && (i >= sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size()) && ((i - sistema_a->termeletricasVtr.size() - 3*sistema_a->hidreletricasVtr.size()) % (sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + 1) == 0) )
		return( false );
	else
		return( true );
}

LMNum Spcdec3Fi::GetUB( cIndex i )
{
	//if ( (flag3 == 0) && (i >= sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size()) && ((i - sistema_a->termeletricasVtr.size() - 3*sistema_a->hidreletricasVtr.size()) % (sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + 1) == 0) )
	//	return( 10000 );
	//else
		return( Inf<LMNum>() );
}

// somente quando tem-se easy components

Index Spcdec3Fi::GetBNC( cIndex wFi )
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

Index Spcdec3Fi::GetBNR( cIndex wFi )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetNRestricoes(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetNRestricoes(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	throw( NDOException( "Spcdec3Fi::GetBNR: this component is not easy" ) );
	//return( 0 );
}

Index Spcdec3Fi::GetBNZ( cIndex wFi )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetNNZsRest(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetNNZsRest(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	throw( NDOException( "Spcdec3Fi::GetBNZ: this component is not easy" ) );
	//return( 0 );
}

void Spcdec3Fi::GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd )
{
	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			subproblemaHA->GetBDesc( wFi - 1, Bbeg, Bind, Bval, lhs, rhs, cst, lbd, ubd);
		else if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			subproblemaD->GetBDesc(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()), Bbeg, Bind, Bval, lhs, rhs, cst, lbd, ubd);
		else	// difficult component
			throw( NDOException( "Spcdec3Fi::GetBDesc: this component is not easy" ) );
	}
}

Index Spcdec3Fi::GetANZ( cIndex wFi, cIndex strt, Index stp )
{
	// Considerando que strt = 0 e stp = Inf<Index>() sempre.
	if ( strt != 0 || stp != Inf<Index>() )
		throw( NDOException( "Spcdec3Fi::GetANZ: dynamic generation of multipliers not implemented" ) );

	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			return ( subproblemaHA->GetANZ(wFi - 1) );
		if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			return ( subproblemaD->GetANZ(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size())) );
	}
	// difficult component
	throw( NDOException( "Spcdec3Fi::GetANZ: this component is not easy" ) );
	//return( 0 );
}

void Spcdec3Fi::GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval , cIndex strt, Index stp )
{
	// Considerando que strt = 0 e stp = Inf<Index>() sempre.
	if ( strt != 0 || stp != Inf<Index>() )
		throw( NDOException( "Spcdec3Fi::GetANZ: dynamic generation of multipliers not implemented" ) );

	if ( EasyComp != 0 )	// easy component
	{
		if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico: Confere se o subp. é easy e se o componente wFi corresponde a ele
			subproblemaHA->GetADesc(wFi - 1, Abeg, Aind, Aval);
		else if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp) )		//Demanda
			subproblemaD->GetADesc(wFi - 1 - (resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()), Abeg, Aind, Aval);
		else	// difficult component
			throw( NDOException( "Spcdec3Fi::GetANZ: this component is not easy" ) );
	}
}

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/

void Spcdec3Fi::SetLambda( cLMRow Lmbd )
{
	Lam1 = Lmbd;
	for( Index i = 0 ; i < NumHCmp ; i++ )
		New_wFi[ i ] = true;
 }

/*--------------------------------------------------------------------------*/

void Spcdec3Fi::SetLamBase( cIndex_Set LmbdB , cIndex LmbdBD )
{
	LamB = LmbdB;
	LamBd = LmbdBD;
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

HpNum Spcdec3Fi::Fi( cIndex wFi )
{
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

	// lambda [pt ph v d (phmax)]
	// Subproblemas HA, HE, T, D
	if( wFi == 0 )
	{
		if (flag3 == 1)
		{
			if( Fit )
				Fit->Stop();
			return( 0 );
		}
		else
		{
			double comp0 = 0;
			for (int t = 0; t < T; t++)
				comp0 += lambda[sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + t*(sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + 1)] * sistema_a->GetReserva(t);
			if( Fit )
				Fit->Stop();
			return( - comp0 );
		}
	}
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
		if (SHAt)
			SHAt->Start();
		subproblemaHA->ResolverProblemaRL(resultadosGurobi, &lambda, 0, wFi - 1);
		if (SHAt)
			SHAt->Stop();
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else if ( (wFi > resultadosGurobi->GetCH()) && (wFi <= resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size()) )		//Hidreletrico
	{
		int r = (wFi - resultadosGurobi->GetCH() - 1) % int (sistema_a->hidreletricasVtr.size());
		int no = (wFi - resultadosGurobi->GetCH() - 1) / int (sistema_a->hidreletricasVtr.size());
		if (SHEt)
			SHEt->Start();
		subproblemaHE->ResolverProblemaRL(resultadosGurobi, &lambda, 0, r, no);
		if (SHEt)
			SHEt->Stop();
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else if ( (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size()) && (wFi <= resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) )		//Termeletrico
	{
		int i = wFi - resultadosGurobi->GetCH() - resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() - 1;
		if (STt)
			STt->Start();
		subproblemaT->ResolverProblemaRL(resultadosGurobi, &lambda, 0, i);
		if (STt)
			STt->Stop();
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else if ( (EasyComp == 0 || EasyComp == 1) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp ) )		//Demanda
	{
		int no = wFi - resultadosGurobi->GetCH() - resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() - sistema_a->termeletricasVtr.size() - 1;
		if (SDt)
			SDt->Start();
		subproblemaD->ResolverProblemaRL(resultadosGurobi, &lambda, 0, no);
		if (SDt)
			SDt->Stop();
		temp = - resultadosGurobi->GetobjComp(wFi - 1);
	}
	else
	{
		throw NDOException( "Spcdec3Fi: Fi called fon an easy component." );	// Essa funçao nao pode ser chamada para easy components
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

bool Spcdec3Fi::NewGi( cIndex wFi )
{
	last_wFi = wFi;

	if( wFi == 0 )
		//if (flag3 == 1)
		throw NDOException( "NewGi: wFi = 0" );
		//else
		//	return( false );

	if ( (EasyComp == 1 || EasyComp == 3) && (wFi > 0) && (wFi <= resultadosGurobi->GetCH()) )	//Hidraulico ( 1 <= comp <= CH )
		throw NDOException( "Spcdec3Fi: NewGi called fon an easy component." );
	if ( (EasyComp == 2 || EasyComp == 3) && (wFi > resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size()) && (wFi <= NumPCmp ) )		//Demanda
		throw NDOException( "Spcdec3Fi: NewGi called fon an easy component." );
	
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

Index Spcdec3Fi::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name , cIndex strt , Index stp )
{
	// timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if( Fit )
		Fit->Start();

	if( Name < MBN )
		throw NDOException( "GetGi: past information is not recorded" );

	if( stp > NumVar )
		stp = NumVar;

	Index temp;
	if( Name == MBN )
	{
		// Name == MaxName. The required information is about the constant
		// subgradient of the linear 0-th component of Fi  - - - - - - - - - - - - -
		// like last_wFi == 0
		if (flag3 == 1)
		{
			SGBse = &InINF;
			if( Fit )
				Fit->Stop();
			return( 0 );
		}
		else
		{
			for (int i = 0; i < nd; i++)
				SubG[i] = 0;
			for (int t = 0; t < T; t++)
				SubG[sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + t*(sistema_a->termeletricasVtr.size() + 3*sistema_a->hidreletricasVtr.size() + 1)] = - sistema_a->GetReserva(t);
			// turn SubG from a "dense" NumVar-vector to a "sparse" one - - - - - - - - -
			Index_Set SGBse2 = Sparsify( SubG , SGBse1 , NumVar );
			*SGBse2 = Inf<Index>();
			temp = SGBse2 - SGBse1;
			SGBse = SGBse1;
		}
	}
	else
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

HpNum Spcdec3Fi::GetVal( cIndex Name )
{
 if( Name < MBN )
  throw( NDOException( "GetVal: past information is not recorded" ) );

 return( 0 );
 }

/*--------------------------------------------------------------------------*/

void Spcdec3Fi::SetGiName( cIndex Name )
{
	if ( !ptr_x_names[Name].empty() )
		ptr_x_names[Name].clear();
		//vetorfloat().swap(ptr_x_names[Name]);

	ptr_x_names[Name] = resultadosGurobi->GetCompX(last_wFi - 1);
	//if ( last_wFi - 1 == resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + 6)
	//	cout << " Name =  " << Name << endl;
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

HpNum Spcdec3Fi::GetLowerBound( cIndex wFi )
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

FiOracle::FiStatus Spcdec3Fi::GetFiStatus( Index wFi )
{
	// Fazer calculo da medida de inviabilidade e imprimir em arquivo auxiliar...
	// <iter> <Fi> <Sigma> <t> <n1> <n2> <nInf>
	if (wFi > GetNrFi())
	{
		// Calcula a solução convexificada
		//cout << "GetFiStatus iter= " << Slvr->NrIter() << endl;
		if (Slvr->NrIter() != 0)
			ConvSol();
		///resultadosGurobi->ArmazenarXtilhat();
		// Escreve cada iteração no log auxiliar
		// <iter> <Fi> <Sigma> <t> <n1> <n2> <nInf>
		// <iter> <Fi> <Sigma> <t> <ñ1> <ñ2> <ñInf> <hrst> <n1> <n2> <nInf>
		// ñ é a norma da solução convexificada e n é a norma da solução da heuristica

		vetorfloat norma(3,0);	// n1, n2 e nINF para todas as restrições (independente do tipo)
		resultadosGurobi->CalcularNormaSubgXtil(norma);

		// Dependendo do modelo dual a norma (e o desvio) utilizada como criterio de para é diferente
		double * norm;
		double desvio;
		if (Aggrgtd)
		{
			norm = &norma[1];
			desvio = 0.1*pow(NumVar, 0.5);
		}
		else
		{
			norm = &norma[0];
			HpNum epsB;
			Slvr->GetPar(Bundle::kEpsLin, epsB);
			desvio = epsB*1e3*NumVar;	// com Epslin = 1e-4 tem-se 0.1*NumVar
			//desvio = 0.1*NumVar;
		}

		//Escrever Log
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
		if (iter_count == 1)
		{
			// na interação 0 os suproblemas não são resolvido, portanto x = 0. Assim t_init deve ser calculado na iteração 1
			//t_init = (1 + abs(fdual_cont))/(5*norma[0]);		// fdual_cont é o valor do prob. extendido continuo
			t_init = (1 + abs(-Slvr->ReadFiVal()))/(5*(*norm));	// aqui então utilização o valor da função dual avaliada para lambda da relax. continua

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
			///tol = desvio/(1 + abs(fdual_cont));
			//tol = 1e-5;		// justificada pelo erro aceitável em sigma e tb em t*g, na casa das dezenas
			if (t_init < 1)
				tol = t_init*desvio/(1 + abs(fdual_cont));
			else
				tol = desvio/(1 + abs(fdual_cont));

			Slvr->SetPar(Bundle::ktStar, 50*t_init);		// altera tStar mas no log fica o valor inicial, tStar = [10,100]*tinit

			Slvr->SetPar(Bundle::ktCurr, t_init);
			//Slvr->SetPar(Bundle::ktInit, t_init);
			//Slvr->SetPar(NDOSolver::ktStar, t_init*10);

			//*log_auxiliar << "Tol:" << char(9) << tol << endl;

			iter_count++;
			return ( kFiChgd );
		}
		// Calcula o t inicial da primeira iteração
		if (Slvr->NrIter() == 0)
			iter_count++;
			
		// Criterio de parada pelo valor relativo do sigma-subgradente e do sigma
		else if ( (abs(BSlvr->ReadSigma()/(1 - Slvr->ReadFiVal())) <= tol) && (abs(BSlvr->Readt()*(*norm)/(1 - Slvr->ReadFiVal())) <= tol) )
		//else if ( (abs(BSlvr->ReadSigma()/(1 - Slvr->ReadFiVal())) <= tol) && (abs((*norm)/(1 - Slvr->ReadFiVal())) <= tol) )
			return( kFiStop );

		//// Criterio de parada pelo gap entre Fi e o valor ótimo (comentada em 10/01/2018)
		//else if (((-LowerB + Slvr->ReadFiVal()) / (-Slvr->ReadFiVal()) <= 1) && (heuristica != NULL))
		//	return(kFiStop);

		else
			return( kFiCont );
			//return( kFiNorm );	//Se usar isso ele vai verificar a condição do Feixes???


		//return( kFiNorm );


		// If you don't use a long-term t-strat *and* you have an alternative stopping rule based on the norm of violation, tStar is ineffective. Well, almost: you have to be sure it's large enough so that the standard stopping condition does not stop you early on (but even for this, you can still use kContAnyway in GetFiStatus() to avoid the problem).
		// posso conferir se kFiNorm convergiria e colocar uma mensagem, porém continuar iterando...


		//  Relaxing Ax >= b then the unfeasibility is [ Ax - b ]_-, meaning that each positive entry of the vector Ax - b is set to zero (in D^*). In practice, the unfeasibility is zero if all contraints are satisfied (obviously).
	}

	// Chamada para wFi = Inf e depois para cada componente wFi = 1,...,GetNrFi()
	// called at every iteration after that the master problem has been solved but before SetLambda() and Fi() are
	return( kFiNorm );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
	
void Spcdec3Fi::Deleted( cIndex Name )
{
	if (Name < MBN)
		ptr_x_names[Name].clear();
	//vetorfloat().swap(ptr_x_names[Name]);
}

void Spcdec3Fi::Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm )
{
	throw NDOException( "Aggregate: Solucao convexificada nao agregada!!" );
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

Spcdec3Fi::~Spcdec3Fi()
{
	BSlvr = NULL;
	heuristica = NULL;
	log_auxiliar = NULL;
	resultadosGurobi = NULL;
	sistema_a = NULL;

	delete[] New_wFi;
	delete[] SGBse1;

	delete subproblemaHA;
	delete subproblemaHE;
	delete subproblemaT;
	delete subproblemaD;

	delete Hrstt;
	//delete BSlvr;
	ambGRB.~GRBEnv();
 }
/*--------------------------------------------------------------------------*/
/*--------------------------------- OTHER ----------------------------------*/
/*--------------------------------------------------------------------------*/

void Spcdec3Fi::ConvSol(  )
{
	// Calcula a soluçao convexificada e armazena na classe Spcdec3Results (ptr_x_til)
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
		else
			if (Slvr->ReadLBMult() != 1)
				throw( NDOException( "Spcdec3Fi::ConvSol, agg: D == 0" ) );
		// Since the sum of that and of all the others has to be 1 as well, it means that all the others are 0. Hence, D = 0.
		// If this is what happens, it is good news: it means that solution is integer optimal.
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
				else
					if (Slvr->ReadLBMult() != 1)
						throw( NDOException( "Spcdec3Fi::ConvSol, HA: D == 0" ) );
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
				else
					if (Slvr->ReadLBMult() != 1)
						throw( NDOException( "Spcdec3Fi::ConvSol, HAe: D == 0" ) );
					//cout << "ReadLBMult() = " << Slvr->ReadLBMult() << endl;
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
			else
				if (Slvr->ReadLBMult() != 1)
					throw( NDOException( "Spcdec3Fi::ConvSol, HE: D == 0" ) );
		}
		for (int comp = resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size(); comp < resultadosGurobi->GetCH() + resultadosGurobi->GetN() * sistema_a->hidreletricasVtr.size() + sistema_a->termeletricasVtr.size(); comp++)
		{//Termeletrico
			mult = Slvr->ReadMult( I, D, comp + 1 );
			if (D != 0)
			{
				if( I )
				{
					for (size_t j = 0; j < ptr_x_names[ I[0] ].size(); j++)		// loop nas variáveis de cada componente
					{
						ptr_x_a[comp][j] = 0;
						for( Index i = 0; i < D; i++ )		// loop nos multiplicadores (names)
							ptr_x_a[comp][j] += mult[i] * ptr_x_names[ I[ i ] ][j];
					}
				}
				else
				{
					for (size_t j = 0; j < ptr_x_names[ 0 ].size(); j++)
					{
						ptr_x_a[comp][j] = 0;
						for( Index i = 0; i < D; i++ )
							ptr_x_a[comp][j] += mult[i] * ptr_x_names[ i ][j];
					}
				}
			}
			else
				if (Slvr->ReadLBMult() != 1)
					throw( NDOException( "Spcdec3Fi::ConvSol, T: D == 0" ) );
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
				else
					if (Slvr->ReadLBMult() != 1)
						throw( NDOException( "Spcdec3Fi::ConvSol, D: D == 0" ) );
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
				else
					if (Slvr->ReadLBMult() != 1)
						throw( NDOException( "Spcdec3Fi::ConvSol, De: D == 0" ) );
			}
		}
	}
	
	// Adicionar componentes do Lower Bound!!
	if (Slvr->ReadLBMult() != 0)
	{
		HpNum multiLB = Slvr->ReadLBMult();
		resultadosGurobi->AdicionarLBXtil(heuristica->GetX(), multiLB);
	}

	ptr_x_a = NULL;
	ptr_x_a_hat = NULL;
	delete ptr_x_a;
	delete ptr_x_a_hat;
 }

HpNum Spcdec3Fi::GetHeuristicSolution(  )
{
	return heuristica->GetFO();
}

vetorfloat Spcdec3Fi::GenPreconditioners()
{
	vetorfloat precond;
	// L = [Lpt Lph Lv Ld (Lphmax/Lreserva)]
	int flag3 = int (sistema_a->GetFlagPhmax());
	//int deltaa = 0;
	// Kind of Preconditioner
	switch (sistema_a->GetFlagPrecond())
	{
	case 1:		//probability of the node
		{
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size() + (1 - flag3); i++)
					precond.push_back(1);
			}
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (int t = 0; t < sistema_a->GetTt2() - sistema_a->GetTt1(); t++)
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size() + (1 - flag3); i++)
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));
				}
			}
			break;
		}
	case 2:		//square root of the probability of the node
		{
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
					precond.push_back(1);
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
					precond.push_back(1);
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					precond.push_back(1);
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
					precond.push_back(1);
				if (flag3 == 1)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
						precond.push_back(1);
				}
				else
					precond.push_back(1);								// reserva
			}
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (int t = 0; t < sistema_a->GetTt2() - sistema_a->GetTt1(); t++)
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));
					if (flag3 == 1)
					{
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
							precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));
					}
					else
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)));	// reserva
				}
			}
			break;
		}
	case 3:		//maximum value of the relaxed constraint
		{
			double phmax;
			double qhmax;
			for (int t = 0; t < T; t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				{
					if (sistema_a->GetFlagTbinaryModel() == 1)
						precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					else
						precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
				{
					phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					precond.push_back(1 / (phmax - 0));
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					precond.push_back(1 / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
				{
					qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						qhmax += (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
					precond.push_back(1 / (qhmax + sistema_a->hidreletricasVtr[r].GetSmax() - 0));
				}
				if (flag3 == 1)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
					{
						phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(1 / (phmax - 0));
					}
				}
				else
				{
					phmax = 0;
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//reserva
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					precond.push_back(1 / (phmax - sistema_a->GetReserva(t)));
				}
			}
			break;
		}
	case 4:		//combinaion between 1 and 3
		{
			double phmax;
			double qhmax;
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				{
					if (sistema_a->GetFlagTbinaryModel() == 1)
						precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					else
						precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
				{
					phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					precond.push_back(1 / (phmax - 0));
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					precond.push_back(1 / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
				{
					qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						qhmax += (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
					precond.push_back(1 / (qhmax + sistema_a->hidreletricasVtr[r].GetSmax() - 0));
				}
				if (flag3 == 1)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
					{
						phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(1 / (phmax - 0));
					}
				}
				else
				{
					phmax = 0;
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//reserva
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					precond.push_back(1 / (phmax - sistema_a->GetReserva(t)));
				}
			}
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (int t = 0; t < sistema_a->GetTt2() - sistema_a->GetTt1(); t++)
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
					{
						if (sistema_a->GetFlagTbinaryModel() == 1)
							precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
						else
							precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
					}
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
					{
						phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (phmax - 0));
					}
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
					{
						qhmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							qhmax += (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (qhmax + sistema_a->hidreletricasVtr[r].GetSmax() - 0));
					}
					if (flag3 == 1)
					{
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
						{
							phmax = 0;
							for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
								phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
							precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (phmax - 0));
						}
					}
					else
					{
						phmax = 0;
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//reserva
							for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
								phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) / (phmax - sistema_a->GetReserva(t)));
					}
				}
			}
			break;
		}
	case 5:		//combination between 2 and 3
		{
			double phmax;
			double qhmax;
			for (int t = 0; t < sistema_a->GetTt1(); t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
				{
					if (sistema_a->GetFlagTbinaryModel() == 1)
						precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
					else
						precond.push_back(1 / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
				{
					phmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					precond.push_back(1 / (phmax - 0));
				}
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
					precond.push_back(1 / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
				for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
				{
					qhmax = 0;
					for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
						qhmax += (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
					precond.push_back(1 / (qhmax + sistema_a->hidreletricasVtr[r].GetSmax() - 0));
				}
				if (flag3 == 1)
				{
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
					{
						phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(1 / (phmax - 0));
					}
				}
				else
				{
					phmax = 0;
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//reserva
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
					precond.push_back(1 / (phmax - sistema_a->GetReserva(t)));
				}
			}
			for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			{
				for (int t = 0; t < sistema_a->GetTt2() - sistema_a->GetTt1(); t++)
				{
					for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)	//pt
					{
						if (sistema_a->GetFlagTbinaryModel() == 1)
							precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (sistema_a->termeletricasVtr[i].GetPmax() - sistema_a->termeletricasVtr[i].GetPmin()));
						else
							precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (sistema_a->termeletricasVtr[i].GetPmax() - 0));
					}
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//ph
					{
						phmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (phmax - 0));
					}
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//v
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (sistema_a->hidreletricasVtr[r].GetVmax() - sistema_a->hidreletricasVtr[r].GetVmin()));
					for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//d
					{
						qhmax = 0;
						for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
							qhmax += (sistema_a->hidreletricasVtr[r].grupoVtr[j].GetQmax() * sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (qhmax + sistema_a->hidreletricasVtr[r].GetSmax() - 0));
					}
					if (flag3 == 1)
					{
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//phmax
						{
							phmax = 0;
							for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
								phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
							precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (phmax - 0));
						}
					}
					else
					{
						phmax = 0;
						for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)	//reserva
							for (int j = 0; j < sistema_a->hidreletricasVtr[r].GetNGrupos(); j++)
								phmax += sistema_a->hidreletricasVtr[r].grupoVtr[j].GetPmax()*sistema_a->hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
						precond.push_back(sqrt(sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen)) / (phmax - sistema_a->GetReserva(t)));
					}
				}
			}
			break;
		}
	default:	//none of them
		{
			for (int t = 0; t < T; t++)
			{
				for (size_t i = 0; i < sistema_a->termeletricasVtr.size() + (3+flag3)*sistema_a->hidreletricasVtr.size() + (1 - flag3); i++)
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
/*------------------------ End File Spcdec3Fi.cpp ------------------------*/
/*--------------------------------------------------------------------------*/



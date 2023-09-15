/*--------------------------------------------------------------------------*/
/*---------------------------- File Scndec2Fi.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Oracle para a decomposiçao por cenários 2 (somente algumas váriaveis chave
 * são consideradas nas restrições de não antecipatividade: pt, u, phg, d)
 * interface defined by the abstract base class FiOracle.
 *
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// order the functions are called:
//	Solve Master Problem
//	GetFiStatus()
//	SetLambda()
//	Fi()
//	all the GetGi() [NewGi(), GetGi(), GetVal() and SetGiName(), respectively]
//	GetLowerBound()

#pragma once

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

#include "Scndec2SPC.h"
//#include "Scndec2Hrstc.h"

#include "Bundle.h"

typedef vector<Scndec2SPC *> vetorSPC;

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

//using namespace met_DecEsp;

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Scndec2Fi --------------------------*/
/*--------------------------------------------------------------------------*/

class Scndec2Fi : public FiOracle
{
	CSistema * sistema_a;
	Scndec2Results * resultadosGurobi;

	GRBEnv ambGRB;										// Ambiente Gurobi
	vetorSPC vtr_subproblemaC;							// Vetor com o subproblema de cada cenário
	//Scndec2Hrstc * heuristica;
	ofstream * log_auxiliar;

	double fdual_RL, norma_sg, tol;
	int n, na;											// numero de variaveis primais do problema e numero de var. duais (multiplicadores) por cenário
	vetorfloat x; //, Lambda, subgrad;
	int flag1, flag2, flag4, flag7, JJ;
	int T, nt;
	int iter_count;
	HpNum fdual_cont, t_init;

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

	Scndec2Fi( CSistema * const sistema_end , Scndec2Results * const resultadosGurobi_end, int tipo_dual);
   
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

	void SetData ( cLMRow UB  );
	
	void SetMaxName( cIndex MxNme = 0 )
	{
		MBN = MxNme;
		ptr_x_names.resize(MBN);
		epsilon_names.resize(MBN);
		comp_names.resize(MBN);
    }

	void CriarModelos ( void );

	//void SetHeuristic ( Scndec2Hrstc *heuristic );

	void SetFiIterLog( ofstream *outs = 0 )
	{
		log_auxiliar = outs;		// Criar arquivo de log auxiliar
    }

	void SetSubpPrecision( double precision );		// Define precisão dos subproblemas

	void SetNDOSolver( NDOSolver *NwSlvr = 0 )
	{
		Slvr = NwSlvr;
		BSlvr = dynamic_cast<Bundle*>( Slvr );
    }

	void SetFiCont( HpNum Fic = 0 )
	{
		fdual_cont = Fic;
    }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

	Index GetNumVar( void ) const
	{
		return( NumVar );
	}

	Index GetNrFi( void ) const
	{
		return( NumQCmp + NumPCmp );	   
	}

	bool GetUC( cIndex i );

	LMNum GetUB( cIndex i );

	cLMRow GetLambda()
	{
		return Lam1;
	}

	// Easy component only:
	//Index GetBNC( cIndex wFi );

	//Index GetBNR( cIndex wFi );

	//Index GetBNZ( cIndex wFi );

	//void GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval , double *lhs , double *rhs , double *cst , double *lbd , double *ubd );
   
	//Index GetANZ( cIndex wFi, cIndex strt = 0 , Index stp = Inf<Index>() );
   
	//void GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval , cIndex strt = 0 , Index stp = Inf<Index>() );
   
/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/

	void SetLambda( cLMRow Lmbd = 0 );

	void SetLamBase( cIndex_Set LmbdB = 0 , cIndex LmbdBD = 0 );

	void SetLowerBound( cHpNum FiLower = - Inf<HpNum>() )
	{
		LowerB = FiLower;
	}

	bool SetPrecision( HpNum Eps );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   HpNum Fi( cIndex wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   bool NewGi( cIndex wFi = Inf<Index>() );

   Index GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name = Inf<Index>() ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

   HpNum GetVal( cIndex Name = Inf<Index>() );

   void SetGiName( cIndex Name );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

	HpNum GetLowerBound( cIndex wFi = Inf<Index>() );

	FiStatus GetFiStatus( Index wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
	
	void Deleted( cIndex Name = Inf<Index>() );

	void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~Scndec2Fi();

/* Destructor of the class: it must be virtual. */

/*--------------------------------------------------------------------------*/
/*--------------------------------- OTHER ----------------------------------*/
/*--------------------------------------------------------------------------*/

	void ConvSol( );
	/* Compute the convexified solution. */

	HpNum GetHeuristicSolution( );		// return the optimal value found in the heuristic

	void SetHrstTime( const bool TimeIt = true )
	{
		if( TimeIt )
			if( Hrstt )
				Hrstt->ReSet();
			else
				Hrstt = new OPTtimers();
		else
		{
			delete Hrstt;
			Hrstt = 0;
		}
	}

	void HrstTime( double &t_us , double &t_ss )
	{
		t_us = t_ss = 0;
		if( Hrstt )
			Hrstt->Read( t_us , t_ss ); 
    }

	double HrstTime( void )
	{
		return( Hrstt ? Hrstt->Read() : 0 );
	}

	vetorfloat GenPreconditioners();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to extend the code by deriving a new class may use these   --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  Index NumVar;      // the number of variables

  HpNum LowerB;     // A (coarse) lower bound on the value of Fi.
  cLMRow UpBnd;      // the upper bounds
  //cBool_Vec UnCnstr; // the variable is unconstrained
  bool Aggrgtd;    /**< Agreggation flag.
					 Variable that defines if outputs aggregated or disaggregated information*/
  Index EasyComp;    //< Easy components type.

  Index_Set SGBse1;       // format of Subgradient
  
  Index NumQVar;     // Number of variable for the quadratic part
  Index NumQCmp;     // Number of quadratic component 
  Index NumPCmp;     // Number of component in the polyedral part
  Index NumECmp;     // Number of easy component
  Index NumHCmp;     // Number of hard component
  //Index_Set VarDiv;  // Index of the group of variables of each component
	 // 				 // component wFi involves group of variables j with
  //                   // VarDiv[ wFi - 1 ] <= j < VarDiv[ wFi ]
  //Index_Set VarI;	 // Index of the first variable of each group of variables
  //                   // VarI[ VarDiv[] ]
  //Index_Set VarF;	 // Index of the last variable of each group of variables
  //                   // VarF[ VarDiv[] ]
  
  Index MBN;         // max n. of names

  cLMRow Lam1;       // point at which to evaluate the function
  cIndex_Set LamB;   // vector of the indices of nonzeroes in Lam1
  Index LamBd;       // lenght of LamB

  Index last_wFi;    // remeber the argument of the last call to NewGi
  bool *New_wFi;     // New_wFi[ n ] is true if it is available a new
                     // subgradient for the n+1 component ? 
  					 // (should be n < NumQCmp)
  vetorfloat2 ptr_x_names;		// Primal solution for each "item" or name
  //vetorfloat * ptr_x_hat;		// Primal solution out of LR for each component
  //vetorfloat * ptr_x_til;		// Convexified primal solution for each component
  vetorfloat epsilon_names;		// Epsilon for each "item" or name
  vetorfloat comp_names;		// identifica a qual componente pertence o item name

  OPTtimers *Hrstt;      ///< OPTtimer for timing purposes

  Bundle *BSlvr;		// pointer to the Bundle

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline SgNum sign( SgNum x )
{
 if( x > 0 )
  return( 1 );
 else
  return( x < 0 ? -1 : 0 );
 }
 
/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 };  // end( class TestFi )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File Scndec2Fi.h -------------------------*/
/*--------------------------------------------------------------------------*/

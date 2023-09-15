/*--------------------------------------------------------------------------*/
/*------------------------- File OSIMPSolver.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the OSIMPSolver class, which implements a generic Master
 * Problem Solver for Bundle algorithms, using a generic OSISolverInterface
 * object. This class conforms to the interface defined by the class MPSolver
 * [see MPSolver.h].
 *
 * \version 1.30
 *
 * \date 21 - 11 - 2015
 *
 * \author Antonio Frangioni (original idea & implementation)\n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone (implementation)\n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' della Calabria \n
 *
 * \author Andrea Nerli (implementation)\n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2001 - 2015 by Antonio Frangioni.
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OSIMPSolver
 #define __OSIMPSolver /* self-identification: #endif at the end of the file*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define PRESERVE_OSI_SOLS 1

/* As the name implies, OSIMPSolver uses an object derived from
   OsiSolverInterface to solve the MP. Certain usage patterns cause the MP to
   be changed, and *after that* the previous optimal primal and dual solution
   to be accessed. Certain OsiSolvers graciously allow that, by keeping the
   previous solution even in face of changes until the MP is explicitly
   re-solved, while others are more strict and delete any previous solution
   as soon as anything changes in the MP.

   This macro, if set to 1, dictates that the OSIMPSolver object stores the
   solution information of the OsiSolver in its own data structures, so that
   it remains available even after changes in the MP. Doing so is safe, but
   it may be useless (and therefore better avoided for higher efficiency in
   both space and time) with certain OsiSolvers. */

/*--------------------------------------------------------------------------*/

#define OSIMPSOLVERLOG 0

/* If OSIMPSOLVERLOG > 0, the OSIMPSolver class produces a log of its
   activities on the ostream object and at the "level of verbosity" set with
   the method SetMPLog() [see below]. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MPSolver.h"

#include "OsiSolverInterface.hpp"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The OSIMPSolver class implements a generic Master Problem Solver for
    Bundle algorithms, using a generic OSISolverInterface object. This class
    conforms to the interface defined by the class MPSolver [see MPSolver.h].
  */

class OSIMPSolver : public MPSolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
   @{ */

  enum StabFun { unset = 0 , none , boxstep , quadratic };

/** Public enum for handling the Bundle stabilization term [see below]. */

  enum OsiAlg { kAuto = 0 , kPrim , kDual , kNet , kBar , kSif , kCon };

/** Public enum which is used to define the algorithm used in Solve() [see
    below] for solving either LP or QP problem:

    - kAuto:  automatic
    - kPrim:  Primal simplex
    - kDual:  Dual simplex
    - kNet:   Network simplex
    - kBar:   Barrier
    - kSif:   Sifting
    - kCon:   Concurrent (Dual, Barrier, and Primal) */

  enum OsiRed { rNo = 0 , rPrim , rDual , rBoth };

/** Public enum which is used to tell whether or not the reduction is used
    [see below]:

    - rNo:    no reduction
    - rPrim:  primal reduction
    - rDual:  dual reduction
    - rBoth:  primal and dual reduction */

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/

  OSIMPSolver( istream *iStrm = 0 );

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- methods from the base class ----------------------*/
/*--------------------------------------------------------------------------*/

  void SetDim( cIndex MxBSz = 0 , FiOracle *Oracle = 0 ,
	       const bool UsAvSt = false );

/*--------------------------------------------------------------------------*/

  void Sett( cHpNum tt = 1 );

/*--------------------------------------------------------------------------*/

   inline void SetPar( const int wp , cHpNum value );

/*--------------------------------------------------------------------------*/

  void SetLowerBound( cHpNum LwBnd = - Inf<HpNum>() ,
		      cIndex wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/

  void SetMPLog( ostream *outs = 0 , const char lvl = 0 );

/**< lvl controls the "level of verbosity" of the code. The first four bits
   of lvl have the following meaning:

    0  =>  no log at all (also assumed if log = 0);

    1  =>  some messages are printed;

    2  =>  as, 1 plus the current Master Problem is also saved as a MPS file
           just prior to being solved;

    3  ==>  as 2, plus very verbose information is also printed

    4 .. 15 unused, available to derived classes. */

/*--------------------------------------------------------------------------*/

  void SetAlgo( const OsiAlg algo = kAuto , const OsiRed reduc = rNo );

/**< Change "int" algorithmic parameters of the NDO solver. This method is
   used to set the method for solving either LP or QP. */

/*--------------------------------------------------------------------------*/
/*--------------------- derived-class-specific methods ---------------------*/
/*--------------------------------------------------------------------------*/

  void SetOsi( OsiSolverInterface *osi = 0 );

/**< Provides OSIMPSolver with an object of any class derived from
   OsiSolverInterface, which is used to actually solve the Master Problem.
   Note that the object becomes "property" of the OSIMPSolver, which will
   delete it when it is deleted. */

/*--------------------------------------------------------------------------*/

  void SetStabType( const StabFun sf = none );

/**< Sets the stabilizing term in the Master Problem. Possible values are:

   - none: no stabilization;

   - boxstep: a box in the infinity-norm in the primal, a 1-norm penalty
              in the dual;

   - quadratic: a 2-norm penalty term both in the primal and in the dual. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

  MPStatus SolveMP( void );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

  HpNum ReadFiBLambda( cIndex wFi = Inf<Index>() );

  HpNum ReadDt( cHpNum tt = 1 );

  HpNum ReadSigma( cIndex wFi = Inf<Index>() );

  HpNum ReadDStart( cHpNum tt = 1 );

  cLMRow Readd( bool Fulld = false );

  void ReadZ( LMRow tz , cIndex_Set &I , Index &D ,
	      cIndex wFi = Inf<Index>() );

  cHpRow ReadMult( cIndex_Set &I , Index &D , cIndex wFi = Inf<Index>() ,
		   const bool IncldCnst = true );

  HpNum ReadLBMult( cIndex wFi = Inf<Index>() );

  HpNum ReadGid( cIndex Nm = Inf<Index>() );

  void MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau );

  void SensitAnals( HpNum &lp , HpNum &cp );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

  Index BSize( cIndex wFi = Inf<Index>() );

  Index BCSize( cIndex wFi = Inf<Index>() );

  Index MaxName( cIndex wFi = Inf<Index>() );

  Index WComponent( cIndex i );

  bool IsSubG( cIndex i );

  Index NumNNVars( void );

  Index NumBxdVars( void );

  bool IsNN( cIndex i );

  void CheckIdentical( const bool Chk = true );
		
  cHpRow ReadLinErr( void );

  HpNum EpsilonD( void );

  bool FiBLambdaIsExact( cIndex wFi );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

  SgRow GetItem( cIndex wFi = Inf<Index>() );

  void SetItemBse( cIndex_Set SGBse = 0 , cIndex SGBDm = 0 );

  Index CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai , HpNum &ScPri );

  Index CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt );

  bool ChangesMPSol( void );

  void SetItem( cIndex Nm = Inf<Index>() );

  void SubstItem( cIndex Nm );

/*--------------------------------------------------------------------------*/

  void RmvItem( cIndex i );

  void RmvItems( void );

/*--------------------------------------------------------------------------*/

  void SetActvSt( cIndex_Set AVrs = 0 , cIndex AVDm = 0 );

  void AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs );

  void RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs );

/*--------------------------------------------------------------------------*/

  void AddVars( cIndex NNwVrs );

  void RmvVars( cIndex_Set whch = 0 , Index hwmny = 0 );

/*--------------------------------------------------------------------------*/

  void ChgAlfa( cHpRow DeltaAlfa );

  void ChgAlfa( cHpRow NewAlfa , cIndex wFi );

/*--------------------------------------------------------------------------*/

  void ChangeCurrPoint( cLMRow DLambda , cHpRow DFi );

  void ChangeCurrPoint( cHpNum Tau , cHpRow DFi );

/*--------------------------------------------------------------------------*/

  void ChgSubG( cIndex strt , Index stp , cIndex wFi );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  ~OSIMPSolver();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  void cleanup( void );

/*--------------------------------------------------------------------------*/

  void tUpdatePrices( cIndex strt = 0 , Index stp = Inf<Index>() );

  void ptUpdatePrices( cIndex strt = 0 , Index stp = Inf<Index>() );

  void ptUpdatePricesInPlace( void );

/* Update the prices after either t has changed or if the current point
   has been changed. The last form changes the prices in tmpHP[], assumed to
   hold the whole vector of cost coefficients of the master problem. */

/*--------------------------------------------------------------------------*/

  void switchToQP( void );

/* Private methods for handling the quadratic stabilization. */

/*--------------------------------------------------------------------------*/

//  inline bool isactive( Index i );
		
  inline void activate( Index i );

  inline void deactivate( Index i );

  inline void resizeHP( Index i );

  inline void resizeI( Index i );

/*--------------------------------------------------------------------------*/

  inline void CheckDS( void );

/* Placeholder method: it should do nothing, but depending on the value of
   the macro CHECK_DS [see OSIMPSolver.C] some extra checking is done in
   there for debugging purposes. */

/*--------------------------------------------------------------------------*/

  Index CheckBCopy( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  OsiSolverInterface* osiSlvr;  ///< pointer to the OsiSolver
  FiOracle *FIO;	        ///< pointer to the FiOracle

  OsiAlg algorithm;          ///< algorithm used in the solver

  bool useactiveset;	     ///< true if the active set are taken in account
  bool first;		     /**< true if it's the first master problem to
				solve */
  bool checkID;
  Bool_Vec weasy;            ///< which Fi-component is easy

  StabFun stab;	             ///< tell the stabilization term adopted
  HpNum t;		     ///< proximal parameter term

  double MaxTime;            ///< max running time for each call to SolveMP
  double OptEps;             ///< optimality tolerance
  double FsbEps;             ///< feasibility tolerance

  double *RhoCol;            ///< entries and vocabulary for rho
  int  *RhoColBse;
  Index	RhoColBDm;
		
  SgRow NewItem;	     ///< the new item
  Index NewItemFi;	     ///< the function component relative to new item
  cIndex_Set NewItemBse;     ///< index vector of the new item
  Index NewItemBDm;	     ///< dimension of new item

  Index MaxNZ;		     ///< max dimension of the item

  bool NewItemisSG;	     ///< subgradient vs constraint
  double NewItemprice;	     ///< price of the item
  double NewItemScPri;       ///< Gid product

  Index_Set NSubG;           ///< number of subgradients (per component)
  Index_Set NConst;          ///< number of constraints (per component)

  Index_Set comp_row;	     ///< row vocabulary
  Index_Set dict_item;	     ///< item vocabulary
  Index_Set dict_slack;	     ///< slack vocabulary
  Index_Set dict_stab;       ///< stabilization variables vocabulary

  Index_Set wcomp;           /**< Has several different uses:
				- if InINF, it means the item is not there
				- the value *masked by the two most
				  siggnificant bits* is the component of
				  the function that the item belongs to
				- the value of the most significant bit
				  is 1 if the item is a subgradient, 0 if
				  it is a constraint
				- the value of the second-most significant
				  bit is 1 if the item was *not* present at
				  the time that the last Master Problem has
				  been solved, 0 otherwise */

  Index_Set comp_col;	     ///< column vocabulary

  Index dict_item_maxname;   ///< current number of items

  HpRow tempHP;		     ///< buffer structures
  int*  tempI;
  Index tempHP_size;
  Index tempI_size;

  LMRow Upper;		     ///< upper bounds on the primal variable
  LMRow Lower;		     ///< lower bounds on the primal variable

  Index NNVars;              ///< variables with >= 0 constraints
  Index BxdVars;             ///< variables with >= 0 or <= UB constraints

  cIndex_Set Aset;           ///< active set structure
  Index Asetdim;	     ///< size of active set

  HpRow FiLambda;            ///< Fi[ k ]( Lambda ) for easy components
  HpRow FiLambda1;           ///< Fi[ k ]( Lambda1 ) for easy components

  bool HasLwrBnd;            ///< true if the global lower bound has been set
  HpRow LwrBnds;             ///< values of lower bounds for each component
  
  HpRow GiPerd;              /**< keep the scalar products between each
				item and the direction: while it is not used
				at all iterations, its computation uses
				information that may change (e.g. if Alphas
				change) so it's better to save it */
  #if( PRESERVE_OSI_SOLS )
   double *csol;             ///< column (primal) solution
   int csols;                ///< lenght of csol
   double *rsol;             ///< row (dual) solution
   int rsols;                ///< lenght of rsol
   double *rcst;             ///< row reduced cost
  #endif

  CoinMessageHandler *derhand;

/*--------------------------------------------------------------------------*/

 };  // end( class OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* OSIMPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OSIMPSolver.h ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- File OSIMPSolver.C ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Implementation of the OSIMPSolver class, which implements a generic  --*/
/*-- Master Problem Solver for Bundle algorithms, using a generic         --*/
/*-- OSISolverInterface object. This class conforms to the interface      --*/
/*-- defined by the class MPSolver [see MPSolver.h].                      --*/
/*--                                                                      --*/
/*--                            VERSION 1.56                              --*/
/*--                           05 - 10 - 2014                             --*/
/*--                                                                      --*/
/*--                            Original Idea:                            --*/
/*--                                                                      --*/
/*--                           Antonio Frangioni                          --*/
/*--                                                                      --*/
/*--                        Operations Research Group                     --*/
/*--                       Dipartimento di Informatica                    --*/
/*--                           Universita' di Pisa                        --*/
/*--                                                                      --*/
/*--                            Implementation:                           --*/
/*--                                                                      --*/
/*--                           Antonio Frangioni                          --*/
/*--                             Andrea Nerli                             --*/
/*--                                                                      --*/
/*--                        Operations Research Group                     --*/
/*--                       Dipartimento di Informatica                    --*/
/*--                           Universita' di Pisa                        --*/
/*--                                                                      --*/
/*--                             Enrico Gorgone                           --*/
/*--                                                                      --*/
/*--                        Operations Research Group                     --*/
/*--                       Dipartimento di Informatica                    --*/
/*--                        Universita' della Calabria                    --*/
/*--                                                                      --*/
/*--             Copyright (C) 2001 - 2014 by Antonio Frangioni           --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

#include "OSIMPSolver.h"

#include "OPTvect.h"

#include "NDOSlver.h"

#if( OSIMPSOLVERLOG )
 #include <strstream>
#endif

#if CPLEX == 1
 #include "OsiCpxSolverInterface.hpp"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OSIMPSOLVERLOG )
 #define MSG( l , m ) if( MPLLvl > l ) (*MPLog) << m << flush
#else
 #define MSG( l , m )
#endif

/*--------------------------------------------------------------------------*/

#define CHECK_DS 0

/* If CHECK_DS > 0, various data structures are checked for correctness
   during the run of the algorithm, tremendously slowing down the algorithm
   but allowing to debug the thing.
   What data structures are checked is coded bit-wise in CHECK_DS:

    bit 0 (+ 1)  =>  the various dict_* are tested
    bit 1 (+ 2)  =>  
    bit 2 (+ 4)  =>  
    */

/*--------------------------------------------------------------------------*/

// true if i was in the last Master Problem
#define WasInMP( i ) ( ! ( wcomp[ i ] & LLBIndex ) )

// the component of i, assuming wcomp[ i ] < InINF
#define WComp( i ) ( wcomp[ i ] & ~IMask )

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if OPT_USE_NAMESPACES
 using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

static cIndex InINF = Inf<Index>();
static cLMNum LMINF = Inf<LMNum>();
static cHpNum HpINF = Inf<HpNum>();

// an Index with the most significant bit only == 1
static cIndex LBIndex = 1 << ( numeric_limits<Index>::digits - 1 );
// an Index with the second-most significant bit only == 1
static cIndex LLBIndex = LBIndex >> 1;
// masking the two most significant bits in an Index
static cIndex IMask = LBIndex + LLBIndex;

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::OSIMPSolver( istream *iStrm )
    :
    MPSolver( )
{
 // reset the pointers to OsiSolver and FiOracle object  - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr = 0;
 FIO = 0;

 // set the default algorithm parameters - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 useactiveset = false;  // no active set
 first = true;          // no problem has been solved so far
 checkID = false;       // by default don't check for identical items

 stab = unset;	        // the stabilization term is unknown
 t = 1;			// the proximity parameter is set to a standard value

 // the master is empty- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MaxBSize = 0;		// max bundle is zero
 NSubG = NConst = 0;    // bundle is empty

 dict_item = 0;	        // no items are in the master problem
 dict_item_maxname = 0;
 dict_slack = 0;
 dict_stab = 0;

 wcomp = 0;
 weasy = 0;

 NewItem = 0;	        // new item information is empty
 NewItemBse = 0;
 NewItemBDm = 0;

 NNVars = BxdVars = 0;
 Upper = Lower = 0;
 FiLambda = FiLambda1 = 0;

 tempHP = 0;             // temporary buffers of HpNum
 tempI = 0;              // and Index
 tempHP_size = 0;        // initially they are empty
 tempI_size = 0;

 Aset = 0;                // no active sets
 Asetdim = 0;

 // no solutions is given  - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 comp_row = 0;
 comp_col = 0;

 RhoCol = 0;
 RhoColBse = 0;
 RhoColBDm = 0;

 FsbEps = Eps<HpNum>();
 HasLwrBnd = false;

 GiPerd = 0;

 #if( PRESERVE_OSI_SOLS )
  csol = rsol = rcst = 0;
  csols = rsols = 0;
 #endif

 derhand = new DerivedHandler;
 }  // end( OSIMPSolver::OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetDim( cIndex MxBSz , FiOracle *Oracle ,
			  const bool UsAvSt )
{
 MSG( 0 , "OSIMPSolver::SetDim()\n" );

 if( ! MxBSz ) {  // deallocate all its memory and quietly wait
  cleanup();      // for new instructions- - - - - - - - - - - - - - - - - - -
  return;
  }

 if( Oracle ) {  // allocate the memory for solving the problem- - - - - - - -
  if( ! osiSlvr )
   throw( NDOException( "OSIMPSolver::SetDim: osiSlvr must be set" ) );

  // discards all the previous settings and deallocate all the memory- - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cleanup();

  // set the problem as a maximum problem  - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->setObjSense( -1.0 );

  // ask FiOracle the instance dimension - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  MaxBSize = MxBSz;
  FIO = Oracle;
  MaxSGLen = Oracle->GetMaxNumVar();
  CrrSGLen = Oracle->GetNumVar();
  NrFi = Oracle->GetNrFi();
  MaxNZ = Oracle->GetMaxNZ();
  useactiveset = UsAvSt;

  // allocate the memory - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  weasy = new bool[ NrFi + 1 ];
  weasy[ 0 ] = true;  // the 0-th is very easy

  NewItem = new SgNum[ MaxNZ ];

  NSubG = new Index[ NrFi + 1 ];
  VectAssign( NSubG , Index( 0 ) , NrFi + 1 );
  NConst = new Index[ NrFi + 1 ];
  VectAssign( NConst , Index( 0 ) , NrFi + 1 );

  comp_row = new Index[ NrFi + 1 ];
  comp_col = new Index[ NrFi + 1 ];
  dict_item = new Index[ MaxBSize ];
  VectAssign( dict_item , InINF , MaxBSize );
  dict_slack = new Index[ MaxSGLen ];
  VectAssign( dict_slack , InINF , MaxSGLen );
  wcomp = new Index[ MaxBSize ];
  VectAssign( wcomp , InINF , MaxBSize );

  Upper = new LMNum[ MaxSGLen ];
  Lower = new LMNum[ MaxSGLen ];

  dict_item_maxname = 0;

  // create all the rows for static constraints and (partial) set of rows
  // for the dynamic constraints. The first constraint is relative to rho.
  // Note that variable r has been deleted since it's redundant. So the first
  // constraint is 1 - rho >= 0. All in all in the absence of total lower
  // bound case this constraint should be rho = 1 (the default case).
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->addRow( 0 , 0 , 0 , 1.0 , 1.0 );

  // for each easy component generate (2 * (BNR + BNC)) structured constraints
  // while for the difficult ones create just the simplex constraint
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index BNC , BNR , BNZ , ANZ;
  Index count, start;

  for( Index i = 0 ; i++ < NrFi ; ) {
   comp_row[ i - 1 ] = osiSlvr->getNumRows();
   BNC = FIO->GetBNC( i );
   if( BNC ) {
    weasy[ i ] = true;

    BNR = FIO->GetBNR( i );  // number of rows of matrix B[ i ]

    double *lhs;
    double *rhs;
    double *lbd;
    double *ubd;

    lbd = new double[ BNC ];
    ubd = new double[ BNC ];
    lhs = new double[ BNR ];
    rhs = new double[ BNR ];

    if( BNR )
     FIO->GetBDesc( i , 0 , 0 , 0 , lhs , rhs , 0 , lbd , ubd );
    else
     FIO->GetBDesc( i , 0 , 0 , 0 , 0 , 0 , 0 , lbd , ubd );

    // allocate the rows having a finite bound  - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    for( Index j = 0 ; j < BNC ; j++ ) // if lbd is well defined
     if( lbd[ j ] > -Inf<double>() )
      osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0.0 );

    for( Index j = 0 ; j < BNC ; j++ ) // if ubd is well defined
     if( ubd[ j ] < Inf<double>() )
      osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0.0 );

    if( BNR ) { // check whether or not some rows exist
     for( Index j = 0 ; j < BNR ; j++ ) // if lhs is well defined
      if( lhs[ j ] > -Inf<double>() )
       osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0.0 );

     for( Index j = 0 ; j < BNR ; j++ ) // if rhs is well defined
      if( rhs[ j ] < Inf<double>() )
       osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0.0 );
     }

    delete[] lbd;
    delete[] ubd;
    delete[] lhs;
    delete[] rhs;

    MSG( 0 , "Structured constraints for the easy component " << i
         << " have been created starting from the row "  << comp_row[ i - 1 ]
	 << endl );
    }
   else {
    weasy[ i ] = false;
    osiSlvr->addRow( 0 , 0 , 0 , 0.0 , 0.0 );

    MSG( 0 , "Structured constraints for the difficult component  " << i
	  << " have been created into the row " << comp_row[ i - 1 ] << endl
	 );
    }
   }

  // if the initial dimension of the items is greater than zero, generate the
  // (partial) rows of the the dynamic part which describes the form of the
  // dual variables z  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_row[ NrFi ] = osiSlvr->getNumRows();
  if( useactiveset ) { // all the constraints are inactive so far  - - - - - -
    MSG( 0 , "The dynamic inactive constraints have been created starting "
         "from the row " << comp_row[ NrFi ] << endl );

    for( Index j = CrrSGLen ; j > 0 ; j-- )
     osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() ,
		                    osiSlvr->getInfinity() );
   }
  else {     // all the constraints are active (now and for ever)  - - - - - -
    MSG( 0 , "The dynamic constraints have been created starting from the row "
	 << comp_row[ NrFi ] << endl );
    for( Index j = CrrSGLen ; j > 0 ; j-- )
     osiSlvr->addRow( 0 , 0 , 0 , 0.0 , 0.0 );
   }

  // allocate and initialize dictionary and values of the column
  // relative to the special variable rho. Again, if a certain total lower
  // bound exists, the coefficient of rho in the objective function is *-*
  // total lower bound. Else this coefficient is equal to 0 and rho is to set
  // to 1. Note also that rho is the multiplier of the subgradient of 0-th
  // component, that is, the b vector; this vector will be provided later.

  RhoCol = new double[ comp_row[ NrFi ] + MaxSGLen ];
  RhoColBse = new int[ comp_row[ NrFi ] + MaxSGLen ];

  RhoColBse[ 0 ] = 0;
  RhoCol[ 0 ] = 1.0;
  RhoColBDm = 1;

  // allocate the memory for the temporary vectors tempI and tempHp
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  resizeI( osiSlvr->getNumRows() );
  resizeHP( osiSlvr->getNumRows() );

  // allocate the memory for the objective function values
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // however, note that if the easy components are defined, then the master
  // problem provides the Lagrangian object values of them

  FiLambda = new HpNum[ NrFi ];
  FiLambda1 = new HpNum[ NrFi ];
  VectAssign( FiLambda , HpINF , NrFi );  // Fi( Lambda ) == ?

  // allocate the memory for Gi * d
  GiPerd = new HpNum[ MaxBSize ];

  // recover the static information of the easy component  - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for ( Index i = 1 ; i <= NrFi ; i++ ) {
   BNC = FIO->GetBNC( i );
   if( BNC ) { // FiOracle returns all the information
	    // in a sparse format, so the columns of the matrices A and B are
	    // in sparse format- - - - - - - - - - - - - - - - - - - - - - - -

    BNR = FIO->GetBNR( i );  // number of rows of matrix B[ i ]
    BNZ = FIO->GetBNZ( i );  // number of non-zeroes of matrix B[ i ]
    ANZ = FIO->GetANZ( i );  // number of non-zeroes of matrix A[ i ]

    // allocate the memory to describe the matrix B[ i ]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    double *cst = new double[ BNC ];

    // allocate the memory for the matrix B[ i ] excluding the infinity bounds
    // relative to variables/constraints - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    int* Bbeg;
    int* Bind;
    double *Bval;

    double *lhs;
    double *rhs;

    double *lbd = new double[ BNC ];
    double *ubd = new double[ BNC ];

    if( BNR ) { // this is equivalent to the condition: if( BndForm & ( ~3 ) )
     Bbeg = new int[ BNC + 1 ];
     Bind = new int[ BNZ ];
     Bval = new double[ BNZ ];
     lhs = new double[ BNR ];
     rhs = new double[ BNR ];
     }
    else {
     Bbeg = Bind = 0;
     Bval = 0;
     lhs = rhs = 0;
     }

    // allocate the memory to describe the matrix A[ i ]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    int *Abeg = new int[ BNC + 1 ];
    int *Aind = new int[ ANZ ];
    double *Aval = new double[ ANZ ];

    // ask FiOracle for the matrices A[ i ] and B[ i ] (and not only them...)
    // Note that the FiOracle may give a partial information of the rows
    // of matrix A[ i ] and precisely the first CrrSGLen rows.
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    FIO->GetBDesc( i , Bbeg , Bind , Bval , lhs , rhs , cst , lbd , ubd );
    FIO->GetADesc( i , Abeg , Aind , Aval );

    // mark the start of the variables x[ i ] of the i-th easy component
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    comp_col[ i ] = osiSlvr->getNumCols();

    // assign a name to the constraints, taking into account only those whose
    // the corresponding bounds are finite - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Index_Set NmRowsRhs = 0;
    Index_Set NmRowsLhs = 0;

    start = 0;
    if( BNR ) { // if some rows exist, give a name to each row excluding
    	   // all of them which have an infinite bound

     count = 0;
     NmRowsRhs = new Index[ BNR ];
     for( Index j = 0 ; j < BNR ; j++ )   // if rhs is well defined,
      if( rhs[ j ] < Inf<double>() )      // assign a name
       NmRowsRhs[ j ] = j - count;
      else {                              // count the infinite values
       NmRowsRhs[ j ] = Inf<Index>();
       count++;
       }

     start = BNR - count;  // start from here ...
     count = 0;
     NmRowsLhs = new Index[ BNR ];
     for( Index j = 0 ; j < BNR ; j++ )   // if lhs is well defined,
      if( lhs[ j ] > - Inf<double>() )      // assign a name
       NmRowsLhs[ j ] = j - count + start;
      else {
       NmRowsLhs[ j ] = Inf<Index>();
       count++;
       }

     start += ( BNR - count );
     } // end assignment - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    count = 0;
    Index_Set NmRowsUbd = new Index[ BNC ];
    for( Index j = 0 ; j < BNC ; j++ )   // if ubd is well defined,
     if( ubd[ j ] < Inf<double>() )      // assign a name
      NmRowsUbd[ j ] = j - count + start;
     else {
      NmRowsUbd[ j ] = Inf<Index>();
      count++;
      }

    start += ( BNC - count );
    count = 0;
    Index_Set NmRowsLbd = new Index[ BNC ];
    for( Index j = 0 ; j < BNC ; j++ )   // if lbd is well defined,
     if( lbd[ j ] > - Inf<double>() )    // assign a name
      NmRowsLbd[ j ] = j - count + start;
     else {
      NmRowsLbd[ j ] = Inf<Index>();
      count++;
      }

    // write the rows relative to the easy components  - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    for ( Index j = 0 ; j < BNC ; j++ ) {  // for each variable j of the
     count = 0;                            // i-th easy component do ...

     // write the static part: B's columns and the box constraints
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     if( BNR )
      for( int k = Bbeg[ j ] ; k < Bbeg[ j + 1 ] ; k++ ) {
       if( rhs[ Bind[ k ] ] < Inf<double>() ) {
        tempI[ count ] = NmRowsRhs[ Bind[ k ] ] + comp_row[ i - 1 ];
        tempHP[ count++ ] = Bval[ k ];   // B[ i ] x[ i ] - e[ i ] <= 0
        }
       if( lhs[ Bind[ k ] ] > - Inf<double>() ) {
        tempI[ count ] = NmRowsLhs[ Bind[ k ] ] + comp_row[ i - 1 ];
        tempHP[ count++ ] = -Bval[ k ];   // -B[ i ] x[ i ] + d[ i ] <= 0
        }
      }

     if( ubd[ j ] < Inf<double>() ) {
      tempI[ count ] = comp_row[ i - 1 ] + NmRowsUbd[ j ];
      tempHP[ count++ ] = 1.0;          // x[ i ] - u[ i ] <= 0
      }

     if( lbd[ j ] > -Inf<double>() ) {
      tempI[ count ] = comp_row[ i - 1 ] + NmRowsLbd[ j ];
      tempHP[ count++ ] = -1.0;         // -x[ i ] + l[ i ] <= 0
      }

     // write the dynamic part: A's columns  - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     for( int k = Abeg[ j ] ; k < Abeg[ j + 1 ] ; k++ ) {
      tempI[ count ] = Aind[ k ] + comp_row[ NrFi ];
      tempHP[ count++ ] = Aval[ k ];
      }

     // add the j-th column of the component i
     // Note that objective function of the easy component has the form
     // as c[ i ] - Lambda * A[ i ] ) x[ i ] - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     osiSlvr->addCol( count , tempI , tempHP , - osiSlvr->getInfinity() ,
		      osiSlvr->getInfinity() , cst[ j ] );

     MSG( 0 , "The variables "<< j << "of the easy component "<< i
          << "has been added" << endl );

     }  // end adding of j-th column of the component i- - - - - - - - - - - -
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // put in RhoCol the vectors e[ i ], d[ i ], u[ i ] and l[ i ] [see above]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( BNR )
     for( Index j = 0 ; j < BNR ; j++ ) {
      if( rhs[ j ] < Inf<double>() ) {
       RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + NmRowsRhs[ j ];
       RhoCol[ RhoColBDm++ ] = - rhs[ j ];
       }
      if( lhs[ j ] > - Inf<double>() ) {
       RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + NmRowsLhs[ j ];
       RhoCol[ RhoColBDm++ ] = lhs[ j ];
       }
      }

    for( Index j = 0 ; j < BNC ; j++ ) {
     if( ubd[ j ] < Inf<double>() ) {
      RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + NmRowsUbd[ j ];
      RhoCol[ RhoColBDm++ ] = - ubd[ j ];
      }
     if( lbd[ j ] > - Inf<double>()) {
      RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + NmRowsLbd[ j ];
      RhoCol[ RhoColBDm++ ] = lbd[ j ];
      }
     }

    // deallocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    delete[] cst;
    delete[] lbd;
    delete[] ubd;
    delete[] Bbeg;
    delete[] Bind;
    delete[] Bval;
    delete[] lhs;
    delete[] rhs;
    delete[] Abeg;
    delete[] Aind;
    delete[] Aval;
    delete[] NmRowsLhs;
    delete[] NmRowsRhs;
    delete[] NmRowsLbd;
    delete[] NmRowsUbd;

    }
   }  // end easy components part description- - - - - - - - - - - - - - - - -
  //- -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // add the description of the difficult components - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 1 ; i <= NrFi ; i++ ) {
   BNC = FIO->GetBNC( i );
   if( BNC == 0 ) { // for each component "i" we must create a
    // variable \f$gamma^i\f$. The coefficient of the \f$gamma^i\f$ into
    // objective is - lower bound of the i-th component. The coefficient
	// is 0, if no lower bound is given as well as rho case but unlike the
    // rho case, \f$gamma^i\f$ must to be set to 0 whenever some item exist.
    // In addition, it provides the simplex set
    // \f$\sum_{k} \theta_k^i + \gamma^i = \rho \f$
    // Note that at beginning no theta multipliers are known.

    comp_col[ i ] = osiSlvr->getNumCols();
    tempI[ 0 ] = comp_row[ i - 1 ];
    tempHP[ 0 ] = 1.0;
    osiSlvr->addCol( 1 , tempI , tempHP , 0.0 , osiSlvr->getInfinity() ,
		     0.0 );
    MSG( 0 , "The variable gamma_" << i << " has been added and its name is "
         << comp_col[ i ] << endl );

    RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ];
    RhoCol[ RhoColBDm ] = - 1.0;
    RhoColBDm++;
    }
   }  // end difficult components part description - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // define the slack and stabilization variables- - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( stab == quadratic ) {
   switchToQP();
   dict_stab = new Index[ MaxSGLen ];
   VectAssign( dict_stab , InINF , MaxSGLen );
   }
  else
   dict_stab = 0;

  for( Index i = 0 ; i < CrrSGLen ; i++) { // ask FiOracle if there are some
                               // lower and(/or) upper bounds on the primal
   bool slack_p = false;       // variables  - - - - - - - - - - - - - - - - -
   bool slack_m = false;

   if( ( Upper[ i ] = Oracle->GetUB( i ) ) < LMINF )
    slack_m = true;  // there is some upper bound on variable i
   if( Oracle->GetUC( i ) )
    Lower[ i ] = -LMINF;
   else {            // the primal variable is constrained to be nonnegative
    Lower[ i ] = 0;
    slack_p = true;
    }

   if( slack_p )
    NNVars++;

   if( slack_p || slack_m )
    BxdVars++;

   // make decision on stabilization - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   switch( stab ) {
    case none:
    case boxstep:
     slack_p = slack_m = true;
     break;
    case quadratic:
     tempI[ 0 ] = comp_row[ NrFi ] + i;
     tempHP[ 0 ] = 1.0;
     dict_stab[ i ] = osiSlvr->getNumCols();
     osiSlvr->addCol( 1 , tempI , tempHP , - osiSlvr->getInfinity() ,
		      osiSlvr->getInfinity() , 0.0 );
     MSG( 0 , "The stabilization variable z_"<< i
	  <<" has been created and its name is " << dict_stab[ i ] << endl );
     break;
    default:
     throw( NDOException( "OSIMPSolver::SetDim: undecided stabilization" ) );
    }

   // define the slack variables and put in the vocabulary, one of each pair
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( slack_p ) {
    tempI[ 0 ] = comp_row[ NrFi ] + i;
    tempHP[ 0 ] = 1.0;
    osiSlvr->addCol( 1 , tempI , tempHP , 0.0 , osiSlvr->getInfinity() ,
		     0.0 );
    dict_slack[ i ] = osiSlvr->getNumCols () - 1;
    MSG( 0 , "The slack s_"<< i <<"+ has been created and its name is"
    	<< dict_slack[ i ] << endl );
    }

   if( slack_m ) {
    tempI[ 0 ] = comp_row[ NrFi ] + i;
    tempHP[ 0 ] = -1.0;
    osiSlvr->addCol( 1 , tempI , tempHP , 0.0 , osiSlvr->getInfinity() ,
		     0.0 );
    if( ! slack_p )
     dict_slack[ i ] = osiSlvr->getNumCols() - 1;
     MSG( 0 , "The slack s_" << i << "- has been created and its name is "<<
    	( slack_p? dict_slack[ i ] + 1 : dict_slack[ i ] ) << endl );
     }

   }  // end stabilization and slack definition- - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // add the column of rho in the master - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_col[ 0 ] = osiSlvr->getNumCols();
  osiSlvr->addCol( RhoColBDm , RhoColBse , RhoCol , 0.0 ,
		   osiSlvr->getInfinity() , 0.0 );
  MSG( 0 , "rho has been created and its name is " << comp_col[ 0 ] << endl );

  // change the price of the slack variables - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tUpdatePrices();       // change the prices for the stabilization
  if( stab != boxstep )  // also change the slack prices (if they have not
   ptUpdatePrices();     // already been changed in tUpdatePrices())

  }  // end( if(Oracle) )
 else {  // sets the max bundle size to n and activate/deactivate the Active
         // Set Mechanism without changing anything else - - - - - - - - - - -

  // deletes the item in excess  - - - - - - - - - - - - - - - - - - - - - - -

  for( ; dict_item_maxname-- > MxBSz ; )
   if( dict_item[ dict_item_maxname ] != InINF )
    RmvItem( dict_item_maxname );

  if( MxBSz != MaxBSize ) {  // the existing items in the bundle (if any) are
   // all kept if MxBSz is larger than MaxBSize, but a smaller value will
   // force deletion of all the items with "name"  >= MxBSz

   Index_Set olddict_item = dict_item;
   Index_Set oldwcomp = wcomp;

   dict_item = new Index[ MxBSz ];
   wcomp = new Index[ MxBSz ];

   VectAssign( dict_item , InINF , MxBSz );
   VectAssign( wcomp , InINF , MxBSz );

   cIndex copydim = ( MxBSz < MaxBSize ) ? MxBSz : MaxBSize;

   if( olddict_item ) {
    VectAssign( dict_item , olddict_item , copydim );
    delete[] olddict_item;
    }

   if( oldwcomp ) {
    VectAssign( wcomp , oldwcomp , copydim );
    delete[] oldwcomp;
    }

   MaxBSize = MxBSz;  // new dimension of the bundle

   // the initial "active set" of variables is *empty* and so all the
   // variables must be non active - - - - - - - - - - - - - - - - - - - - - -

   if( ! useactiveset && UsAvSt )
    for( Index i = 0 ; i < CrrSGLen ; i++ )
     deactivate( i );

   // all the  variables are considered to be always "active"
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( useactiveset && ! UsAvSt )
    for( Index i = 0 ; i < CrrSGLen ; i++ )
     activate( i );

   useactiveset = UsAvSt;
   }
  }  // end( else(Oracle) )- - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 }  // end( OSIMPSolver::SetDim )

/*--------------------------------------------------------------------------*/

inline void OSIMPSolver::Sett( cHpNum tt )
{
 if( t != tt ) {  // actually update t- - -  - - - - - - - - - - - - - - - - -
  t = tt;
  MSG( 0 , "OSIMPSolver::Sett(): t = " << t << endl );
  tUpdatePrices();
  }
 }  // end( OSIMPSolver::Sett )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetOptPrcsn( HpNum OEps )
{
 MSG( 0 , "OSIMPSolver::SetOptPrcsn()\n" );
 osiSlvr->setDblParam( OsiDualTolerance , double( OEps ) );

 #if CPLEX
  if( algorithm == kBar ) {
   OsiCpxSolverInterface *osiCpx = 
                              dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
   if( ! osiCpx )
    throw( NDOException( 
                 "OSIMPSolver::SetOptPrcsn: the OSI solver is not Cplex" ) );

   CPXENVptr env = osiCpx->getEnvironmentPtr();
   CPXsetdblparam( env , CPX_PARAM_BAREPCOMP , OEps );
  }
 #endif

 }  // end( OSIMPSolver::SetOptPrcsn )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetFsbPrcsn( HpNum FEps )
{
 MSG( 0 , "OSIMPSolver::SetFsbPrcsn()\n" );
 FsbEps = FEps;
 osiSlvr->setDblParam( OsiPrimalTolerance , double( FEps ) );

 }  // end( OSIMPSolver::SSetFsbPrcsn )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetLowerBound( cHpNum LwBnd , cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::SetLowerBound()\n" );

 if( ( wFi < InINF ) && weasy[ wFi ] )
  throw( NDOException(
    "OSIMPSolver::SetLowerBound: Setting LowerBound for a easy component" ) );

 if( LwBnd > -HpINF ) {
  if( wFi < InINF ) {  // insert the constraint gamma_i <= 1 and
	               // add LwBnd * gamma_i to the objective function
   osiSlvr->setColUpper( comp_col[ wFi ] , 1.0 );
   osiSlvr->setObjCoeff( comp_col[ wFi ] , LwBnd );
   MSG( 0 , "LB( "<< wFi << " ) = " << LwBnd << endl );
   }
  else {   // set the total lower bound- - - - - - - - - - - - - - - - - - - -
   // Note that it needs change the constraint of rho:
   // r + rho = 1 ----> r = 1 - rho => 0 [ rho <= 1 ]
   // and in the objective function LwBnd * r = LwBnd - LwBnd * rho
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   HasLwrBnd = true;
   osiSlvr->setRowBounds( 0 , - osiSlvr->getInfinity() , 1.0 );
   osiSlvr->setObjCoeff( comp_col[ 0 ] , - LwBnd );
   MSG( 0 , "Total LB = " << LwBnd << endl );
   }
  }
 else {  // if LwBnd = -HpINF reset the previous lower bound ...
  if( wFi < InINF ) {  // if some items exist set gamma_i = 0
   MSG( 0 , "Lower bound of "<< wFi << " has been deleted " << endl );

   for( Index i = 0 ; i < dict_item_maxname ; i++ )
    if( WComponent( i ) == wFi ) {
     osiSlvr->setObjCoeff( comp_col[ wFi ] , 0.0 );
     if( dict_item_maxname )
      osiSlvr->setColUpper( comp_col[ wFi ] , 0.0 );
     break;
     }
   }
  else {  // change the rho constraint --> rho = 1 and put the coefficient
	  // of rho equal to zero in the objective function

   HasLwrBnd = false;
   MSG( 0 , "Total lower bound is - Inf<double>()\n" );
   osiSlvr->setObjCoeff( comp_col[ 0 ] , 0.0 );
   osiSlvr->setRowBounds( 0 ,  1.0 , 1.0 );
   }
  }
 }  // end( OSIMPSolver::SetLowerBound )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetMPLog( ostream *outs , const char lvl )
{
 MPSolver::SetMPLog( outs , lvl );

 #if( OSIMPSOLVERLOG )

  /* if( osiSlvr )
    osiSlvr->messageHandler()->setLogLevel( int( MPLLvl ) ); */

  #if( CPLEX )
   OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>(
								 osiSlvr );
   if( ! osiCpx )
    throw( NDOException( "OSIMPSolver::SetMPLog: the OSI solver is not Cplex" ) );

   CPXENVptr env = osiCpx->getEnvironmentPtr();
   CPXsetintparam( env , CPX_PARAM_MIPDISPLAY , MPLLvl );

   logfile = CPXfopen( "cplex.log" , "w" );
   if( CPXsetlogfile( env , logfile ) )
    throw( NDOException( "OSIMPSolver::SetMPLog: error " ) );
  #endif
 #endif

 } // end( OSIMPSolver::SetMPLog )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetAlgo( const OsiAlg algo, const OsiRed reduc )
{
 algorithm = algo;

 #if CPLEX
 OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>(
								  osiSlvr );
 if( ! osiCpx )
  throw( NDOException( "OSIMPSolver::SetAlgo: the OSI solver is not Cplex" )
	 );

 CPXENVptr env = osiCpx->getEnvironmentPtr();

 // primal and dual reduction  - - - - - - - - - - - - - - - - - - - - - -
 switch ( reduc ) {
  case ( rNo ):
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_NOPRIMALORDUAL );
   break;
  case ( rPrim ):
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALONLY );
   break;
  case ( rDual ):
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_DUALONLY );
   break;
  case ( rBoth ):
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALANDDUAL );
   break;
  default:
   throw( NDOException(
                "OSIMPSolver::SetAlgo: reduction hasn't been selected" ) );
  }

 if( stab != quadratic )  // choose CPLEX algorithm   - - - - - - - -
  switch( algorithm ) {
   case( kAuto ): CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_AUTOMATIC );;
                  break;
   case( kPrim ): CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_PRIMAL );
                  break;
   case( kDual ): CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_DUAL );
                  break;
   case( kNet ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_NET );
                  break;
   case( kBar ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_BARRIER );
                  break;
   case( kSif ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_SIFTING );
                  break;
   case( kCon ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_CONCURRENT );
                  break;
   default:
    throw( NDOException(
            "OSIMPSolver::SetAlgo: no algorithms has been selected" ) );

   } // end switch
 else
  switch( algorithm ) {
   case( kAuto ): CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_AUTOMATIC );
                  break;
   case( kPrim ): CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_PRIMAL );
                  break;
   case( kDual ): CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_DUAL );
                  break;
   case( kNet ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_NET );
                  break;
   case( kBar ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_BARRIER );
                  break;
   case( kSif ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_SIFTING );
                  break;
   case( kCon ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_CONCURRENT );
                  break;
   default:
     throw( NDOException( "OSIMPSolver::SetAlgo: wrong algorithm" ) );

   } // end switch
 #endif

 } // end ( OSIMPSolver::SetAlgo )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetOsi( OsiSolverInterface *osi )
{
 if( osiSlvr ) {
  if( FIO )
   throw( NDOException(
          "OSIMPSolver::SetosiSlvr: cannot change OSI solver in flight" ) );
  else
   delete osiSlvr;  // the MP is empty, just switching solver
  }

 // to dealt with the quadratic master problem it needs to use the specialized
 // solver because of the weakness of the interface of OsiSolver. In fact, it
 // doesn't support the quadratic case - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( osiSlvr = osi ) )
  throw( NDOException( "OSIMPSolver::SetOsi: wrong OSI solver" ) );

 // set to 0 the log: do that to avoid some warnings in cleanup
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // osiSlvr->messageHandler()->setLogLevel( int( MPLLvl ) );
 osiSlvr->passInMessageHandler( derhand );

 } // end( OSIMPSolver::SetOsi )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetStabType( const StabFun sf )
{
 MSG( 0 , "OSIMPSolver::SetStabType(()\n" );
 if( ( ( stab = sf ) != unset ) && FIO )
  throw( NDOException( 
    "OSIMPSolver::SetStabType: cannot change D_t in flight" ) );

 switch( stab ) {
  case none:
  case boxstep:
  case quadratic:
	   break;
  case unset:
  default: throw( NDOException(
	"OSIMPSolver::SetDim: stabilization function not specified yet" ) );
  }
 }  // end( OSIMPSolver::SetStabType )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::MPStatus OSIMPSolver::SolveMP( void )
{
 MSG( 0 , "OSIMPSolver::SolveMP()\n" );

 if( MPt )
  MPt->Start();

 // get the dimension of the Master Problem  - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int nc = osiSlvr->getNumCols();
 int nr = osiSlvr->getNumRows();

 #if( OSIMPSOLVERLOG )
  if( osiSlvr && MPLLvl > 1 && nc && nr ) {
   std::strstream filename;
   filename << "MP_" << FIO->GetNDOSolver()->NrIter() << std::ends;
   osiSlvr->writeMps( filename.str() );
   }
 #endif

 // solve the Master Problem - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( stab == quadratic ) {
  #if CPLEX == 1
   OsiCpxSolverInterface *osiCpx =
                            dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
   if( ! osiCpx )
    throw( NDOException(
                    "OSIMPSolver::SolveMP: the OSI solver is not Cplex" ) );

   CPXENVptr env = osiCpx->getEnvironmentPtr();
   CPXLPptr qp = osiCpx->getLpPtr();

   CPXqpopt( env , qp );  // call the QP solver- - - - - - - - - - - - - - -
  #else
   throw( NDOException( "OSIMPSolver::SolveMP: not implemented yet" ) );
  #endif
  }
 else
  if( first ) {
   osiSlvr->initialSolve();
   first = false;
   }
  else
   osiSlvr->resolve();

 // check the status of the Solver and take the necessary action - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( osiSlvr->isProvenDualInfeasible() ) {  // the feasible set is empty
  MSG( 0 , "MP is dual infeasible" << endl );
  return( kUnfsbl );
  }

 if( osiSlvr->isProvenPrimalInfeasible() ) {
  MSG( 0 , "MP is primal infeasible" << endl );
  return( kUnbndd );
  }

 if( ! osiSlvr->isProvenOptimal() )
  if( osiSlvr->isAbandoned() )
   MSG( 0 , "Warning: numerical difficulties in the solver, ignoring them"
	<< endl );
  else {
   MSG( 0 , "Some unknown error happened" << endl );
   return( kError );
   }

 // record solution information in its data structures - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( PRESERVE_OSI_SOLS )
  // note: for the quadratic MP one can still call getColSolution(),
  // getReducedCost() and getRowPrice(), as Cplex uses the same functions
  // as in the LP case (and precisely CPXgetx, CPXgetdj and CPXgetpi)

  if( nc > csols ) {
   delete[] csol;
   delete[] rcst;
   csols = nc;
   csol = new double[ csols ];
   rcst = new double[ csols ];
   }

  VectAssign( csol , osiSlvr->getColSolution() , nc );
  VectAssign( rcst , osiSlvr->getReducedCost() , nc );

  if( nr > rsols ) {
   delete[] rsol;
   rsols = nr;
   rsol = new double[ rsols ];
   }

  VectAssign( rsol , osiSlvr->getRowPrice() , nr );
 #else
  const double *rsol = osiSlvr->getRowPrice();
  const double *rcst = osiSlvr->getReducedCost();
 #endif

 // to compute Gid use the reduced cost rc_i, in particular we have:
 // g^{top} * delta = rc_i + v_i + \alpha_i

 const double *objcoeff = osiSlvr->getObjCoefficients();
 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( dict_item[ i ] < InINF ) {
   GiPerd[ i ] = rcst[ dict_item[ i ] ] - objcoeff[ dict_item[ i ] ];
   if( IsSubG( i ) )
    GiPerd[ i ] += rsol[ comp_row[ WComp( i ) - 1 ] ];
   }

 // compute the actual function value of the easy components by using the
 // row price: Fi[ wFi ]( Lambda + d* ) =
 // \pi_e[ wFi ] * e[ wFi ] - \pi_d[ wFi ] * d[ wFi ] +
 // \pi_u[ wFi ] * u[ wFi ] - pi_l[ wFi ] * l[ wFi ]
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( FiLambda1 , HpNum( 0 ) , NrFi );
 for( Index wFi = 1 ; wFi <= NrFi; wFi++ )
  if( ( wFi <= NrFi ) && weasy[ wFi ] ) {
   for( Index i = 0 ; i < RhoColBDm ;  i++)
    if( ( Index( RhoColBse[ i ] ) >= comp_row[ wFi - 1 ] ) &&
	( Index( RhoColBse[ i ] ) < comp_row[ wFi ] ) )
     FiLambda1[ wFi - 1 ] -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];

   MSG( 0 , "FiBLambda(" << wFi << ") = " << FiLambda1[ wFi - 1 ] << endl );
   }

 // print the primal/dual solution  - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( OSIMPSOLVERLOG )
  if( MPLLvl > 2 ) {
   #if( PRESERVE_OSI_SOLS == 0 )
    const double *csol = osiSlvr->getColSolution();
   #endif

   *MPLog << "x =" << endl;
   for( int i = 0 ; i < nc ; i++ )
    *MPLog << csol[ i ]  << " ~ ";
   *MPLog << endl;

   *MPLog << "Pi =" << endl;
   for( int i = 0 ; i < nr ; i++ )
    *MPLog << rsol[ i ] << " ~ ";
   *MPLog << endl;
   }
 #endif

 // mark all items as read  - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( wcomp[ i ] < InINF )
   wcomp[ i ] &= ~LLBIndex;

 if( MPt )
  MPt->Stop();

 return( kOK );

 }  // end( OSIMPSolver::SolveMP )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR READING RESULTS ------------------------*/
/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadFiBLambda( cIndex wFi )
{
 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif

 HpNum FiVal = 0;
 if( ! wFi ) {  // the 0-th component- - - - - - - - - - - - - - - - - - - -
  // it is given by Fi[ 0 ]( d^* ) = d^* b. Note that d^* is the dual
  // solution of the dynamic part and it represents the search direction

  for( Index i = comp_row[ NrFi ] ; i < RhoColBDm  ; i++ )
   FiVal -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];

  MSG( 0 , "FiBLambda( 0 ) = " << FiVal << endl );
  }
 else
  if( wFi <= NrFi ) {  // a specific component - - - - - - - - - - - - - - -
   // the easy case is quite different from the difficult one: in fact, in
   // the easy case it computes the actual value of the function while in
   // the other case just an approximation is provided
   if( weasy[ wFi ] )
    FiVal = FiLambda1[ wFi - 1 ];
   else
    FiVal = rsol[ comp_row[ wFi - 1 ] ];

   MSG( 0 , "FiBLambda(" << wFi << ") = " << FiVal << endl );
   }
  else { // all the components - - - - - - - - - - - - - - - - - - - - - - -

   for( Index i = 1 ; i <= NrFi ; i++ )
    if( weasy[ i ] )
     FiVal += FiLambda1[ i - 1 ];
    else
     FiVal += rsol[ comp_row[ i - 1 ] ];

   if( wFi == InINF )   // ... comprised the 0-th
    for( Index i = comp_row[ NrFi ] ; i < RhoColBDm ; i++ )
     FiVal -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];

   MSG( 0 , "FiBLambda( INF ) = " << FiVal << endl );

   }  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( FiVal );

 }  // end( OSIMPSolver::ReadFiBLambda )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadDt( cHpNum tt )
{
 HpNum value = 0;
 HpNum mytt;

 #if( PRESERVE_OSI_SOLS == 1 )
  cLMRow primal = rsol + comp_row[ NrFi ];
 #else
  cLMRow primal = osiSlvr->getRowPrice() + comp_row[ NrFi ];
  const double *csol = osiSlvr->getColSolution();
 #endif

 switch( stab ) {
  case none: break;
  case boxstep:
  mytt = tt * ( 1 + FsbEps );  // enlarge tt a little to take in account
                               // for numerical errors (use the tolerance)
  if( mytt < t ) {  // check whether d lies in the ball of radious tt
   for( Index n = 0 ; n < CrrSGLen ; n++ )
    if( ABS( primal[ n ] ) > mytt ) {
     MSG( 0 , "Dt = +INF" << endl );
     return( HpINF );
     }
    }
   break;
  case quadratic:
   for( Index n = 0 ; n < CrrSGLen ; n++ )
	value += csol[ dict_stab[ n ] ] * csol[ dict_stab[ n ] ];
   value *= ( ( tt == t ? tt : ( t * t ) / tt ) / 2 );
   MSG( 0 , "Dt = " << value << endl );
   break;
  default:
   throw( NDOException( "OSIMPSolver::ReadDt: undecided stabilization" ) );
   }

 MSG( 0 , "Dt = " << value << endl );
 return( value );

 }  // end( OSIMPSolver::ReadDt )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadSigma( cIndex wFi )
{
 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif
 const double *obj = osiSlvr->getObjCoefficients();

 // compute sigma  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum value = 0;
 if( wFi > NrFi ) {  // aggregate the linearization error and the rhs of the
                     // constraints of the difficult components

  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( ( dict_item[ i ] < InINF ) && WasInMP( i ) )
    value -= csol[ dict_item[ i ] ] * obj[ dict_item[ i ] ];

  // add the easy components cost part and the lower bounds contribution of
  // the difficult components  - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 1 ; i <= NrFi ; i++ )
   if( weasy[ i ] ) {
    if( FiLambda[ i - 1 ] < HpINF ) {
     for( Index j = 0; j < FIO->GetBNC( i ) ; j++ )
      value -= csol[ comp_col[ i ] + j ] * obj[ comp_col[ i ] + j ];
     value += FiLambda[ i - 1 ];
     }
    else
     value = HpINF;
    }
   else
    value -= csol[ comp_col[ i ] ] * obj[ comp_col[ i ] ];

  // divide by rho - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( csol[ comp_col[ 0 ] ] != 0 )
   value /= csol[ comp_col[ 0 ] ];

  // If a global LB is present, we have LB * r in the objective function; but
  // we don't actually have r in the problem, rather rho = 1 - r. Therefore
  // we have LB * r = LB ( 1 - rho ) = - LB * rho + LB, and that's why in
  // the objective function we have - LB. But this also means that
  // *a fixed term LB has to be added to Sigma* if the LB is present

  if( HasLwrBnd )
   value += obj[ comp_col[ 0 ] ] * ( 1 - csol[ comp_col[ 0 ] ] );

  // finally add the contribution of the slack variables - - - - - - - - - -
  // - - - - - -  - - - - - - -  - - - - - - - - - - - - - - - - - - - - - -

  if( stab != boxstep )
   for( Index i = 0; i < CrrSGLen ; i++ ) {
    int offset = 0;  // if the lower bound exists, there is a slack s_i^+
    if( Lower[ i ] > -LMINF ) {
     value -= Lower[ i ] * csol[ dict_slack[ i ] ];
     offset = 1;
     }

    if( Upper[ i ] < LMINF )
     value += Upper[ i ] * csol[ dict_slack[ i ] + offset ];
    }
  else
   for( Index i = 0; i < CrrSGLen ; i++ ) {
    // if -t < Lower[ i ], the stabilization constraint -t <= d[ i ] is
    // redundant, and the opposite for Upper[ i ]
    if( -t < Lower[ i ] )
     value -= Lower[ i ] * csol[ dict_slack[ i ] ];
    if( t > Upper[ i ] )
     value += Upper[ i ] * csol[ dict_slack[ i ] + 1 ];
    }

  MSG( 0 , "Sigma* + Sigma_L = " << value << endl );
  }
 else {
  if( weasy[ wFi ] ) {
   if( FiLambda[ wFi - 1 ] < HpINF ) {
    for( Index j = 0; j < FIO->GetBNC( wFi ) ; j++ )
     value -= csol[ comp_col[ wFi ] + j ] * obj[ comp_col[ wFi ] + j ];
    value += FiLambda[ wFi - 1];
    }
   else
    value = HpINF;
  }
  else {
   for( Index i = 0 ; i < dict_item_maxname ; i++ )
    if( ( WComponent( i ) == wFi ) && WasInMP( i ) )
     value -= csol[ dict_item[ i ] ] * obj[ dict_item[ i ] ];
   value -= csol[ comp_col[ wFi ] ] * obj[ comp_col[ wFi ] ];
   }

  if( csol[ comp_col[ 0 ] ] != 0 )
   value /= csol[ comp_col[ 0 ] ];
  MSG( 0 , "Sigma*(" << wFi << ") = " << value << endl );
  }

 return( value );

 }  // end( OSIMPSolver::ReadSigma )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadDStart( cHpNum tt )
{
 MSG( 0 , "OSIMPSolver::ReadDStart()" << endl );

 HpNum total = 0;
 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 switch( stab ) {
  case none: break;  // note: we assume the solution to be feasible
  case boxstep:
   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( -t >= Lower[ i ] )
     total += tt * csol[ dict_slack[ i ] ];

    if( t <= Upper[ i ] )
     total += tt * csol[ dict_slack[ i ] + 1 ];
    }
   break;
  case quadratic:
   for( Index i = 0 ; i < CrrSGLen ; i++ )
    total += csol[ dict_stab[ i ] ] * csol[ dict_stab[ i ] ] ;
   total *= (tt/2);
   break;
  default:
   throw( NDOException( "OSIMPSolver::ReadDStart: undecided stabilization" )
	  );
  }

 MSG( 0 , "D* = " << total << endl );
 return( total );

 }  // end( OSIMPSolver::ReadDStart )

/*--------------------------------------------------------------------------*/

cLMRow OSIMPSolver::Readd( bool Fulld )
{
 MSG( 0 , "OSIMPSolver::Readd()\n" );

 if( Fulld ) {
  #if( PRESERVE_OSI_SOLS == 0 )
   const double *rsol = osiSlvr->getRowPrice();
  #endif

  return( rsol + comp_row[ NrFi ] );
  }
 else
  throw( NDOException( "OSIMPSolver::Readd: Fulld( false )" ) );

 }  // end( OSIMPSolver::Readd )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ReadZ( LMRow tz , cIndex_Set &I , Index &D , cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::ReadZ() \n");

 // give it out in "dense" format  - - - - - - - - - - - - - - - - - - - - -

 D = CrrSGLen;
 I = 0;

 if( wFi == 0 ) { // return the gradient of 0-th component - - - - - - - - -
  VectAssign( tz , LMRow(0) , CrrSGLen );
  for( Index i = 0 ; i < RhoColBDm  ; i++ )
   if( Index( RhoColBse[ i ] ) >= comp_row[ NrFi ] )
    tz[ RhoColBse[ i ] - comp_row[ NrFi ] ] = -RhoCol[ i ];
  return;
  }

 // make a convex combination of subgradient  - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 // use CoinPackedVector to simplify the work - - - - - - - - - - - - - - - -

 CoinPackedVector temp;
 if( wFi == 0 )    // return the gradient of 0-th component - - - - - - - - -
  temp = ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
 else
  if( wFi > NrFi ) {  // take in account all the subgradients - - - - - - - -

   for( Index i = 0 ; i <= dict_item_maxname ; i++ )
    if( ( IsSubG( i ) ) && WasInMP( i ) )
     temp = temp + csol[ dict_item[ i ] ] *
                  ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

   // add the sugradient of the easy components - - - - - - - - - - - - - - -

   for( Index i = 1 ; i <= NrFi ; i++ )
    if( weasy[ i ] )
     for( int j = comp_col[ i ] ; j < comp_col[ i ] + FIO->GetBNC( i ) ; j++ )
      temp = temp + csol[ j ] * ( osiSlvr->getMatrixByCol() )->getVector( j );

   if( wFi == InINF )   // add the linear part
    temp = temp + ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
   }
  else { // disaggregated case  - - - - - - - - - - - - - - - - - - - - - - -
   if( weasy[ wFi ] )
    throw( NDOException( "OSIMPSolver::ReadZ: calling for easy component" ) );

   for( Index i = 0 ; i < dict_item_maxname ; ++i )
    if( IsSubG( i ) && ( WComp( i ) == wFi ) && WasInMP( i ) )
     temp = temp + csol[ dict_item[ i ] ] *
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );
   }

 // assign tz the subgradient  - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( tz , LMNum( 0 ) , CrrSGLen );

 for( int j = 0 ; j < temp.getNumElements() ; j++ )
  if( temp.getIndices()[ j ] >= int( comp_row[ NrFi ] ) )
   tz[ temp.getIndices()[ j ] - comp_row[ NrFi ] ] = -temp.getElements()[ j ];

 }  // end( OSIMPSolver::ReadZ )

/*--------------------------------------------------------------------------*/

cHpRow OSIMPSolver::ReadMult( cIndex_Set &I , Index &D , cIndex wFi ,
			      const bool IncldCnst )
{
 MSG( 0 , "OSIMPSolver::ReadMult()\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ReadMult for 0-th component" ) );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 HpNum rho = csol[ comp_col[ 0 ] ];

 if( rho == 0 ) {  // in the case rho = 0, return a null vector- - - - - - -
  resizeI( 1 );
  tempI[ D = 0 ] = InINF;
  I = reinterpret_cast<cIndex_Set>( tempI );
  return( tempHP );
  }

 if( wFi <= NrFi )
  if( weasy[ wFi ] ) {  // easy component part- - - - - - - - - - - - - - -
   I = 0;
   D = FIO->GetBNC( wFi );
   resizeHP( D );
   VectAssign( tempHP , csol + comp_col[ wFi ] , 1 / rho , D );
   }
  else {                // difficult component part - - - - - - - - - - - -
   resizeHP( MaxBSize );
   resizeI( MaxBSize + 1 );
   I = reinterpret_cast<cIndex_Set>( tempI );
   D = 0;
   if( IncldCnst ) {
    for( Index name = 0 ; name < dict_item_maxname ; name++ )
     if( ( WComponent( name ) == wFi ) && WasInMP( name ) )
      if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
       tempI[ D++ ] = name;
    }
   else
    for( Index name = 0 ; name < dict_item_maxname ; name++ )
     if( IsSubG( name ) && ( WComp( name ) == wFi ) && WasInMP( name ) )
      if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
       tempI[ D++ ] = name;

   tempI[ D ] = InINF;
   }
 else {  // wFi == InINF- - - - - - - - - - - - - - - - - - - - - - - - - -
  resizeHP( MaxBSize );
  resizeI( MaxBSize + 1 );
  I = reinterpret_cast<cIndex_Set>( tempI );
  D = 0;
  if( IncldCnst ) {
   for( Index name = 0 ; name < dict_item_maxname ; name++ )
    if( ( dict_item[ name ] < InINF ) && WasInMP( name ) )
     if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
      tempI[ D++ ] = name;
   }
  else
   for( Index name = 0 ; name < dict_item_maxname ; name++ )
    if( IsSubG( name ) && WasInMP( name ) )
     if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
      tempI[ D++ ] = name;

  tempI[ D ] = InINF;
  }

 return( tempHP );

 }  // end( OSIMPSolver::ReadMult )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadLBMult( cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::ReadLBMult()\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ReadLBMult for 0-th component" ) );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 HpNum value;
 if( wFi <= NrFi ) {  // return gamma_i
  if( weasy[ wFi ] )
   throw( NDOException( "OSIMPSolver::ReadLBMult for an easy component" ) );

  value = csol[ comp_col[ wFi ] ];
  }
 else  // return r = ( 1 - r ) - - - - - - - - - - - - - - - - - - - - - - -
  value = 1 - csol[ comp_col[ 0 ] ];

 return( value );
 }

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadGid( cIndex Nm )
{
 MSG( 0 , "OSIMPSolver::ReadGid() called\n" );

 HpNum value = 0;
 if( Nm >= MaxBSize ) {  // 0-th component: return delta * b
  #if( PRESERVE_OSI_SOLS == 0 )
   const double *rsol = osiSlvr->getRowPrice();
  #endif
  for( Index i = comp_row[ NrFi ] ; i < RhoColBDm ; i++ )
   value -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];
  MSG( 0 , "b * delta = " << value << endl );
  }
 else {
  if( dict_item[ Nm ] == InINF )
   throw( NDOException( "OSIMPSolver::ReadGid: unused item name" ) );

  // note: the scalar product is available for all items, even those that
  //       were not present in the last master problem (i.e., such that
  //       wcomp[ Nm ] & LLBIndex == true) because it is computed in
  //       Check[SubG/Const]() and saved in SetItem().

  value = GiPerd[ Nm ];

  MSG( 0 , "Item[" << Nm << "] * d = " << value << endl );
  }

 return( value );

 }  // end( OSIMPSolver::ReadGid )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau )
{
 MSG( 0, "OSIMPSolver::MakeLambda1()\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 if( useactiveset )   // Lambda and Lambda1 are in sparse format- - - - - -
  for( Index i = 0 ; i < Asetdim ; i++ ) {
   HpNum step = ( Tau / t ) * primal[ Aset[ i ] ];
   if( step > Upper[ i ] ) step = Upper[ Aset[ i ] ];
   if( step < Lower[ i ] ) step = Lower[ Aset[ i ] ];
   Lmbd1[ i ] = Lmbd[ i ] + step;
   }
 else   // Lambda and Lambda1 are in dense format- - - - - - - - - - - - - -
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   HpNum step = ( Tau / t ) * primal[ i ];
   if( step > Upper[ i ] ) step = Upper[ i ];
   if( step < Lower[ i ] ) step = Lower[ i ];
   Lmbd1[ i ] = Lmbd[ i ] + step;
   }

 }  // end( OSIMPSolver::MakeLambda1 )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SensitAnals( HpNum &lp , HpNum &cp )
{
 throw( NDOException( "OSIMPSolver::SensitAnals not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE PROBLEM ----------------*/
/*--------------------------------------------------------------------------*/

Index OSIMPSolver::BSize( cIndex wFi )
{
 if( ( wFi > 0 ) && ( wFi <= NrFi ) )
  return( NSubG[ wFi ] + NConst[ wFi ] );
 else
  return( NSubG[ 0 ] + NConst[ 0 ] );

 }  // end( OSIMPSolver::BSize )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::BCSize( cIndex wFi )
{
 if( ( wFi > 0 ) && ( wFi <= NrFi ) )
  return( NConst[ wFi ] );
 else
  return( NConst[ 0 ] );

 }  // end( OSIMPSolver::BCSize )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::MaxName( cIndex wFi )
{
 Index i = dict_item_maxname;
 if( wFi != InINF )
  for( ; i > 0 ; i-- )
   if( WComponent( i - 1 ) == wFi )
    break;

 return( i );

 }  // end( OSIMPSolver::MaxName )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::WComponent( cIndex i )
{
 return( wcomp[ i ] < InINF ? WComp( i ) : InINF );

 }  // end( OSIMPSolver::WComponent )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::IsSubG( cIndex i )
{
 return( wcomp[ i ] < InINF ? bool( wcomp[ i ] & LBIndex ) : false );

 }  // end( OSIMPSolver::IsSubG )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::NumNNVars( void )
{
 return( NNVars );

 }  // end( OSIMPSolver::NumNNVars )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::NumBxdVars( void )
{
 return( BxdVars );

 }  // end( OSIMPSolver::NumBxdVars )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::IsNN( cIndex i )
{
 return( Lower[ i ] > -LMINF );

 }  // end( OSIMPSolver::IsNN )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::CheckIdentical( const bool Chk )
{
 checkID = Chk;

 }  // end( OSIMPSolver::CheckIdentical )

/*--------------------------------------------------------------------------*/

cHpRow OSIMPSolver::ReadLinErr( void )
{
 resizeHP( MaxBSize );
 VectAssign( tempHP , 0.0 , MaxBSize );

 const double *linerr = osiSlvr->getObjCoefficients();
 for( Index i = 0 ; i < MaxBSize ; i++ )
  if( dict_item[ i ] < InINF )
   tempHP[ i ] = - linerr[ dict_item[ i ] ];  // change sign!

 return( tempHP );

 }  // end( OSIMPSolver::ReadLinErr )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::EpsilonD( void )
{
 return( FsbEps );

 }  // end( OSIMPSolver::EpsilonD )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::FiBLambdaIsExact( cIndex wFi )
{
 return( weasy[ wFi ] );

 }  // end( OSIMPSolver::FiBLambdaIsExact )

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR ADDING / REMOVING / CHANGING DATA ---------------*/
/*--------------------------------------------------------------------------*/

SgRow OSIMPSolver::GetItem( cIndex wFi )
{
 NewItemFi = wFi;  // mark the name of the component
 return( NewItem );

 }  // end( OSIMPSolver::GetItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetItemBse( cIndex_Set SGBse , cIndex SGBDm )
{
 if( ( NewItemBse = SGBse ) )
  NewItemBDm = SGBDm;    // sparse format
 else
  NewItemBDm = CrrSGLen; // dense format

 }  // end( OSIMPSolver::SetItemBse )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai ,
			      HpNum &ScPri )
{
 MSG( 0 , "OSIMPSolver::CheckSubG() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 // check if the subgradient is identical to any other in the bundle - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = CheckBCopy();

 // compute the scalar product Gi^{top} delta- - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NewItemBse )  // sparse format
  ScPri = ScalarProduct( NewItem , primal , NewItemBse );
 else              // dense format
  ScPri = ScalarProduct( NewItem , primal , CrrSGLen );

 // compute the linearization error of the new subgradient - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Ai = Ai - ( DFi - ( Tau / t ) * ScPri );

 NewItemprice = Ai;
 NewItemScPri = ScPri;
 NewItemisSG = true;   // this item is a subgradient ...

 return( IsIde );

 } // end( OSIMPSolver::CheckSubG )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt )
{
 MSG( 0 , "OSIMPSolver::CheckCnst() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 // check if the constraint is identical to any other in the bundle- - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = CheckBCopy();

 // compute ScPri and right hand side of the constraint- - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NewItemBse ) {  // sparse format- - - - - - - - - - - - - - - - - - -
  ScPri = ScalarProduct( NewItem , primal , NewItemBse );

  // compute the right hand side of the constraint
  if( Aset )         // newitem is sparse & is using the active set
   Ai -= ScalarProduct( NewItem , NewItemBse , CrrPnt , Aset );
  else               // newitem is sparse & is not using the active set
   Ai -= ScalarProduct( NewItem , CrrPnt , NewItemBse );
  }
 else {              // dense format - - - - - - - - - - - - - - - - - - -
  ScPri = ScalarProduct( NewItem , primal , CrrSGLen );

  // compute the right hand side of the constraint
  if( Aset )         // newitem is dense & activeset
   Ai -= ScalarProduct( CrrPnt , NewItem , Aset );
  else               // newitem is dense & not activeset
   Ai -= ScalarProduct( NewItem , CrrPnt , CrrSGLen );
  }

 NewItemprice = Ai;
 NewItemScPri = ScPri;
 NewItemisSG = false;  // this item is a constraint ...

 return( IsIde );

 }  // end( OSIMPSolver::CheckCnst )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::ChangesMPSol( void )
{
 MSG( 0 , "OSIMPSolver::ChangesMPSol() called\n" );

 if(! NewItemFi )
  throw( NDOException(
   "OSIMPSolver::ChangesMPSol: a calling to 0-th is not allowed" ) );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 bool ChgMP = true;

 if( NewItemisSG ) {
  /* check whether or not the new subgradient provides a cut in Master
     Problem, that is verify v >= G1 * dir - Alfa1 is violated */

  if( NewItemFi <= NrFi ) {
   if( weasy[ NewItemFi ] )
    throw( NDOException(
     "OSIMPSolver::ChangesMPSol: check for an easy comp. is not allowed" )
     );
   }
  else
   throw( NDOException( "OSIMPSolver::ChangesMPSol: "
		  "check for an aggregated subgradient is not allowed" ) );

  if( rsol[ comp_row[ NewItemFi - 1 ] ] + NewItemprice - NewItemScPri >=
	 - FsbEps * max( ABS( NewItemprice ) , double(1) ) )
   ChgMP = false;

  }
 else // check whether or not the new constraint provides a cut
      // in Master Problem, that is verify Alfa1 >= G1 * dir is violated
  if( NewItemprice - NewItemScPri >=
      - FsbEps * max( ABS( NewItemprice ) , double( 1 ) ) )
   ChgMP = false;

 return( ChgMP );

 } // end( OSIMPSolver::ChangesMPSol )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetItem( cIndex Nm )
{
 // note: the columns of the coefficient matrix (for the relevant part)
 //       contain the *opposite* of the subgradient

 MSG( 0 , "OSIMPSolver::SetItem() \n");

 if( ( ! NewItemFi ) && ( Nm == InINF ) ) {  // the 0th component
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // we have to update the rho column, in particular the part corresponding
  // to the subgradient
  RhoColBDm = comp_row[ NrFi ];

  if( NewItemBse )  // sparse format - - - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < NewItemBDm ; i++ ) {
    RhoColBse[ comp_row[ NrFi ] + i ] = comp_row[ NrFi ] + NewItemBse[ i ];
    RhoCol[ comp_row[ NrFi ] + i ] = - NewItem[ i ];
    RhoColBDm++;
    }
  else 	            // dense format  - - - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    RhoColBse[ comp_row[ NrFi ] + i ] = comp_row[ NrFi ] + i;
    RhoCol[ comp_row[ NrFi ] + i ] = - NewItem[ i ];
    RhoColBDm++;
    }

  int col = comp_col[ 0 ];  // get the name of the column of rho - - - - - -
  HpNum coeff = osiSlvr->getObjCoefficients()[ col ]; // and the lower bound

  osiSlvr->deleteCols( 1 , &col ); // delete the previous rho column

  // update the dictionaries - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( ( dict_item[ i ] < InINF ) &&
       ( dict_item[ i ] > Index( col ) ) )
    dict_item[ i ]--;

  for( Index i = 0 ; i < CrrSGLen ; i++ )
   if( ( dict_slack[ i ] < InINF ) &&
       ( dict_slack[ i ] > Index( col ) ) )
    dict_slack[ i ]--;

  if( dict_stab ) {
   for( Index i = 0 ; i < CrrSGLen ; i++ )
    if( ( dict_stab[ i ] < InINF ) &&
	( dict_stab[ i ] > Index( col ) ) )
     dict_stab[ i ]--;
   }

  // add the rho column to the master problem  - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_col[ 0 ] = osiSlvr->getNumCols();
  osiSlvr->addCol( RhoColBDm , RhoColBse , RhoCol , 0.0 ,
		   osiSlvr->getInfinity() , coeff );
  }
 else {  // it's an actual item- - - - - - - - - - - - - - - - - - - - - - -

  resizeI( NewItemBDm + 1 );   // note: NewItemBD == CrrSGLen in the dense
  resizeHP( NewItemBDm + 1 );  // case

  // copy the item in the temporary structure  - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index i = 0;
  if( NewItemBse )
   for( ; i < NewItemBDm ; i++ ) {
    tempI[ i ] = NewItemBse[ i ] + comp_row[ NrFi ];
    tempHP[ i ] = - NewItem[ i ];
    }
  else
   for( ; i < CrrSGLen ; i++ ) {
    tempI[ i ] = i + comp_row[ NrFi ];
    tempHP[ i ] = - NewItem[ i ];
    }

  if( NewItemisSG ) {  // add the multiplier to the simplex set

   tempI[ i ] = comp_row[ NewItemFi - 1 ];
   tempHP[ i++ ] = 1.0;
   NSubG[ NewItemFi ]++;
   (*NSubG)++;

   // and possibly fix gamma[ NewItemFi ] = 0- - - - - - - - - - - - - - - -

   if( ( osiSlvr->getColUpper()[ comp_col[ NewItemFi ] ] > 0.0 ) &&
       ( osiSlvr->getObjCoefficients()[ comp_col[ NewItemFi ] ] == 0.0 ) )
    osiSlvr->setColUpper( comp_col[ NewItemFi ] , 0.0 );

   wcomp[ Nm ] = NewItemFi | IMask;  // set name, plus it's a "new" item,
                                     // plus is a subgradient
   }
  else {
   NConst[ NewItemFi ]++;
   (*NConst)++;

   wcomp[ Nm ] = NewItemFi | LLBIndex;  // set name, plus it's a "new" item
   }

  // also record the scalar product for future use (ChangeCurrPoint) - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  GiPerd[ Nm ] = NewItemScPri;

  // add the column to the master probem - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->addCol( int( i ) , tempI , tempHP , 0.0 , osiSlvr->getInfinity() ,
		   - NewItemprice );

  // put the item in the dictionary  - - - - - - - - - - - - - - - - - - - -

  dict_item[ Nm ] = osiSlvr->getNumCols() - 1;

  MSG( 0 , "Item " << Nm << " of the component " << NewItemFi
	<< " has the entry " << dict_item[ Nm ] << endl);

  if( dict_item_maxname < Nm + 1 )
   dict_item_maxname = Nm + 1;
  }
 }  // end( OSIMPSolver::SetItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SubstItem( cIndex Nm )
{
 MSG( 0 , "OSIMPSolver::SubstItem()\n" );

 //if( - osiSlvr->getObjCoefficients()[ dict_item[ Nm ] ] > NewItemprice )
 osiSlvr->setObjCoeff( dict_item[ Nm ] , - NewItemprice );

 } // end( OSIMPSolver::SubstItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvItem( cIndex i )
{
 MSG( 0 , "OSIMPSolver::RmvItem()\n" );

 const int index = dict_item[ i ];  // mark the column deleted - - - - - - -
 osiSlvr->deleteCols( 1 , &index );

 #if( PRESERVE_OSI_SOLS )
  // having deleted the column 'index' from the osiSlvr, all the stored
  // solution information corresponding to columns (csol[] and rcst[]) must
  // be shifted left by one from position index onwards

  ShiftVect( csol + index , csols - index - 1 );
  ShiftVect( rcst + index , csols - index - 1 );
  csols--;
 #endif

 Index wFi = WComp( i );

 if( IsSubG( i ) ) {
  NSubG[ wFi ]--;
  (*NSubG)--;
  }
 else {
  NConst[ wFi ]--;
  (*NConst)--;
  }

 // update the item's vocabulary - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 dict_item[ i ] = wcomp[ i ] = InINF;

 if( i + 1 == dict_item_maxname )
  dict_item_maxname--;

 // if the deleted item is the last one for its component, the convexity
 // constraint for that component has to be deactivated  - - - - - - - - - -

 if( ( ! ( NSubG[ wFi ] + NConst[ wFi ] ) ) &&
     ( osiSlvr->getColUpper()[ comp_col[ wFi ] ] == 0.0 ) ) {
  osiSlvr->setColUpper( comp_col[ wFi ] , osiSlvr->getInfinity() );
  MSG( 0 , "removed item " << i << " relative to the component " << wFi << endl );
  }

 // update the dictionaries - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
  if( ( dict_item[ Nm ] != InINF ) &&
      ( dict_item[ Nm ] > Index( index ) ) )
   dict_item[ Nm ]--;

 for( Index Nm = 0 ; Nm < CrrSGLen ; Nm++ ) {
  if( ( dict_slack[ Nm ] != InINF ) &&
      ( dict_slack[ Nm ] > Index( index ) ) )
   dict_slack[ Nm ]--;
  if( stab == quadratic )
   if( ( dict_stab[ Nm ] != InINF ) &&
       ( dict_stab[ Nm ] > Index( index ) ) )
    dict_stab[ Nm ]--;
  }

 if( comp_col[ 0 ] > Index( index ) )
  comp_col[ 0 ]--;

 }  // end( OSIMPSolver::RmvItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvItems( void )
{
 MSG( 0 , "OSIMPSolver::RmvItems()\n" );

 if( ! dict_item_maxname )  // there is nothing to do ...
  return;

 VectAssign( NSubG , Index( 0  ) , NrFi + 1 );
 VectAssign( NConst , Index( 0 ) , NrFi + 1 );

 int count = 0;
 resizeI( dict_item_maxname );

 for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
  if( dict_item[ Nm ] != InINF )
   tempI[ count++ ] = dict_item[ Nm ];

 // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->deleteCols( count , tempI );

 // fill the vectors for the shifting of the vocabularies  - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index elem;
 for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
  if( dict_item[ Nm ] != InINF ) {

   elem = dict_item[ Nm ];

   // update the left part of tempI  - - - - - - - - - - - - - - - - - - - -

   for( Index i = Nm + 1; i < dict_item_maxname ; i++ )
    if( ( dict_item[ i ] > elem ) && ( dict_item[ i ] != InINF ) )
     dict_item[ i ]--;

   // update the position of rho column  - - - - - - - - - - - - - - - - - -

   if( comp_col[ 0 ] > elem )
    comp_col[ 0 ]--;

   // update the slack and stabilization dictionaries  - - - - - - - - - - -

   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( ( dict_slack[ i ] != InINF ) && ( dict_slack[ i ] > elem ) )
     dict_slack[ i ]--;
    if( stab == quadratic )
     if( ( dict_stab[ i ] != InINF ) && ( dict_stab[ i ] > elem ) )
      dict_stab[ i ]--;
    }

   } // end check  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < NrFi ; i++ )
  if( osiSlvr->getColUpper()[ comp_col[ i ] ] == 0.0 )
   osiSlvr->setColUpper( comp_col[ i ] , osiSlvr->getInfinity() );

 VectAssign( dict_item , InINF , dict_item_maxname );
 VectAssign( wcomp , InINF , dict_item_maxname );
 dict_item_maxname = 0;

 }  // end( OSIMPSolver::RmvItems )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetActvSt( cIndex_Set AVrs , cIndex AVDm )
{
 MSG( 0 , "OSIMPSolver::SetActvSt()\n" );

 if( ! useactiveset )
  throw( NDOException( "OSIMPSolver::SetActvSt(): Active set not declared"
		       ) );

 if( ! AVrs ) {	 //  all the variables will be not active- - - - - - - - - -
  while( Asetdim > 0 )
   deactivate( Aset[ Asetdim-- ] );
  }
 else
  if( ! Aset ) {  // if Aset is empty, set all the variables in AVrs as active
    while ( Asetdim < AVDm )
     activate( AVrs[ Asetdim++ ] );
   }
  else
   for( cIndex *newA = AVrs ; ( *newA < InINF ) ||
	                      ( *Aset < InINF ) ; ) {
    if( *newA < *Aset )	{
     activate( *newA );
     newA++;
     }
    else
     if( *newA > *Aset ) {
      deactivate( *Aset );
      Aset++;
      }
     else {
      newA++;
      Aset++;
      }
    }  // end( for )

 Aset = AVrs;
 Asetdim = AVDm;

 }  // end( OSIMPSolver::SetActvSt )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs )
{
 MSG( 0 , "OSIMPSolver::AddActvSt() called\n" );

 if( ! useactiveset )
  throw( NDOException( "OSIMPSolver::AddActvSt(): Active set not declared"
		       ) );
 Aset = AVrs;
 Asetdim += AdDm;

 for( ; *Addd < InINF ; Addd++ )
  activate( *Addd );

 }  // end( OSIMPSolver::AddActvSt )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs )
{
 MSG( 0 , "OSIMPSolver::RmvActvSt() called\n" );

 Aset = AVrs;
 Asetdim = Asetdim - RmDm;

 for( ; *Rmvd < InINF ; Rmvd++ )
  deactivate( *Rmvd );

 } // end( OSIMPSolver::RmvActvSt )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::AddVars( cIndex NNwVrs )
{
 MSG( 0 , "OSIMPSolver::AddVars() called\n" );

 if( ! NNwVrs )  // adding 0 new variables
  return;        // all is done

 cIndex NewSGLen = CrrSGLen + NNwVrs;
 if( NewSGLen > MaxSGLen )
  throw( NDOException( "OSIMPSolver::AddVars: too many variables" ) );

 // allocate the memory for new rows - - - - - - - - - - - - - - - - - - - -
 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index size = 0;
 int * rowStarts = new int[ NNwVrs + 1 ];      // start point of the rows
 HpRow rowlb = new HpNum[ NNwVrs ];            // lower bound of constraints
 HpRow rowub = new HpNum[ NNwVrs ];            // upper bound of constraints

 // recovery the description of  matrix A[ wFi ] of the "easy" - - - - - - -
 // components - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SIndex_Mat Abeg = new SIndex_Set[ NrFi ];
 SIndex_Mat Aind = new SIndex_Set[ NrFi ];
 Mat Aval = new Row[ NrFi ];

 Index_Set BNC = new Index[ NrFi ];
 Index_Set ANZ = new Index[ NrFi ];

 for( Index j = 1 ; j <= NrFi ; j++)
  if( ( BNC[ j - 1 ] = FIO->GetBNC( j ) ) ) {
   ANZ[ j - 1 ] = FIO->GetANZ( j , CrrSGLen , NewSGLen );
   if( ANZ[ j - 1 ] ) {
	size += ANZ[ j - 1 ];
    Abeg[ j - 1 ] = new SIndex[ BNC[ j - 1 ] + 1 ];
    Aind[ j - 1 ] = new SIndex[ ANZ[ j - 1 ] ];
    Aval[ j - 1 ] = new Number[ ANZ[ j - 1 ] ];
    FIO->GetADesc( j , Abeg[ j - 1 ] , Aind[ j - 1 ] , Aval[ j - 1 ] ,
		   CrrSGLen , NewSGLen );
    }
   }
  else {
   Abeg[ j - 1 ] = Aind[ j - 1 ] = 0;
   Aval[ j - 1 ] = 0;
   ANZ[ j - 1 ] = 0;
   }

 // recovery all part of sub-gradients corresponding to variables in Lambda
 // from index CrrSGLen to (NewSGLen - 1)   - - - - - - - - - - - - - - - - -

 cIndex_Set SGBse;
 SgMat tempG = new SgRow[ dict_item_maxname ];

 for( Index j = 0 ; j < dict_item_maxname ; j++ )
  if( dict_item[ j ] < InINF ) {  // ask the components [ CrrSGLen ,
   tempG[ j ] = new SgNum[ NNwVrs ];     // NewSGLen ) for item j
   cIndex SGBDim = FIO->GetGi( tempG[ j ] , SGBse , j , CrrSGLen , NewSGLen );
   if( SGBse )
    Densify( tempG[ j ] - CrrSGLen , SGBse , SGBDim , NewSGLen , CrrSGLen );

   size += NNwVrs;
   }
  else
   tempG[ j ] = 0;

 // recovery the subgradient of the linear 0-th component of Fi  - - - - - -
 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SgRow tempG0 = new SgNum[ NNwVrs ];
 cIndex SGBDim = FIO->GetGi( tempG0 , SGBse , FIO->GetMaxName() , CrrSGLen ,
			     NewSGLen );
 if( SGBse )
  Densify( tempG0 - CrrSGLen , SGBse , SGBDim , NewSGLen , CrrSGLen );

 size += NNwVrs;

 // construct the rows to the problem, one after another   - - - - - - - - -
 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 resizeI( size );
 resizeHP( size );

 Index count = 0;
 for( Index i = CrrSGLen ; i < NewSGLen ; i++ ) {

  rowStarts[ i - CrrSGLen ] = count;

  // write the matrices A[ wFi ] in the problem, one for each easy component

  for( Index j = 1 ; j <= NrFi ; j++)
   if( ANZ[ j - 1 ] ) {
    for( Index k = 0 ; k < BNC[ j - 1 ] ; k++ )
     for( SIndex l =  Abeg[ j - 1 ][ k ]; l < Abeg[ j - 1 ][ k + 1 ] ; l++ )
      if( Index( Aind[ j - 1 ][ l ] ) == i ) {
       tempI[ count ] = k + comp_col[ j ];
       tempHP[ count++ ] = Aval[ j - 1 ][ l ];
       }
    }

  // write in both tempI and tempHP the subgradient part - - - - - - - - - -

  for( Index j = 0 ; j < dict_item_maxname ; j++ )
   if( dict_item[ j ] < InINF ) {
    tempI[ count ] = dict_item[ j ];
    tempHP[ count++ ] = - tempG[ j ][ i - CrrSGLen ];
    }

  // write in both tempI and tempHP the subgradient of the linear 0-th
  // component

  RhoColBse[ RhoColBDm ] = comp_row[ NrFi ] + i;     // update the column
  RhoCol[ RhoColBDm++ ] = - tempG0[ i - CrrSGLen ];  // of rho

  tempI[ count ] = comp_col[ 0 ];
  tempHP[ count++ ] = -tempG0[ i - CrrSGLen ];
  }

 rowStarts[ NNwVrs ] = count;

 // add the rows to problem - - - - - - - - - - - - - - - - - - - - - - - -

 if( useactiveset ) {
  MSG( 0 , "Add constraints for the primal variables as inactive ones\n" );
  for( Index j = 0 ; j < NNwVrs ;  j++ ) {
   rowlb[ j ] = - osiSlvr->getInfinity();
   rowub[ j ] = osiSlvr->getInfinity();
   }
  osiSlvr->addRows( NNwVrs , rowStarts, tempI , tempHP , rowlb , rowub );
  }
 else {
  MSG( 0 , "Add constraints for the primal variables \n" );
  for( Index j = 0 ; j < NNwVrs ;  j++ ) {
   rowlb[ j ] = 0;
   rowub[ j ] = 0;
   }
  osiSlvr->addRows( NNwVrs , rowStarts, tempI , tempHP , rowlb , rowub );
  }

 // deallocate the memory- - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] tempG0;
 for( Index j = 0 ; j < dict_item_maxname ; j++ )
   delete[] tempG[ j ];
 delete[] tempG;

 delete[] BNC;
 delete[] ANZ;

 for( Index j = 1 ; j <= NrFi ; j++ ) {
  delete[] Abeg[ j - 1 ];
  delete[] Aind[ j - 1 ];
  delete[] Aval[ j - 1 ];
  }

 delete[] Abeg;
 delete[] Aind;
 delete[] Aval;

 delete[] rowStarts;
 delete[] rowlb;
 delete[] rowub;

 // add the slacks to the problem  - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index nc = osiSlvr->getNumCols();  // number of columns - 1
 HpRow obj;                         // coefficients of the objective
                                    // function relative to slacks
 int *columnStarts;
 if( dict_stab ) {
  columnStarts = new int[ ( 3 * NNwVrs ) + 1 ];
  resizeI( 3 * NNwVrs );
  resizeHP( 3 * NNwVrs );
  rowlb = new HpNum[ 3 * NNwVrs ];
  }
 else {
  columnStarts = new int[ ( 2 * NNwVrs ) + 1 ];
  resizeI( 2 * NNwVrs );
  resizeHP( 2 * NNwVrs );
  rowlb = new HpNum[ 2 * NNwVrs ];
  }

 // redefine the size both of tempHp and tempI - - - - - - - - - - - - - - -

 bool slack_p;
 bool slack_m;

 count = 0;
 for( Index i = CrrSGLen ; i < NewSGLen ; i++ ) {
  slack_p = false;
  slack_m = false;

  if( ( Upper[ i ] = FIO->GetUB( i ) ) != LMINF )
   slack_m = true;

  if( FIO->GetUC( i ) )
   Lower[ i ] = -LMINF;
  else {
   Lower[ i ] = 0;
   slack_p = true;
   }

  if( slack_p )
   NNVars++;

  if( slack_p || slack_m )
   BxdVars++;

  // on the basis of the adopted stabilization type we have
  // either one or two slacks  - - - - - - - - - - - - - - - - - - - - - - -

  switch ( stab ) {
   case none:    break;
   case boxstep: slack_p = slack_m = true;
                 break;
   case quadratic:
	columnStarts[ count ] = count;
	dict_stab[ i ] = count + nc;
	rowlb[ count ] = - osiSlvr->getInfinity();
	tempI[ count ] = comp_row[ NrFi ] + i;
	tempHP[ count++ ] = 1.0;
	MSG( 0 , "Add z_"<< i << endl );
	break;
   default:
    throw( NDOException( "OSIMPSolver::AddVars: undecided stabilization" ) );
   }

  if( slack_p ) {
   columnStarts[ count ] = count;
   dict_slack[ i ] = count + nc;
   rowlb[ count ] = 0;
   tempI[ count ] = comp_row[ NrFi ] + i;
   tempHP[ count++ ] = 1.0;
   MSG( 0 , "Add the slack s_"<< i <<"+\n" );
   }

  if( slack_m )	{
   columnStarts[ count ] = count;

   // if slack_p is present, slack_m is not indicated
   if( ! slack_p )
    dict_slack[ i ] = count + nc;

   rowlb[ count ] = 0;
   tempI[ count ] = comp_row[ NrFi ] + i;
   tempHP[ count++ ] = -1.0;
   MSG( 0 , "Add the slack s_" << i << "-\n" );
   }
  }

 columnStarts[ count ] = count;  // set the end both of tempHp and tempI
 rowub = new HpNum[ count ];
 obj = new HpNum[ count ];
 for( Index j = 0 ; j < count ;  j++ ) {
  rowub[ j ] = osiSlvr->getInfinity();
  obj[ j ] = 0;
  }

 osiSlvr->addCols( count , columnStarts , tempI , tempHP , rowlb , rowub ,
		   obj );

 // deallocate the memory - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] columnStarts;
 delete[] rowlb;
 delete[] rowub;
 delete[] obj;

 // update the current items length - - - - - - - - - - - - - - - - - - - - -

 CrrSGLen = NewSGLen;

 // change the price of the slack variables - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // change the prices for the stabilization
 tUpdatePrices( CrrSGLen - NNwVrs , CrrSGLen ); 
 if( stab != boxstep )  // change the slack prices, unless they have just
  ptUpdatePrices( CrrSGLen - NNwVrs , CrrSGLen );  // been changed 

 }  // end( OSIMPSolver::AddVars )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvVars( cIndex_Set whch , Index hwmny )
{
 MSG( 0 , "OSIMPSolver::RmvVars() called\n" );

 int count;
 if( whch ) {  // removing a specific subset of variables - - - - - - - - - -
               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Index newCrrSGLen = CrrSGLen - hwmny;

  if( dict_stab )
   resizeI( 3 * hwmny );
  else
   resizeI( 2 * hwmny );

  // remove the rows in the dynamic part - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < hwmny ; i++ )
   tempI[ i ] = whch[ i ] + comp_row[ NrFi ];

  osiSlvr->deleteRows( hwmny , tempI );

  // update rho column  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  CoinPackedVector rhoCpy = 
                   ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
  for( Index j = rhoCpy.getNumElements() ; j-- > 0 ; ) {
   RhoCol[ j ] = rhoCpy.getElements()[ j ];
   RhoColBse[ j ] = rhoCpy.getIndices()[ j ];
   }
  RhoColBDm = rhoCpy.getNumElements();

  // mark the columns to delete   - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  count = 0;
  for( Index i = 0 ; i < hwmny ; i++ ) {
   cIndex wi = whch[ i ];
   bool islwr =  ( Lower[ wi ] > -LMINF );
   bool isupr = ( Upper[ wi ] < LMINF );

   if( islwr )
    NNVars--;

   if( islwr || isupr )
    BxdVars--;

   if( dict_slack[ wi ] != InINF ) {
    tempI[ count++ ] = dict_slack[ wi ];
    if( ( islwr && isupr ) || ( stab == boxstep ) )
     tempI[ count++ ] = dict_slack[ wi ] + 1;
    }
   if( dict_stab )
    tempI[ count++ ] = dict_stab[ wi ];
   }

  // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->deleteCols( count , tempI );

  // shift the vocabularies - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index elem;
  for( int Nm = 0 ; Nm < count ; Nm++ ) {
   elem = tempI[ Nm ];

   // update the left part of tempI  - - - - - - - - - - - - - - - - - - - -

   for( int i = Nm + 1; i < count ; i++ )
    if( tempI[ i ] > elem )
     tempI[ i ]--;

   // update the position of rho column  - - - - - - - - - - - - - - - - - -

   if( comp_col[ 0 ] > elem )
    comp_col[ 0 ]--;

   // update the items' dictionary  - - - - - - - - - -  - - - - - - - - - -

   for( Index i = 0 ; i < dict_item_maxname ; i++ )
    if( ( dict_item[ i ] != InINF ) && ( dict_item[ i ] > elem ) )
     dict_item[ i ]--;

   // update the slack and stabilization dictionaries  - - - - - - - - - - -

   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( ( dict_slack[ i ] != InINF ) && ( dict_slack[ i ] > elem ) )
     dict_slack[ i ]--;
    if( stab == quadratic )
     if( ( dict_stab[ i ] != InINF ) && ( dict_stab[ i ] > elem ) )
      dict_stab[ i ]--;
    }
   }  // end for- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // compact the vocabulary - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Compact( Upper , whch , CrrSGLen );
  Compact( Lower , whch , CrrSGLen );
  Compact( dict_slack , whch , CrrSGLen );
  if( dict_stab )
   Compact( dict_stab , whch , CrrSGLen );

  CrrSGLen = newCrrSGLen;
  }
 else {  // removing all variables- - - - - - - - - - - - - - - - - - - - - -
         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // remove all rows in the dynamic part- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < CrrSGLen ; i++ )
   tempI[ i ] = comp_row[ NrFi ] + i;

  osiSlvr->deleteRows( CrrSGLen , tempI );

  // mark the columns to delete - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( dict_stab )
   resizeI( 3 * CrrSGLen + dict_item_maxname );
  else
   resizeI( 2 * CrrSGLen + dict_item_maxname );

  count = 0;
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   if( dict_slack[ i ] != InINF ) {
	tempI[ count++ ] = dict_slack[ i ];
    if( ( Lower[ i ] > -LMINF ) &&  ( Upper[ i ] < LMINF ) )
     tempI[ count++ ] = dict_slack[ i ] + 1;
    }
    if( dict_stab )
     tempI[ count++ ] = dict_stab[ i ];
   }

  for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
   if( dict_item[ Nm ] != InINF )
    tempI[ count++ ] = dict_item[ Nm ];

  // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->deleteCols( count , tempI );

  // update the 0-th component- - - - - - - - - - - - - - - - - - - - - - - -

  RhoColBDm = comp_row[ NrFi ];

  for( Index Nm = 0 ; Nm < Index(count) ; Nm++ )
    if( comp_col[ 0 ] > tempI[ count ] )
     comp_col[ 0 ]--;

  // delete all the items- - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < NrFi ; i++ )
   if( osiSlvr->getColUpper()[ comp_col[ i ] ] == 0.0 )
    osiSlvr->setColUpper( comp_col[ i ] , osiSlvr->getInfinity() );

  VectAssign( dict_item , InINF , dict_item_maxname );
  VectAssign( wcomp , InINF , dict_item_maxname );
  VectAssign( dict_slack , InINF , MaxSGLen );

  dict_item_maxname = 0;

  CrrSGLen = NNVars = BxdVars = 0;

  } // end removing all elements - - - - - - - - - - - - - - - - - - - - - -
 }  // end( OSIMPSolver::RmvVars )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgAlfa( cHpRow DeltaAlfa )
{
 // note: the obj coefficients of each item column are *the opposite* of
 //       the corresponding linearization error / RHS

 MSG( 0 , "OSIMPSolver::ChgAlfa( DeltaAlfa ) called\n" );

 int numc = osiSlvr->getNumCols();
 resizeHP( numc );
 VectAssign( tempHP , osiSlvr->getObjCoefficients() , numc );

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( IsSubG( i ) )
   tempHP[ dict_item[ i ] ] -= DeltaAlfa[ WComp( i ) ];
   // note the sign: the linearization error has to be *increased* by
   // DeltaAlfa[ WComp( i ) ], but since the objective function coefficients
   // are the *opposite* of the linearization error, one has to subtract

 osiSlvr->setObjective( tempHP );

 }  // end( OSIMPSolver::ChgAlfa( cHpRow ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgAlfa( cHpRow NewAlfa , cIndex wFi )
{
 // note: the obj coefficients of each item column are *the opposite* of
 //       the corresponding linearization error / RHS

 MSG( 0 , "OSIMPSolver::ChgAlfa( NewAlfa ) called\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ChgAlfa( * , 0 ) called" ) );

 if( wFi > NrFi ) {
  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( dict_item[ i ] < InINF )
    osiSlvr->setObjCoeff( dict_item[ i ] , - NewAlfa[ i ] );
  }
 else
  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( WComponent( i ) == wFi )
    osiSlvr->setObjCoeff( dict_item[ i ] , - NewAlfa[ i ] );

 }  // end( OSIMPSolver::ChgAlfa( cHpRow , cIndex ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChangeCurrPoint( cLMRow DLambda , cHpRow DFi )
{
 /* Note: the formula for updating the linearization error of subgradient
          i belonging to component k is

     newalfa[ i ] = alfa[ i ] + DFi[ k ] - DLambda * z[ i ] .

    (and that of a constraint is the same with DFi[ k ] replaced with 0).
    However, in OSIMPSolver:

    - the linearization error / RHS of the column corresponding to z[ i ]
      is *the opposite* of the corresponding objective function coefficient

    - the column itself contains - z[ i ]

    so beware of the sign. */

 MSG( 0 , "OSIMPSolver::ChangeCurrPoint( DLabda ) called\n" );

 // change the coefficients of the objective function- - - - - - - - - - - -

 cIndex numcols = osiSlvr->getNumCols();
 resizeHP( numcols );
 VectAssign( tempHP , osiSlvr->getObjCoefficients() , numcols );

 // update the linearization error/rhs of the items- - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( dict_item[ i ] < InINF ) {
   double temp = - tempHP[ dict_item[ i ] ];
   CoinShallowPackedVector item =
     ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

   int k;
   for( Index j = item.getNumElements() ; j-- > 0 ; ) {
    if( ( k = item.getIndices()[ j ] - int( comp_row[ NrFi ] ) ) >= 0 )
     temp += item.getElements()[ j ] * DLambda[ k ];
    }

   if( IsSubG( i ) )
    temp += DFi[ WComp( i ) ];
   else {
    // it is a constraint: ensure its RHS remains non-negative, as small
    // negative RHS which may crop up by numerical errors would make
    // the primal master problem unfeasible (the dual unbounded)

    if( temp < 0 )
     temp = 0;
    }

   tempHP[ dict_item[ i ] ] = - temp;
   }

 // change the cost of the easy components and the lower bound of the- - - -
 // difficult ones - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i++ < NrFi ; )
  if( weasy[ i ] ) {
   cIndex BNC = FIO->GetBNC( i );
   for( Index j = 0 ; j < BNC ; j++ ) {
    double temp = tempHP[ comp_col[ i ] + j ];
    CoinShallowPackedVector item =
     ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ i ] + j );

    for( Index l = item.getNumElements() ; l > 0 ; l-- ) {
     int k;
     if( ( k = item.getIndices()[ l - 1 ] - int( comp_row[ NrFi ] ) ) >= 0 )
      temp -= item.getElements()[ l - 1 ] * DLambda[ k ];
     }

    tempHP[ comp_col[ i ] + j ] = temp;
    }
   }
  else
   if( osiSlvr->getColUpper()[ comp_col[ i ] ] > 0 )
    tempHP[ comp_col[ i ] ] -= DFi[ i ];

 // update the global lower bound- - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( HasLwrBnd )  // ... if any
  tempHP[ comp_col[ 0 ] ] += DFi[ 0 ];

 // change the bounds of the variables - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < CrrSGLen ; i++ ) {
  if( Upper[ i ] < LMINF )
   Upper[ i ] -= DLambda[ i ];
  if( Lower[ i ] > -LMINF )
   Lower[ i ] -= DLambda[ i ];
  }

 // update the slack prices - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ptUpdatePricesInPlace();

 // finally, change all costs in one blow - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->setObjective( tempHP );

 // set the value of the easy components in the current point - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( FiLambda , HpINF , NrFi );

 }  // end( OSIMPSolver::ChangeCurrPoint( DLambda , DFi ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChangeCurrPoint( cHpNum Tau , cHpRow DFi )
{
 MSG( 0 , "OSIMPSolver::ChangeCurrPoint( Tau ) called\n" );

 cIndex numcols = osiSlvr->getNumCols();
 resizeHP( numcols );
 VectAssign( tempHP , osiSlvr->getObjCoefficients() , numcols );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *delta = rsol + comp_row[ NrFi ];

 // update the linearization error/rhs of the items- - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( dict_item[ i ] < InINF ) {
   double temp = - tempHP[ dict_item[ i ] ];
   temp -= ( Tau / t ) * GiPerd[ i ];

   if( IsSubG( i ) )
    temp += DFi[ WComp( i ) ];
   else {  // it is a constrtaint: ensure its RHS remains non-negative, as
    // small negative RHS which may crop up by numerical errors may make
    // the primal master problem unfeasible (the dual unbounded)

    if( temp < 0 )
     temp = 0;
    }

   tempHP[ dict_item[ i ] ] = - temp;
   }

 // change the cost of the easy components and the lower bound of the- - - -
 // difficult ones - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i++ < NrFi ; )
  if( weasy[ i ] ) {
   cIndex BNC = FIO->GetBNC( i );
   for( Index j = 0 ; j < BNC ; j++ ) {
    double temp = tempHP[ comp_col[ i ] + j ];
    CoinShallowPackedVector item =
     ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ i ] + j );

    for( Index l = item.getNumElements() ; l > 0 ; l-- ) {
     int k;
     if( ( k = item.getIndices()[ l - 1 ] - int( comp_row[ NrFi ] ) ) >= 0 )
      temp -= item.getElements()[ l - 1 ] * delta[ k ] * (Tau / t);
     }

    tempHP[ comp_col[ i ] + j ] = temp;
    }
   }
  else
   if( osiSlvr->getColUpper()[ comp_col[ i ] ] > 0 )
    tempHP[ comp_col[ i ] ] -= DFi[ i ];

 // update the global lower bound- - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( HasLwrBnd )  // ... if any
  tempHP[ comp_col[ 0 ] ] += DFi[ 0 ];

 // change the bounds of the variables - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cLMRow d = Readd( true );
 for( Index i = 0 ; i < CrrSGLen ; i++ ) {
  if( Upper[ i ] < LMINF )
   Upper[ i ] -= d[ i ] * ( Tau / t );
  if( Lower[ i ] > -LMINF )
   Lower[ i ] -= d[ i ] * ( Tau / t );
  }

 // update the slack prices - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ptUpdatePricesInPlace();

 // finally, change all costs in one blow - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->setObjective( tempHP );

 // set the value of the easy components in the current point - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Tau != t )
  VectAssign( FiLambda , HpINF , NrFi );
 else
  VectAssign( FiLambda , FiLambda1 , NrFi);

 }  // end( OSIMPSolver::ChangeCurrPoint( Tau , DFi ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgSubG( cIndex strt , Index stp , cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::ChgSubG() called\n" );

 // Le modifiche devono essere fatte su una striscia
 // continua di variabili e per un insieme di colonne
 // facilmente ricavabile. Ci sono due vie
 //
 // A) nel caso lavori su una componente difficile o su
 // un numero di colonne minore di  CrrSGLen-strt :
 // Si estraggono dalla matrice le colonne in questione
 // si modificano e si reinseriscono in coda modificando
 // opportunamente tutti i dizionari
 //
 // B) altrimenti si estraggono le ultime CrrSGLen-strt
 // righe sistemandole il una matrice di appoggio che
 // viene modificata per colonne e reinserita

 throw( NDOException( "OSIMPSolver::ChgSubG not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::~OSIMPSolver()
{
 MSG( 0 , "OSIMPSolver::~OSIMPSolver() called\n" );

 cleanup();

 #if( CPLEX && OSIMPSOLVERLOG )
  CPXfclose( logfile );
 #endif

 delete derhand;
 delete osiSlvr;

 }  // end( OSIMPSolver::~OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void OSIMPSolver::cleanup( void )
{
 if( osiSlvr )
  osiSlvr->reset();

 // delete the dictionaries everything

 delete[] tempHP;
 tempHP = 0;
 tempHP_size = 0;
 delete[] tempI;
 tempI = 0;
 tempI_size = 0;

 delete[] dict_stab;
 dict_stab = 0;

 delete[] GiPerd;
 GiPerd = 0;

 delete[] FiLambda1;
 FiLambda1 = 0;
 delete[] FiLambda;
 FiLambda = 0;

 delete[] RhoCol;
 RhoCol = 0;
 delete[] RhoColBse;
 RhoColBse = 0;
 RhoColBDm = 0;

 delete[] Lower;
 Lower = 0;
 delete[] Upper;
 Upper = 0;

 delete[] wcomp;
 wcomp = 0;

 delete[] dict_slack;
 dict_slack = 0;
 delete[] dict_item;
 dict_item = 0;

 delete[] comp_col;
 comp_col = 0;
 delete[] comp_row;
 comp_row = 0;
 dict_item_maxname = 0;

 delete[] NConst;
 NConst = 0;
 delete[] NSubG;
 NSubG = 0;

 delete[] NewItem;
 NewItem = 0;

 delete[] weasy;
 weasy = 0;

 #if( PRESERVE_OSI_SOLS )
  delete[] rcst;
  rcst = 0;
  delete[] rsol;
  rsol = 0;
  delete[] csol;
  csol = 0;
 #endif

 t = 1;
 first = true;
 MaxBSize = 0;

 } // end( OSIMPSolver::cleanup )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::tUpdatePrices( cIndex strt , Index stp )
{
 if( ! dict_slack )  // no slacks
  return;            // nothing to do

 if( stp > CrrSGLen )
  stp = CrrSGLen;

 if( stab == boxstep ) {
   resizeI( 2 * ( stp - strt ) );
   resizeHP( 2 * ( stp - strt ) );
   int *ttI = tempI;
   HpRow ttHP = tempHP;

   for( Index i = strt ; i < stp ; i++ ) {
    *(ttI++) = dict_slack[ i ];
    *(ttHP++) = ( Lower[ i ] > -t ) ? Lower[ i ] : -t;
    *(ttI++) = dict_slack[ i ] + 1;
    *(ttHP++) = ( t > Upper[ i ] ) ? -Upper[ i ] : -t;
    }

   osiSlvr->setObjCoeffSet( tempI , ttI , tempHP );
   }
 else
  if( stab == quadratic ) {

   #if CPLEX == 1
    OsiCpxSolverInterface *osiCpx =
                           dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
    if( ! osiCpx )
     throw( NDOException( "OSIMPSolver::tUpdatePrices: "
    	"the OSI solver is not Cplex" ) );

    resizeHP( osiSlvr->getNumCols() );
    VectAssign( tempHP , 0.0 , osiSlvr->getNumCols() );

    CPXENVptr env = osiCpx->getEnvironmentPtr();
    CPXLPptr qp = osiCpx->getLpPtr();

    for( Index i = 0 ; i < CrrSGLen ; i++ )
     tempHP[ dict_stab[ i ] ] = -t;

    CPXcopyqpsep( env , qp , tempHP );
   #else
    throw( NDOException( "OSIMPSolver::tUpdatePrices: not implemented yet" ) );
   #endif

   }
 }  // end( OSIMPSolver::tUpdatePrices )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ptUpdatePrices( cIndex strt , Index stp )
{
 // distinguish box stabilization case because it depends on t too- - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! dict_slack )  // no slacks
  return;            // nothing to do

 if( stp > CrrSGLen )
  stp = CrrSGLen;

 resizeI( 2 * ( stp - strt ) );
 resizeHP( 2 * ( stp - strt ) );
 int *ttI = tempI;
 HpRow ttHP = tempHP;

 if( stab != boxstep )  // anything but boxstep
  for( Index i = strt ; i < stp ; i++ ) {
   int offset = 0;  // if any lower bound exists, there is a slack s_i^+
   if( Lower[ i ] > -LMINF ) {
    *(ttI++) = dict_slack[ i ];
    *(ttHP++) = Lower[ i ];
    offset = 1;
    }
   if( Upper[ i ] < LMINF ) {
    *(ttI++) = dict_slack[ i ] + offset;
    *(ttHP++) = -Upper[ i ];
    }
   }
 else                   // boxstep
  for( Index i = strt ; i < stp ; i++ ) {
   *(ttI++) = dict_slack[ i ];
   *(ttHP++) = ( Lower[ i ] > -t ) ? Lower[ i ] : -t;
   *(ttI++) = dict_slack[ i ] + 1;
   *(ttHP++) = ( t > Upper[ i ] ) ? -Upper[ i ] : -t;
   }

 osiSlvr->setObjCoeffSet( tempI , ttI , tempHP );

 }  // end( OSIMPSolver::ptUpdatePrices )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ptUpdatePricesInPlace( void )
{
 if( ! dict_slack )  // no slacks
  return;            // nothing to do

 if( stab != boxstep )  // anything but boxstep
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   int offset = 0;  // if any lower bound exists, there is a slack s_i^+
   if( Lower[ i ] > -LMINF ) {
    tempHP[ dict_slack[ i ] ] = Lower[ i ];
    offset = 1;
    }
   if( Upper[ i ] < LMINF )
    tempHP[ dict_slack[ i ] + offset ] = -Upper[ i ];
   }
 else                   // boxstep
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   tempHP[ dict_slack[ i ] ] = ( Lower[ i ] > -t ) ? Lower[ i ] : -t;
   tempHP[ dict_slack[ i ] + 1 ] = ( t > Upper[ i ] ) ? -Upper[ i ] : -t;
   }

 }  // end( OSIMPSolver::ptUpdatePricesInPlace )

/*------------------------------------------------------------------------*/

void OSIMPSolver::switchToQP( void )
{
 #if CPLEX == 1
  OsiCpxSolverInterface *osiCpx =
                           dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
  if( ! osiCpx )
   throw( NDOException( "OSIMPSolver::switchToQP: the OSI solver is not Cplex"
			) );

  CPXENVptr env = osiCpx->getEnvironmentPtr();
  CPXLPptr qp = osiCpx->getLpPtr();

  if( CPXchgprobtype( env, qp, CPXPROB_QP ) )
   throw( NDOException( "OSIMPSolver::switchToQP: can't turn to QP problem"
			) );
 #else
  throw( NDOException( "OSIMPSolver::switchToQP: not implemented yet" ) );
 #endif

 } // end( OSIMPSolver::switchToQP )

/*------------------------------------------------------------------------*/

/* inline bool OSIMPSolver::isactive( Index i )
{
 return( osiSlvr->getRowSense()[ comp_row[ NrFi ] + i ] == 'E' );
 } */

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::activate( Index i )
{
 osiSlvr->setRowType( comp_row[ NrFi ] + i , 'E' , 0.0 , 0.0 );
 }

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::deactivate( Index i )
{
 osiSlvr->setRowType( comp_row[ NrFi ] + i , 'N' , 0.0 , 0.0 );
 }

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::resizeHP( Index i )
{
 if( i > tempHP_size ) {
  delete[] tempHP;
  tempHP = new HpNum[ tempHP_size = i ];
  }
 }

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::resizeI( Index i )
{
 if( i > tempI_size ) {
  delete[] tempI;
  tempI = new int[ tempI_size = i ];
  }
 }

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::CheckDS( void )
{
 #if CHECK_DS & 1
 // obvious sanity check in dict_item[]
 if( osiSlvr ) {
  Index numcols = osiSlvr->getNumCols();
  for( Index name = 0 ; name < dict_item_maxname ; name++ )
   if( dict_item[ name ] < InINF )
    if( dict_item[ name ] > numcols )
     cout << "dict_item[ " << name << " ] = " << dict_item[ name ]
	  << " out of range (> " << numcols << ")" << endl;
  }
 #endif
 }

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::CheckBCopy( void )
{
 // check if the item is identical to any other in the bundle - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = InINF;

 if( checkID )
  if( NewItemBse ) {  // sparse format for the new item - - - - - - - - - - -
                      //- - - - - - - - - - - - - - - - - - - - - - - - - - -
   int h;
   Index NumElem;
   for( Index i = 0 ; i < dict_item_maxname ; i++ )
    if( dict_item[ i ] < InINF  )  {
     // the items can't be of different type- - - - - - - - - - - - - - - - -
     if( NewItemisSG ^ IsSubG( i ) )
      continue;

     // the subgradients must be relative to the same component - - - - - - -
     if( NewItemisSG && ( WComp( i ) != NewItemFi ) )
      continue;

     // get the item- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CoinPackedVector item =
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

     item.sortIncrIndex();
     NumElem = item.getNumElements();

     const int *B1 = item.getIndices();
     cIndex_Set B2 = NewItemBse;

     cSgRow g1 = item.getElements();
     cSgRow g2 = NewItem;

     // take off the static part of the item - - - - - - - - - - - - - - - -
     for( ; ( ( h = ( int( *B1 ) - int( comp_row[ NrFi ] ) ) ) < 0 ) &&
	    NumElem ; ) {
      B1++;
      g1++;
      NumElem--;
      }

     // checks whether or not the two items are element-wise identical - - -
     if( NumElem != NewItemBDm )
      continue;

     for( ; NumElem ; NumElem-- ) {
      if( ( *(B1++) - comp_row[ NrFi ] ) != *(B2++) )
       break;

      if( ( - *(g1++) ) != *(g2++) )
       break;
      }

     if( NumElem == 0 ) {
      IsIde = i;
      break;
      }
     }
   }
  else {  // dense format for the new item- - - - - - - - - - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   int h;
   Index l;
   Index NumElem;

   for( Index i = 0 ; i < dict_item_maxname ; i++ )
    if( dict_item[ i ] < InINF )  {
     // the items can't be of different type- - - - - - - - - - - - - - - - -
     if( NewItemisSG ^ IsSubG( i ) )
      continue;

     // the subgradients must be relative to the same component - - - - - - -
     if( NewItemisSG && ( WComp( i ) != NewItemFi ) )
      continue;

     // get the item- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CoinPackedVector item =
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

     item.sortIncrIndex();
     NumElem = item.getNumElements();
     l = 0;

     const int *B = item.getIndices();

     cSgRow g1 = NewItem;
     cSgRow g2 = item.getElements();

     // take off the static part of the item- - - - - - - - - - - - - - - - -
     for( ; ( ( h = ( int( *B ) - int( comp_row[ NrFi ] ) ) ) < 0 )
	    && NumElem ; ) {
      B++;
      g2++;
      NumElem--;
      }

     // checks whether or not the two items are element-wise identical- - - -
     for( ; NumElem ; NumElem-- , l++ ) {
      h = *(B++) - comp_row[ NrFi ];

      for( ; l < h ; l++ )
       if( *(g1++) )
	break;

      if( l < h )
       break;

      if( *(g1++) != ( -*(g2++) ) )
       break;
      }

     // the last part of g1 should be zero- - - - - - - - - - - - - - - - - -

     if( NumElem )
      continue;

     for( ; l < CrrSGLen ; l++ )
      if( *(g1++) ) {
       NumElem = InINF;
       break;
       }

     if( NumElem == 0 ) {
      IsIde = i;
      break;
      }
     }
   }

 return( IsIde );

 }  // end( OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*------------------------- End File OSIMPSolver.C -------------------------*/
/*--------------------------------------------------------------------------*/

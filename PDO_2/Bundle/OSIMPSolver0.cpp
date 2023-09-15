/*--------------------------------------------------------------------------*/
/*--------------------------- File OSIMPSolver.C ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Implementation of the OSIMPSolver class, which implements a generic  --*/
/*-- Master Problem Solver for Bundle algorithms, using a generic         --*/
/*-- OSISolverInterface object. This class conforms to the interface      --*/
/*-- defined by the class MPSolver [see MPSolver.h].                      --*/
/*--                                                                      --*/
/*--                            VERSION 1.34                              --*/
/*--                           09 - 03 - 2014                             --*/
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
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if OSIMPSOLVERLOG > 0
 #define MSG( m ) if( MPLog ) (*MPLog) << m << flush
#else
 #define MSG( m )
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

#include "OSIMPSolver.h"

#include "OPTvect.h"

#include "NDOSlver.h"

#if OSIMPSOLVERLOG > 1
 #include <strstream>
#endif

#if CPLEX == 1
 #include "OsiCpxSolverInterface.hpp"
 #include "cplex.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if OPT_USE_NAMESPACES
 using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::OSIMPSolver( istream *iStrm )
{
 // reset the pointers to OsiSolver and FiOracle object  - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 osiSlvr = 0;
 FIO = 0;

 // set the default algorithm parameters - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 useactiveset = false;  // no active set
 first = true;          // no problem has been solved so far
 checkID = false;       // setta default checkID

 stab = unset;	        // the stabilization term is unknown
 t = 1;			// the proximity parameter is set to a standard value

 // the master will be empty ...
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MaxBSize = 0;		// max bundle is zero
 dict_item = 0;	        // no items are in the master problem
 dict_item_maxname = 0;
 dict_slack = 0;
 dict_stab = 0;
 newly_inserted = 0;	// after last call to either CheckLinErr or
                        // ChangeCurrPoint, no items have been added
 just_inserted = 0;	// after the last call to SolveMP, no items have
                        // been added
 wcomp = 0;

 NewItem = 0;	        // new item information is empty
 NewItemBse = 0;
 NewItemBDm = 0;

 Upper = Lower = 0;
 FiLambda = FiLambda1 = 0;

 // the temporary buffers of HpNum and Index are used to save the information
 // before passing it to the OSISolver; by default they are empty
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 tempHP = 0;
 tempHP_size = 0;
 tempI = 0;
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

 }  // end( OSIMPSolver::OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetDim( cIndex MxBSz , FiOracle *Oracle ,
			  const bool UsAvSt )
{
 MSG( "OSIMPSolver::SetDim()\n" );

 if( ! MxBSz ) {  // deallocate all its memory and quietly wait
  cleanup();      // for new instructions- - - - - - - - - - - - - - - - - - -
  return;
  }

 if( Oracle ) {  // allocate the memory for solving the problem- - - - - - - -
  if( ! osiSlvr )
   throw( NDOException( "OSIMPSolver::SetDim: osiSlvr must be set" ) );

  // discards all the previous settings and deallocate the whole memory  - - -
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

  NewItem = new SgNum[ MaxNZ ];

  comp_row = new Index[ NrFi + 1 ];
  comp_col = new Index[ NrFi + 1 ];
  dict_item = new Index[ MaxBSize ];
  dict_slack = new Index[ MaxSGLen ];
  wcomp = new Index[ MaxBSize ];
  newly_inserted = new bool[ MaxBSize ];
  just_inserted = new bool[ MaxBSize ];

  Upper = new LMNum[ MaxSGLen ];
  Lower = new LMNum[ MaxSGLen ];

  // no item has just been added  - - - - - -  - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VectAssign( newly_inserted , false , MaxBSize );
  VectAssign( just_inserted , false , MaxBSize );

  // dictionaries are empty  - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VectAssign( dict_item , Index( Inf<Index>() ) , MaxBSize );
  VectAssign( dict_slack , Index( Inf<Index>() ) , MaxSGLen );

  VectAssign( wcomp , Index( Inf<Index>() ) , MaxBSize );
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

  for( Index i = 1 ; i <= NrFi ; i++ ) {
   comp_row[ i - 1 ] = osiSlvr->getNumRows();
   Index BNC = FIO->GetBNC( i );
   if( BNC ) {
    MSG( "Structured constraints for the easy component " << i
         << " have been created starting from the row "  << comp_row[ i - 1 ]
	 << endl );

    for( Index j =  2 * (BNC + FIO->GetBNR( i ) ) ; j > 0 ; j-- )
     osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0.0 );
    }
   else {
    MSG( "Structured constraints for the difficult component  " << i
	  << " have been created into the row " << comp_row[ i - 1 ] << endl
	 );
    osiSlvr->addRow( 0 , 0 , 0 , 0.0 , 0.0 );
    }
   }

  // if the initial dimension of the items is greater than zero, generate the
  // (partial) rows of the the dynamic part which describes the form of the
  // dual variables z  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_row[ NrFi ] = osiSlvr->getNumRows();
  if( useactiveset ) { // all the constraints are inactive so far  - - - - - -
    MSG( "The dynamic inactive constraints have been created starting "
         "from the row " << comp_row[ NrFi ] << endl );

    for( Index j = CrrSGLen ; j > 0 ; j-- )
     osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() ,
		                    osiSlvr->getInfinity() );
   }
  else {     // all the constraints are active (now and for ever)  - - - - - -
    MSG( "The dynamic constraints have been created starting from the row "
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
  VectAssign( FiLambda , HpNum( Inf<HpNum>() ), NrFi );  // Fi( Lambda ) == ?

  // allocate the memory for Gi * d
  GiPerd = new HpNum[ MaxBSize ];

  // recover the static information of the easy component  - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for ( Index i = 1 ; i <= NrFi ; i++ ) {
   Index BNC , BNR , BNZ , ANZ;
   if( ( BNC = FIO->GetBNC( i ) ) ) { // FiOracle returns all the information
	    // in a sparse format, so the columns of the matrices A and B are
	    // in sparse format- - - - - - - - - - - - - - - - - - - - - - - -

    BNR = FIO->GetBNR( i );  // number of rows of matrix B[ i ]
    BNZ = FIO->GetBNZ( i );  // number of non-zeroes of matrix B[ i ]
    ANZ = FIO->GetANZ( i );  // number of non-zeroes of matrix A[ i ]

    // allocate the memory to describe the matrix B[ i ]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    int* Bbeg = new int[ BNC + 1 ];
    int* Bind = new int[ BNZ ];
    double *Bval = new double[ BNZ ];
    double *lhs = new double[ BNR ];
    double *rhs = new double[ BNR ];
    double *cst = new double[ BNC ];
    double *lbd = new double[ BNC ];
    double *ubd = new double[ BNC ];

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

    // write the rows relative to the easy components  - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    for ( Index j = 0 ; j < BNC ; j++ ) {  // for each variable j of the
     Index count = 0;                      // i-th easy component do ...

     // write the static part: B's columns and the box constraints
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     for( int k = Bbeg[ j ] ; k < Bbeg[ j + 1 ] ; k++ ) {
      tempI[ count ] = Bind[ k ] + comp_row[ i - 1 ];
      tempHP[ count++ ] = Bval[ k ];   // B[ i ] x[ i ] - e[ i ] <= 0
      tempI[ count ] = Bind[ k ] + comp_row[ i - 1 ] + BNR;
      tempHP[ count++ ] = -Bval[ k ];  // -B[ i ] x[ i ] + d[ i ] <= 0
      }

     tempI[ count ] = comp_row[ i - 1 ] + ( 2 * BNR ) + j;
     tempHP[ count++ ] = 1.0;          // x[ i ] - u[ i ] <= 0
     tempI[ count ] = comp_row[ i - 1 ] + ( 2 * BNR ) + BNC + j;
     tempHP[ count++ ] = -1.0;         // -x[ i ] + l[ i ] <= 0

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

     MSG( "The variables "<< j << "of the easy component "<< i
          << "has been added" << endl );

     }  // end adding of j-th column of the component i- - - - - - - - - - - -
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // put in RhoCol the vectors: e[ i ] , d[ i ] , u[ i ] and l[ i ]
    // [see above]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    for( Index j = 0 ; j < BNR ; j++ ) {
     RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + j;
     RhoCol[ RhoColBDm++ ] = - rhs[ j ];
     RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + BNR + j;
     RhoCol[ RhoColBDm++ ] = lhs[ j ];
     }

    for( Index j = 0 ; j < BNC ; j++ ) {
     RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + 2 * BNR + j;
     RhoCol[ RhoColBDm++ ] = - ubd[ j ];
     RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ] + 2 * BNR + BNC + j;
     RhoCol[ RhoColBDm++ ] = lbd[ j ];
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
    }
   }  // end easy components part description- - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // add the description of the difficult components - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 1 ; i <= NrFi ; i++ ) {
   if( FIO->GetBNC( i ) == 0 ) { // for each component "i" we must create a
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
    MSG( "The variable gamma " << i << " has been added and its name is "
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
   VectAssign( dict_stab , Index( Inf<Index>() ) , MaxSGLen );
   }
  else
   dict_stab = 0;

  for( Index i = 0 ; i < CrrSGLen ; i++) { // ask FiOracle if there are some
                               // lower and(/or) upper bounds on the primal
   bool slack_p = false;       // variables  - - - - - - - - - - - - - - - - -
   bool slack_m = false;

   if( ( Upper[ i ] = Oracle->GetUB( i ) ) < Inf<LMNum>() )
    slack_m = true;  // there is some upper bound on variable i
   if( Oracle->GetUC( i ) )
    Lower[ i ] = -Inf<LMNum>();
   else {            // the primal variable is constrained to be nonnegative
    Lower[ i ] = 0;
    slack_p = true;
    }

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
     MSG( "The stabilization variable z_"<< i
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
    MSG( "The slack s_"<< i <<"+ has been created and its name is"
    	<< dict_slack[ i ] << endl );
    }

   if( slack_m ) {
    tempI[ 0 ] = comp_row[ NrFi ] + i;
    tempHP[ 0 ] = -1.0;
    osiSlvr->addCol( 1 , tempI , tempHP , 0.0 , osiSlvr->getInfinity() ,
		     0.0 );
    if( ! slack_p )
     dict_slack[ i ] = osiSlvr->getNumCols() - 1;
     MSG( "The slack s_"<< i << "- has been created and its name is "
      << slack_p? dict_slack[ i ] + 1 : dict_slack[ i ] << endl );
     }

   }  // end stabilization and slack definition- - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // add the column of rho in the master - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_col[ 0 ] = osiSlvr->getNumCols();
  osiSlvr->addCol( RhoColBDm , RhoColBse , RhoCol , 0.0 ,
		   osiSlvr->getInfinity() , 0.0 );
  MSG( "The column of rho has been created and its name is " << comp_col[ 0 ]
       << endl );

  // change the price of the slack variables - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tUpdatePrices();       // change the prices for the stabilization

  if( stab != boxstep )  // change the slack prices
   ptUpdatePrices();

  }  // end( if(Oracle) )
 else {  // sets the max bundle size to n and activate/deactivate the Active
         // Set Mechanism without changing anything else - - - - - - - - - - -

  // deletes the item in excess  - - - - - - - - - - - - - - - - - - - - - - -

  for( ; dict_item_maxname-- > MxBSz ; )
   if( dict_item[ dict_item_maxname ] != Inf<Index>() )
    RmvItem( dict_item_maxname );

  if( MxBSz != MaxBSize ) {  // the existing items in the bundle (if any) are
   // all kept if MxBSz is larger than MaxBSize, but a smaller value will
   // force deletion of all the items with "name"  >= MxBSz

   Index_Set olddict_item = dict_item;
   Index_Set oldwcomp = wcomp;
   Bool_Vec oldnewly_inserted = newly_inserted;
   Bool_Vec oldjust_inserted = just_inserted;

   dict_item = new Index[ MxBSz ];
   wcomp = new Index[ MxBSz ];
   newly_inserted = new bool[ MxBSz ];
   just_inserted = new bool[ MxBSz ];

   VectAssign( newly_inserted , false , MxBSz );
   VectAssign( just_inserted , false , MxBSz );
   VectAssign( dict_item , Index( Inf<Index>() ) , MxBSz );
   VectAssign( wcomp , Index( Inf<Index>() ) , MxBSz );

   Index copydim = ( MxBSz < MaxBSize ) ? MxBSz : MaxBSize;

   if( olddict_item ) {
    VectAssign( dict_item , olddict_item , copydim );
    delete[] olddict_item;
    }

   if( oldwcomp ) {
    VectAssign( wcomp , oldwcomp , copydim );
    delete[] oldwcomp;
    }

   if( oldnewly_inserted ) {
    VectAssign( newly_inserted , oldnewly_inserted , copydim );
    delete[] oldnewly_inserted;
    }

   if( oldjust_inserted ) {
    VectAssign( just_inserted , oldjust_inserted , copydim );
    delete[] oldjust_inserted;
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
  MSG( "OSIMPSolver::Sett(): t = " << t << endl );
  tUpdatePrices();
  }
 }  // end( OSIMPSolver::Sett )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetOptPrcsn( HpNum OEps )
{
 MSG( "OSIMPSolver::SetOptPrcsn()\n" );
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
 MSG( "OSIMPSolver::SetFsbPrcsn()\n" );
 FsbEps = FEps;
 osiSlvr->setDblParam( OsiPrimalTolerance , double( FEps ) );

 }  // end( OSIMPSolver::SSetFsbPrcsn )

/*------------------------------------------------------------------------*/

void OSIMPSolver::SetLowerBound( cHpNum LwBnd , cIndex wFi )
{
 MSG( "OSIMPSolver::SetLowerBound()\n" );

 if( ( wFi < Inf<Index>() ) && iseasy( wFi ) )
  throw( NDOException(
    "OSIMPSolver::SetLowerBound: Setting LowerBound for a easy component" ) );

 if( LwBnd > -Inf<HpNum>() ) {
  if( wFi < Inf<Index>() ) {  // insert the constraint gamma_i <= 1 and
	                      // add LwBnd * gamma_i to the objective function
   osiSlvr->setColUpper( comp_col[ wFi ] , 1.0 );
   osiSlvr->setObjCoeff( comp_col[ wFi ] , LwBnd );
   MSG( "LB( "<< wFi << " ) = " << LwBnd << endl );
   }
  else {   // set the total lower bound- - - - - - - - - - - - - - - - - - - -
   // Note that it needs change the constraint of rho:
   // r + rho = 1 ----> r = 1 - rho => 0 [ rho <= 1 ]
   // and in the objective function LwBnd * r = LwBnd - LwBnd * rho
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   HasLwrBnd = true;
   osiSlvr->setRowBounds( 0 , - osiSlvr->getInfinity() , 1.0 );
   osiSlvr->setObjCoeff( comp_col[ 0 ] , - LwBnd );
   MSG( "Total LB = " << LwBnd << endl );
   }
  }
 else {  // if LwBnd = -Inf<HpNum>() reset the previous lower bound ...
  if( wFi < Inf<Index>() ) { // if some items exist set gamma_i = 0
   MSG( "Lower bound of "<< wFi << " has been deleted " << endl );
   osiSlvr->setObjCoeff( comp_col[ wFi ] , 0.0 );
   if( dict_item_maxname )
    osiSlvr->setColUpper( comp_col[ wFi ] , 0.0 );
   }
  else {  // change the rho constraint --> rho = 1 and put the coefficient
	  // of rho equal to zero in the objective function

   HasLwrBnd = false;
   MSG( "Total lower bound is - Inf<double>()\n" );
   osiSlvr->setObjCoeff( comp_col[ 0 ] , 0.0 );
   osiSlvr->setRowBounds( 0 ,  1.0 , 1.0 );
   }
  }
 }  // end( OSIMPSolver::SetLowerBound )

/*------------------------------------------------------------------------*/

void OSIMPSolver::SetMPLog( ostream *outs , const char lvl )
{
 MPSolver::SetMPLog( outs , lvl );

 if( osiSlvr )
  osiSlvr->messageHandler()->setLogLevel( int( MPLLvl ) );

 #if CPLEX
  OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>(
								 osiSlvr );
  if( ! osiCpx )
   throw( NDOException( "OSIMPSolver::SetAlgo: the OSI solver is not Cplex"
			) );

  CPXENVptr env = osiCpx->getEnvironmentPtr();

  CPXsetintparam( env , CPX_PARAM_MIPDISPLAY , MPLLvl );

  CPXFILEptr logfile = CPXfopen( "cplex.log" , "w" );
  int status = CPXsetlogfile( env , logfile );
  assert( ! status );
 #endif

 } // end( OSIMPSolver::SetMPLog )

/*------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------*/

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

 #if CPLEX == 1
  OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>( osi );
  if( ! osiCpx )
   throw( NDOException(
                   "OSIMPSolver::SetosiSlvr: the OSI solver is not Cplex" ) );
 #endif

 if( ! ( osiSlvr = osi ) )
  throw( NDOException( "OSIMPSolver::SetosiSlvr: wrong OSI solver" ) );

 osiSlvr->messageHandler()->setLogLevel( int( MPLLvl ) );

 } // end( OSIMPSolver::SetOsi )

/*------------------------------------------------------------------------*/

void OSIMPSolver::SetStabType( const StabFun sf )
{
 MSG( "OSIMPSolver::SetStabType(()\n" );
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

/*------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*------------------------------------------------------------------------*/

OSIMPSolver::MPStatus OSIMPSolver::SolveMP( void )
{
 MSG( "OSIMPSolver::SolveMP()\n" );

 if( MPt )
  MPt->Start();

 // get the dimension of the Master Problem  - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int nc = osiSlvr->getNumCols();
 int nr = osiSlvr->getNumRows();

 #if OSIMPSOLVERLOG > 1
  if( osiSlvr && nc && nr ) {
   std::strstream filename;
   filename << "MP_" << FIO->GetNDOSolver()->NrIter() << ".mps" << std::ends;
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
  MSG( "MP is dual infeasible" << endl );
  return( kUnfsbl );
  }

 if( osiSlvr->isProvenPrimalInfeasible() ) {
  MSG( "MP is primal infeasible" << endl );
  return( kUnbndd );
  }

 if( ! osiSlvr->isProvenOptimal() )
  if( osiSlvr->isAbandoned() )
   MSG( "Warning: numerical difficulties in the solver, ignoring them"
	<< endl );
  else {
   MSG( "Some unknown error happened" << endl );
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
  if( wcomp[ i ] <= NrFi ) {
   GiPerd[ i ] = rcst[ dict_item[ i ] ] - objcoeff[ dict_item[ i ] ];
   if( IsSubG( i ) )
    GiPerd[ i ] += rsol[ comp_row[ wcomp[ i ] - 1 ] ];
   }

 // compute the actual function value of the easy components by using the row
 // price: Fi[ wFi ]( Lambda + d* ) = \pi_e[ wFi ] * e[ wFi ] - \pi_d[ wFi ] *
 // d[ wFi ] + \pi_u[ wFi ] * u[ wFi ] - pi_l[ wFi ] * l[ wFi ]
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( FiLambda1 , HpNum( 0 ) , NrFi );
 for( Index wFi = 1; wFi <= NrFi; wFi++ )
  if( wFi <= NrFi && iseasy( wFi ) ) {
   for( Index i = 0 ; i < RhoColBDm ;  i++)
    if( ( Index( RhoColBse[ i ] ) >= comp_row[ wFi - 1 ] ) &&
	( Index( RhoColBse[ i ] ) < comp_row[ wFi ] ) )
     FiLambda1[ wFi - 1 ] -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];

   MSG( "FiBLambda(" << wFi << ") = " << FiLambda1[ wFi - 1 ] << endl );
   }

 // print the primal/dual solution  - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( OSIMPSOLVERLOG > 2 )
  if( MPLog ) {
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

 VectAssign( just_inserted , false , dict_item_maxname );

 if( MPt )
  MPt->Stop();

 MSG( "MP solved with success" << endl );
 return( kOK );

 }  // end( OSIMPSolver::SolveMP )

/*------------------------------------------------------------------------*/
/*--------------------- METHODS FOR READING RESULTS ----------------------*/
/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadFiBLambda( cIndex wFi )
{
 MSG( "OSIMPSolver::ReadFiBLambda()\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif

 HpNum FiVal = 0;
 if( ! wFi ) { // compute the function value of the 0-th component - - - - -
  // it is given by Fi[ 0 ]( d^* ) = d^* b. Note that d^* is
  // the dual solution of the dynamic part and it  represents the search
  // direction - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = comp_row[ NrFi ] ; i < RhoColBDm  ; i++ )
   FiVal -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];
  MSG( "FiBLambda( 0 ) = " << total << endl );
  }
 else
  if( wFi <= NrFi ) { // Note that the difficult easy case is quite
	  // different from the difficult one. In fact, in the easy case
	  // it computes the actual value of the function while in the other case
	  // just an approximation is provided
   if( iseasy( wFi ) )
	FiVal = FiLambda1[ wFi -1 ];
   else
	FiVal = rsol[ comp_row[ wFi - 1 ] ];

   MSG( "FiBLambda(" << wFi << ") = " << total << endl );
   }
  else { // return all the components of the function  - - - - - - - - - - -
         //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index i = 1 ; i <= NrFi ; i++ )
	if( iseasy( i ) )
     FiVal += FiLambda1[ i -1 ];
	else
     FiVal += rsol[ comp_row[ i - 1 ] ];

   if( wFi == Inf<Index>() )   // add the linear part  - - - - - - - - - - -
    for( Index i = comp_row[ NrFi ] ; i < RhoColBDm ; i++ )
     FiVal -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];

   MSG( "FiBLambda( INF ) = " << total << endl );

   }  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( FiVal );

 }  // end( OSIMPSolver::ReadFiBLambda )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadDt( cHpNum tt )
{
 MSG( "OSIMPSolver::ReadDt()\n" );

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
     MSG( "Dt = +Inf<double>()" << endl );
     return( Inf<HpNum>() );
     }
    }
   break;
  case quadratic:
   for( Index n = 0 ; n < CrrSGLen ; n++ )
	value += csol[ dict_stab[ n ] ] * csol[ dict_stab[ n ] ];
   value *= ( ( tt == t ? tt : ( t * t ) / tt ) / 2 );
   MSG( "Dt = " << value << endl );
   break;
  default:
   throw( NDOException( "OSIMPSolver::ReadDt: undecided stabilization" ) );
   }

 MSG( "Dt = " << value << endl );
 return( value );

 }  // end( OSIMPSolver::ReadDt )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadSigma( cIndex wFi )
{
 MSG( "OSIMPSolver::ReadSigma()\n" );

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
   if( ( wcomp[ i ] <= NrFi ) && ( ! just_inserted[ i ] ) )
    value -= csol[ dict_item[ i ] ] * obj[ dict_item[ i ] ];

  // add the easy components cost part and the lower bounds contribution of
  // the difficult components  - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 1 ; i <= NrFi ; i++ )
   if( iseasy( i ) ) {
    if( FiLambda[ i - 1] < Inf<HpNum>() ) {
     for( Index j = 0; j < FIO->GetBNC( i ) ; j++ )
      value -= csol[ comp_col[ i ] + j ] * obj[ comp_col[ i ] + j ];
     value += FiLambda[ i - 1];
     }
    else
     value = Inf<HpNum>();
    }
   else
    value -= csol[ comp_col[ i ] ] * obj[ comp_col[ i ] ];

  // divide per rho  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  // finally add the contribution of the slack variables
  // - - - - - -  - - - - - - -  - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0; i < CrrSGLen ; i++ ) {
   // if -t < Lower[ i ], the stabilization constraint -t <= d[ i ] is
   // redundant, and the opposite for Upper[ i ]
   if( -t < Lower[ i ] )
    value -= Lower[ i ] * csol[ dict_slack[ i ] ];
   if( t > Upper[ i ] )
    value += Upper[ i ] * csol[ dict_slack[ i ] + 1 ];
   }

  MSG( "Sigma* + Sigma_L = " << value << endl );
  }
 else {
  if( iseasy( wFi ) ) {
   if( FiLambda[ wFi - 1] < Inf<HpNum>() ) {
    for( Index j = 0; j < FIO->GetBNC( wFi ) ; j++ )
     value -= csol[ comp_col[ wFi ] + j ] * obj[ comp_col[ wFi ] + j ];
    value += FiLambda[ wFi - 1];
    }
   else
    value = Inf<HpNum>();
  }
  else {
   for( Index i = 0 ; i < dict_item_maxname ; i++ )
    if( ( wcomp[ i ] == wFi ) && ( ! just_inserted[ i ] ) )
     value -= csol[ dict_item[ i ] ] * obj[ dict_item[ i ] ];
   value -= csol[ comp_col[ wFi ] ] * obj[ comp_col[ wFi ] ];
   }

  if( csol[ comp_col[ 0 ] ] != 0 )
   value /= csol[ comp_col[ 0 ] ];
  MSG( "Sigma*(" << wFi << ") = " << value << endl );
  }

 return( value );

 }  // end( OSIMPSolver::ReadSigma )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadDStart( cHpNum tt )
{
 MSG( "OSIMPSolver::ReadDStart()" << endl );

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

 MSG( "D* = " << total << endl );
 return( total );

 }  // end( OSIMPSolver::ReadDStart )

/*------------------------------------------------------------------------*/

cLMRow OSIMPSolver::Readd( bool Fulld )
{
 MSG( "OSIMPSolver::Readd()\n" );

 if( Fulld ) {
  #if( PRESERVE_OSI_SOLS == 0 )
   const double *rsol = osiSlvr->getRowPrice();
  #endif

  return( rsol + comp_row[ NrFi ] );
  }
 else
  throw( NDOException( "OSIMPSolver::Readd: Fulld( false )" ) );

 }  // end( OSIMPSolver::Readd )

/*------------------------------------------------------------------------*/

void OSIMPSolver::ReadZ( LMRow tz , cIndex_Set &I , Index &D , cIndex wFi )
{
 MSG("OSIMPSolver::ReadZ() \n");

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
    if( ( dict_item[ i ] < Inf<Index>() ) && ( IsSubG( i ) ) &&
	( ! just_inserted[ i ] ) )
     temp = temp + csol[ dict_item[ i ] ] *
                  ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

   // add the sugradient of the easy components - - - - - - - - - - - - - - -

   for( Index i = 1 ; i <= NrFi ; i++ )
    if( iseasy( i ) )
     for( int j = comp_col[ i ] ; j < comp_col[ i ] + FIO->GetBNC( i ) ; j++ )
      temp = temp + csol[ j ] * ( osiSlvr->getMatrixByCol() )->getVector( j );

   if( wFi == Inf<Index>() )   // add the linear part   - - - - - - - - - - -
    temp = temp + ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
   }
  else { // disaggregated case  - - - - - - - - - - - - - - - - - - - - - - -
   if( iseasy( wFi ) )
    throw( NDOException( "OSIMPSolver::ReadZ: calling for easy component" ) );

   for( Index i = 0 ; i < dict_item_maxname ; ++i )
    if( ( wcomp[ i ] == wFi ) && ( IsSubG( i ) ) && ( ! just_inserted[ i ] ) )
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

/*------------------------------------------------------------------------*/

cHpRow OSIMPSolver::ReadMult( cIndex_Set &I , Index &D , cIndex wFi ,
			      const bool IncldCnst )
{
 MSG( "OSIMPSolver::ReadMult()\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ReadMult for 0-th component" ) );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 HpNum rho = csol[ comp_col[ 0 ] ];

 if( rho == 0 ) {  // in the case rho = 0, return a null vector- - - - - - -
  resizeI( 1 );
  tempI[ D = 0 ] = Inf<Index>();
  I = reinterpret_cast<cIndex_Set>( tempI );
  return( tempHP );
  }

 if( wFi <= NrFi )
  if( iseasy( wFi ) ) {  // easy component part- - - - - - - - - - - - - - -
   I = 0;
   D = FIO->GetBNC( wFi );
   resizeHP( D );
   VectAssign( tempHP , csol + comp_col[ wFi ] , 1 / rho , D );
   }
  else {                 // difficult component part - - - - - - - - - - - -
   resizeHP( MaxBSize );
   resizeI( MaxBSize + 1 );
   I = reinterpret_cast<cIndex_Set>( tempI );
   D = 0;
   for( Index name = 0 ; name < dict_item_maxname ; name++ ) {
    if( ( dict_item[ name ] < Inf<Index>() ) && ( wcomp[ name ] == wFi ) &&
	( IncldCnst || IsSubG( name ) ) )
     if( ( ! just_inserted[ name ] ) &&
	 ( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 ) )
      tempI[ D++ ] = name;
    }
   tempI[ D ] = Inf<Index>();
   }
 else {  // wFi == Inf<Index>()- - - - - - - - - - - - - - - - - - - - - - -
  resizeHP( MaxBSize );
  resizeI( MaxBSize + 1 );
  I = reinterpret_cast<cIndex_Set>( tempI );
  D = 0;
  for( Index name = 0 ; name < dict_item_maxname ; name++ ) {
   if( ( dict_item[ name ] < Inf<Index>() ) &&
       ( IncldCnst || IsSubG( name ) ) )
    if( ( ! just_inserted[ name ] ) &&
	( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
     tempI[ D++ ] = name;
   }

  tempI[ D ] = Inf<Index>();
  }

 return( tempHP );

 }  // end( OSIMPSolver::ReadMult )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadLBMult( cIndex wFi )
{
 MSG( "OSIMPSolver::ReadLBMult()\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ReadLBMult for 0-th component" ) );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 HpNum value;
 if( wFi <= NrFi ) {  // return gamma_i
  if( iseasy( wFi ) )
   throw( NDOException( "OSIMPSolver::ReadLBMult for an easy component" ) );

  value = csol[ comp_col[ wFi ] ];
  }
 else  // return r = ( 1 - r ) - - - - - - - - - - - - - - - - - - - - - - -
  value = 1 - csol[ comp_col[ 0 ] ];

 return( value );
 }

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadGid( cIndex Nm )
{
 MSG( "OSIMPSolver::ReadGid() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif

 HpNum value = 0;
 if( Nm >= MaxBSize ) {  // zero component: return delta * b
  for( Index i = comp_row[ NrFi ] ; i < RhoColBDm ; i++ )
   value -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];
  MSG( "b * delta = " << value << endl );
  }
 else {
  if( dict_item[ Nm ] == Inf<Index>() )
   throw( NDOException( "OSIMPSolver::ReadGid: unused item name" ) );

  if( just_inserted[ Nm ] ) {  // compute the product directly, since it's
                               // not possible to use the reduced cost

   CoinShallowPackedVector item =
                ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ Nm ] );
   const int *j = item.getIndices();
   for( int i = item.getNumElements() ; i-- > 0 ; ) {
    const Index js = static_cast<Index>( *(j++) );
    if( js >= comp_row[ NrFi ] )
     value -= item[ js ] * rsol[ js ];
    }
   }
  else
   value = GiPerd[ Nm ];

  MSG( "Item[" << Nm << "] * delta = " << value << endl );
  }

 return( value );

 }  // end( OSIMPSolver::ReadGid )

/*------------------------------------------------------------------------*/

void OSIMPSolver::MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau )
{
 MSG( "OSIMPSolver::MakeLambda1()\n" );

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

/*------------------------------------------------------------------------*/

void OSIMPSolver::SensitAnals( HpNum &lp , HpNum &cp )
{
 throw( NDOException( "OSIMPSolver::SensitAnals not implemented yet" ) );
 }

/*------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE PROBLEM --------------*/
/*------------------------------------------------------------------------*/

Index OSIMPSolver::BSize( cIndex wFi )
{
 MSG( "OSIMPSolver::BSize()\n" );

 Index count = 0;
 if( wFi == Inf<Index>() ) {
  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( dict_item[ i ] < Inf<Index>() )
    count++;
  }
 else
  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( wcomp[ i ] == wFi )
    count++;

 return( count );

 }  // end( OSIMPSolver::BSize )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::BCSize(cIndex wFi)
{
 MSG( "OSIMPSolver::BCSize()\n" );

 Index count = 0;
 if( wFi == Inf<Index>() ) {
  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( ( dict_item[ i ] < Inf<Index>() ) && ( ! IsSubG( i ) ) )
    count++;
  }
 else
  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( ( wcomp[ i ] == wFi ) && ( ! IsSubG( i ) ) )
    count++;

 return( count );

 }  // end( OSIMPSolver::BCSize )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::MaxName( cIndex wFi )
{
 MSG( "OSIMPSolver::MaxName()\n" );

 Index i = dict_item_maxname;
 if( wFi != Inf<Index>() )
  for( ; i > 0 ; i-- )
   if( wcomp[ i - 1 ] == wFi )
    break;

 return( i );

 }  // end( OSIMPSolver::MaxName )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::WComponent( cIndex i )
{
 MSG( "OSIMPSolver::WComponent()\n" );

 return( wcomp[ i ] );

 }  // end( OSIMPSolver::WComponent )

/*------------------------------------------------------------------------*/

bool OSIMPSolver::IsSubG( cIndex i )
{
 MSG( "OSIMPSolver::IsSubG()\n" );

 if( dict_item[ i ] != Inf<Index>() ) {
  const CoinShallowPackedVector v =
   osiSlvr->getMatrixByCol()->getVector( dict_item[ i ] );
  if( v[ comp_row[ wcomp[ i ] - 1 ] ] == 1 )
    return( true );
  }

 return( false );

 }  // end( OSIMPSolver::IsSubG )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::NumNNVars( void )
{
 MSG( "OSIMPSolver::NumNNVars()\n" );

 Index count = 0;
 for( Index i = 0 ; i < CrrSGLen ; i++ )
  if( Lower[ i ] > -Inf<LMNum>() )
   count++;

 return( count );

 } // end( OSIMPSolver::NumNNVars )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::NumBxdVars( void )
{
 MSG( "OSIMPSolver::NumBxdVars()\n" );

 Index count = 0;
 for( Index i = 0 ; i < CrrSGLen ; i++ )
  if( ( Lower[ i ] > -Inf<LMNum>() ) || ( Upper[ i ] < Inf<LMNum>() ) )
   count++;

 return( count );

 } // end( OSIMPSolver::NumBxdVars )

/*------------------------------------------------------------------------*/

bool OSIMPSolver::IsNN( cIndex i )
{
 MSG( "OSIMPSolver::IsNN()\n" );

 return( Lower[ i ] > -Inf<LMNum>() );

 } // end( OSIMPSolver::IsNN )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::CheckLinErr( cHpNum AEps , const bool All , cIndex wFi )
{
 // Checks the existence of linearization errors that be strictly less than
 // - AEps  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MSG( "OSIMPSolver::CheckLinErr()\n" );

 const double *linerr = osiSlvr->getObjCoefficients();
 resizeI( dict_item_maxname );

 // look for the smallest linearization error   - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double smallest = 0;
 for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ ) {
  if( ( dict_item[ Nm ] < Inf<Index>() ) && 
      ( wcomp[ Nm ] == wFi || wFi == Inf<Index>() ) &&
      ( newly_inserted[ Nm ] || All ) ) {
   if( ( linerr[ dict_item[ Nm ] ] > AEps ) &&
       ( - linerr[ dict_item[ Nm ] ] < smallest ) )
    smallest = - linerr[ dict_item[ Nm ] ];
   newly_inserted[ Nm ] = false;
   }
  }

 // change all the alphas  - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( smallest < 0 )
  for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
   if( ( dict_item[ Nm ] < Inf<Index>() ) &&
       ( wcomp[ Nm ] == wFi || wFi == Inf<Index>() ) )
    osiSlvr->setObjCoeff( dict_item[ Nm ] ,
			  linerr[ dict_item[ Nm ] ] + smallest );
 return( smallest );

 }  // end( OSIMPSolver::CheckLinErr )

/*------------------------------------------------------------------------*/

void OSIMPSolver::CheckIdentical( const bool Chk )
{
 MSG( "OSIMPSolver::CheckIdentical()\n" );
 checkID = Chk;

 }  // end( OSIMPSolver::CheckIdentical )

/*------------------------------------------------------------------------*/

cHpRow OSIMPSolver::ReadLinErr( void )
{
 MSG( "OSIMPSolver::ReadLinErr()\n" );
 resizeHP( MaxBSize );
 VectAssign( tempHP , 0.0 , MaxBSize );

 const double *linerr = osiSlvr->getObjCoefficients();
 for( Index i = 0 ; i < MaxBSize ; i++ )
  if( dict_item[ i ] < Inf<Index>() )
   tempHP[ i ] = - linerr[ dict_item[ i ] ];  // change sign!

 return( tempHP );

 }  // end( OSIMPSolver::ReadLinErr )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::EpsilonD( void )
{
 MSG( "OSIMPSolver::EpsilonD()\n" );

 return( FsbEps );

 } // end( OSIMPSolver::EpsilonD )

/*------------------------------------------------------------------------*/

bool OSIMPSolver::FiBLambdaIsExact( cIndex wFi )
{
 MSG( "OSIMPSolver::FiBLambdaIsExact()\n" );

 return( iseasy( wFi ) );

 } // end( OSIMPSolver::FiBLambdaIsExact )

/*------------------------------------------------------------------------*/
/*------------ METHODS FOR ADDING / REMOVING / CHANGING DATA -------------*/
/*------------------------------------------------------------------------*/

SgRow OSIMPSolver::GetItem( cIndex wFi )
{
 MSG( "OSIMPSolver::GetItem()\n" );
 NewItemFi = wFi;  // mark the name of the component

 if( NewItem )
  return( NewItem );
 else
  throw( NDOException( "OSIMPSolver::GetItem allocation error" ) );

 } // end( OSIMPSolver::GetItem )

/*------------------------------------------------------------------------*/

void OSIMPSolver::SetItemBse( cIndex_Set SGBse , cIndex SGBDm )
{
 MSG( "OSIMPSolver::SetItemBse()\n" );
 if( ( NewItemBse = SGBse ) )
  NewItemBDm = SGBDm;    // sparse format
 else
  NewItemBDm = CrrSGLen; // dense format

 } // end( OSIMPSolver::SetItemBse )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai ,
			      HpNum &ScPri )
{
 MSG( "OSIMPSolver::CheckSubG() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 // compute the scalar product Gi^{top} delta  - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NewItemBse )  // sparse format
  ScPri = ScalarProduct( NewItem , primal , NewItemBse );
 else              // dense format
  ScPri = ScalarProduct( NewItem , primal , CrrSGLen );

 // compute the linearization error of the new subgradient
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Ai = Ai - ( DFi - ( Tau / t ) * ScPri );

 NewItemprice = Ai;
 NewItemScPri = ScPri;
 NewItemisSG = true;   // this item is a subgradient ...

 // Note!!!! It doesn't check whether the item is or not identical to some
 // other items

 return( Inf<Index>() );

 } // end( OSIMPSolver::CheckSubG )

/*------------------------------------------------------------------------*/

Index OSIMPSolver::CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt )
{
 MSG( "OSIMPSolver::CheckCnst() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 // compute ScPri and right hand side of the constraint- - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NewItemBse ) {  // sparse format
  ScPri = ScalarProduct( NewItem , primal , NewItemBse );

  // compute the right hand side of the constraint - - - - - - - - - - - - -

  if( Aset )         // newitem is sparse & is using the active set
   Ai -= ScalarProduct( NewItem , NewItemBse , CrrPnt , Aset );
  else               // newitem is sparse & is not using the active set
   Ai -= ScalarProduct( NewItem , CrrPnt , NewItemBse );
  }
 else {              // dense format
  ScPri = ScalarProduct( NewItem , primal , CrrSGLen );

  // compute the right hand side of the constraint - - - - - - - - - - - - -

  if( Aset )         // newitem is dense & activeset
   Ai -= ScalarProduct( CrrPnt , NewItem , Aset );
  else               // newitem is dense & not activeset
   Ai -= ScalarProduct( NewItem , CrrPnt , CrrSGLen );
  }

 NewItemprice = Ai;
 NewItemScPri = ScPri;
 NewItemisSG = false;  // this item is a constraint ...

 // Note!!!! It doesn't check whether the item is or not identical to some
 // other items

 return( Inf<Index>() );

 }  // end( OSIMPSolver::CheckCnst )

/*------------------------------------------------------------------------*/

bool OSIMPSolver::ChangesMPSol( void )
{
 MSG( "OSIMPSolver::ChangesMPSol() called\n" );

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
   if( iseasy( NewItemFi ) )
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

/*------------------------------------------------------------------------*/

void OSIMPSolver::SetItem( cIndex Nm )
{
 // note: the columns of the coefficient matrix (for the relevant part)
 //       contain the *opposite* of the subgradient

 MSG("OSIMPSolver::SetItem() \n");

 if( ( ! NewItemFi ) && ( Nm == Inf<Index>() ) ) {  // the 0th component
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

  //update the dictionaries  - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < dict_item_maxname ; i++ )
   if( ( dict_item[ i ] < Inf<Index>() ) &&
       ( dict_item[ i ] > Index( col ) ) )
    dict_item[ i ]--;

  for( Index i = 0 ; i < CrrSGLen ; i++ )
   if( ( dict_slack[ i ] < Inf<Index>() ) &&
       ( dict_slack[ i ] > Index( col ) ) )
    dict_slack[ i ]--;

  if( dict_stab ) {
   for( Index i = 0 ; i < CrrSGLen ; i++ )
    if( ( dict_stab[ i ] < Inf<Index>() ) &&
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

   // and possibly fix gamma[ NewItemFi ] = 0- - - - - - - - - - - - - - - -

   if( ( osiSlvr->getColUpper()[ comp_col[ NewItemFi ] ] > 0.0 ) &&
       ( osiSlvr->getObjCoefficients()[ comp_col[ NewItemFi ] ] == 0.0 ) )
    osiSlvr->setColUpper( comp_col[ NewItemFi ] , 0.0 );
   }

  // set the name and mark it as a new item- - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  wcomp[ Nm ] = NewItemFi;
  newly_inserted[ Nm ] = just_inserted[ Nm ] = true;

  // add the column to the master probem - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->addCol( int( i ) , tempI , tempHP , 0.0 , osiSlvr->getInfinity() ,
	       - NewItemprice );

  // put the item in the dictionary  - - - - - - - - - - - - - - - - - - - -

  dict_item[ Nm ] = osiSlvr->getNumCols () - 1;

  MSG( "Item " << Nm << " of the component " << NewItemFi
	<< " has the entry " << dict_item[ Nm ] << endl);

  if( dict_item_maxname < Nm + 1 )
   dict_item_maxname = Nm + 1;
  }
 }  // end( OSIMPSolver::SetItem )

/*------------------------------------------------------------------------*/

void OSIMPSolver::SubstItem( cIndex Nm )
{
 MSG( "OSIMPSolver::SubstItem()\n" );

 if( - osiSlvr->getObjCoefficients()[ dict_item[ Nm ] ] > NewItemprice )
  osiSlvr->setObjCoeff( dict_item[ Nm ] , - NewItemprice );

 } // end( OSIMPSolver::SubstItem )

/*------------------------------------------------------------------------*/

void OSIMPSolver::RmvItem( cIndex i )
{
 MSG( "OSIMPSolver::RmvItem()\n" );

 Index Nm;
 const int index = dict_item[ i ];  // mark the column deleted - - - - - - -
 osiSlvr->deleteCols( 1 , &index );

 #if( PRESERVE_OSI_SOLS )
  // having deleted the column 'index' from the osiSlvr, all the stored solution
  // information corresponding to columns (csol[] and rcst[]) must be shifted
  // left by one from position index onwards

  ShiftVect( csol + index , csols - index - 1 );
  ShiftVect( rcst + index , csols - index - 1 );
  csols--;
 #endif

 Index wFi = wcomp[ i ];

 // update the item's vocabulary - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 dict_item[ i ] = Inf<Index>();
 wcomp[ i ] = Inf<Index>();

 if( i + 1 == dict_item_maxname )
  dict_item_maxname--;

 // if the deleted item is the last one for its component, the convexity
 // constraint for that component has to be deactivated  - - - - - - - - - -

 if( ( ! ( MaxName( wFi ) ) ) &&
       ( osiSlvr->getColUpper()[ comp_col[ wFi ] ] == 0.0 ) ) {
  osiSlvr->setColUpper( comp_col[ wFi ] , osiSlvr->getInfinity() );
  MSG( "removed item " << i << " relative to the component " << wFi << endl );
  }

 // update the dictionaries - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
  if( ( dict_item[ Nm ] != Inf<Index>() ) &&
      ( dict_item[ Nm ] > Index( index ) ) )
   dict_item[ Nm ]--;

 for( Nm = 0 ; Nm < CrrSGLen ; Nm++ ) {
  if( ( dict_slack[ Nm ] != Inf<Index>() ) &&
      ( dict_slack[ Nm ] > Index( index ) ) )
   dict_slack[ Nm ]--;
  if( stab == quadratic )
   if( ( dict_stab[ Nm ] != Inf<Index>() ) &&
       ( dict_stab[ Nm ] > Index( index ) ) )
    dict_stab[ Nm ]--;
  }

 if( comp_col[ 0 ] > Index( index ) )
  comp_col[ 0 ]--;

 }  // end( OSIMPSolver::RmvItem )

/*------------------------------------------------------------------------*/

void OSIMPSolver::RmvItems( void )
{
 MSG( "OSIMPSolver::RmvItems()\n" );

 if( ! dict_item_maxname )  // there is nothing to do ...
  return;

 int count = 0;
 resizeI( dict_item_maxname );

 for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
  if( dict_item[ Nm ] != Inf<Index>() )
   tempI[ count++ ] = dict_item[ Nm ];

 // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->deleteCols( count , tempI );

 // fill the vectors for the shifting of the vocabularies  - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index elem;
 for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
  if( dict_item[ Nm ] != Inf<Index>() ) {

   elem = dict_item[ Nm ];

   // update the left part of tempI  - - - - - - - - - - - - - - - - - - - -

   for( Index i = Nm + 1; i < dict_item_maxname ; i++ )
    if( ( dict_item[ i ] > elem ) && ( dict_item[ i ] != Inf<Index>() ) )
     dict_item[ i ]--;

   // update the position of rho column  - - - - - - - - - - - - - - - - - -

   if( comp_col[ 0 ] > elem )
	comp_col[ 0 ]--;

   // update the slack and stabilization dictionaries  - - - - - - - - - - -

   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( ( dict_slack[ i ] != Inf<Index>() ) && ( dict_slack[ i ] > elem ) )
     dict_slack[ i ]--;
    if( stab == quadratic )
     if( ( dict_stab[ i ] != Inf<Index>() ) && ( dict_stab[ i ] > elem ) )
      dict_stab[ i ]--;
    }

   } // end check  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < NrFi ; i++ )
  if( osiSlvr->getColUpper()[ comp_col[ wcomp[ i ] ] ] == 0.0 )
   osiSlvr->setColUpper( comp_col[ wcomp[ i ] ] , osiSlvr->getInfinity() );

 VectAssign( dict_item , Index( Inf<Index>() ) , dict_item_maxname );
 VectAssign( wcomp , Index( Inf<Index>() ) , dict_item_maxname );
 dict_item_maxname = 0;

 }  // end( OSIMPSolver::RmvItems )

/*------------------------------------------------------------------------*/

void OSIMPSolver::SetActvSt( cIndex_Set AVrs , cIndex AVDm )
{
 MSG("OSIMPSolver::SetActvSt()\n" );

 if( ! useactiveset )
  throw( NDOException( "OSIMPSolver::SetActvSt(): Active set not declared" ) );

 if( ! AVrs ) {	 //  all the variables will be not active  - - - - - - - - -
  while( Asetdim > 0 )
   deactivate( Aset[ Asetdim-- ] );
  }
 else
  if( ! Aset ) {  // if Aset is empty, define all the variables in AVrs as active
    while ( Asetdim < AVDm )
     activate( AVrs[ Asetdim++ ] );
   }
  else
   for( cIndex *newA = AVrs ; ( *newA < Inf<Index>() ) ||
	                      ( *Aset < Inf<Index>() ) ; ) {
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

/*------------------------------------------------------------------------*/

void OSIMPSolver::AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs )
{
 MSG( "OSIMPSolver::AddActvSt() called\n" );

 if( ! useactiveset )
  throw( NDOException( "OSIMPSolver::AddActvSt(): Active set not declared" ) );

 Aset = AVrs;
 Asetdim += AdDm;

 for( ; *Addd < Inf<Index>() ; Addd++ )
  activate( *Addd );

 }  // end( OSIMPSolver::AddActvSt )

/*------------------------------------------------------------------------*/

void OSIMPSolver::RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs )
{
 MSG( "OSIMPSolver::RmvActvSt() called\n" );

 Aset = AVrs;
 Asetdim = Asetdim - RmDm;

 for( ; *Rmvd < Inf<Index>() ; Rmvd++ )
  deactivate( *Rmvd );

 } // end( OSIMPSolver::RmvActvSt )

/*------------------------------------------------------------------------*/

void OSIMPSolver::AddVars( cIndex NNwVrs )
{
 MSG( "OSIMPSolver::AddVars() called\n" );

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
  if( dict_item[ j ] < Inf<Index>() ) {  // ask the components [ CrrSGLen ,
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
   if( dict_item[ j ] < Inf<Index>() ) {
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
  MSG( "Add constraints for the primal variables as inactive ones\n" );
  for( Index j = 0 ; j < NNwVrs ;  j++ ) {
   rowlb[ j ] = - osiSlvr->getInfinity();
   rowub[ j ] = osiSlvr->getInfinity();
   }
  osiSlvr->addRows( NNwVrs , rowStarts, tempI , tempHP , rowlb , rowub );
  }
 else {
  MSG( "Add constraints for the primal variables \n" );
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

  if( ( Upper[ i ] = FIO->GetUB( i ) ) != Inf<LMNum>() )
   slack_m = true;

  if( FIO->GetUC( i ) )
   Lower[ i ] = -Inf<LMNum>();
  else {
   Lower[ i ] = 0;
   slack_p = true;
   }

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
    MSG( "Add z_"<< i << endl );
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
   MSG( "Add the slack s_"<< i <<"+\n" );
   }

  if( slack_m )	{
   columnStarts[ count ] = count;

   // if slack_p is present, slack_m is not indicated
   if( ! slack_p )
    dict_slack[ i ] = count + nc;

   rowlb[ count ] = 0;
   tempI[ count ] = comp_row[ NrFi ] + i;
   tempHP[ count++ ] = -1.0;
   MSG( "Add the slack s_" << i << "-\n" );
   }
  }

 columnStarts[ count ] = count;  // set the end both of vector tempHp and tempI
 rowub = new HpNum[ count ];
 obj = new HpNum[ count ];
 for( Index j = 0 ; j < count ;  j++ ) {
  rowub[ j ] = osiSlvr->getInfinity();
  obj[ j ] = 0;
  }

 osiSlvr->addCols( count , columnStarts , tempI , tempHP , rowlb , rowub , obj );

 // deallocate the memory - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] columnStarts;
 delete[] rowlb;
 delete[] rowub;
 delete[] obj;

 //update the current items length   - - - - - - - - - - - - - - - - - - - - -

 CrrSGLen = NewSGLen;

 // change the price of the slack variables - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 tUpdatePrices( CrrSGLen - NNwVrs , CrrSGLen );
 // change the prices for the stabilization

 if( stab != boxstep )  // change the slack prices
  ptUpdatePrices( CrrSGLen - NNwVrs , CrrSGLen );

 }  // end( OSIMPSolver::AddVars )

/*------------------------------------------------------------------------*/

void OSIMPSolver::RmvVars( cIndex_Set whch , Index hwmny )
{
 MSG( "OSIMPSolver::RmvVars() called\n" );

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
   if( dict_slack[ whch[ i ] ] != Inf<Index>() ) {
	tempI[ count++ ] = dict_slack[ whch[ i ] ];
    if( ( ( Lower[ whch[ i ] ] > -Inf<LMNum>() ) &&
	  ( Upper[ whch[ i ] ] < +Inf<LMNum>() ) ) || ( stab == boxstep ) )
     tempI[ count++ ] = dict_slack[ whch[ i ] ] + 1;
    }
   if( dict_stab )
    tempI[ count++ ] = dict_stab[ whch[ i ] ];
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
    if( ( dict_item[ i ] != Inf<Index>() ) && ( dict_item[ i ] > elem ) )
     dict_item[ i ]--;

   // update the slack and stabilization dictionaries  - - - - - - - - - - -

   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( ( dict_slack[ i ] != Inf<Index>() ) && ( dict_slack[ i ] > elem ) )
     dict_slack[ i ]--;
    if( stab == quadratic )
     if( ( dict_stab[ i ] != Inf<Index>() ) && ( dict_stab[ i ] > elem ) )
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

  // remove all rows in the dynamic part  - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < CrrSGLen ; i++ )
   tempI[ i ] = comp_row[ NrFi ] + i;

  osiSlvr->deleteRows( CrrSGLen , tempI );

  // mark the columns to delete   - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( dict_stab )
   resizeI( 3 * CrrSGLen + dict_item_maxname );
  else
   resizeI( 2 * CrrSGLen + dict_item_maxname );

  count = 0;
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   if( dict_slack[ i ] != Inf<Index>() ) {
	tempI[ count++ ] = dict_slack[ i ];
    if( ( Lower[ i ] > -Inf<LMNum>() ) &&  ( Upper[i] < +Inf<LMNum>() ) )
     tempI[ count++ ] = dict_slack[ i ] + 1;
    }
    if( dict_stab )
     tempI[ count++ ] = dict_stab[ i ];
   }

  for( Index Nm = 0 ; Nm < dict_item_maxname ; Nm++ )
   if( dict_item[ Nm ] != Inf<Index>() )
    tempI[ count++ ] = dict_item[ Nm ];

  // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->deleteCols( count , tempI );

  // update the 0-th component  - - - - - - - - - - - - - - - - - - - - - - -

  RhoColBDm = comp_row[ NrFi ];

  for( Index Nm = 0 ; Nm < Index(count) ; Nm++ )
    if( comp_col[ 0 ] > tempI[ count ] )
     comp_col[ 0 ] --;

  // delete all the items  - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < NrFi ; i++ )
   if( osiSlvr->getColUpper()[ comp_col[ wcomp[ i ] ] ] == 0.0 )
    osiSlvr->setColUpper( comp_col[ wcomp[ i ] ] , osiSlvr->getInfinity() );

  VectAssign( dict_item , Index( Inf<Index>() ) , dict_item_maxname );
  VectAssign( wcomp , Index( Inf<Index>() ) , dict_item_maxname );
  VectAssign( dict_slack , Index( Inf<Index>() ) , MaxSGLen );

  dict_item_maxname = 0;

  CrrSGLen = 0;

  } // end removing all elements - - - - - - - - - - - - - - - - - - - - - -
 }  // end( OSIMPSolver::RmvVars )

/*------------------------------------------------------------------------*/

void OSIMPSolver::ChgAlfa( cHpRow NewAlfa , cIndex wFi )
{
 // note: the obj coefficients of each item column are *the opposite* of
 //       the corresponding linearization error / RHS

 MSG( "OSIMPSolver::ChgAlfa() called\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ChgAlfa( * , 0 ) called" ) );

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( ( dict_item[ i ] < Inf<Index>() ) &&
      ( ( wFi > NrFi ) || ( wFi == wcomp[ i ] ) ) )
    osiSlvr->setObjCoeff( dict_item[ i ] , - NewAlfa[ i ] );

 } // end( OSIMPSolver::ChgAlfa )

/*------------------------------------------------------------------------*/

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

 MSG( "OSIMPSolver::ChangeCurrPoint( DLabda ) called\n" );

 // no items are new - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( newly_inserted , false , dict_item_maxname );

 // change the coefficients of the objective function- - - - - - - - - - - -

 const double* old = osiSlvr->getObjCoefficients();

 // update the either the linearization error or the rhs of the item - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( dict_item[ i ] < Inf<Index>() ) {
   double temp = - old[ dict_item[ i ] ];
   CoinShallowPackedVector item =
     ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

   int k;
   for( Index j = item.getNumElements() ; j-- > 0 ; ) {
    if( ( k = item.getIndices()[ j ] - comp_row[ NrFi ] ) >= 0 )
     temp += item.getElements()[ j ] * DLambda[ k ];
    }

   if( IsSubG( i ) )
    temp += DFi[ wcomp[ i ] ];
   else {
    // it is a constraint: ensure its RHS remains non-negative, as small
    // negative RHS which may crop up by numerical errors would make
    // the primal master problem unfeasible (the dual unbounded)

    if( temp < 0 )
     temp = 0;
    }

   osiSlvr->setObjCoeff( dict_item[ i ] , - temp );
   }

 // change the cost of the easy components - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 1 ; i <= NrFi ; i++ ) {
  Index BNC;
  if( ( BNC = FIO->GetBNC( i ) ) )
   for( Index j = 0 ; j < BNC ; j++ ) {
    double temp = old[ comp_col[ i ] + j ];
    CoinShallowPackedVector item =
      ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ i ] + j );

    for( Index l = item.getNumElements() ; l > 0 ; l-- ) {
     int k;
     if( ( k = item.getIndices()[ l - 1 ] - comp_row[ NrFi ] ) >= 0 )
      temp -= item.getElements()[ l - 1 ] * DLambda[k];
     }

    osiSlvr->setObjCoeff( comp_col[ i ] + j , temp );
    }
  }

 // change the bounds of the variables - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < CrrSGLen ; i++ ) {
  if( Upper[ i ] < Inf<LMNum>() )
   Upper[ i ] -= DLambda[ i ];
  if( Lower[ i ] > -Inf<LMNum>() )
   Lower[ i ] -= DLambda[ i ];
  }

 // change lower bound of the difficult components - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 1 ; i <= NrFi ; i++ )
  if( ! iseasy( i ) ) {
   double temp = old[ comp_col[ i ] ];
   if( osiSlvr->getColUpper()[ comp_col[ i ] ] > 0 ) {
    temp -= DFi[ i ];   // ????? sign - ????
    osiSlvr->setObjCoeff( comp_col[ i ] , temp );
    }
   }

 // change the total lower bound - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( HasLwrBnd ) {  // ... if it exists
  double temp = old[ comp_col[ 0 ] ];
  temp += DFi[ 0 ];
  for( Index i = 1 ; i <= NrFi ; i++ )
   if( ! iseasy( i ) )
    temp += DFi[ i ];  // ????? sign + ????

  osiSlvr->setObjCoeff( comp_col[ 0 ] , temp );
  }

 // update the price of slack - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ptUpdatePrices();

 // set the value of the easy components in the current point - - - - - - -

 VectAssign( FiLambda , HpNum( Inf<HpNum>() ) , NrFi );

 }  // end( OSIMPSolver::ChangeCurrPoint( DLambda , DFi ) )

/*------------------------------------------------------------------------*/

void OSIMPSolver::ChangeCurrPoint( cHpNum Tau , cHpRow DFi )
{
 MSG( "OSIMPSolver::ChangeCurrPoint( Tau ) called\n" );

 VectAssign( newly_inserted , false , dict_item_maxname );
 const double *old = osiSlvr->getObjCoefficients();

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *delta = rsol + comp_row[ NrFi ];

 for( Index i = 0 ; i < dict_item_maxname ; i++ )
  if( dict_item[ i ] < Inf<Index>() ) {
   double temp = - old[ dict_item[ i ] ];
   temp -= ( Tau / t ) * ReadGid( i );

   if( IsSubG( i ) )
    temp += DFi[ wcomp[ i ] ];
   else { // it is a constrtaint: ensure its RHS remains non-negative, as
    // small negative RHS which may crop up by numerical errors would make
    // the primal master problem unfeasible (the dual unbounded)

    if( temp < 0 )
     temp = 0;
    }

   osiSlvr->setObjCoeff( dict_item[ i ] , - temp );
   }

 // change the cost of the easy components - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 1 ; i <= NrFi ; i++ ) {
  Index BNC;
  if( ( BNC = FIO->GetBNC( i ) ) )
   for( Index j = 0 ; j < BNC ; j++ ) {
    double temp = old[ comp_col[ i ] + j ];
    CoinShallowPackedVector item =
      ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ i ] + j );

    for( Index l = item.getNumElements() ; l > 0  ; l-- ) {
     int k;
     if( ( k = item.getIndices()[ l - 1 ] - comp_row[ NrFi] ) >= 0 )
      temp -= item.getElements()[ l - 1 ] * delta[ k ] * (Tau / t);
     }

    osiSlvr->setObjCoeff( comp_col[ i ] + j , temp );
    }
  }

 // change the bounds of the variables - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cLMRow d = Readd( true );
 for( Index i = 0 ; i < CrrSGLen ; i++ ) {
  if( Upper[ i ] < Inf<LMNum>() )
   Upper[ i ] -= d[ i ] * ( Tau / t );
  if( Lower[ i ] > -Inf<LMNum>() )
   Lower[ i ] -= d[ i ] * ( Tau / t );
  }

 // change lower bound of the difficult components - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 1 ; i <= NrFi ; i++ )
  if( ! iseasy( i ) ) {
   double temp = old[ comp_col[ i ] ];
   if( osiSlvr->getColUpper()[ comp_col[ i ] ] > 0 ) {
    temp -= DFi[ i ]; // ????? sign - ????
    osiSlvr->setObjCoeff( comp_col[ i ] , temp );
    }
   }

 // change the total lower bound - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( HasLwrBnd ) {  // ... if it exists
  double temp = old[ comp_col[ 0 ] ];
  temp += DFi[ 0 ];
  for( Index i = 1 ; i <= NrFi ; i++ )
   if( ! iseasy( i ) )
    temp += DFi[ i ]; // ????? sign - ????
  osiSlvr->setObjCoeff( comp_col[ 0 ] , temp );
  }

 // update the price of slack - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ptUpdatePrices();

 // set the value of the easy components in the current point
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Tau != t )
  VectAssign( FiLambda , HpNum( Inf<HpNum>() ) , NrFi );
 else
  VectAssign( FiLambda , FiLambda1 , NrFi);
 }  // end( OSIMPSolver::ChangeCurrPoint( Tau , DFi ) )

/*------------------------------------------------------------------------*/

void OSIMPSolver::ChgSubG( cIndex strt , Index stp , cIndex wFi )
{
 MSG( "OSIMPSolver::ChgSubG() called\n" );

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
/*------------------------------------------------------------------------*/

OSIMPSolver::~OSIMPSolver()
{
 MSG( "OSIMPSolver::~OSIMPSolver() called\n" );

 cleanup();

 delete osiSlvr;

 }  // end( OSIMPSolver::~OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void OSIMPSolver::cleanup( void )
{
 if( osiSlvr ) {
  osiSlvr->reset();
  osiSlvr->messageHandler()->setLogLevel( int( MPLLvl ) );
  }

 // delete the dictionaries everything

 delete[] dict_item;
 dict_item = 0;
 delete[] dict_slack;
 dict_slack = 0;
 delete[] dict_stab;
 dict_stab = 0;
 delete[] wcomp;
 wcomp = 0;
 delete[] newly_inserted;
 newly_inserted = 0;
 delete[] just_inserted;
 just_inserted = 0;
 delete[] comp_row;
 comp_row = 0;
 delete[] comp_col;
 comp_col = 0;
 dict_item_maxname = 0;

 delete[] RhoCol;
 RhoCol = 0;
 delete[] RhoColBse;
 RhoColBse = 0;
 RhoColBDm = 0;

 delete[] NewItem;
 delete[] tempHP;
 tempHP_size = 0;
 delete[] tempI;
 tempI_size = 0;
 delete[] Lower;
 Lower = 0;
 delete[] Upper;
 Upper = 0;

 delete[] GiPerd;
 GiPerd = 0;
 delete[] FiLambda1;
 FiLambda1 = 0;
 delete[] FiLambda;
 FiLambda = 0;

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
 if( stp > CrrSGLen )
  stp = CrrSGLen;

 switch( stab ) {
  case none: break;
  case boxstep:
   if( dict_slack )
    for( Index i = strt ; i < stp ; i++ ) {
     osiSlvr->setObjCoeff( dict_slack[ i ] ,
			   ( Lower[ i ] > -t ) ? Lower[ i ] : -t );
     osiSlvr->setObjCoeff( dict_slack[ i ] + 1 ,
			   ( t > Upper[ i ] ) ? -Upper[ i ] : -t );
     }
   break;
  case quadratic:
   if( dict_stab )
    for( Index i = strt ; i < stp ; i++ )
     ChgQCoef( - t , dict_stab[ i ] );
   break;
  default:
   throw( NDOException( "OSIMPSolver: undecided stabilization" ) );
  }
 }  // end( OSIMPSolver::tUpdatePrices )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ptUpdatePrices( cIndex strt , Index stp )
{
 // distinguish box stabilization case because it depends on t too- - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( stp > CrrSGLen )
  stp = CrrSGLen;

 if( stab != boxstep ) {
  if( dict_slack )
   for( Index i = strt ; i < stp ; i++ ) {
    int offset = 0;  // if any lower bound exists, there is a slack s_i^+
    if( Lower[ i ] > -Inf<LMNum>() ) {
     osiSlvr->setObjCoeff( dict_slack[ i ] , Lower[ i ] );
     offset = 1;
     }
    if( Upper[i] < +Inf<LMNum>() )
     osiSlvr->setObjCoeff( dict_slack[ i ] + offset , -Upper[ i ] );
    }
  }
 else                 // boxstep case ...
  if( dict_slack )
   for( Index i = strt ; i < stp ; i++ ) {
    osiSlvr->setObjCoeff( dict_slack[ i ] ,
			  ( Lower[ i ] > -t ) ? Lower[ i ] : -t );
    osiSlvr->setObjCoeff( dict_slack[ i ] + 1 ,
			  ( t > Upper[ i ] ) ? -Upper[ i ] : -t );
    }

 }  // end( OSIMPSolver::ptUpdatePrices )

/*------------------------------------------------------------------------*/

void OSIMPSolver::switchToQP( void )
{
 #if CPLEX == 1
  OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
  if( ! osiCpx )
   throw( NDOException( "OSIMPSolver::switchToQP: the OSI solver is not Cplex" ) );

  CPXENVptr env = osiCpx->getEnvironmentPtr();
  CPXLPptr qp = osiCpx->getLpPtr();

  if( CPXchgprobtype( env, qp, CPXPROB_QP ) )
   throw( NDOException( "OSIMPSolver::switchToQP: can't turn to QP problem" ) );
 #else
  throw( NDOException( "OSIMPSolver::switchToQP: not implemented yet" ) );
 #endif

 } // end( OSIMPSolver::switchToQP )

/*------------------------------------------------------------------------*/

void OSIMPSolver::ChgQCoef( HpNum value , Index nm )
{
 #if CPLEX == 1
  OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
  if( ! osiCpx )
   throw( NDOException( "OSIMPSolver::ChgQCoef: the OSI solver is not Cplex" ) );

  CPXENVptr env = osiCpx->getEnvironmentPtr();
  CPXLPptr qp = osiCpx->getLpPtr();

  CPXchgqpcoef( env , qp , nm , nm , value );
 #else
  throw( NDOException( "OSIMPSolver::ChgQCoef: not implemented yet" ) );
 #endif

 } // end( OSIMPSolver::ChgQCoef )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::GetQCoef( Index nm )
{
 #if CPLEX == 1
  HpNum qcoef;
  OsiCpxSolverInterface *osiCpx = dynamic_cast<OsiCpxSolverInterface*>( osiSlvr );
  if( ! osiCpx )
   throw( NDOException( "OSIMPSolver::GetQCoef: the OSI solver is not Cplex" ) );

  CPXENVptr env = osiCpx->getEnvironmentPtr();
  CPXLPptr qp = osiCpx->getLpPtr();

  CPXgetqpcoef( env , qp , nm , nm , &qcoef );
  return( qcoef );

 #else
  throw( NDOException( "OSIMPSolver::GetQCoef: not implemented yet" ) );
 #endif

 }  // end( OSIMPSolver::GetQCoef )

/*------------------------------------------------------------------------*/

inline bool OSIMPSolver::isactive( Index i )
{
 return( osiSlvr->getRowSense()[ comp_row[ NrFi ] + i ] == 'E' );
 }

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
 }  // end( OSIMPSolver::resizeHP )

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::resizeI( Index i )
{
 if( i > tempI_size ) {
  delete[] tempI;
  tempI = new int[ tempI_size = i ];
  }
 }  // end( OSIMPSolver::resizeI )

/*------------------------------------------------------------------------*/

inline bool OSIMPSolver::iseasy( Index i )
{
 return( osiSlvr->getColLower()[ comp_col[ i ] ] < 0.0 );
 }

/*------------------------------------------------------------------------*/

inline void OSIMPSolver::CheckDS( void )
{
 #if CHECK_DS & 1
 // obvious sanity check in dict_item[]
 if( osiSlvr ) {
  Index numcols = osiSlvr->getNumCols();
  for( Index name = 0 ; name < dict_item_maxname ; name++ )
   if( dict_item[ name ] < Inf<Index>() )
    if( dict_item[ name ] > numcols )
     cout << "dict_item[ " << name << " ] = " << dict_item[ name ]
	  << " out of range (> " << numcols << ")" << endl;
  }
 #endif
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File OSIMPSolver.C -------------------------*/
/*--------------------------------------------------------------------------*/

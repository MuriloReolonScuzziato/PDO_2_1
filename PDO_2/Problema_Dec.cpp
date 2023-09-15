/*--------------------------------------------------------------------------*/
/*---------------------------- File Problema_Dec.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para resolver o problema por decomposição
 *
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Problema_Dec.h"

// NDOSolver header(s) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#if( WHICH_NDOSOLVER == 0 )
	#include "Bundle.h"
#elif( WHICH_NDOSOLVER == 1 )
	#include "Bundle.h"
	#include "Volume.h"
#elif( WHICH_NDOSOLVER == 2 )
	#include "Volume.h"
#else
	#error "Unable to find the solver"
#endif

// MPSolver header(s), if used - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
#if WHICH_NDOSOLVER <= 1
	#if WHICH_MPSOLVER == 0
		#include "QPPnltMP.h"
	//#elif WHICH_MPSOLVER == 1
	//	#include "OSIMPSolver.h"
	//	#include "OsiCpxSolverInterface.hpp"
	#elif WHICH_MPSOLVER >= 1
		#include "OSIMPSolver.h"
		#if WHICH_RMPSOLVER == 0
			#include "OsiCpxSolverInterface.hpp"
			//typedef OsiCpxSolverInterface RealSolverInterface;
		#elif WHICH_RMPSOLVER == 1
			#include "OsiClpSolverInterface.hpp"
			//typedef OsiClpSolverInterface RealSolverInterface;
		#else
			#include "OsiGrbSolverInterface.hpp"
			//typedef OsiGrbSolverInterface RealSolverInterface;
		#endif
	#else
		#error "Unable to find the MPsolver."
	#endif
#endif

#include <fstream>
#include <sstream>

#include <windows.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
	using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// parameter file name

#if( WHICH_NDOSOLVER == 0 )
	//const char *const logF = "log.bn";
	//const char *const logMP = "MPlog.bn";
	#if( WHICH_MPSOLVER )
		const char *const ParF = "ParValue.osi";
	#else
		const char *const ParF = "ParValue.qp";
	#endif
#else
	const char *const ParF = "ParValue.vol";
#endif

// Output file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

const char *const out = "report.txt";

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// set algorithmic parameters for Cplex, if adopted- - - - - - - - - - - - - -

#if WHICH_NDOSOLVER <= 1
 #if (WHICH_MPSOLVER)
  #if WHICH_RMPSOLVER == 0

void SetCplexA( OsiCpxSolverInterface * osiCpx , const int algorithm = 0 ,
	  const int reduction = 3 )
{

 CPXENVptr env = osiCpx->getEnvironmentPtr ();  // Cplex environment

 // choose CPLEX algorithm - - - - - - - - - - - - - - - - - - - - - - - - - -

 switch ( algorithm ) {
  case ( 0 ): //Automatic: let CPLEX choose (default)
	  CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_AUTOMATIC );
    break;
  case ( 1 ): // Primal simplex
	CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_PRIMAL );
    break;
  case ( 2 ): // Dual simplex
    CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_DUAL );
    break;
  case ( 3 ): // Network simplex
    CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_NET );
    break;
  case ( 4 ): // Barrier
    CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_BARRIER );
    break;
  case ( 5 ): // Sifting
   CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_SIFTING );
   break;
  case ( 6 ): // Concurrent (Dual, Barrier, and Primal)
   CPXsetintparam( env , CPX_PARAM_QPMETHOD  , CPX_ALG_CONCURRENT );
   break;
  default:
   throw( NDOException( "Main: does not select any algorithm" ) );
  }

 // Primal and dual reduction type - - - - - - - - - - - - - - - - - - - - - -

 switch ( reduction ) {
  case ( 0 ): // No primal or dual reductions
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_NOPRIMALORDUAL );
   break;
  case ( 1 ): // Only primal reductions
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALONLY );
   break;
  case ( 2 ): // CPX_PREREDUCE_DUALONLY
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_DUALONLY );
   break;
  case ( 3 ): // Both primal and dual reductions (default)
   CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALANDDUAL );
   break;
  default:
   throw( NDOException( "Main: does not select any reduction" ) );
  }

 } // end ( SetCplexA )

  #elif WHICH_RMPSOLVER == 2

void SetGrbA( OsiGrbSolverInterface * osiGrb , const int algorithm = 0 ,
	  const int reduction = 3 )
{

	GRBenv * env = osiGrb->getEnvironmentPtr();  // Gurobi environment

 // choose Gurobi algorithm - - - - - - - - - - - - - - - - - - - - - - - - - -

 switch ( algorithm ) {
  case ( 0 ): //Automatic: let Gurobi choose (default)
	GRBsetintparam( env, "Method", GRB_METHOD_AUTO);
    break;
  case ( 1 ): // Primal simplex
	GRBsetintparam( env, "Method", GRB_METHOD_PRIMAL);
    break;
  case ( 2 ): // Dual simplex
	GRBsetintparam( env, "Method", GRB_METHOD_DUAL);
    break;
  case ( 3 ): // Barrier
	GRBsetintparam( env, "Method", GRB_METHOD_BARRIER);
    break;
  case ( 4 ): // Concurrent (Dual, Barrier, and Primal)
   GRBsetintparam( env, "Method", GRB_METHOD_CONCURRENT);
   break;
  default:
   throw( NDOException( "Main: does not select any algorithm" ) );
  }

 // Primal and dual reduction type - - - - - - - - - - - - - - - - - - - - - -

 //switch ( reduction ) {
 // case ( 0 ): // No primal or dual reductions
 //  CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_NOPRIMALORDUAL );
 //  break;
 // case ( 1 ): // Only primal reductions
 //  CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALONLY );
 //  break;
 // case ( 2 ): // CPX_PREREDUCE_DUALONLY
 //  CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_DUALONLY );
 //  break;
 // case ( 3 ): // Both primal and dual reductions (default)
 //  CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALANDDUAL );
 //  break;
 // default:
 //  throw( NDOException( "Main: does not select any reduction" ) );
 // }

 } // end ( SetGrbA )

  #endif
 #endif
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

Problema_Dec::Problema_Dec(CSistema * const sistema_end, string ID)
{
	// Fazer essa classe generica (interface para a decomposiçao espacial e por cenários) -> usar typedef
	sistema_a = sistema_end;
	resultadosGurobi = new Resultados(sistema_a);
	IDpar = ID;
	startt = NULL;
}
Problema_Dec::~Problema_Dec(void)
{
	delete startt;
	delete resultadosGurobi;
}

int Problema_Dec::ResolverRL(string arv, string real)
{
	// Log files
	string log = "log_" + arv + "_" + real + "_";
	string MPlog = "MPlog";

	//// Log files
	//string log = "log_" + arv + "_" + real;
	//string MPlog = "MPlog_" + arv + "_" + real;
	//#if (WHICH_MPSOLVER == 0)
	//	log.append("_QPP");
	//	MPlog.append("_QPP");
	//#elif (WHICH_MPSOLVER == 1)
	//	log.append("_CPLEXQ");
	//	MPlog.append("_CPLEXQ");
	//#else
	//	#if( WHICH_RMPSOLVER == 0 )
	//		log.append("_CPLEXL");
	//		MPlog.append("_CPLEXL");
	//	#elif ( WHICH_RMPSOLVER == 1 )
	//		log.append("_CLP");
	//		MPlog.append("_CLP");
	//	#elif ( WHICH_RMPSOLVER == 2 )
	//		log.append("_GRB");
	//		MPlog.append("_GRB");
	//	#endif
	//#endif
	log.append(IDpar);
	log.append(".bn");
	MPlog.append(".bn");
	const char *const logF = log.c_str();
	const char *const logMP = MPlog.c_str();

	// enter the try-block - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	try 
	{
		// open parameters file - - - - - - - - - - - - - - - - - - - - - - - - - -

		ofstream Out( out , ofstream::app );
		if( ! Out.is_open() )
			cerr << "Warning: cannot open log file """ << out << """" << endl;
		// construct the FiOracle object- - - - - - - - - - - - - - - - - - - - - -
		#if( WHICH_DEC == 1 )
			Spcdec1Fi *Fi = new Spcdec1Fi( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#elif( WHICH_DEC == 2 )
			Spcdec2Fi *Fi = new Spcdec2Fi( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#elif( WHICH_DEC == 3 )
			Spcdec3Fi *Fi = new Spcdec3Fi( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#elif( WHICH_DEC == 4 )
			Scndec1Fi *Fi = new Scndec1Fi( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#elif( WHICH_DEC == 5 )
			Scndec2Fi *Fi = new Scndec2Fi( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#elif( WHICH_DEC == 6 )
			Scndec3Fi *Fi = new Scndec3Fi( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#elif( WHICH_DEC == 0 )
			EspFiOracle *Fi = new EspFiOracle( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#endif
//		cout << "Oracle criado!" << endl;
		// Criar modelos do oracle - - - - - - - - - - - - - - - - - - - - - - - - - - 

		//Fi->SetData(cntr , UB);	// Definir limites para os Lambdas e the minimum of the function (L0)
		Fi->CriarModelos();
		//cout << "Modelos do Oracle criados!" << endl;
		
		// Add the identification to the parameters file - - - - - - - - - - - - - 
		string ParFN;
		int ip = 0;
		while (ParF[ip] != '\0')
		{
			if (ParF[ip] == '.')
				ParFN.append(IDpar);
			ParFN.push_back(ParF[ip]);
			ip++;
		}

		// open parameters file - - - - - - - - - - - - - - - - - - - - - - - - - -
		//ifstream ParFile( ParF );
		ifstream ParFile( ParFN );
		if( ! ParFile.is_open() )
			cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;

		// construct the NDOSolver object - - - - - - - - - - - - - - - - - - - - -

		#if( WHICH_NDOSOLVER == 0 )
			Bundle *s = new Bundle( &ParFile );
			//cout << "Bundle criado!" << endl;
			// select the verbosity of log - - - - - - - - - - - - - - - - - - - - -

			int NDOlvl;
			DfltdSfInpt( &ParFile , NDOlvl , int( 0 ) );

		#elif( WHICH_NDOSOLVER == 1 )

			Volume *sv = new Volume( &ParFile );
			Bundle *s = new Bundle( &ParFile );

			// select the verbosity of log - - - - - - - - - - - - - - - - - - - - - -

			int Vollvl, Bndlvl;
			DfltdSfInpt( &ParFile , Vollvl , int( 0 ) );
			DfltdSfInpt( &ParFile , Bndlvl , int( 0 ) );

		#else

			Volume *s = new Volume( &ParFile );

			// select the verbosity of log - - - - - - - - - - - - - - - - - - - - - -

			int NDOlvl;
			DfltdSfInpt( &ParFile , NDOlvl , int( 0 ) );

		#endif

		// select the verbosity of FiOracle log - - - - - - - - - - - - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		int Filvl;
		DfltdSfInpt( &ParFile , Filvl , int( 0 ) );

		//construct the MPSolver object- - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		#if( WHICH_NDOSOLVER <= 1 )

			// select the verbosity of MPlog - - - - - - - - - - - - - - - - - - - - -

			int MPlvl;
			DfltdSfInpt( &ParFile , MPlvl , int( 0 ) );

			#if( WHICH_MPSOLVER == 0)

				QPPenaltyMP *MP = new QPPenaltyMP( &ParFile );

			#elif( WHICH_MPSOLVER == 1 )

				OSIMPSolver *MP = new OSIMPSolver( );
				#if( WHICH_RMPSOLVER == 0 )
					OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
					int algorithm, reduction;
					DfltdSfInpt( &ParFile , algorithm , int(4) );
					DfltdSfInpt( &ParFile , reduction , int(3) );
					SetCplexA( osiSlvr , algorithm , reduction );
				#elif( WHICH_RMPSOLVER == 2 )
					OsiGrbSolverInterface *osiSlvr = new OsiGrbSolverInterface();
					GRBenv * env = osiSlvr->getEnvironmentPtr();
					GRBsetintparam( env, "LogToConsole", 0);
				#else
					cout << "CLP não implementado com estabilização quadratica!" << endl;
				#endif

				MP->SetOsi( osiSlvr );						// pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::quadratic );	// choose the stabilization typedef

			#else

				OSIMPSolver *MP = new OSIMPSolver( );
				
				#if( WHICH_RMPSOLVER == 0 )
					OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
					if (WHICH_DEC >= 4)
					{
						// alterar Markowitz tolerance para sobrepujar o erro no OsiCpxSolverInterface
						CPXENVptr env = osiSlvr->getEnvironmentPtr ();
						CPXsetdblparam( env, CPX_PARAM_EPMRK, 0.1);
					}

					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetCplexA( osiSlvr , algorithm , reduction );
				#elif ( WHICH_RMPSOLVER == 1 )
					OsiClpSolverInterface *osiSlvr = new OsiClpSolverInterface();
					dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->getModelPtr()->setLogLevel(0);
					//dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->messageHandler()->setLogLevel(0);
					//dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->setHintParam(OsiDoReducePrint);

				#else
					OsiGrbSolverInterface *osiSlvr = new OsiGrbSolverInterface();
					GRBenv * env = osiSlvr->getEnvironmentPtr();
					GRBsetintparam( env, "LogToConsole", 0);
					
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetGurobiA( osiSlvr , algorithm , reduction );
				#endif

				MP->SetOsi( osiSlvr );					 // pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::boxstep ); // choose the stabilization typedef

			#endif
			//cout << "MP criado!" << endl;

			// set the verbosity of MPlog - - - - - - - - - - - - - - - - - - - - - -
			ofstream LOGMP( logMP , ofstream::out );
			if( ! LOGMP.is_open() )
				cerr << "Warning: cannot open log file """ << logMP << """" << endl;
			else
				MP->SetMPLog( &LOGMP , MPlvl );

			// activate or not the check for identical subgradients (columns in the Master Problem)
			MP->CheckIdentical(true);

			// pass the MPSolver to the Bundle - - - - - - - - - - - - - - - - - - - -
			s->SetMPSolver( MP );
			MP->SetMPTime();

		#endif  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		// Criar Heuristicas - - - - - - - - - - - - - - - - - - - - - - - - - - 
		#if ( HEURISTIC )
			HpNum kMaxTme_;
			s->GetPar(NDOSolver::kMaxTme, kMaxTme_);
			Heuristicas *heuristica = new Heuristicas( sistema_a, resultadosGurobi, kMaxTme_, &ParFile);
			Fi->SetHeuristic(heuristica);
			//cout << "Heuristica criada!" << endl;
		#endif
		// Definir aqui os ajustes das heuristicas (Fi->SetHeuristic( Heur ) ) em que Heur é um ponteiro para a classe de heuristicas!!!

		// Imprimir tempo subproblemas
		#if ( SUBP_TEMPO )
			Fi->SetSubpTimes();
		#endif

		// set the verbosity of FiOracle's log  - - - - - - - - - - - - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			
		ofstream LOGFile( logF , ofstream::out );
		LOGFile << std::setprecision(15);			// set the precision of log file
		if( ! LOGFile.is_open() )
			cerr << "Warning: cannot open log file """ << logF << """" << endl;
		else
			Fi->SetFiLog( &LOGFile , Filvl );

		// set the verbosity of FiOracle's iteration log  - - - - - - - - - - - - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		string logIF = logF;
		int pos = logIF.length() - 3;
		logIF.erase(pos, 3);
		logIF.append(".it");
		ofstream LOGFileIter( logIF , ofstream::out );
		LOGFileIter << std::setprecision(15);			// set the precision of log file
		if( ! LOGFileIter.is_open() )
			cerr << "Warning: cannot open log file """ << logIF << """" << endl;
		else
			Fi->SetFiIterLog( &LOGFileIter );

		// do timing

		Fi->SetFiTime();
		#if ( HEURISTIC )
			Fi->SetHrstTime();
		#endif

		double tndo;                     // time for NDO
		double torcl;                    // time for Oracle
		double tmp;                      // time for MPSolver
		double tstart;					 // time for smart start
		double theur;					 // time for heuristics
		Index NrIter;                    // number of iterations

		#if( WHICH_NDOSOLVER == 0 )  // - - - - - - - - - - - - - - - - - - - - - -

			// do timing
			s->SetNDOTime();
			//startt = new OPTtimers();

			// pass the FiOracle to the NDOSolver and vice-versa - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			s->SetFiOracle( Fi );
			Fi->SetNDOSolver( s );
			
			// Gerar Lambda inicial e definir lower bound - - - - - - - - - - - - - - - - - - - - - - - - - -

			HpNum FoIn;	
			startt = new OPTtimers();
			startt->Start();
			Index na;
			LMRow Lmbd;
			na = Fi->GetNumVar();
			Lmbd = new LMNum[na];
			ExtDE * model_ext;
			bool flag_cont;
			int flag_rede;
			flag_cont = sistema_a->GetFlagVarBin();
			flag_rede = sistema_a->GetFlagModeloRede();
			sistema_a->SetFlagVarBin( 0 );
			#if ( WARM_START == 2)
				sistema_a->SetFlagBarraUnica( 1 );
			#endif
			model_ext = new ExtDE(sistema_a);
			//cout << "Modelo estendido criado!" << endl;
			FoIn = model_ext->ResolverProblema(Lmbd, resultadosGurobi);
			sistema_a->SetFlagVarBin( flag_cont );
			sistema_a->SetFlagModeloRede( flag_rede );
			delete model_ext;
			Fi->SetFiCont(FoIn);

			#if( WARM_START )
				//// Zerar lambda : mesma coisa que não usar smart start!
				//for (int i = 0; i < na; i++)
				//	Lmbd[i] = 0;
				s->SetLambda(Lmbd);
				
				//Fi->SetLowerBound( - FoIn );		// The continuous solution is not a valid UB on the optimum of the Lagrangian dual!!
				startt->Stop();
			#else
				startt->Stop();
			#endif
			delete Lmbd;

			// Caso de ajuste do Bundle
			//s->SetPar(NDOSolver::kEpsLin, 1e-4);
			HpNum tInit;
			s->GetPar(Bundle::ktMinor, tInit);
			s->SetPar(Bundle::ktInit, tInit);

			HpNum tMax;
			s->GetPar(Bundle::ktStar, tMax);
			s->SetPar(Bundle::ktMaior, tMax*100);		//Let me mention that with t* = 0.01, tMax = 100000000 is clearly ludicrous. You don't expect t to ever go much above t^*, which is "large" already.

			//s->SetPar(Bundle::ktInit, 100.0);
			//s->SetPar(NDOSolver::ktStar, 1.0e+8);
			//s->SetPar(Bundle::ktSPar1, 0);
			//Fi->SetSubpPrecision( 1e-9 );
			
			// defining the value of t*
			//HpNum tstar_calc = 0;
			//if (WHICH_DEC < 4)		// Decomposição espacial
			//{
			//	if (WHICH_MPSOLVER < 2)
			//		tstar_calc = pow((Fi->GetNumVar()), 0.25);
			//	else
			//		tstar_calc = pow((Fi->GetNumVar()), 0.5);
			//}
			//else					// Decomposição por cenários
			//{
			//	if (WHICH_MPSOLVER < 2)
			//		tstar_calc = pow((Fi->GetNumVar()), 0.25)/1;
			//	else
			//		tstar_calc = (std::floor(pow((Fi->GetNumVar()), 0.5)*1e3)/1e3)*100;	//ou dividir por 100	// truncar o numero (com pelo menos 3 casas decimais, pois depende do valor da raiz)
			//		//tstar_calc = pow((Fi->GetNumVar()), 0.5)/100;
			//}
			// dependendo do num de var. Poderia tb considerar FoIn no tstar
			//s->SetPar(NDOSolver::ktStar, tstar_calc);
			///s->SetPar(Bundle::ktInit, tstar_calc/10);
			//s->SetPar(Bundle::ktInit, 0.0001);
			//s->SetPar(NDOSolver::ktStar, 0.1);
			//s->SetPar(Bundle::ktInit, 0.1);

			// set the verbosity of log    - - - - - - - - - - - - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			if( LOGFile.is_open() )
				s->SetNDOLog( &LOGFile , NDOlvl );
			
			// minimize the function - - - - - - - - - - - - - - - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			
			s->KeepBestLambda();
			#if ( WHICH_RMPSOLVER == 1 )
				dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->getModelPtr()->setLogLevel(0);
			#endif
			NDOSolver::NDOStatus Status = s->Solve();
			//s->Solve();
			
			//// Imprimir Lambdas solução
			//Index tam = Fi->GetNumVar();
			//cLMRow Ls = Fi->GetLambda();
			////for (int j = 0; j < tam; j++)
			////	cout << Ls[j] << endl;
			//
			//ofstream * inFile;
			//inFile = new ofstream( "LM.txt", ios::out );
			//if ( inFile->is_open() )
			//{
			//	for (int j = 0; j < tam; j++)
			//		*inFile << Ls[j] << endl;
			//	inFile->close();
			//}
			//else
			//	cout << "Unable to open file";
			///resultadosGurobi->ImprimirArmazenarXtilhat();

			// get the running time - - - - - - - - - - - - - - - - - - - - - - - - - -

			tndo = s->NDOTime();						// time for NDO
			torcl = Fi->FiTime();						// time for Oracle
			tmp = MP->MPTime();							// time for MPSolver
			tstart = ( startt ? startt->Read() : 0 );	// time for smart start
			theur = Fi->HrstTime();						// time for heuristics
			NrIter = s->NrIter();						// number of iterations
			
			// save and print the heuristic primal results
			#if ( HEURISTIC )
				resultadosGurobi->GravarSolHeuristica(heuristica->GetFO(), heuristica->GetX(), - s->ReadBestFiVal());
				//string saida = "Hrst_" + arv + "_" + IDpar + ".txt";
				string saida = "log_" + arv + "_" + real + "_" + IDpar + ".hr";
				resultadosGurobi->EscreverArquivo(saida, 0, tndo, "");
			#endif

		#elif( WHICH_NDOSOLVER == 1 )
			// implementar (Bundle + Volume)
		#else
			// implementar (Volume)
		#endif

		HpNum OV = - s->ReadBestFiVal();            // get the optimal value
		HpNum OpV = Inf<HpNum>();
		#if ( HEURISTIC )
			OpV = Fi->GetHeuristicSolution();		// get the optimal primal value
		#endif

		// print results- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		if( s->IsOptimal() )
			cout << " (optimal value)" << endl;
		else
			cout << " (not provably optimal)" << endl;
		cout << "Bundle status: " << Status << endl;

		cout << endl << "Iter = " << s->NrIter() << ", time = " << s->NDOTime()	<< endl;

		// deallocate everything- - - - - - - - - - - - - - - - - - - - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		#if WHICH_NDOSOLVER <= 1

			s->SetMPSolver();
			delete MP;

			#if( MPSOLVER )
				delete osiSolver;
			#endif

		#endif
		
		// Imprimir tempo subproblemas
		#if ( SUBP_TEMPO )
			vetorfloat timeSbp = Fi->SbPTime();
		for (size_t nsbp = 0; nsbp < timeSbp.size(); nsbp++)
			LOGFile << "Subproblem " << nsbp << " time (s): " << timeSbp[nsbp] << "\n";
		#endif

		//s->SetFiOracle( NULL );
		delete s;
		#if ( HEURISTIC )
			delete heuristica;
		#endif
		delete Fi;
				
		#if( WARM_START )
			LOGFile << "Heuristic time (s): " << theur << "\n" << "SmartStart time (s): " << tstart << "\n" << "Cont. Relaxation: " << FoIn;
			LOGFile.close();
		#else
			LOGFile << "Heuristic time (s): " << theur << "\n" << "Cont. Relaxation: (not used as lambda0) " << FoIn;;
			LOGFile.close();
		#endif

		


		LOGFileIter.close();

		// output the results- - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		#if( WARM_START )
			Out << setprecision(10) << logF << " : " << "Cont. Relaxation:\t" << FoIn << "\t" << "Time_NDO:\t" << tndo << "\t" << "Time_SmartStart:\t" << tstart << "\t" << "Time_Fi:\t" << torcl << "\t" << "Time_Hrtcs:\t" << theur << "\t" <<  "Time_MP:\t" << tmp << "\t" << "Niter:\t" << NrIter << "\t" << "Optimal Heuristic:\t" << OpV << "\t";
			//Out << logF << "\t" << "Cont. Relaxation: " << "\t" << FoIn << "\t" << "Time_NDO: " << "\t" << tndo << "\t" << "Time_SmartStart: " << "\t" << tstart << "\t" << "Time_Fi: " << "\t" << torcl << "\t" << "Time_Hrtcs: " << "\t" << theur << "\t" <<  "Time_MP: " << "\t" << tmp << "\t" << "Niter: " << "\t" << NrIter << "\t" << "Optimal Heuristic: " << "\t" << OpV << "\t";
		#else
			Out << setprecision(10) << logF << " : " << "Time_NDO:\t" << tndo << "\t" << "Time_SmartStart:\t" << tstart << "\t" << "Time_Fi:\t" << torcl << "\t" << "Time_Hrtcs:\t" << theur << "\t" <<  "Time_MP:\t" << tmp << "\t" << "Niter:\t" << NrIter << "\t" << "Optimal Heuristic:\t" << OpV << "\t";
		#endif

		switch( Status )
		{
		case( NDOSolver::kOK ) :
			Out.precision( 16 );
			Out << "Status: OK, Value:\t" << OV << endl;
			break;
		case( NDOSolver::kStopped ) :
			Out.precision( 16 );
			Out << "Status: Stopped, Value:\t" << OV << endl;
			break;
		case( NDOSolver::kUnfsbl ) :
			Out << "Status: Unfeas." << endl;
			break;
		case( NDOSolver::kUnbndd ) :
			Out << "Status: Unbound." << endl;
			break;
		case( NDOSolver::kStpIter ) :
			Out.precision( 16 );
			Out << "Status: kStpIter, Value:\t" << OV << endl;
			break;
		case( NDOSolver::kStpTime ) :
			Out.precision( 16 );
			Out << "Status: kStpTime, Value:\t" << OV << endl;
			break;
		default :
			Out << "Status: Error" << endl;
		}

		Out.close();

	}  // end( try-block )

	// managing exceptions - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	catch( exception &e )
	{
		cerr << e.what() << endl;
		return( 1 );
	}
	catch(...)
	{
		cerr << "Error: unknown exception thrown" << endl;
		return( 1 );
	}

	// the end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	return( 0 );

};  // end( ResolverRL )

int Problema_Dec::ResolverHrst()
{
	// deixar essa função como resolver heuristica sozinha!!
	// para qualquer ajsute de alfa, beta, gama...
	// imprimir log de relatório...
	#if ( HEURISTIC )
	{
		ofstream LOG( "Hrst_result.it" , ofstream::out );
		// consideram-se os ajustes de parametros padrao do constructor da heuristica
		Heuristicas *heuristica = new Heuristicas( sistema_a, resultadosGurobi, 10);
		LOG << setprecision(15);
		heuristica->SetLog(&LOG);

		double resultado = heuristica->ResolverHeuristica();
		cout << setprecision(15) << "Heur. = " << resultado <<endl;
		resultadosGurobi->GravarSolHeuristica(resultado, heuristica->GetX(), resultado);
		resultadosGurobi->EscreverArquivo("Hrst_result.hr", 0, 0, "C:/Users/Murilo/Documents/Visual Studio 2010/Projects/PDO_2_1/PDO_2/");
	
		delete heuristica;
	}
	#else
		cout << "Definir heuristica!" << endl;
	#endif	

	return( 0 );

}

int Problema_Dec::Teste()		// testar resoluçao dos subproblemas (depois que estiver td ok alterar o Oracle...)
{
	////-------------------------------------------------------------------------
	//// Testar resoluçao de cada subproblema separadamente
	//GRBEnv ambGRB;
	//// criar lambda aleatório para testar os resultados...
	//vetorfloat lambda, x, xa;
	//double fdual;
	//int NumVar = (sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())) * (sistema_a->termeletricasVtr.size() + (4+flag3)*sistema_a->hidreletricasVtr.size());
	//lambda.resize(NumVar);
	//
	//srand( ii );
	//for (int i = 0; i < NumVar; i++)
	//	lambda[i] = double((rand() % 100) - (rand() % 100)) / pow(10.0,(rand() % 5));

	//// Tempo
	//LARGE_INTEGER clockPerSec;
	//LARGE_INTEGER inicio, fim;
	//QueryPerformanceFrequency(&clockPerSec);

	//CResultados * resultadosGurobi0;
	//resultadosGurobi0 = new CResultados(sistema_a);

	//// Subproblema Hidreletrico
	//CSubProblemaHE * subproblemaHE;
	//subproblemaHE = new CSubProblemaHE(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//subproblemaHE->ResolverProblemaRL(resultadosGurobi0, &lambda, 0);	
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//cout << setprecision(12) << - resultadosGurobi0->GetFobj() << endl;
	//
	//ConjSPHE * conjspHE;
	//conjspHE = new ConjSPHE(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//for (int n = 0; n < (sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())); n++)
	//	for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
	//		conjspHE->ResolverProblemaRL(resultadosGurobi, &lambda, 0, r, n);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//cout << setprecision(12) << - resultadosGurobi->GetFobj() << endl << endl;
	
	//// Subproblema Hidraulico
	//CSubProblemaHA * subproblemaHA;
	//subproblemaHA = new CSubProblemaHA(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//subproblemaHA->ResolverProblemaRL(resultadosGurobi0, &lambda, 0);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//fdual = - resultadosGurobi0->GetFobj();
	//cout << - resultadosGurobi0->GetFobj() << endl;

	//ConjSPHA * conjspHA;
	//conjspHA = new ConjSPHA(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//for (int n = 0; n < conjspHA->GetNCascatas(); n++)
	//	conjspHA->ResolverProblemaRL(resultadosGurobi, &lambda, 0, n);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//cout << - resultadosGurobi->GetFobj() << endl << endl;

	//// Subproblema Demanda
	//CSubProblemaD * subproblemaD;
	//subproblemaD = new CSubProblemaD(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//subproblemaD->ResolverProblemaRL(resultadosGurobi0, &lambda, 0);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//fdual = - resultadosGurobi0->GetFobj();
	//cout << - resultadosGurobi0->GetFobj() << endl;

	//ConjSPD * conjspD;
	//conjspD = new ConjSPD(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//for (int n = 0; n < (sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())); n++)
	//	conjspD->ResolverProblemaRL(resultadosGurobi, &lambda, 0, n);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//fdual = - resultadosGurobi->GetFobj();
	//cout << - resultadosGurobi->GetFobj() << endl << endl;

	//// Subproblema Termelétrico
	//CSubProblemaT * subproblemaT;
	//subproblemaT = new CSubProblemaT(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//subproblemaT->ResolverProblemaRL(resultadosGurobi0, &lambda, 0);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//fdual = - resultadosGurobi0->GetFobj();
	//cout << - resultadosGurobi0->GetFobj() << endl;

	//ConjSPT * conjspT;
	//conjspT = new ConjSPT(sistema_a, ambGRB);
	//QueryPerformanceCounter(&inicio);
	//for (int n = 0; n < sistema_a->termeletricasVtr.size(); n++)
	//	conjspT->ResolverProblemaRL(resultadosGurobi, &lambda, 0, n);
	//QueryPerformanceCounter(&fim);
	//cout << double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart) << "  :  ";
	//fdual = - resultadosGurobi->GetFobj();
	//cout << - resultadosGurobi->GetFobj() << endl;

	//resultadosGurobi->ExportarX("x.txt");
	//resultadosGurobi->ExportarXA("xa.txt");
	//resultadosGurobi0->ExportarX("x0.txt");
	//resultadosGurobi0->ExportarXA("xa0.txt");

	//SgRow SubG;
	//SubG = new double[resultadosGurobi0->GetNA()];
	//resultadosGurobi0->GetSubGradiente(5, SubG);

	//vetorfloat subgradiente;
	//subgradiente.resize(resultadosGurobi0->GetNA());

	//for (int i = 0; i < resultadosGurobi0->GetNA(); i++)
	//	subgradiente[i] = SubG[i];

	
	////-------------------------------------------------------------------------
	// Teste com a Heuristica

	//// Carregar soluçao do ED
	//vetorfloat x1, x2;
	//x1.resize(resultadosGurobi->Getn());
	//x2.resize(resultadosGurobi->Getn());
	////ifstream inFile( "x_opt.txt", ios::in );
	//ifstream inFile( "x_hat.txt", ios::in ); 
	//if ( !inFile )                                                            
	//{                                                                               
	//	cout << "File x.txt could not be opened" << endl;
	//	exit( 1 );
	//}
	//while ( ! inFile.eof() )
	//{
	//	for (size_t i = 0; i < x1.size(); i++)
	//		inFile >> x1[i];
	//}
	//inFile.close();

	////ifstream inFile2( "x_opt.txt", ios::in );
	//ifstream inFile2( "x_til.txt", ios::in ); 
	//if ( !inFile2 )                                                            
	//{                                                                               
	//	cout << "File x.txt could not be opened" << endl;
	//	exit( 1 );
	//}
	//while ( ! inFile2.eof() )
	//{
	//	for (size_t i = 0; i < x2.size(); i++)
	//		inFile2 >> x2[i];
	//}
	//inFile2.close();


	// deixar essa função como resolver heuristica sozinha!!
	// para qualquer ajsute de alfa, beta, gama...
	// imprimir log de relatório...

	//Heuristicas *heuristica = new Heuristicas( sistema_a, resultadosGurobi, 600 );
	//cout << setprecision(15) << "Heur. = " << heuristica->ResolverHeuristica(350) <<endl;

	return( 0 );
}
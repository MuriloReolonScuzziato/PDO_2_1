/*--------------------------------------------------------------------------*/
/*---------------------- File Problema_Dec_Teste.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para resolver o problema por decomposição
 * Para varios parametros (t^* e tinit) e OSIsolvers diferentes
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

#include "Problema_Dec_Teste.h"

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
		#define slvr_I 0
		#define slvr_F 1
	#elif WHICH_MPSOLVER == 1
		#include "OSIMPSolver.h"
		#include "OsiCpxSolverInterface.hpp"
		#define slvr_I 0
		#define slvr_F 1
	#elif WHICH_MPSOLVER == 2
		#include "OSIMPSolver.h"
		#if WHICH_RMPSOLVER == 0
			#include "OsiCpxSolverInterface.hpp"
			typedef OsiCpxSolverInterface RealSolverInterface;
			#define slvr_I 0
			#define slvr_F 1
		#elif WHICH_RMPSOLVER == 1
			#include "OsiClpSolverInterface.hpp"
			typedef OsiClpSolverInterface RealSolverInterface;
			#define slvr_I 1
			#define slvr_F 2
		#elif WHICH_MPSOLVER == 2
			#include "OsiGrbSolverInterface.hpp"
			typedef OsiGrbSolverInterface RealSolverInterface;
			#define slvr_I 2
			#define slvr_F 3
		#else WHICH_MPSOLVER == 3
			#include "OSIMPSolver.h"
			#include "OsiCpxSolverInterface.hpp"
			#include "OsiClpSolverInterface.hpp"
			#include "OsiGrbSolverInterface.hpp"
			#define slvr_I 0
			#define slvr_F 3
		#endif
	#else
		#error "Unable to find the MPsolver."
	#endif
#endif

#include <fstream>
#include <sstream>

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

void SetCplexAA( OsiCpxSolverInterface * osiCpx , const int algorithm = 0 ,
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

void SetGrbAA( OsiGrbSolverInterface * osiGrb , const int algorithm = 0 ,
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

Problema_Dec_Teste::Problema_Dec_Teste(CSistema * const sistema_end , int tipo)
{
	// Fazer essa classe generica (interface para a decomposiçao espacial e por cenários) -> usar typedef
	sistema_a = sistema_end;
	resultadosGurobi = new Resultados(sistema_a);
	tipo_param = tipo;
	startt = 0;
}
Problema_Dec_Teste::~Problema_Dec_Teste(void)
{
	delete startt;
}

int Problema_Dec_Teste::ResolverRL(string arv, string real)
{
	if (tipo_param == 1)		// Uma bateria de teste com parametros e solvers (do bundle) diferentes para um mesmo problema
		return ResolverRL1(arv, real);
	if (tipo_param == 2)		// Uma bateria de teste com parametros da heuristica diferentes para um mesmo problema
		return ResolverRL2(arv, real);
	return ( 1 );
}
int Problema_Dec_Teste::ResolverRL1(string arv, string real)
{
	// Observações sobre o t* adjustment test:
	// Para a decomposição spacial...
	// Para a decomposição por cenários: t* máximo de 10/100, alterar Mipgap do primeiro teste (além do Eps), t^min deve ser de 1e-2/1e-3
	//

	// Note that increasing t is something that the Bundle does much more rapidly and happily than decreasing it


	/// se n converge no tempo limite para 1e8 reduzir para 1e7 e tentar, até ter uma solução ótima sem esgotar o tempo limite, só depois entçao iniciar de valores baixos para encontrar o t* ideal
	/// tempo máximo? 300000

	// Log files
	string log = "log_" + arv + "_" + real;
	//const char *const logF = log.c_str();
	string MPlog = "MPlog_" + arv + "_" + real;
	//const char *const logMP = MPlog.c_str();
	size_t log_tam = log.size();
	size_t MPlog_tam = MPlog.size();

	// Define valores de t^* (tinit = t^* e t^* / 10)
	vetorfloat t_s;
	t_s.push_back( 1e7 );		// t* para obter o valor de referencia (função dual ótima)
	for (int i = 0; i <= 10; i++)		// [0 11]
		t_s.push_back( pow(10, double(i - 4)) );
	
	ostringstream convert;
	// Loop para resolver o problema com um ou mais solvers e diferentes parametros (t^* e tinit)
	int stop_cond = 0;
	for (int solver = slvr_I; solver < slvr_F; solver++)
	{
		for (size_t i = 0; i < t_s.size(); i++)
		{
			for (int ii = 1; ii < 2; ii++)
			{
				#if (WHICH_MPSOLVER == 0)
					log.append("_QPP");
					MPlog.append("_QPP");
				#elif (WHICH_MPSOLVER == 1)
					log.append("_CPLEXQ");
					MPlog.append("_CPLEXQ");
				#else
				if (solver == 0)
				{
					log.append("_CPLEXL");
					MPlog.append("_CPLEXL");
				}
				else if (solver == 1)
				{
					log.append("_CLP");
					MPlog.append("_CLP");
				}
				else if (solver == 2)
				{
					log.append("_GRP");
					MPlog.append("_GRP");
				}
				#endif
				convert << "_" << t_s[i] << "_" << t_s[i] / pow(10, double(ii)) << ".bn";
				log.append(convert.str());
				MPlog.append(convert.str());
				convert.str("");

				const char *const logF = log.c_str();
				const char *const logMP = MPlog.c_str();
				resultadosGurobi->ZerarSolucoes();
				
				stop_cond = ChamarBundle1(logF, logMP, t_s[i], t_s[i] / pow(10, double(ii)), solver, i);
				if (stop_cond == 2)
					break;		// um bom valor para o t* já foi atingido!
				log.erase(log_tam, log.size());
				MPlog.erase(MPlog_tam, MPlog.size());
			}
			if (stop_cond == 2)
				break;
		}
		if (stop_cond == 2)
			break;
	}

	return ( 0 );
}
int Problema_Dec_Teste::ResolverRL2(string arv, string real)
{
	// Log files
	string log = "log_" + arv + "_" + real;
	//const char *const logF = log.c_str();
	string MPlog = "MPlog_" + arv + "_" + real;
	//const char *const logMP = MPlog.c_str();
	size_t log_tam = log.size();
	size_t MPlog_tam = MPlog.size();

	// Define os casos
	int nc = 11;		// number of cases (different parameters adjusts)
	int nw = 12;		// number of options for the window size
	int BVMW[] = {0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2};
	int BVFW[] = {0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2};
	double alfa[] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0};
	double beta[] = {0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5};
	double gama[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.5};
	int MWS[] = {1, 1, 1, 1, 6, 6, 6, 6, 12, 12, 12, 24};
	int FWS[] = {1, 6, 12, 24, 1, 6, 12, 24, 1, 6, 12, 0};

	ostringstream convert;
	// Loop para resolver usar diferentes tamanhos de janelas nas heuristicas
	for (int i = 0; i < nw; i++)
		for (int ii = 0; ii < nc; ii++)		// Loop para usar diferentes parametros nas heuristicas
		{
			#if (WHICH_MPSOLVER == 0)
				log.append("_QPP");
				MPlog.append("_QPP");
			#elif (WHICH_MPSOLVER == 1)
				log.append("_CPLEXQ");
				MPlog.append("_CPLEXQ");
			#else
				#if( WHICH_RMPSOLVER == 0 )
					log.append("_CPLEXL");
					MPlog.append("_CPLEXL");
				#elif ( WHICH_RMPSOLVER == 1 )
					log.append("_CLP");
					MPlog.append("_CLP");
				#elif ( WHICH_RMPSOLVER == 2 )
					log.append("_GRB");
					MPlog.append("_GRB");
				#endif
			#endif
			
			convert << "_case" << ii + 1 << "_" << MWS[i] << "_" << FWS[i] << ".bn";
			log.append(convert.str());
			MPlog.append(".bn");
			//MPlog.append(convert.str());
			convert.str("");

			const char *const logF = log.c_str();
			const char *const logMP = MPlog.c_str();
			resultadosGurobi->ZerarSolucoes();
			ChamarBundle2(logF, logMP, MWS[i], FWS[i], BVMW[ii], BVFW[ii], alfa[ii], beta[ii], gama[ii]);
			log.erase(log_tam, log.size());
			MPlog.erase(MPlog_tam, MPlog.size());
		}
	return ( 0 );
}

int Problema_Dec_Teste::ChamarBundle1(const char *const logF, const char *const logMP, const HpNum tstar_, const HpNum tinit_, const int solver, int iter)
{
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
		#elif( WHICH_DEC == 0 )
			EspFiOracle *Fi = new EspFiOracle( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#endif

		// Criar modelos do oracle - - - - - - - - - - - - - - - - - - - - - - - - - - 

		//Fi->SetData(cntr , UB);	// Definir limites para os Lambdas e the minimum of the function (L0)
		Fi->CriarModelos();
		
		// open parameters file - - - - - - - - - - - - - - - - - - - - - - - - - -

		ifstream ParFile( ParF );
		if( ! ParFile.is_open() )
			cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;

		// construct the NDOSolver object - - - - - - - - - - - - - - - - - - - - -

		#if( WHICH_NDOSOLVER == 0 )
			Bundle *s = new Bundle( &ParFile );
			
			// Define novos valores para os parametros t^* e tinit
			s->SetPar(NDOSolver::ktStar, tstar_);
			s->SetPar(Bundle::ktInit, tinit_);
			if (iter == 0)		// primeira iteração para obter o valor de referencia
			{
				s->SetPar(NDOSolver::kEpsLin, 1e-6);
				s->SetPar(Bundle::ktInit, 1.0);
				s->SetPar(Bundle::ktSPar1, 0);
				Fi->SetSubpPrecision( 1e-9 );
			}
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

		// Criar Heuristicas - - - - - - - - - - - - - - - - - - - - - - - - - - 
		#if ( HEURISTIC )
			HpNum kMaxTme_;
			s->GetPar(NDOSolver::kMaxTme, kMaxTme_);
			Heuristicas *heuristica = new Heuristicas( sistema_a, resultadosGurobi, kMaxTme_ );
			Fi->SetHeuristic(heuristica);
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

				OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
				int algorithm, reduction;
				DfltdSfInpt( &ParFile , algorithm , int(4) );
				DfltdSfInpt( &ParFile , reduction , int(3) );
				SetCplexA( osiSlvr , algorithm , reduction );

				MP->SetOsi( osiSlvr );					 // pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::quadratic ); // choose the stabilization typedef

			#elif ( WHICH_MPSOLVER == 2 )

				OSIMPSolver *MP = new OSIMPSolver( );

				#if( WHICH_RMPSOLVER == 0 )
					OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetCplexA( osiSlvr , algorithm , reduction );
				#elif ( WHICH_RMPSOLVER == 1 )
					OsiClpSolverInterface *osiSlvr = new OsiClpSolverInterface();
					dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->getModelPtr()->setLogLevel(0);
				#else
					OsiGrbSolverInterface *osiSlvr = new OsiGrbSolverInterface();
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetGurobiA( osiSlvr , algorithm , reduction );
				#endif

				MP->SetOsi( osiSlvr );					 // pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::boxstep ); // choose the stabilization typedef
			#else

				OSIMPSolver *MP = new OSIMPSolver( );

				if ( solver == 0 )
				{
					OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetCplexAA( osiSlvr , algorithm , reduction );
				}
				else if ( solver == 1 )
				{
					OsiClpSolverInterface *osiSlvr = new OsiClpSolverInterface();
					dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->getModelPtr()->setLogLevel(0);
				}
				else if ( solver == 2 )
				{
					OsiGrbSolverInterface *osiSlvr = new OsiGrbSolverInterface();
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetGurobiAA( osiSlvr , algorithm , reduction );
				}

				MP->SetOsi( osiSlvr );					 // pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::boxstep ); // choose the stabilization typedef

			#endif

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

		#else
			#error "Codigo implementado apenas para o NDOsolver Bundle."
		#endif  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		//// select main things in the FiOracle - - - - - - - - - - - - - - - - - - -
		////- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		//bool yie , agg;

		//Fi->GetPar( FlwFiOrcl::kYiE , yie );    // yie: easy components
		//Fi->GetPar( FlwFiOrcl::kAgg , agg );    // agg: aggregation

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
		logIF.append("_iter.txt");
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


		HpNum OV;
		#if( WHICH_NDOSOLVER == 0 )  // - - - - - - - - - - - - - - - - - - - - - -

			// do timing

			s->SetNDOTime();

			// pass the FiOracle to the NDOSolver and vice-versa - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			s->SetFiOracle( Fi );
			Fi->SetNDOSolver( s );

			// Gerar Lambda inicial e definir lower bound - - - - - - - - - - - - - - - - - - - - - - - - - - 
			#if( WARM_START )
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
				FoIn = model_ext->ResolverProblema(Lmbd);
				sistema_a->SetFlagVarBin( flag_cont );
				sistema_a->SetFlagModeloRede( flag_rede );
				delete model_ext;
				s->SetLambda(Lmbd);
				delete Lmbd;

				//Fi->SetLowerBound( - FoIn );		// The continuous solution is not a valid UB on the optimum of the Lagrangian dual!!
				startt->Stop();
			#endif
				
			// set the verbosity of log    - - - - - - - - - - - - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			if( LOGFile.is_open() )
				s->SetNDOLog( &LOGFile , NDOlvl );

			// minimize the function - - - - - - - - - - - - - - - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			NDOSolver::NDOStatus Status;
			if ((iter > 0) || (sistema_a->GetFlagVarBin() == 1))	// problema inteiro
			{
				s->KeepBestLambda();
				Status = s->Solve();
				//s->Solve();
			}
			else
			{
				Index na;
				LMRow Lmbd;
				na = Fi->GetNumVar();
				Lmbd = new LMNum[na];
				ExtDE * model_ext;
				//bool flag_cont;
				//int flag_rede;
				model_ext = new ExtDE(sistema_a);
				OV = model_ext->ResolverProblema(Lmbd);
				Status = NDOSolver::kOK;
			}

			// get the running time - - - - - - - - - - - - - - - - - - - - - - - - - -

			tndo = s->NDOTime();						// time for NDO
			torcl = Fi->FiTime();						// time for Oracle
			tmp = MP->MPTime();							// time for MPSolver
			tstart = ( startt ? startt->Read() : 0 );	// time for smart start
			theur = Fi->HrstTime();						// time for heuristics
			NrIter = s->NrIter();						// number of iterations
			

		#elif( WHICH_NDOSOLVER == 1 )
			// implementar (Bundle + Volume)
		#else
			// implementar (Volume)
		#endif

		if ((iter > 0) || (sistema_a->GetFlagVarBin() == 1))	// problema inteiro
			OV = - s->ReadBestFiVal();            // get the optimal value
		HpNum OpV = Inf<HpNum>();
		#if ( HEURISTIC )
			OpV = Fi->GetHeuristicSolution();		// get the optimal primal value
		#endif

		// print results- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		if( s->IsOptimal() )
			cout << " (optimal value)" << endl;
		else
			cout << " (not provably optimal)" << endl;

		cout << endl << "Iter = " << s->NrIter() << ", time = " << s->NDOTime()	<< endl;

		// test the conditions of stopping (reference value)- - - - - - - - -
		HpNum eU;
		if (iter == 0)
			refere_value = OV;
		else
			s->GetPar(NDOSolver::kEpsLin, eU);

		// deallocate everything- - - - - - - - - - - - - - - - - - - - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		#if WHICH_NDOSOLVER <= 1

			s->SetMPSolver();
			delete MP;

			#if( MPSOLVER )
				delete osiSolver;
			#endif

		#endif

		delete s;
		#if ( HEURISTIC )
			delete heuristica;
		#endif
		delete Fi;

		#if( WARM_START )
			LOGFile << "Heuristic time (s): " << theur << "\n" << "SmartStart time (s): " << tstart << "\n" << "Cont. Relaxation: " << FoIn;
			LOGFile.close();
		#endif
  
		LOGFileIter.close();

		// output the results- - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		#if( WARM_START )
			Out << logF << " : " << "Cont. Relaxation: " << FoIn << "\t" << "Time_NDO: " << tndo << "\t" << "Time_SmartStart: " << tstart << "\t" << "Time_Fi: " << torcl << "\t" << "Time_Hrtcs: " << theur << "\t" <<  "Time_MP: " << tmp << "\t" << "Niter: " << NrIter << "\t" << "Optimal Heuristic: " << OpV << "\t";
			//Out << logF << "\t" << "Cont. Relaxation: " << "\t" << FoIn << "\t" << "Time_NDO: " << "\t" << tndo << "\t" << "Time_SmartStart: " << "\t" << tstart << "\t" << "Time_Fi: " << "\t" << torcl << "\t" << "Time_Hrtcs: " << "\t" << theur << "\t" <<  "Time_MP: " << "\t" << tmp << "\t" << "Niter: " << "\t" << NrIter << "\t" << "Optimal Heuristic: " << "\t" << OpV << "\t";
		#else
			Out << logF << " : " << "Time_NDO: " << tndo << "\t" << "Time_SmartStart: " << tstart << "\t" << "Time_Fi: " << torcl << "\t" << "Time_Hrtcs: " << theur << "\t" <<  "Time_MP: " << tmp << "\t" << "Niter: " << NrIter << "\t" << "Optimal Heuristic: " << OpV << "\t";
		#endif

		switch( Status )
		{
		case( NDOSolver::kOK ) :
			Out.precision( 16 );
			Out << "Status: OK, Value: " << OV << endl;
			break;
		case( NDOSolver::kStopped ) :
			Out.precision( 16 );
			Out << "Status: Stopped, Value: " << OV << endl;
			break;
		case( NDOSolver::kUnfsbl ) :
			Out << "Status: Unfeas." << endl;
			break;
		case( NDOSolver::kUnbndd ) :
			Out << "Status: Unbound." << endl;
			break;
		case( NDOSolver::kStpIter ) :
			Out.precision( 16 );
			Out << "Status: kStpIter, Value: " << OV << endl;
			break;
		case( NDOSolver::kStpTime ) :
			Out.precision( 16 );
			Out << "Status: kStpTime, Value: " << OV << endl;
			break;
		default :
			Out << "Status: Error" << endl;
		}

		Out.close();

		// test the conditions of stopping (reference value)- - - - - - - - -
		if (iter > 0)
		{
			double real_error;
			real_error = (refere_value - OV) / OV;
			if (real_error <= eU)
			{
				cout << "Valor de t* determinado! t* = " << tstar_ << endl;
				return ( 2 );		// condição para indicar que o for seja parado
			}
		}

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
int Problema_Dec_Teste::ChamarBundle2(const char *const logF, const char *const logMP, int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a )
{
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
		#elif( WHICH_DEC == 0 )
			EspFiOracle *Fi = new EspFiOracle( sistema_a, resultadosGurobi, MODEL_CHRCTR);
		#endif

		// Criar modelos do oracle - - - - - - - - - - - - - - - - - - - - - - - - - - 

		//Fi->SetData(cntr , UB);	// Definir limites para os Lambdas e the minimum of the function (L0)
		Fi->CriarModelos();
			
		// open parameters file - - - - - - - - - - - - - - - - - - - - - - - - - -

		ifstream ParFile( ParF );
		if( ! ParFile.is_open() )
			cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;

		// construct the NDOSolver object - - - - - - - - - - - - - - - - - - - - -

		#if( WHICH_NDOSOLVER == 0 )
			Bundle *s = new Bundle( &ParFile );

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

		// Criar Heuristicas - - - - - - - - - - - - - - - - - - - - - - - - - - 
		#if ( HEURISTIC )
			HpNum kMaxTme_;
			s->GetPar(NDOSolver::kMaxTme, kMaxTme_);
			Heuristicas *heuristica = new Heuristicas( sistema_a, resultadosGurobi, ws_a, wlas_a, bvmw, bvfw, alfa_a, beta_a, gama_a, kMaxTme_ );
			Fi->SetHeuristic(heuristica);
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

				OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
				int algorithm, reduction;
				DfltdSfInpt( &ParFile , algorithm , int(4) );
				DfltdSfInpt( &ParFile , reduction , int(3) );
				SetCplexA( osiSlvr , algorithm , reduction );

				MP->SetOsi( osiSlvr );					 // pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::quadratic ); // choose the stabilization typedef

			#else

				OSIMPSolver *MP = new OSIMPSolver( );

				#if( WHICH_RMPSOLVER == 0 )
					OsiCpxSolverInterface *osiSlvr = new OsiCpxSolverInterface();
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetCplexA( osiSlvr , algorithm , reduction );
				#elif ( WHICH_RMPSOLVER == 1 )
					OsiClpSolverInterface *osiSlvr = new OsiClpSolverInterface();
					dynamic_cast<OsiClpSolverInterface*>( osiSlvr )->getModelPtr()->setLogLevel(0);
				#else
					OsiGrbSolverInterface *osiSlvr = new OsiGrbSolverInterface();
					//int algorithm, reduction;
					//DfltdSfInpt( &ParFile , algorithm , int(0) );
					//DfltdSfInpt( &ParFile , reduction , int(3) );
					//SetGurobiA( osiSlvr , algorithm , reduction );
				#endif

				MP->SetOsi( osiSlvr );					 // pass the RealSolver to MPSolver
				MP->SetStabType( OSIMPSolver::boxstep ); // choose the stabilization typedef

			#endif

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

		//// select main things in the FiOracle - - - - - - - - - - - - - - - - - - -
		////- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		//bool yie , agg;

		//Fi->GetPar( FlwFiOrcl::kYiE , yie );    // yie: easy components
		//Fi->GetPar( FlwFiOrcl::kAgg , agg );    // agg: aggregation

		// set the verbosity of FiOracle's log  - - - - - - - - - - - - - - - - - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			
		ofstream LOGFile( logF , ofstream::out );
		LOGFile << std::setprecision(15);			// set the precision of log file
		if( ! LOGFile.is_open() )
			cerr << "Warning: cannot open log file """ << logF << """" << endl;
		else
			Fi->SetFiLog( &LOGFile , Filvl );

		// do timing

		Fi->SetFiTime();
		Fi->SetHrstTime();

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

			// set the verbosity of log    - - - - - - - - - - - - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			if( LOGFile.is_open() )
				s->SetNDOLog( &LOGFile , NDOlvl );


			// Gerar Lambda inicial e definir lower bound - - - - - - - - - - - - - - - - - - - - - - - - - -
			#if( WARM_START )
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
				FoIn = model_ext->ResolverProblema(Lmbd);
				sistema_a->SetFlagVarBin( flag_cont );
				sistema_a->SetFlagModeloRede( flag_rede );
				delete model_ext;
				s->SetLambda(Lmbd);
				delete Lmbd;

				//Fi->SetLowerBound( - FoIn );		// The continuous solution is not a valid UB on the optimum of the Lagrangian dual!!
				startt->Stop();
			#endif

			// minimize the function - - - - - - - - - - - - - - - - - - - - - - - - -
			// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			s->KeepBestLambda();
			NDOSolver::NDOStatus Status = s->Solve();
			//s->Solve();

			// get the running time - - - - - - - - - - - - - - - - - - - - - - - - - -

			tndo = s->NDOTime();						// time for NDO
			torcl = Fi->FiTime();						// time for Oracle
			tmp = MP->MPTime();							// time for MPSolver
			tstart = ( startt ? startt->Read() : 0 );	// time for smart start
			theur = Fi->HrstTime();						// time for heuristics
			NrIter = s->NrIter();						// number of iterations
			

		#elif( WHICH_NDOSOLVER == 1 )
			// implementar (Bundle + Volume)
		#else
			// implementar (Volume)
		#endif

		HpNum OV = - s->ReadBestFiVal();            // get the optimal value
		HpNum OpV = Fi->GetHeuristicSolution();		// get the optimal primal value

		// print results- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		if( s->IsOptimal() )
			cout << " (optimal value)" << endl;
		else
			cout << " (not provably optimal)" << endl;

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

		delete s;
		#if ( HEURISTIC )
			delete heuristica;
		#endif
		delete Fi;

		#if( WARM_START )
			LOGFile << "Heuristic time (s): " << theur << "\n" << "SmartStart time (s): " << tstart << "\n" << "Cont. Relaxation: " << FoIn;
			LOGFile.close();
		#endif
  
		// output the results- - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		#if( WARM_START )
			Out << logF << " : " << "Cont. Relaxation: " << FoIn << "\t" << "Time_NDO: " << tndo << "\t" << "Time_SmartStart: " << tstart << "\t" << "Time_Fi: " << torcl << "\t" << "Time_Hrtcs: " << theur << "\t" <<  "Time_MP: " << tmp << "\t" << "Niter: " << NrIter << "\t" << "Optimal Heuristic: " << OpV << "\t";
		#else
			Out << logF << " : " << "Time_NDO: " << tndo << "\t" << "Time_SmartStart: " << tstart << "\t" << "Time_Fi: " << torcl << "\t" << "Time_Hrtcs: " << theur << "\t" <<  "Time_MP: " << tmp << "\t" << "Niter: " << NrIter << "\t" << "Optimal Heuristic: " << OpV << "\t";
		#endif

		switch( Status )
		{
		case( NDOSolver::kOK ) :
			Out.precision( 16 );
			Out << "Status: OK, Value: " << OV << endl;
			break;
		case( NDOSolver::kStopped ) :
			Out.precision( 16 );
			Out << "Status: Stopped, Value: " << OV << endl;
			break;
		case( NDOSolver::kUnfsbl ) :
			Out << "Status: Unfeas." << endl;
			break;
		case( NDOSolver::kUnbndd ) :
			Out << "Status: Unbound." << endl;
			break;
		case( NDOSolver::kStpIter ) :
			Out.precision( 16 );
			Out << "Status: kStpIter, Value: " << OV << endl;
			break;
		case( NDOSolver::kStpTime ) :
			Out.precision( 16 );
			Out << "Status: kStpTime, Value: " << OV << endl;
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
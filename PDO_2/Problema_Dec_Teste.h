/*--------------------------------------------------------------------------*/
/*------------------------ File Problema_Dec_Teste.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Interface para resolver o problema por decomposição
 * Para varios parametros (t^* e tinit) e OSIsolvers diferentes
 * Definição dos valores dos parametros na funçao ResolverRL()
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#pragma once

#define WHICH_NDOSOLVER 0
/**< Select what NDOSolver is used to solve the problem.
   2 solvers is available:
    - 0 ==> Bundle
    - 1 ==> Volume + Bundle
    - 2 ==> Volume */

#if WHICH_NDOSOLVER <= 1
	#define WHICH_MPSOLVER 2
	/**< Select what MPSolver solver is used in the Bundle class.
	- 0 ==> QPPnltyMO
	- 1 ==> OSIMPSolver with quadratic stabilizing term (CPLEX)
	- 2 ==> OSIMPSolver with linear (box) stabilizing term */

 #if WHICH_MPSOLVER == 2
	#define WHICH_RMPSOLVER 0
	/**< Select what LP solver is used in MPSolver class.
	3 of the OsiSolverInterface compatible solvers are available:
	- 0 ==> CPLEX
	- 1 ==> CLP (COIN-OR)
	- 2 ==> Gurobi
	- 3 ==> Varios solvers */
		// Loop para vários OSIsolvers:
		// 0 ==> OSIMPSolver with linear (box) stabilizing term (CPLEX)
		// 1 ==> OSIMPSolver with linear (box) stabilizing term (CLP)
		// 2 ==> OSIMPSolver with linear (box) stabilizing term (Gurobi)

 #endif
#endif

#define WHICH_DEC 3		// Decomposition strategy  (colocar as demais e arrumar calculo x_med para cada uma na classe de resultados)
	// 0 ==> Teste (antiga) -> Space Decomposition: Spcdec2Fi
	// 1 ==> Space Decomposition: SpcDec_1
	// 2 ==> Space Decomposition: SpcDec_2
	// 3 ==> Space Decomposition: SpcDec_3
	// 4 ==> Scenario Decomposition: ScnDec_1
	// 5 ==> Scenario Decomposition: ScnDec_2
	// 6 ==> Hybrid Decomposition: HbrDec_1
	// 7 ==> Hybrid Decomposition: HbrDec_2

#define WARM_START 1	// fazer um modelo para a decomposição por cenários tb e otimizar código para somente entregar a f.o. ótima e os multiplicadores
	// 0 ==> No
	// 1 ==> Yes	-> usa a solução dual do modelo estendido resolvido no Gurobi
	// 2 ==> Yes	-> usa a solução dual do modelo estendido resolvido no Gurobi desconsiderando a rede

#define HEURISTIC 0		// Heuristic
	// 0 ==> No
	// 1 ==> Yes

#define MODEL_CHRCTR 4	// Kind of the Dual problem solved
	// 0 ==> Agreggated
	// 1 ==> Disaggregated
	// 2 ==> Disaggregated with Easy Components_1		(para a decomposiçao Space, tem 3 variantes possíveis, considerar isso tudo em WHICH_DEC???, com as outras possíveis decomp. espaciais...)
	// 3 ==> Disaggregated with Easy Components_2
	// 4 ==> Disaggregated with Easy Components_3

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( WHICH_DEC == 1 )
	#include "Spcdec1Fi.h"
	typedef Spcdec1Results Resultados;
#elif( WHICH_DEC == 2 )
	#include "Spcdec2Fi.h"
	typedef Spcdec2Results Resultados;
#elif( WHICH_DEC == 3 )
	#include "Spcdec3Fi.h"
	typedef Spcdec3Results Resultados;
	typedef Hrstc Heuristicas;
#elif( WHICH_DEC == 4 )
	#include "ScnDec1Fi.h"
	typedef Scndec1Results Resultados;
#elif( WHICH_DEC == 5 )
	#include "ScnDec2Fi.h"
	typedef Scndec2Results Resultados;
#elif( WHICH_DEC == 6 )
	#include "HbrDec1Fi.h"
	typedef Hbrdec1Results Resultados;
#elif( WHICH_DEC == 7 )
	#include "HbrDec2Fi.h"
	typedef Hbrdec2Results Resultados;
#else( WHICH_DEC == 0 )
	#include "EspFiOracle.h"
	typedef ResultadosConj Resultados;		// define um tipo genérico (Resultados) para o ponteiro, que vai depender de qual decomposição é usada
#endif

#if( WARM_START )
	#if( WHICH_DEC == 1 )
		#include "Spcdec1ExtDE.h"
		typedef Spcdec1ExtDE ExtDE;
	#elif( WHICH_DEC == 2 )
		#include "Spcdec2ExtDE.h"
		typedef Spcdec2ExtDE ExtDE;
	#elif( WHICH_DEC == 3 )
		#include "Spcdec3ExtDE.h"
		typedef Spcdec3ExtDE ExtDE;
	#elif( WHICH_DEC == 4 )
		#include "Scndec1ExtDE.h"
		typedef Scndec1ExtDE ExtDE;
	#elif( WHICH_DEC == 5 )
		#include "Scndec2ExtDE.h"
	#elif( WHICH_DEC == 6 )
		#include "HbrDec1ExtDE.h"
	#elif( WHICH_DEC == 7 )
		#include "HbrDec2ExtDE.h"
	#else( WHICH_DEC == 0 )
		#include "Problema_ED_ext.h"
		typedef Problema_ED_ext ExtDE;
	#endif
#else
	#include "OPTtypes.h"
	using namespace OPTtypes_di_unipi_it;
#endif

#include <ctime>

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Problema_Dec -------------------------*/
/*--------------------------------------------------------------------------*/

class Problema_Dec_Teste
{
	CSistema * sistema_a;
	int tipo_param;			// parametros do bundle = 1, parametros da heuristica = 2
	Resultados * resultadosGurobi;
	OPTtimers *startt;
	double refere_value;

	int ResolverRL1(string arv, string real);
	int ChamarBundle1(const char *const logF, const char *const logMP, const HpNum tstar_, const HpNum tinit_, const int solver, int iter);
	int ResolverRL2(string arv, string real);
	int ChamarBundle2(const char *const logF, const char *const logMP, int ws_a, int wlas_a, int bvmw, int bvfw, double alfa_a, double beta_a, double gama_a );

public:

	Problema_Dec_Teste(CSistema * const sistema_end, int tipo);
	~Problema_Dec_Teste(void);

	int ResolverRL(string arv, string real);

};

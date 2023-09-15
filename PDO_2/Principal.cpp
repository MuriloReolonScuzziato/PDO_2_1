/*--------------------------------------------------------------------------*/
/*----------------------------- File Principal.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * main() to solve de UUC problema.
 *
 * \version 1
 *
 * \date 29 - 09 - 2023
 *
 * \author Murilo Reolon Scuzziato \n
 *         Universidade Federal de Santa Catarina \n
 *
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

//#define MULT_ENTRADAS 0
//// 0 ==> Rodar problema para um unico dado de entrada padrao
//// N ==> Rodar problema para varios (N) dados de entrada diferentes (diferentes sorteios das realizaçoes)
//#if ( MULT_ENTRADAS )
//	#define TAM_ARVi 1
//	#define TAM_ARVf 1
//	// Tamanho da arvore de cenarios
//	// 1: 1 cenário para demanda e 1 para as afluencias
//	// 2: 2 cenário para demanda e 2 para as afluencias = 4 cenários
//	// ...
//#endif

#define WHICH_METHOD 0
// 0 ==> Deterministic Equivalent
// 1 ==> Decomposition	-> tipo de decomposicao definida na classe Problema_Dec
// 2 ==> Only the Heuristic

//#if( WHICH_METHOD == 0 ) -> implementar, resolução através da interface OSI (não permite funçao objetivo quadrática) !!!
//	#define WHICH_SOLVER 0
//	// 0 ==> Gurobi
//	// 1 ==> Cplex
//#endif

#if (WHICH_METHOD)
	#define MULT_PAR 0
	// 0 ==> Um caso unico com os parametros dos arquivos
	// 1 ==> Uma bateria de teste com parametros e solvers (do bundle) diferentes para um mesmo problema
	// 2 ==> Uma bateria de teste com parametros da heuristica diferentes para um mesmo problema
#endif

#define WHICH_RUN 0
// 0 ==> Debug local (executed with Visual, the Debug or Release Configuration)
// 1 ==> Release local (.exe)
// 2 ==> Release remote	(.exe)

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iomanip>   
using std::setw;
using std::setfill;
using std::right;

#include <time.h>

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( WHICH_METHOD == 0 )
	#include "Problema_ED.h"
	#include <windows.h>				// Ver biblioteca para contar tempo das classes do Frangioni, usar funçoes da bib. ctime
#else
	#if ( !MULT_PAR)
		#include "Problema_Dec.h"
	#else
		#include "Problema_Dec_Teste.h"
	#endif
#endif

//#include <vld.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
	using namespace NDO_di_unipi_it;
#endif
	
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// Arquivos de entrada (dados do sistema modelado)

#if (WHICH_RUN < 2)
const string pasta_dados = "C:/Dados_PDO/Entrada/Sistema_82bar/";
const string pasta_resultados = "C:/Dados_PDO/Resultados/Sistema_82bar/";
	
//const string pasta_dados = "C:/Dados_PDO/Entrada/System_46buses/";
//const string pasta_resultados = "C:/Dados_PDO/Resultados/System_46buses/";

//const string pasta_dados = "C:/Dados_PDO/Entrada/Sistema_4583bar/";
//const string pasta_resultados = "C:/Dados_PDO/Resultados/Sistema_4583bar/";

//const string pasta_dados = "D:/Dados_PDO/Entrada/Sistema_pequeno/";
//const string pasta_resultados = "D:/Dados_PDO/Resultados/Sistema_pequeno/";

//const string pasta_dados = "D:/Dados_PDO/Entrada/Sistema_teste/";
//const string pasta_resultados = "D:/Dados_PDO/Resultados/Sistema_teste/";
#endif

#if (WHICH_RUN == 2)
//const string pasta_dados = "D:/Murilo/PDO/Dados/Sistema_teste/";
//const string pasta_resultados = "D:/Murilo/PDO/Resultados/Sistema_teste/";

//Sistema_82bar
const string pasta_dados = "D:/Murilo/PDO2/Dados/";
const string pasta_resultados = "D:/Murilo/PDO2/Resultados/";

//const string pasta_dados = "D:/Murilo/PDO/Dados/Sistema_46bar/";
//const string pasta_resultados = "D:/Murilo/PDO/Resultados/Sistema_46bar/";

//const string pasta_dados = "D:/Murilo/PDO/Dados/Sistema_4583bar/";
//const string pasta_resultados = "D:/Murilo/PDO/Resultados/Sistema_4583bar/";
#endif

//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>


/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void str2val( const char* const str , T &sthg )
{
	istringstream( str ) >> sthg;
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/
#if (WHICH_RUN > 0)
int main(int argc, char *argv[])
#else
int main()
#endif
{
	//_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

#if (WHICH_RUN == 0)

	int argc = 14;
	char *argv[14];
	char *nome = "pdo_02";

	char *entr1 = "0";		// flag_modelo_rede					// modelo de rede de transmissão= 0: barra unica; 1: rede completa; 2: modelo compacto; 3: modelo compacto com lazy constraints
	char *entr2 = "0";		// flag_vfol (antigo v_meta)		// em alguns casos com =0 o método de resolução do ed pode influencia no resultado: casosgrb_intparam_method in problema_ed_ext
	char *entr3 = "0";		// flag_phmax (antigo valor_agua)	// se =0 para algumas decomposições entao relaxa-se a restrição de reserva girante
	char *entr4 = "10";		// number of initial linear cuts for the thermal cost function
	char *entr5 = "1";		// binary variables: 0 continuous; 1 binary;.
	char *entr6 = "1";		// 0 use least-square approximantion, 1 use perspective-cut formulation, 2 perspective-cut formulation dynamically (entr4 is the initial number of cuts (2 <= entr4 <= entr6), and entr6 is the maximum number of cuts that can be added)
	char *entr7 = "1";		// the binary variables of the thermal units: 0 old model; 1 modern model (tight and compact milp formulation for the thermal uc problem). 
	char *entr8 = "";		// caso (pasta) dos dados e resultados: "A/", "B/",...
	char *entr9 = "0";		// kind of preconditioner: 0 none of them; 1 probability of the node; 2 square root of the probability of the node; 3 maximum value of the relaxed constraint; 4 combinaion between 1 and 3; 5 combination between 2 and 3
	char *entr10 = "1";		// identification of the parameters file
	char *entr11 = "1";		// rodar problema para varios (n) dados de entrada diferentes (diferentes sorteios das realizaçoes)
	char *entr12 = "1";		// tamanho inicial da arvore de cenarios: =1, 1 cenário para demanda e 1 para as afluencias
	char *entr13 = "1";		// tamanho final da arvore de cenarios: =2, 2 cenário para demanda e 2 para as afluencias = 4 cenários

	// entr9 é só para decomposição
	// entr6 dynamically implementada somente para o ed
	argv[0] = nome;
	argv[1] = entr1;
	argv[2] = entr2;
	argv[3] = entr3;
	argv[4] = entr4;
	argv[5] = entr5;
	argv[6] = entr6;
	argv[7] = entr7;
	argv[8] = entr8;
	argv[9] = entr9;
	argv[10] = entr10;
	argv[11] = entr11;
	argv[12] = entr12;
	argv[13] = entr13;

#endif		
	// read the command-line parameters- - - - - - - - - - - - - - - - - - - - -

	if(argc != 14)
	{
		cout << "Numero de parametros incorreto!!" << endl;
		return ( 0 );
	}
	// flag_barra_unica, flag_volume_meta, flag_valor_agua, flag_aprox_custoT, flag_var_bin, flag_thermal_cost_approximation, flag_thermal_binary_model
	int flag1;
	str2val( argv[ 1 ] , flag1 );
	bool flag2;
	str2val( argv[ 2 ] , flag2 );
	bool flag3;
	str2val( argv[ 3 ] , flag3 );
	int flag4;
	str2val( argv[ 4 ] , flag4 );	// Numero de partes para linearizar a função de custo de produção das termicas (0:quadratica, 1:uma reta, 2:duas retas,...)
	bool flag5;
	str2val( argv[ 5 ] , flag5 );
	int flag6;
	str2val( argv[ 6 ] , flag6 );
	int flag7;
	str2val( argv[ 7 ] , flag7 );
	int flag8;
	str2val( argv[ 9 ] , flag8 );
	int mult_ent;
	str2val( argv[ 11 ] , mult_ent );
	int tam_ini;
	str2val( argv[ 12 ] , tam_ini );
	int tam_fin;
	str2val( argv[ 13 ] , tam_fin );

	// Ajustar entradas - - - - - - - - - - - - - - - - -
	vetorstring demanda, afluencias;
	ostringstream convert, convert2;
	if (!mult_ent)
	{
		demanda.push_back("Demanda.txt");
		afluencias.push_back("Afluencias.txt");
		convert << 0;
		convert2 << 0;
	}
	else
	{
		for (int j = tam_ini; j <= tam_fin; j++)
		{
			convert2 << pow(double(j), 2);
			for (int i = 1; i <= mult_ent; i++)
			{
				convert << i;
				demanda.push_back("Demanda" + convert2.str() + "_" + convert.str() + ".txt");
				afluencias.push_back("Afluencias" + convert2.str() + "_" + convert.str() + ".txt");
				//demanda.push_back("Demanda" + convert2.str() + ".txt");
				//afluencias.push_back("Afluencias" + convert2.str() + ".txt");
				convert.str("");
			}
			convert2.str("");
		}
	}

	// Loop para resolver um problema com multiplas entradas
	CSistema * sistema;
	for (size_t n_ent = 0; n_ent < demanda.size(); n_ent++)
	{
		// Criar objeto Sistema e ajustar flags - - - - - - - - - - - - - - - - -
		sistema = new CSistema( pasta_dados + argv[8], "Hidreletricas.txt", "UnidadesH.txt", "Termeletricas.txt", "Parametros.txt", demanda[n_ent], "Barras.txt", "Linhas.txt", afluencias[n_ent], "Cond_Inic_H.txt", "Cond_Inic_T.txt");
		sistema->SetFlags(flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8);
		//cout << "modelo carregado!" << endl;
		// Resolver o problema de UUC- - - - - - - - - - - - - - - - - - - - - - - -
		#if( WHICH_METHOD == 0 )	// ED - arquivos: Problema_ED.cpp(h)
		{
			CProblema_ED * modeloED;
			double tempoModelo, tempoResolucao;
			LARGE_INTEGER clockPerSec;
			LARGE_INTEGER inicio, fim;
			QueryPerformanceFrequency(&clockPerSec);
			QueryPerformanceCounter(&inicio);
			modeloED = new CProblema_ED(sistema);
			QueryPerformanceCounter(&fim);
			tempoModelo = double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart);
			//modeloED->ConferirCoerenciaDados();
			cout << "ED: ";
			QueryPerformanceCounter(&inicio);
			modeloED->ResolverProblema();
			QueryPerformanceCounter(&fim);
			tempoResolucao = double(fim.QuadPart - inicio.QuadPart)/double(clockPerSec.QuadPart);
			if (mult_ent != 0)
			{
				convert << n_ent%mult_ent + 1;
				convert2 << pow(double(n_ent/mult_ent + tam_ini), 2);
			}
			modeloED->EscreverResultados("out_ED_" + convert2.str() + "_" + convert.str() + ".txt", tempoModelo, tempoResolucao, pasta_resultados + argv[8]);
			convert.str("");
			convert2.str("");
			delete modeloED;
		}
		#elif ( WHICH_METHOD == 1 )		// Decomposiçao - arquivos: Problema_Dec.cpp(h) + Oracles + Subproblems
		{
			// Problema_Dec que decide qual decomposição vai utilizar e chama o bundle
			// Dec. Espacial - arquivos: Spcdec2Fi.cpp(h); SubProblemaHA.cpp(h); SubProblemaHE.cpp(h); SubProblemaT.cpp(h); SubProblemaD.cpp(h)
			// Dec. por Cenarios - arquivos: CenFiOracle.cpp(h); SubProblemaCen.cpp(h)

			#if ( !MULT_PAR)
				if (mult_ent != 0)
				{
					convert << n_ent%mult_ent + 1;
					convert2 << pow(double(n_ent/mult_ent + tam_ini), 2);
				}
				//Problema_Dec modeloDec(sistema, &resultadosDecEsp);
				Problema_Dec * modeloDec;
				modeloDec = new Problema_Dec(sistema, argv[10]);
//				cout << "modelo decomp. criado!" << endl;
			#else
				// nessa classe resolvem-se vários problemas com ajustes diferentes
				if (mult_ent != 0)
				{
					convert << n_ent%mult_ent + 1;
					convert2 << pow(double(n_ent/mult_ent + tam_ini), 2);
				}
				Problema_Dec_Teste * modeloDec;
				modeloDec = new Problema_Dec_Teste(sistema, MULT_PAR);
			#endif
			modeloDec->ResolverRL(convert2.str(), convert.str());
//			cout << "RL resolvida!" << endl;
			//modeloDec->Teste();
			convert.str("");
			convert2.str("");
			delete modeloDec;
			
			// Colocar todos flags de ajust no arquivo de saida!!!
			
			// Conferir soluções
			//resultadosDecEsp.AlocarX_med();
			//resultadosDecEsp.ExportarXmed("x.txt");

		}
		#else
		{
			if (mult_ent != 0)
			{
				convert << n_ent%mult_ent + 1;
				convert2 << pow(double(n_ent/mult_ent + tam_ini), 2);
			}
			//Problema_Dec modeloDec(sistema, &resultadosDecEsp);
			Problema_Dec * modeloDec;
			modeloDec = new Problema_Dec(sistema, argv[10]);

			modeloDec->ResolverHrst();
			convert.str("");
			convert2.str("");
			delete modeloDec;
		}
		#endif

		delete sistema;
	}

	return ( 0 );
}

// contar tempo
// #include <time.h>
// clock_t tempo;
// tempo = clock();
// ...
// double(clock() - tempo)/double(CLOCKS_PER_SEC)

// pause
//std::cin.ignore();

// parei aqui
// deixar como opcao resolver o PL da heuristica tb!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! precisa dele: qual seria a solucao sem:

// as informações estao sendo armazenadas corretamente no ptr_x_til e no ptr_x_subp dessa classe
// mas eles tem a precisao do double, entao 1512 pode ficar 1511.99999999999, e isso n pode ser usado como entrada de outro subproblema ou de RHS das restrições proximais???
// alem disso tem inviabilida 9 no servidor 127, tratar essa questão numerica na heuristica!!! conversar com frangioni

// servidor da resultados diferentes do q esse pc?!?

// nao precisa deixar redonda a solucao, tem q ver o q ta causando o x_til = 0. provavelmente ALFA = 0... verificar onde ele eh alterado, n deixar isso, nem q ocorra erro... alias verificar esse erro!!!

// conferir os vetores x_til e x_hat da classe heuristica!!

// parei aqui
// criar string com endereço de saida na classe de resultados!!! pra escrever resultados a partir do FiOracle!

// deixar tempo limite dos subproblemas de acordo com o tempo total limite, uma porcentagem dele... (levando em consideração um numero de possiveis iterações...)
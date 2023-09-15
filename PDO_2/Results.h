#pragma once

// Definition of the abstract base class FiOracle, which sets the interface
// for the individual (derived) results classes

#include "Sistema.h"

#include "MatrizEsparsa.h"

#include "gurobi_c++.h"
using namespace std;

typedef vector<string> vetorstring;
typedef vector<GRBModel *> vetorGRBModel;
typedef vector<GRBVar *> vetorGRBVar;
typedef vector<GRBConstr *> vetorGRBCons;

class Results
{
public:
	// Atributos gerais para todas as decomposições
	CSistema * sistema_a;
	int flag1, flag2, flag3, flag4, flag7;
	vetorfloat x, xa, obj_subp, x_med, lambda;
	int n, na, nd, CH, R, N, I, status;
	double obj, CFO, CIO, DEF, CTL, CTQ, VFOL, mipgap;

	vetorstring stringVariaveis;
	ofstream * inFile;

	vetorfloat2 X_til;				// armazenar x_til de todas as iterações para conferencia
	vetorfloat2 X_hat;				// armazenar x_til de todas as iterações para conferencia

	Results(void)
	{
		sistema_a = NULL;
		flag1 = 0;flag2 = 0;flag3 = 0;flag4 = 0;flag7 = 0;
		x.resize(0);xa.resize(0);obj_subp.resize(0);x_med.resize(0);lambda.resize(0);
		n = 0; na = 0; nd = 0; CH = 0; R = 0; N = 0; I = 0; status = 0;
		obj = 0; CFO = 0; CIO = 0; DEF = 0; CTL = 0; CTQ = 0; VFOL = 0; mipgap = 0;
		stringVariaveis.resize(0);
		inFile = NULL;
	}
	virtual ~Results(void) {}

	// funçoes q sao usadas na heuristica, devem ser declaradas virtuais aqui, pois a heuristica "ve" um objeto da classe Results
	virtual vetorfloat GetX_med(bool hat_or_til) {cout << "Virtual" << endl; return vetorfloat();}
	double GetLambda(int posicao) {return lambda[posicao];}

	virtual void AlocarXmed(CMatrizEsparsa * x_spr) {cout << "Virtual" << endl; }

	// funçoes genéricas para todas as decomposições
	inline void EscreverArquivo(string nome_arquivo, double tempo_modelo, double tempo_resol, string saida);
	inline void GravarSolHeuristica(double obj_a, CMatrizEsparsa * x_spr, double LB);
	inline void GravarSolExtDE(GRBVar * vars, GRBConstr * constrs, int nvar, int nvar_a, int nconstr, bool ScnDec = false);
	inline void GravarSolExtDE(int nvar, int nvar_a, int nconstr, bool ScnDec = false);
	inline CMatrizEsparsa MatrizCalcFluxo(int n_a);
	inline CMatrizEsparsa MatrizCalcFluxo_cen(int n_a);
	inline void ArmazenarXtilhat();
	inline void ImprimirArmazenarXtilhat();
};

/* @} end( Hrst ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void Results::EscreverArquivo(string nome_arquivo, double tempo_modelo, double tempo_resol, string saida)
{
	int flag2 = int (sistema_a->GetFlagVfol());
	size_t n = x.size();
	size_t nt = stringVariaveis.size() - flag2*sistema_a->hidreletricasVtr.size();
	int num_cen = sistema_a->GetNCenarios();
	//num_cen = 1;		// para escrever resultado de um cenario
	int T1 = sistema_a->GetTt1();
	int T2 = sistema_a->GetTt2();
	//CreateDirectoryA(varString.c_str(), NULL);		// criar pasta
	inFile = new ofstream( saida + nome_arquivo, ios::out );
	if ( inFile->is_open() )
	{
		if (status != 2)
		{
			*inFile << "Solução não ótima! " << "Status final = " << status << char(9) << "Gap = " << mipgap << endl;
			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagModeloRede() << char(9) << sistema_a->GetFlagVfol() << char(9) << sistema_a->GetFlagPhmax() << char(9) << sistema_a->GetFlagInitAproxCT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
		}
		else
		{
			*inFile << "Solução ótima! " << "Status final = " << status << char(9) << "Gap = " << mipgap << endl;
			*inFile << "Funcao objetivo :" << char(9) << std::scientific << setprecision(10) << obj << char(9) << "DEF (MW): " << DEF << char(9) << "VFol (hm3): " << VFOL << char(9) << "CIO : " << CIO << char(9) << "CFO : " << CFO << endl;
			*inFile << "CTL :" << CTL << char(9) << "CTQ :" << CTQ << endl;
			*inFile << "Tempo de resolucao (s):" << char(9) << std::scientific << setprecision(4) << tempo_resol << char(9) << "Criar modelo (s):" << tempo_modelo << endl;
			*inFile	<< "flags =" << char(9) << sistema_a->GetFlagModeloRede() << char(9) << sistema_a->GetFlagVfol() << char(9) << sistema_a->GetFlagPhmax() << char(9) << sistema_a->GetFlagInitAproxCT() << char(9) << sistema_a->GetFlagVarBin() << char(9) << sistema_a->GetFlagMaxAproxCT() << char(9) << sistema_a->GetFlagTbinaryModel() << endl;
		}
		*inFile << endl;
		// Estágio 1
		*inFile << setw(10) << left << setprecision(6) << "Estágio 1" << endl;
		*inFile << setw(10) << left << setprecision(6) << "Periodo";
		for (int t = 0; t < T1; t++)
			*inFile << setw(14) << left << char(9) << t + 1;
		for (size_t i = 0; i < nt; i++)
		{
			*inFile << endl;
			*inFile << setw(10) << left << stringVariaveis[i];
			for (int t = 0; t < T1; t++)
				*inFile << char(9) << setw(15) << right << double(x[i + t*nt] > 0 ? double (x[i + t*nt] <= 1e-10 ? 0 : x[i + t*nt]) : x[i + t*nt]);
		}
		// Estágio 2
		for (int n_c = 0; n_c < num_cen; n_c++)
		{
			*inFile << endl << endl;
			*inFile << setw(10) << left << setprecision(6) << "Estágio 2 - Cenário "<< n_c + 1 << endl;
			*inFile << setw(10) << left << setprecision(6) << "Periodo";
			for (int t = T1; t < T2; t++)
				*inFile << setw(14) << left << char(9) << t + 1;
			for (size_t i = 0; i < nt; i++)
			{
				*inFile << endl;
				*inFile << setw(10) << left << stringVariaveis[i];
				for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1); t++)
				{
					*inFile << char(9) << setw(15) << right << double(x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] > 0 ? double (x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] <= 1e-10 ? 0 : x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]) : x[i + t*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]);
				}
			}
			if (sistema_a->GetFlagVfol() == 1)
				for (size_t i = 0; i < sistema_a->hidreletricasVtr.size(); i++)
				{
					*inFile << endl;
					*inFile << setw(10) << left << stringVariaveis[i + nt];
					for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1) - 1; t++)
						*inFile << char(9) << setw(15) << right << " - ";
					*inFile << char(9) << setw(15) << right << double (x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] > 0 ? double (x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()] <= 1e-10 ? 0 : x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]) : x[i + (T2 + n_c* (T2-T1))*nt + flag2*n_c*sistema_a->hidreletricasVtr.size()]);
				}
		}
		// Fluxo nas linhas
		vetorfloat fluxo;
		fluxo.resize(sistema_a->linhasVtr.size()*(T1 + num_cen * (T2 - T1)));
		if (sistema_a->GetFlagModeloRede() == 1)
		{
			vetorfloat teta;
			CMatrizEsparsa M(0);
			M = MatrizCalcFluxo(int (x.size()) );		//Multiplicar matriz MatrizLimFluxo pelo vetor x
			fluxo = M.MultiplicarPorVetorDenso(&x);
		}
		else	// sistema_a->GetFlagModeloRede() != 1
		{
			int delta = 0;
			int cen = 0;
			double D = 0;
			int JJ = 0;
			for (size_t i = 0; i < R; i++)
				JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
			MatrixXd Agh = MatrixXd::Zero(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
					for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
						if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
							Agh(b, r) = 1;
			MatrixXd Agt = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
						if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
							Agt(b, i) = 1;
			MatrixXd Agtu = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
			if (sistema_a->GetFlagTbinaryModel() == 1)
			{
				for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
					for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
						for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
							if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
								Agtu(b, i) = sistema_a->termeletricasVtr[i].GetPmin();
			}
			MatrixXd Adef = MatrixXd::Identity(sistema_a->barrasVtr.size(),sistema_a->barrasVtr.size());
			Agh = - sistema_a->Beta * Agh;
			Agt = - sistema_a->Beta * Agt;
			Agtu = - sistema_a->Beta * Agtu;
			Adef = - sistema_a->Beta * Adef;
			VectorXd Ad(sistema_a->barrasVtr.size());
			VectorXd gg(sistema_a->linhasVtr.size());
			VectorXd g, Rhs;
			VectorXd flow(sistema_a->linhasVtr.size());
			for (int t = 0; t < N; t++)
			{
				Ad.resize(sistema_a->barrasVtr.size());
				if (sistema_a->GetTt1() <= t)
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				// Determinar o vetor de demandas
				for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				{
					D = 0;
					for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
						D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
					Ad(b) = D;
				}
				Rhs = - sistema_a->Beta * ( Ad );
				// Determinar o vetor de gerações
				gg = VectorXd::Zero(sistema_a->linhasVtr.size());
				g.resize(R);
				for (size_t r = 0; r < R; r++)
					g(r) = x[r + delta + (3+flag4+flag7)*I];
				gg = gg + Agh * g;
				g.resize(I);
				for (size_t i = 0; i < I; i++)
					g(i) = x[i + delta];
				gg = gg + Agt * g;
				g.resize(I);
				if (sistema_a->GetFlagTbinaryModel() == 1)
					for (size_t i = 0; i < I; i++)
						g(i) = x[i + delta + I];
				gg = gg + Agtu * g;
				g.resize(sistema_a->barrasVtr.size());
				if ( sistema_a->GetFlagModeloRede() == 0 )		// modelo de barra unica, deficit alocado igualmente entre todas barras, o que não é verdade!
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						g(b) = x[delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ] / sistema_a->barrasVtr.size();
				}
				else
				{
					for (size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
						g(b) = x[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
				}
				gg = gg + Adef * g;
				// Conferir limites das linhas extrapolados para cada periodo
				flow = - gg + Rhs;		// Beta ( d - g): fluxo nas linhas para o periodo t
				for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
					fluxo[l + t*sistema_a->linhasVtr.size()] = flow(l);
				if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
					delta += nt + flag2*R;
				else
					delta += nt;
			}
		}
		// Escrever fluxo
		*inFile << endl << endl;
		if ( sistema_a->GetFlagModeloRede() == 0 )	
		{
			*inFile << "Atenção: Modelo de barra única, caso exista déficit ele será alocado igualmente entre as barras! Portanto o resultado pode ser diferente do modelo de rede.";
			*inFile << endl << endl;
		}
		*inFile << setw(10) << left << setprecision(5) << "Fluxo nas linhas (MW)" << endl;
		*inFile << setw(10) << left << setprecision(5) << "Estágio 1" << endl;
		for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		{
			*inFile << setprecision(5) << "L" ;
			*inFile << setprecision(5) << l + 1;
			*inFile << setprecision(1) << std::fixed << " (" << sistema_a->linhasVtr[l].GetCapacidade() << ")" << setw(8);
			*inFile << setprecision(5) << std::scientific;
			for (int t = 0; t < T1; t++)
			{
				if ((fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade()) ||
					(fluxo[l + t*sistema_a->linhasVtr.size()] <= - 0.9999*sistema_a->linhasVtr[l].GetCapacidade()))
					*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()] << "*";
				else
					*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()];

			}
			*inFile << endl;
		}
		for (int n_c = 0; n_c < num_cen; n_c++)
		{
			*inFile << setw(10) << left << setprecision(6) << "Estágio 2 - Cenário "<< n_c + 1 << endl;
			for (size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
			{
				*inFile << setprecision(5) << "L" ;
				*inFile << setprecision(5) << l + 1;
				*inFile << setprecision(1) << std::fixed << " (" << sistema_a->linhasVtr[l].GetCapacidade() << ")" << setw(8);
				*inFile << setprecision(5) << std::scientific;
				for (int t = T1 + n_c * (T2-T1); t < T2 + n_c * (T2-T1); t++)
				{
					if ((fluxo[l + t*sistema_a->linhasVtr.size()] >= 0.9999*sistema_a->linhasVtr[l].GetCapacidade()) ||
						(fluxo[l + t*sistema_a->linhasVtr.size()] <= - 0.9999*sistema_a->linhasVtr[l].GetCapacidade()))
						*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()] << "*";
					else
						*inFile << char(9) << setw(15) << left << fluxo[l + t*sistema_a->linhasVtr.size()];
				}
				*inFile << endl;
			}
		}
		inFile->close();
	}
	else
		cout << "Unable to open file";
}
inline void Results::GravarSolHeuristica(double obj_a, CMatrizEsparsa * x_spr, double LB)
{
	obj = obj_a;
	mipgap = (obj_a - LB) / (LB + 1);

	for (size_t i = 0; i < n; i++)
		x[i] = x_spr->GetElemento(i, 0);

	CIO = 0;		// Custo imediato de operação (custo do primeiro estágio), considerando o custo de operação e partida das termelétricas
	CTL = 0;		// Custo de operação das termelétricas linear
	CTQ = 0;		// Custo de operação das termelétricas quadrático
	CFO = 0;		// Custo futuro de operação ( valor incremental da água, da função de custo futuro)
	DEF = 0;		// Somatorio dos déficits
	VFOL = 0;		// Somatorio dos vfolga

	int nt;
	//n = x.size();
	int n_a = n - flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size();
	double deltaT;
	int delta = 0;
	int flag1a = 0;		// referente à var. teta
	if (flag1 == 1)
		flag1a = 1;
	int flag1d = 1;		// referente à var. def
	if (flag1 == 0)
		flag1d = 0;
	int cen;
	for (int t = 0; t < N; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			nt = (n_a / N) + flag2*sistema_a->hidreletricasVtr.size();
		else
			nt = (n_a / N);
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			deltaT = sistema_a->GetDeltaT1();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				if (sistema_a->GetFlagTbinaryModel() == 0)
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT;		// Custo de operaçao quadrático
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						CTL += x[i + 3*sistema_a->termeletricasVtr.size() + delta];		// Custo de operaçao linear
						CIO += x[i + 2*sistema_a->termeletricasVtr.size() + delta] + x[i + 3*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
					}
					else
						CIO += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
				}
				else
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT;		// Custo de operaçao quadrático
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta];		// Custo de operaçao linear
						CIO += x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta];		// Custo de 1 estagio
					}
					else
						CIO += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida();		// Custo de 1 estagio
				}
			}
		}
		else
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			deltaT = sistema_a->GetDeltaT2();
			for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			{
				if (sistema_a->GetFlagTbinaryModel() == 0)
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta])*deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						CTL += x[i + 3*sistema_a->termeletricasVtr.size() + delta] * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
						CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] + x[i + 3*sistema_a->termeletricasVtr.size() + delta]) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
					}
					else
						CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) + x[i + 2*sistema_a->termeletricasVtr.size() + delta] ) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
				}
				else
				{
					CTQ += sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao quadrático
					if (sistema_a->GetFlagInitAproxCT() > 1)
					{
						CTL += x[i + 4*sistema_a->termeletricasVtr.size() + delta]* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de operaçao linear
						CFO += (x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida() + x[i + 4*sistema_a->termeletricasVtr.size() + delta])* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
					}
					else
						CFO += (sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta] + sistema_a->termeletricasVtr[i].GetPmin()*x[i + sistema_a->termeletricasVtr.size() + delta], x[i + sistema_a->termeletricasVtr.size() + delta]) * deltaT + x[i + 2*sistema_a->termeletricasVtr.size() + delta] * sistema_a->termeletricasVtr[i].GetCoefCustoPartida())* sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		// Custo de 2 estagio
				}
			}
		}
		delta = delta + nt;
	}

	delta = (n_a / N) - (1 - flag1d) - flag1d*sistema_a->barrasVtr.size(); //(n_a / N) = 34 52
	for (int t = 0; t < N; t++)
	{
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			nt = (n_a / N) + flag2*sistema_a->hidreletricasVtr.size();
		else
			nt = (n_a / N);
		if ((0 <= t) && (sistema_a->GetTt1() > t))
		{
			deltaT = sistema_a->GetDeltaT1();
			for (size_t i = 0; i < (1 - flag1d) + flag1d*sistema_a->barrasVtr.size(); i++)
				DEF += x[i + delta];		// Def
		}
		else
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			deltaT = sistema_a->GetDeltaT2();
			for (size_t i = 0; i < (1 - flag1d) + flag1d*sistema_a->barrasVtr.size(); i++)
				DEF += x[i + delta];					// Def
		}
		delta = delta + nt;
	}

	if (sistema_a->GetFlagVfol() == 1)
	{
		delta = (n_a / N) * sistema_a->GetTt2();
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				VFOL += x[delta++];					// vfol
			delta += (n_a / N) * (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
	}
	else
		VFOL = 0;
}
inline CMatrizEsparsa Results::MatrizCalcFluxo(int n_a)
{
	CMatrizEsparsa Alb(int (sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
	{
		Alb.InserirElemento(int (l), sistema_a->linhasVtr[l].de_barra->GetIdentBarra(), 1);
		Alb.InserirElemento(int (l), sistema_a->linhasVtr[l].para_barra->GetIdentBarra(), -1);
	}
	Alb.RemoverColuna(sistema_a->GetBarraRef() - 1);
	CMatrizEsparsa TT(int (sistema_a->linhasVtr.size()), int (sistema_a->linhasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		TT.InserirElemento(int(l), int(l), 100 / (sistema_a->linhasVtr[l].GetReatancia()));
	CMatrizEsparsa Alin(int (N * sistema_a->linhasVtr.size()), n_a);
	n_a = n_a - int(flag2*sistema_a->GetNCenarios()*sistema_a->hidreletricasVtr.size());
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < N; t++)
	{
		Alin.InserirMatriz(ll,c + int ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),ll + int (sistema_a->linhasVtr.size()) - 1,c + int((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
		ll = ll + int(sistema_a->linhasVtr.size());
		if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
			c += (n_a / N) + int(flag2*sistema_a->hidreletricasVtr.size());
		else
			c += (n_a / N);
	}
	return Alin;
}
inline CMatrizEsparsa Results::MatrizCalcFluxo_cen(int n_a)
{
	CMatrizEsparsa Alb(int(sistema_a->linhasVtr.size()), int(sistema_a->barrasVtr.size()));
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
	{
		Alb.InserirElemento(int (l), sistema_a->linhasVtr[l].de_barra->GetIdentBarra(), 1);
		Alb.InserirElemento(int (l), sistema_a->linhasVtr[l].para_barra->GetIdentBarra(), -1);
	}
	Alb.RemoverColuna(sistema_a->GetBarraRef() - 1);
	CMatrizEsparsa TT(sistema_a->linhasVtr.size(),sistema_a->linhasVtr.size());
	for ( size_t l = 0; l < sistema_a->linhasVtr.size(); l++)
		TT.InserirElemento(l, l, 100 / (sistema_a->linhasVtr[l].GetReatancia()));
	CMatrizEsparsa Alin(N * sistema_a->linhasVtr.size(),n_a);
	n_a = n_a - sistema_a->hidreletricasVtr.size();
	int ll = 0;
	int c = 0;
	TT.MultiplicarPorMatriz(&Alb);
	for ( int t = 0; t < N; t++)
	{
		Alin.InserirMatriz(ll,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size()),ll + sistema_a->linhasVtr.size() - 1,c + ((3+flag4+flag7)*sistema_a->termeletricasVtr.size() + sistema_a->barrasVtr.size()) - 2, &TT, 0, 0);
		ll = ll + sistema_a->linhasVtr.size();
		c += (n_a / N);
	}
	return Alin;
}
inline void Results::GravarSolExtDE(GRBVar * vars, GRBConstr * constrs, int nvar, int nvar_a, int nconstr, bool ScnDec)
{
	// gravar em x, xa e L os resultados do smart start
	for  (int i = 0; i < nvar; i++)
		x[i] = vars[i].get(GRB_DoubleAttr_X);
	for  (int i = 0; i < nvar_a; i++)
		xa[i] = vars[i + nvar].get(GRB_DoubleAttr_X);
	if (ScnDec)
	{
		int n_rest_per = (nconstr - sistema_a->GetNCenarios()*R)/(sistema_a->GetTt2()*sistema_a->GetNCenarios());
		lambda.resize(n_rest_per*N);
		double valor = 0;
		for (int t = 0; t < sistema_a->GetTt1(); t++)
		{
			for  (int i = 0; i < n_rest_per; i++)
			{
				valor = 0;
				for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					valor = constrs[i + n_rest_per*(t + cen*sistema_a->GetTt2())].get(GRB_DoubleAttr_Pi) * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
				lambda[i + n_rest_per*t] = valor;
			}
		}
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			for (int t = sistema_a->GetTt1(); t < sistema_a->GetTt2(); t++)
				for  (int i = 0; i < n_rest_per; i++)
					lambda[i + n_rest_per*(t + cen*(sistema_a->GetTt2()-sistema_a->GetTt1()))] = constrs[i + n_rest_per*(t + cen*(sistema_a->GetTt2()))].get(GRB_DoubleAttr_Pi);
	}
	else
	{
		lambda.resize(nconstr);
		for  (int i = 0; i < nconstr; i++)
			lambda[i] = constrs[i].get(GRB_DoubleAttr_Pi);
	}
	// lambda: restrições devem estar na mesma ordem no Ext DE e na heuristica
	///cout << "Results: conferir se restricoes do modelo Ext estao na mesma ordem que na Heuristica!" << endl;
}
inline void Results::GravarSolExtDE(int nvar, int nvar_a, int nconstr, bool ScnDec)
{
	// gravar em x, xa e L os resultados do smart start
	for  (int i = 0; i < nvar; i++)
		x[i] = 0;
	for  (int i = 0; i < nvar_a; i++)
		xa[i] = 0;
	if (ScnDec)
	{
		int n_rest_per = (nconstr - sistema_a->GetNCenarios()*R)/(sistema_a->GetTt2()*sistema_a->GetNCenarios());
		lambda.resize(n_rest_per*N);
		for  (int i = 0; i < n_rest_per*N; i++)
			lambda[i] = 0;
	}
	else
	{
		lambda.resize(nconstr);
		for  (int i = 0; i < nconstr; i++)
			lambda[i] = 0;
	}
}
inline void Results::ArmazenarXtilhat()
{
	vetorfloat temp = GetX_med(false);
	X_hat.push_back(temp);
	temp = GetX_med(true);
	X_til.push_back(temp);
}
inline void Results::ImprimirArmazenarXtilhat()
{

	// Escrever x_hat e x_til
	ofstream * inFile;
	inFile = new ofstream( "x_hat.txt", ios::out );
	if ( inFile->is_open() )                                                            
	{
		*inFile << std::scientific << setprecision(16);
		for (size_t i = 0; i < X_hat[0].size(); i++)
		{
			for (size_t ii = 0; ii < X_hat.size(); ii++)
				*inFile << X_hat[ii][i] << char(9);
			*inFile << endl;
		}
	}
	else
		cout << "Unable to open file";
	inFile->close();

	// Ler solução x
	inFile = new ofstream( "x_til.txt", ios::out );
	if ( inFile->is_open() )                                                            
	{
		*inFile << std::scientific << setprecision(16);
		for (size_t i = 0; i < X_til[0].size(); i++)
		{
			for (size_t ii = 0; ii < X_til.size(); ii++)
				*inFile << X_til[ii][i] << char(9);
			*inFile << endl;
		}
	}
	else
		cout << "Unable to open file";
	inFile->close();
}

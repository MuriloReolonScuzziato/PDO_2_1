#pragma once

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;
using std::setw;
using std::left;

#include <fstream>
using std::ifstream;
using std::ios;
using std::ofstream;

#include <vector>
using std::vector;

#include <eigen>

typedef vector<double> vetorfloat;
typedef vector<vetorfloat> vetorfloat2;
typedef vector<vetorfloat2> vetorfloat3;
typedef vector<int> vetorint;
typedef vector<vetorint> vetorint2;
typedef vector<vetorint2> vetorint3;
typedef vector<int>::const_iterator myiter;

// Algoritmo de ordenção (Gnome sort) - modificado
// ---------------------------------------------------------------------------------------------------
template<class T>
void gnome_sort( vector<T> &lista, vector<double> &lista2 )
{
	vector<T>::size_type i = 1;
 
	while( i < lista.size() )
	{
		if( i == 0 || lista.at( i-1 ) <= lista.at( i ) )
		{
			i++;
		}
		else
		{
			std::swap( lista[ i - 1 ], lista [ i ] );std::swap( lista2[ i - 1 ], lista2 [ i ] );
			--i;
		}
	}
}
// ---------------------------------------------------------------------------------------------------

// Algoritmo de ordenção por referencia
// ---------------------------------------------------------------------------------------------------
struct ordering 
{
	bool operator ()(std::pair<size_t, myiter> const& a, std::pair<size_t, myiter> const& b) 
	{
		return *(a.second) < *(b.second);
    }
};
template <typename T> void sort_from_ref( vector<T> & in, vector<std::pair<size_t, myiter> > const& reference)
{
    vector<T> ret(in.size());

    size_t const size = in.size();
    for (size_t i = 0; i < size; ++i)
        ret[i] = in[reference[i].first];
	in = ret;
}

struct ordering2
{
	bool operator ()(std::pair<size_t, std::pair<myiter, myiter>> const& a, std::pair<size_t, std::pair<myiter, myiter>> const& b) 
	{
		if (*(a.second.first) < *(b.second.first)) return true;
		if (*(a.second.first) == *(b.second.first)) return *(a.second.second) < *(b.second.second);
		return false;
    }
};
template <typename T> void sort_from_ref2( vector<T> & in, vector<std::pair<size_t, std::pair<myiter, myiter>> > const& reference)
{
    vector<T> ret(in.size());

    size_t const size = in.size();
    for (size_t i = 0; i < size; ++i)
        ret[i] = in[reference[i].first];
	in = ret;
}
// ---------------------------------------------------------------------------------------------------

class CMatrizEsparsa
{
	vetorfloat val;
	vetorint lprim, col, lprox;
	int nlin, ncol, nnz;

public:
	CMatrizEsparsa(void);
	CMatrizEsparsa(int linha, int coluna);
	CMatrizEsparsa(int linha_coluna);		// Cria uma matriz identidade quadrada
	CMatrizEsparsa(CMatrizEsparsa &matriz);		// Cria matriz copia
	~CMatrizEsparsa(void);

	void InserirElemento(int linhas, int colunas, double valor);
	void SubstituirElemento(int linhas, int colunas, double valor);
	void ImprimirLinha(int linha);
	void ImprimirMatriz();
	void ImprimirMatrizVisual();
	void ImprimirMatrizArquivo();
	void JuntarColuna(CMatrizEsparsa * matriz, bool teste = false);
	void RemoverColuna(int coluna);
	//void Transpor();
	void OrdenarLista();
	void RemoverElemento(int linha, int coluna);
	void MultiplicarPorMatriz(CMatrizEsparsa * m2);
	void MultiplicarPorMatriz(Eigen::MatrixXd &m2);
	void MultiplicarPorEscalar(double valor);
	void InserirMatriz(int linha_i, int coluna_i, int linha_f, int coluna_f, CMatrizEsparsa * m2, int linha2_i, int coluna2_i);
	void InserirMatriz(int linha_i, int coluna_i, int linha_f, int coluna_f, Eigen::MatrixXd &m2, int linha2_i, int coluna2_i);
	//CMatrizEsparsa MultiplicarMatrizes(CMatrizEsparsa * m1, CMatrizEsparsa * m2);
	void ZerarMatriz() {val.clear(); lprim.clear(); col.clear(); lprox.clear(); nlin = 0; ncol = 0; nnz = 0;}
	void ZerarMatriz(int n_linhas, int n_colunas);
	void RemoverTodosElementos() {val.clear(); col.clear(); lprox.clear(); nnz = 0; for (int k = 0; k < nlin; k++) lprim[k] = - 1;}
	void SomarComMatriz(CMatrizEsparsa * m2);
	vetorfloat MultiplicarPorVetorDenso(vetorfloat * vetor);

	int GetNnz() {return nnz;}
	int GetNlin() {return nlin;}
	int GetNcol() {return ncol;}
	double GetElemento(int linha, int coluna);
	int GetValorCol(int l) {return col[l];}
	double GetValorVal(int l) {return val[l];}
	int GetValorLprim(int l) {return lprim[l];}
	int GetValorLprox(int l) {return lprox[l];}
	void SparseMatriz(vetorfloat &nf, vetorint &nr, vetorint &nc );
	void SparseMatriz(double *Bval, int *Bind, int *Bbeg );
	//CMatrizEsparsa operator = (CMatrizEsparsa * m2);
};


// Funções auxiliares para manipulação das matrizes
// ------------------------------------------------
void InserirMatriz(vetorfloat2 * matriz1, int linhaI1, int colunaI1, int linhaF1, int colunaF1, vetorfloat2 * matriz2, int linhaI2, int colunaI2);
void JuntarColunas(vetorfloat2 * matriz1, vetorfloat2 * matriz2);
void IniciaMatriz(vetorfloat2 * matriz1, double valor);
void EyeMatriz(vetorfloat2 * matriz1);
void DimensionarMatriz(vetorfloat2 * matriz1, int linhas, int colunas);
void RemoverColuna(vetorfloat2 * matriz1, int coluna);
vetorfloat2 MultiplciarMatrizes(vetorfloat2 * matrizdiagonal, vetorfloat2 * matriz1);		// Somente para quando a primeira matriz for diagonal
vetorfloat2 SomarMatrizes(vetorfloat2 * matriz1, vetorfloat2 * matriz2);
void MultPorEscalar(vetorfloat2 * matriz1, double escalar);
void SparseMatriz(CMatrizEsparsa * matrizcompleta, vetorint * indexLinha, vetorint * indexColuna, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count );
void SparseMatriz(CMatrizEsparsa * matrizcompleta, vetorint * indexLinha, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count );
void SparseMatriz(vetorfloat2 * matrizcompleta, vetorint * indexLinha, vetorint * indexColuna, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count );
void AlocarLimites(vetorfloat2 * L, vetorint * LimTipo, vetorfloat * LimValor);
// ------------------------------------------------
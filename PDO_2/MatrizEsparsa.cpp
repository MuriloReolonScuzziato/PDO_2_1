#include "MatrizEsparsa.h"

// Funções auxiliares para manipulação das matrizes
// ------------------------------------------------
void InserirMatriz(vetorfloat2 * matriz1, int linhaI1, int colunaI1, int linhaF1, int colunaF1, vetorfloat2 * matriz2, int linhaI2, int colunaI2)
{
	// parte da matriz1 recebe parte da matriz2
	int linhaF2 = linhaF1 - linhaI1 + linhaI2;
	int colunaF2 = colunaF1 - colunaI1 + colunaI2;
	for (int i = 0; i <= linhaF1 - linhaI1; i++)
	{
		for (int j = 0; j <= colunaF1 - colunaI1; j++)
		{
			matriz1->at(i + linhaI1).at(j + colunaI1) = matriz2->at(i + linhaI2).at(j + colunaI2);
		}
	}
}
void JuntarColunas(vetorfloat2 * matriz1, vetorfloat2 * matriz2)
{
	for (size_t i = 0; i < matriz2->size(); i++)
	{
		matriz1->push_back(matriz2->at(i));
	}
}
void IniciaMatriz(vetorfloat2 * matriz1, double valor)
{
	for (size_t ii = 0; ii < matriz1->size(); ii++)
			for (size_t iii = 0; iii < matriz1->at(ii).size(); iii++)
				matriz1->at(ii).at(iii) = valor;
}
void EyeMatriz(vetorfloat2 * matriz1)
{
	for (size_t ii = 0; ii < matriz1->size(); ii++)
			for (size_t iii = 0; iii < matriz1->at(ii).size(); iii++)
				if (ii == iii)
					matriz1->at(ii).at(iii) = 1;
}
void DimensionarMatriz(vetorfloat2 * matriz1, int linhas, int colunas)
{
	matriz1->resize(linhas);
	for (size_t i = 0; i < matriz1->size(); i++)
	{
		matriz1->at(i).resize(colunas);
	}
}
void RemoverColuna(vetorfloat2 * matriz1, int coluna)
{
	for (size_t i = 0; i < matriz1->size(); i++)
	{
		matriz1->at(i).erase(matriz1->at(i).begin() + coluna );
	}
}
vetorfloat2 MultiplciarMatrizes(vetorfloat2 * matrizdiagonal, vetorfloat2 * matriz1)		// Somente para quando a primeira matriz for diagonal
{
	vetorfloat2 resultado;
	DimensionarMatriz(&resultado, int(matrizdiagonal->size()), int(matriz1->at(0).size()));
	for (size_t i = 0; i < matrizdiagonal->size(); i++)
	{
		for (size_t j = 0; j < matriz1->at(0).size(); j++)
		{
			resultado[i][j] = matrizdiagonal->at(i).at(i) * matriz1->at(i).at(j);
		}
	}
	return resultado;
}
vetorfloat2 SomarMatrizes(vetorfloat2 * matriz1, vetorfloat2 * matriz2)
{
	vetorfloat2 resultado;
	DimensionarMatriz(&resultado, int(matriz1->size()), int(matriz1->at(0).size()));
	for (size_t i = 0; i < matriz1->size(); i++)
	{
		for (size_t j = 0; j < matriz1->at(i).size(); j++)
		{
			resultado[i][j] = matriz1->at(i).at(j) + matriz2->at(i).at(j);
		}
	}
	return resultado;
}
void MultPorEscalar(vetorfloat2 * matriz1, double escalar)
{
	for (size_t i = 0; i < matriz1->size(); i++)
	{
		for (size_t j = 0; j < matriz1->at(i).size(); j++)
			matriz1->at(i).at(j) = matriz1->at(i).at(j)*escalar;
	}
}
void SparseMatriz(CMatrizEsparsa * matrizcompleta, vetorint * indexLinha, vetorint * indexColuna, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count )
{
	int l;
	for (int i = 0; i < matrizcompleta->GetNlin(); i++)
	{
		l = matrizcompleta->GetValorLprim(i);
		while ( l != -1 )
		{
			indexLinha->push_back(i + *n_restricoes);
			indexColuna->push_back(matrizcompleta->GetValorCol(l));
			indexValor->push_back(matrizcompleta->GetValorVal(l));
			*count = *count + 1;
			l = matrizcompleta->GetValorLprox(l);
		}
		nnZ->push_back(*count);
	}
}
void SparseMatriz(CMatrizEsparsa * matrizcompleta, vetorint * indexLinha, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count )
{
	int flag = false;
	if (nnZ->size() == 0)
		for (size_t i = 0; i < matrizcompleta->GetNcol() + 1; i++)
			nnZ->push_back(0);
	//if (indexLinha->size() == 0)
	//{
	//	indexLinha->push_back(0);
	//	indexValor->push_back(0);
	//	flag = true;
	//}

	std::vector<int>::iterator it;
	std::vector<double>::iterator itd;
	int l;
	for (int i = 0; i < matrizcompleta->GetNlin(); i++)
	{
		l = matrizcompleta->GetValorLprim(i);
		while ( l != -1 )
		{
			nnZ->at(matrizcompleta->GetValorCol(l) + 1)++;
			// deve-se tb atualizar os elementos de nnZ das colunas a direita desta
			for (int j = matrizcompleta->GetValorCol(l) + 2; j < matrizcompleta->GetNcol() + 1; j++)
				nnZ->at(j)++;
			
			it = indexLinha->begin() + nnZ->at(matrizcompleta->GetValorCol(l) + 1) - 1;
			indexLinha->insert(it, i + *n_restricoes);
			itd = indexValor->begin() + nnZ->at(matrizcompleta->GetValorCol(l) + 1) - 1;
			indexValor->insert(itd, matrizcompleta->GetValorVal(l));
			l = matrizcompleta->GetValorLprox(l);
		}
	}
	//if (flag)		// remover ultimo elemento (adicionado artificialmente para tero iterator begin())
	//{
	//	indexLinha->pop_back();
	//	indexValor->pop_back();
	//}
}
void SparseMatriz(vetorfloat2 * matrizcompleta, vetorint * indexLinha, vetorint * indexColuna, vetorfloat * indexValor, vetorint * nnZ, int * n_restricoes, int * count )
{
	for (size_t i = 0; i < matrizcompleta->size(); i++)
	{
		for (size_t j = 0; j < matrizcompleta->at(i).size(); j++)
		{
			if (matrizcompleta->at(i).at(j) != 0)
			{
				indexLinha->push_back(int(i + *n_restricoes));
				indexColuna->push_back(int(j));
				indexValor->push_back(matrizcompleta->at(i).at(j));
				*count = *count + 1;
			}
		}
		nnZ->push_back(*count);
	}
}
void AlocarLimites(vetorfloat2 * L, vetorint * LimTipo, vetorfloat * LimValor)
{
	for (size_t i = 0; i < L->size(); i++)
	{
		LimTipo->push_back(int(L->at(i).at(0)));
		LimValor->push_back(double(L->at(i).at(1)));
	}
}
// ------------------------------------------------


CMatrizEsparsa::CMatrizEsparsa(void)
{
	val.clear();
	lprim.clear();
	col.clear();
	lprox.clear();
	nlin = 0;
	ncol = 0;
	nnz = 0;
}
CMatrizEsparsa::CMatrizEsparsa(int linhas, int colunas)
{
	nlin = linhas;
	ncol = colunas;
	nnz = 0;
	lprim.resize(linhas);
	for (int k = 0; k < linhas; k++)
		lprim[k] = - 1;
	col.resize(nnz);
	val.resize(nnz);
	lprox.resize(nnz);
}
CMatrizEsparsa::CMatrizEsparsa(int linhas_colunas)
{
	nlin = linhas_colunas;
	ncol = linhas_colunas;
	nnz = 0;
	lprim.resize(linhas_colunas);
	for (int k = 0; k < linhas_colunas; k++)
		lprim[k] = - 1;
	col.resize(nnz);
	val.resize(nnz);
	lprox.resize(nnz);
	for (int k = 0; k < linhas_colunas; k++)
		InserirElemento(k, k, 1);
}
CMatrizEsparsa::CMatrizEsparsa(CMatrizEsparsa &matriz)		// Criar matriz copia
{
	val = matriz.val;
	lprim = matriz.lprim;
	col = matriz.col;
	lprox = matriz.lprox;
	nlin = matriz.nlin;
	ncol = matriz.ncol;
	nnz = matriz.nnz;
}
CMatrizEsparsa::~CMatrizEsparsa(void)
{
	lprim.clear();
	col.clear();
	val.clear();
	lprox.clear();
}
void CMatrizEsparsa::ZerarMatriz(int n_linhas, int n_colunas)
{
	val.clear();
	lprim.resize(n_linhas);
	for (int k = 0; k < n_linhas; k++)
		lprim[k] = - 1;
	col.clear();
	lprox.clear();
	nlin = n_linhas;
	ncol = n_colunas;
	nnz = 0;
}
void CMatrizEsparsa::InserirElemento(int linha, int coluna, double valor)
{
	//if (valor != 0)
	if ( (valor < -1e-8) || (valor > 1e-8) )	// ex.: se for um valor de 1e-10 ou -1e-10 é considerado 0! (zeroTol) - depende da FeasibilityTol do solver
	{
		int temp;
		temp = lprim[linha];
		lprim[linha] = int (col.size());
		col.push_back(coluna);
		val.push_back(valor);
		lprox.push_back(temp);
		nnz++;
	}
}
void CMatrizEsparsa::SubstituirElemento(int linha, int coluna, double valor)		// Insere elemento caso ele ainda não exista
{
	if ( (valor < -1e-8) || (valor > 1e-8) )		// ex.: se for um valor de 1e-10 ou -1e-10 é considerado 0! (zeroTol) - depende da FeasibilityTol do solver
	//if (valor != 0)
	{
		int l;
		l = lprim[linha];
		bool elem_existe = false;
		if (l != -1)		// Se a linha já existe
		{
			while ( l != -1 )		// Percorrer colunas da linha
			{
				if (col[l] == coluna)		// Sobreescrever elemento (caso a coluna já exista)
				{
					val[l] = valor;
					elem_existe = true;
					break;
				}
				l = lprox[l];
			}
			// Inserir elemento inexistente
			if (elem_existe == false)
				InserirElemento(linha, coluna, valor);
			//

		}
		// Inserir elemento inexistente
		else
		{
			InserirElemento(linha, coluna, valor);
		}
		//
	}
	else
		RemoverElemento(linha, coluna);
}
void CMatrizEsparsa::ImprimirLinha(int linha)
{
	int l, viz;
	double valor;
	l = lprim[linha];
	while ( l != -1 )
	{
		viz = col[l];
		valor = val[l];
		cout << viz << " " << valor << "\n";
		l = lprox[l];
	}
}
void CMatrizEsparsa::ImprimirMatriz()
{
	int l, viz;
	double valor;
	for (int k = 0; k < nlin; k++)
	{
		l = lprim[k];
		while ( l != -1 )
		{
			viz = col[l];
			valor = val[l];
			cout << viz << " " << valor << "\n";
			l = lprox[l];
		}
	}
}
void CMatrizEsparsa::ImprimirMatrizVisual()
{
	int l;
	vetorfloat temp_valor;
	vetorint temp_col;
	temp_valor.resize(ncol);
	temp_col.resize(ncol);
	for (int k = 0; k < nlin; k++)
	{
		l = lprim[k];
		while ( l != -1 )
		{
			temp_valor[col[l]] = val[l];
			temp_col[col[l]] = col[l];
			l = lprox[l];
		}
		for (int kk = 0; kk < ncol; kk++)
		{
			if ( kk == temp_col[kk])
				cout << setw(3) << left << temp_valor[kk] << " ";
			else
				cout << setw(3) << left << "0" << " ";
		}
		temp_valor.clear();temp_valor.resize(ncol);
		temp_col.clear();temp_col.resize(ncol);
		cout << "\n";
	}
}
void CMatrizEsparsa::ImprimirMatrizArquivo()
{
	int l;
	vetorfloat temp_valor;
	vetorint temp_col;
	temp_valor.resize(ncol);
	temp_col.resize(ncol);
	ofstream * inFile;
	//inFile = new ofstream( "matriz.txt", ios::app );
	inFile = new ofstream( "matriz.txt", ios::out );
	if ( inFile->is_open() )
	{
		for (int k = 0; k < nlin; k++)
		{
			l = lprim[k];
			while ( l != -1 )
			{
				temp_valor[col[l]] = val[l];
				temp_col[col[l]] = col[l];
				l = lprox[l];
			}
			for (int kk = 0; kk < ncol; kk++)
			{
				if ( kk == temp_col[kk])
					*inFile << setw(15) << std::scientific << setprecision(6) << left << temp_valor[kk] << " " << char(9);
				else
					*inFile << setw(15) << std::scientific << setprecision(6) << left << "0" << " " << char(9);
			}
			temp_valor.clear();temp_valor.resize(ncol);
			temp_col.clear();temp_col.resize(ncol);
			*inFile << "\n";
		}
		inFile->close();
	}
	else
		cout << "Unable to open file";
	delete inFile;
}
void CMatrizEsparsa::JuntarColuna(CMatrizEsparsa * matriz, bool teste)
{
	//if (teste)
	//{
	//	cout << "JuntarColuna()" << endl;
	//	cout << "matriz:" << endl;
	//	cout << "col=" << matriz->ncol << " lin=" << matriz->nlin << " nnz=" << matriz->nnz << " lprim=" << matriz->lprim.size() << endl;
	//	cout << "val=" << matriz->val.size() << " col=" << matriz->col.size() << " lprox=" << matriz->lprox.size() << endl;

	//	cout << "this:" << endl;
	//	cout << "col=" << ncol << " lin=" << nlin << " nnz=" << nnz << " lprim=" << lprim.size() << endl;
	//	cout << "val=" << val.size() << " col=" << col.size() << " lprox=" << lprox.size() << endl;
	//}

	if (ncol != matriz->ncol)
		cout << "Impossivel, matrizes com numero de colunas diferente \n";
	else
	{
		int tamanho_anterior = int(val.size());
		for (size_t k = 0; k < matriz->val.size(); k++)
		{
			val.push_back(matriz->val[k]);
			col.push_back(matriz->col[k]);
			if (matriz->lprox[k] == -1)
				lprox.push_back(-1);
			else
				lprox.push_back(matriz->lprox[k] + tamanho_anterior);
		}
		for (int k = 0; k < matriz->nlin; k++)
		{
			if (matriz->lprim[k] == -1)
				lprim.push_back(-1);
			else
				lprim.push_back(matriz->lprim[k] + nnz);
		}
		nlin += matriz->nlin;
		nnz += matriz->nnz;
	}
}
//void CMatrizEsparsa::RemoverColuna(int coluna)
//{
//	for (int k = 0; k < nlin; k++)
//	{
//		RemoverElemento(k, coluna);
//	}
//	ncol--;
//}
void CMatrizEsparsa::RemoverColuna(int coluna)
{
	int l;
	for (int k = 0; k < nlin; k++)
	{
		RemoverElemento(k, coluna);
		
		l = lprim[k];
		while ( l != -1 )
		{
			if ( col[l] > coluna )
			{
				col[l] -= 1;
			}
			l = lprox[l];
		}
	}
	ncol--;
	// Ineficiente fazer busca exaustiva!!
}
void CMatrizEsparsa::OrdenarLista()
{
	vetorfloat temp_val, temp_val2;
	vetorint temp_lprim, temp_lprox, temp_col, temp_col2;
	int l;
	int count, temp_idisp;
	temp_idisp = 0;
	for (int k = 0; k < nlin; k++)
	{
		l = lprim[k];
		count = 0;
		while ( l != -1)
		{
			temp_col2.push_back(col[l]);
			temp_val2.push_back(val[l]);
			l = lprox[l];
			count++;
		}
		gnome_sort(temp_col2, temp_val2);
		for (size_t kk = 0; kk < temp_col2.size(); kk++)
		{
			temp_col.push_back(temp_col2[kk]);
			temp_val.push_back(temp_val2[kk]);
		}
		temp_col2.clear();
		temp_val2.clear();
		for (int kk = 0; kk < count; kk++)
		{
			if (kk == 0)
				temp_lprox.push_back( - 1 );
			else
				temp_lprox.push_back( kk - 1 + temp_idisp);
		}
		temp_idisp += count;
		if (count == 0)
			temp_lprim.push_back( -1);
		else
			temp_lprim.push_back(temp_idisp - 1);
	}
	val = temp_val;
	lprim = temp_lprim;
	lprox = temp_lprox;
	col = temp_col;
}
void CMatrizEsparsa::RemoverElemento(int linha, int coluna)
{
	// quando remove-se um elemento nnz é decrementado, porém os vetores col, val,... continuam com o valor, porem ele n é acessado!
	// se o elemento desses vetores fosse removido os indices teriam q ser atualizados
	int l, * anterior;
	l = lprim[linha];
	anterior = &lprim[linha];
	while ( l != -1 )
	{
		if ( col[l] == coluna )
		{
			col[l] = -1;
			val[l] = -1;
			*anterior = lprox[l];
			lprox[l] = -1;
			nnz--;
			break;
		}
		anterior = &lprox[l];
		l = lprox[l];
	}
}
double CMatrizEsparsa::GetElemento(int linha, int coluna)
{
	int l;
	double elemento = 0;
	l = lprim[linha];
	while ( l != -1 )
	{
		if ( col[l] == coluna )
		{
			elemento = val[l];
			break;
		}
		l = lprox[l];
	}
	return elemento;
}
//CMatrizEsparsa CMatrizEsparsa::MultiplicarMatrizes(CMatrizEsparsa * m1, CMatrizEsparsa * m2)
//{
//	m1->OrdenarLista();
//	m2->OrdenarLista();
//	CMatrizEsparsa mr(m1->nlin, m2->ncol);
//	if (m1->GetNcol() != m2->nlin)
//		cout << "Matrizes com dimensões incompatíveis!! \n";
//	else
//	{
//		vetorint temp_posicoes;
//		double elemento;
//		int l;
//		for (int k = 0; k < m1->nlin; k++)
//		{
//			temp_posicoes.clear();
//			l = m1->lprim[k];
//			while ( l != -1 )
//			{
//				temp_posicoes.push_back(m1->col[l]);
//				l = m1->lprox[l];
//			}
//			for (int kk = 0; kk < m2->ncol; kk++)
//			{
//				elemento = 0;
//				for (size_t kkk = 0; kkk < temp_posicoes.size(); kkk++)
//				{
//					if (m2->GetElemento(temp_posicoes[kkk], kk) != 0)
//					{
//						elemento += m1->GetElemento(k, temp_posicoes[kkk]) * m2->GetElemento(temp_posicoes[kkk], kk);
//					}
//					else
//					{
//						elemento += 0;
//					}
//				}
//				if (elemento != 0)
//					mr.InserirElemento(k, kk, elemento);
//			}
//		}
//	}
//	return mr;
//}
void CMatrizEsparsa::MultiplicarPorMatriz(CMatrizEsparsa * m2)
{
	CMatrizEsparsa mr(nlin, m2->ncol);
	// Multiplicação simbólica
	vetorint ptr;
	ptr.resize(m2->ncol);
	for (int i = 0; i < m2->ncol; i++)
		ptr[i] = -1;
	int nloc_C = -1;
	int jp, j, kp, k, l;
	for (int i = 0; i < nlin; i++)
	{
		jp = lprim[i];
		while (jp != -1)
		{
			j = col[jp];
			kp = m2->lprim[j];
			while (kp != -1 )
			{
				k = m2->col[kp];
				if ( ptr[k] == -1 )
				{
					nloc_C++;
					ptr[k] = nloc_C;
					l = mr.lprim[i];
					mr.lprim[i] = nloc_C;
					mr.lprox.push_back(l);
					mr.col.push_back(k);
					//mr.lprox[nloc_C] = l;
					//mr.col[nloc_C] = k;
				}
				kp = m2->lprox[kp];
			}
			jp = lprox[jp];
		}
		for (int ii = 0; ii < m2->ncol; ii++)
			ptr[ii] = -1;
		kp = mr.lprim[i];
		while (kp != -1)
		{
			ptr[mr.col[kp]] = -1;
			kp = mr.lprox[kp];
		}
	}
	// Multiplicação numérica
	double alfa;
	int cp;
	for (int ii = 0; ii < m2->ncol; ii++)
			ptr[ii] = -1;
	for (int i = 0; i < nlin; i++)
	{
		jp = mr.lprim[i];
		while (jp != -1 )
		{
			mr.val.push_back(0);
			ptr[mr.col[jp]] = jp;
			jp = mr.lprox[jp];
		}
		jp = lprim[i];
		while (jp != -1 )
		{
			j = col[jp];
			alfa = val[jp];
			kp = m2->lprim[j];
			while (kp != -1)
			{
				k = m2->col[kp];
				cp = ptr[k];
				mr.val[cp] = mr.val[cp] + alfa * m2->val[kp];
				kp = m2->lprox[kp];
			}
			jp = lprox[jp];
		}
	}
	// Receber resultado
	val = mr.val;
	lprim = mr.lprim;
	col = mr.col;
	lprox = mr.lprox;
	nlin = mr.nlin;
	ncol = mr.ncol;
	nnz = mr.val.size();
	//nnz = mr.nnz;
}
void CMatrizEsparsa::MultiplicarPorMatriz(Eigen::MatrixXd &m2)
{
	CMatrizEsparsa mr(nlin, m2.cols());
	// Multiplicação numérica
	double valor = 0;
	int l;
	for (int k = 0; k < nlin; k++)
	{
		for (int kk = 0; kk < ncol; kk++)
		{
			valor = 0;
			l = lprim[k];
			while ( l != -1 )
			{
				valor += val[l] * m2(col[l], kk);
				l = lprox[l];
			}
			mr.SubstituirElemento(k, kk, valor); 
		}
	}
	// Receber resultado
	val = mr.val;
	lprim = mr.lprim;
	col = mr.col;
	lprox = mr.lprox;
	nlin = mr.nlin;
	ncol = mr.ncol;
	nnz = mr.val.size();
}
void CMatrizEsparsa::MultiplicarPorEscalar(double valor)
{
	for (size_t k = 0; k < val.size(); k++)
		val[k] *= valor;
}
void CMatrizEsparsa::InserirMatriz(int linha_i, int coluna_i, int linha_f, int coluna_f, CMatrizEsparsa * m2, int linha2_i, int coluna2_i)
{
	int delta_linha = linha_f - linha_i;
	int delta_coluna = coluna_f - coluna_i;
	int l, l2;
	if ( (delta_linha + 1 > nlin) || (delta_coluna + 1> ncol) || ((linha2_i + delta_linha + 1) > m2->nlin) || ((coluna2_i + delta_coluna + 1) > m2->ncol) )
		cout << " Valores de entrada incompatíveis!! \n";
	else
	{
		for (int k = linha_i; k <= linha_f; k++)
		{
			l = lprim[k];
			while ( l != -1 )
			{
				if ( (col[l] >= coluna_i) && (col[l] <= coluna_f) )
				{
					l2 = l;
					l = lprox[l];
					RemoverElemento(k, col[l2]);
				}
				else
					l = lprox[l];
			}
		}
		for (int k = linha2_i; k <= linha2_i + delta_linha; k++)
		{
			l = m2->lprim[k];
			while ( l != -1 )
			{
				if ( (m2->col[l] <= coluna2_i + delta_coluna) && (m2->col[l] >= coluna2_i) )
					InserirElemento(k - linha2_i + linha_i, m2->col[l] - coluna2_i + coluna_i, m2->val[l]);
				l = m2->lprox[l];
			}
		}
	}
}
void CMatrizEsparsa::InserirMatriz(int linha_i, int coluna_i, int linha_f, int coluna_f, Eigen::MatrixXd &m2, int linha2_i, int coluna2_i)
{
	int delta_linha = linha_f - linha_i;
	int delta_coluna = coluna_f - coluna_i;
	int l, l2;
	if ( (delta_linha + 1 > nlin) || (delta_coluna + 1> ncol) || ((linha2_i + delta_linha + 1) > m2.rows()) || ((coluna2_i + delta_coluna + 1) > m2.cols()) )
		cout << " Valores de entrada incompatíveis!! \n";
	else
	{
		//
		for (int k = linha_i; k <= linha_f; k++)
		{
			l = lprim[k];
			while ( l != -1 )
			{
				if ( (col[l] >= coluna_i) && (col[l] <= coluna_f) )
				{
					l2 = l;
					l = lprox[l];
					RemoverElemento(k, col[l2]);
				}
				else
					l = lprox[l];
			}
		}
		for (int k = linha2_i; k <= linha2_i + delta_linha; k++)
			for (int kk = coluna2_i; kk <= coluna2_i + delta_coluna; kk++)
				InserirElemento(k - linha2_i + linha_i, kk - coluna2_i + coluna_i, m2(k,kk));
		//
		// ou
		//for (int k = linha2_i; k <= linha2_i + delta_linha; k++)
		//	for (int kk = coluna2_i; kk <= coluna2_i + delta_coluna; kk++)
		//		SubstituirElemento(k - linha2_i + linha_i, kk - coluna2_i + coluna_i, m2(k,kk));
	}
}
//void CMatrizEsparsa::Transpor()
//{
//	CMatrizEsparsa temp_m1;
//	temp_m1.ncol = nlin;
//	temp_m1.nlin = ncol;
//	temp_m1.nnz = nnz;
//
//
//	//lprim.clear();
//	//lprim.resize(ncol);
//	//for (int k = 0; k < ncol; k++)
//	//	lprim[k] = 0;
//}
//CMatrizEsparsa CMatrizEsparsa::operator = (CMatrizEsparsa * m2)
//{
//	CMatrizEsparsa temp();
//
//	//temp.x = x + param.x;
//	//temp.y = y + param.y;
//	return (temp);
//}

void CMatrizEsparsa::SomarComMatriz(CMatrizEsparsa * m2)
{
	if ( (ncol != m2->ncol) || (nlin != m2->nlin) ) 
		cout << "Impossivel, matrizes com dimensões diferentes \n";
	else
	{
		int l;
		for (int k = 0; k < nlin; k++)
		{
			l = m2->lprim[k];
			while ( l != -1 )
			{
				if ( GetElemento(k, m2->col[l]) != 0)
					SubstituirElemento(k, m2->col[l], GetElemento(k, m2->col[l]) + m2->val[l]);
					//InserirElemento(k, m2->col[l], GetElemento(k, m2->col[l]) + m2->val[l]);
				else
					//SubstituirElemento(k, m2->col[l], m2->val[l]);
					InserirElemento(k, m2->col[l], m2->val[l]);
				l = m2->lprox[l];
			}
		}
	}
}
vetorfloat CMatrizEsparsa::MultiplicarPorVetorDenso(vetorfloat * vetor)
{
	vetorfloat resultado;
	if ( ncol != vetor->size())
		cout << " Vetores de tamanhos incompatíveis!" << endl;
	else
	{
		resultado.resize(nlin);
		int l, j;
		for (int k = 0; k < nlin; k++)
		{
			resultado[k] = 0;
			l = lprim[k];
			while ( l != -1)
			{
				j = col[l];
				resultado[k] += val[l] * vetor->at(j);
				l = lprox[l];
			}
		}
	}
	return resultado;
}
void CMatrizEsparsa::SparseMatriz(vetorfloat &nf, vetorint &nr, vetorint &nc )
{
	nf.resize(nnz);
	nc.resize(ncol + 1);
	nr.resize(nnz);

	vetorint rr;
	rr.resize(nnz);
	int l;
	for (int k = 0, i = 0; i < nlin; i++)
	{
		l = lprim[i];
		while ( l != -1 )
		{
			rr[l] = i;
			l = lprox[l];
		}
	}

	for (int i = 0; i < nnz; i++)
		nc[col[i]+1]++;
	for (int i = 1; i <= ncol; i++)
		nc[i] += nc[i-1];
	vetorint nn = nc;

	for (int i = 0; i < nnz; i++)
	{
		int x = nn[col[i]]++;
		nf[x] = val[i];
		nr[x] = rr[i];
	}
}
void CMatrizEsparsa::SparseMatriz(double *Bval, int *Bind, int *Bbeg )
{
	for (int i = 0; i < ncol + 1; i++)
		Bbeg[i] = 0;

	int * rr;
	int * nn;
	rr = new int[nnz];
	nn = new int[ncol+1];
	int l;
	for (int k = 0, i = 0; i < nlin; i++)
	{
		l = lprim[i];
		while ( l != -1 )
		{
			rr[l] = i;
			l = lprox[l];
		}
	}

	for (int i = 0; i < nnz; i++)
		Bbeg[col[i]+1]++;
	for (int i = 1; i <= ncol; i++)
		Bbeg[i] += Bbeg[i-1];

	for (int i = 0; i < ncol+1; i++)
		nn[i] = Bbeg[i];

	for (int i = 0; i < nnz; i++)
	{
		int x = nn[col[i]]++;
		Bval[x] = val[i];
		Bind[x] = rr[i];
	}
	delete [] rr;
	delete [] nn;
}
// verificar o somar, multiplicar e inserir matriz
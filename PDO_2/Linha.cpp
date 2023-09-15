#include "Linha.h"

CLinha::CLinha(string nome_linha_a, CBarra * de_barra_a, CBarra * para_barra_a, double reatancia_a, double capacidade_a, int ind_de_barra_a, int ind_para_barra_a)
{
	nome_linha = nome_linha_a;
	de_barra = de_barra_a;
	para_barra = para_barra_a;
	ind_de_barra = ind_de_barra_a;
	ind_para_barra = ind_para_barra_a;
	reatancia = reatancia_a;
	capacidade = capacidade_a;
}


CLinha::~CLinha(void)
{

}

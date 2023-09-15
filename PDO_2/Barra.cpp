#include "Barra.h"

CBarra::CBarra(	string nome_barra_a, int ident)
{
	nome_barra = nome_barra_a;
	ident_barra = ident;
}
CBarra::~CBarra(void)
{
	demandasPtr.clear();
	hidrosPtr.clear();
	termosPtr.clear();
}
void CBarra::AddDemanda(CDemanda * demanda_end)
{
	demandasPtr.push_back(demanda_end);
}
void CBarra::AddHidro(CHidreletrica * hidro_end)
{
	hidrosPtr.push_back(hidro_end);
}
void CBarra::AddTermo(CTermeletrica * termo_end)
{
	termosPtr.push_back(termo_end);
}
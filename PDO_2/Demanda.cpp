#include "Demanda.h"

CDemanda::CDemanda(int n_demanda_a)
{
	ostringstream convert;   // stream used for the conversion
	convert << n_demanda_a;
	string inicial = "D";
	nome_demanda = inicial+convert.str();

	//D_1 = D_1_a;
	//D_2 = D_2_a;
}
CDemanda::~CDemanda(void)
{
}
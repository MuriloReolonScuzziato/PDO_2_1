#include "Termeletrica.h"


CTermeletrica::CTermeletrica(string nome_usina_a, int ident, double pmin_a, double pmax_a, int t_down_a, int t_up_a, double rampa_down_a, double rampa_up_a, vetorfloat coef_custo_oper_a, vetorfloat coef_custo_partida_a)
{
	nome_usina = nome_usina_a;
	ident_usina = ident;
	pmin = pmin_a;
	pmax = pmax_a;
	rampa_down = rampa_down_a;
	rampa_up = rampa_up_a;
	t_down = t_down_a;
	t_up = t_up_a;
	coef_custo_oper = coef_custo_oper_a;
	coef_custo_partida = coef_custo_partida_a;
}
CTermeletrica::~CTermeletrica(void)
{
}
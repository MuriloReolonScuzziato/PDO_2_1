// Define caracteristicas da classe usinas (grupos)
#include "Unidades.h"

CUnidades::CUnidades(int n_usina_a, int grupo_a, int n_unidades_a, double qmin_a, double qmax_a, double pmin_a, double pmax_a, double coef_perdas_hid_a, vetorfloat coef_rend_a)
{
	n_usina = n_usina_a;
	grupo = grupo_a;
	n_unidades = n_unidades_a;
	qmin = qmin_a;
	qmax = qmax_a;
	pmin = pmin_a;
	pmax = pmax_a;
	coef_perdas_hid = coef_perdas_hid_a;
	coef_rend = coef_rend_a;
}
CUnidades::~CUnidades(void)
{
}
double CUnidades::PotenciaGerada(double v_a, double q_a, double d_s_a)
{
	double ph = 1e10;
	for (size_t ii = 0; ii < coef_phg.size(); ii++)
		ph = std::fmin(ph, coef_phg[ii][0] + coef_phg[ii][1]*v_a + coef_phg[ii][2]*q_a + coef_phg[ii][3]*d_s_a);
	return ( ph );
}
//double CUnidades::GradPotenciaGeradaV(double v_a, double q_a, double d_a)
//{
//	return ((coef_phg.at(1) + coef_phg.at(4)*q_a + coef_phg.at(5)*d_a + 2*coef_phg.at(7)*v_a));
//}
//double CUnidades::HessPotenciaGeradaVV()
//{
//	return (2*coef_phg.at(7));
//}
//double CUnidades::HessPotenciaGeradaVQ()
//{
//	return (coef_phg.at(4));
//}
//double CUnidades::HessPotenciaGeradaVD()
//{
//	return (coef_phg.at(5));
//}
//double CUnidades::GradPotenciaGeradaQ(double v_a, double q_a, double d_a)
//{
//	return ((coef_phg.at(2) + coef_phg.at(4)*v_a + coef_phg.at(6)*d_a + 2*coef_phg.at(8)*q_a));
//}
//double CUnidades::HessPotenciaGeradaQV()
//{
//	return ((coef_phg.at(4)));
//}
//double CUnidades::HessPotenciaGeradaQQ()
//{
//	return ((2*coef_phg.at(8)));
//}
//double CUnidades::HessPotenciaGeradaQD()
//{
//	return ((coef_phg.at(6)));
//}
//double CUnidades::GradPotenciaGeradaD(double v_a, double q_a, double d_a)
//{
//	return ((coef_phg.at(3) + coef_phg.at(5)*v_a + coef_phg.at(6)*q_a + 2*coef_phg.at(9)*d_a));
//}
//double CUnidades::HessPotenciaGeradaDV()
//{
//	return ((coef_phg.at(5)));
//}
//double CUnidades::HessPotenciaGeradaDQ()
//{
//	return ((coef_phg.at(6)));
//}
//double CUnidades::HessPotenciaGeradaDD()
//{
//	return ((2*coef_phg.at(9)));
//}


//// fazer hess daqui pra baixo
//double CUnidades::PotenciaGerada(double v_a, double q_a, double d_a, double s_a)
//{
//	return ((coef_phg.at(0) + coef_phg.at(1)*v_a + coef_phg.at(2)*q_a + coef_phg.at(3)*(d_a - s_a) + coef_phg.at(4)*v_a*q_a + coef_phg.at(5)*v_a*(d_a - s_a) + coef_phg.at(6)*q_a*(d_a - s_a) + coef_phg.at(7)*pow(v_a,2) + coef_phg.at(8)*pow(q_a,2) + coef_phg.at(9)*pow((d_a - s_a),2)));
//}
//double CUnidades::GradPotenciaGeradaV(double v_a, double q_a, double d_a, double s_a)
//{
//	return ((coef_phg.at(1) + coef_phg.at(4)*q_a + coef_phg.at(5)*(d_a - s_a) + 2*coef_phg.at(7)*v_a));
//}
//double CUnidades::GradPotenciaGeradaQ(double v_a, double q_a, double d_a, double s_a)
//{
//	return ((coef_phg.at(2) + coef_phg.at(4)*v_a + coef_phg.at(6)*(d_a - s_a) + 2*coef_phg.at(8)*q_a));
//}
//double CUnidades::GradPotenciaGeradaD(double v_a, double q_a, double d_a, double s_a)
//{
//	return ((coef_phg.at(3) + coef_phg.at(5)*v_a + coef_phg.at(6)*q_a + 2*coef_phg.at(9)*(d_a - s_a)));
//}
//double CUnidades::GradPotenciaGeradaS(double v_a, double q_a, double d_a, double s_a)
//{
//	return ((- coef_phg.at(3) - coef_phg.at(5)*v_a - coef_phg.at(6)*q_a - 2*coef_phg.at(9)*(d_a - s_a)));
//}
//double CUnidades::PotenciaGerada(double v_a, double q_a, double d_a)
//{
//	if (q_a >= qmin)	// Existe a possibilidade de q_a estar fora dos limites da aproximação da função
//		return (coef_phg.at(0) + coef_phg.at(1)*v_a + coef_phg.at(2)*q_a + coef_phg.at(3)*d_a + coef_phg.at(4)*v_a*q_a + coef_phg.at(5)*v_a*d_a + coef_phg.at(6)*q_a*d_a + coef_phg.at(7)*pow(v_a,2) + coef_phg.at(8)*pow(q_a,2) + coef_phg.at(9)*pow(d_a,2));
//	else
//		return pmin;
//}
//double CUnidades::GradPotenciaGeradaV(double v_a, double q_a, double d_a)
//{
//	if (q_a >= qmin)
//		return (coef_phg.at(1) + coef_phg.at(4)*q_a + coef_phg.at(5)*d_a + 2*coef_phg.at(7)*v_a);
//	else
//		return 0;	// Retornar 0 para o gradiente quando a unidade está desligada
//}
//double CUnidades::GradPotenciaGeradaQ(double v_a, double q_a, double d_a)
//{
//	if (q_a >= qmin)
//		return (coef_phg.at(2) + coef_phg.at(4)*v_a + coef_phg.at(6)*d_a + 2*coef_phg.at(8)*q_a);
//	else
//		return 0;
//}
//double CUnidades::GradPotenciaGeradaD(double v_a, double q_a, double d_a)
//{
//	if (q_a >= qmin)
//		return (coef_phg.at(3) + coef_phg.at(5)*v_a + coef_phg.at(6)*q_a + 2*coef_phg.at(9)*d_a);
//	else
//		return 0;
//}
//double CUnidades::PotenciaGerada(double v_a, double q_a, double d_a, double s_a)
//{
//	if (q_a >= qmin)	// Existe a possibilidade de q_a estar fora dos limites da aproximação da função
//		return (coef_phg.at(0) + coef_phg.at(1)*v_a + coef_phg.at(2)*q_a + coef_phg.at(3)*(d_a - s_a) + coef_phg.at(4)*v_a*q_a + coef_phg.at(5)*v_a*(d_a - s_a) + coef_phg.at(6)*q_a*(d_a - s_a) + coef_phg.at(7)*pow(v_a,2) + coef_phg.at(8)*pow(q_a,2) + coef_phg.at(9)*pow((d_a - s_a),2));
//	else
//		return pmax;
//}
//double CUnidades::GradPotenciaGeradaV(double v_a, double q_a, double d_a, double s_a)
//{
//	if (q_a >= qmin)
//		return (coef_phg.at(1) + coef_phg.at(4)*q_a + coef_phg.at(5)*(d_a - s_a) + 2*coef_phg.at(7)*v_a);
//	else
//		return 0;
//}
//double CUnidades::GradPotenciaGeradaQ(double v_a, double q_a, double d_a, double s_a)
//{
//	if (q_a >= qmin)
//		return (coef_phg.at(2) + coef_phg.at(4)*v_a + coef_phg.at(6)*(d_a - s_a) + 2*coef_phg.at(8)*q_a);
//	else
//		return 0;
//}
//double CUnidades::GradPotenciaGeradaD(double v_a, double q_a, double d_a, double s_a)
//{
//	if (q_a >= qmin)
//		return (coef_phg.at(3) + coef_phg.at(5)*v_a + coef_phg.at(6)*q_a + 2*coef_phg.at(9)*(d_a - s_a));
//	else
//		return 0;
//}
//double CUnidades::GradPotenciaGeradaS(double v_a, double q_a, double d_a, double s_a)
//{
//	if (q_a >= qmin)
//		return (- coef_phg.at(3) - coef_phg.at(5)*v_a - coef_phg.at(6)*q_a - 2*coef_phg.at(9)*(d_a - s_a));
//	else
//		return 0;
//}
//double CUnidades::PotenciaGeradaUnid(double queda_liquida_a, double q_a, double rendimento_hidraulico_a)
//{
//	return (9.81e-3*rendimento_hidraulico_a*queda_liquida_a*q_a)*n_unidades;
//}
double CUnidades::QuedaLiquida(double queda_bruta_a, double q_a)
{
	return queda_bruta_a - pow(q_a,2);
}
double CUnidades::RendimentoHidraulico(double queda_liquida_a, double q_a)
{
	return (coef_rend.at(0)+coef_rend.at(1)*q_a+coef_rend.at(2)*queda_liquida_a+coef_rend.at(3)*q_a*queda_liquida_a+coef_rend.at(4)*pow(q_a,2)+coef_rend.at(5)*pow(queda_liquida_a,2));
}

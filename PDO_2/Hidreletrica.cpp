// Define caracteristicas da classe hidreletrica
#include "Hidreletrica.h"

//int CHidreletrica::n_hidreletricas = 0;

CHidreletrica::CHidreletrica(string nome_usina_a, int ident, double vmin_a, double vmax_a, double smin_a, double smax_a, int influe_vert_a, int n_grupos_a, int usina_jusante_a, int tempo_viagem_a, vetorfloat coef_mont_a, vetorfloat coef_jus_a)
{
	nome_usina = nome_usina_a;
	ident_usina = ident;
	vmin = vmin_a;
	vmax = vmax_a;
	smin = smin_a;
	smax = smax_a;
	influe_vert = influe_vert_a;
	n_grupos = n_grupos_a;
	usina_jusante = usina_jusante_a;
	tempo_viagem = tempo_viagem_a;
	coef_mont = coef_mont_a;
	coef_jus = coef_jus_a;
	//n_hidreletricas++;
	Vmin.resize(0);
	Vmax.resize(0);
	usinas_montante.resize(0);
}
CHidreletrica::~CHidreletrica(void)
{
	//n_hidreletricas--;
	grupoVtr.clear();
}
//void CHidreletrica::add_grupo(int usina_a, int grupo_a, int n_unidades_a, double qmin_a, double qmax_a, double pmin_a, double pmax_a, double coef_perdas_hid_a, vetorfloat coef_rend_a)
//{
//	grupo.push_back(new CUnidades(usina_a, grupo_a, n_unidades_a, qmin_a, qmax_a, pmin_a, pmax_a, coef_perdas_hid_a, coef_rend_a));
//}
double CHidreletrica::NivelMontante(double v_a)
{
	double nivel_montante = 0;
	for (size_t i = 0; i < coef_mont.size(); i++)
		 nivel_montante = nivel_montante + coef_mont.at(i)*pow(v_a,int(i));
	return nivel_montante;
}
double CHidreletrica::NivelJusante(double d_a)
{
	double nivel_jusante = 0;
	for (size_t i = 0; i < coef_jus.size(); i++)
		 nivel_jusante = nivel_jusante + coef_jus.at(i)*pow(d_a,int(i));
	return nivel_jusante;
}
double CHidreletrica::QuedaBruta(double v_a, double d_a)
{
	return (NivelMontante(v_a) - NivelJusante(d_a));
}
double CHidreletrica::GetDmax()
{
	double dmax = 0;
	for (int j = 0; j < n_grupos; j++)
		dmax += grupoVtr[j].GetQmax() * grupoVtr[j].GetNUnidades();
	dmax += smax;
	return dmax;
}
bool CHidreletrica::FcmIsCte()
{
	double soma = 0;
	for (size_t i = 1; i < coef_mont.size(); i++)
		soma += coef_mont[i];
	if (soma != 0)
		return false;
	else
		return true;
}
bool CHidreletrica::FcjIsCte()
{
	double soma = 0;
	for (size_t i = 1; i < coef_jus.size(); i++)
		soma += coef_jus[i];
	if (soma != 0)
		return false;
	else
		return true;
}
double CHidreletrica::GetGeracaoOtima(double v_a, double Q_a, double de_a, int grupo)
{
	// Essa função retorna o valor 'real' da função, sem linearização, portanto se a função fosse representada
	// por um polinomio de ordem alta a única diferença seria que
	// essa aproximação não representa o 'gap' de geração quando se tem numero de unidades diferentes operando
	// por exemplo, qmin = 80 e qmax = 100. Entre 100 e 160 não teria-se geração para esse usina, porem aqui a sol~ução é operar com 2 unidades em baixa eficiencia.

	double ph = -1e9;
	double cm, cj, h, ro;
	Q_a = floor(Q_a*1e9)/1e9;		// para evitar imprecisão numérica no valor de Q_a (problema ao comparar com Qmin e Qmax)
	for (int n = 1; n <= grupoVtr[grupo].GetNUnidades(); n++)
	{
		if ((Q_a >= grupoVtr[grupo].GetQmin()) && (Q_a <= grupoVtr[grupo].GetQmax()*n))		// a regiao que cada combinação pode gerar vai de qmin até n*qmax, isso exlcui regiões proibidas entre as unidades (aproximação!)
		{
			cm = NivelMontante( v_a );
			if (cm < 0)
				cm = 0;
			cj = NivelJusante ( Q_a + de_a );
			h = cm - cj - grupoVtr[grupo].GetCoefPerdasHid()*pow( Q_a/n, 2);
			ro = grupoVtr[grupo].RendimentoHidraulico( h, Q_a/n);
			ph = std::fmax(ph, 0.00981*ro*h*Q_a);
		}
	}
	return ( ph );
}
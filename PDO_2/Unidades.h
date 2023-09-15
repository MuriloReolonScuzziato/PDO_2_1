#pragma once
#include <string>
using std::string;
//#define CRH 6
//#define CPHG 10

#include <cmath>
using std::floor;
//using std::abs;

#include <vector>
using std::vector;

typedef vector<double> vetorfloat;
typedef vector<vetorfloat> vetorfloat2;
typedef vector<int> vetorint;

class CUnidades
{
	int n_usina, grupo, n_unidades;
	double qmin, qmax, pmin, pmax, coef_perdas_hid;
	//double coef_rend [CRH], coef_phg [CPHG];
	vetorfloat coef_rend;
	vetorfloat2 coef_phg;

public:
	CUnidades(int, int, int, double, double, double, double, double, vetorfloat);
	~CUnidades(void);
	
	
	double PotenciaGerada(double v_a, double q_a, double d_s_a); //{return (coef_phg.at(0) + coef_phg.at(1)*v_a + coef_phg.at(2)*q_a + coef_phg.at(3)*d_a);}
	//double PotenciaGerada(double v_a, double q_a, double d_a) {return (coef_phg.at(2)*q_a);}		// FPH constante
	double CoefPhg(int appr) {return coef_phg[appr][0];}
	double CoefPhgV(int appr) {return coef_phg[appr][1];}
	double CoefPhgQ(int appr) {return coef_phg[appr][2];}
	double CoefPhgD(int appr) {return coef_phg[appr][3];}
	double CoefPhgS(int appr) {return - coef_phg[appr][3];}
	int GetNappFPH() {return int (coef_phg.size());}
	
	// A função abaixo n precisa?!?! pois o s já está considerado no d?!?
	//double PotenciaGerada(double v_a, double q_a, double de_a, double s_a) {return (coef_phg.at(0) + coef_phg.at(1)*v_a + coef_phg.at(2)*q_a + coef_phg.at(3)*(d_a - s_a));}
	
	//double PotenciaGerada(double, double, double);
	//double GradPotenciaGeradaV(double, double, double);
	//double HessPotenciaGeradaVV();
	//double HessPotenciaGeradaVQ();
	//double HessPotenciaGeradaVD();
	//double GradPotenciaGeradaQ(double, double, double);
	//double HessPotenciaGeradaQV();
	//double HessPotenciaGeradaQQ();
	//double HessPotenciaGeradaQD();
	//double GradPotenciaGeradaD(double, double, double);
	//double HessPotenciaGeradaDV();
	//double HessPotenciaGeradaDQ();
	//double HessPotenciaGeradaDD();

	//double PotenciaGerada(double, double, double, double);			// Para usinas em que o vertimento não influencia no canal de fuga
	//double GradPotenciaGeradaV(double, double, double, double);
	//double GradPotenciaGeradaQ(double, double, double, double);
	//double GradPotenciaGeradaD(double, double, double, double);
	//double GradPotenciaGeradaS(double, double, double, double);

	double QuedaLiquida(double , double);
	double RendimentoHidraulico(double, double);

	int GetNUsina() {return n_usina;}
	int GetGrupo() {return grupo;}
	int GetNUnidades() {return n_unidades;}
	double GetQmin() {return qmin;}
	double GetQmax() {return qmax;}
	double GetPmin() {return pmin;}
	double GetPmax() {return pmax;}
	double GetCoefPerdasHid() {return coef_perdas_hid;}

	void SetCoefPhg(vetorfloat2 coef_phg_a) {coef_phg = coef_phg_a;};
};

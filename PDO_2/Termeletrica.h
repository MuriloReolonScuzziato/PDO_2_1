#pragma once
#include <string>
using std::string;

#include <vector>
using std::vector;
typedef vector<double> vetorfloat;

class CTermeletrica
{
	string nome_usina;
	int ident_usina;
	double pmin, pmax, rampa_down, rampa_up, pt0;
	int t_down, t_up, u0, x0;
	vetorfloat coef_custo_oper, coef_custo_partida;
	vetorfloat coef_a0_aprox_linear, coef_a1_aprox_linear;

	vetorfloat ptc;		// 1xnumero_cortes_total vetor com o valor da potencia de corte de cada usina (pontos onde são criados os PC)

public:
	CTermeletrica(string nome_usina_a, int ident, double pmin_a, double pmax_a, int t_down_a, int t_up_a, double rampa_down_a, double rampa_up_a, vetorfloat coef_custo_oper_a, vetorfloat coef_custo_partida_a);
	~CTermeletrica(void);

	double CustoOperacao(double pt_a, double u_a) {return (coef_custo_oper[0]*u_a + coef_custo_oper[1]*pt_a + coef_custo_oper[2]*pow(pt_a,2));}
	double GradCustoOperacao(double pt_a) {return (coef_custo_oper[1] + 2*coef_custo_oper[2]*pt_a);}
	double GradCustoOperacao() {return (coef_custo_oper[0]);}
	double HessCustoOperacao() {return (2*coef_custo_oper[2]);}
	//double CustoPartida(double u_a, double u_a0) {return (coef_custo_partida[1]*u_a*(1 - u_a0));}
	double GetCoefCustoPartida() {return (coef_custo_partida[1]);}
	//double CustoPartida(double x_a0, double u_a, double u_a0, double delta_t) {return (coef_custo_partida[0]*u_a*(1 - exp(-(x_a0*delta_t)/coef_custo_partida[2])) + coef_custo_partida[1]*u_a*(1 - u_a0));}
	//double GradCustoPartidaX0(double x_a0, double u_a, double delta_t) {return (coef_custo_partida[0]*u_a*(-(x_a0*delta_t)/coef_custo_partida[2])*(- exp(-(x_a0*delta_t)/coef_custo_partida[2])));}
	//double GradCustoPartidaU(double x_a0, double u_a0, double delta_t) {return (coef_custo_partida[0]*(1 - exp(-(x_a0*delta_t)/coef_custo_partida[2])) + coef_custo_partida[1]*(1 - u_a0));}
	//double GradCustoPartidaU0(double u_a) {return (- coef_custo_partida[1]*u_a);}
	//double CustoPartida(int x_a0, int u_a, int u_a0, double delta_t) {return (0);}
	//double GradCustoPartidaX0(int x_a0, int u_a, double delta_t) {return (0);}
	//double GradCustoPartidaU(int x_a0, int u_a0, double delta_t) {return (0);}
	//double GradCustoPartidaU0(int u_a) {return (0);}

	string GetNomeUsina() {return nome_usina;}
	int GetIdentUsina() {return ident_usina;}
	double GetPmin() {return pmin;}
	double GetPmax() {return pmax;}
	double GetRampaDown() {return rampa_down;}
	double GetRampaUp() {return rampa_up;}
	int GetTDown() {return t_down;}
	int GetTUp() {return t_up;}
	double GetPt0() {return pt0;}
	int GetU0() {return u0;}
	int GetX0() {return x0;}
	double GetCoefCustoOper(int a) {return coef_custo_oper[a];}
	void SetPt0(double pt0_a) {pt0 = pt0_a;}
	void SetU0(int u0_a) {u0 = u0_a;}
	void SetX0(int x0_a) {x0 = x0_a;}
	void AdicionarReta(double a0, double a1) {coef_a0_aprox_linear.push_back(a0); coef_a1_aprox_linear.push_back(a1);}
	void SetCoefCustoOperacao(double a0, double a1, double a2) {coef_custo_oper[0] = a0; coef_custo_oper[1] = a1; coef_custo_oper[2] = a2;}
	double GetCoefA0(int pos) {return coef_a0_aprox_linear[pos];}
	double GetCoefA1(int pos) {return coef_a1_aprox_linear[pos];}
	void AdicionarPontoPt(double ptc_valor) {ptc.push_back(ptc_valor);}
	double GetPtc(int i) {return ptc[i];}
};
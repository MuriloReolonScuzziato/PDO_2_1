#pragma once
#include "Unidades.h"

typedef vector<CUnidades> vetorgrupo;	// o vetor grupo tem dados de cada unidade do grupo

class CHidreletrica
{
	string nome_usina;
	int ident_usina;
	double vmin, vmax, smin, smax, v0, vmeta, custo_agua, dmin, dmax;
	int influe_vert, n_grupos, usina_jusante, tempo_viagem;
	vetorfloat coef_mont, coef_jus, probabilidade;
	vetorfloat2 cen_afluencia;
	vetorfloat Vmin, Vmax;
	vetorint usinas_montante;
public:
	//static int n_hidreletricas;
	vetorgrupo grupoVtr;
	CHidreletrica(string nome_usina_a, int ident, double vmin_a, double vmax_a, double smin_a, double smax_a, int influe_vert_a, int n_grupos_a, int usina_jusante_a, int tempo_viagem_a, vetorfloat coef_mont_a, vetorfloat coef_jus_a);
	~CHidreletrica(void);
	//void add_grupo(int, int, int, double, double, double, double, double, vetorfloat);
	double NivelMontante(double);
	double NivelJusante(double);
	double QuedaBruta(double, double);
	double GetGeracaoOtima(double v_a, double Q_a, double de_a, int grupo);

	string GetNomeUsina() {return nome_usina;}
	int GetIdentUsina() {return ident_usina;}
	double GetVmin() {return vmin;}
	double GetVmax() {return vmax;}
	double GetVmin(int pos) {return Vmin[pos];}
	double GetVmax(int pos) {return Vmax[pos];}
	double GetSmin() {return smin;}
	double GetSmax() {return smax;}
	double GetDmax();
	int GetTempoViagem() {return tempo_viagem;}
	double GetV0() {return v0;}
	double GetVMeta() {return vmeta;}
	double GetAfluencia(int cenario, int periodo) {return cen_afluencia[cenario][periodo];}
	double GetProbAfluencia(int cenario) {return probabilidade[cenario];}
	int GetInflueVert() {return influe_vert;}
	int GetNGrupos() {return n_grupos;}
	int GetUsinaJusante() {return usina_jusante;}
	void SetAfluencia(vetorfloat afluencia_a) {cen_afluencia.push_back(afluencia_a);}
	void SetProbAfluencia(double prob) {probabilidade.push_back(prob);}
	void SetV0(double v0_a) {v0 = v0_a;}
	void SetVMeta(double vmeta_a) {vmeta = vmeta_a;}
	void SetCustoAgua(double custo) {custo_agua = custo;}
	double GetCustoAgua() {return custo_agua;}
	bool FcmIsCte();
	bool FcjIsCte();
	void SetVminSize(int tam) {Vmin.resize(tam);}
	void SetVmaxSize(int tam) {Vmax.resize(tam);}
	void SetVolElem(int pos, double valormin, double valormax) {Vmin[pos] = valormin; Vmax[pos] = valormax;}
	void SetUsinaMont(int ident_usina) {usinas_montante.push_back(ident_usina);}
	int GetNUsinaM() {return int(usinas_montante.size());}
	int GetUsinaMont(int pos) {return usinas_montante[pos];}
};
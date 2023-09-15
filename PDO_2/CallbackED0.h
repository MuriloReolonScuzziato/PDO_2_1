#pragma once

#include "gurobi_c++.h"
#include "Sistema.h"
#include <cmath>
using namespace std;

typedef vector<vetorint> vetorint2;

class CallbackED : public GRBCallback
{
public:
	GRBVar* vars;
    int numvars;
    double lastnode;
	double lastiter;
	int flag7;
	CSistema * sistema_a;
	int T, T1, T2, I, R, n_var_por_periodo;
	double tol;
	//vetorint2 count_persp_cuts;		// numero de usinas por numero de periodos (I x T), conta o numero de cortes por periodo e térmica
	int count_persp_cuts;

	double total_cuts;
	vetorfloat max_difer;
	int delta;
	double difer, Foriginal, p_bar, coeficiente;
	GRBLinExpr restricao;
	GRBVar variavel;

	CallbackED(GRBVar* xvars, int xnumvars, CSistema * sistema_end);
	~CallbackED(void);
	
protected:
    void callback();
};
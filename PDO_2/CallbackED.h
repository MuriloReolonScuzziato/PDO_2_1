#pragma once

#include "gurobi_c++.h"
#include "Sistema.h"
#include <cmath>
using namespace std;

typedef vector<vetorint> vetorint2;
typedef vector<vetorint2> vetorint3;

class CallbackED : public GRBCallback
{
public:
	GRBVar* vars;
    int numvars;
    double lastnode;
	double lastiter;
	double indice;
	int flag7, flag4, flag3;
	CSistema * sistema_a;
	int T, T1, T2, I, L, R, B, n_var_por_periodo;
	double tol;
	int count_persp_cuts;
	MatrixXd Agh, Agt, Agtu, Adef;

	double total_cuts;
	int delta;
	double difer, Foriginal, p_bar, coeficiente;
	GRBLinExpr restricao;
	GRBVar variavel;

	CallbackED(GRBVar* xvars, int xnumvars, CSistema * sistema_end);
	~CallbackED(void);
	
protected:
    void callback();
};
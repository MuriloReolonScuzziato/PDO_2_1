#include "CallbackED.h"

CallbackED::CallbackED(GRBVar* xvars, int xnumvars, CSistema * sistema_end)
{
	vars = xvars;
    numvars = xnumvars;
    lastnode = lastiter = -1;
	sistema_a = sistema_end;
	T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	T1 = sistema_a->GetTt1();
	T2 = sistema_a->GetTt2();
	R = sistema_a->hidreletricasVtr.size();
	I = sistema_a->termeletricasVtr.size();
	n_var_por_periodo = (numvars - R*sistema_a->GetNCenarios()) / T;
	tol = 0.10;	//0.25;
	flag7 = sistema_a->GetFlagTbinaryModel();
	//count_persp_cuts.resize(sistema_a->termeletricasVtr.size());
	//vetorint cuts (T, sistema_a->GetFlagAproxCustoT());
	//vetorint cuts (T, 0);
	//for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
	//	count_persp_cuts[i] = cuts;
	count_persp_cuts = 0;

	//// Escrever numero de cortes
	//total_cuts = 0;
	//for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
	//	for (int t = 0; t < T; t++)
	//		total_cuts += count_persp_cuts[i][t];
	//cout << "Numero total de cortes: (Inicial) " << total_cuts << endl;
	////
	cout << "Numero total de cortes: (Inicial) " << T*sistema_a->termeletricasVtr.size()*(sistema_a->GetFlagMaxAproxCT()-sistema_a->GetFlagAproxCustoT()) << endl;

	max_difer.resize(I);
	for (size_t i = 0; i < I; i++)
		max_difer[i] = sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmax(), 1) - sistema_a->termeletricasVtr[i].CustoOperacao(sistema_a->termeletricasVtr[i].GetPmin(), 1);
}

CallbackED::~CallbackED(void)
{
	delete vars;
	delete sistema_a;
}
void CallbackED::callback()
{
	try
	{
		if (where == GRB_CB_MIPSOL)
		{
			delta = 0;
			difer = 0;
			double * x = getSolution(vars, numvars);
			//cout << "GRB_CP_MIPSOL" << endl;
			for (int t = 0; t < T; t++)
			{
				for (size_t i = 0; i < I; i++)
				{
					if ( x[i + I + delta] > 0.5)
					{
						p_bar = (x[i + delta]+flag7*sistema_a->termeletricasVtr[i].GetPmin()) / x[i + I + delta];
						//Foriginal = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*x[i + I + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*x[i + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(x[i + delta],2);
						Foriginal = sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta]+flag7*sistema_a->termeletricasVtr[i].GetPmin(), x[i + I + delta]);
						difer = (Foriginal - x[i + (3+flag7)*I + delta]) / max_difer[i];
						// conferir se a difernça for maior que a tol.
						//if ( (abs(difer) >= tol) &&  count_persp_cuts[i][t] <= (sistema_a->GetFlagMaxAproxCT()-sistema_a->GetFlagAproxCustoT()) )
						if ( (abs(difer) >= tol) &&  count_persp_cuts < T*sistema_a->termeletricasVtr.size()*(sistema_a->GetFlagMaxAproxCT()-sistema_a->GetFlagAproxCustoT()) )
						{
							//cout << vars[i + delta].get(GRB_StringAttr_VarName) << " " << x[i + delta] << endl;
							//cout << vars[i + I + delta].get(GRB_StringAttr_VarName) << " " << x[i + I + delta] << endl;
							//cout << vars[i + 3*I + delta].get(GRB_StringAttr_VarName) << " " << x[i + 3*I + delta] << endl;
							//cout << "Diferença: " << difer << endl;
							//cout << "addLazy (MIPSOL) in T" << i + 1 << " period " << t << endl;
							// F
							coeficiente = 1;
							variavel = vars[i + (3+flag7)*I + delta];
							restricao.addTerms( &coeficiente, &variavel, 1);
							// pt
							coeficiente = - (2*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*p_bar + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1));
							variavel = vars[i + delta];
							restricao.addTerms( &coeficiente, &variavel, 1);
							// u
							if ( flag7 == 0 )
								coeficiente = - (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) - sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(p_bar,2));
							else
								coeficiente = - (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) - sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(p_bar,2)) - (2*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*p_bar + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)) * sistema_a->termeletricasVtr[i].GetPmin();
							variavel = vars[i + I + delta];
							restricao.addTerms( &coeficiente, &variavel, 1);
							// adicionar restrição
							addLazy(restricao, GRB_GREATER_EQUAL, 0);
							restricao.clear();
							//count_persp_cuts[i][t] += 1;
							count_persp_cuts += 1;
						}
					}
				}
				if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
					delta += n_var_por_periodo + R;
				else
					delta += n_var_por_periodo;
			}
			delete[] x;
			//// Escrever numero de cortes
			//total_cuts = 0;
			//for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			//	for (int t = 0; t < T; t++)
			//		total_cuts += count_persp_cuts[i][t];
			//cout << "Numero total de cortes: (GRB_CB_MIPSOL) " << total_cuts << endl;
			////
			cout << "Numero total de cortes: (GRB_CP_MIPSOL) " << count_persp_cuts << endl;

		}
		else if (where == GRB_CB_MIPNODE)
		{
			delta = 0;
			difer = 0;
			double * x;
			if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
				x = getNodeRel(vars, numvars);
			else
			{
				x = NULL;
				delete[] x;
				return;
			}
			//cout << "GRB_CB_MIPNODE" << endl;
			for (int t = 0; t < T; t++)
			{
				for (size_t i = 0; i < I; i++)
				{
					if ( x[i + I + delta] > 0.5)
					{
						p_bar = (x[i + delta]+flag7*sistema_a->termeletricasVtr[i].GetPmin()) / x[i + I + delta];
						//Foriginal = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*x[i + I + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*x[i + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(x[i + delta],2);
						Foriginal = sistema_a->termeletricasVtr[i].CustoOperacao(x[i + delta]+flag7*sistema_a->termeletricasVtr[i].GetPmin(), x[i + I + delta]);
						difer = (Foriginal - x[i + (3+flag7)*I + delta]) / max_difer[i];
						if ( (abs(difer) >= tol) &&  count_persp_cuts < T*sistema_a->termeletricasVtr.size()*(sistema_a->GetFlagMaxAproxCT()-sistema_a->GetFlagAproxCustoT()) )
						{
							//cout << vars[i + delta].get(GRB_StringAttr_VarName) << " " << x[i + delta] << endl;
							//cout << vars[i + I + delta].get(GRB_StringAttr_VarName) << " " << x[i + I + delta] << endl;
							//cout << vars[i + 3*I + delta].get(GRB_StringAttr_VarName) << " " << x[i + 3*I + delta] << endl;
							//cout << "Diferença: " << difer << endl;
							//cout << "addLazy (MIPNODE) in T" << i + 1 << " period " << t << endl;
							// F
							coeficiente = 1;
							variavel = vars[i + (3+flag7)*I + delta];
							restricao.addTerms( &coeficiente, &variavel, 1);
							// pt
							coeficiente = - (2*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*p_bar + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1));
							variavel = vars[i + delta];
							restricao.addTerms( &coeficiente, &variavel, 1);
							// u
							if ( flag7 == 0 )
								coeficiente = - (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) - sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(p_bar,2));
							else
								coeficiente = - (sistema_a->termeletricasVtr[i].GetCoefCustoOper(0) - sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(p_bar,2)) - (2*sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*p_bar + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)) * sistema_a->termeletricasVtr[i].GetPmin();
							variavel = vars[i + I + delta];
							restricao.addTerms( &coeficiente, &variavel, 1);
							// adicionar restrição
							addLazy(restricao, GRB_GREATER_EQUAL, 0);
							//addCut(restricao, GRB_GREATER_EQUAL, 0);
							restricao.clear();
							count_persp_cuts += 1;

							// ver tamanho do problema estatico com 10 cuts e ir comentando callback para ver onde está o problema (linha)
						}
					}
				}
				if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
					delta += n_var_por_periodo + R;
				else
					delta += n_var_por_periodo;
			}
			delete[] x;
			//// Escrever numero de cortes
			//total_cuts = 0;
			//for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
			//	for (int t = 0; t < T; t++)
			//		total_cuts += count_persp_cuts[i][t];
			//cout << "Numero total de cortes: (GRB_CB_MIPNODE) " << total_cuts << endl;
			////
			cout << "Numero total de cortes: (GRB_CB_MIPNODE) " << count_persp_cuts << endl;

			// conferir o valor de F com a função original quadrática
			// se a diferença for maior que a tol. (e o numero de cortes <= numero max. de cortes por periodo e usina) adicionar um corte
			// caso contrário não adicionar
		}
    }
	catch (GRBException e)
	{
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
    }
	catch (...)
	{
		cout << "Error during callback" << endl;
    }
}
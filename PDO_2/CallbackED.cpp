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
	L = sistema_a->linhasVtr.size();
	B = sistema_a->barrasVtr.size();
	n_var_por_periodo = (numvars - R*sistema_a->GetNCenarios()) / T;
	tol = 1e-6;
	flag7 = sistema_a->GetFlagTbinaryModel();
	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag3 = int (sistema_a->GetFlagPhmax());
	count_persp_cuts = T*sistema_a->termeletricasVtr.size()*sistema_a->GetFlagInitAproxCT();
	//cout << "Numero de cortes: " << count_persp_cuts << endl;
	//cout << "Numero max. de cortes: " << T*sistema_a->termeletricasVtr.size()*(sistema_a->GetFlagMaxAproxCT()) << endl;

	if (sistema_a->GetFlagModeloRede() == 3)
	{
		Agh = MatrixXd::Zero(int (sistema_a->barrasVtr.size()), int(sistema_a->hidreletricasVtr.size()));
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t r = 0; r < sistema_a->hidreletricasVtr.size(); r++)
				for ( size_t br = 0; br < sistema_a->barrasVtr[b].hidrosPtr.size(); br++)
					if (sistema_a->barrasVtr[b].hidrosPtr[br] == &sistema_a->hidreletricasVtr[r])
						Agh(b, r) = 1;
		Agt = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
		for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
			for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
				for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
					if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
						Agt(b, i) = 1;
		Agtu = MatrixXd::Zero(int(sistema_a->barrasVtr.size()), int(sistema_a->termeletricasVtr.size()));
		if (sistema_a->GetFlagTbinaryModel() == 1)
		{
			for ( size_t b = 0; b < sistema_a->barrasVtr.size(); b++)
				for ( size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
					for ( size_t bi = 0; bi < sistema_a->barrasVtr[b].termosPtr.size(); bi++)
						if (sistema_a->barrasVtr[b].termosPtr[bi] == &sistema_a->termeletricasVtr[i])
							Agtu(b, i) = sistema_a->termeletricasVtr[i].GetPmin();
		}
		Adef = MatrixXd::Identity(sistema_a->barrasVtr.size(),sistema_a->barrasVtr.size());
		
		Agh = - sistema_a->Beta * Agh;
		Agt = - sistema_a->Beta * Agt;
		Agtu = - sistema_a->Beta * Agtu;
		Adef = - sistema_a->Beta * Adef;
	}
	else
	{
		Agh.resize(0,0);
		Agt.resize(0,0);
		Agtu.resize(0,0);
		Adef.resize(0,0);
	}
}

CallbackED::~CallbackED(void)
{
	delete vars;
}
void CallbackED::callback()
{
	try
	{
		if (sistema_a->GetFlagMaxAproxCT() >= 2)		// callback para função quadrática
		{
			if (where == GRB_CB_MIPSOL)
			{
				delta = 0;
				difer = 0;
				double * x = getSolution(vars, numvars);
				//cout << "GRB_CP_MIPSOL" << endl;

				//if (count_persp_cuts == T*sistema_a->termeletricasVtr.size()*(sistema_a->GetFlagMaxAproxCT()))
				//{
				//	count_persp_cuts += 1;		// adiciona 1 só para n cair nessa condição mais
				//	cout << "GRB_CP_MIPSOL: Numero max. de cortes atingido!" << endl;
				//}

				for (int t = 0; t < T; t++)
				{
					for (size_t i = 0; i < I; i++)
					{
						if ( x[i + I + delta] > 1e-6)
						{
							p_bar = (x[i + delta]+flag7*sistema_a->termeletricasVtr[i].GetPmin()) / x[i + I + delta];
							if ( flag7 == 0 )
								Foriginal = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*x[i + I + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*x[i + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(x[i + delta],2)/x[i + I + delta];
							else
								Foriginal = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*x[i + I + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*(x[i + delta]+x[i + I + delta]*sistema_a->termeletricasVtr[i].GetPmin()) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(x[i + delta]+x[i + I + delta]*sistema_a->termeletricasVtr[i].GetPmin(),2)/x[i + I + delta];
							if (Foriginal > 0)
								difer = (Foriginal - x[i + (3+flag7)*I + delta]) / Foriginal;
							else
								difer = (Foriginal - x[i + (3+flag7)*I + delta]) / (Foriginal + 1e-9);
							// conferir se a difernça for maior que a tol.
							if ( (difer >= tol) &&  (count_persp_cuts < T*sistema_a->termeletricasVtr.size()*sistema_a->GetFlagMaxAproxCT()) )
							{
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
				// Escrever numero de cortes
				//cout << "Numero total de cortes: (GRB_CP_MIPSOL) " << count_persp_cuts << endl;

			}
			else if (where == GRB_CB_MIPNODE)
			{
				delta = 0;
				difer = 0;
				double * x;
				if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
				{
					x = getNodeRel(vars, numvars);
			
					//cout << "GRB_CB_MIPNODE" << endl;

					//if (count_persp_cuts == T*sistema_a->termeletricasVtr.size()*(sistema_a->GetFlagMaxAproxCT()))
					//{
					//	count_persp_cuts += 1;		// adiciona 1 só para n cair nessa condição mais
					//	cout << "GRB_CB_MIPNODE: Numero max. de cortes atingido!" << endl;
					//}

					for (int t = 0; t < T; t++)
					{
						for (size_t i = 0; i < I; i++)
						{
							if ( x[i + I + delta] > 1e-6)
							{
								p_bar = (x[i + delta]+flag7*sistema_a->termeletricasVtr[i].GetPmin()) / x[i + I + delta];
								if ( flag7 == 0 )
									Foriginal = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*x[i + I + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*x[i + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(x[i + delta],2)/x[i + I + delta];
								else
									Foriginal = sistema_a->termeletricasVtr[i].GetCoefCustoOper(0)*x[i + I + delta] + sistema_a->termeletricasVtr[i].GetCoefCustoOper(1)*(x[i + delta]+x[i + I + delta]*sistema_a->termeletricasVtr[i].GetPmin()) + sistema_a->termeletricasVtr[i].GetCoefCustoOper(2)*pow(x[i + delta]+x[i + I + delta]*sistema_a->termeletricasVtr[i].GetPmin(),2)/x[i + I + delta];
								if (Foriginal > 0)
									difer = (Foriginal - x[i + (3+flag7)*I + delta]) / Foriginal;
								else
									difer = (Foriginal - x[i + (3+flag7)*I + delta]) / (Foriginal + 1e-9);
								if ( (difer >= tol) &&  (count_persp_cuts < T*sistema_a->termeletricasVtr.size()*sistema_a->GetFlagMaxAproxCT()) )
								{
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
								}
							}
						}
						if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
							delta += n_var_por_periodo + R;
						else
							delta += n_var_por_periodo;
					}
					delete[] x;
					// Escrever numero de cortes
					//cout << "Numero total de cortes: (GRB_CB_MIPNODE) " << count_persp_cuts << endl;

					// conferir o valor de F com a função original quadrática
					// se a diferença for maior que a tol. (e o numero de cortes <= numero max. de cortes por periodo e usina) adicionar um corte
					// caso contrário não adicionar
				}
				else
				{
					x = NULL;
					delete[] x;
				}
			}
		}
		if (sistema_a->GetFlagModeloRede() == 3)		// callback para os limites de transmissão do modelo compacto
		{
			if (where == GRB_CB_MIPSOL)
			{
				delta = 0;
				difer = 0;
				double * x = getSolution(vars, numvars);
				//
				int cen = 0;
				double rhs = 0;
				double D = 0;
				int JJ = 0;
				for (size_t i = 0; i < R; i++)
					JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
				VectorXd Ad(B);
				VectorXd gg(L);
				VectorXd g, Rhs;
				VectorXd flow(L);
				for (int t = 0; t < T; t++)
				{
					Ad.resize(sistema_a->barrasVtr.size());
					if (sistema_a->GetTt1() <= t)
						cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					// Determinar o vetor de demandas
					for (size_t b = 0; b < B; b++)
					{
						D = 0;
						for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
							D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
						Ad(b) = D;
					}
					Rhs = - sistema_a->Beta * ( Ad );
					// Determinar o vetor de gerações
					gg = VectorXd::Zero(L);
					g.resize(R);
					for (size_t r = 0; r < R; r++)
						g(r) = x[r + delta + (3+flag4+flag7)*I];
					gg = gg + Agh * g;
					g.resize(I);
					for (size_t i = 0; i < I; i++)
						g(i) = x[i + delta];
					gg = gg + Agt * g;
					g.resize(I);
					if (sistema_a->GetFlagTbinaryModel() == 1)
						for (size_t i = 0; i < I; i++)
							g(i) = x[i + delta + I];
					gg = gg + Agtu * g;
					g.resize(B);
					for (size_t b = 0; b < B; b++)
						g(b) = x[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
					gg = gg + Adef * g;
					// Conferir limites das linhas extrapolados para cada periodo
					flow = gg - Rhs;		// Beta ( d - g): fluxo nas linhas para o periodo t
					for (size_t l = 0; l < L; l++)
					{
						// Adicionar restrição se limites forem violados (individualmente para limite negativo e positivo)
						if ( flow(l) > sistema_a->linhasVtr[l].GetCapacidade() || flow(l) < - sistema_a->linhasVtr[l].GetCapacidade())
						{
							// ph
							for (size_t r = 0; r < R; r++)
							{
								coeficiente = Agh(l, r);
								variavel = vars[r + delta + (3+flag4+flag7)*I];
								restricao.addTerms( &coeficiente, &variavel, 1);
							}
							// pt
							for (size_t i = 0; i < I; i++)
							{
								coeficiente = Agt(l, i);
								variavel = vars[i + delta];
								restricao.addTerms( &coeficiente, &variavel, 1);
							}
							// u
							if (sistema_a->GetFlagTbinaryModel() == 1)
							{
								for (size_t i = 0; i < I; i++)
								{
									coeficiente = Agtu(l, i);
									variavel = vars[i + delta + I];
									restricao.addTerms( &coeficiente, &variavel, 1);
								}
							}
							// def
							for (size_t b = 0; b < B; b++)
							{
								coeficiente = Adef(l, b);
								variavel = vars[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
								restricao.addTerms( &coeficiente, &variavel, 1);
							}
							// adicionar restrição
							//cout << "Teste para conferir rhs" << endl;
							//cout << "periodo: " << t << " linha: " << l << " -> R" << t*L + l + 10 << " e R" << T*L + t*L + l + 10 << endl;
							//cout << "rhs: " << sistema_a->linhasVtr[l].GetCapacidade() + Rhs(l) << endl;
							//cout << "lhs: " << - sistema_a->linhasVtr[l].GetCapacidade() + Rhs(l) << endl;
							if ( flow(l) > sistema_a->linhasVtr[l].GetCapacidade() )
							{
								rhs = sistema_a->linhasVtr[l].GetCapacidade() + Rhs(l);
								addLazy(restricao, GRB_LESS_EQUAL, rhs);
							}
							else
							{
								rhs = - sistema_a->linhasVtr[l].GetCapacidade() + Rhs(l);
								addLazy(restricao, GRB_GREATER_EQUAL, rhs);
							}
							restricao.clear();
						}
					}
					if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
						delta += n_var_por_periodo + R;
					else
						delta += n_var_por_periodo;
				}
				delete[] x;
			}
			//else if (where == GRB_CB_MIPNODE)
			//{
			//	delta = 0;
			//	difer = 0;
			//	double * x;
			//	if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
			//	{
			//		x = getNodeRel(vars, numvars);
			//	
			//		int cen = 0;
			//		double rhs = 0;
			//		double D = 0;
			//		int JJ = 0;
			//		for (size_t i = 0; i < R; i++)
			//			JJ += sistema_a->hidreletricasVtr[i].GetNGrupos();
			//		VectorXd Ad(B);
			//		VectorXd gg(L);
			//		VectorXd g, Rhs;
			//		VectorXd flow(L);
			//		for (int t = 0; t < T; t++)
			//		{
			//			Ad.resize(sistema_a->barrasVtr.size());
			//			if (sistema_a->GetTt1() <= t)
			//				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			//			// Determinar o vetor de demandas
			//			for (size_t b = 0; b < B; b++)
			//			{
			//				D = 0;
			//				for (size_t d = 0; d < sistema_a->barrasVtr[b].demandasPtr.size(); d++)
			//					D += sistema_a->barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(sistema_a->GetTt2() - sistema_a->GetTt1()));
			//				Ad(b) = D;
			//			}
			//			Rhs = - sistema_a->Beta * ( Ad );
			//			// Determinar o vetor de gerações
			//			gg = VectorXd::Zero(L);
			//			g.resize(R);
			//			for (size_t r = 0; r < R; r++)
			//				g(r) = x[r + delta + (3+flag4+flag7)*I];
			//			gg = gg + Agh * g;
			//			g.resize(I);
			//			for (size_t i = 0; i < I; i++)
			//				g(i) = x[i + delta];
			//			gg = gg + Agt * g;
			//			g.resize(I);
			//			if (sistema_a->GetFlagTbinaryModel() == 1)
			//				for (size_t i = 0; i < I; i++)
			//					g(i) = x[i + delta + I];
			//			gg = gg + Agtu * g;
			//			g.resize(B);
			//			for (size_t b = 0; b < B; b++)
			//				g(b) = x[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
			//			gg = gg + Adef * g;
			//			// Conferir limites das linhas extrapolados para cada periodo
			//			flow = gg - Rhs;		// Beta ( d - g): fluxo nas linhas para o periodo t
			//			for (size_t l = 0; l < L; l++)
			//			{
			//				// Adicionar restrição se limites forem violados (individualmente para limite negativo e positivo)
			//				if ( flow(l) > sistema_a->linhasVtr[l].GetCapacidade() || flow(l) < - sistema_a->linhasVtr[l].GetCapacidade())
			//				{
			//					// ph
			//					for (size_t r = 0; r < R; r++)
			//					{
			//						coeficiente = Agh(l, r);
			//						variavel = vars[r + delta + (3+flag4+flag7)*I];
			//						restricao.addTerms( &coeficiente, &variavel, 1);
			//					}
			//					// pt
			//					for (size_t i = 0; i < I; i++)
			//					{
			//						coeficiente = Agt(l, i);
			//						variavel = vars[i + delta];
			//						restricao.addTerms( &coeficiente, &variavel, 1);
			//					}
			//					// u
			//					if (sistema_a->GetFlagTbinaryModel() == 1)
			//					{
			//						for (size_t i = 0; i < I; i++)
			//						{
			//							coeficiente = Agtu(l, i);
			//							variavel = vars[i + delta + I];
			//							restricao.addTerms( &coeficiente, &variavel, 1);
			//						}
			//					}
			//					// def
			//					for (size_t b = 0; b < B; b++)
			//					{
			//						coeficiente = Adef(l, b);
			//						variavel = vars[b + delta + (3+flag4+flag7)*(sistema_a->termeletricasVtr.size()) + (4+flag3)*(sistema_a->hidreletricasVtr.size()) + 3*JJ];
			//						restricao.addTerms( &coeficiente, &variavel, 1);
			//					}
			//					// adicionar restrição
			//					if ( flow(l) > sistema_a->linhasVtr[l].GetCapacidade() )
			//					{
			//						rhs = sistema_a->linhasVtr[l].GetCapacidade() + Rhs(l);
			//						addLazy(restricao, GRB_LESS_EQUAL, rhs);
			//					}
			//					else
			//					{
			//						rhs = - sistema_a->linhasVtr[l].GetCapacidade() + Rhs(l);
			//						addLazy(restricao, GRB_GREATER_EQUAL, rhs);
			//					}
			//					restricao.clear();
			//				}
			//			}
			//			if (((t + 1 - T1)%(T2 - T1) == 0) && ((t + 1 - T1) > 0))
			//				delta += n_var_por_periodo + R;
			//			else
			//				delta += n_var_por_periodo;
			//		}
			//		delete[] x;
			//	}
			//	else
			//	{
			//		x = NULL;
			//		delete[] x;
			//	}
			//}
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
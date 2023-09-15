#include "Spcdec3SPT.h"

// Criar restrições
// ------------------------------------------------

//Funções para montar o lado esquerdo das restrições lineares
CMatrizEsparsa Spcdec3SPT::MatrizTup(int n_a, int nu)
{
	// Tup é o tempo que a unidade deve ficar ligada após ser ligada, portanto para ficar ligada por 3 periodos, Tup = 2;
	if (flag7 == 0)
	{
		// T * Tup
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tup1, Tup2, Tup, cen, t_a;
		int c = 1;
		int ca = 0;
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT1())
			Tup1 = 0;
		else
			Tup1 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT2())
			Tup2 = 0;
		else
			Tup2 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT2());

		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}
			if (t < sistema_a->GetTt1())
				Tup = Tup1;
			else
				Tup = Tup2;

			if (Tup > 0)
			{
				for (int tt = 0; tt < Tup; tt++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c, 1);
					if (t_a - tt - 1 >= 0)		// cada if adiciona um elemento dos dois elementos do lado direito da equação
						if (t_a - tt - 1 >= sistema_a->GetTt1())
							a.InserirElemento(0, c - tt*(n_a / T) - 1*(n_a / T), -1);
						else
							a.InserirElemento(0, ca - tt*(n_a / T) - 1*(n_a / T), -1);
					if (t_a - tt - 2 >= 0)
						if (t_a - tt - 2 >= sistema_a->GetTt1())
							a.InserirElemento(0, c - tt*(n_a / T) - 2*(n_a / T), 1);
						else
							a.InserirElemento(0, ca - tt*(n_a / T) - 2*(n_a / T), 1);
					
					Alin.JuntarColuna(&a);
				}
			}
			else {}
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T (para usinas com Tup > 0)
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tup1, Tup2, Tup, cen, t_a;
		int c = 1;
		int ca = 0;
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT1())
			Tup1 = 0;
		else
			Tup1 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT2())
			Tup2 = 0;
		else
			Tup2 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT2());

		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}
			if (t < sistema_a->GetTt1())
				Tup = Tup1;
			else
				Tup = Tup2;

			if (Tup > 0)
			{
				a.RemoverTodosElementos();
				if  (t_a >= Tup)
				{
					a.InserirElemento(0, c, -1);		// adiciona termo de u
					for (int tt = 0; tt <= Tup; tt++)
					{
						if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
							a.InserirElemento(0, c + 1 - tt*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + 1 - tt*(n_a / T), 1);
					}
				}
				else if ((Tup >= sistema_a->GetTt2()) && (t_a == sistema_a->GetTt2() - 1))
				{
					a.InserirElemento(0, c, -1);		// adiciona termo de u
					for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
					{
						if (t_a - tt >= sistema_a->GetTt1())		// adiciona cada termo do somatório (up_t)
							a.InserirElemento(0, c + 1 - tt*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + 1 - tt*(n_a / T), 1);
					}
				}
				Alin.JuntarColuna(&a);
			}
			else {}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Spcdec3SPT::MatrizTdown(int n_a, int nu)
{
	// O valor de Tdown tem significado similar ao Tup
	if (flag7 == 0)
	{
		// T * Tdown
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tdown1, Tdown2, Tdown, cen, t_a;
		int c = 1;
		int ca = 0;
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT1())
			Tdown1 = 0;
		else
			Tdown1 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT2())
			Tdown2 = 0;
		else
			Tdown2 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT2());

		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}
		
			if (t < sistema_a->GetTt1())
				Tdown = Tdown1;
			else
				Tdown = Tdown2;

			if (Tdown > 0)
			{
				for (int tt = 0; tt < Tdown; tt++)
				{
					a.RemoverTodosElementos();
					a.InserirElemento(0, c , 1);
					if (t_a - tt - 1 >= 0)		// cada if adiciona um elemento dos dois elementos do lado direito da equação
						if (t_a - tt - 1 >= sistema_a->GetTt1())
							a.InserirElemento(0, c - tt*(n_a / T) - 1*(n_a / T), -1);
						else
							a.InserirElemento(0, ca - tt*(n_a / T) - 1*(n_a / T), -1);
					if (t_a - tt - 2 >= 0)
						if (t_a - tt - 2 >= sistema_a->GetTt1())
							a.InserirElemento(0, c - tt*(n_a / T) - 2*(n_a / T), 1);
						else
							a.InserirElemento(0, ca - tt*(n_a / T) - 2*(n_a / T), 1);
					
					Alin.JuntarColuna(&a);
				}
			}
			else {}
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T (para usinas com Tdown > 0)
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int Tdown1, Tdown2, Tdown, cen, t_a;
		int c = 1;
		int ca = 0;
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT1())
			Tdown1 = 0;
		else
			Tdown1 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT2())
			Tdown2 = 0;
		else
			Tdown2 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT2());

		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}
		
			if (t < sistema_a->GetTt1())
				Tdown = Tdown1;
			else
				Tdown = Tdown2;

			if (Tdown > 0)
			{
				a.RemoverTodosElementos();
				if  (t_a >= Tdown)
				{
					a.InserirElemento(0, c, 1);
					for (int tt = 0; tt <= Tdown; tt++)
					{
						if (t_a - tt >= sistema_a->GetTt1())
							a.InserirElemento(0, c + 2 - tt*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + 2 - tt*(n_a / T), 1);
					}
				}
				else if ((Tdown >= sistema_a->GetTt2()) && (t_a == sistema_a->GetTt2() - 1))
				{
					a.InserirElemento(0, c, 1);
					for (int tt = 0; tt <= sistema_a->GetTt2() - 1; tt++)
					{
						if (t_a - tt >= sistema_a->GetTt1())
							a.InserirElemento(0, c + 2 - tt*(n_a / T), 1);
						else
							a.InserirElemento(0, ca + 2 - tt*(n_a / T), 1);
					}
				}
				Alin.JuntarColuna(&a);
			}
			else {}
			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Spcdec3SPT::MatrizTupDown(int n_a, int nu)
{
	// T (nos primeiros periodos tem-se linhas com 0, pois não são consideradas determinadas restrições para t < 2 ou t < Tup + 2 ou t < Tdown + 2
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int cen, t_a;
	int c = 1;
	int ca = 0;
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}
		
		a.RemoverTodosElementos();
		if (t >= 1)
		{
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + 1, -1);
			a.InserirElemento(0, c + 2, 1);
			if (t_a - 1 >= sistema_a->GetTt1())
				a.InserirElemento(0, c - 1*(n_a / T), -1);
			else
				a.InserirElemento(0, ca - 1*(n_a / T), -1);
		}
		else
		{
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + 1, -1);
			a.InserirElemento(0, c + 2, 1);
		}
		Alin.JuntarColuna(&a);
		
		c += (n_a / T);
	}
	return Alin;
}
CMatrizEsparsa Spcdec3SPT::MatrizRampaUp(int n_a, int nu)
{
	if (flag7 == 0)
	{
		// T
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int c = 0;
		int ca, t_a, cen;
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}

			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
				{
					a.InserirElemento(0, c - (n_a / T), -1);
					a.InserirElemento(0, c + 1 - (n_a / T), sistema_a->termeletricasVtr[nu].GetPmin() - sistema_a->termeletricasVtr[nu].GetRampaUp());
				}
				else
				{
					a.InserirElemento(0, ca - (n_a / T), -1);
					a.InserirElemento(0, ca + 1 - (n_a / T), sistema_a->termeletricasVtr[nu].GetPmin() - sistema_a->termeletricasVtr[nu].GetRampaUp());
				}
			}
			Alin.JuntarColuna(&a);

			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		// T
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int c = 0;
		int ca, t_a, cen;
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}

			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
					a.InserirElemento(0, c - (n_a / T), -1);
				else
					a.InserirElemento(0, ca - (n_a / T), -1);
			}
			Alin.JuntarColuna(&a);

			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Spcdec3SPT::MatrizRampaDown(int n_a, int nu)
{
	// T
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int c = 0;
		int ca, t_a, cen;
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}

			a.RemoverTodosElementos();
			a.InserirElemento(0, c , -1);
			a.InserirElemento(0, c + 1, sistema_a->termeletricasVtr[nu].GetPmin() - sistema_a->termeletricasVtr[nu].GetRampaDown());
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
					a.InserirElemento(0, c - (n_a / T), 1);
				else
					a.InserirElemento(0, ca - (n_a / T), 1);
			}
			Alin.JuntarColuna(&a);

			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int c = 0;
		int ca, t_a, cen;
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}

			a.RemoverTodosElementos();
			a.InserirElemento(0, c , -1);
			if (t_a > 0)
			{
				if (t_a - 1 >= sistema_a->GetTt1())
					a.InserirElemento(0, c - (n_a / T), 1);
				else
					a.InserirElemento(0, ca - (n_a / T), 1);
			}
			Alin.JuntarColuna(&a);

			c += (n_a / T);
		}
		return Alin;
	}
}
CMatrizEsparsa Spcdec3SPT::MatrizLimPtMin(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + 1, - sistema_a->termeletricasVtr[nu].GetPmin());
			Alin.JuntarColuna(&a);
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int ca, t_a, cen;
		int c = 0;
		c += (n_a / T);						// referente à t = 0
		// t = 0
		a.RemoverTodosElementos();
		a.InserirElemento(0, 3, (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));		// ud
		Alin.JuntarColuna(&a);
		for (int t = 1; t < T; t++)			// nó inicial, já adicionadas restrições, por isso t começa em 1
		{
			if (t >= sistema_a->GetTt2())
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			}
			else
			{
				ca = c;
				t_a = t;
			}

			if (sistema_a->termeletricasVtr[nu].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + 3, (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));		// ud
				if (t_a - 1 >= sistema_a->GetTt1())		// todos periodos menos do inicio do segundo estágio (incluindo o de final de horizonte)
				{
					a.InserirElemento(0, c - (n_a / T), 1);		// pt-1
					a.InserirElemento(0, c - (n_a / T) + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));		// u-1
				}
				else									// periodos no inicio do segundo estágio
				{
					a.InserirElemento(0, ca - (n_a / T), 1);
					a.InserirElemento(0, ca - (n_a / T) + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
				}
				Alin.JuntarColuna(&a);
			}
			else
			{
				a.RemoverTodosElementos();
				a.InserirElemento(0, c + 3, (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));		// ud
				if (t_a - 1 >= sistema_a->GetTt1())		// todos periodos menos do inicio do segundo estágio (incluindo o de final de horizonte)
				{
					a.InserirElemento(0, c - (n_a / T), 1);		// pt-1
					a.InserirElemento(0, c - (n_a / T) + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));		// u-1
					a.InserirElemento(0, c - (n_a / T) + 2, (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));		// up-1
				}
				else									// periodos no inicio do segundo estágio
				{
					a.InserirElemento(0, ca - (n_a / T), 1);
					a.InserirElemento(0, ca - (n_a / T) + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
					a.InserirElemento(0, ca - (n_a / T) + 2, (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
				}
				Alin.JuntarColuna(&a);
			}
			c += (n_a / T);
		}
		return Alin;


		// Código antigo
		//int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		//CMatrizEsparsa Alin(0, n_a);
		//CMatrizEsparsa a(1, n_a);
		//int c = 0;
		//for (int t = 0; t < T; t++)
		//{
		//	if (sistema_a->termeletricasVtr[nu].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
		//	{
		//		if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//		{
		//			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c, 1);
		//				a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//				a.InserirElemento(0, c + 3 + (n_a / T)*(1 + n_cen*(sistema_a->GetTt2() - sistema_a->GetTt1())), sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//				Alin.JuntarColuna(&a);
		//			}
		//		}
		//		else if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//		{
		//			a.RemoverTodosElementos();
		//			a.InserirElemento(0, c, 1);
		//			a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//			Alin.JuntarColuna(&a);
		//		}
		//		else
		//		{
		//			a.RemoverTodosElementos();
		//			a.InserirElemento(0, c, 1);
		//			a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//			a.InserirElemento(0, c + 3 + (n_a / T), sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//			Alin.JuntarColuna(&a);
		//		}
		//	}
		//	else
		//	{
		//		if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//		{
		//			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//			{
		//				a.RemoverTodosElementos();
		//				a.InserirElemento(0, c, 1);
		//				a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//				a.InserirElemento(0, c + 2, sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//				a.InserirElemento(0, c + 3 + (n_a / T)*(1 + n_cen*(sistema_a->GetTt2() - sistema_a->GetTt1())), sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//				Alin.JuntarColuna(&a);
		//			}
		//		}
		//		else if (((t + 1 - sistema_a->GetTt1())%(sistema_a->GetTt2() - sistema_a->GetTt1()) == 0) && ((t + 1 - sistema_a->GetTt1()) > 0))
		//		{
		//			a.RemoverTodosElementos();
		//			a.InserirElemento(0, c, 1);
		//			a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//			a.InserirElemento(0, c + 2, sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//			Alin.JuntarColuna(&a);
		//		}
		//		else
		//		{
		//			a.RemoverTodosElementos();
		//			a.InserirElemento(0, c, 1);
		//			a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//			a.InserirElemento(0, c + 2, sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//			a.InserirElemento(0, c + 3 + (n_a / T), sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//			Alin.JuntarColuna(&a);
		//		}
		//	}
		//	c += (n_a / T);
		//}
		//return Alin;
	}
}
CMatrizEsparsa Spcdec3SPT::MatrizLimPtMax(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0,n_a);
		CMatrizEsparsa a(1,n_a);
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + 1, - sistema_a->termeletricasVtr[nu].GetPmax());
			Alin.JuntarColuna(&a);
			c += (n_a / T);
		}
		return Alin;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		CMatrizEsparsa Alin(0, n_a);
		CMatrizEsparsa a(1, n_a);
		int c = 0;
		for (int t = 0; t < T; t++)
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c, 1);
			a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
			a.InserirElemento(0, c + 2, sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
			Alin.JuntarColuna(&a);
			c += (n_a / T);
		}
		return Alin;

		// código antigo
		//int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		//CMatrizEsparsa Alin(0,n_a);
		//CMatrizEsparsa a(1,n_a);
		//int c = 0;
		//for (int t = 0; t < T; t++)
		//{
		//	if (sistema_a->termeletricasVtr[nu].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
		//	{
		//		a.RemoverTodosElementos();
		//		a.InserirElemento(0, c, 1);
		//		a.InserirElemento(0, c + 1, - (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()));
		//		a.InserirElemento(0, c + 2, sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin());
		//		Alin.JuntarColuna(&a);
		//	}
		//	else
		//	{
		//		a.RemoverTodosElementos();		// restrição já adicionada na função ptmin
		//		Alin.JuntarColuna(&a);
		//	}
		//	c += (n_a / T);
		//}
		//return Alin;
	}
}
CMatrizEsparsa Spcdec3SPT::MatrizRestCP(int n_a, int nu)
{
	// T * I
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int c = 1;
	int ca, t_a, cen;
	for (int t = 0; t < T; t++)
	{
		if (t >= sistema_a->GetTt2())
		{
			cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
			ca = c - (n_a / T)*(cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) + 0*sistema_a->GetTt1());
			t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
		}
		else
		{
			ca = c;
			t_a = t;
		}

		a.RemoverTodosElementos();
		a.InserirElemento(0, c + 1, 1);
		a.InserirElemento(0, c, - sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());
		if (t_a > 0)
		{
			if (t_a - 1 >= sistema_a->GetTt1())
				a.InserirElemento(0, c - (n_a / T), sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());
			else
				a.InserirElemento(0, ca - (n_a / T), sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());
		}
		Alin.JuntarColuna(&a);

		c += (n_a / T);
	}
	return Alin;

	//int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	//CMatrizEsparsa Alin(0,n_a);
	//CMatrizEsparsa a(1,n_a);
	//CMatrizEsparsa Alin3(0,n_a);
	//int c = 1;
	//Alin.ZerarMatriz(0, n_a);
	//for (int t = 0; t < T; t++)
	//{
	//	a.RemoverTodosElementos();
	//	if (t == 0)
	//	{
	//		a.InserirElemento(0, c + 1, 1);
	//		a.InserirElemento(0, c, - sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());
	//	}
	//	else
	//	{
	//		a.InserirElemento(0, c + 1, 1);
	//		a.InserirElemento(0, c, - sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());
	//		a.InserirElemento(0, c - (n_a / T), sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());
	//	}
	//	Alin.JuntarColuna(&a);
	//	c += (n_a / T);
	//}
	//CMatrizEsparsa Alin2((sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1())), n_a);
	//CMatrizEsparsa A1(sistema_a->GetTt2() - sistema_a->GetTt1(), (n_a / T)*sistema_a->GetTt1());
	//CMatrizEsparsa A2((sistema_a->GetTt2() - sistema_a->GetTt1()), (n_a / T)*sistema_a->GetTt2() - (n_a / T)*sistema_a->GetTt1());
	//Alin2.InserirMatriz(0, 0, sistema_a->GetTt1() - 1, (n_a / T)*sistema_a->GetTt1() - 1, &Alin, 0, 0);
	//A1.InserirMatriz(0, 0, sistema_a->GetTt2() - sistema_a->GetTt1() - 1, (n_a / T)*sistema_a->GetTt1() - 1, &Alin, sistema_a->GetTt1(), 0);
	//A2.InserirMatriz(0, 0, sistema_a->GetTt2() - sistema_a->GetTt1() - 1, (n_a / T)*sistema_a->GetTt2() - (n_a / T)*sistema_a->GetTt1() - 1, &Alin, sistema_a->GetTt1(), (n_a / T)*sistema_a->GetTt1());
	//int linha = 0;
	//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
	//{
	//	linha = sistema_a->GetTt1() + (cen)*(sistema_a->GetTt2() - sistema_a->GetTt1());
	//	Alin2.InserirMatriz(linha, 0, linha + (sistema_a->GetTt2() - sistema_a->GetTt1()) - 1, (n_a / T)*sistema_a->GetTt1() - 1, &A1, 0, 0);
	//	Alin2.InserirMatriz(linha, (n_a / T)*(sistema_a->GetTt1()+cen*(sistema_a->GetTt2() - sistema_a->GetTt1())), linha + (sistema_a->GetTt2() - sistema_a->GetTt1()) - 1, (n_a / T)*(sistema_a->GetTt2() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1())) - 1, &A2, 0, 0);
	//}
	//Alin3.JuntarColuna(&Alin2);
	//return Alin3;
}
CMatrizEsparsa Spcdec3SPT::MatrizCortesF(int n_a, int nu)
{
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	CMatrizEsparsa Alin(0, n_a);
	CMatrizEsparsa a(1, n_a);
	int c = 0;
	
	for (int t = 0; t < T; t++)
	{
		for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)		// Adiciona uma restrição por corte
		{
			a.RemoverTodosElementos();
			a.InserirElemento(0, c, - sistema_a->termeletricasVtr[nu].GetCoefA1(n_cort));
			if ( flag7 == 0 )
				a.InserirElemento(0, c + 1, - sistema_a->termeletricasVtr[nu].GetCoefA0(n_cort));
			else
				a.InserirElemento(0, c + 1, - sistema_a->termeletricasVtr[nu].GetCoefA0(n_cort) - sistema_a->termeletricasVtr[nu].GetCoefA1(n_cort) * sistema_a->termeletricasVtr[nu].GetPmin());
			a.InserirElemento(0, c + 3 + flag7, 1);
			Alin.JuntarColuna(&a);
		}
		c += (n_a / T);
	}
	return Alin;
}

// Funções para montar limites das restrições lineares (A*x (tipo) Limite) -> Lim = [tipo Limite] -> tipo = (0,1,2):(GRB_LESS_EQUAL, GRB_EQUAL, or GRB_GREATER_EQUAL) 
vetorfloat2 Spcdec3SPT::LimTup(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int Tup1, Tup2, Tup, t_a;
		//
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT1())
			Tup1 = 0;
		else
			Tup1 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT2())
			Tup2 = 0;
		else
			Tup2 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT2());
		for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
		{
			if (t >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = t;

			if (t < sistema_a->GetTt1())
				Tup = Tup1;
			else
				Tup = Tup2;

			if (Tup > 0)		// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
			{
				for (int tt = 0; tt < Tup; tt++)
				{
					if (t_a - tt <= 0)		// if referente ao elemento tt (se ele for constante entra nesse loop)
					{
						if (t_a - tt > - sistema_a->termeletricasVtr[nu].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] = double (sistema_a->termeletricasVtr[nu].GetU0());
							L[0] = 2;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] = double((1 - sistema_a->termeletricasVtr[nu].GetU0()));
							L[0] = 2;
						}
					}
					else
					{
						L[1] = 0;
						L[0] = 2;
					}
					if (t_a - tt - 1 <= 0)	// if referente ao elemento tt - 1 (se ele for constante entra nesse loop)
					{
						if (t_a - tt - 1 > - sistema_a->termeletricasVtr[nu].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] += - sistema_a->termeletricasVtr[nu].GetU0();
							//L[0] = 2;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] += - (1 - sistema_a->termeletricasVtr[nu].GetU0());
							//L[0] = 2;
						}
					}
					else
					{
						//L[1] = 0;
						//L[0] = 2;
					}
					Lim.push_back(L);
					L.clear();L.resize(2);
				}
			}
		}
		return Lim;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int Tup1, Tup2, Tup;
		//
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0; L[1] = 0;
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT1())
			Tup1 = 0;
		else
			Tup1 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTUp() < sistema_a->GetDeltaT2())
			Tup2 = 0;
		else
			Tup2 = int (sistema_a->termeletricasVtr[nu].GetTUp() / sistema_a->GetDeltaT2());
		for (int t = 0; t < T; t++)
		{
			if (t < sistema_a->GetTt1())
				Tup = Tup1;
			else
				Tup = Tup2;

			if (Tup > 0)		// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
				Lim.push_back(L);
		}
		return Lim;
	}
}
vetorfloat2 Spcdec3SPT::LimTdown(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int Tdown1, Tdown2, Tdown, t_a;
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0;
		L[1] = 1;
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT1())
			Tdown1 = 0;
		else
			Tdown1 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT2())
			Tdown2 = 0;
		else
			Tdown2 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT2());

		for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
		{
			if (t >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = t;

			if (t < sistema_a->GetTt1())
				Tdown = Tdown1;
			else
				Tdown = Tdown2;

			if (Tdown > 0)		// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
			{
				for (int tt = 0; tt < Tdown; tt++)
				{
					if (t_a - tt <= 0)		// if referente ao elemento tt (se ele for constante entra nesse loop)
					{
						if (t_a - tt > - sistema_a->termeletricasVtr[nu].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] += double (sistema_a->termeletricasVtr[nu].GetU0());
							//L[0] = 0;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] += double((1 - sistema_a->termeletricasVtr[nu].GetU0()));
							//L[0] = 0;
						}
					}
					else
					{
						//L[1] = 0;
						//L[0] = 0;
					}
					if (t_a - tt - 1 <= 0)	// if referente ao elemento tt - 1 (se ele for constante entra nesse loop)
					{
						if (t_a - tt - 1 > - sistema_a->termeletricasVtr[nu].GetX0())	// estado da usina x0 periodos antes (de t(-x0) a t(0))
						{
							L[1] += - sistema_a->termeletricasVtr[nu].GetU0();
							//L[0] = 0;
						}
						else	// estado da usina antes do perido -x0
						{	
							L[1] += - (1 - sistema_a->termeletricasVtr[nu].GetU0());
							//L[0] = 0;
						}
					}
					else
					{
						//L[1] = 0;
						//L[0] = 0;
					}
					Lim.push_back(L);
					L.clear();L.resize(2);L[0] = 0;L[1] = 1;
				}
			}
		}
		return Lim;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int Tdown1, Tdown2, Tdown;
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0;L[1] = 1;
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT1())
			Tdown1 = 0;
		else
			Tdown1 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT1());
		if (sistema_a->termeletricasVtr[nu].GetTDown() < sistema_a->GetDeltaT2())
			Tdown2 = 0;
		else
			Tdown2 = int (sistema_a->termeletricasVtr[nu].GetTDown() / sistema_a->GetDeltaT2());

		for (int t = 0; t < T; t++)		// Em tese só precisaria rodar até T2 e para todos nos do cen > 0 o RHS seria igual
		{
			if (t < sistema_a->GetTt1())
				Tdown = Tdown1;
			else
				Tdown = Tdown2;

			if (Tdown > 0)		// Se tup = 0, não adiciona-se restrição (a usina fica no min. um período em operação)
				Lim.push_back(L);
		}
		return Lim;
	}
}
vetorfloat2 Spcdec3SPT::LimTupDown(int n_a, int nu)
{
	// T
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	vetorfloat2 Lim;
	vetorfloat L, L0;
	L.resize(2); L0.resize(2);
	L0[0] = 1;
	L[0] = 1;L[1] = 0;
	for (int t = 0; t < T; t++)
	{
		if (t >= 1)
			Lim.push_back(L);
		else
		{
			L0[1] = sistema_a->termeletricasVtr[nu].GetU0();
			Lim.push_back(L0);
		}
	}
	return Lim;
}
vetorfloat2 Spcdec3SPT::LimRampaUp(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetPmin() + sistema_a->termeletricasVtr[nu].GetPt0() + sistema_a->termeletricasVtr[nu].GetU0() * (sistema_a->termeletricasVtr[nu].GetRampaUp() - sistema_a->termeletricasVtr[nu].GetPmin());
			else
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetPmin();
		}
		return Lim;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetRampaUp() + (sistema_a->termeletricasVtr[nu].GetPt0() - sistema_a->termeletricasVtr[nu].GetPmin()*sistema_a->termeletricasVtr[nu].GetU0());
			else
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetRampaUp();
		}
		return Lim;
	}
}
vetorfloat2 Spcdec3SPT::LimRampaDown(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetPmin() - sistema_a->termeletricasVtr[nu].GetPt0();
			else
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetPmin();
		}
		return Lim;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
		{
			if (t == 0)
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetRampaDown() - (sistema_a->termeletricasVtr[nu].GetPt0() - sistema_a->termeletricasVtr[nu].GetPmin()*sistema_a->termeletricasVtr[nu].GetU0());
			else
				Lim[t][1] = sistema_a->termeletricasVtr[nu].GetRampaDown();
		}
		return Lim;
	}
}
vetorfloat2 Spcdec3SPT::LimPtMin(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		for (int t = 0; t < T; t++)
			Lim[t][0] = 2;
		return Lim;
	}
	else
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		vetorfloat L; L.resize(2);
		L[0] = 0; L[1] = 0;
		// t = 0
		if (sistema_a->termeletricasVtr[nu].GetU0() == 1)
			if (sistema_a->termeletricasVtr[nu].GetX0() > 1)  // indica que a usina foi ligada antes do periodo t = -1, portanto up-1 = 0;
				L[1] = sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPt0();		// valor de pt0 está em valor de ptmin a ptmax e n de 0 a ptmax-ptmin
		Lim.push_back(L);
		L[0] = 0; L[1] = 0;
		for (int t = 1; t < T; t++)
			Lim.push_back(L);
		return Lim;

		// código antigo
		//// sum(t*n_aberturas)
		//int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		//vetorfloat2 Lim;
		//vetorfloat L; L.resize(2);
		//L[0] = 0; L[1] = 0;
		//for (int t = 0; t < T; t++)
		//{
		//	if (sistema_a->termeletricasVtr[nu].GetTUp() == 0)		// units that can be online for just one period, se Tup > 0 adiciona-se somente uma restrição para os limites de pt (ptmin restrição e ptmax linha de zeros)
		//	{
		//		if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//		{
		//			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//				Lim.push_back(L);
		//		}
		//		else
		//			Lim.push_back(L);
		//	}
		//	else
		//	{
		//		if (t == sistema_a->GetTt1() - 1)		// nó com múltiplas aberturas (não somente uma)
		//		{
		//			for (int n_cen = 0; n_cen < sistema_a->GetNCenarios(); n_cen++)
		//				Lim.push_back(L);
		//		}
		//		else
		//			Lim.push_back(L);
		//	}
		//}
		//return Lim;
	}
}
vetorfloat2 Spcdec3SPT::LimPtMax(int n_a, int nu)
{
	if (flag7 == 0)
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		return Lim;
	}
	else
	{
		// T
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		vetorfloat2 Lim;
		DimensionarMatriz(&Lim, T, 2);
		IniciaMatriz(&Lim, 0);
		return Lim;
	}
}
vetorfloat2 Spcdec3SPT::LimRestCP(int n_a, int nu)
{
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, T, 2);
	IniciaMatriz(&Lim, 0);
	for (int t = 0; t < T; t++)
	{
		if (t == 0)
		{
			Lim[t][1] = - sistema_a->termeletricasVtr[nu].GetCoefCustoPartida()*sistema_a->termeletricasVtr[nu].GetU0();
			Lim[t][0] = 2;
		}
		else
		{
			Lim[t][1] = 0;
			Lim[t][0] = 2;
		}
	}
	return Lim;
}
vetorfloat2 Spcdec3SPT::LimCortesF(int n_a, int nu)
{
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	vetorfloat2 Lim;
	DimensionarMatriz(&Lim, T * sistema_a->GetFlagInitAproxCT(), 2);
	IniciaMatriz(&Lim, 0);
	int c = 0;
	for (int t = 0; t < T; t++)
	{
		for (int n_cort = 0; n_cort < sistema_a->GetFlagInitAproxCT(); n_cort++)
		{
			Lim[c][0] = 2;
			c++;
		}
	}
	return Lim;
}

// Matriz dos coeficientes e limites das restrições lineares
void Spcdec3SPT::MatrizRestricoesLineares(int n_a, CMatrizEsparsa &MM, int nu)		// Matriz M é composta de 3 vetores, cada coluna (vetor) é : [indexL indexC Valor]
{
	CMatrizEsparsa M(0);
	M = MatrizTup(n_a, nu);									
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizTdown(n_a, nu);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if ( sistema_a->GetFlagTbinaryModel() == 1 )
	{
		M = MatrizTupDown(n_a, nu);
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
	M = MatrizRampaUp(n_a, nu);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizRampaDown(n_a, nu);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPtMin(n_a, nu);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	M = MatrizLimPtMax(n_a, nu);
	MM.JuntarColuna(&M); M.ZerarMatriz();
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		M = MatrizRestCP(n_a, nu);
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		M = MatrizCortesF(n_a, nu);
		MM.JuntarColuna(&M); M.ZerarMatriz();
	}
}
void Spcdec3SPT::MatrizLimitesLineares(int n_a, vetorint * LimTipo, vetorfloat * LimValor, int nu)
{
	vetorfloat2 L;
	L = LimTup(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimTdown(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	if ( sistema_a->GetFlagTbinaryModel() == 1 )
	{
		L = LimTupDown(n_a, nu);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	L = LimRampaUp(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimRampaDown(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPtMin(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	L = LimPtMax(n_a, nu);
	AlocarLimites(&L, LimTipo, LimValor);
	if ( sistema_a->GetFlagTbinaryModel() == 0 )
	{
		L = LimRestCP(n_a, nu);
		AlocarLimites(&L, LimTipo, LimValor);
	}
	if (sistema_a->GetFlagInitAproxCT() > 1)
	{
		L = LimCortesF(n_a, nu);
		AlocarLimites(&L, LimTipo, LimValor);
	}
}
// ------------------------------------------------

Spcdec3SPT::Spcdec3SPT(CSistema * const sistema_end, const GRBEnv &ambiente_gurobi)
{
	sistema_a = sistema_end;

	if (sistema_a->GetFlagInitAproxCT() > 1)
		flag4 = 1;
	else
		flag4 = 0;
	flag7 = int (sistema_a->GetFlagTbinaryModel());
	n_usinas = sistema_a->termeletricasVtr.size();
	modelosGRB.resize(n_usinas);
	vars.resize(n_usinas);
	for (int nu = 0; nu < n_usinas; nu++)
		modelosGRB[nu] = new GRBModel(ambiente_gurobi);		// um modelo de otimização para cada usina
	n.resize(n_usinas);
	for (size_t i = 0; i < sistema_a->termeletricasVtr.size(); i++)
		n[i] = (sistema_a->GetTt1() + (sistema_a->GetTt2() - sistema_a->GetTt1()) * sistema_a->GetNCenarios() ) * (3 + flag4 + flag7);

	// Criar variáveis
	for (int nu = 0; nu < n_usinas; nu++)
	{
		CriarVariaveis(nu);
		modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.
		vars[nu] = modelosGRB[nu]->getVars();
	}
	
	// Adicionar restrições
	for (int nu = 0; nu < n_usinas; nu++)
	{
		CriarRestricoes(nu);
		modelosGRB[nu]->update();
	}
	
	// Fixar condições iniciais
	for (int nu = 0; nu < n_usinas; nu++)
	{
		FixarCondIniciais(nu);
		modelosGRB[nu]->update();
	}

	// Definir ajustes do solver
	for (int nu = 0; nu < n_usinas; nu++)
	{
		//modelosGRB[nu]->write("probT.lp");	// Escreve modelo em arquivo

		// Nao escrever detalhes
		modelosGRB[nu]->getEnv().set(GRB_IntParam_OutputFlag, 0);
		
		//modelosGRB[nu].getEnv().set(GRB_IntParam_LogToConsole, 0);		// Nao escrever detalhes na tela
		//modelosGRB[nu].getEnv().set(GRB_StringParam_LogFile, "log.txt");		// Escrever log em arquivo

		modelosGRB[nu]->getEnv().set(GRB_DoubleParam_MIPGap, 1e-6);	// Define o gap de tolerancia
		//modelosGRB[nu]->getEnv().set(GRB_IntParam_ScaleFlag, 0);		// Desabilita o model scaling
		//modelosGRB[nu]->getEnv().set(GRB_DoubleParam_MarkowitzTol, 0.01);		// The Markowitz tolerance is used to limit numerical error in the simplex algorithm. Specifically, larger values reduce the error introduced in the simplex basis factorization. A larger value may avoid numerical problems in rare situations, but it will also harm performance
		//modelosGRB[nu]->getEnv().set(GRB_IntParam_Threads, 4);
		if (sistema_a->GetFlagInitAproxCT() == 0)
			modelosGRB[nu]->getEnv().set(GRB_IntParam_Presolve, 0);	//Desabilitar o presolve
		//modelosGRB[nu].getEnv().set(GRB_IntParam_PreSparsify, 1);
		//modelosGRB[nu].getEnv().set(GRB_IntParam_Method, 0);

		//modelosGRB[nu].getEnv().set(GRB_IntParam_MIPFocus, 2);
		//modelosGRB[nu].getEnv().set(GRB_DoubleParam_ImproveStartGap, 0.10);

		modelosGRB[nu]->getEnv().set(GRB_DoubleParam_TimeLimit, 600);		// Limita o tempo de resolução do problema
		// modelosGRB[nu]->getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-2);	// Define tolerancia de otimalidade
	}
}
Spcdec3SPT::~Spcdec3SPT(void)
{
	for (size_t i = 0; i < modelosGRB.size(); i++)
		delete modelosGRB[i];
	modelosGRB.clear();
	for (size_t i = 0; i < vars.size(); i++)
		delete vars[i];
	vars.clear();
}

void Spcdec3SPT::SetPrecision( double precision )
{
	for (int nu = 0; nu < n_usinas; nu++)
	{
		if (sistema_a->GetFlagVarBin() == true)
			modelosGRB[nu]->getEnv().set(GRB_DoubleParam_MIPGap, precision);
		else
			modelosGRB[nu]->getEnv().set(GRB_DoubleParam_OptimalityTol, precision);
	}
}
void Spcdec3SPT::CriarVariaveis(int nu)
{
	try 
	{
		// variáveis em cada problema somente para uma usina termeletrica
		// Estagio 1
		//x = [pt u cp F]
		//x = [pt u up ud F]	"modern" model
		int T = sistema_a->GetTt1();
		for (int t = 0; t < T; t++)
		{
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				modelosGRB[nu]->addVar(0, double (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()), 0.0, GRB_CONTINUOUS, "");		//pt
			else
				modelosGRB[nu]->addVar(0, double (sistema_a->termeletricasVtr[nu].GetPmax()), 0.0, GRB_CONTINUOUS, "");		//pt
			if (sistema_a->GetFlagVarBin() == true)
				modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");		//u
			else
				modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");	//u
			if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
			{
				if (sistema_a->GetFlagVarBin() == true)	
				{
					modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");		//up
					modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");		//ud
				}
				else
				{
					modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");		//up
					modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");		//ud
				}
			}
			else
				modelosGRB[nu]->addVar(0, sistema_a->termeletricasVtr[nu].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");		//cp
			if (sistema_a->GetFlagInitAproxCT() > 1)	//F
				modelosGRB[nu]->addVar(0, double (sistema_a->termeletricasVtr[nu].CustoOperacao(sistema_a->termeletricasVtr[nu].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
		}
		// Estagio 2
		//x = [pt u cp F]
		T = sistema_a->GetTt2() - sistema_a->GetTt1();
		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
		{
			for (int t = 0; t < T; t++)
			{
				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
					modelosGRB[nu]->addVar(0, double (sistema_a->termeletricasVtr[nu].GetPmax() - sistema_a->termeletricasVtr[nu].GetPmin()), 0.0, GRB_CONTINUOUS, "");		//pt
				else
					modelosGRB[nu]->addVar(0, double (sistema_a->termeletricasVtr[nu].GetPmax()), 0.0, GRB_CONTINUOUS, "");		//pt
				if (sistema_a->GetFlagVarBin() == true)
					modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");		//u
				else
					modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");		//u

				if (sistema_a->GetFlagTbinaryModel() == 1)		// "modern" model
				{
					if (sistema_a->GetFlagVarBin() == true)
					{
						modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");		//up
						modelosGRB[nu]->addVar(0, 1, 0.0, GRB_BINARY, "");		//ud
					}
					else
					{
						modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");		//up
						modelosGRB[nu]->addVar(0, 1, 0.0, GRB_CONTINUOUS, "");		//ud
					}
				}
				else
					modelosGRB[nu]->addVar(0, sistema_a->termeletricasVtr[nu].GetCoefCustoPartida(), 0.0, GRB_CONTINUOUS, "");		//cp
				if (sistema_a->GetFlagInitAproxCT() > 1)	//F
					modelosGRB[nu]->addVar(0, double (sistema_a->termeletricasVtr[nu].CustoOperacao(sistema_a->termeletricasVtr[nu].GetPmax(), 1)), 0.0, GRB_CONTINUOUS, "");
			}
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
void Spcdec3SPT::CriarRestricoes(int nu)
{
	try
	{
		vetorint LimTipo;
		vetorfloat LimValor;
		LimTipo.resize(0);LimValor.resize(0);
		CMatrizEsparsa MM(0, n[nu]);

		MatrizRestricoesLineares(n[nu], MM, nu);
		MatrizLimitesLineares(n[nu], &LimTipo, &LimValor, nu);

		GRBLinExpr restricao;
		double coeficiente;
		GRBVar variavel;
		int l;
		for (int i = 0; i < MM.GetNlin(); i++)		// loop no numero de restrições
		{
			l = MM.GetValorLprim(i);
			while ( l != -1 )
			{
				coeficiente = MM.GetValorVal(l);
				variavel = vars[nu][MM.GetValorCol(l)];
				restricao.addTerms( &coeficiente, &variavel, 1);
				l = MM.GetValorLprox(l);
			}
			switch (LimTipo[i]) {
			case 0:
				modelosGRB[nu]->addConstr(restricao, GRB_LESS_EQUAL, LimValor[i], "");
				break;
			case 1:
				modelosGRB[nu]->addConstr(restricao, GRB_EQUAL, LimValor[i], "");
				break;
			case 2:
				modelosGRB[nu]->addConstr(restricao, GRB_GREATER_EQUAL, LimValor[i], "");
				break;
			default:
				cout << "Tipo inválido de restrição adicionada" << endl;
			}
			restricao.clear();
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
void Spcdec3SPT::CriarFuncaoObjetivoRL(const vetorfloat * const lambda, int nu)
{
	try 
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int nt = (n[nu] / T);
		int nt_dual = lambda->size() / T;
		int delta_dual = 0;
		double deltaT;
		int delta = 0;
		int cen = 0;
		
		// Termos lineares (Custo de operação e partida + termos da RL)
		// Termos da RL
		// lambda [pt ph v d phmax]
		for (int t = 0; t < T; t++)
		{
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				deltaT = sistema_a->GetDeltaT1();
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					vars[nu][3 + delta].set(GRB_DoubleAttr_Obj, 1*deltaT);						//F
					vars[nu][delta].set(GRB_DoubleAttr_Obj, - lambda->at(nu + delta_dual));		//pt
				}
				else
				{
					vars[nu][delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1)*deltaT - lambda->at(nu + delta_dual));			//pt
					vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(0)*deltaT);		//u
				}
				vars[nu][2 + delta].set(GRB_DoubleAttr_Obj, 1);		//cp
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				//	if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				//	{
				deltaT = sistema_a->GetDeltaT2();
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					vars[nu][3 + delta].set(GRB_DoubleAttr_Obj, 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//F
					vars[nu][delta].set(GRB_DoubleAttr_Obj, - lambda->at(nu + delta_dual));		//pt
				}
				else
				{
					vars[nu][delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) - lambda->at(nu + delta_dual));		//pt
					vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//u
				}
				vars[nu][2 + delta].set(GRB_DoubleAttr_Obj, 1 * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//cp
			}
			delta += nt;
			delta_dual += nt_dual;
		}
		modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.

		// Termos quadráticos (tem-se que usar GRBQuadExpr, que é mais devagar, mas é o único jeito de incluir termos quadráticos)
		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			GRBQuadExpr fo;
			//x = [pt u cp F]
			double * coeffs;
			coeffs = new double[n[nu]];
			for (int i = 0; i < n[nu]; i++)
				coeffs[i] = vars[nu][i].get(GRB_DoubleAttr_Obj);
			fo.addTerms(coeffs, vars[nu], n[nu]);
			delete coeffs;
			double coeficiente;
			GRBVar variavel;

			delta = 0;
			for (int t = 0; t < T; t++)
			{
				if ((0 <= t) && (sistema_a->GetTt1() > t))
				{
					deltaT = sistema_a->GetDeltaT1();
					coeficiente = sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*deltaT;
					variavel = vars[nu][delta];
					fo.addTerms(&coeficiente, &variavel, &variavel, 1);
				}
				else
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
					//	if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
					//	{
					deltaT = sistema_a->GetDeltaT2();
					coeficiente = sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);
					variavel = vars[nu][delta];
					fo.addTerms(&coeficiente, &variavel, &variavel, 1);
				}
				delta += nt;
			}
			modelosGRB[nu]->setObjective(fo, GRB_MINIMIZE);
			modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.

			// o algoritmo acima pode alternativamente ser dado por:
			//double * coeffs;
			//coeffs = new double[n[nu]];
			//for (int i = 0; i < n[nu]; i++)
			//	coeffs[i] = vars[nu][i].get(GRB_DoubleAttr_Obj);
			//fo.addTerms(coeffs, vars[nu], n[nu]);
			//for (int i = 0; i < n[nu]; i++)
			//	coeffs[i] = 0;

			//delta = 0;
			//for (int t = 0; t < T; t++)
			//{
			//	if ((0 <= t) && (sistema_a->GetTt1() > t))
			//	{
			//		deltaT = sistema_a->GetDeltaT1();
			//		coeffs[delta] = sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*deltaT;		//pt^2
			//	}
			//	else
			//		for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
			//			if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
			//			{
			//				deltaT = sistema_a->GetDeltaT2();
			//				coeffs[delta] = sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);		//pt^2
			//			}
			//	delta += nt;
			//}
			//fo.addTerms(coeffs, vars[nu], vars[nu], n[nu]);
			//delete coeffs;
			//modelosGRB[nu]->setObjective(fo, GRB_MINIMIZE);
		}
		
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
void Spcdec3SPT::CriarFuncaoObjetivoRL2(const vetorfloat * const lambda, int nu)
{
	try 
	{
		int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
		int nt = (n[nu] / T);
		int nt_dual = lambda->size() / T;
		int delta_dual = 0;
		double deltaT;
		int delta = 0;
		int cen = 0;
		
		// Termos lineares (Custo de operação e partida + termos da RL)
		// Termos da RL
		// lambda [pt ph v d phmax]
		for (int t = 0; t < T; t++)
		{
			if ((0 <= t) && (sistema_a->GetTt1() > t))
			{
				deltaT = sistema_a->GetDeltaT1();
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					vars[nu][4 + delta].set(GRB_DoubleAttr_Obj, 1*deltaT);						//F
					vars[nu][delta].set(GRB_DoubleAttr_Obj, - lambda->at(nu + delta_dual));		//pt
					vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, - lambda->at(nu + delta_dual)*sistema_a->termeletricasVtr[nu].GetPmin());		//u
				}
				else
				{
					vars[nu][delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1)*deltaT - lambda->at(nu + delta_dual));			//pt
					if (sistema_a->GetFlagInitAproxCT() == 0)
						vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[nu].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[nu].GetPmin()*(sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1) + sistema_a->termeletricasVtr[nu].GetPmin()*sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)))*deltaT - lambda->at(nu + delta_dual)*sistema_a->termeletricasVtr[nu].GetPmin());		//u
					else
						vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(0)*deltaT + sistema_a->termeletricasVtr[nu].GetPmin()*sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1)*deltaT - lambda->at(nu + delta_dual)*sistema_a->termeletricasVtr[nu].GetPmin());		//u
				}
				vars[nu][2 + delta].set(GRB_DoubleAttr_Obj,  sistema_a->termeletricasVtr[nu].GetCoefCustoPartida());		//up
			}
			else
			{
				cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
				//for (int cen = 0; cen < sistema_a->GetNCenarios(); cen++)
				//	if ((sistema_a->GetTt1() + cen*(sistema_a->GetTt2() - sistema_a->GetTt1()) <= t) && (sistema_a->GetTt1() + (cen + 1)*(sistema_a->GetTt2() - sistema_a->GetTt1()) > t))
				//	{
				deltaT = sistema_a->GetDeltaT2();
				if (sistema_a->GetFlagInitAproxCT() > 1)
				{
					vars[nu][4 + delta].set(GRB_DoubleAttr_Obj, 1*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//F
					vars[nu][delta].set(GRB_DoubleAttr_Obj, - lambda->at(nu + delta_dual));		//pt
					vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, - lambda->at(nu + delta_dual)*sistema_a->termeletricasVtr[nu].GetPmin());		//u
				}
				else
				{
					vars[nu][delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) - lambda->at(nu + delta_dual));		//pt
					if (sistema_a->GetFlagInitAproxCT() == 0)
						vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, (sistema_a->termeletricasVtr[nu].GetCoefCustoOper(0) + sistema_a->termeletricasVtr[nu].GetPmin()*(sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1) + sistema_a->termeletricasVtr[nu].GetPmin()*sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)))*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) - lambda->at(nu + delta_dual)*sistema_a->termeletricasVtr[nu].GetPmin());		//u
					else
						vars[nu][1 + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoOper(0)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) + sistema_a->termeletricasVtr[nu].GetPmin()*sistema_a->termeletricasVtr[nu].GetCoefCustoOper(1)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen) - lambda->at(nu + delta_dual)*sistema_a->termeletricasVtr[nu].GetPmin());		//u
				}
				vars[nu][2 + delta].set(GRB_DoubleAttr_Obj, sistema_a->termeletricasVtr[nu].GetCoefCustoPartida() * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen));		//up
			}
			delta += nt;
			delta_dual += nt_dual;
		}
		modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.


		if (sistema_a->GetFlagInitAproxCT() == 0)
		{
			GRBQuadExpr fo;
			//x = [pt u cp F]
			double * coeffs;
			coeffs = new double[n[nu]];
			for (int i = 0; i < n[nu]; i++)
				coeffs[i] = vars[nu][i].get(GRB_DoubleAttr_Obj);
			fo.addTerms(coeffs, vars[nu], n[nu]);
			delete coeffs;
			double coeficiente;
			GRBVar variavel, variavel2;

			delta = 0;
			for (int t = 0; t < T; t++)
			{
				if ((0 <= t) && (sistema_a->GetTt1() > t))
				{
					deltaT = sistema_a->GetDeltaT1();
					coeficiente = sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*deltaT;	//pt^2
					variavel = vars[nu][delta];
					fo.addTerms(&coeficiente, &variavel, &variavel, 1);
					coeficiente = 2*sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*sistema_a->termeletricasVtr[nu].GetPmin()*deltaT;	//pt*u
					variavel = vars[nu][delta];
					variavel2 = vars[nu][1 + delta];
					fo.addTerms(&coeficiente, &variavel, &variavel2, 1);
				}
				else
				{
					cen = (t - sistema_a->GetTt1()) / (sistema_a->GetTt2() - sistema_a->GetTt1());
					deltaT = sistema_a->GetDeltaT2();
					coeficiente = sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);	//pt^2
					variavel = vars[nu][delta];
					fo.addTerms(&coeficiente, &variavel, &variavel, 1);
					coeficiente = 2*sistema_a->termeletricasVtr[nu].GetCoefCustoOper(2)*sistema_a->termeletricasVtr[nu].GetPmin()*deltaT * sistema_a->hidreletricasVtr[0].GetProbAfluencia(cen);	//pt*u
					variavel = vars[nu][delta];
					variavel2 = vars[nu][1 + delta];
					fo.addTerms(&coeficiente, &variavel, &variavel2, 1);
				}
				delta += nt;
			}
			modelosGRB[nu]->setObjective(fo, GRB_MINIMIZE);
			modelosGRB[nu]->update();	// Atualiza o modelo Gurobi.
		}
	} 
	catch(GRBException e) 
	{
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
	} 
	catch(...) 
	{
			cout << "Exception during optimization" << endl;	
	}
}
void Spcdec3SPT::FixarCondIniciais(int nu)
{
	int T = sistema_a->GetTt1() + sistema_a->GetNCenarios() * (sistema_a->GetTt2() - sistema_a->GetTt1());
	int c = 0;
	int TUr = 0;
	int t_a = 0;
	if (sistema_a->GetFlagTbinaryModel() == 1)	// "modern" model
	{
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = t;
			
			// Verifica condições iniciais e fixa valores das variáveis caso necessário
			if (sistema_a->termeletricasVtr[nu].GetU0() == 0)
				TUr = max(0, sistema_a->termeletricasVtr[nu].GetTDown() + 1 - sistema_a->termeletricasVtr[nu].GetX0());
			else
				TUr = max(0, sistema_a->termeletricasVtr[nu].GetTUp() + 1 - sistema_a->termeletricasVtr[nu].GetX0());
			if ((TUr >= 1) && (t_a + 1 <= TUr))
			{
				//u
				vars[nu][1 + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[nu].GetU0());
				vars[nu][1 + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[nu].GetU0());
				//up
				vars[nu][2 + c].set(GRB_DoubleAttr_LB, 0);
				vars[nu][2 + c].set(GRB_DoubleAttr_UB, 0);
				//ud
				vars[nu][3 + c].set(GRB_DoubleAttr_LB, 0);
				vars[nu][3 + c].set(GRB_DoubleAttr_UB, 0);
			}
			c += (n[nu] / T);
		}
	}
	else
	{
		for (int t = 0; t < T; t++)
		{
			if (t >= sistema_a->GetTt2())
				t_a = sistema_a->GetTt1() + (t - sistema_a->GetTt1()) % (sistema_a->GetTt2() - sistema_a->GetTt1());
			else
				t_a = t;

			// Verifica condições iniciais e fixa valores das variáveis caso necessário
			if (sistema_a->termeletricasVtr[nu].GetU0() == 0)
				TUr = max(0, sistema_a->termeletricasVtr[nu].GetTDown() + 1 - sistema_a->termeletricasVtr[nu].GetX0());
			else
				TUr = max(0, sistema_a->termeletricasVtr[nu].GetTUp() + 1 - sistema_a->termeletricasVtr[nu].GetX0());
			if ((TUr >= 1) && (t_a + 1 <= TUr))
			{
				//u
				vars[nu][1 + c].set(GRB_DoubleAttr_LB, sistema_a->termeletricasVtr[nu].GetU0());
				vars[nu][1 + c].set(GRB_DoubleAttr_UB, sistema_a->termeletricasVtr[nu].GetU0());
				//cp
				vars[nu][2 + c].set(GRB_DoubleAttr_LB, 0);
				vars[nu][2 + c].set(GRB_DoubleAttr_UB, 0);
			}
			c += (n[nu] / T);
		}
	}
}
int Spcdec3SPT::ResolverProblemaRL(Spcdec3Results * resultadoGurobi, const vetorfloat * const lambda, const int iter, int usina)
{
	// Resolve o subproblema de uma usina especifica
	try
	{
		// Criar função objetivo, a seleçao dos lambdas de interesse é feita dentro da funçao abaixo
		if (flag7 == 0)
			CriarFuncaoObjetivoRL(lambda, usina);
		else
			CriarFuncaoObjetivoRL2(lambda, usina);
		
		modelosGRB[usina]->reset();

		// Otimizar
		modelosGRB[usina]->optimize();
		//if (usina == 6)
		//	cout << "Sol.SPT, t=1, i=" << usina << ": pt=" << double(vars[usina][0].get(GRB_DoubleAttr_X)) << ": u=" << double(vars[usina][1].get(GRB_DoubleAttr_X)) << ": F=" << double(vars[usina][4].get(GRB_DoubleAttr_X)) << endl;
		// ???


		// Salvar resultados
		x_i.clear();
		x_i.resize(n[usina]);
		int nStatus = modelosGRB[usina]->get(GRB_IntAttr_Status);
		if (nStatus != 2)
			cout << "The T subproblem was not solved until optimality! " << nStatus << endl;
		if ((nStatus == 2) || (nStatus == 9))
		{
			fo_i = double(modelosGRB[usina]->get(GRB_DoubleAttr_ObjVal));
			for (int i = 0; i < n[usina]; i++)
			{	
				x_i[i] = double(vars[usina][i].get(GRB_DoubleAttr_X));
			}
		}
		else
		{
			for (int i = 0; i < n[usina]; i++)
			{	
				x_i[i] = 0;
			}
			fo_i = 0;
		}
		resultadoGurobi->GravarSolucao(fo_i, x_i, nStatus, resultadoGurobi->GetCH() + resultadoGurobi->GetN() * resultadoGurobi->GetR() + usina);

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		resultadoGurobi->GravarSolucao(fo_i, x_i, e.getErrorCode(), resultadoGurobi->GetCH() + resultadoGurobi->GetN() * resultadoGurobi->GetR() + usina);
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
return 0;
}
// Define caracteristicas da classe sistema
#include "Sistema.h"
#include <ppl.h>
using namespace Concurrency;

MatrixXd CombinaVetores( MatrixXd A_a, VectorXd b_a)
{
	int linhasA = int (A_a.rows());
	int colunasA = int (A_a.cols());
	int linhasb = int (b_a.size());
	MatrixXd M = MatrixXd::Zero(linhasb*linhasA, 1 + colunasA);
	for (int i = 0; i < linhasb*linhasA; i++)
	{
		int l1 = int(floor(double((i)/linhasb)));
		int l2 = int(fmod(double(i), double(linhasb)));
		int j = 0;
		for ( j = 0; j < colunasA; j++)
			M(i,j) = A_a(l1,j);
		M(i,j) = b_a(l2);
	}
	return M;
}

vetorfloat2 CSistema::GerarCoeficientes(CHidreletrica * hidreletricaPtr_a, CUnidades * grupoVtr_a)		// Gera coeficientes para o grupo!!
{
	// Gerar varios planos para cada grupo, aproximados com Taylor de primeira ordem
	// Para considerar (aproximar) para uma função quadrática, ela deve ser SDP ou DP, o que não acontece com a função de produção!
	// Modelo da função

	// Esse algoritmo gera coeficientes para N planos lineares para aproximar a função de produção dos grupos de unidades geradoras
	// Coeficientes gerados para tres variáveis, v (volume armazenado), q (vazao turbinada no grupo) e de (defluencia extra, que é toda vazao defluente, além de 'q' que influencia no canal de fuga)
	// Se 's' influencia no canal de fuga: de = d - q;
	// Se 's' nao influencia no canal de fuga: de = d - q - s;
	// Os planos gerados são:
	// phg <= fphg_n(v, q, de) = a0_n + a1_n * v + a2_n * q + a3_n * de; para n = 1,...,N
	// No modelo esses planos são aplicados como:
	// phg <= fphg_n(v, q, d, s, z) = b0_n + b1_n * v + b2_n * q + b3_n * ds + phg_max * (1 - z); para n = 1,...,N
	// Em que:
	// b0_n = a0_n;
	// b1_n = a1_n;
	// b2_n = a2_n - a3_n;
	// b3_n = a3_n;
	// ds = d - s (para usinas q o vert. n influecia a jusante) ou ds = d
	// Para as usinas em que o vertimento não influencia a cota de jusante phg <= fphg_n(v, q, d, z), e b4_n = 0;
	// O termo phg_max * (1 - z) deve ser adicionado aos planos para evitar inviabilidades (fphg_n < 0), quando o grupo está desligado
	// ou seja, b0_n + b1_n * v pode ser menor que 0, assim soma-se o termo adicional para evitar fphg_n negativos qdo z = 0
	// Assim esse termo só é adicionado quando (b0_n + b1_n * v_min < 0) ou (b0_n + b1_n * v_max < 0)

	// Ajustes
	// -----------------------
	double var_v0 = 10;		// (na tese = 10) Variação (%) em torno do volume inicial para a construção dos planos (menor e maior volumes considerados para aproximação), poderia levar em conta o menor Vmin e o maior Vmax determinados em DeterminarLimitesV()
	int n_pontos_v = 3;		// numero de pontos considerados na aproximação, (n_pontos_ - 1) é o numero de quadrantes considerados
	int n_pontos_q = 5;
	int n_pontos_de = 2;
	int res_planos = 20;	// numero de pontos considerados em cada dimensão para construção do plano (talvez pudesse considerar um numero diferente para cada tipo de variavel, maior para a vazao por exemplo)
	// -----------------------

	// Limites de v
	// Variação de volume em +/- var_v0 % do volume inicial
	double vmin = double (hidreletricaPtr_a->GetV0() * (1 - var_v0/100));
	if (vmin < hidreletricaPtr_a->GetVmin())
		vmin = hidreletricaPtr_a->GetVmin();
	double vmax = double (hidreletricaPtr_a->GetV0() * (1 + var_v0/100));
	if (vmax > hidreletricaPtr_a->GetVmax())
		vmax = hidreletricaPtr_a->GetVmax();
	// Limites de q
	double qmin = grupoVtr_a->GetQmin();
	double qmax = grupoVtr_a->GetQmax() * grupoVtr_a->GetNUnidades();
	// Limites de de
	double demin = 0;
	double demax = 0;	// Considerar em demax somente os grupos diferentes do atual e o vertimento (qdo aplicável)
	for (int j = 0; j < hidreletricaPtr_a->GetNGrupos(); j++)
	{
		if (&hidreletricaPtr_a->grupoVtr[j] == grupoVtr_a)
			continue;
		demax += hidreletricaPtr_a->grupoVtr[j].GetQmax() * hidreletricaPtr_a->grupoVtr[j].GetNUnidades();
	}
	if (hidreletricaPtr_a->GetInflueVert() == 1)
		demax += hidreletricaPtr_a->GetSmax();		// Se demax = 0, a função não depende de 'de'!
	// 
	int nplanos = (n_pontos_v-1)*(n_pontos_q-1)*(n_pontos_de-1);
	double delta_v = (vmax - vmin)/(n_pontos_v-1);
	double delta_q = (qmax - qmin)/(n_pontos_q-1);
	double delta_de = (demax - demin)/(n_pontos_de-1);
	// Construção do planos
	VectorXd vm, qm, dem;
	vetorfloat2 coef;
	MatrixXd X;
	#if (FPH_APP == 1)		// Minimos quadrados (usa-se inversão de matrizes e calculam-se todos os pontos do vetor de pg)
		for (int iv = 0; iv < (n_pontos_v - 1); iv++)
		{
			for (int iq = 0; iq < (n_pontos_q - 1); iq++)
			{
				for (int ide = 0; ide < (n_pontos_de - 1); ide++)
				{
					vm.setLinSpaced(res_planos, vmin+iv*delta_v, vmin+(iv+1)*delta_v);
					qm.setLinSpaced(2*res_planos, qmin+iq*delta_q, qmin+(iq+1)*delta_q);		// como a nao linearidade em q é maior, utilizam-se 2x mais pontos (além do vetor de vazoes já possuir mais divisões, n_pontos_q maior)
					dem.setLinSpaced(res_planos, demin+ide*delta_de, demin+(ide+1)*delta_de);
					
					X = CombinaVetores(vm, qm);
					X = CombinaVetores(X, dem);
					VectorXd v, q, de;
					v = X.col(0);
					q = X.col(1);
					de = X.col(2);		// aqui de é a defluencia total que influencia no canal de fuga (na cota de jusante)

					// Calcular pg para todos os pontos de v, q e de
					int i_max = int (X.rows());
					VectorXd pg = VectorXd::Zero(i_max);
					for (int i = 0; i < i_max; i++)
						pg(i) = hidreletricaPtr_a->GetGeracaoOtima( v(i), q(i), de(i), grupoVtr_a->GetGrupo() - 1);		// retorna a geração ótima para diferentes numeros de unidades, qdo for o caso

					VectorXd I;
					I = VectorXd::Ones(i_max);
					// Algumas usinas tem a função de cota montante e/ou cota de jusante constantes, neste caso alguns coeficientes são nulos
					// Se fcm = constante, coef_v = 0	(vmin == vmax)
					// Se fcj = constante, coef_d = 0
					// Além disso, se a função não depender da defluencia extra a3_n = 0, n precisa incluir no sistema linear;
					// Se demax = 0, a função não depende de 'de'!
					VectorXd coef_v;
					vetorfloat coef_a;
					if ((hidreletricaPtr_a->FcmIsCte() == true) && (hidreletricaPtr_a->FcjIsCte() == true))
					{
						X.resize(v.size(), 2);
						X.col(0) = I;
						X.col(1) = q;
						coef_v = X.householderQr().solve(pg);
						coef_a.push_back(coef_v(0));	// a0
						coef_a.push_back(0);			// a1
						coef_a.push_back(coef_v(1));	// a2
						coef_a.push_back(0);			// a3
					}
					else if ((hidreletricaPtr_a->FcmIsCte() == true) && (hidreletricaPtr_a->FcjIsCte() == false))
					{
						if ( demax == 0 )
						{
							X.resize(v.size(), 2);
							X.col(0) = I;
							X.col(1) = q;
							coef_v = X.householderQr().solve(pg);
							coef_a.push_back(coef_v(0));	// a0
							coef_a.push_back(0);			// a1
							coef_a.push_back(coef_v(1));	// a2
							coef_a.push_back(0);			// a3
						}
						else
						{
							X.resize(v.size(), 3);
							X.col(0) = I;
							X.col(1) = q;
							X.col(2) = de;
							coef_v = X.householderQr().solve(pg);
							coef_a.push_back(coef_v(0));				// a0
							coef_a.push_back(0);						// a1
							coef_a.push_back(coef_v(1) - coef_v(2));	// a2
							coef_a.push_back(coef_v(2));				// a3
						}
					}
					else if ((hidreletricaPtr_a->FcmIsCte() == false) && (hidreletricaPtr_a->FcjIsCte() == true))
					{
						X.resize(v.size(), 3);
						X.col(0) = I;
						X.col(1) = v;
						X.col(2) = q;
						coef_v = X.householderQr().solve(pg);
						coef_a.push_back(coef_v(0));	// a0
						coef_a.push_back(coef_v(1));	// a1
						coef_a.push_back(coef_v(2));	// a2
						coef_a.push_back(0);			// a3
					}
					else
					{
						if ( demax == 0 )
						{
							X.resize(v.size(), 3);
							X.col(0) = I;
							X.col(1) = v;
							X.col(2) = q;
							coef_v = X.householderQr().solve(pg);
							coef_a.push_back(coef_v(0));	// a0
							coef_a.push_back(coef_v(1));	// a1
							coef_a.push_back(coef_v(2));	// a2
							coef_a.push_back(0);			// a3
						}
						else
						{
							X.resize(v.size(), 4);
							X.col(0) = I;
							X.col(1) = v;
							X.col(2) = q;
							X.col(3) = de;
							coef_v = X.householderQr().solve(pg);
							//coef_v = X.colPivHouseholderQr().solve(pg);
							coef_a.push_back(coef_v(0));				// a0
							coef_a.push_back(coef_v(1));				// a1
							coef_a.push_back(coef_v(2)-coef_v(3));		// a2
							coef_a.push_back(coef_v(3));				// a3
						}
					}
					coef.push_back(coef_a);
				}
			}
		}
	#else		// Taylor, só calcula-se pg para os pontos de interesse, portanto res_planos pode ser maior
		double der_v, der_q, der_de, pgf, pgi, der_q_a;
		for (int iv = 0; iv < (n_pontos_v - 1); iv++)
		{
			for (int ide = 0; ide < (n_pontos_de - 1); ide++)
			{
				der_q_a = 0;
				//for (int iq = 0; iq < (n_pontos_q - 1); iq++)
				for (int iq = (n_pontos_q - 2); iq >= 0; iq--)
				{
					vm.setLinSpaced(res_planos, vmin+iv*delta_v, vmin+(iv+1)*delta_v);
					qm.setLinSpaced(2*res_planos, qmin+iq*delta_q, qmin+(iq+1)*delta_q);		// como a nao linearidade em q é maior, utilizam-se 2x mais pontos (além do vetor de vazoes já possuir mais divisões, n_pontos_q maior)
					dem.setLinSpaced(res_planos, demin+ide*delta_de, demin+(ide+1)*delta_de);
					// aqui não precisa-se calcular pg para todos os pontos, somente para os pontos em que se calcula a derivada
					pgi = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos - 1), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
					//cout << vm(res_planos/2 - 1) << ":" << qm(res_planos - 1) << ":"  << dem(res_planos/2 - 1) << ":" << pgi << endl;
					vetorfloat coef_a;
					if ((hidreletricaPtr_a->FcmIsCte() == true) && (hidreletricaPtr_a->FcjIsCte() == true))
					{
						der_v = 0;
						pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
						der_q = (pgf - pgi)/(qm(res_planos) - qm(res_planos - 1));
						der_de = 0;
						coef_a.push_back(pgi - der_v*vm(res_planos/2 - 1) - der_q*qm(res_planos - 1) - der_de*dem(res_planos/2 - 1));				// a0
						coef_a.push_back(der_v);				// a1
						coef_a.push_back(der_q - der_de);		// a2
						coef_a.push_back(der_de);				// a3
					}
					else if ((hidreletricaPtr_a->FcmIsCte() == true) && (hidreletricaPtr_a->FcjIsCte() == false))
					{
						if ( demax == 0 )
						{
							der_v = 0;
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
							der_q = (pgf - pgi)/(qm(res_planos) - qm(res_planos - 1));
							der_de = 0;
							coef_a.push_back(pgi - der_v*vm(res_planos/2 - 1) - der_q*qm(res_planos - 1) - der_de*dem(res_planos/2 - 1));				// a0
							coef_a.push_back(der_v);				// a1
							coef_a.push_back(der_q - der_de);		// a2
							coef_a.push_back(der_de);				// a3
						}
						else
						{
							der_v = 0;
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
							der_q = (pgf - pgi)/(qm(res_planos) - qm(res_planos - 1));
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos - 1), dem(res_planos/2), grupoVtr_a->GetGrupo() - 1);
							der_de = (pgf - pgi)/(dem(res_planos/2) - dem(res_planos/2 - 1));
							coef_a.push_back(pgi - der_v*vm(res_planos/2 - 1) - der_q*qm(res_planos - 1) - der_de*dem(res_planos/2 - 1));				// a0
							coef_a.push_back(der_v);				// a1
							coef_a.push_back(der_q - der_de);		// a2
							coef_a.push_back(der_de);				// a3
						}
					}
					else if ((hidreletricaPtr_a->FcmIsCte() == false) && (hidreletricaPtr_a->FcjIsCte() == true))
					{
						pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2), qm(res_planos - 1), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
						der_v = (pgf - pgi)/(vm(res_planos/2) - vm(res_planos/2 - 1));
						pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
						der_q = (pgf - pgi)/(qm(res_planos) - qm(res_planos - 1));
						der_de = 0;
						coef_a.push_back(pgi - der_v*vm(res_planos/2 - 1) - der_q*qm(res_planos - 1) - der_de*dem(res_planos/2 - 1));				// a0
						coef_a.push_back(der_v);				// a1
						coef_a.push_back(der_q - der_de);		// a2
						coef_a.push_back(der_de);				// a3
					}
					else
					{
						if ( demax == 0 )
						{
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2), qm(res_planos - 1), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
							der_v = (pgf - pgi)/(vm(res_planos/2) - vm(res_planos/2 - 1));
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
							der_q = (pgf - pgi)/(qm(res_planos) - qm(res_planos - 1));
							//cout << qm(res_planos) << ":" << qm(res_planos - 1) << ":" << pgf << ":" << pgi << ":" << der_q << endl;
							der_de = 0;
							coef_a.push_back(pgi - der_v*vm(res_planos/2 - 1) - der_q*qm(res_planos - 1) - der_de*dem(res_planos/2 - 1));				// a0
							coef_a.push_back(der_v);				// a1
							coef_a.push_back(der_q - der_de);		// a2
							coef_a.push_back(der_de);				// a3
						}
						else
						{
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2), qm(res_planos - 1), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
							der_v = (pgf - pgi)/(vm(res_planos/2) - vm(res_planos/2 - 1));
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos), dem(res_planos/2 - 1), grupoVtr_a->GetGrupo() - 1);
							der_q = (pgf - pgi)/(qm(res_planos) - qm(res_planos - 1));
							pgf = hidreletricaPtr_a->GetGeracaoOtima( vm(res_planos/2 - 1), qm(res_planos - 1), dem(res_planos/2), grupoVtr_a->GetGrupo() - 1);
							der_de = (pgf - pgi)/(dem(res_planos/2) - dem(res_planos/2 - 1));
							coef_a.push_back(pgi - der_v*vm(res_planos/2 - 1) - der_q*qm(res_planos - 1) - der_de*dem(res_planos/2 - 1));				// a0
							coef_a.push_back(der_v);				// a1
							coef_a.push_back(der_q - der_de);		// a2
							coef_a.push_back(der_de);				// a3
						}
					}
					#if (FPH_REM)
						if (der_q >= der_q_a)
							coef.push_back(coef_a);
						else
							nplanos--;
					#else
						coef.push_back(coef_a);
					#endif
					der_q_a = der_q;
				}
			}
		}
	#endif

	// Calculo do erro
	double fun, erro;
	vm.setLinSpaced(2*res_planos, vmin, vmax);
	qm.setLinSpaced(4*res_planos, qmin, qmax);		// como a nao linearidade em q é maior, utilizam-se 2x mais pontos (além do vetor de vazoes já possuir mais divisões, n_pontos_q maior)
	dem.setLinSpaced(2*res_planos, demin, demax);
	
	X = CombinaVetores(vm, qm);
	if (demax != 0)
		X = CombinaVetores(X, dem);
	erro = 0;
	for (int i = 0; i < X.rows(); i++)
	{
		// parei aqui
		// avaliar o código abaixo, deve ter condição se o canal de fuga depende de 's'?!?
		// na classe unidade a função PotenciaGerada() tb deve ter outra considerando o 's'???

		// na classe grupos é armazenado o coeficiente referente à cada variavel do modelo, se o 's' n influencia naquela usina, seu coeficiente já vai ser 0
		// entao soh precisa uma função q calcula a potencia gerada (na classe unidades), considerando todas variáveis do modelo (v, q, d, s, z)

		//fun = 1e10;
		//for (int ii = 0; ii < nplanos; ii++)
		//	fun = std::min(fun, coef[ii][0] + coef[ii][1]*X(i,0) + coef[ii][2]*X(i,1) + coef[ii][3]*X(i,2));
		// de = d - q - s (qdo s n influencia na jusante) ou = d - q; pois d eh sempre = Q + s
		// portanto 'de + q' vai ser 'd - s' qdo 's' n influencia na jusante e
		// de + q = d qdo o 's' influencia na jusante
		//fun = grupoVtr_a->PotenciaGerada(X(i,0), X(i,1), X(i,2)+X(i,1));		// essa funcao n pode ser usada aqui pois os coefiecientes da phg ainda n foram setados
		fun = 1e10;

		// abaixo, como estamos trablhando com as var. v,q e de, o coeficiente de q deve ser a soma b2_n + b3_n

		if (demax != 0)
		{
			for (size_t ii = 0; ii < coef.size(); ii++)
				fun = std::min(fun, coef[ii][0] + coef[ii][1]*X(i,0) + (coef[ii][2]+coef[ii][3])*X(i,1) + coef[ii][3]*X(i,2));
			erro += pow((hidreletricaPtr_a->GetGeracaoOtima(X(i,0), X(i,1), X(i,2), grupoVtr_a->GetGrupo() - 1) - fun), 2);
			// debug
			//if (hidreletricaPtr_a->GetIdentUsina() + 1 == 4)
			//	cout << "v = " << X(i,0) << " q = " << X(i,1) << " d = " << X(i,2)+X(i,1) << " pg_ap = " << fun << " pg = " << hidreletricaPtr_a->GetGeracaoOtima(X(i,0), X(i,1), X(i,2), grupoVtr_a->GetGrupo() - 1) << endl;
		}
		else
		{
			for (size_t ii = 0; ii < coef.size(); ii++)
				fun = std::min(fun, coef[ii][0] + coef[ii][1]*X(i,0) + coef[ii][2]*X(i,1));
				//fun = std::min(fun, coef[ii][0] + coef[ii][1]*X(i,0) + coef[ii][2]*X(i,1) + coef[ii][3]*X(i,1));
			erro += pow((hidreletricaPtr_a->GetGeracaoOtima(X(i,0), X(i,1), 0, grupoVtr_a->GetGrupo() - 1) - fun), 2);
			// debug
			//cout << "v = " << X(i,0) << " q = " << X(i,1) << " d = " << 0+X(i,1) << " pg_ap = " << fun << " pg = " << hidreletricaPtr_a->GetGeracaoOtima(X(i,0), X(i,1), 0, grupoVtr_a->GetGrupo() - 1) << endl;
		}
		//cout << endl;
	}
	erro = sqrt(erro/X.rows());
	erro = 100*erro / (grupoVtr_a->GetPmax() * grupoVtr_a->GetNUnidades());

	//// debug
	//for (size_t ii = 0; ii < coef.size(); ii++)
	//	cout << coef[ii][0] << " : " << coef[ii][1] << " : " << coef[ii][2] << " : " << coef[ii][3] << endl;
	//cout << erro / 100 * (grupoVtr_a->GetPmax() * grupoVtr_a->GetNUnidades()) << " MW : " << erro << " % " << endl;
	//cout << endl;

	if (erro > 10)		// erro maior que 10 %
		cout << "Atencao, erro de "<< erro << "% na linerizacao da funcao de producao do grupo " << grupoVtr_a->GetGrupo() << " usina H" << hidreletricaPtr_a->GetIdentUsina() + 1 << endl;
	
	return coef;
}
vetorfloat CSistema::GerarCoeficientes2(CHidreletrica * hidreletricaPtr_a, CUnidades * grupoVtr_a, CSistema * sistema_a)		// Gera coeficientes para o grupo!!
{
	// Código antigo
	// arrumado:
	// GetPmax para GetQmax

	MatrixXd X;
	VectorXd coef_v;
	double fcm, fcj, h, ro;
	srand(1);

	// Variação de volume em +/- 10% do volume inicial
	double vmin = double (hidreletricaPtr_a->GetV0() * 0.9);
	if (vmin < hidreletricaPtr_a->GetVmin())
		vmin = hidreletricaPtr_a->GetVmin();
	double vmax = double (hidreletricaPtr_a->GetV0() * 1.1);
	if (vmax > hidreletricaPtr_a->GetVmax())
		vmax = hidreletricaPtr_a->GetVmax();
	VectorXd vm;
	vm.setLinSpaced(10, vmin, vmax);
	//vm += VectorXd::Random(10)*(vmax-vmin)/1000;
	
	//fmat vm = linspace<fmat>(double(hidreletricaPtr_a->GetVmin()), double(hidreletricaPtr_a->GetVmax()), 10) + randu<fmat>(10,1);
	VectorXd qm;
	qm.setLinSpaced(100, grupoVtr_a->GetQmin(), grupoVtr_a->GetQmax());
	//qm += VectorXd::Random(100)*(grupoVtr_a->GetQmax()-grupoVtr_a->GetQmin())/1000;
	double dm_max = 0;
	int tam_dm = 0;
	for (int j = 0; j < hidreletricaPtr_a->GetNGrupos(); j++)
	{
		//if (&hidreletricaPtr_a->grupoVtr[j] == grupoVtr_a)
		//	continue;
		dm_max = dm_max + hidreletricaPtr_a->grupoVtr[j].GetQmax() * hidreletricaPtr_a->grupoVtr[j].GetNUnidades();
		tam_dm = tam_dm + hidreletricaPtr_a->grupoVtr[j].GetNUnidades();
	}

	VectorXd dm, d1;
	d1.setLinSpaced(2*tam_dm, 0, dm_max);
	//d1 += VectorXd::Random(2*tam_dm)*dm_max/1000;
	if (hidreletricaPtr_a->GetInflueVert() == 1)
	{
		VectorXd d2;
		d2.setLinSpaced(5, dm_max, dm_max + hidreletricaPtr_a->GetSmax());
		//d2 += VectorXd::Random(5)*hidreletricaPtr_a->GetSmax()/1000;
		dm.resize(d1.size() + d2.size());		
		dm.segment(0, d1.size()) = d1;
		dm.segment(d1.size(), d2.size()) = d2;
	}
	else
	{
		dm.resize(d1.size());
		dm = d1;
	}

	X = CombinaVetores(vm, qm);
	X = CombinaVetores(X, dm);
	VectorXd v, q, d;
	v = X.col(0);
	q = X.col(1);
	d = X.col(2);		// aqui d é a defluencia total que influencia no canal de fuga (na cota de jusante)

	VectorXd vv, qq, dd, pgg, I;
	int i_max = int (X.rows());
	VectorXd pg = VectorXd::Zero(i_max);
	//double Q;

	for (int nn = 0; nn < grupoVtr_a->GetNUnidades(); nn++)
	{
		for (int i = 0; i < i_max; i++)
		{
			//Q = q(i)*(nn+1);


			fcm = hidreletricaPtr_a->NivelMontante( v(i) );
			fcj = hidreletricaPtr_a->NivelJusante ( q(i)*(nn+1) + d(i) );
			h = fcm - fcj - grupoVtr_a->GetCoefPerdasHid()*pow( q(i), 2);
			ro = grupoVtr_a->RendimentoHidraulico( h, q(i));
			pg(i) = double (0.00981*ro*h*q(i)*(nn+1));
			//if (pg(i) < 0)		// nao considerar geração negativa na regressão linear
			//	pg(i) = 0;
		}
		vv.conservativeResize(vv.size() + v.size());
		vv.segment(nn*i_max, i_max) = v;
		qq.conservativeResize(qq.size() + q.size());
		qq.segment(nn*i_max, i_max) = q*double(nn + 1);
		dd.conservativeResize(dd.size() + d.size());
		dd.segment(nn*i_max, i_max) = d;
		pgg.conservativeResize(pgg.size() + pg.size());
		pgg.segment(nn*i_max, i_max) = pg;
	}
	v = vv;
	q = qq;
	d = dd;
	pg = pgg;
	
	// aqui n pega-se o máximo da geração apra duas combinações diferentes???
	// implementar isso na classe de unidades... retorna a maior potencia para determinada vazao turbinada (otimização interna)
	// como se fosse entregue para a usina uma meta de geração ou vazao, e ela própria otimiza sua operação

	I = VectorXd::Ones(i_max*(grupoVtr_a->GetNUnidades()));
	
	// Algumas usinas tem a função de cota montante e/ou cota de jusante constantes, neste caso alguns coeficientes são nulos
	// Se fcm = constante, coef_v = 0	(vmin == vmax)
	// Se fcj = constante, coef_d = 0
	vetorfloat coef;
	if ((hidreletricaPtr_a->FcmIsCte() == true) && (hidreletricaPtr_a->FcjIsCte() == true))
	{
		X.resize(v.size(), 2);
		X.col(0) = I;
		X.col(1) = q;
		coef_v = X.householderQr().solve(pg);
		coef.push_back(coef_v(0));	// a0
		coef.push_back(0);			// a1
		coef.push_back(coef_v(1));	// a2
		coef.push_back(0);			// a3
	}
	else if ((hidreletricaPtr_a->FcmIsCte() == true) && (hidreletricaPtr_a->FcjIsCte() == false))
	{
		X.resize(v.size(), 3);
		X.col(0) = I;
		X.col(1) = q;
		X.col(2) = d;
		coef_v = X.householderQr().solve(pg);
		coef.push_back(coef_v(0));	// a0
		coef.push_back(0);			// a1
		coef.push_back(coef_v(1));	// a2
		coef.push_back(coef_v(2));	// a3
	}
	else if ((hidreletricaPtr_a->FcmIsCte() == false) && (hidreletricaPtr_a->FcjIsCte() == true))
	{
		X.resize(v.size(), 3);
		X.col(0) = I;
		X.col(1) = v;
		X.col(2) = q;
		coef_v = X.householderQr().solve(pg);
		coef.push_back(coef_v(0));	// a0
		coef.push_back(coef_v(1));	// a1
		coef.push_back(coef_v(2));	// a2
		coef.push_back(0);			// a3
	}
	else
	{
		X.resize(v.size(), 4);
		X.col(0) = I;
		X.col(1) = v;
		X.col(2) = q;
		X.col(3) = d;
		coef_v = X.householderQr().solve(pg);
		//coef_v = X.colPivHouseholderQr().solve(pg);
		for (int i = 0; i < int(coef_v.rows()); i++ )
		{
			//if (abs(coef_v(i)) <= tol)
			//	coef.push_back(0);
			//else
			coef.push_back(coef_v(i));
			//cout << coef_v(i) << endl;
			//cout << coef[i] << endl;
		}
	}

	// Calculo erro
	double fun, erro;
	erro = 0;
	for (int i = 0; i < i_max; i++)
	{
		//fun = coef[0] + coef[1]*v(i) + coef[2]*q(i) + coef[3]*d(i) + coef[4]*v(i)*q(i) + coef[5]*v(i)*d(i) + coef[6]*q(i)*d(i) + coef[7]*pow(v(i),2) + coef[8]*pow(q(i),2) + coef[9]*pow(d(i),2);
		fun = coef[0] + coef[1]*v(i) + coef[2]*q(i) + coef[3]*d(i);
		erro = erro + pow((pg(i) - fun),2);
	}
	erro = sqrt(erro/i_max);
	//cout << coef[0] + coef[1]*hidreletricaPtr_a->GetVmin() + coef[2]*grupoVtr_a->GetQmin() + coef[3]*grupoVtr_a->GetQmin() + coef[4]*hidreletricaPtr_a->GetVmin()*grupoVtr_a->GetQmin() + coef[5]*hidreletricaPtr_a->GetVmin()*grupoVtr_a->GetQmin() + coef[6]*grupoVtr_a->GetQmin()*grupoVtr_a->GetQmin() + coef[7]*pow(hidreletricaPtr_a->GetVmin(),2) + coef[8]*pow(grupoVtr_a->GetQmin(),2) + coef[9]*pow(grupoVtr_a->GetQmin(),2) << endl;
	//cout << "Erro : " << erro << endl;

	return coef;
}

CSistema::CSistema(string diretorio_a, string arq_hidreletricas_a, string arq_unidades_a, string arq_termeletricas_a, string arq_parametros_a, string arq_demanda_a, string arq_barras_a, string arq_linhas_a, string arq_afluencia_a, string arq_cond_inic_h_a, string arq_cond_inic_t_a)
{
	diretorio = diretorio_a;
	arq_parametros = diretorio_a+arq_parametros_a;
	CarregarParametros();
	arq_hidreletricas = diretorio_a+arq_hidreletricas_a;
	CriarHidreletricas();
	arq_unidades = diretorio_a+arq_unidades_a;
	CriarGrupos();
	arq_termeletricas = diretorio_a+arq_termeletricas_a;
	CriarTermeletricas();
	arq_demanda = diretorio_a+arq_demanda_a;
	CriarDemandas();
	arq_barras = diretorio_a+arq_barras_a;
	CriarBarras();
	arq_linhas = diretorio_a+arq_linhas_a;
	CriarLinhas();
	arq_afluencia = diretorio_a+arq_afluencia_a;
	CarregarAfluencias();
	arq_cond_inic_h = diretorio_a+arq_cond_inic_h_a;
	CarregarCIH();
	arq_cond_inic_t = diretorio_a+arq_cond_inic_t_a;
	CarregarCIT();
	// Adequação do numero de períodos considerados
	Tt1 = Tt1 / int(delta_t1);
	Tt2 = Tt2 / int(delta_t1);
}
CSistema::~CSistema(void)
{
	hidreletricasVtr.clear();
	termeletricasVtr.clear();
	demandasVtr.clear();
	barrasVtr.clear();
	linhasVtr.clear();
}
void CSistema::CriarHidreletricas()
{
	string nome_usina_a;
	double vmin_a, vmax_a, smin_a, smax_a;
	int n_grupos_a, usina_jusante_a, vert_canal_fuga_a, tempo_viagem_a;
	vetorfloat coef_mont_a, coef_jus_a;
	int ident_H = 0;

	ifstream inFile( arq_hidreletricas, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_hidreletricas << " could not be opened" << endl;
		exit( 1 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> nome_usina_a >> vmin_a >> vmax_a >> smin_a >> smax_a >> vert_canal_fuga_a >> n_grupos_a >> usina_jusante_a >> tempo_viagem_a;
		// Adequar tempo de viagem ao delta_t1 (= delta_t2), transformar de horas para periodos
		tempo_viagem_a = tempo_viagem_a / int (delta_t1);
		coef_mont_a.resize(CNM);
		coef_jus_a.resize(CNJ);
		for (int j=0;j<CNM;j++)
			inFile >> coef_mont_a[j];
		for (int j=0;j<CNJ;j++)
			inFile >> coef_jus_a[j];
		hidreletricasVtr.push_back(CHidreletrica(nome_usina_a, ident_H, vmin_a, vmax_a, smin_a, smax_a, vert_canal_fuga_a, n_grupos_a, usina_jusante_a, tempo_viagem_a, coef_mont_a, coef_jus_a));
		ident_H++;
	}
	// usina a jusante ou montante = 0 significa nenhuma, diferente da identificacao da usina q vai de 0 a H
	for (size_t h = 0; h < hidreletricasVtr.size(); h++)
		if (hidreletricasVtr[h].GetUsinaJusante() != 0)
			hidreletricasVtr[hidreletricasVtr[h].GetUsinaJusante() - 1].SetUsinaMont(h + 1);

	cout << "-> " << ident_H << " hidreletricas adicionadas!" << endl;
	inFile.close();
}
void CSistema::CriarGrupos()
{
	int n_usina_a, grupo_a, n_unidades_a;
	double qmin_a, qmax_a, pmin_a, pmax_a, coef_perdas_hid_a;
	vetorfloat coef_rend_a;
	
	ifstream inFile( arq_unidades, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_unidades << " could not be opened" << endl;
		exit( 2 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> n_usina_a >> grupo_a >> n_unidades_a >> qmin_a >> qmax_a >> pmin_a >> pmax_a >> coef_perdas_hid_a;
		coef_rend_a.resize(CRH);
		for (int j=0;j<CRH;j++)
			inFile >> coef_rend_a[j];
		//hidreletricasVtr[usina_a]->add_grupo(usina_a, grupo_a, n_unidades_a, qmin_a, qmax_a, pmin_a, pmax_a, coef_perdas_hid_a, coef_rend_a);
		hidreletricasVtr[n_usina_a - 1].grupoVtr.push_back(CUnidades(n_usina_a, grupo_a, n_unidades_a, qmin_a, qmax_a, pmin_a, pmax_a, coef_perdas_hid_a, coef_rend_a));
	}
	inFile.close();
}
void CSistema::GerarFpgh()
{
	// for paralelo, fonte: https://msdn.microsoft.com/en-us/library/dd728073(v=vs.100).aspx
	parallel_for (size_t(0), hidreletricasVtr.size(), [&](size_t h)
	{
		for (size_t j = 0; j < hidreletricasVtr[h].grupoVtr.size(); j++)
		{
			hidreletricasVtr[h].grupoVtr[j].SetCoefPhg(GerarCoeficientes(&hidreletricasVtr[h], &hidreletricasVtr[h].grupoVtr[j]));
		}
	});

	//// para debugar precisa-se do código em série
	//for (size_t h = 0; h < hidreletricasVtr.size(); h++)
	//{
	//	for (size_t j = 0; j < hidreletricasVtr[h].grupoVtr.size(); j++)
	//	{
	//		hidreletricasVtr[h].grupoVtr[j].SetCoefPhg(GerarCoeficientes(&hidreletricasVtr[h], &hidreletricasVtr[h].grupoVtr[j]));
	//	}
	//}
	
	////// Exportar coeficientes para plotar e conferir no matlab
	//// Imprimir coeficientes
	//// ------------------------------------------------
	//ofstream * inFile;
	//inFile = new ofstream( "CoeficientesFphg.txt", ios::out );
	//if ( inFile->is_open() )
	//{
	//	for (size_t r = 0; r < hidreletricasVtr.size(); r++)
	//	{
	//		for (size_t j = 0; j < hidreletricasVtr[r].grupoVtr.size(); j++)
	//		{
	//			*inFile << "H" << r+1 << " grupo " << j+1 << endl;
	//			for (int nap = 0; nap < hidreletricasVtr[r].grupoVtr[j].GetNappFPH() ; nap++)
	//				*inFile << std::scientific << setprecision(10) << setw(15) << right << hidreletricasVtr[r].grupoVtr[j].CoefPhg(nap) << char(9) << hidreletricasVtr[r].grupoVtr[j].CoefPhgV(nap) << char(9) << hidreletricasVtr[r].grupoVtr[j].CoefPhgQ(nap) << char(9) << hidreletricasVtr[r].grupoVtr[j].CoefPhgD(nap) << endl;
	//		}
	//	}
	//	inFile->close();
	//}
	//else
	//	cout << "Unable to open file";
	//delete inFile;
	//// ------------------------------------------------
}
void CSistema::CriarTermeletricas()
{
	string nome_usina_a;
	double pmin_a, pmax_a, rampa_down_a, rampa_up_a;
	int t_down_a, t_up_a;
	vetorfloat coef_custo_oper_a, coef_custo_partida_a;
	int ident_T = 0;

	ifstream inFile( arq_termeletricas, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_termeletricas << " could not be opened" << endl;
		exit( 3 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> nome_usina_a >> pmin_a >> pmax_a >> t_down_a >> t_up_a >> rampa_down_a >> rampa_up_a;
		// Adequar tup/tdown ao delta_t1 (= delta_t2), transformar de horas para periodos
		t_down_a = t_down_a / int (delta_t1);
		t_up_a = t_up_a / int (delta_t1);
		// Adequar rampas ao delta_t1 (= delta_t2), transformar de horas para periodos
		rampa_down_a = rampa_down_a * int (delta_t1);
		rampa_up_a = rampa_up_a * int (delta_t1);
		coef_custo_oper_a.resize(COT);
		coef_custo_partida_a.resize(CPT);
		for (int j=0;j<COT;j++)
			inFile >> coef_custo_oper_a[j];
		for (int j=0;j<CPT;j++)
			inFile >> coef_custo_partida_a[j];
		termeletricasVtr.push_back(CTermeletrica( nome_usina_a, ident_T, pmin_a, pmax_a, t_down_a, t_up_a, rampa_down_a, rampa_up_a, coef_custo_oper_a, coef_custo_partida_a));
		ident_T++;
	}
	cout << "-> " << ident_T << " termeletricas adicionadas!" << endl;
	inFile.close();
}
void CSistema::AproximarFCT()
{
	if (flag_init_custoT > 0)
	{
		if (flag_appr_custo >= 2)
			if (flag_init_custoT < 2)
			{
				cout << "PCFD: numero de aproximações para a Função de custo das T é menor que 2, ajustado para 2!!!" << endl;
				flag_init_custoT = 2;
			}

		double deltax;
		int n_pontos;
		double discretiz_x = 0.1;
		VectorXd x, y, I, coef;
		MatrixXd X;
		
		double coeff[2];
		double p_bar;
		double aprox_custo;
		//if (flag_appr_custo >= 2)		// perspective-cut formulation dynamically: criam-se flag_appr_custo cortes, mas inicialmente só se utilizam flag_aprox_custoT.
		//	aprox_custo = flag_appr_custo;
		//else							// perspective-cut formulation statical
		//	aprox_custo = flag_init_custoT;

		aprox_custo = flag_init_custoT;
		for (size_t t = 0; t < termeletricasVtr.size(); t++)
		{
			deltax = (termeletricasVtr[t].GetPmax() - termeletricasVtr[t].GetPmin()) / aprox_custo;
			/// arrumado e ainda não testado, assim como deve ser visto quais cortes iniciais devem ser adicionados!!!
			///deltax = (termeletricasVtr[t].GetPmax() - termeletricasVtr[t].GetPmin()) / (aprox_custo - 1);


			if (flag_appr_custo == 0)
				n_pontos = int (ceil( deltax / discretiz_x));
			for (int i = 0; i < aprox_custo; i++)
			{
				if (flag_appr_custo == 0)		// least-square approximation
				{
					x.resize(n_pontos);
					x.setLinSpaced(n_pontos,termeletricasVtr[t].GetPmin() + i*deltax,termeletricasVtr[t].GetPmin() + (i + 1)*deltax);
					y.resize(n_pontos);
					for (int n = 0; n < x.size(); n++)
						y[n] = termeletricasVtr[t].CustoOperacao(x[n], 1);
					I = VectorXd::Ones(n_pontos);
					X.resize(n_pontos, 2);
					X.col(0) = I;
					X.col(1) = x;
					coef = X.householderQr().solve(y);
					if (flag_init_custoT == 1)
						termeletricasVtr[t].SetCoefCustoOperacao(coef(0), coef(1), 0);
					else
						termeletricasVtr[t].AdicionarReta(coef(0), coef(1));
				}
				else		// Perspective-cut formulation
				{
					p_bar = ( termeletricasVtr[t].GetPmin() + i*deltax + termeletricasVtr[t].GetPmin() + (i + 1)*deltax ) / 2;
					termeletricasVtr[t].AdicionarPontoPt(p_bar);
					coeff[0] = 2*termeletricasVtr[t].GetCoefCustoOper(2)*p_bar + termeletricasVtr[t].GetCoefCustoOper(1);
					coeff[1] = termeletricasVtr[t].GetCoefCustoOper(0) - termeletricasVtr[t].GetCoefCustoOper(2)*pow(p_bar,2);
					termeletricasVtr[t].AdicionarReta(coeff[1], coeff[0]);
				}

				// np é o numero de partes da função de custo
				// se np for 0 ou 1 o coeficiente das térmicas não muda, no caso de 1 o termo quadrático é zerado e o linear e constante atualizados.
				// se np > 1 deve-se criar mais uma variável, F por térmica. A variavel F vai na f.o. e cria-se mais np restrições lineares por térmica.
			}
		}
	}
}
void CSistema::CarregarParametros()
{
	ifstream inFile( arq_parametros, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_parametros << " could not be opened" << endl;
		exit( 4 );
	}
	while ( ! inFile.eof() )
		inFile >> Tt1 >> Tt2 >> delta_t1 >> delta_t2 >> fator_reserva >> barra_ref >> custo_def;
	custo_vfol = 20*custo_def;
	//custo_vfol = (20/0.0036)*custo_def;

	// A implementação só suporta deltas iguais e menores que T!
	if (delta_t1 > Tt1)
	{
		cout << "Atenção: delta_t1 > Tt1, delta_t1 será alterado!" << endl;
		delta_t1 = Tt1;
	}
	if (delta_t2 > (Tt2 - Tt1))
	{
		cout << "Atenção: delta_t2 > Tt2, delta_t2 será alterado!" << endl;
		delta_t2 = (Tt2 - Tt1);
	}
	if (delta_t1 != delta_t2)
	{
		cout << "Atenção: delta_t1 != delta_t2, delta_t1 será usado!" << endl;
		delta_t2 = delta_t1;
	}
	cout << "-> " << "T_1 = " << Tt1 << " ; T_2 = " << Tt2 << endl;
	inFile.close();
}
void CSistema::CriarDemandas()
{
	int n_demanda_a;
	int ultima_demanda = 0;
	//vetorfloat D_1_a, D_2_a;
	vetorfloat D_a;
	double L, L_a;

	ifstream inFile( arq_demanda, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_demanda << " could not be opened" << endl;
		exit( 5 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> n_demanda_a;
		D_a.resize(Tt2 / int(delta_t1));
		for (int j=0;j<Tt2/int(delta_t1);j++)
		{
			L = 0;
			for (int jj=0;jj<int(delta_t1);jj++)
			{
				inFile >> L_a;
				L += L_a;
			}
			D_a[j] = L / int(delta_t1);
		}
		if (n_demanda_a > ultima_demanda)
		{
			ultima_demanda = n_demanda_a;
			demandasVtr.push_back(CDemanda(n_demanda_a));
		}
		demandasVtr[n_demanda_a - 1].SetDemanda(D_a);

		//D_a.resize(Tt2);
		//for (int j=0;j<Tt2;j++)
		//	inFile >> D_a[j];
		//if (n_demanda_a > ultima_demanda)
		//{
		//	ultima_demanda = n_demanda_a;
		//	demandasVtr.push_back(CDemanda(n_demanda_a));
		//}
		//demandasVtr[n_demanda_a - 1].SetDemanda(D_a);
	}
	inFile.close();
}
void CSistema::CriarBarras()
{
	string nome_barra_a;
	string nome_barra_anterior = "AAAA";
	int cont_barras = 0;
	int i_demanda, i_hidro, i_termo, cont_barras_iguais;
	int ident_B = 0;

	ifstream inFile( arq_barras, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_barras << " could not be opened" << endl;
		exit( 6 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> nome_barra_a >> i_demanda >> i_hidro >> i_termo;
		if (nome_barra_anterior == nome_barra_a)
		{
			if (i_demanda != 0)
				barrasVtr[cont_barras - cont_barras_iguais].AddDemanda(&demandasVtr[i_demanda - 1]);
			if (i_hidro != 0)
				barrasVtr[cont_barras - cont_barras_iguais].AddHidro(&hidreletricasVtr[i_hidro - 1]);
			if (i_termo != 0)
				barrasVtr[cont_barras - cont_barras_iguais].AddTermo(&termeletricasVtr[i_termo - 1]);
			cont_barras_iguais++;
			cont_barras--;
		}
		else
		{
			barrasVtr.push_back(CBarra(nome_barra_a, ident_B));
			ident_B++;
			cont_barras_iguais = 1;
			if (i_demanda != 0)
				barrasVtr[cont_barras].AddDemanda(&demandasVtr[i_demanda - 1]);
			if (i_hidro != 0)
				barrasVtr[cont_barras].AddHidro(&hidreletricasVtr[i_hidro - 1]);
			if (i_termo != 0)
				barrasVtr[cont_barras].AddTermo(&termeletricasVtr[i_termo - 1]);
		}
		cont_barras++;
		nome_barra_anterior = nome_barra_a;
	}
	cout << "-> " << ident_B << " barras adicionadas!" << endl;
	inFile.close();
}
void CSistema::CriarLinhas()
{
	string nome_linha_a;
	double reatancia_a, capacidade_a;
	int de_barra, para_barra;
	int ident_L = 0;
	ifstream inFile( arq_linhas, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_linhas << " could not be opened" << endl;
		exit( 7 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> nome_linha_a >> de_barra >> para_barra >> reatancia_a >> capacidade_a;
		linhasVtr.push_back(CLinha(nome_linha_a, &barrasVtr[de_barra - 1], &barrasVtr[para_barra - 1], reatancia_a, capacidade_a, de_barra, para_barra));
		ident_L++;
	}
	cout << "-> " << ident_L << " linhas adicionadas!" << endl;
	inFile.close();
}
void CSistema::CarregarAfluencias()
{
	int n_usina;
	int cenario;
	vetorfloat afluencia_a;
	double prob, y, y_a;
	n_cenarios = 0;
	ifstream inFile( arq_afluencia, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_afluencia << " could not be opened" << endl;
		exit( 8 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> n_usina >> cenario >> prob;
		if (cenario >= n_cenarios)
			n_cenarios = cenario;
		afluencia_a.resize(Tt2/int(delta_t1));
		for (int j = 0; j < Tt2/int(delta_t1); j++)
		{
			y = 0;
			for (int jj=0;jj<int(delta_t1);jj++)
			{
				inFile >> y_a;
				y += y_a;
			}
			afluencia_a[j] = y / int(delta_t1);
		}
		hidreletricasVtr[n_usina - 1].SetAfluencia(afluencia_a);
		hidreletricasVtr[n_usina - 1].SetProbAfluencia(prob/100);

		//inFile >> n_usina >> cenario >> prob;
		//if (cenario >= n_cenarios)
		//	n_cenarios = cenario;
		//afluencia_a.resize(Tt2);
		//for (int j = 0; j < Tt2; j++)
		//	inFile >> afluencia_a[j];
		//hidreletricasVtr[n_usina - 1].SetAfluencia(afluencia_a);
		//hidreletricasVtr[n_usina - 1].SetProbAfluencia(prob/100);
	}
	inFile.close();
}
void CSistema::CarregarCIH()
{
	int n_usina;
	double v0_a, vmeta_a, custo_agua_a;

	ifstream inFile( arq_cond_inic_h, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_cond_inic_h << " could not be opened" << endl;
		exit( 9 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> n_usina >> v0_a >> vmeta_a >> custo_agua_a;
		hidreletricasVtr[n_usina - 1].SetV0(v0_a);
		hidreletricasVtr[n_usina - 1].SetVMeta(vmeta_a);
		hidreletricasVtr[n_usina - 1].SetCustoAgua(custo_agua_a);
	}
	inFile.close();
}
void CSistema::CarregarCIT()
{
	int n_usina;
	double pt0_a;
	int u0_a, x0_a;

	ifstream inFile( arq_cond_inic_t, ios::in );   
	if ( !inFile )                                                            
	{                                                                               
		cout << "File "<< arq_cond_inic_t << " could not be opened" << endl;
		exit( 10 );
	}
	while ( ! inFile.eof() )
	{
		inFile >> n_usina >> pt0_a >> u0_a >> x0_a;
		termeletricasVtr[n_usina - 1].SetPt0(pt0_a);
		termeletricasVtr[n_usina - 1].SetU0(u0_a);
		termeletricasVtr[n_usina - 1].SetX0(x0_a);
	}
	inFile.close();
}
void CSistema::SetFlags(int mod_rede, bool vfol, bool phmax, int init_custoT, bool var_bin, int appr_custo, int Tbinary_model, int preconditioner)
{
	flag_modelo_rede = mod_rede;
	flag_phmax = phmax;
	flag_vfol = vfol;
	flag_init_custoT = init_custoT;
	flag_var_bin = var_bin;
	flag_appr_custo = appr_custo;
	flag_Tbinary_model = Tbinary_model;
	flag_preconditioner = preconditioner;
	AproximarFCT();
	GerarFpgh();
	CalcularReserva();
	if (flag_modelo_rede != 1)	// calcular matriz beta no modelo de barra unica somente para calcular fluxos de potencias dos resultados
	{
		Beta.resize(int(linhasVtr.size()), int(barrasVtr.size()));		// matriz usada no modelo compacto
		CalcularMatrizBeta();

		//BetaS.ZerarMatriz(int(linhasVtr.size()), int(linhasVtr.size()));
		//CalcularMatrizBetaS();
	}
	else
		Beta.resize(0,0);
	DeterminarLimitesV();
}
void CSistema::SetPreconditioners(vetorfloat precond_a)
{
	precond = precond_a;
}
void CSistema::CalcularReserva()
{
	res_gir.resize(Tt1 + n_cenarios*(Tt2 - Tt1));
	size_t B = barrasVtr.size();
	double D;
	int cen;
	for (int t = 0; t < Tt1 + n_cenarios*(Tt2 - Tt1); t++)
	{
		D = 0;
		for (size_t b = 0; b < B; b++)
		{
			if ((0 <= t) && (Tt1 > t))
				for (size_t d = 0; d < barrasVtr[b].demandasPtr.size(); d++)
					D += barrasVtr[b].demandasPtr[d]->GetD(0, t);
			else
			{
				cen = (t - Tt1) / (Tt2 - Tt1);		// retorna um inteiro 0,1,...
				for (size_t d = 0; d < barrasVtr[b].demandasPtr.size(); d++)
					D += barrasVtr[b].demandasPtr[d]->GetD(cen, t - cen*(Tt2 - Tt1));
			}
		}
		res_gir[t] = D*fator_reserva/100;
	}
}
void CSistema::CalcularMatrizBeta()
{
	///clock_t tempo;
	///tempo = clock();

	//double zerotol = 1e-8;
	//double valor;
	// essas matrizes não são tão grandes para serem tratas como sparsas, porém Alin é
	MatrixXd B = MatrixXd::Zero(int(barrasVtr.size()), int(barrasVtr.size()));
	MatrixXd St = MatrixXd::Zero(int(linhasVtr.size()), int(barrasVtr.size()));
	MatrixXd TT = MatrixXd::Zero(int(linhasVtr.size()), int(linhasVtr.size()));

	//SparseMatrix<double> St(int(linhasVtr.size()), int(barrasVtr.size()));
	// algoritmo para formar a matriz de incidencia linhaxbarra e Ybarra matriz B
	int ifr, ito;
	for ( size_t l = 0; l < linhasVtr.size(); l++)
	{
		ifr = linhasVtr[l].GetIndDeBarra() - 1;
		ito = linhasVtr[l].GetIndParaBarra() - 1;
		St(l, ifr) = 1;
		St(l, ito) = -1;
		B(ifr, ifr) = B(ifr, ifr) + 1/(linhasVtr[l].GetReatancia());
		B(ifr, ito) = B(ifr, ito) - 1/(linhasVtr[l].GetReatancia());
		B(ito, ifr) = B(ito, ifr) - 1/(linhasVtr[l].GetReatancia());
		B(ito, ito) = B(ito, ito) + 1/(linhasVtr[l].GetReatancia());
		TT(l, l) = 1 / (linhasVtr[l].GetReatancia());
		// nao precisa multiplicar por 100 pois ao calcular beta fica 100 * 1/100
	}
	//
	//MatrixXd Beta(int(linhasVtr.size()), int(barrasVtr.size()));
	Beta = TT*St;

	// Eliminar linha e coluna de referencia de B antes de invertê-la (em seguida incluí-se uma linha e uma coluna com 0's)
	int refbus = GetBarraRef() - 1;
	// Remover linha
	__int64 numRows = B.rows() - 1;
	__int64 numCols = B.cols();
	B.block(refbus, 0, numRows - refbus, numCols) = B.block(refbus + 1, 0, numRows - refbus, numCols);
	B.conservativeResize(numRows,numCols);
	// Remover coluna
	numRows = B.rows();
	numCols = B.cols() - 1;
	B.block(0, refbus, numRows, numCols - refbus) = B.block(0, refbus + 1, numRows, numCols - refbus);
	B.conservativeResize(numRows,numCols);
	// Inverter matriz
	
	MatrixXd Ieye = MatrixXd::Identity(B.rows(), B.rows());
	//B = B.partialPivLu().solve(Ieye);
	B = B.fullPivLu().solve(Ieye);
	//B = B.inverse();

	// Incluir linha de 0's na posição da barra de referencia
	numRows = B.rows() + 1;
	numCols = B.cols();
	B.conservativeResize(numRows,numCols);
	MatrixXd temp = B.block(refbus, 0, numRows - refbus - 1, numCols);
	B.block(refbus + 1, 0, numRows - refbus - 1, numCols) = temp;
	for (int i = 0; i < B.cols(); i++)
		B(refbus, i) = 0;
	// Incluir coluna de 0's na posição da barra de referencia
	numRows = B.rows();
	numCols = B.cols() + 1;
	B.conservativeResize(numRows,numCols);
	temp = B.block(0, refbus, numRows, numCols - refbus - 1);
	B.block(0, refbus + 1, numRows, numCols - refbus - 1) = temp;
	for (int i = 0; i < B.rows(); i++)
		B(i, refbus) = 0;
	Beta = Beta*B;		// Beta é uma matriz densa!!
	///cout << "Calc. Beta com Eigen: " << double(clock() - tempo)/double(CLOCKS_PER_SEC) << endl;

	//cout << "B, linha 0 : " << endl;
	//for (int l = 0; l < 10; l++)
	//	cout << Beta(0 , l) << ";";
	//cout << "Hrst, Beta final" << endl;
	//std::cin.ignore();

	//// Imprimir Beta em arquivo
	//ofstream * inFile;
	//inFile = new ofstream( "Beta.txt", ios::out );
	//if ( inFile->is_open() )                                                            
	//{
	//	*inFile << std::scientific << setprecision(10);
	//	for (size_t i = 0; i < Beta.rows(); i++)
	//	{
	//		for (size_t j = 0; j < Beta.cols(); j++)
	//			*inFile << Beta(i,j) << char(9);
	//		*inFile << endl;
	//	}
	//	*inFile << endl;
	//}
	//else
	//	cout << "Unable to open file";
	//inFile->close();
	//delete inFile;

}
void CSistema::CalcularMatrizBetaS()
{
	clock_t tempo;
	tempo = clock();

	// essas matrizes não são tão grandes para serem tratas como sparsas, porém Alin é
	MatrixXd B = MatrixXd::Zero(int(barrasVtr.size()), int(barrasVtr.size()));
	CMatrizEsparsa St(int(linhasVtr.size()), int(barrasVtr.size()));
	// algoritmo para formar a matriz de incidencia linhaxbarra e Ybarra matriz B
	int ifr, ito;
	for ( size_t l = 0; l < linhasVtr.size(); l++)
	{
		ifr = linhasVtr[l].GetIndDeBarra() - 1;
		ito = linhasVtr[l].GetIndParaBarra() - 1;
		St.SubstituirElemento(l, ifr, 1);
		St.SubstituirElemento(l, ito, -1);
		B(ifr, ifr) = B(ifr, ifr) + 1/(linhasVtr[l].GetReatancia());
		B(ifr, ito) = B(ifr, ito) - 1/(linhasVtr[l].GetReatancia());
		B(ito, ifr) = B(ito, ifr) - 1/(linhasVtr[l].GetReatancia());
		B(ito, ito) = B(ito, ito) + 1/(linhasVtr[l].GetReatancia());
		BetaS.SubstituirElemento(l, l, 1 / (linhasVtr[l].GetReatancia()));

		//B(ifr, ifr) = B(ifr, ifr) + 100/(linhasVtr[l].GetReatancia());
		//B(ifr, ito) = B(ifr, ito) - 100/(linhasVtr[l].GetReatancia());
		//B(ito, ifr) = B(ito, ifr) - 100/(linhasVtr[l].GetReatancia());
		//B(ito, ito) = B(ito, ito) + 100/(linhasVtr[l].GetReatancia());
		//BetaS.SubstituirElemento(l, l, 100 / (linhasVtr[l].GetReatancia()));
	}
	//
	//MatrixXd Beta(int(linhasVtr.size()), int(barrasVtr.size()));
	BetaS.MultiplicarPorMatriz(&St);
	// Eliminar linha e coluna de referencia de B antes de invertê-la (em seguida incluí-se uma linha e uma coluna com 0's)
	int refbus = GetBarraRef() - 1;
	// Remover linha
	int numRows = B.rows() - 1;
	int numCols = B.cols();
	B.block(refbus, 0, numRows - refbus, numCols) = B.block(refbus + 1, 0, numRows - refbus, numCols);
	B.conservativeResize(numRows,numCols);
	// Remover coluna
	numRows = B.rows();
	numCols = B.cols() - 1;
	B.block(0, refbus, numRows, numCols - refbus) = B.block(0, refbus + 1, numRows, numCols - refbus);
	B.conservativeResize(numRows,numCols);
	// Inverter matriz
	MatrixXd Ieye = MatrixXd::Identity(B.rows(), B.rows());
	//B = B.partialPivLu().solve(Ieye);
	B = B.fullPivLu().solve(Ieye);
	//B = B.inverse();

	// Incluir linha de 0's na posição da barra de referencia
	numRows = B.rows() + 1;
	numCols = B.cols();
	B.conservativeResize(numRows,numCols);
	MatrixXd temp = B.block(refbus, 0, numRows - refbus - 1, numCols);
	B.block(refbus + 1, 0, numRows - refbus - 1, numCols) = temp;
	for (int i = 0; i < B.cols(); i++)
		B(refbus, i) = 0;
	// Incluir coluna de 0's na posição da barra de referencia
	numRows = B.rows();
	numCols = B.cols() + 1;
	B.conservativeResize(numRows,numCols);
	temp = B.block(0, refbus, numRows, numCols - refbus - 1);
	B.block(0, refbus + 1, numRows, numCols - refbus - 1) = temp;
	for (int i = 0; i < B.rows(); i++)
		B(i, refbus) = 0;
	BetaS.MultiplicarPorMatriz(B);
	cout << "Calc. Beta com Sparse: " << double(clock() - tempo)/double(CLOCKS_PER_SEC) << endl;

	cout << "BetaS, linha 0 : " << endl;
	for (int l = 0; l < 10; l++)
		cout << BetaS.GetElemento(0 , l) << ";";
	cout << "Hrst, BetaS final" << endl;
	std::cin.ignore();

	//// Imprimir Beta em arquivo
	//ofstream * inFile;
	//inFile = new ofstream( "BetaS.txt", ios::out );
	//if ( inFile->is_open() )                                                            
	//{
	//	*inFile << std::scientific << setprecision(10);
	//	for (size_t i = 0; i < BetaS.GetNlin(); i++)
	//	{
	//		for (size_t j = 0; j < BetaS.GetNcol(); j++)
	//			*inFile << BetaS.GetElemento(i,j) << char(9);
	//		*inFile << endl;
	//	}
	//	*inFile << endl;
	//}
	//else
	//	cout << "Unable to open file";
	//inFile->close();
	//delete inFile;
}
void CSistema::DeterminarLimitesV()
{
	// para o loop paralelo as variaveis utilizadas devem ser declaradas dentro do loop, pois são individuais para cada laço
	parallel_for (size_t(0), hidreletricasVtr.size(), [&](size_t h)
	//for (size_t h = 0; h < hidreletricasVtr.size(); h++)
	{
		double dmontmax0 = 0;
		double dmax = 0;
		double dmontmax = 0;
		double vmin, vmax, vminbase, vmaxbase, vminbaseE, vmaxbaseE;
		hidreletricasVtr[h].SetVminSize(Tt1+n_cenarios*(Tt2-Tt1));
		hidreletricasVtr[h].SetVmaxSize(Tt1+n_cenarios*(Tt2-Tt1));
		// encontrar usinas a montante e definir o dmontmax (igual para qq periodo), dmontmax é a soma das defluencias máximas das usinas a montante
		for (int hm = 0; hm < hidreletricasVtr[h].GetNUsinaM(); hm++)
		{
			dmontmax0 += hidreletricasVtr[hidreletricasVtr[h].GetUsinaMont(hm) - 1].GetSmax();
			for (int j = 0; j < hidreletricasVtr[hidreletricasVtr[h].GetUsinaMont(hm) - 1].GetNGrupos(); j++)
				dmontmax0 += hidreletricasVtr[hidreletricasVtr[h].GetUsinaMont(hm) - 1].grupoVtr[j].GetQmax() * hidreletricasVtr[hidreletricasVtr[h].GetUsinaMont(hm) - 1].grupoVtr[j].GetNUnidades();
		}
		// definir dmax (igual para qq period), defluencia maxima da usina h
		dmax += hidreletricasVtr[h].GetSmax();		// apesar de ser possível, nao faz muito sentido, pois não é comum verter com reservatório baixo
		for (int j = 0; j < hidreletricasVtr[h].GetNGrupos(); j++)
			dmax += hidreletricasVtr[h].grupoVtr[j].GetQmax() * hidreletricasVtr[h].grupoVtr[j].GetNUnidades();
		//
		vminbase = hidreletricasVtr[h].GetV0();
		vmaxbase = hidreletricasVtr[h].GetV0();
		for (int t = 0; t < Tt1; t++)
		{
			// somar y (que varia com o periodo) ao dmotmax
			dmontmax = dmontmax0 + hidreletricasVtr[h].GetAfluencia(0, t);
			vmin = std::max(vminbase - dmax * 0.0036, hidreletricasVtr[h].GetVmin());		// vmin turbinar e verter tudo sem chegar nada de agua
			vmax = std::min(vmaxbase + dmontmax * 0.0036, hidreletricasVtr[h].GetVmax());	// vmax não turinar nem verter e chegar o máximo de afluencia total (incremental + usinas a montante)
			hidreletricasVtr[h].SetVolElem(t, vmin, vmax);

			// vbase atualizado para o próximo periodo
			vminbase = vmin;
			vmaxbase = vmax;
		}
		vminbaseE = vminbase;
		vmaxbaseE = vmaxbase;
		for (int cen  = 0; cen < n_cenarios;cen++)
		{
			vminbase = vminbaseE;
			vmaxbase = vmaxbaseE;
			for (int t = Tt1; t < Tt2; t++)
			{
				// somar y (que varia com o periodo) ao dmotmax
				dmontmax = dmontmax0 + hidreletricasVtr[h].GetAfluencia(0, t);
				vmin = std::max(vminbase - dmax * 0.0036, hidreletricasVtr[h].GetVmin());
				vmax = std::min(vmaxbase + dmontmax * 0.0036, hidreletricasVtr[h].GetVmax());
				hidreletricasVtr[h].SetVolElem(t + cen*(Tt2-Tt1), vmin, vmax);

				// vbase atualizado para o próximo periodo
				vminbase = vmin;
				vmaxbase = vmax;
			}
		}
	//}
	});
	// Conferir viabilidade do volume meta
	for (size_t h = 0; h < hidreletricasVtr.size(); h++)
		if (hidreletricasVtr[h].GetVmax(Tt1+n_cenarios*(Tt2-Tt1) - 1) < hidreletricasVtr[h].GetVMeta())
			cout << "Volume meta da hidrelétrica " << h << " é inviável!" << endl;
}

//void CSistema::AjustarPreconditioners()
//{
//
//	// arrumar essa função para cada tipo de decomposição!!!!
//
//	// L = [Lpt Lph Lv Ld flag_phmax*(Lphmax)]
//	int flag3 = int (flag_phmax);
//	//int deltaa = 0;
//	// Kind of Preconditioner
//	switch (flag_preconditioner)
//	{
//	case 1:		//probability of the node
//		{
//			for (int t = 0; t < Tt1; t++)
//			{
//				for (size_t i = 0; i < termeletricasVtr.size() + (3+flag3)*hidreletricasVtr.size() + (1 - flag3); i++)
//					precond.push_back(1);
//			}
//			for (int cen = 0; cen < n_cenarios; cen++)
//			{
//				for (int t = 0; t < Tt2 - Tt1; t++)
//				{
//					for (size_t i = 0; i < termeletricasVtr.size() + (3+flag3)*hidreletricasVtr.size() + (1 - flag3); i++)
//						precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen));
//				}
//			}
//			break;
//		}
//	case 2:		//square root of the probability of the node
//		{
//			for (int t = 0; t < Tt1; t++)
//			{
//				for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//					precond.push_back(1);
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//					precond.push_back(1);
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//					precond.push_back(1);
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//					precond.push_back(1);
//				if (flag3 == 1)
//				{
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//						precond.push_back(1);
//				}
//				else
//					precond.push_back(1);								// reserva
//			}
//			for (int cen = 0; cen < n_cenarios; cen++)
//			{
//				for (int t = 0; t < Tt2 - Tt1; t++)
//				{
//					for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)));
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)));
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)));
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)));
//					if (flag3 == 1)
//					{
//						for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//							precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)));
//					}
//					else
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)));	// reserva
//				}
//			}
//			break;
//		}
//	case 3:		//maximum value of the relaxed constraint
//		{
//			double phmax;
//			double qhmax;
//			for (int t = 0; t < Tt1 + n_cenarios*(Tt2 - Tt1); t++)
//			{
//				for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//				{
//					if (GetFlagTbinaryModel() == 1)
//						precond.push_back(1 / (termeletricasVtr[i].GetPmax() - termeletricasVtr[i].GetPmin()));
//					else
//						precond.push_back(1 / (termeletricasVtr[i].GetPmax() - 0));
//				}
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//				{
//					phmax = 0;
//					for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//						phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//					precond.push_back(1 / (phmax - 0));
//				}
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//					precond.push_back(1 / (hidreletricasVtr[r].GetVmax() - hidreletricasVtr[r].GetVmin()));
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//				{
//					qhmax = 0;
//					for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//						qhmax += (hidreletricasVtr[r].grupoVtr[j].GetQmax() * hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
//					precond.push_back(1 / (qhmax + hidreletricasVtr[r].GetSmax() - 0));
//				}
//				if (flag3 == 1)
//				{
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//					{
//						phmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(1 / (phmax - 0));
//					}
//				}
//				else
//				{
//					phmax = 0;
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//reserva
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//					precond.push_back(1 / (phmax - res_gir[t]));
//				}
//			}
//			break;
//		}
//	case 4:		//combinaion between 1 and 3
//		{
//			double phmax;
//			double qhmax;
//			for (int t = 0; t < Tt1; t++)
//			{
//				for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//				{
//					if (GetFlagTbinaryModel() == 1)
//						precond.push_back(1 / (termeletricasVtr[i].GetPmax() - termeletricasVtr[i].GetPmin()));
//					else
//						precond.push_back(1 / (termeletricasVtr[i].GetPmax() - 0));
//				}
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//				{
//					phmax = 0;
//					for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//						phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//					precond.push_back(1 / (phmax - 0));
//				}
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//					precond.push_back(1 / (hidreletricasVtr[r].GetVmax() - hidreletricasVtr[r].GetVmin()));
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//				{
//					qhmax = 0;
//					for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//						qhmax += (hidreletricasVtr[r].grupoVtr[j].GetQmax() * hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
//					precond.push_back(1 / (qhmax + hidreletricasVtr[r].GetSmax() - 0));
//				}
//				if (flag3 == 1)
//				{
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//					{
//						phmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(1 / (phmax - 0));
//					}
//				}
//				else
//				{
//					phmax = 0;
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//reserva
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//					precond.push_back(1 / (phmax - res_gir[t]));
//				}
//			}
//			for (int cen = 0; cen < n_cenarios; cen++)
//			{
//				for (int t = 0; t < Tt2 - Tt1; t++)
//				{
//					for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//					{
//						if (GetFlagTbinaryModel() == 1)
//							precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (termeletricasVtr[i].GetPmax() - termeletricasVtr[i].GetPmin()));
//						else
//							precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (termeletricasVtr[i].GetPmax() - 0));
//					}
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//					{
//						phmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (phmax - 0));
//					}
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//						precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (hidreletricasVtr[r].GetVmax() - hidreletricasVtr[r].GetVmin()));
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//					{
//						qhmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							qhmax += (hidreletricasVtr[r].grupoVtr[j].GetQmax() * hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
//						precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (qhmax + hidreletricasVtr[r].GetSmax() - 0));
//					}
//					if (flag3 == 1)
//					{
//						for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//						{
//							phmax = 0;
//							for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//								phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//							precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (phmax - 0));
//						}
//					}
//					else
//					{
//						phmax = 0;
//						for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//reserva
//							for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//								phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(hidreletricasVtr[0].GetProbAfluencia(cen) / (phmax - res_gir[t]));
//					}
//				}
//			}
//			break;
//		}
//	case 5:		//combination between 2 and 3
//		{
//			double phmax;
//			double qhmax;
//			for (int t = 0; t < Tt1; t++)
//			{
//				for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//				{
//					if (GetFlagTbinaryModel() == 1)
//						precond.push_back(1 / (termeletricasVtr[i].GetPmax() - termeletricasVtr[i].GetPmin()));
//					else
//						precond.push_back(1 / (termeletricasVtr[i].GetPmax() - 0));
//				}
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//				{
//					phmax = 0;
//					for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//						phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//					precond.push_back(1 / (phmax - 0));
//				}
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//					precond.push_back(1 / (hidreletricasVtr[r].GetVmax() - hidreletricasVtr[r].GetVmin()));
//				for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//				{
//					qhmax = 0;
//					for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//						qhmax += (hidreletricasVtr[r].grupoVtr[j].GetQmax() * hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
//					precond.push_back(1 / (qhmax + hidreletricasVtr[r].GetSmax() - 0));
//				}
//				if (flag3 == 1)
//				{
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//					{
//						phmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(1 / (phmax - 0));
//					}
//				}
//				else
//				{
//					phmax = 0;
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//reserva
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//					precond.push_back(1 / (phmax - res_gir[t]));
//				}
//			}
//			for (int cen = 0; cen < n_cenarios; cen++)
//			{
//				for (int t = 0; t < Tt2 - Tt1; t++)
//				{
//					for (size_t i = 0; i < termeletricasVtr.size(); i++)	//pt
//					{
//						if (GetFlagTbinaryModel() == 1)
//							precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (termeletricasVtr[i].GetPmax() - termeletricasVtr[i].GetPmin()));
//						else
//							precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (termeletricasVtr[i].GetPmax() - 0));
//					}
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//ph
//					{
//						phmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (phmax - 0));
//					}
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//v
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (hidreletricasVtr[r].GetVmax() - hidreletricasVtr[r].GetVmin()));
//					for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//d
//					{
//						qhmax = 0;
//						for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//							qhmax += (hidreletricasVtr[r].grupoVtr[j].GetQmax() * hidreletricasVtr[r].grupoVtr[j].GetNUnidades());
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (qhmax + hidreletricasVtr[r].GetSmax() - 0));
//					}
//					if (flag3 == 1)
//					{
//						for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//phmax
//						{
//							phmax = 0;
//							for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//								phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//							precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (phmax - 0));
//						}
//					}
//					else
//					{
//						phmax = 0;
//						for (size_t r = 0; r < hidreletricasVtr.size(); r++)	//reserva
//							for (int j = 0; j < hidreletricasVtr[r].GetNGrupos(); j++)
//								phmax += hidreletricasVtr[r].grupoVtr[j].GetPmax()*hidreletricasVtr[r].grupoVtr[j].GetNUnidades();
//						precond.push_back(sqrt(hidreletricasVtr[0].GetProbAfluencia(cen)) / (phmax - res_gir[t]));
//					}
//				}
//			}
//			break;
//		}
//	default:	//none of them
//		{
//			for (int t = 0; t < Tt1 + n_cenarios*(Tt2 - Tt1); t++)
//			{
//				for (size_t i = 0; i < termeletricasVtr.size() + (3+flag3)*hidreletricasVtr.size() + (1 - flag3); i++)
//					precond.push_back(1);
//			}
//			break;
//		}
//	}
//}


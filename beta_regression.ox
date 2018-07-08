/*******************************************************************************
 * PROGRAMA: beta_regression.ox                                                *
 *                                                                             *
 * USO:      Esse script faz simulações de Monte Carlo para estimar parâmetros *
 *           da Regressão Beta proposta por Ferrari e Cribari-Neto (2004) com  *
 *           parâmetro de precisão constante.                                  *
 *                                                                             *
 *           A estimação é feita por Estimação Pontual, Estimação Intervalar e *
 *           por Testes de Hipóteses, esta última fazendo simulação de tamanho *
 *           (impondo H0 verdadeira) e simulação de poder (impondo H0 falsa).  *
 *           Para o menor dos tamanhos de amostra (25 observações), é realiza- *
 *           da uma correção de viés (BC1) da Estimação Pontual.               *
 *                                                                             *
 *           A geração de ocorrências de uma Distribuição Beta pode ser feita  *
 *           usando a função ranbeta do Ox, ou usando o método de rejeição ca- *
 *           nônico com uma distribuição uniforme como distribuição instrumen- *
 *           tal e mantendo reproducibidade dos dados gerados em outras lin-   *
 *           guagens.                                                          *
 *                                                                             *
 * AUTOR:    Leonardo Melo                                                     *
 *                                                                             *
 * VERSÃO:   1.0                                                               *
 *                                                                             *
 * DATA:     25 de Maio de 2018                                                *
 *                                                                             *
 * NOTAS:    Ao executar esse programa, a diretiva de pré-processamento condi- *
 *           cional DONT_PANIC em "beta_fun.ox" não será executada, ocultando  *
 *           a impressão de "beta_fun.ox" e deixando apenas a impressão desse  *
 *           script.                                                           *
 *                                                                             *
 ******************************************************************************/

// Inclusão de arquivos cabeçalhos, Definições e Importação de arquivos //

#import  <maximize>
#define  DONT_PANIC
#include "beta_fun.ox"		

// Escolha de efetuar simulações adicionais para gerar dados do Poder dos Testes //
static decl s_bPlots = 0; 

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * s_bPlots == 1 :  Calcula e armazena em dados.txt o poder dos 3 testes do  *
   *                  programa para diferentes valores de beta1 na em H0.      *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Protótipo de Função //
montecarlo(const N, const R, const B, const rbeta_usr, const SEED, const verbose, const file);

// Main //

main()
{
	// Iniciando o cronômetro
	decl dtime = timer();
	decl R = 1000, B = 250;
	decl k, file = fopen("dados.txt", "w");

	// Imprimindo informações na tela
	println("\t");												      	  
	println("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]");
	println("[]                                                                          []");
	println("[]                     RESULTADOS PARA A REGRESSÃO BETA                     []");
	println("[]                      -  SIMULAÇÃO DE MONTE CARLO -                       []");
	println("[]                                                                          []");
	println("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]");

	decl SEED = {1802558351, -1609231228};	// semente 3112 do R 
	decl vfalse = <0.6; 0.7; 0.8; 0.85; 0.9; 1.0; 1.1>;
	
	for	(k = 0; k < 7; k++)
	{
		if(s_bPlots != 1)              // Se não quiser gerar o poder para diferentes
			dfalse = vfalse[1];    // valores de beta1, então gere só para 1 deles
		else {                         // e faça a simulação de poder apenas alterando 
			dfalse = vfalse[k];    // o tamanho amostral
		}

		montecarlo(25, R,B,SEED,1,0,file);
		montecarlo(50, R,0,SEED,1,0,file);
		montecarlo(100,R,0,SEED,1,0,file);
		montecarlo(150,R,0,SEED,1,0,file);
		montecarlo(200,R,0,SEED,1,0,file);
		montecarlo(250,R,0,SEED,1,0,file);

		if(s_bPlots != 1)				 
			k = 6;
	}

	// Imprimir data e tempo de execução
	println("\n");
	println( " Data: ", date());
	println( " Tempo de execução: ", timespan(dtime) , "\n\n");
	println( " Fim. \n");

	println("==============================================================================");
	println("\n\n");

    fclose(file);  
}	

// Função Monte Carlo //
montecarlo(const N, const R, const B, const SEED, const rbeta_usr, const verbose, const file)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *  Uso                                                                    *
	 *                                                                         *
	 *    N             Tamanho amostral.                                      *
	 *    R             Número de Réplicas de Monte Carlo.                     *
	 *    B             Número de Réplicas de Bootstrap.                       *
	 *    F_H0false     Valor de H0 quando falsa.                              *
	 *    SEED          Semente para garantir reprodutibilidade com R.         *
	 *                                                                         *
	 *    rbeta_usr     Boolean (default é TRUE). Geração da distribuição Beta *
	 *                  pelo método da Rejeição usando a distribuição Uniforme *
	 *                  como distribuição instrumental. Aviso: menos eficiente!*
	 *                  Se FALSE, a função ranbeta é utilizada.                *
	 *                                                                         *
	 *    verbose       Boolean (default é FALSE). Se TRUE, visualiza o anda-  *
	 *                  mento do laço com a mensagem: Réplica 'j' de 'R'.      *
	 *                                                                         *
	 *    file          Arquivo que armazenará os dados para o plot do Poder   *
	 *                  do teste.                                              *
	 *                                                                         *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	// Declaração de parâmetros da simulação e de outras variáveis
	decl vparams = <0.85; -0.93; -1.22; 40.5>;  // vetor de parâmetros verdadeiros
	decl dc = 0.1;                              // Incremento na razão de chances
	decl vparamsODDS = exp(dc*vparams[:2]);     // Razão de chances dos parâmetros Beta
	decl vestimators = <0;0;0;1>;               // chute inicial
	decl vest_bs = <0;0;0;1>;
	
	decl j, i, ir;
	decl fail, failB, failT, failF;
	fail = failB = failT = failF = 0;
		// Contador de falhas de convergência da estimação:
		// 1) da regressão                    2) da reg. nas réplicas de Bootstrap
		// 3) do teste quando H0 verdadeira   4) do teste quando H0 falso
	
	decl dval, dval0;
          // valor da LogVeros. com os Estimadores e da Restrita
	
	// Variáveis da Estimação Pontual //
	decl mestimators, mestimators_bs;
	mestimators = mestimators_bs = zeros(R, 4);	
	  // matrizes de estimativas Monte Carlo, de Erros Padrão e de Bootstrap
	decl mbootstrap = zeros(B, 4);
	
	// Variáveis da Estimação Intervalar //
	decl vcIC99, vcIC95, vcIC90;
	vcIC99 = vcIC95 = vcIC90 = <0;0;0;0>;
          // contadores de que pertencem ao IC

	decl vcIC99ODDS, vcIC95ODDS, vcIC90ODDS;
	vcIC99ODDS = vcIC95ODDS = vcIC90ODDS = <0;0;0>;					
          // contadores de que pertencem ao IC de Razão de Chances			 
	
	// Variáveis do Teste de Hipóteses //											   
	decl vrest, mK_1;
	decl vdados = zeros(1,5);
	
	decl dVC10_2 = quanchi(0.90, 2);	// Tamanho												
	decl dVC5_2 = quanchi(0.95, 2);		// val. críticos a 10%, 5% e 1% da Qui-Quadrado
	decl dVC1_2 = quanchi(0.99, 2);		// com 2 graus de liberdade (2 restrições no TH)

	decl dVC10_1 = quanchi(0.90, 1);	// Poder												
	decl dVC5_1 = quanchi(0.95, 1);		// val. críticos a 10%, 5% e 1% da Qui-Quadrado
	decl dVC1_1 = quanchi(0.99, 1);		// com 1 graus de liberdade (1 restrição no TH)

	decl vLR1, vSc1, vWa1, vLR2, vSc2, vWa2;
	vLR1 = vSc1 = vWa1 = vLR2 = vSc2 = vWa2 = zeros(R, 1);
          // Vetores estatísticas de teste (1 para poder e 2 para tamanho)
	
	// Variáveis Independentes
	ranseed("GM");                // Gerador George Marsaglia
	ranseed(SEED);							
	g_vx2 = rann(N, 1);					  

	ranseed(SEED);                // Fixar semente (Uniforme) 
	g_vx3 = ranu(N, 1); 
	g_mX = 1~g_vx2~g_vx3;         // Matrix X

	
	// \eta e \mu (usando a Logit como funcao de ligação)
	decl veta = g_mX * vparams[:2];
	decl vmu = (exp(veta)) ./ (1 + exp(veta));


	ranseed(SEED);                // Fixar semente (Beta) 
	MaxControl(50, -1);           // Limitando o número de iterações

	// Declaração de Variáveis que serão usadas no Laço
	decl mtemp, mFishInv, vSE;
	decl vIC99min, vIC99max, vIC95min;
	decl vIC95max, vIC90min, vIC90max;          // IC dos parâmetros
	decl vIC99minODDS, vIC99maxODDS, vIC95minODDS;
	decl vIC95maxODDS, vIC90minODDS, vIC90maxODDS;     // IC da Razão de Chances
	decl vU, mFishInv_rest;                     // Para a estatística do teste Escore			
	decl veta_hat, vmu_hat, dphi_hat;           // Para o Bootstrap Paramétrico
	
	// Laço Monte Carlo
	for (j = 0; j < R; j++)
	{
		decl iSEED = {j+1,j+1};
		ranseed(iSEED);
		// Gerando N ocorrências da distribuição Beta 
		if (rbeta_usr)
			g_vy = vfBeta_AccRej(N, vmu .* vparams[3], (1 - vmu) .* vparams[3]);
		else
		 	g_vy = ranbeta(N, 1, vmu .* vparams[3], (1 - vmu) .* vparams[3]);

		// ESQUECENDO OS VALORES DOS PARÂMETROS //

		// Encontrando os estimadores por Máxima Verossimilhança (MV)
		ir = MaxBFGS(dfRegBeta, &vestimators, &dval, 0, FALSE); 	

		if (ir == MAX_CONV || ir == MAX_WEAK_CONV)
		{
			/* * * ESTIMAÇÃO PONTUAL * * */
		
			// Armazenando o estimador da MV na matriz de Monte Carlo 
			mestimators[j][] = vestimators';

			/* * * ESTIMAÇÃO INTERVALAR * * */
			decl mKest = mfKmatrix(vestimators);       // Informação de Fisher
			mFishInv = invertsym(mKest);               // Inversa da Informação de Fisher
			vSE = sqrt(diagonal(mFishInv));            // Erro padrão
			
			// Intervalos de Confiança dos Parâmetros
			vIC99min = (vestimators - quann(1 - 0.01/2) * vSE');					
			vIC99max = (vestimators + quann(1 - 0.01/2) * vSE');			
			
			vIC95min = (vestimators - quann(1 - 0.05/2) * vSE');
			vIC95max = (vestimators + quann(1 - 0.05/2) * vSE');
			
			vIC90min = (vestimators - quann(1 - 0.10/2) * vSE');
			vIC90max = (vestimators + quann(1 - 0.10/2) * vSE');
			
			// IC da Razão de Chances (apenas Beta_1, Beta_2 e Beta_3)
			vIC99minODDS = exp(dc*(vestimators[:2] - quann(1 - 0.01/2) * (vSE[:2])'));	 
			vIC99maxODDS = exp(dc*(vestimators[:2] + quann(1 - 0.01/2) * (vSE[:2])'));
			
			vIC95minODDS = exp(dc*(vestimators[:2] - quann(1 - 0.05/2) * (vSE[:2])'));
			vIC95maxODDS = exp(dc*(vestimators[:2] + quann(1 - 0.05/2) * (vSE[:2])'));
			
			vIC90minODDS = exp(dc*(vestimators[:2] - quann(1 - 0.10/2) * (vSE[:2])'));
			vIC90maxODDS = exp(dc*(vestimators[:2] + quann(1 - 0.10/2) * (vSE[:2])'));

			// Incrementando contadores dos IC's
			vcIC99 = vcIC99 + (vparams .> vIC99min .&& vparams .< vIC99max); 					
			vcIC95 = vcIC95 + (vparams .> vIC95min .&& vparams .< vIC95max); 					
			vcIC90 = vcIC90 + (vparams .> vIC90min .&& vparams .< vIC90max); 					

			vcIC99ODDS = vcIC99ODDS +
                                    (vparamsODDS .> vIC99minODDS .&& vparamsODDS .< vIC99maxODDS);
			vcIC95ODDS = vcIC95ODDS +
                                    (vparamsODDS .> vIC95minODDS .&& vparamsODDS .< vIC95maxODDS);
			vcIC90ODDS = vcIC90ODDS +
                                    (vparamsODDS .> vIC90minODDS .&& vparamsODDS .< vIC90maxODDS);
   																					   
			/* * * TESTE DE HIPÓTESES * * */
			
			// ======================= //
			//    Com H0 verdadeira    //
			// ======================= //
			
			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 *   As Hipóteses do Teste são (2 restrições):                           *
			 *                                                                       *
			 *      H0:      beta1 =  0.85      vs.     H1:       beta1 !=  0.85     *
			 *          (e)  beta2 = -0.93      vs.         (ou)  beta2 != -0.93     *
			 *                                                                       *
			 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */	
			 
			// Calculando o valor da Log-Verossimilhança Restrita
			vrest = <0;1>;			
			MaxBFGS(dfRegBetaTRUE, &vrest, &dval0, 0, FALSE);	
	
			// Armazenando a estatística de teste da Razão de Verossimilhança
			vLR2[j] = 2*(dval - dval0);
			

			// Armazenando a estatística de teste Escore
			vrest = (0.85)|(-0.93)|vrest;
			Num1Derivative(dfRegBeta, vrest, &vU);
			mtemp = mfKmatrix(vrest);
			mFishInv_rest = invertsym(mtemp);

			vSc2[j] = (vU[:1])' * mFishInv_rest[:1][:1] * (vU[:1]);


			// Armazenando a estatística de teste de Wald
			mK_1 = invert(mFishInv[:1][:1]);
			vWa2[j] = (vestimators[:1] - vrest[:1])' * mK_1 * (vestimators[:1] - vrest[:1]);

			
			// ==================== //
			//     Com H0 falsa     //
			// ==================== //
	
			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 *   A Hipótese do Teste é (1 restrição):                                *
			 *                                                                       *
			 *      H0:      beta1 = dfalse     vs.     H1:       beta1 != dfalse    *
			 *                                                                       *
			 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
		 
			// Armazenando a estatística de teste da Razão de Verossimilhança
			vrest = <0;0;1>;			
			MaxBFGS(dfRegBetaFALSE, &vrest, &dval0, 0, FALSE);	
			
			vLR1[j] = 2*(dval - dval0);
			

			// Armazenando a estatística de teste Escore
			vrest = dfalse|vrest;
			Num1Derivative(dfRegBeta, vrest, &vU);
			mtemp = mfKmatrix(vrest);
			mFishInv_rest = invertsym(mtemp);
			
			vSc1[j] = (vU[0])' * mFishInv_rest[0][0] * (vU[0]);

			
			// Armazenando a estatística de teste de Wald
			mK_1 = invert(mFishInv[0][0]);
			vWa1[j] = (vestimators[0] - vrest[0])' * mK_1 * (vestimators[0] - vrest[0]);


			/* * * BOOTSTRAP PARAMÉTRICO * * */
		  	if (B) {				 
				veta_hat = g_mX * vestimators[:2];
				vmu_hat = (exp(veta_hat)) ./ (1 + exp(veta_hat));
				dphi_hat = vestimators[3];
				ranseed(SEED);				
				for (i = 0; i < B; i++)
				{
					if (rbeta_usr)
						s_y_bs = vfBeta_AccRej(N, vmu_hat.*dphi_hat, (1-vmu_hat).*dphi_hat);
					else
					 	s_y_bs = ranbeta(N, 1, vmu_hat.*dphi_hat, (1-vmu_hat).*dphi_hat);
						
					ir = MaxBFGS(dfRegBetaBS, &vest_bs, &dval, 0, FALSE);
	
					if (ir == MAX_CONV || ir == MAX_WEAK_CONV)
						mbootstrap[i][] = vest_bs';
					else {
						println("Falha de Convergência na Réplica de Bootstrap");
						++failB;
						--i; // Contando falhas de convergência
					}
				}
	
				// Calculando o estimador corrigido de Bootstrap (BC1)
				mestimators_bs[j][] = ((2 * vestimators - meanc(mbootstrap)'))';
		   	}

			if (verbose)
				println("Réplica: ", (j+1), " de ", R);		
		}
		
		else { 
			++fail;
			--j;
		} // Contando falhas de convergência
	}

	// Arrumando dados para impressão
	
	// Média dos Estimadores de MV
	decl vest_mean = meanc(mestimators)';
	decl vest_mean_bs = meanc(mestimators_bs)';

	// Variância dos Estimadores de MV
	decl vest_var = varc(mestimators)';
	decl vest_var_bs = varc(mestimators_bs)';
		
	// Viés e viés relativo dos EMV
	decl vest_bias = vest_mean - vparams;		
	decl vest_bias_bs = vest_mean_bs - vparams;		
	decl vest_pbias = 100*vest_bias./vparams;
	decl vest_pbias_bs = 100*vest_bias_bs./vparams;
	
	// EQM dos EMV
  	decl vEQM = vest_bias .^2 + vest_var;
  	decl vEQM_bs = vest_bias_bs .^2 + vest_var_bs;

	// Salvando dados para Plot (se s_bPlots == 1)
	if(s_bPlots)
	{
		decl u;
		vdados[0] = N;
		vdados[1] = dfalse;
		vdados[2] = 100 * (sumc(vLR1 .> dVC5_1) / R);
		vdados[3] = 100 * (sumc(vSc1 .> dVC5_1) / R);
		vdados[4] = 100 * (sumc(vWa1 .> dVC5_1) / R);
	
		if (isfile(file))
	    {
	        for (u = 0; u < 5; u++)
				fprint(file, vdados[u], ";");
			fprintln(file, "");
	    }
	}

	// Imprimindo Informações na tela
	println("\n==============================================================================\n\n");
	println("\t TAMANHO DA AMOSTRA: ", "                ", "% d", N, "\n\n");
	println("\t Método:                              ", "BFGS");
	println("\t Função de Geração Beta:              ", rbeta_usr?"Rejeição\n":"ranbeta\n");
	println("\t Algoritmo de Geração Uniforme:       ", ranseed(""));
	println("\t Sementes no Ox:                     ", "% d", SEED[0], " e ", "% d", SEED[1]);
	println("\t Número de Réplicas de Monte Carlo:  ", "% d", R);
	println("\t Número de Réplicas de Bootstrap:    ", "% d", B, "\n\n\n");
	println(" ESTIMAÇÃO PONTUAL");
	println("==============================================================================");
	print("% 13.6f", "%c", {"    Parâmetro", "    Média EMVs", "Viés", "    Viés Rel %"},
			vparams~vest_mean~vest_bias~vest_pbias);
	println("% 13.6f", "%c", {"    Variância", "EQM  "},
			vest_var~vEQM);
	if(B) {
		println("  -------------   Correção de Bootstrap   -------------");
		print("% 13.6f", "%c", {"    Parâmetro", "    Média EMVs", "Viés", "    Viés Rel %"},
				vparams~vest_mean_bs~vest_bias_bs~vest_pbias_bs);
		println("% 13.6f", "%c", {"    Variância", "EQM  "},
				vest_var_bs~vEQM_bs);
	}
	println("\t\n");
	println(" ESTIMAÇÃO INTERVALAR");
	println("==============================================================================\n\n");
	println("   ------- Parâmetros Beta1, Beta2, Beta3 e Phi  -------");
	println("% 13.2f", "%c", {"    Parâmetro", "    [ 90 % ]", "  [ 95 % ]", "    [ 99 % ]"},
			vparams~(100 * (vcIC90) / R)~(100 * (vcIC95) / R)~(100 * (vcIC99) / R));
	println("\t\n");
	println("   ----- Razão de Chances para Beta1, Beta2 e Beta3 ----");
	println("% 13.2f", "%c", {"   Verdadeira", "    [ 90 % ]", "  [ 95 % ]", "    [ 99 % ]"},
			vparamsODDS~(100 * (vcIC90ODDS) / R)~(100 * (vcIC95ODDS) / R)~(100 * (vcIC99ODDS) / R));
	println("\t\n");
	println(" TESTE DE HIPÓTESES");
	println("==============================================================================\n\n");
	println("   ---------------  Simulação de Tamanho  ---------------","\n\n");
	println("        Valor crítico 10%:                ","%10.6f", dVC10_2);
	println("        Valor crítico 5%:                 ","%10.6f", dVC5_2);
	println("        Valor crítico 1%:                 ","%10.6f", dVC1_2, "\n\n");
	println("      ---  Probabilidade de cometer o erro tipo I  ---   ", "\n\n");
	println("        Razão de Verossimilhança a 10%: ","%10.2f", double (100 * (sumc(vLR2 .> dVC10_2) / R)), " %");
	println("        Razão de Verossimilhança a  5%: ","%10.2f", double (100 * (sumc(vLR2 .> dVC5_2) / R)), " %");
	println("        Razão de Verossimilhança a  1%: ","%10.2f", double (100 * (sumc(vLR2 .> dVC1_2) / R)), " %\n\n");
	println("        Teste Escore a 10%:             ","%10.2f", double (100 * (sumc(vSc2 .> dVC10_2) / R)), " %");
	println("        Teste Escore a  5%:             ","%10.2f", double (100 * (sumc(vSc2 .> dVC5_2) / R)), " %");
	println("        Teste Escore a  1%:             ","%10.2f", double (100 * (sumc(vSc2 .> dVC1_2) / R)), " %\n\n");
	println("        Teste Wald a 10%:               ","%10.2f", double (100 * (sumc(vWa2 .> dVC10_2) / R)), " %");
	println("        Teste Wald a  5%:               ","%10.2f", double (100 * (sumc(vWa2 .> dVC5_2) / R)), " %");
	println("        Teste Wald a  1%:               ","%10.2f", double (100 * (sumc(vWa2 .> dVC1_2) / R)), " %\n\n\n");
	println("   ----------------  Simulação de Poder  ----------------","\n\n");
	println("        Valor crítico 10%:                ","%10.6f", dVC10_1);
	println("        Valor crítico 5%:                 ","%10.6f", dVC5_1);
	println("        Valor crítico 1%:                 ","%10.6f", dVC1_1, "\n\n");
	println("      ---    Probabilidade de rejeitar H0 falsa    ---   ", "\n\n");
	println("           ---  Hipótese Nula:   Beta1 = ", dfalse, "   ---  \n\n");
	println("        Razão de Verossimilhança a 10%: ","%10.2f", double (100 * (sumc(vLR1 .> dVC10_1) / R)), " %");
	println("        Razão de Verossimilhança a  5%: ","%10.2f", double (100 * (sumc(vLR1 .> dVC5_1) / R)), " %");
	println("        Razão de Verossimilhança a  1%: ","%10.2f", double (100 * (sumc(vLR1 .> dVC1_1) / R)), " %\n\n");
	println("        Teste Escore a 10%:             ","%10.2f", double (100 * (sumc(vSc1 .> dVC10_1) / R)), " %");
	println("        Teste Escore a  5%:             ","%10.2f", double (100 * (sumc(vSc1 .> dVC5_1) / R)), " %");
	println("        Teste Escore a  1%:             ","%10.2f", double (100 * (sumc(vSc1 .> dVC1_1) / R)), " %\n\n");
	println("        Teste Wald a 10%:               ","%10.2f", double (100 * (sumc(vWa1 .> dVC10_1) / R)), " %");
	println("        Teste Wald a  5%:               ","%10.2f", double (100 * (sumc(vWa1 .> dVC5_1) / R)), " %");
	println("        Teste Wald a  1%:               ","%10.2f", double (100 * (sumc(vWa1 .> dVC1_1) / R)), " %\n\n");
	println("==============================================================================\n\n");
	println(" Num de falhas de Convergência da Maximização:         ", fail, " de ", R);
	if(B)
		println(" Num de falhas de Convergência do Bootstrap:           ", failB, " de ", (R*B));
	println(" Num de falhas de Convergência do Teste H0 Verdadeira: ", failT, " de ", R);
	println(" Num de falhas de Convergência do Teste H0 Falsa:      ", failF, " de ", R, "\n\n");
	println("==============================================================================\n\n");

}

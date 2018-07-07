/*******************************************************************************
 * PROGRAMA: beta_fun.ox                                                       *
 *                                                                             *
 * USO:      Esse script contém declarações de funções importantes ao programa *
 *           principal "beta_regression.ox", entre elas: funções de log-veros- *
 *           similhança da Regressao Beta (com função de ligação Logística) –  *
 *           inclusive para testes de hipóteses com H0 verdadeira e H0 falsa e *
 *           bootstrap –, os respectivos gradientes, uma funcão de cálculo do  *
 *           valor da log-verossimilhança, uma função para geração de ocorrên- *
 *           cias da Distribuição Beta pelo método da rejeição, usando como    *
 *           distribuição instrumental uma Uniforme, e uma função para calcu-  *
 *           lar analiticamente a inversa da informação de Fisher.             *
 *                                                                             *
 * AUTOR:    Leonardo Melo                                                     *
 *                                                                             *
 * VERSÃO:   1.0                                                               *
 *                                                                             *
 * DATA:     25 de Maio de 2018                                                *
 *                                                                             *
 ******************************************************************************/

// Inclusão de arquivos cabeçalhos //

#include <oxstd.h>
#include <oxprob.h>

// Declaração de Variáveis Globais com Escopo de Arquivo e Ligação Externa //

decl g_vy;
decl s_y_bs;
decl g_mX;
decl g_vx2;
decl g_vx3;
decl dfalse;                    // Valor de beta1 sob H0 quando H0 é falsa
	 
// ==================================================================== //
//   Declaração da Função de Log-Verossimilhança que entra em MaxBFGS   //
// ==================================================================== //

dfRegBeta(const vP, const adFunc, const avScore, const amHess)
{
	// Declarando variáveis	da Log-Verossimilhança //
	decl k = columns(g_mX) + 1;
	decl dphi = vP[k-1];
	decl veta = g_mX * vP[0:(k-2)];
	
	// Usando a Logit como função de ligação //
	decl vmu = (exp(veta)) ./ (1 + exp(veta)); 

	// Criando o vetor de Log-Verossimilhanças para cada observação t //
	decl vLogLike_t = loggamma(dphi) - loggamma(vmu .* dphi) 
	                - loggamma((1-vmu) .* dphi) 
			+ ((vmu .* dphi) - 1).*log(g_vy) 
	                + (((1-vmu).*dphi) - 1).*log(1-g_vy); 
					
	// Somando o vetor para obter a Log-Verossimilhança //
	adFunc[0] = sumc(vLogLike_t);

	// Declarando variáveis do gradiente //
	if(avScore)
	{
		decl vystar = log( (g_vy) ./ (1-g_vy) );
		// y* = log{y/(1-y)}
		decl vmustar = polygamma(vmu .* dphi, 0)
					 - polygamma((1-vmu).*dphi, 0);  
		// mu* = psi(mu*phi) - psi((1-mu)*phi)
	 	decl mT = diag( exp(veta) ./ (1 + exp(veta)) .^2 );  
	 	// T = diag{1/g'(mu_1), ..., 1/g'(mu_n)}
		
		// Derivadas com relação ao vetor de parâmetros beta (U_beta)
		(avScore[0])[0:(k-2)] = dphi * g_mX' * mT * (vystar - vmustar);
		
		// Derivada com relação a phi (U_phi)
		decl dU_phi = polygamma(dphi, 0) - vmu .* polygamma(vmu .* dphi, 0)
						  - (1-vmu) .* polygamma( (1-vmu).*dphi, 0)
						  + vmu .* log(g_vy) + (1-vmu).*log(1-g_vy);
		(avScore[0])[(k-1)] = sumc(dU_phi); 		 	
	}	     	  

	// Retorno da função //

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ){
		return 0;  // indica falha
	} else{
		return 1;  // indica sucesso
	}
}


// ==================================================================== //
//   Declaração da Função de Log-Verossimilhança que entra em MaxBFGS   //
//                             no Bootstrap                             //
// ==================================================================== //

dfRegBetaBS(const vP, const adFunc, const avScore, const amHess)
{
	// Declarando variáveis	da Log-Verossimilhança //
	decl k = columns(g_mX) + 1;
	decl dphi = vP[k-1];
	decl veta = g_mX * vP[0:(k-2)];
	
	// Usando a Logit como função de ligação //
	decl vmu = (exp(veta)) ./ (1 + exp(veta)); 
 	
	// Criando o vetor de Log-Verossimilhanças para cada observação t //
	decl vLogLike_t = loggamma(dphi) - loggamma(vmu .* dphi) 
	                - loggamma((1-vmu) .* dphi)
	                + ((vmu .* dphi) - 1).*log(s_y_bs) 
	                + (((1-vmu).*dphi) - 1).*log(1-s_y_bs);
 					
	// Somando o vetor para obter a Log-Verossimilhança //
	adFunc[0] = sumc(vLogLike_t);

	// Declarando variáveis do gradiente //
	if(avScore)
	{
		decl vystar = log( (s_y_bs) ./ (1-s_y_bs) );       
		// y* = log{y/(1-y)}
		decl vmustar = polygamma(vmu .* dphi, 0)
					 - polygamma((1-vmu).*dphi, 0);         
		// mu* = psi(mu*phi) - psi((1-mu)*phi)
	 	decl mT = diag( exp(veta) ./ (1 + exp(veta)) .^2 );  
	 	// T = diag{1/g'(mu_1), ..., 1/g'(mu_n)}
		
		// Derivadas com relação ao vetor de parâmetros beta (U_beta)
		(avScore[0])[0:(k-2)] = dphi * g_mX' * mT * (vystar - vmustar);

		// Derivada com relação a phi (U_phi)
		decl dU_phi = polygamma(dphi, 0) - vmu .* polygamma(vmu .* dphi, 0)
						  - (1-vmu) .* polygamma( (1-vmu).*dphi, 0)
						  + vmu .* log(s_y_bs) + (1-vmu).*log(1-s_y_bs);
		(avScore[0])[(k-1)] = sumc(dU_phi); 		 	
	}	     	  

	// Retorno da função //

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ){
		return 0;  // indica falha
	} else{
		return 1;  // indica sucesso
	}
}	


// ==================================================================== //
//   Declaração da Função de Log-Verossimilhança que entra em MaxBFGS   //
//             no Teste de Hipóteses quando H0 é VERDADEIRA             //
// ==================================================================== //

dfRegBetaTRUE(const vP, const adFunc, const avScore, const amHess)
{
	// Declarando variáveis	da Log-Verossimilhança //
	decl dphi = vP[1];
	decl veta = g_mX * ((0.85)|(-0.93)|vP[0]);

	// Usando a Logit como função de ligação //
	decl vmu = (exp(veta)) ./ (1 + exp(veta)); 

	// Criando o vetor de Log-Verossimilhanças para cada observação t //
	decl vLogLike_t = loggamma(dphi) - loggamma(vmu .* dphi) 
	                - loggamma((1-vmu) .* dphi) 
			+ ((vmu .* dphi) - 1).*log(g_vy) 
	                + (((1-vmu).*dphi) - 1).*log(1-g_vy);
					
	// Somando o vetor para obter a Log-Verossimilhança //
	adFunc[0] = sumc(vLogLike_t);
	
	// Declarando variáveis do gradiente //
	if(avScore)
	{
		decl vystar = log( (g_vy) ./ (1-g_vy) );
		// y* = log{y/(1-y)}
		decl vmustar = polygamma(vmu .* dphi, 0)
					 - polygamma((1-vmu).*dphi, 0);          
	        // mu* = psi(mu*phi) - psi((1-mu)*phi)
	 	decl mT = diag( exp(veta) ./ (1 + exp(veta)) .^2 );  
	 	// T = diag{1/g'(mu_1), ..., 1/g'(mu_n)}
		
		// Derivadas com relação ao vetor de parâmetros beta (U_beta)
		(avScore[0])[0] = (dphi * g_mX' * mT * (vystar - vmustar))[2];

		// Derivada com relação a phi (U_phi)
		decl dU_phi = polygamma(dphi, 0) - vmu .* polygamma(vmu .* dphi, 0)
						  - (1-vmu) .* polygamma( (1-vmu).*dphi, 0)
						  + vmu .* log(g_vy) + (1-vmu).*log(1-g_vy);
		(avScore[0])[1] = sumc(dU_phi); 		 	
	}	     	  

	// Retorno da função //

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ){
		return 0;  // indica falha
	} else{
		return 1;  // indica sucesso
	}
}


// ==================================================================== //
//   Declaração da Função de Log-Verossimilhança que entra em MaxBFGS   //
//               no Teste de Hipóteses quando H0 é FALSA                //
// ==================================================================== //

dfRegBetaFALSE(const vP, const adFunc, const avScore, const amHess)
{
	// Declarando variáveis	da Log-Verossimilhança //
	decl dphi = vP[2];
	decl veta = g_mX * (dfalse|vP[:1]);

	// Usando a Logit como função de ligação //
	decl vmu = (exp(veta)) ./ (1 + exp(veta)); 

	// Criando o vetor de Log-Verossimilhanças para cada observação t //
	decl vLogLike_t = loggamma(dphi) - loggamma(vmu .* dphi) 
	                - loggamma((1-vmu) .* dphi) 
			+ ((vmu .* dphi) - 1).*log(g_vy) 
			+ (((1-vmu).*dphi) - 1).*log(1-g_vy);

	// Somando o vetor para obter a Log-Verossimilhança //
	adFunc[0] = sumc(vLogLike_t);

	// Declarando variáveis do gradiente //
	if(avScore)
	{
		decl vystar = log( (g_vy) ./ (1-g_vy) );            
		// y* = log{y/(1-y)}
		decl vmustar = polygamma(vmu .* dphi, 0)
					 - polygamma((1-vmu).*dphi, 0);          
		// mu* = psi(mu*phi) - psi((1-mu)*phi)
	 	decl mT = diag( exp(veta) ./ (1 + exp(veta)) .^2 );  
	 	// T = diag{1/g'(mu_1), ..., 1/g'(mu_n)}
		
		// Derivadas com relação ao vetor de parâmetros beta2 e beta3 
		// (lembrando que Beta1 == dfalse é dado)
		(avScore[0])[0:1] = (dphi * g_mX' * mT * (vystar - vmustar))[1:2];

		// Derivada com relação a phi
		decl dU_phi = polygamma(dphi, 0) - vmu .* polygamma(vmu .* dphi, 0)
						  - (1-vmu) .* polygamma( (1-vmu).*dphi, 0)
						  + vmu .* log(g_vy) + (1-vmu).*log(1-g_vy);
		(avScore[0])[2] = sumc(dU_phi); 		 	

	}	     	  

	// Retorno da função //

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ){
		return 0;  // indica falha
	} else{
		return 1;  // indica sucesso
	}
}	


// ================================================================== //
//   Declaração da Função que gera ocorrências da Distribuição Beta   //
//                  pelo Método da Aceitação-Rejeição                 //
// ================================================================== //

vfBeta_AccRej(const ilen, const vshape1, const vshape2)
{
	/* Notas explicativas:

	-> Y ~ Beta(dshape1, dshape2), X ~ Unif(0,1), U ~ Unif(0,1)
	-> Y tem PDF 'f' da distribuição Beta
	-> X tem PDF 'g' da distribuição Uniforme, i.e., g = 1.
	-> Como M = max (f/g), então M é o valor máximo da distribuição Beta
	-> O valor máximo de uma distribuição ocorre quando x = moda.		 */
	
	decl vY = zeros(ilen, 1);   // Vetor que receberá as ocorrências de Y
	decl dX;                    // Variável que receberá as ocorrências de X
	decl dU;                    // Variável que receberá as ocorrências de U
	decl dfX;                   // PDF de Beta em X
	decl i = 0;                 // Inteiro para indexar e contar

	while (i < ilen)
	{
		// Moda da distribuição Beta //
		decl dmode = (vshape1[i] - 1) / (vshape1[i] + vshape2[i] - 2);
	
		// Valor de Gamma(a+b)/(Gamma(a)*Gamma(b)) //
		decl dgammaconst = gammafact(vshape1[i] + vshape2[i]) / 
		                  (gammafact(vshape1[i]) * gammafact(vshape2[i]));
		
		// O valor de M é o valor de f(moda) //
		decl dM = dgammaconst * dmode ^ (vshape1[i] - 1) * 
		          (1 - dmode) ^ (vshape2[i] - 1);
		
		// Gerando ocorrência de X ~ Unif(0,1) //
		dX = ranu(1,1);
		
		// Avaliando o valor de f(X) (Que é o limite f(X)/g(X), pois g = 1) //
		dfX = dgammaconst * dX ^ (vshape1[i] - 1) * (1 - dX) ^ (vshape2[i] - 1);
		
		// Gerando ocorrência de U ~ Unif(0,1) de forma independente //
		dU = dM * ranu(1,1);

		if (ismissing(dU))
		{
			println("\nNão é possível gerar a Beta pelo algoritmo selecionado. \n",
				"Por favor, escolha outro algoritmo ou mude os parâmetros.");
			exit(1);
		}
		
		// Testando //
		if (dU < dfX)
		{
			// Se U menor que o limite, então aceita-se X como ocorrência de Y //
			vY[i] = dX;  
			i += 1;
		}
	}
	return vY;
}


// ========================================================================== //
//   Declaração da Função que calcula analiticamente a Informação de Fisher   //
// ========================================================================== //

mfKmatrix(const vP)
{
	decl k = 3;
	decl dphik = vP[k]; 
	decl vetak = g_mX * vP[0:(k-1)];
	decl vmuk = exp(vetak) ./ (1 + exp(vetak)); 

	// Trigamma
	decl vtrigam1k = polygamma(vmuk * dphik, 1);	
	decl vtrigam2k = polygamma((1 - vmuk) * dphik, 1);
	decl dtrigam3k = polygamma(dphik, 1);			
	
	// Matrizes T, W, D e vetor c
	decl mT = diag( exp(vetak) ./ (1 + exp(vetak)) .^2 );
	decl mW = diag( dphik * (vtrigam1k + vtrigam2k) ) * mT .^2; 					
	decl mD = diag(vtrigam1k .* (vmuk .^ 2) 
	        + vtrigam2k .* (1 - vmuk) .^2 - dtrigam3k);	
	decl vc = dphik * (vtrigam1k .* vmuk - vtrigam2k .* (1 - vmuk));					

	// Matriz K 
	decl mKbb	= dphik * g_mX' * mW * g_mX;					
	decl mKbp	= g_mX' * mT * vc; 								
	decl dKpp	= trace(mD);									
	decl mK = (mKbb ~ mKbp) | (mKbp' ~ dKpp);					

	return(mK);
}


// Pré-Processador Condicional //

#ifndef DONT_PANIC

main()
{
	// Impressão inicial //
	println("");
	println("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]");
	println("[]                                                        []");
	println("[] Esse arquivo não tem execução/impressão importante.    []");
	println("[] Sua finalidade é fornecer funções ao programa:         []");
	println("[]                                                        []");
	println("[]     beta_regression.ox                                 []");
	println("[]                                                        []");
	println("[] Tente executá-lo.                                      []");
	println("[]                                                        []");
	println("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]");
	println("");
}

#endif

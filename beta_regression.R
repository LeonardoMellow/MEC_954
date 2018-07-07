#!/usr/bin/Rscript

# ============================================================================ #
#   PROGRAMA:   beta_regression.R                                              #
#                                                                              #
#   USO:        Esse script faz simulações de Monte Carlo para estimar parâme- #
#               tros da Regressão Beta proposta por Ferrari e Cribari-Neto     #
#               (2004) com parâmetro de precisão constante.                    #
#                                                                              #
#               A estimação é feita por Estimação Pontual, Estimação Interva-  #
#               lar e por Testes de Hipóteses, esta última fazendo simulação   #
#               de tamanho (impondo H0 verdadeira) e simulação de poder (im-   #
#               pondo H0 falsa e trocando os valores dos parâmetros.           #
#                                                                              #
#               A geração de ocorrências de uma Distribuição Beta pode ser     #
#               feita usando a função rbeta da base do R, ou usando o método   #
#               de rejeição canônico com uma distribuição uniforme como dis-   #
#               tribuição instrumental e mantendo reproducibilidade dos dados  #
#               gerados em outras linguagens.                                  #
#                                                                              #
#   AUTOR:      Leonardo Melo                                                  #
#                                                                              #
#   VERSÃO:     1.0                                                            #
#                                                                              #
#   DATA:       28 de Maio de 2018                                             #
#                                                                              #
#   NOTAS:      Ao executar esse programa, a diretiva de pré-processamento     #
#               condicional DONT_PANIC em "beta_fun.R" não será executada, o-  #
#               cultando a impressão de "beta_fun.R" e deixando apenas a im-   #
#               pressão desse script.                                          #
#                                                                              #
#   COMPILAÇÃO: Rscript beta_regression.R                                      #
#                                                                              #
# ============================================================================ #

# Iniciando o cronômetro
time <- proc.time()

# "Diretiva Condicional" para ocultar impressões de outros programas
DONT_PANIC <- 1

# ================================== #
#    CARREGANDO PACOTES E FUNÇÕES    #
# ================================== #
source("beta_fun.R", chdir = TRUE)
# O arquivo beta_fun.R deve estar no mesmo diretório de beta_regression.R

# =============== #
#    SIMULAÇÃO    #
# =============== #
# Imprimindo Informações na tela
{ cat(" \n")												      	  
  cat("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n")
  cat("[]                                                              []\n")
  cat("[]               RESULTADOS PARA A REGRESSÃO BETA               []\n")
  cat("[]                -  SIMULAÇÃO DE MONTE CARLO -                 []\n")
  cat("[]                                                              []\n")
  cat("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n\n")
}

# # Função Monte Carlo # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
montecarlo <- function(N = 25, R = 1000, B = 250, F_H0Value = 0.7, 
                       SEED = 3112, rbeta_usr = 1, print = 1) {
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
  # Uso #
  # # # # 
  #
  #   N            Tamanho amostral (default é 25).
  #   R            Número de réplicas de Monte Carlo (default é 1000).
  #   B            Número de réplicas de Bootstrap (default é 250). 
  #   F_H0Value    Valor de H0 quando falsa (default é 0.7).
  #
  #   SEED         Semente para garantir reprodutibilidade com Ox
  #
  #   rbeta_usr    Boolean (default é TRUE). Geração da distribuição Beta pelo 
  #                método da Rejeição usando a distribuição Uniforme como dis-
  #                tribuição instrumental. Aviso: menos eficiente!
  #                Se FALSE, a função rbeta é utilizada.
  #
  #   print        Boolean (default é TRUE). Imprime resultados como no Ox.
  #                
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
  
  # # PARÂMETROS GERAIS # #
  vparams <- c(0.85, -0.93, -1.22, 40.5)    # vetor de parâmetros verdadeiros
  ODDS_c <- 0.1                             # incremento na razão de chances 
  vparamsODDS <- exp(ODDS_c*vparams[1:3])   # Razão de chances dos parâm. Beta
  
  fail <- failT <- failF <- 0                
  # Contadores de Falhas de convergência da estimação:
  # 1) da regressão                    2) da regressão nas réplicas de bootstrap
  # 3) do teste quando H0 verdadeira   4) do teste quando H0 falsa
  
  # # VARIÁVEIS DA ESTIMAÇÃO PONTUAL # #
  mestimators <- mestimatorsSE <- mestimators_bs <- matrix(NA, R, 4)
  # mbootstrap <- matrix(NA, B, 4)
  
  # # VARIÁVEIS DA ESTIMAÇÃO INTERVALAR # #
  a10 <- qnorm(1 - 0.10/2); a5 <- qnorm(1 - 0.05/2); a1 <- qnorm(1 - 0.01/2)
  
  # # VALORES CRÍTICOS PARA TH A 10%, 5% e 1% # #              
  VC10_1 <- qchisq(0.90,1); VC5_1 <- qchisq(0.95,1); VC1_1 <- qchisq(0.99,1)
  # 1 restrição => Qui-Quadrado com 1 df (poder)                            
  VC10_2 <- qchisq(0.90,2); VC5_2 <- qchisq(0.95,2); VC1_2 <- qchisq(0.99,2)
  # 2 restrições => Qui-Quadrado com 1 df (tamanho)                         
  
  
  # # VARIÁVEIS INDEPENDENTES # #
  vx2 <- rann(N, SEED); vcoef <- rep(1, N)
  RNGkind("Marsaglia"); set.seed(SEED)
  vx3 <- runif(N)
  s_mX <<- cbind(vcoef, vx2, vx3)  # <<- para ser armazenada na memória
  
  # # PARÂMETROS DA DISTRIBUIÇÃO # #
  veta <- c(s_mX %*% vparams[1:3])
  vmu <- (exp(veta))/(1 + exp(veta))
  dphi <- vparams[4]
  
  
  # # LAÇO MONTE CARLO # #
  nb_cores <- detectCores() - 1
  cl <- makeCluster(nb_cores, type = "FORK")
  registerDoParallel(cl)
  
  par_seed <- matrix(c(rep(401,R), 1:R, 1:R), ncol = R, byrow = TRUE)
  loop <- foreach(j=1:R, .options.RNG=par_seed, .combine = 'c') %dorng% { 
    
    if (rbeta_usr) 
      s_vy <<- beta_accept_reject(N, vmu*dphi, (1 - vmu)*dphi)
    else 
      s_vy <<- rbeta(N, vmu*dphi, (1 - vmu)*dphi)
    
    # Ajustar modelo nesses dados simulados
    tryCatch( {
      fit.sim <- optim(c(0,0,0,1), optLogLike, grLogLike, 
                       method = "L-BFGS-B", lower = c(-5,-5,-5,0.1))
    }, error = function(cond) {
      return(fit.sim <- list(convergence = 42))
    } )
    
    ## PONTUAL ##
    vestimators <- fit.sim$par
    
    ## INTERVALAR ##
    EI <- EInt(vestimators, c(a10, a5, a1), ODDS_c)
    
    # Contando se o parâmetro verdadeiro está dentro do intervalo
    cIC90 <- (vparams < EI$vIC90max & vparams > EI$vIC90min)
    cIC95 <- (vparams < EI$vIC90max & vparams > EI$vIC95min)
    cIC99 <- (vparams < EI$vIC90max & vparams > EI$vIC99min)
    
    cIC90ODDS <- (vparamsODDS < EI$vIC90maxODDS & vparamsODDS > EI$vIC90minODDS)
    cIC95ODDS <- (vparamsODDS < EI$vIC90maxODDS & vparamsODDS > EI$vIC95minODDS)
    cIC99ODDS <- (vparamsODDS < EI$vIC90maxODDS & vparamsODDS > EI$vIC99minODDS)
    
    ## TESTE DE HIPÓTESES ##
    TH <- TH_parallel(fit.sim$value, EI$mFishInv, vestimators, F_H0Value)
    
    if (B) { ## BOOTSTRAP PARAMÉTRICO ##
      mbootstrap <- param_BS(vestimators, SEED, N, B, rbeta_usr)
      
      # Calculando o estimador corrigido de Bootstrap (BC1)
      vestimators_bs <- 2*vestimators - apply(mbootstrap, 2, mean)
    } else {
      vestimators_bs <- NULL
    }
    
    ## RESULTADO ##
    res <- list()
    res[[paste("MC", j, sep = "_")]] <-   
      list(y = s_vy, conv = fit.sim$convergence, 
           convTH_T = TH$convergenceT, convTH_F = TH$convergenceF,
           theta = vestimators, boot = vestimators_bs,
           cIC90 = cIC90, cIC95 = cIC95, cIC99 = cIC99,
           cIC90ODDS = cIC90ODDS, cIC95ODDS = cIC95ODDS, cIC99ODDS = cIC99ODDS,
           vLR1 = TH$vLR1, vSc1 = TH$vSc1, vWa1 = TH$vWa1, 
           vLR2 = TH$vLR2, vSc2 = TH$vSc2, vWa2 = TH$vWa2 )
    res
  }
  mbinded <- matrix(NA, R, 8)
  conv <- convT <- convF <- c()
  
  vcIC90 <- vcIC95 <- vcIC99 <- integer(4)
  vcIC90ODDS <- vcIC95ODDS <- vcIC99ODDS <- integer(3)
  
  vLR1 <- vSc1 <-  vWa1 <- vLR2 <- vSc2 <-  vWa2 <- c()
  for (i in 1:R) {
    conv <- c(conv, loop[[paste0("MC_",i)]]$conv) 
    convT <- c(convT, loop[[paste0("MC_",i)]]$convTH_T)
    convF <- c(convF, loop[[paste0("MC_",i)]]$convTH_F)
    
    if ( conv[i] == 0 & convT[i] == 0 & convF[i] == 0 ) {
      mbinded[i,] <- c(loop[[paste0("MC_",i)]]$theta, loop[[paste0("MC_",i)]]$boot)
      
      vcIC90 <- vcIC90 + loop[[paste0("MC_",i)]]$cIC90
      vcIC95 <- vcIC95 + loop[[paste0("MC_",i)]]$cIC95
      vcIC99 <- vcIC99 + loop[[paste0("MC_",i)]]$cIC99
      
      vcIC90ODDS <- vcIC90ODDS + loop[[paste0("MC_",i)]]$cIC90ODDS
      vcIC95ODDS <- vcIC95ODDS + loop[[paste0("MC_",i)]]$cIC95ODDS
      vcIC99ODDS <- vcIC99ODDS + loop[[paste0("MC_",i)]]$cIC99ODDS
      
      vLR1 <- c(vLR1, loop[[paste0("MC_",i)]]$vLR1)
      vSc1 <- c(vSc1, loop[[paste0("MC_",i)]]$vSc1)
      vWa1 <- c(vWa1, loop[[paste0("MC_",i)]]$vWa1)
      
      vLR2 <- c(vLR2, loop[[paste0("MC_",i)]]$vLR2)
      vSc2 <- c(vSc2, loop[[paste0("MC_",i)]]$vSc2)
      vWa2 <- c(vWa2, loop[[paste0("MC_",i)]]$vWa2)
    }
  }
  
  fail <- sum(conv != 0); failT <- sum(convT != 0); failF <- sum(convF != 0)

  # # ARRUMANDO DADOS PARA IMPRESSÃO # #
  
  # PONTUAL #
  vparam <- t(t(vparams))
  rownames(vparam) <- c("beta1", "beta2", "beta3", "phi")
  colnames(vparam) <- "     Parâm."
  
  vmean <- apply(mbinded, 2, function(x) mean(x, na.rm = T)) 
  vvar <- apply(mbinded, 2, 
                function(x) var(x, na.rm = T)*(sum(!is.na(x)) - 1)/sum(!is.na(x))) # R^{-1} \sum, como no Ox
  
  vest_mean <- t(t(vmean[1:4])); vest_mean_bs <- t(t(vmean[5:8]))
  vest_var <- t(t(vvar[1:4])); vest_var_bs <- t(t(vvar[5:8]))
  vest_bias <- vest_mean - vparam; vest_bias_bs <- vest_mean_bs - vparam
  vest_pbias <- 100*vest_bias/vparam; vest_pbias_bs <- 100*vest_bias_bs/vparam 
  veqm <- vest_bias^2 + vest_var; veqm_bs <- vest_bias_bs^2 + vest_var_bs
  
  pontual <- cbind(round(vest_mean,6), round(vest_bias,6), 
                   round(vest_pbias,6), round(veqm,6))
  rownames(pontual) <- c("beta1", "beta2", "beta3", "phi")
  colnames(pontual) <- c("    Média EMVs", "        Viés ", 
                         " Viés Rel. (%)", "         EQM ")
  if (B) {
    pontual_bs <- cbind(round(vest_mean_bs,6), round(vest_bias_bs,6), 
                        round(vest_pbias_bs,6), round(veqm_bs,6))
    rownames(pontual_bs) <- c("beta1", "beta2", "beta3", "phi")
    colnames(pontual_bs) <- c("    Média EMVs", "        Viés ", 
                              " Viés Rel. (%)", "         EQM ")
  } else {
    pontual_bs <- NULL
  }
  
  # INTERVALAR #
  interv <- cbind( (100*vcIC90/R) %>% round(2) %>% t() %>% t(),
                   (100*vcIC95/R) %>% round(2) %>% t() %>% t(),
                   (100*vcIC99/R) %>% round(2) %>% t() %>% t() )
  rownames(interv) <- c("beta1", "beta2", "beta3", "phi")
  colnames(interv) <- c("      [ 90 % ]","      [ 95 % ]","      [ 99 % ]")
  
  interv2 <- cbind( (100*vcIC90ODDS/R) %>% round(2) %>% t() %>% t(),
                    (100*vcIC95ODDS/R) %>% round(2) %>% t() %>% t(),
                    (100*vcIC99ODDS/R) %>% round(2) %>% t() %>% t() )
  rownames(interv2) <- c("beta1", "beta2", "beta3")
  colnames(interv2) <- c("      [ 90 % ]","      [ 95 % ]","      [ 99 % ]")
  
  ## TESTES DE HIPÓTESES
  vTH10pow <- apply(cbind(vLR1, vSc1, vWa1), 2, 
                    function(x) round(100*sum(x > VC10_1)/R,2))
  vTH5pow <- apply(cbind(vLR1, vSc1, vWa1), 2, 
                   function(x) round(100*sum(x > VC5_1)/R,2))
  vTH1pow <- apply(cbind(vLR1, vSc1, vWa1), 2,
                   function(x) round(100*sum(x > VC1_1)/R,2))
  vTH10siz <- apply(cbind(vLR2, vSc2, vWa2), 2, 
                    function(x) round(100*sum(x > VC10_2)/R,2))
  vTH5siz <- apply(cbind(vLR2, vSc2, vWa2), 2, 
                   function(x) round(100*sum(x > VC5_2)/R,2))
  vTH1siz <- apply(cbind(vLR2, vSc2, vWa2), 2, 
                   function(x) round(100*sum(x > VC1_2)/R,2))
  
  
  # Imprimindo Informações na tela como no Ox
  if (print) {
    cat("==================================================================\n")
    cat("\n "); cat("   TAMANHO DA AMOSTRA: ", "             ", N, "\n")
    cat("   \n "); cat("   Método: ", "                         ", "BFGS \n")
    cat("    Função de Geração Beta: ", "         ", 
        ifelse(rbeta_usr, "Rejeição\n", "rbeta\n"))
    cat("    Algoritmo Geração Uniforme: ", "     ", RNGkind()[1], "\n")
    cat("    Sementes no R: ", "                  ", SEED, "\n")
    set.seed(SEED)
    cat("    Sementes correspondentes no Ox: ", " ", .Random.seed[2], "e",
        .Random.seed[3], "\n")
    cat("   \n "); cat("   Número de Réplicas de MC: ", "       ", R, "\n")
    cat("    Número de Réplicas de Bootstrap:  ", B, "\n")
    cat("\n \n "); cat(" ESTIMAÇÃO PONTUAL\n")              #
    cat("==================================================================\n")
    cat("\n "); print(vparam); cat("\n "); print(pontual)
    if (B) {
      cat("\n "); cat(" Bootstrap \n "); print(pontual_bs)
    }
    cat("\n \n \n"); cat(" ESTIMAÇÃO INTERVALAR\n")
    cat("==================================================================\n")
    cat("\n "); cat("------ Parâmetros Beta1, Beta2, Beta3 e Phi ------\n\n")
    print(interv); cat("\n\n")
    cat(" --- Razão de Chances para Beta1, Beta2 e Beta3 ---\n\n")
    print(interv2); cat("\n \n \n"); cat(" TESTE DE HIPÓTESES \n")
    cat("==================================================================\n")
    cat("\n "); cat("-------------  Simulação de Tamanho  -------------\n\n")
    cat("    Valor crítico 10%:                 ", VC10_2, "\n")
    cat("    Valor crítico 5%:                  ", VC5_2, "\n")
    cat("    Valor crítico 1%:                  ", VC1_2, "\n\n")
    cat("  ---  Probabilidade de cometer o erro tipo I  ---\n\n")
    cat("    Razão de Verossimilhança a 10%:    ", vTH10siz[[1]], "%\n")
    cat("    Razão de Verossimilhança a  5%:    ", vTH5siz[[1]], "%\n")
    cat("    Razão de Verossimilhança a  1%:    ", vTH1siz[[1]], "%\n\n")
    cat("    Teste Escore a 10%:                ", vTH10siz[[2]], "%\n")
    cat("    Teste Escore a  5%:                ", vTH5siz[[2]], "%\n")
    cat("    Teste Escore a  1%:                ", vTH1siz[[2]], "%\n\n")
    cat("    Teste Wald a 10%:                  ", vTH10siz[[3]], "%\n")
    cat("    Teste Wald a  5%:                  ", vTH5siz[[3]], "%\n")
    cat("    Teste Wald a  1%:                  ", vTH1siz[[3]], "%\n\n\n\n")
    cat("-------------   Simulação de Poder   -------------\n\n")
    cat("    Valor crítico 10%:                 ", VC10_1, "\n")
    cat("    Valor crítico 5%:                  ", VC5_1, "\n")
    cat("    Valor crítico 1%:                  ", VC1_1, "\n\n")
    cat("  ---    Probabilidade de rejeitar H0 falsa    ---\n\n")
    cat("     ---   Hipótese Nula:     Beta1 =", F_H0Value,"  ---\n\n")
    cat("    Razão de Verossimilhança a 10%:    ", vTH10pow[[1]], "%\n")
    cat("    Razão de Verossimilhança a  5%:    ", vTH5pow[[1]], "%\n")
    cat("    Razão de Verossimilhança a  1%:    ", vTH1pow[[1]], "%\n\n")
    cat("    Teste Escore a 10%:                ", vTH10pow[[2]], "%\n")
    cat("    Teste Escore a  5%:                ", vTH5pow[[2]], "%\n")
    cat("    Teste Escore a  1%:                ", vTH1pow[[2]], "%\n\n")
    cat("    Teste Wald a 10%:                  ", vTH10pow[[3]], "%\n")
    cat("    Teste Wald a  5%:                  ", vTH5pow[[3]], "%\n")
    cat("    Teste Wald a  1%:                  ", vTH1pow[[3]], "%\n\n")
    cat("==================================================================\n\n")
    cat("  Num falhas Convergência: Regressão:        ", fail, "de ", R, "\n")
    cat("  Num falhas Convergência: Teste H0 True:    ", failT, "de ", R, "\n")
    cat("  Num falhas Convergência: Teste H0 False:   ", failF, "de ", R,"\n\n")
    cat("==================================================================\n\n")
  }
  
  stopCluster(cl)
  set.seed(SEED)
  results <- list(
    N = N, R = R, URNG = RNGkind()[1], SEED = SEED, 
    OxSEEDs = c(.Random.seed[2],.Random.seed[3]), 
    EP1 = pontual, EP_bs = pontual_bs, EI = interv, EIodds = interv2, 
    TH10 = matrix(c(vTH10pow,vTH10siz), ncol = 3, byrow = T,
                  dimnames = list(c("Poder","Tamanho"),
                                  c("LR", "Score", "Wald"))),
    TH5 = matrix(c(vTH5pow,vTH5siz), ncol = 3, byrow = T,
                 dimnames = list(c("Poder","Tamanho"),
                                 c("LR", "Score", "Wald"))),
    TH1 = matrix(c(vTH1pow,vTH1siz), ncol = 3, byrow = T,
                 dimnames = list(c("Poder","Tamanho"),
                                 c("LR", "Score", "Wald"))),
    VC = matrix(c(VC10_1,VC10_2,VC5_1,VC5_2,VC1_1,VC1_2), ncol = 3,
                dimnames = list(c("Poder","Tamanho"),
                                c("VC10", "VC5", "VC1"))),
    fail = matrix(c(fail, failT, failF), ncol = 1, 
                  dimnames = list(c("Regressão","Teste TRUE",
                                    "Teste FALSE"),
                                  c("Falhas")))
  )
  return(results)
}


N25  <- montecarlo(R = 1000, B = 250, N = 25)
N50  <- montecarlo(R = 1000, B = 0,   N = 50)
N100 <- montecarlo(R = 1000, B = 0,   N = 100)
N150 <- montecarlo(R = 1000, B = 0,   N = 150)
N200 <- montecarlo(R = 1000, B = 0,   N = 200)
N250 <- montecarlo(R = 1000, B = 0,   N = 250)

## Imprimir data e tempo de execução
end <- proc.time() - time
cat( " Data: ", as.character(Sys.Date()), "\n")
cat( " Tempo de execução: ", round(end[[3]], 2) , "\n\n")
cat( " Fim.\n")
cat("==================================================================\n\n")


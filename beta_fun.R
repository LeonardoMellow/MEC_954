#!/usr/bin/Rscript

# ============================================================================ #
#   PROGRAMA:   beta_fun.R                                                     #
#                                                                              #
#               !!! IMPORTANTE !!!                                             #
#                                                                              #
#               Todos os arquivos de extensão .R devem estar no mesmo diretó-  #
#               rio.                                                           #
#                                                                              #
#   USO:        Esse script serve para carregar funções úteis aos demais pro-  #
#               gramas. Nele constam: uma função para calcular o valor da log- #
#               verossimilhança; uma função para geração de ocorrências da Dis-#
#               tribuição Beta pelo método canônico da rejeição; e uma função  #
#               para geração de ocorrências da Distribuição Normal pelo método #
#               polar, usado em Ox, a fim de garantir a reproducibilidade dos  #
#               dados.                                                         #
#                                                                              #
#   AUTOR:      Leonardo Melo                                                  #
#                                                                              #
#   VERSÃO:     1.0                                                            #
#                                                                              #
#   DATA:       28 de Maio de 2018                                             #
#                                                                              #
#   COMPILAÇÃO: Rscript beta_fun.R                                             #
#                                                                              #
# ============================================================================ #


## Definição de "Diretiva Condicional" para ocultar impressões deste programa
if (!exists("DONT_PANIC")) {
  DONT_PANIC <- 0
}


# ============================================================================ #
#         FUNÇÃO PARA CHECAR EXISTÊNCIA, INSTALAR E CARREGAR PACOTES           # 
# ============================================================================ #

check.packages <- function(pkg) {
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, quietly = TRUE,
         warn.conflicts = FALSE, character.only = TRUE)
}


## Carregando Pacotes
check.packages(c("dplyr", "data.table", "betareg", "foreach",
                 "doRNG", "parallel", "doParallel", "doMC"))


## ==================== ##
##    FUNÇÃO LOGLIKE    ##
## ==================== ##

LogLike <- function(mu, phi, s_vy) {
  
  fn <- lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
        ((mu * phi) - 1) * log(s_vy) + (((1 - mu) * phi) - 1) * log(1 - s_vy) 
  
  return(sum(fn))
}


## ======================================== ##
##    FUNÇÃO LOGLIKE PARA USAR EM optim     ##
## ======================================== ##

 ## OBS: Ela difere da função LogLike(mu, phi, s_vy) por 2 motivos:  # # # # # # 
 #                                                                             #
 #  1º:  Só recebe um chute inicial dos parâmetros como argumento, sendo ne-   #
 #       cessário que a matriz X e o vetor y estejam na memória.
 #  2º:  Como a função optim faz minimização por default, é necessário colocar #
 #       o sinal negativo em seu retorno. Poderia também ser utilizado o parâ- #
 #       metro control(list = c(fnscale = -1)) na função optim.                #
 #                                                                             #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

optLogLike <- function(par) {
  
  dphi <- par[length(par)]
  vbeta <- par[1:(length(par) - 1)] %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  return(-LogLike(vmu, dphi, s_vy))
}

optLogLikeBS <- function(par) { # para uso no Bootstrap
  
  dphi <- par[length(par)]
  vbeta <- par[1:(length(par) - 1)] %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  return(-LogLike(vmu, dphi, y_bs))
}

optLogLikeTRUE <- function(par) { # para uso no Teste de H0 Verdadeira
  
  dphi <- par[length(par)]
  vbeta <- c(0.85, -0.93, par[1:(length(par) - 1)]) %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  return(-LogLike(vmu, dphi, s_vy))
}

optLogLikeFALSE <- function(par) {# para uso no Teste de H0 Falsa
  
  dphi <- par[length(par)]
  vbeta <- c(F_H0Value, par[1:(length(par) - 1)]) %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  return(-LogLike(vmu, dphi, s_vy))
}


## ================================================= ##
##    GRADIENTE DA optLOGLIKE PARA A FUNÇÃO optim    ##
## ================================================= ##

grLogLike <- function(par) {
  
  dphi <- par[length(par)]
  vbeta <- par[1:(length(par) - 1)] %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  vystar <- log( (s_vy)/(1 - s_vy) )
  vmustar <- psigamma(vmu * dphi) - psigamma((1 - vmu) * dphi)
  mT <- diag( exp(veta) / (1 + exp(veta))^2 )
  
  # Derivadas parciais (com relação a beta e a phi)
  dUbeta <- dphi * (t(s_mX) %*% mT %*% (vystar - vmustar))
  dUphi <- sum(psigamma(dphi) - vmu * psigamma(vmu * dphi) -
                  (1 - vmu) * psigamma( (1 - vmu) * dphi ) +
                  vmu * log(s_vy) + (1 - vmu) * log(1 - s_vy))
  
  return(-c(dUbeta, dUphi)) 
}

grLogLikeBS <- function(par) { # para uso no Bootstrap
  
  dphi <- par[length(par)]
  vbeta <- par[1:(length(par) - 1)] %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  vystar <- log( (y_bs)/(1 - y_bs) )
  vmustar <- psigamma(vmu * dphi) - psigamma((1 - vmu) * dphi)
  mT <- diag( exp(veta) / (1 + exp(veta))^2 )
  
  # Derivadas parciais (com relação a beta e a phi)
  dUbeta <- dphi * (t(s_mX) %*% mT %*% (vystar - vmustar))
  dUphi <- sum(psigamma(dphi) - vmu * psigamma(vmu * dphi) -
                 (1 - vmu) * psigamma( (1 - vmu) * dphi ) +
                 vmu * log(y_bs) + (1 - vmu) * log(1 - y_bs))
  
  return(-c(dUbeta, dUphi)) 
}

grLogLikeTRUE <- function(par) { # para uso no Teste de H0 Verdadeira
  
  dphi <- par[length(par)]
  vbeta <- c(0.85, -0.93, par[1:(length(par) - 1)]) %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  vystar <- log( (s_vy)/(1 - s_vy) )
  vmustar <- psigamma(vmu * dphi) - psigamma((1 - vmu) * dphi)
  mT <- diag( exp(veta) / (1 + exp(veta))^2 )
  
  # Derivadas parciais (com relação a beta e a phi)
  dUbeta <- dphi * (t(s_mX) %*% mT %*% (vystar - vmustar))[3]
  dUphi <- sum(psigamma(dphi) - vmu * psigamma(vmu * dphi) -
                 (1 - vmu) * psigamma( (1 - vmu) * dphi ) +
                 vmu * log(s_vy) + (1 - vmu) * log(1 - s_vy))
  
  return(-c(dUbeta, dUphi)) 
}

grLogLikeFALSE <- function(par) { # para uso no Teste de H0 Falsa
  
  dphi <- par[length(par)]
  vbeta <- c(F_H0Value, par[1:(length(par) - 1)]) %>% as.matrix()
  veta <- (s_mX %*% vbeta) %>% as.numeric()
  vmu <- exp(veta) / (1 + exp(veta))
  
  vystar <- log( (s_vy)/(1 - s_vy) )
  vmustar <- psigamma(vmu * dphi) - psigamma((1 - vmu) * dphi)
  mT <- diag( exp(veta) / (1 + exp(veta))^2 )
  
  # Derivadas parciais (com relação a beta e a phi)
  dUbeta <- dphi * (t(s_mX) %*% mT %*% (vystar - vmustar))[2:3]
  dUphi <- sum(psigamma(dphi) - vmu * psigamma(vmu * dphi) -
                 (1 - vmu) * psigamma( (1 - vmu) * dphi ) +
                 vmu * log(s_vy) + (1 - vmu) * log(1 - s_vy))
  
  return(-c(dUbeta, dUphi)) 
}


## ===================================================== ##
##   FUNÇÃO PARA ENCONTRAR A MATRIZ DA INFO DE FISHER    ##
## ===================================================== ##

mfKmatrix <- function(vP)
{
  vetak <- (s_mX %*% vP[1:3]) %>% as.numeric()
  vmuk <- exp(vetak) / (1 + exp(vetak))
  dphik <- vP[4]
  
  vtrigam1 <- trigamma(vmuk*dphik)
  vtrigam2 <- trigamma((1 - vmuk)*dphik)
  vtrigam3 <- trigamma(dphik)
  
  mTk <- diag(c(exp(vetak) / (1 + exp(vetak))^2))
  mWk <- diag(c(dphik * (vtrigam1 + vtrigam2)) ) * mTk^2;
  mDk <- diag(c(vtrigam1 * (vmuk ^ 2) + vtrigam2 * (1 - vmuk)^2 - vtrigam3))
  vck <- dphik * (vtrigam1 * vmuk - vtrigam2 * (1 - vmuk))
  
  mKbbk <- dphik * (t(s_mX) %*% mWk %*% s_mX)
  mKbpk <- t(s_mX) %*% mTk %*% vck
  dKppk <- sum(diag(mDk))
  mKk <- rbind( cbind(mKbbk, mKbpk), cbind(t(mKbpk), dKppk) )

  return(mKk)
}


# ============================================================================ #
#         MÉTODO DA ACEITAÇÃO / REJEIÇÃO CANÔNICO, USANDO UMA UNIFORME         # 
# ============================================================================ #

beta_accept_reject <- function(N, shape1, shape2) {
  # Explicação mais detalhada dessa função encontra-se em beta_fun.ox
  Y <- NULL
  i <- 1
  while (length(Y) < N) {
    Mode <- (shape1[i] - 1)/(shape1[i] + shape2[i] - 2) 
    M <- (1/beta(shape1[i], shape2[i])) * Mode^(shape1[i] - 1) * (1 - Mode)^(shape2[i] - 1)
    X <- runif(1)
    fX <- (1/beta(shape1[i], shape2[i])) * X^(shape1[i] - 1) * (1 - X)^(shape2[i] - 1)
    U <- M * runif(1)
    if (is.na(U)) {
      stop("\t\tNão é possível gerar a Beta pelo algoritmo selecionado. \n
                Por favor, escolha outro algoritmo ou mude os parâmetros.")
    }
    if (U < fX)  {
      Y <- c(Y,X)
      i <- i + 1
    }
  }
  return(Y)
}


## Função análoga à função rann do Ox quando o gerador uniforme é o GM.
## OBS: Método Polar para gerar dados normais
rann <- function(N, SEED) {
  
  RNGkind("Marsaglia")
  set.seed(SEED)
  X <- c()
  while (length(X) < N) {
    U1 <- runif(1)
    U2 <- runif(1)
    V1 <- 2*U1 - 1
    V2 <- 2*U2 - 1
    W <- V1^2 + V2^2
    if (W <= 1) 
      X <- c(X, V1*sqrt(-(2*log(W))/W), V2*sqrt(-(2*log(W))/W))
  }
  return(X[1:N])
  
}


# ============================================================================ #
#              FUNÇÕES PARA DEIXAR O CÓDIGO PRINCIPAL MAIS LIMPO               # 
# ============================================================================ #

param_BS <- function(vP, SEED, N, B, rbeta_usr) {
  veta_hat <- c(s_mX %*% vP[1:3])
  vmu_hat <- (exp(veta_hat))/(1 + exp(veta_hat))
  dphi_hat <- vP[4]
  i <- 0
  
  matrixbootstrap <- matrix(NA, B, 4)
  
  set.seed(SEED)
  while (i < B) {
    if (rbeta_usr)
      y_bs <<- beta_accept_reject(N, vmu_hat*dphi_hat, (1 - vmu_hat)*dphi_hat)
    else
      y_bs <<- rbeta(N, vmu_hat*dphi_hat, (1 - vmu_hat)*dphi_hat)
    
    tryCatch(
      {
        fit.bs <- optim(c(0,0,0,1), optLogLikeBS, grLogLikeBS,
                        method = "L-BFGS-B", lower = c(-5,-5,-5,0.1))
      },
      error = function(cond) {
        return(fit.bs <- list(convergence = 42))
      }
    )
    
    if (fit.bs$convergence == 0) {
      i <- i + 1
      matrixbootstrap[i,] <- fit.bs$par
    } 
  }
  return(matrixbootstrap)
}  

EInt <- function(vP, vVC, ODDS_c) {
  a10 <- vVC[1]; a5 <- vVC[2]; a1 <- vVC[3]
  
  # Calculado o Erro Padrão
  mKest <- mfKmatrix(vP)            # Informação de Fisher estimada
  mFishInv <- solve(mKest)          # Inversa da Informação de Fisher
  vSE <- sqrt(diag(mFishInv))       # Vetor de Erros Padrão
  
  # Intervalos de Confiança dos Parâmetros
  vIC90min <- vP - a10 * vSE; vIC90max <- vP + a10 * vSE
  vIC95min <- vP -  a5 * vSE; vIC95max <- vP +  a5 * vSE
  vIC99min <- vP -  a1 * vSE; vIC99max <- vP +  a1 * vSE
  
  # Intervalos de Confiança da Razão de Chances
  vIC90minODDS <- exp(ODDS_c*(vP[1:3] - a10 * vSE[1:3]))
  vIC90maxODDS <- exp(ODDS_c*(vP[1:3] + a10 * vSE[1:3]))
  
  vIC95minODDS <- exp(ODDS_c*(vP[1:3] - a5 * vSE[1:3]))
  vIC95maxODDS <- exp(ODDS_c*(vP[1:3] + a5 * vSE[1:3]))
  
  vIC99minODDS <- exp(ODDS_c*(vP[1:3] - a1 * vSE[1:3]))
  vIC99maxODDS <- exp(ODDS_c*(vP[1:3] + a1 * vSE[1:3]))
  
  return(list(vIC90min = vIC90min, vIC95min = vIC95min, vIC99min = vIC99min,
              vIC90max = vIC90max, vIC95max = vIC95max, vIC99max = vIC99max,
              vIC90minODDS = vIC90minODDS, vIC90maxODDS = vIC90maxODDS, 
              vIC95minODDS = vIC95minODDS, vIC95maxODDS = vIC95maxODDS, 
              vIC99minODDS = vIC99minODDS, vIC99maxODDS = vIC99maxODDS,
              mFishInv = mFishInv))
}

TH_parallel <- function(LL_Value, mFishInv, vP, F_H0Value) {

    # # # # # # # # # # # # #
  #   Com H0 verdadeira   # # # # # # # # # # # # # # # # # # # # # # # # # 
  #                                                                       #
  #   As hipóteses do Teste são (2 restrições):                           #
  #                                                                       #
  #     H0:     beta1  =  0.85       vs.      H1:      beta1  !=  0.85    #
  #         (e) beta2  = -0.93                    (ou) beta2  != -0.93    #
  #                                                                       #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # # RAZÃO DE VEROSSIMILHANÇA # #
  dloglike_irr <- -LL_Value
  
  tryCatch(
    {
      rest_optimT <- optim(c(0,1), optLogLikeTRUE, grLogLikeTRUE, 
                       method = "L-BFGS-B", lower = c(-5,-5,-5,0.1))
    },
    error = function(cond) {
      return(rest_optimT <- list(convergence = 42))
    }
  )
  
  dloglike_res <- -rest_optimT$value
  vLR2 <- 2*(dloglike_irr - dloglike_res)
  
  # # TESTE ESCORE # #
  vrest <- c(0.85, -0.93, rest_optimT$par)
  vU <- -grLogLike(vrest) 
  mK <- mfKmatrix(vrest)
  mFishInv_rest <- solve(mK)
  
  vSc2 <- c(t(vU[1:2]) %*% mFishInv_rest[1:2,1:2] %*% t(t(vU[1:2])))
  
  # # TESTE WALD # #
  mK_1 <- solve(mFishInv[1:2, 1:2])
  vWa2 <- c(t(vP[1:2] - vrest[1:2]) %*% 
              mK_1 %*% 
              t(t(vP[1:2] - vrest[1:2])))
  
  # # # # # # # # # # #
  #   Com H0 falsa.   # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  #                                                                       #
  #   A hipótese do Teste é (1 restrição):                                #
  #                                                                       #
  #     H0:     beta1  =  F_H0Value    vs.  H1:    beta1  !=  F_H0Value   #
  #                                                                       #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # # RAZÃO DE VEROSSIMILHANÇA # #
  # Guardando na memória
  F_H0Value <<- F_H0Value
  
  tryCatch(
    {
      rest_optimF <- optim(c(0,0,1), optLogLikeFALSE, grLogLikeFALSE, 
                       method = "L-BFGS-B", lower = c(-5,-5,-5,0.1))
    },
    error = function(cond) {
      return(rest_optimF <- list(convergence = 42))
    }
  )
  
  dloglike_res <- -rest_optimF$value
  vLR1 <- 2*(dloglike_irr - dloglike_res)
  
  # # TESTE ESCORE # #
  vrest <- c(F_H0Value, rest_optimF$par)
  vU <- -grLogLike(vrest)
  mK <- mfKmatrix(vrest)
  mFishInv_rest <- solve(mK)
  
  vSc1 <- c(t(vU[1]) %*% mFishInv_rest[1,1] %*% t(t(vU[1])))
  
  # # TESTE WALD # #
  mK_1 <- solve(mFishInv[1,1])
  vWa1 <- c(t(vP[1] - vrest[1]) %*% 
              mK_1 %*% 
              t(t(vP[1] - vrest[1])))
  
  return(list(vLR1 = vLR1, vSc1 = vSc1, vWa1 = vWa1,
              vLR2 = vLR2, vSc2 = vSc2, vWa2 = vWa2,
              convergenceT = rest_optimT$convergence,
              convergenceF = rest_optimF$convergence))

}

## Pré-Processador Condicional

if (!DONT_PANIC) {

  cat("\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n")
  cat("[]                                                        []\n")
  cat("[] Esse arquivo não tem execução/impressão importante.    []\n")
  cat("[] Sua finalidade é fornecer funções aos programas:       []\n")
  cat("[]                                                        []\n")
  cat("[]     beta_regression.R                                  []\n")
  cat("[]     beta_regression_comparison.R                       []\n")
  cat("[]                                                        []\n")
  cat("[] Tente executar um deles.                               []\n")
  cat("[]                                                        []\n")
  cat("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n\n")
    
}



 #### Algoritmo #2 do artigo de Simar e Wilson (2007).
 ## Autor: Lucas Cardoso Queiroz

 #### Estimation and inference in two-stage, semi-parametric models
 ###  of production processes. Journal of Econometrics. Vol.136, pp.31-64

 #### Comparação com a função dea.env.robust do pacote rDEA (Estimates bias-corrected efficiency
 ###  scores in input- and output-oriented DEA models with environmental (exogenous) variables)
 

 rm(list=ls(all=TRUE))

 library(truncnorm)
 library(Benchmarking)

 ## Leitura da base de dados 
 dt01 <- ### Dados das variáveis ambientais
 dt02 <- ### Dados utilizados no modelo DEA

 dados <- merge(dt01, dt02, by=c("DMU","year"))

 
 library(maxLik)
 Z <- as.matrix(dados["variaveis_Z"])
 X <- as.matrix(dados["DEA_inputs"])
 Y <- as.matrix(dados["DEA_outputs"])

 dados$delta <- Benchmarking::dea(X, Y, RTS="crs", ORIENTATION="out")$eff
 modelo.lm   <- lm(delta ~ variaveis_Z, data=dados)
 zi          <- as.matrix(model.matrix(modelo.lm))

 ## Ajuste do modelo linear
 dados2       <- dados[which(dados$delta > 1),]
 modelo.lm2   <- lm(delta ~ vel.vento100, data=dados2)
 zi2          <- as.matrix(model.matrix(modelo.lm2))

 loglik <- function(param, Y, Z) {
   sig.2 <- param[1]
   beta <- as.matrix(param[-1])
  
   erro <- Y - (Z%*%beta)
   t.a  <- 1 - (Z%*%beta)
  
   ll <- sum( log( dtruncnorm(erro, a=t.a, b=Inf, mean=0, sd = sqrt(sig.2)) ) )
   return(ll)
 }

 ## Uso do pacote maxLik
 solucao <- maxLik(loglik, 
                  start=c(anova(modelo.lm2)["Residuals","Mean Sq"], coef(modelo.lm2)), 
                  Y=dados2$delta, Z=zi2)

 beta    <- as.matrix( solucao$estimate[-1] )
 sig2    <- solucao$estimate[1]

 ## Loop L1 vezes:
 L.1           <- 100
 delta_est     <- matrix(NA, 48, L.1)
 delta_est_hat <- matrix(NA, 48, L.1)
 t.a  <- 1 - (zi%*%beta)

 for(i in 1:L.1){
   
   erro <- rtruncnorm(n=1, a=t.a, b=Inf, mean = 0, sd = sqrt(sig2))
  
   delta_est[,i] <- zi %*% beta + erro
   Y2  <- Y * (dados$delta/delta_est[,i])	
   delta_est_hat[,i] <- Benchmarking::dea(X, Y2, RTS="crs", ORIENTATION="out")$eff
 }

 ## Usando o pacote rDEA para comparação
 library(rDEA)
 di_env   <- rDEA::dea.env.robust(X, Y, Z = Z, 
                                 model="output", 
                                 RTS="constant",
                                 L1=L.1, 
                                 L2=2000)

 ## Estimativa do vício (ESTIMADO - REAL)
 Bias <- delta_est_hat - delta_est

 ## summary(apply(Bias, 1, mean) - di_env$bias)
 plot(apply(Bias, 1, mean) ~ di_env$bias, pch=19, col="red")
 abline(a=0, b=1, lwd=2, lty=2, col="blue")

 ## di_env$delta_hat
 ## di_env$bias
 ## di_env$delta_hat - di_env$bias ## Primeiro calcula o "bias" e depois corrige
 ## di_env$delta_hat_hat

 plot(di_env$delta_hat_hat ~ I(dados$delta - apply(Bias, 1, mean)), 
      pch=19, col="red")
 abline(a=0, b=1, lwd=2, lty=2)


 ## plot(1/(dados$delta - apply(Bias, 1, mean)) ~ I(1/dados$delta),
 ##      pch=19, col="red")
 ## abline(a=0, b=1, lwd=2, lty=3)

 ## Passo 5
 modelo.lm3 <- lm(di_env$delta_hat_hat ~ variaveis_Z, data=dados)
 solucao2   <- maxLik(loglik, 
                  start=c(anova(modelo.lm3)["Residuals","Mean Sq"], coef(modelo.lm3)), 
                  Y=di_env$delta_hat_hat, Z=zi)

 beta_hh    <- as.matrix( solucao2$estimate[-1] )
 sig2_hh    <- solucao2$estimate[1]
 sd2_hh     <- sqrt(sig2_hh)

 ## Passo 6
 L.2 <- 2000 ## 2000 no artigo
 beta_hat_est <- matrix(NA, L.2, 2)
 sig_hat_est <- c(NA, L.2)
 t.a2  <- 1 - (zi%*%beta_hh)

 for(i in 1:L.2){
  erro <- rtruncnorm(n=1, a=t.a2, b=Inf, mean = 0, sd = sd2_hh)
  delta_est_est <- zi %*% beta_hh + erro
  
  solucao      <- maxLik(loglik, 
                       start=c(anova(modelo.lm3)["Residuals","Mean Sq"], coef(modelo.lm3)), 
                       Y=delta_est_est, Z=zi)
  
  beta_hat_est[i,]   <- as.matrix( solucao$estimate[-1] )
  sig_hat_est[i]     <- solucao$estimate[1]
}

 colMeans(beta_hat_est)
 di_env$beta_hat_hat_star

 sqrt(mean(sig_hat_est)) 
 di_env$sigma_hat_hat_star

 di_env$beta_ci
 quantile(beta_hat_est[,1], probs = c(0.025,0.975))
 quantile(beta_hat_est[,2], probs = c(0.025,0.975))

 di_env$sigma_ci
 quantile(sqrt(sig_hat_est), probs = c(0.025,0.975))

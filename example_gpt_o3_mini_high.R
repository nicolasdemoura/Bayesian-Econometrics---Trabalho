# -------------------------------
# SCRIPT: Full Bayesian Analysis for Plausibly Exogenous IV Model
# Baseado no Apêndice A.3 de Conley et al. (2012)
# -------------------------------

# Carregar os pacotes necessários
library(MASS)      # Para mvrnorm (amostras da normal multivariada)
library(MCMCpack)  # Para riwish (amostra da Inverse-Wishart)

# -------------------------------
# 1. Função de Simulação dos Dados
# -------------------------------
simulate_data <- function(n, beta, pi, gamma, Sigma) {
  # n: número de observações
  # beta: coeficiente de x na equação estrutural (escalar)
  # pi: coeficiente da primeira etapa (escalar)
  # gamma: vetor de dimensão 2 = (gamma1, gamma2) para o efeito direto dos instrumentos em y 
  #        (note: o primeiro elemento multiplica a constante, o segundo multiplica z)
  # Sigma: matriz 2x2 de covariância dos erros (v e ε)
  
  # Gerar o vetor de instrumentos:
  #   Cada observação tem Z_i = (1, z_i), onde z_i ~ N(0,1)
  z <- rnorm(n)
  Z <- cbind(1, z)
  
  # Gerar os erros bivariados (v, ε) para cada observação
  errors <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  # errors[,1] = v, errors[,2] = ε
  
  # Equação de primeira etapa: x = (z de segunda coluna)*pi + v
  x <- Z[, 2] * pi + errors[, 1]
  
  # Equação estrutural: y = beta * x + (1*gamma1 + z*gamma2) + ε
  y <- beta * x + Z %*% gamma + errors[, 2]
  
  # Retorna um data.frame com as variáveis
  data.frame(y = y, x = x, Z1 = Z[,1], Z2 = Z[,2])
}

# -------------------------------
# 2. Função do Gibbs Sampler (Full Bayesian)
# -------------------------------
mcmc_plausibly_exogenous <- function(data, n_iter = 10000, burn_in = 2000,
                                     mu_pi = 0, V_pi = 100,
                                     mu_theta = rep(0, 3), V_theta = diag(c(100, 1, 1)),
                                     nu0 = 3, S0 = diag(2)) {
  # data: data.frame com colunas y, x, Z1, Z2
  # n_iter: número total de iterações do sampler
  # burn_in: iterações a descartar
  # Priors para pi ~ N(mu_pi, V_pi)
  # Priors para theta = (beta, gamma1, gamma2) ~ N(mu_theta, V_theta)
  # Prior para Sigma ~ Inverse-Wishart(nu0, S0)
  
  n <- nrow(data)
  
  # Extração dos dados:
  y <- data$y
  x <- data$x
  Z1 <- data$Z1  # coluna de 1's
  Z2 <- data$Z2  # variável instrumental relevante para a primeira etapa
  
  # Preparar a matriz do segundo estágio: W = [x, Z1, Z2]
  # Na equação estrutural: y = beta*x + gamma1*Z1 + gamma2*Z2 + ε
  W <- cbind(x, Z1, Z2)
  
  # Número de parâmetros:
  # pi é escalar; theta é 3x1; Sigma é 2x2
  n_pi <- 1
  n_theta <- 3
  
  # Inicializar os parâmetros
  pi_current <- mu_pi            # Inicia com o valor do prior para pi
  theta_current <- mu_theta      # (beta, gamma1, gamma2)
  Sigma_current <- S0            # Inicializar com a escala do prior (pode ser ajustado)
  
  # Armazenar as amostras
  pi_draws     <- numeric(n_iter)
  theta_draws  <- matrix(NA, nrow = n_iter, ncol = n_theta) 
  Sigma_draws  <- vector("list", n_iter)
  
  # Guardar as primeiras iterações
  pi_draws[1] <- pi_current
  theta_draws[1, ] <- theta_current
  Sigma_draws[[1]] <- Sigma_current
  
  # Inverso da variância a priori para pi (escalar)
  V_pi_inv <- 1 / V_pi
  
  # Inverso da matriz de covariância a priori para theta (3x3)
  V_theta_inv <- solve(V_theta)
  
  # Loop do Gibbs Sampler
  for (iter in 2:n_iter) {
    # -------------------------------
    # Passo 1: Atualizar pi (primeira etapa)
    # Modelo: x = Z2 * pi + v,  v ~ N(0, sigma_v^2), onde sigma_v^2 = Sigma_current[1,1]
    sigma_v2 <- Sigma_current[1, 1]
    V_pi_post <- 1 / (V_pi_inv + sum(Z2^2) / sigma_v2)
    mu_pi_post <- V_pi_post * (V_pi_inv * mu_pi + sum(Z2 * x) / sigma_v2)
    pi_current <- rnorm(1, mean = mu_pi_post, sd = sqrt(V_pi_post))
    
    # -------------------------------
    # Passo 2: Atualizar theta = (beta, gamma1, gamma2) (segunda etapa)
    # Modelo: y = x*beta + Z1*gamma1 + Z2*gamma2 + ε,  ε ~ N(0, sigma_eps^2), onde sigma_eps^2 = Sigma_current[2,2]
    sigma_eps2 <- Sigma_current[2, 2]
    # Posterior: theta ~ N(mu_theta_post, V_theta_post), onde
    V_theta_post <- solve(V_theta_inv + t(W) %*% W / sigma_eps2)
    mu_theta_post  <- V_theta_post %*% (V_theta_inv %*% mu_theta + t(W) %*% y / sigma_eps2)
    theta_current <- as.vector(mvrnorm(1, mu = mu_theta_post, Sigma = V_theta_post))
    
    # -------------------------------
    # Passo 3: Atualizar Sigma (erro conjunto de ambas as equações)
    # Para cada observação, calcular resíduos:
    #   Resíduo da primeira equação: r1 = x - Z2*pi
    #   Resíduo da segunda equação: r2 = y - (x*beta + gamma1*Z1 + gamma2*Z2)
    r1 <- x - Z2 * pi_current
    r2 <- y - (x * theta_current[1] + Z1 * theta_current[2] + Z2 * theta_current[3])
    # Matriz soma dos quadrados dos resíduos (SSE)
    SSE <- t(cbind(r1, r2)) %*% cbind(r1, r2)
    # Parâmetros do posterior Inverse-Wishart:
    nu_post <- nu0 + n
    S_post <- S0 + SSE
    # Amostrar Sigma do posterior:
    Sigma_current <- riwish(nu_post, S_post)
    
    # Armazenar as amostras desta iteração
    pi_draws[iter] <- pi_current
    theta_draws[iter, ] <- theta_current
    Sigma_draws[[iter]] <- Sigma_current
  }
  
  # Remover burn-in
  list(pi_draws = pi_draws[(burn_in + 1):n_iter],
       theta_draws = theta_draws[(burn_in + 1):n_iter, ],
       Sigma_draws = Sigma_draws[(burn_in + 1):n_iter])
}

# -------------------------------
# 3. Script Principal: Simulação e Execução do MCMC
# -------------------------------

# Parâmetros "verdadeiros" para simulação
n <- 150
beta_true <- 2            # Efeito causal
pi_true <- 0.25           # Coeficiente da primeira etapa
gamma_true <- c(0, 0.1)   # Efeitos diretos dos instrumentos: (gamma1, gamma2)
Sigma_true <- matrix(c(1, 0.5,
                       0.5, 1), nrow = 2, byrow = TRUE)  # Matriz de covariância dos erros

# Simular os dados
set.seed(123)
data_sim <- simulate_data(n, beta_true, pi_true, gamma_true, Sigma_true)

# Priors para a análise Bayesiana
# Para pi ~ N(mu_pi, V_pi)
mu_pi_prior <- 0
V_pi_prior <- 100
# Para theta = (beta, gamma1, gamma2) ~ N(mu_theta, V_theta)
mu_theta_prior <- c(0, 0, 0)
V_theta_prior <- diag(c(100, 1, 1))  # Priorizando maior incerteza sobre beta e maior precisão para os desvios (γ)
# Para Sigma ~ Inverse-Wishart(nu0, S0)
nu0_prior <- 3         # Mínimo para matriz 2x2 é 3
S0_prior  <- diag(2)    # Escala padrão (pode ser ajustada)

# Definir número de iterações e burn-in
n_iter <- 10000
burn_in <- 2000

# Executar o Gibbs sampler
results <- mcmc_plausibly_exogenous(data_sim, n_iter = n_iter, burn_in = burn_in,
                                    mu_pi = mu_pi_prior, V_pi = V_pi_prior,
                                    mu_theta = mu_theta_prior, V_theta = V_theta_prior,
                                    nu0 = nu0_prior, S0 = S0_prior)

# Exibir resumos dos parâmetros
cat("Posterior para pi:\n")
print(summary(results$pi_draws))

cat("\nPosterior para theta = (beta, gamma1, gamma2):\n")
print(apply(results$theta_draws, 2, summary))

cat("\nPosterior para Sigma (primeira iteração após burn-in):\n")
print(results$Sigma_draws[[1]])





# -----------------------------------------------------------
# 4. Diagnóstico Simples
# -----------------------------------------------------------

# (1) Carregar ou já ter em memória o objeto "results"
#     que contém as amostras do Gibbs:
# - results$pi_draws       -> vetor de amostras para π
# - results$theta_draws    -> matriz de amostras para θ = (β, γ₁, γ₂)
# - results$Sigma_draws    -> lista de matrizes 2×2 para Σ
#
# Exemplo para carregar:
# load("gibbs_results.RData")

# (2) Extrair amostras
pi_draws    <- results$pi_draws
theta_draws <- results$theta_draws
Sigma_draws <- results$Sigma_draws

# (3) Para Σ, extrair as três entradas de interesse
sigma11 <- sapply(Sigma_draws, function(S) S[1,1])
sigma12 <- sapply(Sigma_draws, function(S) S[1,2])
sigma22 <- sapply(Sigma_draws, function(S) S[2,2])

# -----------------------------------------------------------
# A. Trace plots + ACF
# -----------------------------------------------------------

# == Parâmetro π ==
par(mfrow = c(1,2))
plot(pi_draws, type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ pi),
     xlab = "Iteração", ylab = expression(pi))
acf(pi_draws, main = bquote("ACF: " ~ pi))

# == Parâmetro β ==
par(mfrow = c(1,2))
plot(theta_draws[,1], type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ beta),
     xlab = "Iteração", ylab = expression(beta))
acf(theta_draws[,1], main = bquote("ACF: " ~ beta))

# == Parâmetro γ[1] ==
par(mfrow = c(1,2))
plot(theta_draws[,2], type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ gamma[1]),
     xlab = "Iteração", ylab = expression(gamma[1]))
acf(theta_draws[,2], main = bquote("ACF: " ~ gamma[1]))

# == Parâmetro γ[2] ==
par(mfrow = c(1,2))
plot(theta_draws[,3], type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ gamma[2]),
     xlab = "Iteração", ylab = expression(gamma[2]))
acf(theta_draws[,3], main = bquote("ACF: " ~ gamma[2]))

# == Elementos de Σ (Sigma[1,1], Sigma[1,2], Sigma[2,2]) ==
par(mfrow = c(1,2))
plot(sigma11, type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ Sigma[1,1]),
     xlab = "Iteração", ylab = expression(Sigma[1,1]))
acf(sigma11, main = bquote("ACF: " ~ Sigma[1,1]))

par(mfrow = c(1,2))
plot(sigma12, type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ Sigma[1,2]),
     xlab = "Iteração", ylab = expression(Sigma[1,2]))
acf(sigma12, main = bquote("ACF: " ~ Sigma[1,2]))

par(mfrow = c(1,2))
plot(sigma22, type = "l", col = "blue",
     main = bquote("Trace Plot: " ~ Sigma[2,2]),
     xlab = "Iteração", ylab = expression(Sigma[2,2]))
acf(sigma22, main = bquote("ACF: " ~ Sigma[2,2]))

# -----------------------------------------------------------
# B. Histogramas
# -----------------------------------------------------------
par(mfrow = c(3,3))

# π
hist(pi_draws, col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ pi),
     xlab = expression(pi))

# β
hist(theta_draws[,1], col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ beta),
     xlab = expression(beta))

# γ[1]
hist(theta_draws[,2], col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ gamma[1]),
     xlab = expression(gamma[1]))

# γ[2]
hist(theta_draws[,3], col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ gamma[2]),
     xlab = expression(gamma[2]))

# Sigma[1,1]
hist(sigma11, col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ Sigma[1,1]),
     xlab = expression(Sigma[1,1]))

# Sigma[1,2]
hist(sigma12, col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ Sigma[1,2]),
     xlab = expression(Sigma[1,2]))

# Sigma[2,2]
hist(sigma22, col = "lightblue", border = "white",
     main = bquote("Histograma: " ~ Sigma[2,2]),
     xlab = expression(Sigma[2,2]))

# Duas "células" vazias para fechar a grade 3x3
plot.new(); plot.new()

# Fim do script

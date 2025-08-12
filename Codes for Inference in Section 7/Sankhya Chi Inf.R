########################################
#Codes for Sankhya Paper
require(distr)
require(coda)
########################################
#0. Density of chi-square distribution#########################
WS_pdf <- function(nu1, nu2, lam1, lam2, sigma2, y) {
  add <- convpow(2*sigma2*Chisq(df = nu1, ncp = lam1) + 
                   sigma2*Chisq(df = nu2, ncp = lam2), 1)
  as2pdf <- d(add)
  return(as2pdf(y))
}
#1. Read in data and convert to array#######################
#1.1 Cleft#####################################
Cle_fir <- read.csv("Cle_fir.csv", header=FALSE)
#1.2 Control####################################
Con_fir <- read.csv("Con_fir.csv", header=FALSE)
#1.3 Convert the data to array###############################
#1.3.1 Cleft##############################
cle_fir_array <- array(0, dim = c(24, 3, 13))
for (subject in 1:13) {
  cle_fir_array[, , subject] <- matrix(as.matrix(Cle_fir[subject, ]),
            nrow = 24, ncol = 3)
}
#1.3.2 Control###############################
con_fir_array <- array(0, dim = c(24, 3, 12))
for (subject in 1:12) {
  con_fir_array[, , subject] <- matrix(as.matrix(Con_fir[subject, ]),
            nrow = 24, ncol = 3)
}
#1.4 Compute sample mean shapes################################
#1.4.1 Cleft#######################################
Cle_fir_mean <- apply(cle_fir_array, c(1, 2), mean)
#1.4.2 Control#######################################
Con_fir_mean <- apply(con_fir_array, c(1, 2), mean)
###########################################
#2. Compute the AS measure defined in (7.2)#######################
#2.1 Cleft#######################################
cle_fir_L2 <- c()
for (subject in 1:13) {
  pairlm <- sum((cle_fir_array[c(1:6, 20:24), 1, subject] + 
                   cle_fir_array[c(13:8, 18:14), 1, subject])^2 + 
                  (cle_fir_array[c(1:6, 20:24), 2, subject] - 
                     cle_fir_array[c(13:8, 18:14), 2, subject])^2 + 
                  (cle_fir_array[c(1:6, 20:24), 3, subject] - 
                     cle_fir_array[c(13:8, 18:14), 3, subject])^2)
  
  sololm <- cle_fir_array[7, 1, subject]^2 + 
    cle_fir_array[19, 1, subject]^2
  
  cle_fir_L2[subject] <- pairlm + sololm
}
#2.2 Control###################################
con_fir_L2 <- c()
for (subject in 1:12) {
  pairlm <- sum((con_fir_array[c(1:6, 20:24), 1, subject] + 
                   con_fir_array[c(13:8, 18:14), 1, subject])^2 + 
      (con_fir_array[c(1:6, 20:24), 2, subject] - 
         con_fir_array[c(13:8, 18:14), 2, subject])^2 + 
      (con_fir_array[c(1:6, 20:24), 3, subject] - 
         con_fir_array[c(13:8, 18:14), 3, subject])^2)
  
  sololm <- con_fir_array[7, 1, subject]^2 + 
    con_fir_array[19, 1, subject]^2
  
  con_fir_L2[subject] <- pairlm + sololm
}
#3. Inference under Matrix Normality########################
#3.1 Estimate sigma^2##############################
#3.1.1 Cleft#####################################
clefir_subtrmean <- apply(cle_fir_array, 3, 
              function(a) {a - Cle_fir_mean}, simplify = F)
#merge them into a long vector
clefir_subtrmean_vec <- c()
for (i in 1:13) {
  clefir_subtrmean_vec <- append(clefir_subtrmean_vec, 
            matrix(clefir_subtrmean[[i]], ncol = 1))
}

sig2_clefir <- var(clefir_subtrmean_vec)
#3.1.2 Control#####################################
confir_subtrmean <- apply(con_fir_array, 3, 
              function(a) {a - Con_fir_mean}, simplify = F)
#merge them into a long vector
confir_subtrmean_vec <- c()
for (i in 1:12) {
  confir_subtrmean_vec <- append(confir_subtrmean_vec, 
            matrix(confir_subtrmean[[i]], ncol = 1))
}

sig2_confir <- var(confir_subtrmean_vec)
#3.2 MoM estimators given by (6.1)########################
#3.2.1 Cleft##########################################
lam1_cle_fir_MoM <- (4*sig2_clefir*mean(cle_fir_L2) - var(cle_fir_L2) -
        (4*2*sig2_clefir*sig2_clefir - 2*(2*sig2_clefir)^2)*33 - 
  2*sig2_clefir^2*2)/(4*2*sig2_clefir*sig2_clefir - 4*(2*sig2_clefir)^2)

lam2_cle_fir_MoM <- 
  (mean(cle_fir_L2) - 2*sig2_clefir*(33 + lam1_cle_fir_MoM) -
     sig2_clefir*2)/sig2_clefir
#3.2.2 Control##########################################
lam1_con_fir_MoM <- (4*sig2_confir*mean(con_fir_L2) - var(con_fir_L2) -
        (4*2*sig2_confir*sig2_confir - 2*(2*sig2_confir)^2)*33 - 
  2*sig2_confir^2*2)/(4*2*sig2_confir*sig2_confir - 4*(2*sig2_confir)^2)

lam2_con_fir_MoM <- 
  (mean(con_fir_L2) - 2*sig2_confir*(33 + lam1_con_fir_MoM) -
     sig2_confir*2)/sig2_confir
#3.3 MoM estimators given by (6.4)########################
#3.3.1 Cleft#########################################
lambda1_clefir <- sum((Cle_fir_mean[c(1:6, 20:24), 1] + 
    Cle_fir_mean[c(13:8, 18:14), 1])^2 + (Cle_fir_mean[c(1:6, 20:24), 2] - 
    Cle_fir_mean[c(13:8, 18:14), 2])^2 + (Cle_fir_mean[c(1:6, 20:24), 3] - 
    Cle_fir_mean[c(13:8, 18:14), 3])^2)/(2*sig2_clefir)

lambda2_clefir <- (Cle_fir_mean[7, 1]^2 + Cle_fir_mean[19, 1]^2)/sig2_clefir
#3.3.2 Control#############################################
lambda1_confir <- sum((Con_fir_mean[c(1:6, 20:24), 1] + 
    Con_fir_mean[c(13:8, 18:14), 1])^2 + (Con_fir_mean[c(1:6, 20:24), 2] - 
    Con_fir_mean[c(13:8, 18:14), 2])^2 + (Con_fir_mean[c(1:6, 20:24), 3] - 
    Con_fir_mean[c(13:8, 18:14), 3])^2)/(2*sig2_confir)

lambda2_confir <- (Con_fir_mean[7, 1]^2 + Con_fir_mean[19, 1]^2)/sig2_confir
#3.4 Expectation and variance of AS#########################
#3.4.1 Cleft############################################
AS_clefir_E <- 2*28.03*(33 + 0.034) + 28.03*2

AS_clefir_Var <- 8*28.03^2*(33 + 2*0.034) + 2*28.03^2*2
#3.4.2 Control############################################
AS_confir_E <- 2*5.15*(33 + 0.94) + 5.15*2

AS_confir_Var <- 8*5.15^2*(33 + 2*0.94) + 2*5.15^2*2
#3.5 Density plots of AS for the two groups########################
plot(seq(0, 5000, 0.01), WS_pdf(nu1 = 33, nu2 = 2, 
          lam1 = lambda1_confir, lam2 = 0, 
          sigma2 = sig2_confir, y = seq(0, 5000, 0.01)), type = 'l',
     xlab = 'z', ylab = 'density', lwd = 2,
     main = 'Distribution of AS')
lines(seq(0, 5000, 0.01), WS_pdf(nu1 = 33, nu2 = 2, 
          lam1 = lambda1_clefir, lam2 = 0, 
          sigma2 = sig2_clefir, y = seq(0, 5000, 0.01)),
      lty = 2, lwd = 2)
#3.6 The direct measure of asymmetry given in (7.7)##################
#3.6.1 Compute the measure###########################################
ASmu_cle_fir <- 2*0.034*28.03

ASmu_con_fir <- 2*0.94*5.15
#3.6.2 Permutation test####################################
#Test statistic computed on observations
dlam1_clecon <- abs(lambda1_clefir - lambda1_confir)

#Merge the cleft and control data together
cle_con_arr <- array(0, dim = c(24, 3, 25))
cle_con_arr[, , 1:13] <- cle_fir_array
cle_con_arr[, , 14:25] <- con_fir_array

dlam1_permutation <- c()
for (index in 1:10000) {
  #permute the samples
  id <- sample(1:25, 25, replace = F)
  cleft_id <- id[1:13]
  control_id <- id[14:25]
  
  cleft_data <- cle_con_arr[, , cleft_id]
  control_data <- cle_con_arr[, , control_id]
  
  cleft_data_mean <- apply(cleft_data, c(1, 2), mean)
  control_data_mean <- apply(control_data, c(1, 2), mean)
  
  #Compute estimate of sigma^2
  #Cleft
  cle_subtrmean_per <- apply(cleft_data, 3, 
            function(a) {a - cleft_data_mean}, simplify = F)
  #merge them into a long vector
  cle_subtrmean_vec_per <- c()
  for (i in 1:13) {
    cle_subtrmean_vec_per <- append(cle_subtrmean_vec_per, 
          matrix(cle_subtrmean_per[[i]], ncol = 1))
  }
  
  sig2_cle_per <- var(cle_subtrmean_vec_per)
  
  #Control
  con_subtrmean_per <- apply(control_data, 3, 
            function(a) {a - control_data_mean}, simplify = F)
  #merge them into a long vector
  con_subtrmean_vec_per <- c()
  for (i in 1:12) {
    con_subtrmean_vec_per <- append(con_subtrmean_vec_per, 
          matrix(con_subtrmean_per[[i]], ncol = 1))
  }
  
  sig2_con_per <- var(con_subtrmean_vec_per)
  
  #Estimate of lambda1
  #Cleft
  lam1_cle_per <- sum((cleft_data_mean[c(1:6, 20:24), 1] + 
    cleft_data_mean[c(13:8, 18:14), 1])^2 + (cleft_data_mean[c(1:6, 20:24), 2] - 
    cleft_data_mean[c(13:8, 18:14), 2])^2 + (cleft_data_mean[c(1:6, 20:24), 3] - 
    cleft_data_mean[c(13:8, 18:14), 3])^2)/(2*sig2_cle_per)
  
  #Control
  lam1_con_per <- sum((control_data_mean[c(1:6, 20:24), 1] + 
    control_data_mean[c(13:8, 18:14), 1])^2 + (control_data_mean[c(1:6, 20:24), 2] - 
    control_data_mean[c(13:8, 18:14), 2])^2 + (control_data_mean[c(1:6, 20:24), 3] - 
    control_data_mean[c(13:8, 18:14), 3])^2)/(2*sig2_con_per)
  
  #Test statistics
  dlam1_permutation[index] <- abs(lam1_cle_per - lam1_con_per)
}

mean(dlam1_permutation >= dlam1_clecon)
###################################################
#4. Inference under Ex-chisq distribution######################
#4.1 Cleft###########################################
#4.1.1 Simulated annealing###########################
cle_fir_SA <- matrix(0, nrow = 100001, ncol = 2)
beta0 <- 1
alpha0 <- 1.005
cle_fir_SA[1, 1] <- runif(1, 0, 2)
cle_fir_SA[1, 2] <- runif(1, 0, 10)
for (i in 1:100000) {
  betat <- alpha0^i*beta0
  
  para_prop <- c(runif(1, 0, 2), runif(1, 0, 10))
  
  a <- min(1, exp(-betat*
     (-sum(log(WS_pdf(33, 2, para_prop[1], 0, para_prop[2], 
                      cle_fir_L2))) - 
        -sum(log(WS_pdf(33, 2, cle_fir_SA[i, 1],
      0, cle_fir_SA[i, 2], cle_fir_L2))))))
  
  u <- runif(1)
  if (u < a) {
    cle_fir_SA[i+1, ] <- para_prop
  }
  
  else {
    cle_fir_SA[i+1, ] <- cle_fir_SA[i, ]
  }
}
#4.1.2 Check trace plots#########################
traceplot(as.mcmc(cle_fir_SA[, 1]))
traceplot(as.mcmc(cle_fir_SA[, 2]))
#4.1.3 Direct measure of asymmetry given in (7.7)###############
ASmu_cle_fir_Exchisq <- 2*0.58*0.66
#4.1.4 Expectation and variance of AS########################
AS_clefir_E_Exchisq <- 2*0.66*(33 + 0.58) + 0.66*2

AS_clefir_Var_Exchisq <- 8*0.66^2*(33 + 2*0.58) + 2*0.66^2*2
#4.1.5 Log-likelihood################################
sum(log(WS_pdf(33, 2, 0.58, 0, 0.66, cle_fir_L2)))
#4.2 Control###########################################
#4.2.1 Simulated annealing###########################
con_fir_SA <- matrix(0, nrow = 100001, ncol = 2)
beta0 <- 1
alpha0 <- 1.005
con_fir_SA[1, 1] <- runif(1, 0, 1)
con_fir_SA[1, 2] <- runif(1, 0, 2)
for (i in 1:100000) {
  betat <- alpha0^i*beta0
  
  para_prop <- c(runif(1, 0, 1), runif(1, 0, 2))
  
  a <- min(1, exp(-betat*
     (-sum(log(WS_pdf(33, 2, para_prop[1], 0, para_prop[2], 
                      con_fir_L2))) - 
        -sum(log(WS_pdf(33, 2, con_fir_SA[i, 1],
      0, con_fir_SA[i, 2], con_fir_L2))))))
  
  u <- runif(1)
  if (u < a) {
    con_fir_SA[i+1, ] <- para_prop
  }
  
  else {
    con_fir_SA[i+1, ] <- con_fir_SA[i, ]
  }
}
#4.2.2 Check trace plots#########################
traceplot(as.mcmc(con_fir_SA[, 1]))
traceplot(as.mcmc(con_fir_SA[, 2]))
#4.2.3 Direct measure of asymmetry given in (7.7)###############
ASmu_con_fir_Exchisq <- 2*0.23*0.28
#4.2.4 Expectation and variance of AS########################
AS_confir_E_Exchisq <- 2*0.28*(33 + 0.23) + 0.28*2

AS_confir_Var_Exchisq <- 8*0.28^2*(33 + 2*0.23) + 2*0.28^2*2
#4.2.5 Log-likelihood################################
sum(log(WS_pdf(33, 2, 0.23, 0, 0.28, con_fir_L2)))
########################################################
#5. Inference under Ex-Gaussian distribution######################
#5.1 Density function for Ex-Gaussian###################
exGaussian_pdf <- function(z, a, b, mu, sig2, c) {
  #z is the input, a and b are weights,
  #mu and sig2 are the mean and variance for the Gaussian
  #c is the parameter for the exponential distribution
  
  alpha <- a*mu + c*a^2*sig2/b
  beta <- sqrt(2*a^2*sig2)
  
  part1 <- exp(c*a*mu/b + c^2*a^2*sig2/(2*b^2))
  part2 <- dexp(z, c/b)
  part3 <- pnorm(sqrt(2)*(z - alpha)/beta)
  return(part1*part2*part3)
}
#5.2 Cleft###########################################
#5.2.1 Simulated annealing###########################
cle_fir_SA_ExGau <- matrix(0, nrow = 100001, ncol = 2)
beta0 <- 1
alpha0 <- 1.005
cle_fir_SA_ExGau[1, 1] <- runif(1, 0, 1)
cle_fir_SA_ExGau[1, 2] <- rgamma(1, 5, 1)
for (i in 1:100000) {
  betat <- alpha0^i*beta0
  
  para_prop <- c(runif(1, 0, 1), rgamma(1, 5, 1))
  
  a <- min(1, exp(-betat*
     (-sum(log(exGaussian_pdf(cle_fir_L2, 2*para_prop[2],
    para_prop[2], 33 + para_prop[1], 2*(33 + 2*para_prop[1]), 1/2))) - 
      -sum(log(exGaussian_pdf(cle_fir_L2, 2*cle_fir_SA_ExGau[i, 2],
    cle_fir_SA_ExGau[i, 2], 33 + cle_fir_SA_ExGau[i, 1], 
    2*(33 + 2*cle_fir_SA_ExGau[i, 1]), 1/2))))))
  
  u <- runif(1)
  if (u < a) {
    cle_fir_SA_ExGau[i+1, ] <- para_prop
  }
  
  else {
    cle_fir_SA_ExGau[i+1, ] <- cle_fir_SA_ExGau[i, ]
  }
}
#5.2.2 Check trace plots#########################
traceplot(as.mcmc(cle_fir_SA_ExGau[, 1]))
traceplot(as.mcmc(cle_fir_SA_ExGau[, 2]))
#5.2.3 Direct measure of asymmetry given in (7.7)###############
ASmu_cle_fir_ExGau <- 2*0.48*0.75
#5.2.4 Expectation and variance of AS########################
AS_clefir_E_ExGau <- 2*0.75*(33 + 0.48) + 0.75*2

AS_clefir_Var_ExGau <- 8*0.75^2*(33 + 2*0.48) + 2*0.75^2*2
#5.2.5 Log-likelihood################################
sum(log(WS_pdf(33, 2, 0.48, 0, 0.75, cle_fir_L2)))
#5.3 Control###########################################
#5.3.1 Simulated annealing###########################
con_fir_SA_ExGau <- matrix(0, nrow = 100001, ncol = 2)
beta0 <- 1
alpha0 <- 1.005
con_fir_SA_ExGau[1, 1] <- runif(1, 0, 1)
con_fir_SA_ExGau[1, 2] <- rgamma(1, 5, 1)
for (i in 1:100000) {
  betat <- alpha0^i*beta0
  
  para_prop <- c(runif(1, 0, 1), rgamma(1, 5, 1))
  
  a <- min(1, exp(-betat*
     (-sum(log(exGaussian_pdf(con_fir_L2, 2*para_prop[2],
    para_prop[2], 33 + para_prop[1], 2*(33 + 2*para_prop[1]), 1/2))) - 
      -sum(log(exGaussian_pdf(con_fir_L2, 2*con_fir_SA_ExGau[i, 2],
    con_fir_SA_ExGau[i, 2], 33 + con_fir_SA_ExGau[i, 1], 
    2*(33 + 2*con_fir_SA_ExGau[i, 1]), 1/2))))))
  
  u <- runif(1)
  if (u < a) {
    con_fir_SA_ExGau[i+1, ] <- para_prop
  }
  
  else {
    con_fir_SA_ExGau[i+1, ] <- con_fir_SA_ExGau[i, ]
  }
}
#5.3.2 Check trace plots#########################
traceplot(as.mcmc(con_fir_SA_ExGau[, 1]))
traceplot(as.mcmc(con_fir_SA_ExGau[, 2]))
#5.3.3 Direct measure of asymmetry given in (7.7)###############
ASmu_con_fir_ExGau <- 2*0.51*0.35
#5.3.4 Expectation and variance of AS########################
AS_confir_E_ExGau <- 2*0.35*(33 + 0.51) + 0.35*2

AS_confir_Var_ExGau <- 8*0.35^2*(33 + 2*0.51) + 2*0.35^2*2
#5.3.5 Log-likelihood################################
sum(log(WS_pdf(33, 2, 0.51, 0, 0.35, con_fir_L2)))













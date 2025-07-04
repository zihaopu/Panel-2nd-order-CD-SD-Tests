#library(doMC)

library(MASS)
library(ggplot2)
library(reshape2)
library(MASS)
library(Rsolnp)
#install.packages("ggplot2")
#library(ggplot2)
set.seed(2024)

# Y = read.csv("C:\\Users\\puzihao\\Desktop\\empricial\\Y.csv")
# X = read.csv("C:\\Users\\puzihao\\Desktop\\empricial\\X.csv")
# Y = as.matrix(Y)
# X = as.matrix(X)

###### 1. data #########
Y = ymat

N = nrow(Y)

T = ncol(Y)

m = 3 #y_t-1, x_1, x_2

X = array(0,dim = c(N,T,m-1))

X[,,1] = xmat1

X[,,2] = xmat2


######### 2. setttings #########
q1 = 5

q2 = 1.75


######## 3. functions ##############
marchtest  = function (zt, lag) 
{
  if (!is.matrix(zt)) 
    zt = as.matrix(zt)
  nT = dim(zt)[1]
  k = dim(zt)[2]
  C0 = cov(zt)
  zt1 = scale(zt, center = TRUE, scale = FALSE)
  C0iv = ginv(C0)
  wk = zt1 %*% C0iv
  wk = wk * zt1
  rt2 = apply(wk, 1, sum) - k
  m1 = acf(rt2, lag.max = lag, plot = F)
  acf = m1$acf[2:(lag + 1)]
  c1 = c(1:lag)
  deno = rep(nT, lag) - c1
  Q = sum(acf^2/deno) * nT * (nT + 2)
  pv1 = 1 - pchisq(Q, lag)
  
  rk = rank(rt2)
  m2 = acf(rk, lag.max = lag, plot = F)
  acf = m2$acf[2:(lag + 1)]
  mu = -(rep(nT, lag) - c(1:lag))/(nT * (nT - 1))
  v1 = rep(5 * nT^4, lag) - (5 * c(1:lag) + 9) * nT^3 + 9 * 
    (c(1:lag) - 2) * nT^2 + 2 * c(1:lag) * (5 * c(1:lag) + 
                                              8) * nT + 16 * c(1:lag)^2
  v1 = v1/(5 * (nT - 1)^2 * nT^2 * (nT + 1))
  QR = sum((acf - mu)^2/v1)
  pv2 = 1 - pchisq(QR, lag)
  
  ret_march = matrix(0,2,2)
  ret_march[1,] = c(Q,pv1)
  ret_march[2,] = c(QR,pv2)
  return(ret_march)
  
}

lingtest = function(zt,lag){
  # 1 l_1 norm
  N_test = nrow(zt)
  T_test = ncol(zt)
  #norm 1 zt
  zt_norm1 = apply(abs(zt),2,sum) # c(||X1||,||x2||, …… ,||Xt||)
  y_bar =  mean(zt_norm1)
  zt_norm1_ybar = matrix(0,1,T_test)
  for (j in 1:T_test) {
    zt_norm1_ybar[1,j] = zt_norm1[j] - y_bar
  }
  
  gamma_all = matrix(0,lag+1,1)
  for (l in 0:lag) {
    for (t in (l+1):T_test) {
      gamma_all[l+1,1] = gamma_all[l+1,1] + zt_norm1_ybar[1,t] * zt_norm1_ybar[1,t-l]
    }
  }
  gamma_all = gamma_all / T_test
  rou_all = gamma_all[2:(lag+1),] / gamma_all[1,]
  
  # statistic T_n
  T_n = T_test * sum(rou_all^2)
  # statistic W_n
  W_n = 0
  for (l in 1:lag) {
    W_n = W_n + (lag-l+1) / lag * rou_all[l]^2 / (T_test - l)
  }
  W_n = W_n * T_test*(T_test + 2)
  pv1 = 1-pchisq(T_n, lag)
  pv2 = 1 - pgamma(W_n, shape = 3/4*lag*(lag+1)^2/(2*lag^2 + 3*lag+1), scale = 2/3 *(2*lag^2 + 3*lag+1)/(lag*(lag+1)))
  
  ret_ling = matrix(0,2,2)
  ret_ling[1,] = c(T_n,pv1)
  ret_ling[2,] = c(W_n,pv2)
  
  return(ret_ling)
}

selectM_norm = function(lagall,q,K){
  #penalty
  if(max(lagall[-1]) <= sqrt(log(N*T)*q)){
    pi = log(N*T)
  }else{
    pi = 2
  }
  
  #find m
  Q_m = 0
  L_max <- 0
  M_out = 1
  for (m in 1:K) {
    Q_m = Q_m + lagall[m+1]^2
    L_m <- Q_m - pi * m
    if(L_m > L_max){
      L_max <- L_m
      M_out = m
    }
  }
  return(M_out)
}

selectM_chisq = function(lagall,lagtime,lagcross,q1,q2,K){
  #penalty
  if(max(lagtime[-1]) <= sqrt(log(N*T)*q1) && max(lagcross[-1]) <= sqrt(log(N*T)*q2)){
    pi = log(N*T)
  }else{
    pi = 2
  }
  
  #find m
  Q_m = 0
  L_max <- 0
  M_out = 1
  for (m in 1:K) {
    Q_m = Q_m + lagall[m+1]
    L_m <- Q_m - pi * m
    if(L_m > L_max){
      L_max <- L_m
      M_out = m
    }
  }
  return(M_out)
}

lagk = function(matt, k_1){
  matt_tilde = matrix(0,N,k_1)
  matt_t = cbind(matt_tilde,matt[,1:(T-k_1)])
  matt_final = matt %*% t(matt_t)
  autocov_t = sum(diag(matt_final))
  autocov_k = sum(matt_final) - autocov_t
  return(c(autocov_k,autocov_t))
}

lagk_baltagi = function(matt){
  matt_tilde = t(matt)
  matt_final = matt %*% matt_tilde
  diag_matt = sqrt(diag(matt_final))
  matt_final <- t(t(matt_final) / diag_matt)
  matt_final <- matt_final / diag_matt
  matt_final <- T*matt_final^2 - 1
  diag2 = sum(diag(matt_final))
  auto = (sum(matt_final) - diag2) / 2
  LMp = sqrt(1/N/(N-1)) * auto
  LMbc = LMp - N / (2*(T-1))
  LMbc2 = LMbc^2
  pv = 1 - pchisq(LMbc2,1)
  return(c(LMbc2,pv))
}

CDtest = function(zt) {
  N_test = nrow(zt)
  T_test = ncol(zt)
  
  cross_prod = tcrossprod(zt)
  norms = sqrt(rowSums(zt^2))  # Changed to rowSums
  norm_matrix = outer(norms, norms)
  
  rho = cross_prod / norm_matrix
  diag(rho) = 0
  
  CD = sqrt((2 * T_test) / (N_test * (N_test - 1))) * sum(rho[lower.tri(rho)])
  
  
  pv = 2 * pnorm(-abs(CD))
  
  return(c(CD,pv))
}

##### 4. estimation ######
l_N<-matrix(rep(1,N))

l_T<-l<-matrix(rep(1,T))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T<-diag(T) - l_T%*%t(l_T)/T

index_N = c(1:N)

y = array(0,dim = c(N,(T+1),m))

y[,2:(T+1),1] = Y 
y[,1:T,2:m] = X

upper = 0
lower = 0
for(j in 1:N){
  mean_j = matrix(c(mean(y[j,2:(T+1),1]), apply(y[j,1:T,],2,mean)[-1]), m, 1)
  for (t in 3:(T-1)) {
    t_1_ = matrix(y[j,t-1,] - mean_j, m, 1) 
    t_ = matrix(y[j,t,] - mean_j, m, 1) 
    upper = upper + t_1_ %*% t(t_1_)
    lower = lower + t_1_ %*% t(t_)
  }
}

theta_all = solve(upper) %*% lower

est_gamma = theta_all[,1]

est_phi = est_gamma[1]
est_beta = matrix(est_gamma[2:m],m-1,1)

est_V_mat = matrix(0,N,T)
est_V_mat[,1] = Y[,1] - X[,1,]%*%est_beta 
for (t in 2:T) {
  est_V_mat[,t] = Y[,t] - X[,t,]%*%est_beta - Y[,t-1]*est_phi
}
est_mu = apply(est_V_mat,1,mean)
est_mu = matrix(est_mu,N,1)
est_U_mat = matrix(0,N,T)
for (t in 1:T) {
  est_U_mat[,t] = est_V_mat[,t] - est_mu
}

est_Y =matrix(0,N,T)
est_Y[,1] = X[,1,]%*%est_beta + est_mu
for (t in 2:T) {
  est_Y[,t] = X[,t,]%*%est_beta + Y[,t-1]*est_phi + est_mu
}



####### 5 test for CH ###########
####### 5.1 preparation for test########

# est_U_mat = GDP_diff
# N = nrow(GDP)
# T = ncol(GDP)

# est_U_mat = House_diff
# N = nrow(House_diff)
# T = ncol(House_diff)

index_N = c(1:N)

# for time 
est_U_2 = est_U_mat^2
est_U_4 = est_U_mat^4
est_sigma_2 = apply(est_U_2, 1, mean)
est_tau_4 = apply(est_U_4, 1, mean)
est_U_sd = matrix(0,N,T)
for (j in 1:N) {
  est_U_sd[j,] = est_U_2[j,] - est_sigma_2[j] # u_it^2 - est_sigma^2 
}

#sigma_bar_8 = 4*mean(est_sigma_2^4)
#sigma_bar_4 = 2*mean(est_sigma_2^2)
sigma_bar_8 = mean((est_tau_4 - est_sigma_2^2)^2)
sigma_bar_4 = mean(est_tau_4 - est_sigma_2^2)

var_time = sigma_bar_8# / sigma_bar_4^2

# for cross
sigma_4 = est_tau_4 - est_sigma_2^2
#sigma_4 = 2 * est_sigma_2^2
#omega_4 = omega ^ 2
var_cross_up = 0
var_cross_low = 0
for (j1 in 1:N) {
  index_i1 = index_N[-j1]
  for (j2 in index_i1) {
    var_cross_up = var_cross_up + sigma_4[j1] * sigma_4[j2]
    #sigma_all = sigma_all + omega_4[j1] * omega_4[j2]
    #var_cross_low = var_cross_low + sqrt(sigma_4[j1] * sigma_4[j2])
  }
}
var_cross_up = var_cross_up / ( (N-1)*N )
#var_cross_low = var_cross_low / ( (N-1)*N )

var_cross = var_cross_up #/ (var_cross_low)^2

##### 6. Compute rho 1~M #####

M = T-10

C_temp = matrix(0, M+1, 1)
C_time = matrix(0, M+1, 1)
C_time_nobias = matrix(0, M+1, 1)

rho_time = matrix(0, M+1, 1)
rho_time_nobias = matrix(0, M+1, 1)
rho_cross = matrix(0, M+1, 1)
rho_joint = matrix(0, M+1, 1)
rho_joint_nobias = matrix(0, M+1, 1)

for (k in 0:M) {
  
  value_autocov = lagk(est_U_sd, k)
  C_temp[(k+1),1] = value_autocov[1]
  C_time[(k+1),1] = value_autocov[2]
  C_time_nobias[(k+1),1] = value_autocov[2]
  
  #### 6.1 for cross rho ####
  C_temp[k+1, 1] = C_temp[(k+1),1] / sqrt(N*(N-1)*(T-k)) #/ normalize_cross
  
  if (k==0){
    rho_cross[k+1,1] = C_temp[(k+1),1] / sqrt(2*var_cross)
  }else{
    rho_cross[k+1,1] = C_temp[(k+1),1] / sqrt(var_cross)
  }
  
  ####  6.2 for time rho ####
  if(k != 0) {
    C_time[(k+1),1] = C_time[(k+1),1] / sqrt(N*(T-k)) + sqrt(N*(T-k)) / T* sigma_bar_4
    C_time_nobias[(k+1),1] = C_time_nobias[(k+1),1] / sqrt(N*(T-k))  
    
    rho_time[(k+1),1] = C_time[(k+1),1] / sqrt(var_time)
    rho_time_nobias[(k+1),1] = C_time_nobias[(k+1),1] / sqrt(var_time)
  }
  
  
  
  
  #### 6.3 for joint####
  rho_joint[1,1] <- rho_cross[1,1]^2
  rho_joint_nobias[1,1] <- rho_cross[1,1]^2
  
  for (k in 1:M) {
    rho_joint[(k+1),1] = rho_cross[(k+1),1]^2 + rho_time[(k+1),1]^2 
    rho_joint_nobias[(k+1),1] = rho_cross[(k+1),1]^2 + rho_time_nobias[(k+1),1]^2 
  }
}

#####6.4 test ########
M_time = selectM_norm(rho_time, q1, M)

M_cross =  selectM_norm(rho_cross, q2, M)

#M_joint = max(M_time,M_cross)
# 
Q_all = matrix(0,23,3) # Q_c,auto; Q_c,5; Q_c,10; 
                      # Q_r,auto; Q_r_5, Q_r_10; Q_r0_auto; Q_r0_5; Q_r0_10;
                     # Q_cr_5; Q_cr_10 ;Q_cr0_5; Q_cr0_10
                      # Tn, Wn, Q1, Q2  (M=5) (M=10)
                      # LMbc CD

rownames(Q_all) = c("Q_CD,auto", "Q_CD,5", "Q_CD,10", "Q_CDl,auto", "Q_CDl,5", "Q_CDl,10", "Q_SD,auto", "Q_SD,5", "Q_SD,10",
                    "Q_cr,5", "Q_cr,10", "Q_cr0,5", "Q_cr0,10", 
                    "Q1_5", "Q2_5", "Q1_10", "Q2_10","Tn_5", "Wn_5",  "Tn_10", "Wn_10","LMbc","CD")
colnames(Q_all) = c("Statistics","P-value", "Rej")

test_all = function(value,df_){
  test_all_i = c()
  test_all_i[1] = value
  test_all_i[2] = 1 - pchisq(value, df = df_)
  if(test_all_i[2] <= 0.05){
    test_all_i[3] = 1
  }else{
    test_all_i[3] = 0
  }
  return(test_all_i)
}

##### 6.4.1 Q_time #####
# Q_c auto
Q_all[7,] = test_all(sum(rho_time[2:(M_time+1),1]^2) , M_time) 

# Q_c,5
Q_all[8,] = test_all( sum(rho_time[2:(5+1),1]^2) , 5 ) 

#Q_c
Q_all[9,] = test_all(sum(rho_time[2:(10+1),1]^2) , 10)


#### 6.4.2 Q_cross ####
# Q_r_auto
Q_all[4,] =  test_all(sum(rho_cross[2:(M_cross+1),1]^2) , M_cross)

# Q_r_5
Q_all[5,] = test_all(sum(rho_cross[2:(5+1),1]^2) , 5)


# Q_r_10
Q_all[6,] = test_all(sum(rho_cross[2:(10+1),1]^2) , 10)


# Q_r0_auto
Q_all[1,] =  test_all(sum(rho_cross[1:(M_cross+1),1]^2) , M_cross+1)


# Q_r0_5
Q_all[2,] = test_all(sum(rho_cross[1:(5+1),1]^2) , 6)


# Q_r0_10
Q_all[3,] = test_all(sum(rho_cross[1:(10+1),1]^2) , 11)


#### 6.4.3 Q_joint #####
# Q_cr_5
Q_all[10,] = test_all(sum(rho_joint[2:(5+1),1]) , 10)

# Q_cr_10
Q_all[11,] = test_all(sum(rho_joint[2:(10+1),1]) , 20)


# Q_cr0_5
Q_all[12,] = test_all(sum(rho_joint[1:(5+1),1]) , 11)


# Q_cr0_10
Q_all[13,] = test_all(sum(rho_joint[1:(10+1),1]) , 21)

########BOOT #########

est_U_mat_reserve <- est_U_mat


Q_all_Boot = matrix(0,9,1)

M = max(M_time, M_cross, 11)

B = 500

for(boot in 1:B){
  #random_mat = matrix(rnorm(N*T),N,T)
  #est_U_mat_random = random_mat *  est_U_mat_reserve
  
  random_mat = matrix(sample(c(-1, 1), size = N * T, replace = TRUE, prob = c(0.5, 0.5)), N, T)
  
  est_U_2 <- est_U_mat_reserve^2 * random_mat 
  est_U_4 = est_U_mat_reserve^4 * random_mat^2 
  est_sigma_2 = apply(est_U_2, 1, mean)
  est_tau_4 = apply(est_U_4, 1, mean)
  
  
  sigma_bar_8 = mean((est_tau_4 - est_sigma_2^2)^2)
  #sigma_bar_4_ <- sigma_bar_4
  sigma_bar_4 = mean(est_tau_4 - est_sigma_2^2)
  
  
  var_time = sigma_bar_8# / sigma_bar_4^2
  
  # for cross
  sigma_4 = est_tau_4 - est_sigma_2^2
  #sigma_4 = 2 * est_sigma_2^2
  #omega_4 = omega ^ 2
  var_cross_up = 0
  var_cross_low = 0
  for (j1 in 1:N) {
    index_i1 = index_N[-j1]
    for (j2 in index_i1) {
      var_cross_up = var_cross_up + sigma_4[j1] * sigma_4[j2]
      #sigma_all = sigma_all + omega_4[j1] * omega_4[j2]
      #var_cross_low = var_cross_low + sqrt(sigma_4[j1] * sigma_4[j2])
    }
  }
  var_cross_up = var_cross_up / ( (N-1)*N )
  #var_cross_low = var_cross_low / ( (N-1)*N )
  
  var_cross = var_cross_up #/ (var_cross_low)^2
  
  est_U_sd = matrix(0,N,T)
  for (j in 1:N) {
    est_U_sd[j,] = est_U_2[j,] - est_sigma_2[j] # u_it^2 - est_sigma^2
  }
  
  
  ##### boot1. Compute rho 1~M 
  
  #for (m_n in 1:M_N) {
  #m_n = 2
  
  #M = m_n * M_gap
  
  C_temp = matrix(0, M+1, 1)
  C_time = matrix(0, M+1, 1)
  C_time_nobias = matrix(0, M+1, 1)
  
  rho_time = matrix(0, M+1, 1)
  rho_time_nobias = matrix(0, M+1, 1)
  rho_cross = matrix(0, M+1, 1)
  rho_joint = matrix(0, M+1, 1)
  rho_joint_nobias = matrix(0, M+1, 1)
  
  for (k in 0:M) {
    
    value_autocov = lagk(est_U_sd, k)
    C_temp[(k+1),1] = value_autocov[1]
    C_time[(k+1),1] = value_autocov[2]
    C_time_nobias[(k+1),1] = value_autocov[2]
    
    ####  for cross rho 
    C_temp[k+1, 1] = C_temp[(k+1),1] / sqrt(N*(N-1)*(T-k)) #/ normalize_cross
    
    if (k==0){
      rho_cross[k+1,1] = C_temp[(k+1),1] / sqrt(2*var_cross)
    }else{
      rho_cross[k+1,1] = C_temp[(k+1),1] / sqrt(var_cross)
    }
    
    ####  5.2 for time rho 
    if(k != 0) {
      C_time[(k+1),1] = C_time[(k+1),1] / sqrt(N*(T-k)) + sqrt(N*(T-k)) / T* sigma_bar_4 
      C_time_nobias[(k+1),1] = C_time_nobias[(k+1),1] / sqrt(N*(T-k))  
      
      rho_time[(k+1),1] = C_time[(k+1),1] / sqrt(var_time)
      rho_time_nobias[(k+1),1] = C_time_nobias[(k+1),1] / sqrt(var_time)
    }
    
    #### 5.3 for joint
    rho_joint[1,1] <- rho_cross[1,1]^2
    rho_joint_nobias[1,1] <- rho_cross[1,1]^2
    
    for (k in 1:M) {
      rho_joint[(k+1),1] = rho_cross[(k+1),1]^2 + rho_time[(k+1),1]^2 
      rho_joint_nobias[(k+1),1] = rho_cross[(k+1),1]^2 + rho_time_nobias[(k+1),1]^2 
    }
  }
  
  Q_all_bb = matrix(0,9,1)
  
  Q_all_bb[1] = sum(rho_cross[1:(M_cross+1),1]^2) # Q_cross_withR0
  Q_all_bb[2] = sum(rho_cross[1:6,1]^2) # Q_cross_withR0
  Q_all_bb[3] = sum(rho_cross[1:11,1]^2) # Q_cross_withR0
  
  Q_all_bb[4] = sum(rho_cross[2:(M_cross+1),1]^2) # Q_cross_withR0
  Q_all_bb[5] = sum(rho_cross[2:6,1]^2) # Q_cross_withR0
  Q_all_bb[6] = sum(rho_cross[2:11,1]^2) # Q_cross_withR0
 
  Q_all_bb[7] = sum(rho_time[2:(M_time+1),1]^2) # Q_time
  Q_all_bb[8] = sum(rho_time[2:6,1]^2) # Q_cross_withR0
  Q_all_bb[9] = sum(rho_time[2:11,1]^2) # Q_cross_withR0
    
  
  Q_all_Boot = cbind(Q_all_Boot,Q_all_bb)
}

Q_all_Boot = Q_all_Boot[,-1]

for (ii in 1:9) {
  Q_all[ii,2] = sum(Q_all_Boot[ii,] > Q_all[ii,1]) / B
}



######### 6.4.4 benchmark ###########
test_Q5 = marchtest(t(est_U_mat),5)
ling5 = lingtest(est_U_mat,5)
test_Q10 = marchtest(t(est_U_mat),10)
ling10 = lingtest(est_U_mat,10)
Q_all[14:15,1:2] = test_Q5
Q_all[16:17,1:2] = test_Q10
Q_all[18:19,1:2] = ling5
Q_all[20:21,1:2] = ling10

Q_all[22,1:2] = lagk_baltagi(est_U_mat)

Q_all[23,1:2] = CDtest(est_U_mat)


for (ind in 1:23) {
  if(Q_all[ind,2] <= 0.05){
    Q_all[ind,3] = 1
  }else{
    Q_all[ind,3] = 0
  }
}

#Q_all[22,1] = lagk_baltagi(est_U_mat)


Q_all

write.csv(Q_all, file = "/home/zhpu/Panel CH Test/empricial/Q_all_varspatial_westcoast.csv")



#######7 House Modelling##########

######7.1 Pre for GMM#########
# for Pl_b
P1 = W
P2 = W^2 - diag(W^2)

l_N<-matrix(rep(1,N))

l_T<-matrix(rep(1,T))

I_T<-diag(1,T)

I_N<-diag(1,N)

xmat3 = W %*% xmat1

xmat4 = W %*% xmat2


P1_b = matrix(0,T*N, T*N)
P2_b = matrix(0,T*N, T*N)
for(z in 1:T){
  P2_b[(z*N-N+1):(z*N) , (z*N-N+1):(z*N)] = P2
  P1_b[(z*N-N+1):(z*N) , (z*N-N+1):(z*N)] = P1
}

# for M_bar
Z_bar = matrix(0,T,5)
Z_bar[,1] = apply(ymat,2,mean)
Z_bar[,2] = apply(xmat1,2,mean)
Z_bar[,3] = apply(xmat2,2,mean)
Z_bar[,4] = apply(xmat3,2,mean)
Z_bar[,5] = apply(xmat4,2,mean)

M_bar = I_T - Z_bar%*%ginv(t(Z_bar)%*%Z_bar)%*%t(Z_bar)
M_b = kronecker(M_bar,I_N)

Srou = function(rou_){
  return(I_N - rou_ * W)
}

vectorization = function(mat){
  ncol_ = ncol(mat)
  nrow_ = nrow(mat)
  vec_mat = matrix(0,nrow_ * ncol_,1)
  for (i in 1:ncol_) {
    vec_mat[(i*nrow_-nrow_+1):(i*nrow_),1]  = mat[,i]
  }
  return(vec_mat)
}

matrization = function(vec, nrow_ = N, ncol_ = T){
  vec_mat = matrix(0,nrow_, ncol_)
  for (i in 1:ncol_) {
    vec_mat[,i] = vec[(i*nrow_-nrow_+1):(i*nrow_)]
  }
  return(vec_mat)
}

Y_vec = vectorization(Y)

X_vec = matrix(0,T*N,4)



X_vec[,1] = vectorization(xmat1)
X_vec[,2] = vectorization(xmat2)
X_vec[,3] = vectorization(xmat3)
X_vec[,4] = vectorization(xmat4)

q = 12
h = 2
M_k = 2
M_h = 2

Q = matrix(0,N*T, q)
Q[,1:4] = X_vec
Q[,5] = vectorization(W %*% xmat1)
Q[,6] = vectorization(W %*% xmat2)
Q[,7] = vectorization(W %*% xmat3)
Q[,8] = vectorization(W %*% xmat4)

Q[,9] = vectorization(W^2 %*% xmat1)
Q[,10] = vectorization(W^2 %*% xmat2)
Q[,11] = vectorization(W^2 %*% xmat3)
Q[,12] = vectorization(W^2 %*% xmat4)

P_Q <- M_b %*% Q %*% ginv(t(Q) %*% M_b %*% Q) %*% t(Q) %*% M_b

L = cbind(vectorization(W %*% Y), X_vec)

Q_hat = P_Q %*% Q

Q_hat_tensor = array(0,dim = c(N,T,q))
for (i in 1:N) {
  for (t in 1:T) {
    Q_hat_tensor[i,t,] = Q_hat[(t*N-N+i),]
  }
}

#M_b_Q # M_b * Q 


#####7.2 GMM function#####
GMM = function(par){
  
  rou = par[1]
  beta1 = par[2]
  beta2 = par[3]
  theta1 = par[4]
  theta2 = par[5]
  beta = matrix(par[2:5],4,1)
  delta = matrix(par,5,1)
  
  #####7.2.1 g_NT######
  
  ksi = kronecker(I_T, Srou(rou)) %*% Y_vec - X_vec %*% beta
  
  g_NT = matrix(0,h+q,1)
  
  g_NT[1,1] = t(ksi) %*% M_b %*% P1_b %*% M_b %*% ksi
  g_NT[2,1] = t(ksi) %*% M_b %*% P2_b %*% M_b %*% ksi
  g_NT[(h+1):(h+q),1] = t(ksi) %*% M_b %*% Q
  
  g_NT = g_NT / N /T
  
  #####7.2.2 simga_g #####
  # 1. for sigma_p
  
  e_hat = M_b %*% (Y_vec - L %*% delta)
  e_hat = matrization(e_hat)
  
  gamma_e = array(0,dim = c(M_k+1 , N))
  for (k in 0:M_k) {
    e_hat_tilde = matrix(0,N,k)
    e_hat_t = cbind(e_hat_tilde,e_hat[,1:(T-k)])
    e_hat_final = e_hat %*% t(e_hat_t)
    gamma_e[k+1,] = diag(e_hat_final) / T
  }
  
  zeta = matrix(0,N, N)
  for (z in 1:N) {
    for (j in z:N) {
      if(z == j){
        zeta[z,j] = T * gamma_e[1,z]^2
        for (k in 1:M_k) {
          zeta[z,j] = zeta[z,j] + 2* (T-k)*(1-k/(M_k+1))* gamma_e[k+1,z] * gamma_e[k+1,j]
        }
      }else{
        
        for (k in 1:M_k) {
          zeta[z,j] = zeta[z,j] + 2* (T-k)*(1-k/(M_k+1))* gamma_e[k+1,z] * gamma_e[k+1,j]
        }
        zeta[z,j] <- zeta[j,z] <- zeta[z,j] +  T * gamma_e[1,z] * gamma_e[1,j] 
      }
    }
  }
  
  sigma_p = matrix(0,h,h)
  for (z in 1:N) {
    for (j in 1:N) {
      sigma_p[1, 1] = sigma_p[1, 1] + P1[j,z] * (P1[z,j] + P1[j,z]) * zeta[z,j]
      sigma_p[2, 2] = sigma_p[2, 2] + P2[j,z] * (P2[z,j] + P2[j,z]) * zeta[z,j]
      sigma_p[1, 2] = sigma_p[1, 2] + P1[j,z] * (P2[z,j] + P2[j,z]) * zeta[z,j]
      sigma_p[2, 1] = sigma_p[2, 1] + P2[j,z] * (P1[z,j] + P1[j,z]) * zeta[z,j]
    }
  }
  sigma_p = sigma_p / N /T
  
  ####7.2.3 omega_QMe####
  # Omega_iLPe 
  omega_iLPe_h = array(0,dim = c(M_h+1,N,q,q))
  for (h_ in 0:M_h) {
    for(z in 1:N){
      for (t in (h_ +1):T) {
        omega_iLPe_h[h_+1, z,,] = omega_iLPe_h[h_+1, z,,] + e_hat[z,t] * e_hat[z,t-h_] * matrix(Q_hat_tensor[z,t,],q,1) %*% matrix(Q_hat_tensor[z,t-h_,],1,q)
      }
      omega_iLPe_h[h_+1, z,,] = omega_iLPe_h[h_+1, z,,]/T
    }
  }
  
  omega_LPe = 0
  for (z in 1:N) {
    omega_iLPe = 0
    for (h_ in 1:M_h) {
      omega_iLPe = omega_iLPe + (1-h_ / (M_h + 1))*(omega_iLPe_h[h_+1,z,,] + t(omega_iLPe_h[h_+1,z,,]))
    }
    omega_iLPe = omega_iLPe_h[1,z,,] + omega_iLPe
    omega_LPe = omega_LPe + omega_iLPe
  }
  omega_LPe = omega_LPe / N
  
  sigma_g = matrix(0,h+q,h+q)
  sigma_g[1:h,1:h] = sigma_p
  sigma_g[(h+1):(h+q),(h+1):(h+q)] = omega_LPe
  
  return(t(g_NT) %*% ginv(sigma_g) %*% g_NT)
}

# for mean spatial
est_mean<-solnp(c(0.5,0.3,0.1,0.1,0.1), fun = GMM,
                   LB = c(0.001,-100,-100,-100,-100), 
                   UB = c(0.999,100, 100,100,100),
                control = list(itermax = 5))

est_delta = est_mean$pars



e_hat = M_b %*% (Y_vec - L %*% est_delta)

est_U_mat = matrization(e_hat)

est_U_2 = est_U_mat^2

ymat = est_U_mat^2

Y = ymat



#####7.3 for var spatial ######

#pzh = est_U_mat

est_U_mat = pzh

est_U_2 = est_U_mat^2
est_sigma_2 = apply(est_U_2, 1, mean)
est_omega = est_sigma_2

W_est_U_2 = W %*% est_U_2 
W2_est_U_2 = (W^2) %*% est_U_2 

# log input
log_est_U_2 =log(est_U_2) 
W_log = W %*% log_est_U_2
W2_log = (W^2) %*% log_est_U_2

####7.3.1 T+N modelling####
spatial_GARCH<-function(par){
  
  beta0<-par[1]
  beta1<-par[2]
  beta2<-par[3]
  beta3<-par[4]
  beta4<-par[5]
  beta5<-par[6]
  beta6<-par[7]
  beta7<-par[8]
  beta8<-par[9]
  
  value<-0 
  for(i in 1:N){
    h0<- beta0 + beta1 * W_log[i,1] +  beta2 * xmat1[i,1] + beta3 * xmat2[i,1] + beta4 * xmat3[i,1] + beta5 * xmat4[i,1]
    value<-value + log(dnorm(est_U_2[i,t],sd=sqrt(exp(h0))))
    
    for(t in 2:T){
      
      h0= beta0 + beta1 * W_log[i,t] +  beta2 * xmat1[i,t] + beta3 * xmat2[i,t] + beta4 * xmat3[i,t] + beta5 * xmat4[i,t] + beta6*W_log[i,t-1] + beta7 * log_est_U_2[i,t-1] + beta8 * h0
      
      
      value<-value+log(dnorm(est_U_2[i,t],sd=sqrt(exp(h0))))
    }
    
    
  }
  
  return(-value)
}

est_var<-solnp(c(0.5,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1), fun = spatial_GARCH,
                UB = c(Inf,0.999,Inf,Inf,Inf,Inf,0.999,0.999,0.999), 
                LB = c(-Inf,0.001,-Inf,-Inf,-Inf,-Inf,0.001,0.001,0.001),
                control = list(itermax = 5))

est_var_delta = est_var$pars

W_log_0 = matrix(0,N,1)

W_log_1 = cbind(W_log_0, W_log[,1:T-1])

log_1 = cbind(W_log_0, log_est_U_2[,1:T-1])

est_H = est_var_delta[1] + est_var_delta[2] * W_log + est_var_delta[3] * xmat1 +est_var_delta[4] * xmat2 + est_var_delta[5] * xmat3 + est_var_delta[6] * xmat4 + est_var_delta[7] * W_log_1 + est_var_delta[8] * log_1

est_epsilon_mat = est_U_mat / sqrt(exp(est_H))

est_U_mat = est_epsilon_mat


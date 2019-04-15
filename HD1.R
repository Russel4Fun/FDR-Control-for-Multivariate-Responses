library(Matrix)
library(MASS)
library(parallel)
library(doSNOW)
library(foreach)
set.seed(1234)
p = 500 ; q = 20 ; n = 100; m = 100
B1 = rsparsematrix(nrow = p, ncol = q, density = 0.01, rand.x = rnorm)
B2 = rsparsematrix(nrow = p, ncol = q, density = 0.02, rand.x = rnorm)
B3 = rsparsematrix(nrow = p, ncol = q, density = 0.03, rand.x = rnorm)
B4 = rsparsematrix(nrow = p, ncol = q, density = 0.04, rand.x = rnorm)
Y_choice1 = list()
Y_choice2 = list()
Y_choice3 = list()
Y_choice4 = list()
X_choice = list()
for(k in 1:m){
  X = matrix(0, n, p); E = matrix(0, n, q)
  for(i in 1:p){
    X[,i] = mvrnorm(1, mu = rep(10,n), Sigma = diag(n))
  }
  for(j in 1:q){
    E[,j] = mvrnorm(1, mu = rep(0,n), Sigma = diag(n))
  }
  Y1 = X %*% B1 + E 
  Y2 = X %*% B2 + E 
  Y3 = X %*% B3 + E
  Y4 = X %*% B4 + E 
  Y_choice1[[k]] = Y1
  Y_choice2[[k]] = Y2
  Y_choice3[[k]] = Y3
  Y_choice4[[k]] = Y4
  X_choice[[k]] = X
}

FDR_control <- function(Y_list, X_list, B){
  p = 500 ; q = 20 ; n = 100; m = 100
  FDR_control_short <- function(Y,X,B){
    p = 500 ; q = 20 ; n = 100;m = 100
    HC <- function(P){
      p = dim(P)[1]
      q = dim(P)[2]
      p_value = numeric(p)
      selection = numeric(q)
      for(i in 1:p){
        temp = sort(P[i,])
        for(j in 1:q){
          selection[j] = (j/q - temp[j])/sqrt(temp[j] *(1-temp[j]))
        }
        p_value[i] = sqrt(q)*max(selection)
      }
      return(p_value)
    }
    
    C <- function(x, a){
      ((a^2 - a*sqrt(a^2 + 4*(1-x)*x))/2 + x)/(1+a^2)
    }
    
    C1 <- function(x, a){
      part1 = 1/(1+a^2)
      part2 = a*(1-2*x)/((1+a^2)*sqrt(a^2 + 4*x*(1-x)))
      return(part1-part2)
    }
    
    
    D_pvalue <- function(p_value, q){
      p = length(p_value)
      result = numeric(p)
      n = q
      for(i in 1:p){
        sum = 0
        b = p_value[i]
        for(k in 1:n){
          c = C(k/n, b/sqrt(n))
          c1 = C1(k/n, b/sqrt(n))
          first_part = dbeta(c, shape1 = k, shape2 = n+1-k)
          second_part = c/k
          third_part = 1 - (1-k/n)*c1/(1-c)
          sum = sum + first_part*second_part*third_part
        }
        result[i] = sum
      }
      return(result)
    }
    
    p_value <- function(modelobject){
      if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
      f <- summary(modelobject)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      return(p)  
    }
    
    B_hat = matrix(0,p,q)
    P = matrix(0,p,q)
    for(i in 1:p){
      for(j in 1:q){
        fit = lm(Y[,j] ~ X[,i])
        P[i,j] = p_value(fit)
        attributes(fit$coefficients) = NULL
        B_hat[i,j] = fit$coefficients[2]
      }
    }
    
    P_raw = D_pvalue(HC(P), 20)
    P_adjust = p.adjust(P_raw, method = 'hochberg')
    H0_Rejected = which(P_adjust < 0.05)
    B_result = apply(B, 1, max)
    H0_True = which(B_result == 0)
    Rejected_total = length(H0_Rejected)
    V = 0
    for(t in 1:length(Rejected)){
      if(H0_Rejected[t] %in% H0_True){V = V+1}
    }
    return(V/Rejected_total)
  }
  Fdp = numeric(100)
  for(t in 1:100){
    print(t)
    Fdp[t] = FDR_control_short(Y_list[[t]], X_list[[t]], B)
  }
  return(Fdp)
}
  
F1 = FDR_control(Y_choice1, X_choice, B1)
F2 = FDR_control(Y_choice2, X_choice, B2)
F3 = FDR_control(Y_choice3, X_choice, B3)
F4 = FDR_control(Y_choice4, X_choice, B4)
  

X1 = X_choice[[1]] 
Y1 = Y_choice1[[1]] 
B_hat = matrix(0,p,q)
P = matrix(0,p,q)
for(i in 1:p){
  for(j in 1:q){
    fit = lm(Y1[,j] ~ X1[,i])
    P[i,j] = p_value(fit)
    attributes(fit$coefficients) = NULL
    B_hat[i,j] = fit$coefficients[2]
  }
}

P_hat = numeric(500)
for(i in 1:500){
  P_hat[i] = test.hc(P[i,], M=diag(20),k0=1,k1=20)$pvalue
}
rejected = which(P_hat < 0.05)
B1_adjust = apply(B1,1,max)
H0_true = which(B1_adjust>0)

FDR_control_pkg <- function(Y_list, X_list, B){
  p = 500 ; q = 20 ; n = 100; m = 100
  FDR_control_pkg_short <- function(Y, X, B){
    p = 500 ; q = 20 ; n = 100; m = 100
    P = matrix(0,p,q)
    for(i in 1:p){
      for(j in 1:q){
        fit = lm(Y[,j] ~ X[,i])
        P[i,j] = p_value(fit)
        attributes(fit$coefficients) = NULL
        B_hat[i,j] = fit$coefficients[2]
      }
    }
    P_hat = numeric(500)
    for(i in 1:500){
      P_hat[i] = test.hc(P[i,], M=diag(20),k0=1,k1=10)$pvalue
    }
    P_adjust = p.adjust(P_hat, method = 'hochberg')
    H0_Rejected = which(P_adjust < 0.05)
    B_result = apply(B, 1, max)
    H0_True = which(B_result == 0)  
    Rejected_total = length(H0_Rejected)
    V = 0
    for(t in 1:length(Rejected)){
      if(H0_Rejected[t] %in% H0_True){V = V+1}
    }
    return(V/Rejected_total)
  }
  Fdp = numeric(100)
  for(t in 1:100){
    print(t)
    Fdp[t] = FDR_control_pkg_short(Y_list[[t]], X_list[[t]], B)
  }
  return(Fdp)
}

FDR_control_pkg_bj <- function(Y_list, X_list, B){
  p = 500 ; q = 20 ; n = 100; m = 100
  FDR_control_pkg_short_bj <- function(Y, X, B){
    p = 500 ; q = 20 ; n = 100; m = 100
    P = matrix(0,p,q)
    for(i in 1:p){
      for(j in 1:q){
        fit = lm(Y[,j] ~ X[,i])
        P[i,j] = p_value(fit)
        attributes(fit$coefficients) = NULL
        B_hat[i,j] = fit$coefficients[2]
      }
    }
    P_hat = numeric(500)
    for(i in 1:500){
      P_hat[i] = test.bj(P[i,], M=diag(20),k0=1,k1=20)$pvalue
    }
    P_adjust = p.adjust(P_hat, method = 'hochberg')
    H0_Rejected = which(P_adjust < 0.05)
    B_result = apply(B, 1, max)
    H0_True = which(B_result == 0)  
    Rejected_total = length(H0_Rejected)
    V = 0
    for(t in 1:length(Rejected)){
      if(H0_Rejected[t] %in% H0_True){V = V+1}
    }
    return(V/Rejected_total)
  }
  Fdp = numeric(100)
  for(t in 1:100){
    print(t)
    Fdp[t] = FDR_control_pkg_short_bj(Y_list[[t]], X_list[[t]], B)
  }
  return(Fdp)
}

F1_hc = FDR_control_pkg(Y_choice1, X_choice, B1)
F2_hc = FDR_control_pkg(Y_choice2, X_choice, B2)
F3_hc = FDR_control_pkg(Y_choice3, X_choice, B3)
F4_hc = FDR_control_pkg(Y_choice4, X_choice, B4)
F1_bj = FDR_control_pkg_bj(Y_choice1, X_choice, B1)
F2_bj = FDR_control_pkg_bj(Y_choice2, X_choice, B2)
F3_bj = FDR_control_pkg_bj(Y_choice3, X_choice, B3)
F4_bj = FDR_control_pkg_bj(Y_choice4, X_choice, B4)

data1 = as.data.frame(rbind(cbind(F1_hc,rep(1,100)),cbind(F1_bj,rep(2,100))))
data2 = as.data.frame(rbind(cbind(F2_hc,rep(1,100)),cbind(F2_bj,rep(2,100))))
data3 = as.data.frame(rbind(cbind(F3_hc,rep(1,100)),cbind(F3_bj,rep(2,100))))
data4 = as.data.frame(rbind(cbind(F4_hc,rep(1,100)),cbind(F4_bj,rep(2,100))))
colnames(data1) = c('FDR','class')
colnames(data2) = c('FDR','class')
colnames(data3) = c('FDR','class')
colnames(data4) = c('FDR','class')

ggplot(data1, aes(x = factor(class), y = FDR, fill = factor(class))) + 
  geom_boxplot() + xlab('Different Method') +
  scale_fill_discrete(name="Method",
                      labels = c('Higher Criticism','Berk-Jones')) 


ggplot(data2, aes(x = factor(class), y = FDR, fill = factor(class))) + 
  geom_boxplot() + xlab('Different Method') +
  scale_fill_discrete(name="Method",
                      labels = c('Higher Criticism','Berk-Jones'))

ggplot(data3, aes(x = factor(class), y = FDR, fill = factor(class))) + 
  geom_boxplot() + xlab('Different Method') +
  scale_fill_discrete(name="Method",
                      labels = c('Higher Criticism','Berk-Jones')) 

ggplot(data4, aes(x = factor(class), y = FDR, fill = factor(class))) + 
  geom_boxplot() + xlab('Different Method') +
  scale_fill_discrete(name="Method",
                      labels = c('Higher Criticism','Berk-Jones'))



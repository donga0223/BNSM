args <- commandArgs(trailingOnly = TRUE)
myrep <- args[1]
k <- as.numeric(myrep)

source("R_function/thetaprior_function_forpaper.R")

library(sna)
m <-  6  ## number of relation 
V <- 50  ## numoer of node


n.group <- 3
myemp <- c(1,9,1,9,9,1)
myepp <- c(1,9,9,1,1,9)

true_np <- 0.05
mymodel_np <- "obs" 

g10 <- obs.g10 <- list()
trueep10 <- trueem10 <- matrix(NA, nrow = 1, ncol = V)



set.seed(1200+10*k+1)
tmp_ind <- sample(1:V, V*0.1)
tmp1 <- sample(1:length(tmp_ind), round(length(tmp_ind)/2))
d_ind <- tmp_ind[tmp1]
d_ind2 <- tmp_ind[-tmp1]
mis <- sample(V,5)

tmp <- gen.emep.dep.net(V, m, deprho = -0.5, similarity = 0.7, netprob = true_np, d_ind = d_ind, d_ind2 = d_ind2, em = myemp, ep = myepp, myseed = (1200+k))
tmp$obs.g[,mis,] <- NA

g10 <- tmp$g
obs.g10 <- tmp$obs.g
trueem10 <- tmp$trueem
trueep10 <- tmp$trueep


#model_np <- netmean(tmp$obs.g)
if(mymodel_np == "true"){
  model_np <- rep(true_np,m)
}else if(mymodel_np == "obs"){
  model_np <- netmean(tmp$obs.g)
}



a1.em.list <- c(3,1,3,1)
b1.em.list <- c(5,9,5,9)

a1.theta.list <- c(3,3,1,1)
b1.theta.list <- c(5,5,9,9)


theta1.em.density <- theta2.em.density <- theta3.em.density <- theta4.em.density <- theta5.em.density <- list()
theta1.ep.density <- theta2.ep.density <- theta3.ep.density <- theta4.ep.density <- theta5.ep.density <- list()
theta1.em.CI <- theta2.em.CI <- theta3.em.CI <- theta4.em.CI <- theta5.em.CI <- list()
theta1.ep.CI <- theta2.ep.CI <- theta3.ep.CI <- theta4.ep.CI <- theta5.ep.CI <- list()
theta1.b01 <- theta2.b01 <- theta3.b01 <- theta4.b01 <- theta5.b01 <- list()
theta1.pred.y <- theta2.pred.y <- theta3.pred.y <- theta4.pred.y <- theta5.pred.y <- list()
theta1.np <- theta2.np <- theta3.np <- theta4.np <- theta5.np <- list()



for(i in 1:length(a1.em.list)){
  
  a1.em = a1.em.list[i]; b1.em = b1.em.list[i]; a2.em = 1; b2.em = 5
  a1.ep = 1; b1.ep = 1; a2.ep = 1; b2.ep = 30 ### uniform prior
  a1.theta = a1.theta.list[i]; b1.theta = b1.theta.list[i]; a2.theta = 1; b2.theta = 5 
  
  alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
  alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
  alpha.acc.theta = 0; beta.acc.theta = 0; total.alpha.theta = 0; total.beta.theta = 0
  
  assign(paste("theta",i,sep=""), bnsm.hyper.theta.hyper(dat = obs.g10, emsd = list(ema = 10, emb = 10, alpha.sd = .5, beta.sd = .5)
                                                         , epsd = list(epa = 10, epb = 10, alpha.sd = .5, beta.sd = .5)
                                                         , thetasd = list(thetaa = 3, thetab = 9, alpha.sd = .5, beta.sd = .5)
                                                         , diag = FALSE, mode = "graph", model.checking = TRUE
                                                         , reps = 3, draws = 1500, tinning = 100, burntime = 15000
                                                         , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = FALSE))
  
  save.image(paste("simulation/BNSM_simulation_missing_rep_", k, ".RData", sep = ""))
  
  
  
  
  
  
  if(i==1){
    theta1.em.density <- density(theta1$em)
    theta1.ep.density <- density(theta1$ep)
    theta1.em.CI <- apply(theta1$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta1.ep.CI <- apply(theta1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta1.b01 <- apply(theta1$net, c(3,4,5), mean)
    theta1.pred.y <- pred.y.mean(theta1$pred.y, iter = 1500, reps = 3, m = 6)
    theta1.np <- theta1$np
    rm(theta1)
  }else if(i==2){
    theta2.em.density <- density(theta2$em)
    theta2.ep.density <- density(theta2$ep)
    theta2.em.CI <- apply(theta2$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta2.ep.CI <- apply(theta2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta2.b01 <- apply(theta2$net, c(3,4,5), mean)
    theta2.pred.y <- pred.y.mean(theta2$pred.y, iter = 1500, reps = 3, m = 6)
    theta2.np <- theta2$np
    rm(theta2)
  }else if(i==3){
    theta3.em.density <- density(theta3$em)
    theta3.ep.density <- density(theta3$ep)
    theta3.em.CI <- apply(theta3$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta3.ep.CI <- apply(theta3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta3.b01 <- apply(theta3$net, c(3,4,5), mean)
    theta3.pred.y <- pred.y.mean(theta3$pred.y, iter = 1500, reps = 3, m = 6)
    theta3.np <- theta3$np
    rm(theta3)
  }else if(i==4){
    theta4.em.density <- density(theta4$em)
    theta4.ep.density <- density(theta4$ep)
    theta4.em.CI <- apply(theta4$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta4.ep.CI <- apply(theta4$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta4.b01 <- apply(theta4$net, c(3,4,5), mean)
    theta4.pred.y <- pred.y.mean(theta4$pred.y, iter = 1500, reps = 3, m = 6)
    theta4.np <- theta4$np
    rm(theta4)
  }
  
  print(k)
  save.image(paste("simulation/BNSM_simulation_missing_rep_", k, ".RData", sep = ""))
}




args <- commandArgs(trailingOnly = TRUE)
myrep <- args[1]
k <- as.numeric(myrep)
school <- args[2]
print(school)
print(k)

source("thetaprior_function_forpaper.R")

library(sna)
m <-  6
n.group <- 3

n_draws = 1500
n_tinning = 100
n_burntime = 15000

if(school == "ARTS02"){
  V <- 55
  myemp <- c(2,5,10,1,10,1)
  myepp <- c(1,50,1,10,1,10)
  group.p <- 0.1
  true_np <- c(0.05,0.09, 0.01, 0.025, 0.04, 0.13)
  mis <- sample(V, 5)
}else if(school == "ARTS03"){
  V <- 55
  myemp <- c(2,5,5,2,5,2)
  myepp <- c(1,50,1,10,1,10)
  group.p <- 0.1
  true_np <- c(0.04,0.15, 0.01, 0.05, 0.08, 0.20)
  mis <- sample(V, 8)
}else if(school == "PREP02"){
  V <- 105
  myemp <- c(1,1,5,2,5,2)
  myepp <- c(1,100,1,100,1,9)
  group.p <- 0.5
  true_np <- c(0.05,0.15, 0.005, 0.1, 0.12, 0.16)
  mis <- sample(V, 3)
}else if(school == "TECH03"){
  V <- 100
  myemp <- c(5,5,3,9,3,9)
  myepp <- c(1,100,1,100,1,1)
  group.p <- 0.2
  true_np <- c(0.025,0.08, 0.005, 0.01, 0.06, 0.17) 
  mis <- sample(V, 10)
}


g10 <- obs.g10 <- list()
trueep10 <- trueem10 <- matrix(NA, nrow = 1, ncol = V)

set.seed(1200+k)
tmp_ind <- sample(1:V, V*group.p)
tmp1 <- sample(1:length(tmp_ind), round(length(tmp_ind)/2))
d_ind <- tmp_ind[tmp1]
d_ind2 <- tmp_ind[-tmp1]

tmp <- gen.emep.dep.net(V, m, deprho = -0.5, similarity = 0.7, netprob = true_np, d_ind = d_ind, d_ind2 = d_ind2, em = myemp, ep = myepp, myseed = (1200+k))

g10 <- tmp$g
obs.g10 <- tmp$obs.g
obs.g10[,mis,] <- NA
trueem10 <- tmp$trueem
trueep10 <- tmp$trueep

a1.em.list <- c(3,1,1)
b1.em.list <- c(5,9,1)
b2.em.list <- c(5,5,30)

a1.theta.list <- c(3,1,1)
b1.theta.list <- c(5,9,1)
b2.theta.list <- c(5,5,30)

theta1.em.density <- theta2.em.density <- theta3.em.density <- theta4.em.density <- theta5.em.density <- list()
theta1.ep.density <- theta2.ep.density <- theta3.ep.density <- theta4.ep.density <- theta5.ep.density <- list()
theta1.em.CI <- theta2.em.CI <- theta3.em.CI <- theta4.em.CI <- theta5.em.CI <- list()
theta1.ep.CI <- theta2.ep.CI <- theta3.ep.CI <- theta4.ep.CI <- theta5.ep.CI <- list()
theta1.b01 <- theta2.b01 <- theta3.b01 <- theta4.b01 <- theta5.b01 <- list()
theta1.pred.y <- theta2.pred.y <- theta3.pred.y <- theta4.pred.y <- theta5.pred.y <- list()
theta1.np <- theta2.np <- theta3.np <- theta4.np <- theta5.np <- list()


one_theta_all.em.density <- one_theta1.em.density <- one_theta2.em.density <- one_theta3.em.density <- one_theta4.em.density <- one_theta5.em.density <- one_theta6.em.density <- list()
one_theta_all.ep.density <- one_theta1.ep.density <- one_theta2.ep.density <- one_theta3.ep.density <- one_theta4.ep.density <- one_theta5.ep.density <- one_theta6.ep.density <- list()
one_theta_all.em.CI <- one_theta1.em.CI <- one_theta2.em.CI <- one_theta3.em.CI <- one_theta4.em.CI <- one_theta5.em.CI <- one_theta6.em.CI <- list()
one_theta_all.ep.CI <- one_theta1.ep.CI <- one_theta2.ep.CI <- one_theta3.ep.CI <- one_theta4.ep.CI <- one_theta5.ep.CI <- one_theta6.ep.CI <- list()
one_theta_all.b01 <- one_theta1.b01 <- one_theta2.b01 <- one_theta3.b01 <- one_theta4.b01 <- one_theta5.b01 <- one_theta6.b01 <- list()
one_theta_all.pred.y <- one_theta1.pred.y <- one_theta2.pred.y <- one_theta3.pred.y <- one_theta4.pred.y <- one_theta5.pred.y <- one_theta6.pred.y <- list()
one_theta_all.np <- one_theta1.np <- one_theta2.np <- one_theta3.np <- one_theta4.np <- one_theta5.np <- one_theta6.np <- list()


one_theta <- bnsm.one.theta(dat = obs.g10
                                   , emp = c(1,11), epp = c(1,11)
                                   , thetaalpha = 0.5, thetabeta = 0.5
                                   , diag = FALSE, mode = "graph", model.checking = TRUE
                                   , reps = 3, draws = n_draws, tinning = n_tinning, burntime = n_burntime
                                   , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE)

print("done one_theta")
one_theta_all.em.density <- density(one_theta$em)
one_theta_all.ep.density <- density(one_theta$ep)
one_theta_all.em.CI <- apply(one_theta$em, 2, quantile, prob = c(0.025,0.975, 0.5))
one_theta_all.ep.CI <- apply(one_theta$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
one_theta_all.b01 <- apply(one_theta$net, c(3,4,5), mean)
one_theta_all.pred.y <- pred.y.mean(one_theta$pred.y, iter = n_draws, reps = 3, m = m)
one_theta_all.np <- one_theta$np
rm(one_theta)

for(i in 1:6){
  print(i)
  assign(paste("one_theta",i,sep=""), bnsm.one.theta(dat = obs.g10[i,,]
                                   , emp = c(1,11), epp = c(1,11)
                                   , thetaalpha = 0.5, thetabeta = 0.5
                                   , diag = FALSE, mode = "graph", model.checking = TRUE
                                   , reps = 3, draws = n_draws, tinning = n_tinning, burntime = n_burntime
                                   , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE))

  if(i==1){
    one_theta1.em.density <- density(one_theta1$em)
    one_theta1.ep.density <- density(one_theta1$ep)
    one_theta1.em.CI <- apply(one_theta1$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta1.ep.CI <- apply(one_theta1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta1.b01 <- apply(one_theta1$net, c(3,4,5), mean)
    one_theta1.pred.y <- pred.y.mean(one_theta1$pred.y, iter = n_draws, reps = 3, m = 1)
    one_theta1.np <- one_theta1$np
    rm(one_theta1)
  }else if(i==2){
    one_theta2.em.density <- density(one_theta2$em)
    one_theta2.ep.density <- density(one_theta2$ep)
    one_theta2.em.CI <- apply(one_theta2$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta2.ep.CI <- apply(one_theta2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta2.b01 <- apply(one_theta2$net, c(3,4,5), mean)
    one_theta2.pred.y <- pred.y.mean(one_theta2$pred.y, iter = n_draws, reps = 3, m = 1)
    one_theta2.np <- one_theta2$np
    rm(one_theta2)
  }else if(i==3){
    one_theta3.em.density <- density(one_theta3$em)
    one_theta3.ep.density <- density(one_theta3$ep)
    one_theta3.em.CI <- apply(one_theta3$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta3.ep.CI <- apply(one_theta3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta3.b01 <- apply(one_theta3$net, c(3,4,5), mean)
    one_theta3.pred.y <- pred.y.mean(one_theta3$pred.y, iter = n_draws, reps = 3, m = 1)
    one_theta3.np <- one_theta3$np
    rm(one_theta3)
  }else if(i==4){
    one_theta4.em.density <- density(one_theta4$em)
    one_theta4.ep.density <- density(one_theta4$ep)
    one_theta4.em.CI <- apply(one_theta4$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta4.ep.CI <- apply(one_theta4$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta4.b01 <- apply(one_theta4$net, c(3,4,5), mean)
    one_theta4.pred.y <- pred.y.mean(one_theta4$pred.y, iter = n_draws, reps = 3, m = 1)
    one_theta4.np <- one_theta4$np
    rm(one_theta4)
  }else if(i==5){
    one_theta5.em.density <- density(one_theta5$em)
    one_theta5.ep.density <- density(one_theta5$ep)
    one_theta5.em.CI <- apply(one_theta5$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta5.ep.CI <- apply(one_theta5$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta5.b01 <- apply(one_theta5$net, c(3,4,5), mean)
    one_theta5.pred.y <- pred.y.mean(one_theta5$pred.y, iter = n_draws, reps = 3, m = 1)
    one_theta5.np <- one_theta5$np
    rm(one_theta5)
  }else if(i==6){
    one_theta6.em.density <- density(one_theta6$em)
    one_theta6.ep.density <- density(one_theta6$ep)
    one_theta6.em.CI <- apply(one_theta6$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta6.ep.CI <- apply(one_theta6$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    one_theta6.b01 <- apply(one_theta6$net, c(3,4,5), mean)
    one_theta6.pred.y <- pred.y.mean(one_theta6$pred.y, iter = n_draws, reps = 3, m = 1)
    one_theta6.np <- one_theta6$np
    rm(one_theta6)
  }
}
save.image(file = paste("data_matched_sim/BNSM_fit/BNSM_data_matched_simulation_", school,"_", k, ".RData", sep=""))

  
#save.image(file = paste("data_matched_sim/BNSM_fit/BNSM_data_matched_simulation_", school,"_", k, ".RData", sep=""))



for(i in 1:length(b1.em.list)){
  print(i)
  a1.em = a1.em.list[i]; b1.em = b1.em.list[i]; a2.em = 1; b2.em = b2.em.list[i]
  a1.ep = 1; b1.ep = 1; a2.ep = 1; b2.ep = 30 ### uniform prior
  a1.theta = a1.theta.list[i]; b1.theta = b1.theta.list[i]; a2.theta = 1; b2.theta = b2.theta.list[i] 
  
  alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
  alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
  alpha.acc.theta = 0; beta.acc.theta = 0; total.alpha.theta = 0; total.beta.theta = 0
  
  assign(paste("theta",i,sep=""), bnsm.hyper.theta.hyper(dat = obs.g10, emsd = list(ema = 10, emb = 10, alpha.sd = .5, beta.sd = .5)
                                                         , epsd = list(epa = 10, epb = 10, alpha.sd = .5, beta.sd = .5)
                                                         , thetasd = list(thetaa = 3, thetab = 9, alpha.sd = .5, beta.sd = .5)
                                                         , diag = FALSE, mode = "graph", model.checking = TRUE
                                                         , reps = 3, draws = n_draws, tinning = n_tinning, burntime = n_burntime
                                                         , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = FALSE))
  
  
  if(i==1){
    theta1.em.density <- density(theta1$em)
    theta1.ep.density <- density(theta1$ep)
    theta1.em.CI <- apply(theta1$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta1.ep.CI <- apply(theta1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta1.b01 <- apply(theta1$net, c(3,4,5), mean)
    theta1.pred.y <- pred.y.mean(theta1$pred.y, iter = n_draws, reps = 3, m = 6)
    theta1.np <- theta1$np
    rm(theta1)
  }else if(i==2){
    theta2.em.density <- density(theta2$em)
    theta2.ep.density <- density(theta2$ep)
    theta2.em.CI <- apply(theta2$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta2.ep.CI <- apply(theta2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta2.b01 <- apply(theta2$net, c(3,4,5), mean)
    theta2.pred.y <- pred.y.mean(theta2$pred.y, iter = n_draws, reps = 3, m = 6)
    theta2.np <- theta2$np
    rm(theta2)
  }else if(i==3){
    theta3.em.density <- density(theta3$em)
    theta3.ep.density <- density(theta3$ep)
    theta3.em.CI <- apply(theta3$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta3.ep.CI <- apply(theta3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    theta3.b01 <- apply(theta3$net, c(3,4,5), mean)
    theta3.pred.y <- pred.y.mean(theta3$pred.y, iter = n_draws, reps = 3, m = 6)
    theta3.np <- theta3$np
    rm(theta3)
  }
}
save.image(file = paste("data_matched_sim/BNSM_fit/BNSM_data_matched_simulation_", school,"_", k, ".RData", sep=""))



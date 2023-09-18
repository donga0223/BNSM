args <- commandArgs(trailingOnly = TRUE)
myrep <- args[1]
k <- as.numeric(myrep)

source("thetaprior_function_forpaper.R")

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





if(2==3){
  rm(list=ls())
  
  load("/Users/dongahkim/Downloads/paper.simulation.d.3.0.05.theta.hyper.missing.RData")
  setwd("/Users/dongahkim/OneDrive/BNSM/paper simulation/simulation")
  
  myrep <- 1
  mycol <- c("red", "gray", "pink", "yellow", "green")
  mylegend <- c("True" , "M1", "M2", "M3", "M4")
  
  em.reorder <- order(trueem10[myrep,])
  ep.reorder <- order(trueep10[myrep,])
  trueem1 <- trueem10[myrep,em.reorder]
  trueep1 <- trueep10[myrep,ep.reorder]
  
  theta1.em.CI[[2]] <- theta1.em.CI[[1]][,em.reorder]
  theta2.em.CI[[2]] <- theta2.em.CI[[1]][,em.reorder]
  theta3.em.CI[[2]] <- theta3.em.CI[[1]][,em.reorder]
  theta4.em.CI[[2]] <- theta4.em.CI[[1]][,em.reorder]
  
  theta1.ep.CI[[2]] <- theta1.ep.CI[[1]][,ep.reorder]
  theta2.ep.CI[[2]] <- theta2.ep.CI[[1]][,ep.reorder]
  theta3.ep.CI[[2]] <- theta3.ep.CI[[1]][,ep.reorder]
  theta4.ep.CI[[2]] <- theta4.ep.CI[[1]][,ep.reorder]
  
  pdf("paper_simulation_missing.pdf")
  par(mar = c(4.5, 5, 3.5, 1.5))
  
  plot(c(0, 1 + V), c(0,1), type = "n", ylab = "e- Credible Interval",xlab = "node number"
       , cex.lab = 2, cex.main = 2, cex.axis = 2, main = "Posterior for each node")
  for (i in 1:V) lines(rep(i-0.15, 2), theta1.em.CI[[2]][1:2, i], lwd = 3, col = mycol[2])
  for (i in 1:V) lines(rep(i-0.05, 2), theta2.em.CI[[2]][1:2, i], lwd = 3, col = mycol[3])
  for (i in 1:V) lines(rep(i+0.05, 2), theta3.em.CI[[2]][1:2, i], lwd = 3, col = mycol[4])
  for (i in 1:V) lines(rep(i+0.15, 2), theta4.em.CI[[2]][1:2, i], lwd = 3, col = mycol[5])
  
  for (i in 1:V)  lines(c(i-0.2, i+0.2), rep(trueem1[i],2), col = mycol[1], lwd = 2)
  legend("topleft", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)
  
  
  plot(c(0, 1 + V), c(0,1), type = "n", ylab = "e+ Credible Interval",xlab = "node number"
       , cex.lab = 2, cex.main = 2, cex.axis = 2, main = "Posterior for each node")
  for (i in 1:V) lines(rep(i-0.15, 2), theta1.ep.CI[[2]][1:2, i], lwd = 3, col = mycol[2])
  for (i in 1:V) lines(rep(i-0.05, 2), theta2.ep.CI[[2]][1:2, i], lwd = 3, col = mycol[3])
  for (i in 1:V) lines(rep(i+0.05, 2), theta3.ep.CI[[2]][1:2, i], lwd = 3, col = mycol[4])
  for (i in 1:V) lines(rep(i+0.15, 2), theta4.ep.CI[[2]][1:2, i], lwd = 3, col = mycol[5])
  for (i in 1:V)  lines(c(i-0.2, i+0.2), rep(trueep1[i],2), col = "red", lwd = 2)
  legend("topleft", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)
  
  
  xs <- seq(0,1,by = 0.001)
  
  emy <- 0.95*dbeta(xs, 1,9)+0.05*dbeta(xs, 9, 1)
  epy <- 0.95*dbeta(xs, 1,9)+0.05*dbeta(xs, 9, 1)
  plot(xs, emy, type = "l", lwd = 3, col = 2, xlab = "e-", ylab = "Density", main = "Posterior for all observations in 50 reps", cex.main = 2, cex.lab = 2, cex.axis = 2)
  
  #plot(density(trueem10[1,]), col = mycol[1], lwd = 2, xlab = "e-",cex.main = 2, cex.axis = 2, cex.lab = 2, "Posterior for all observation")
  lines(theta1.em.density[[1]], col = mycol[2], lwd = 3)
  lines(theta2.em.density[[1]], col = mycol[3], lwd = 3)
  lines(theta3.em.density[[1]], col = mycol[4], lwd = 3)
  lines(theta4.em.density[[1]], col = mycol[5], lwd = 3)
  legend("topright", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)
  
  plot(xs, epy, type = "l", lwd = 3, col = 2, xlab = "e+", ylab = "Density", main = "Posterior for all observations in 50 reps", cex.main = 2, cex.lab = 2, cex.axis = 2)
  #plot(density(trueep10[1,]), col = mycol[1], lwd = 2, xlab = "e+",cex.main = 2, cex.lab = 2, cex.axis = 2, "Posterior for all observation", ylim = c(0,8))
  lines(theta1.ep.density[[1]], col = mycol[2], lwd = 3)
  lines(theta2.ep.density[[1]], col = mycol[3], lwd = 3)
  lines(theta3.ep.density[[1]], col = mycol[4], lwd = 3)
  lines(theta4.ep.density[[1]], col = mycol[5], lwd = 3)
  legend("topright", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)
  
  
  m <- dim(theta1.np[[1]])[2]
  
  netnp <- array(dim=c(4500, m*4))
  for(i in 1:6){
    for(k in 1:4){
      #print(k+(i-1)*3)
      netnp[,(k+(i-1)*4)] <- get(paste("theta", k, ".np", sep=""))[[1]][,i]
    }
  }
  mn <- 4
  true_np <- netmean(g10[[1]])
  boxplot(netnp, col = mycol[2:5], main = "Posterior for Theta", cex.main = 2, xaxt = "n")
  for (i in 1:V)  lines(c(4*(i-1)+0.5, 4*i+0.5), rep(true_np[i],2), col = 2, lwd = 2)
  abline(v = (1:m)*mn+0.5, lty = 2)
  
  axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  #axis(side=1, at = 1:(mn*m), mgp=c(0,1.5,0), labels = rep(paste0("m",1:mn),6), tick = 0, cex.axis = 0.6)
  legend("topleft", legend = "True mean density", lwd = 3, cex = 1.5, col = 2)
  #legend("topright", legend = mylegend[-1], col = mycol[-1], lwd = 3, cex = 2)
  
  
  union.intersection.net <- function(observed){
    union.net <- intersection.net <- observed
    aa <- which(is.na(observed[1,,1]))
    for(i in 1:dim(observed)[1]){
      observed1 <- ifelse(is.na(observed[i,,]),0,observed[i,,])
      #union.net[i,,] <- ifelse(observed1[i,,]+t(observed1[i,,]) > 0, 1, 0)
      union.net[i,,] <- ifelse(observed1+t(observed1) > 0, 1, 0)
      union.net[i,aa,] <- ifelse( union.net[i,aa,]==1,1,NA)
      
      intersection.net[i,,] <- ifelse(observed[i,,]+t(observed[i,,]) > 1, 1, 0)
    }
    return(list(union.net = union.net, intersection.net = intersection.net))
  }
  
  union.net <- union.intersection.net(obs.g10[[myrep]])$union.net
  intersection.net <- union.intersection.net(obs.g10[[myrep]])$intersection.net
  
  
  roc.table <- function(bnet, g){
    bnet <- diag.remove(bnet)
    g <- diag.remove(g)
    TP <- sum(ifelse(bnet == 1 & g == 1, 1,0), na.rm = TRUE)
    FP <- sum(ifelse(bnet == 1 & g == 0, 1,0), na.rm = TRUE)
    FN <- sum(ifelse(bnet == 0 & g == 1, 1,0), na.rm = TRUE)
    TN <- sum(ifelse(bnet == 0 & g == 0, 1,0), na.rm = TRUE)
    TPR <- TP/(TP+FN)
    FPR <- FP/(FP+TN)
    res <- c(FPR, TPR)
    names(res) <- c("FPR", "TPR")
    return(res)
  }
  
  
  
  
  uni.xy <- roc.table(union.net, g10[[myrep]])
  inter.xy <- roc.table(intersection.net, g10[[myrep]])
  
  
  library(ROCR)
  library(sna)
  
  roc.theta1.b01 <- c(theta1.b01[[1]])
  roc.theta2.b01 <- c(theta2.b01[[1]])
  roc.theta3.b01 <- c(theta3.b01[[1]])
  roc.theta4.b01 <- c(theta4.b01[[1]])
  roc.g10 <- c(g10[[myrep]])
  
  diagna <- which(is.na(roc.theta1.b01))
  roc.theta1.b01 <- roc.theta1.b01[-diagna] 
  roc.theta2.b01 <- roc.theta2.b01[-diagna] 
  roc.theta3.b01 <- roc.theta3.b01[-diagna] 
  roc.theta4.b01 <- roc.theta4.b01[-diagna] 
  roc.g10 <- roc.g10[-diagna]
  
  theta1.pred <- prediction(roc.theta1.b01, roc.g10)
  theta2.pred <- prediction(roc.theta2.b01, roc.g10)
  theta3.pred <- prediction(roc.theta3.b01, roc.g10)
  theta4.pred <- prediction(roc.theta4.b01, roc.g10)
  
  theta1.perf <- performance(theta1.pred,"tpr","fpr")
  theta2.perf <- performance(theta2.pred,"tpr","fpr")
  theta3.perf <- performance(theta3.pred,"tpr","fpr")
  theta4.perf <- performance(theta4.pred,"tpr","fpr")
  
  roc.threshold <- function(theta.b01, g, threshold){
    iter <- length(theta.b01)
    res <- matrix(NA, ncol = 2, nrow = iter)
    for(i in 1:iter){
      aaa <- mean(threshold[[i]],na.rm = TRUE)
      est.b01 <- (theta.b01[[i]]>=aaa)*1
      res[i,] <- roc.table(est.b01, g[[i]])
    }
    return(res)
  }
  
  roc1.obsmean <- roc1.estmean <- matrix(NA, ncol = 2, nrow = mn)
  for(i in 1:mn){
    obsmean <- roc.threshold(get(paste("theta",i,".b01", sep = "")), g10, threshold = obs.g10)
    roc1.obsmean[i,] <- obsmean
    estmean <- roc.threshold(get(paste("theta",i,".b01", sep = "")), g10, threshold = get(paste("theta",i,".b01", sep = "")))
    roc1.estmean[i,] <- estmean
    print(i)
  }
  roc1.obsmean
  roc1.estmean
  
  plot(theta1.perf, col = mycol[2], lwd = 2, main = paste("M.Fle-Network,", "ROC curve", sep = " "), cex.main = 2, cex.lab = 2, cex.axis = 2)
  plot(theta2.perf, col = mycol[3], lwd = 2, add = TRUE)
  plot(theta3.perf, col = mycol[4], lwd = 2, add = TRUE)
  plot(theta4.perf, col = mycol[5], lwd = 2, add = TRUE)
  points(uni.xy[1], uni.xy[2], col = 2, pch = 16, cex = 2)
  points(inter.xy[1], inter.xy[2], col = 2, pch = 15, cex = 2)
  for (i in 1:mn) points(roc1.obsmean[i,1], roc1.obsmean[i,2], pch = 17, col = mycol[(i+1)], cex = 2)
  for (i in 1:mn) points(roc1.estmean[i,1], roc1.estmean[i,2], pch = 8, col = mycol[(i+1)], cex = 2)
  legend("bottomright", legend = c("Union", "Intersection", expression(t == bar(Y)), expression(t == bar(hat(theta))), mylegend[2:5])
         , pch = c(16,15,17,8,20,20,20,20), col = c(2,2,1,1,mycol[2:5]), cex = 1.5)
  
  
  plot(theta1.perf, col = mycol[2], lwd = 2, main = paste("M.Fle-Network,", "ROC curve", sep = " ")
       , cex.main = 2, cex.lab = 2, cex.axis = 2, xlim = c(0, 0.4), ylim = c(0.6,1))
  plot(theta2.perf, col = mycol[3], lwd = 2, add = TRUE)
  plot(theta3.perf, col = mycol[4], lwd = 2, add = TRUE)
  plot(theta4.perf, col = mycol[5], lwd = 2, add = TRUE)
  points(uni.xy[1], uni.xy[2], col = 2, pch = 16, cex = 2)
  points(inter.xy[1], inter.xy[2], col = 2, pch = 15, cex = 2)
  for (i in 1:mn) points(roc1.obsmean[i,1], roc1.obsmean[i,2], pch = 17, col = mycol[(i+1)], cex = 2)
  for (i in 1:mn) points(roc1.estmean[i,1], roc1.estmean[i,2], pch = 8, col = mycol[(i+1)], cex = 2)
  legend("bottomright", legend = c("Union", "Intersection", expression(t == bar(Y)), expression(t == bar(hat(theta))), mylegend[2:5])
         , pch = c(16,15,17,8,20,20,20,20), col = c(2,2,1,1,mycol[2:5]), cex = 1.5)
  
  
  netnp <- array(dim=c(4500, m*4))
  for(i in 1:6){
    for(k in 1:4){
      #print(k+(i-1)*3)
      netnp[,(k+(i-1)*4)] <- get(paste("theta", k, ".pred.y", sep=""))[[1]][,i]
    }
  }
  mn <- 4
  true_np <- netmean(obs.g10[[1]])
  boxplot(netnp, col = mycol[2:5], main = "Posterior predictive check", cex.main = 2, xaxt = "n")
  for (i in 1:V)  lines(c(4*(i-1)+0.5, 4*i+0.5), rep(true_np[i],2), col = 2, lwd = 2)
  abline(v = (1:m)*mn+0.5, lty = 2)
  
  axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  axis(side=1, at = 1:(mn*m), mgp=c(0,1.5,0), labels = rep(paste0("m",1:mn),6), tick = 0, cex.axis = 0.6)
  legend("topleft", legend = "Observed mean density", lwd = 3, cex = 1.5, col = 2)
  
  
  dev.off()
  
}


ppc.pvalue <- function(obs.data, m, mn){
  p.res <- matrix(NA, nrow = mn, ncol = m)
  obs.mean <- netmean(obs.data)
  
  for(k in 1:mn){
    for(i in 1:m){
      aa <- mean((get(paste("theta", k, ".pred.y", sep = ""))[[1]][,i] >= obs.mean[i])*1)
      p.res[k,i] <- ifelse(aa<0.5, aa, 1-aa)
    }
  }
  return(list(p.res = p.res, obs.mean = obs.mean))
}

ppc.pvalue(obs.g10[[1]], m = 6, mn = 4)


sym.rate <- function(school.data,m){
  school.data[which(is.na(school.data))] <- 0
  school.data <- diag.remove(school.data)
  res <- matrix(NA, nrow = m, ncol = 3)
  asym.rate <- c()
  for(i in 1:m){
    res[i,] <- table(school.data[i,,]+t(school.data[i,,]))
    asym.rate[i] <- res[i,2]/(res[i,2]+res[i,3])
  }
  return(list(res = res, asym.rate = asym.rate))
}

sym.rate(obs.g10[[1]],m = 6)

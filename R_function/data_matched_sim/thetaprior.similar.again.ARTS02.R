## Load function
source("thetaprior_function_forpaper.R")
library(sna)

## number of iteration 
myrep <- 50

## number of relation
m <-  6


n.group <- 3

## choose school name from one of those ("ARTS02", "ARTS03", "PREP02", "TECH03")
school <- "ARTS02"
print(school)

n_draws = 1500
n_tinning = 100
n_burntime = 15000

if(1==2){
  set.seed(1)
  data = c(rbeta(80, 5,5), rbeta(20, 1,9))
  plot(density(data), xlim = c(0,1))
  data = c(rbeta(90, 1,50), rbeta(10, 1, 10))
  plot(density(data), xlim = c(0,1))
}



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
trueep10 <- trueem10 <- matrix(NA, nrow = myrep, ncol = V)


for(k in 1:myrep){
  set.seed(1200+k)
  tmp_ind <- sample(1:V, V*group.p)
  tmp1 <- sample(1:length(tmp_ind), round(length(tmp_ind)/2))
  d_ind <- tmp_ind[tmp1]
  d_ind2 <- tmp_ind[-tmp1]
  
  tmp <- gen.emep.dep.net(V, m, deprho = -0.5, similarity = 0.7, netprob = true_np, d_ind = d_ind, d_ind2 = d_ind2, em = myemp, ep = myepp, myseed = (1200+k))
  
  g10[[k]] <- tmp$g
  obs.g10[[k]] <- tmp$obs.g
  obs.g10[[k]][,mis,] <- NA
  trueem10[k,] <- tmp$trueem
  trueep10[k,] <- tmp$trueep
  
  #### missing 
  
}


a1.em.list <- c(3,1,1)
b1.em.list <- c(5,9,1)
b2.em.list <- c(5,5,30)

a1.theta.list <- c(3,1,1)
b1.theta.list <- c(5,9,1)
b2.theta.list <- c(5,5,30)

theta1.em.density <- theta2.em.density <- theta3.em.density <- theta4.em.density <- theta5.em.density <- one_theta.em.density <- list()
theta1.ep.density <- theta2.ep.density <- theta3.ep.density <- theta4.ep.density <- theta5.ep.density <- one_theta.ep.density <- list()
theta1.em.CI <- theta2.em.CI <- theta3.em.CI <- theta4.em.CI <- theta5.em.CI <- one_theta.em.CI <- list()
theta1.ep.CI <- theta2.ep.CI <- theta3.ep.CI <- theta4.ep.CI <- theta5.ep.CI <- one_theta.ep.CI <- list()
theta1.b01 <- theta2.b01 <- theta3.b01 <- theta4.b01 <- theta5.b01 <- one_theta.b01 <- list()
theta1.pred.y <- theta2.pred.y <- theta3.pred.y <- theta4.pred.y <- theta5.pred.y <- one_theta.pred.y <- list()
theta1.np <- theta2.np <- theta3.np <- theta4.np <- theta5.np <- one_theta.np <- list()


for(k in 1:myrep){
  
  assign("one_theta", bnsm.one.theta(dat = obs.g10[[k]]
                                     , emp = c(1,11), epp = c(1,11)
                                     , thetaalpha = 0.5, thetabeta = 0.5
                                     , diag = FALSE, mode = "graph", model.checking = TRUE
                                     , reps = 3, draws = n_draws, tinning = n_draws, burntime = n_burntime
                                     , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE))
  
  one_theta.em.density[[k]] <- density(one_theta$em)
  one_theta.ep.density[[k]] <- density(one_theta$ep)
  one_theta.em.CI[[k]] <- apply(one_theta$em, 2, quantile, prob = c(0.025,0.975, 0.5))
  one_theta.ep.CI[[k]] <- apply(one_theta$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
  one_theta.b01[[k]] <- apply(one_theta$net, c(3,4,5), mean)
  one_theta.pred.y[[k]] <- pred.y.mean(one_theta$pred.y, iter = n_draws, reps = 3, m = 6)
  one_theta.np[[k]] <- one_theta$np
  rm(one_theta)
  
  for(i in 1:length(b1.em.list)){
    a1.em = a1.em.list[i]; b1.em = b1.em.list[i]; a2.em = 1; b2.em = b2.em.list[i]
    a1.ep = 1; b1.ep = 1; a2.ep = 1; b2.ep = 30 ### uniform prior
    a1.theta = a1.theta.list[i]; b1.theta = b1.theta.list[i]; a2.theta = 1; b2.theta = b2.theta.list[i] 
    
    alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
    alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
    alpha.acc.theta = 0; beta.acc.theta = 0; total.alpha.theta = 0; total.beta.theta = 0
    
    assign(paste("theta",i,sep=""), bnsm.hyper.theta.hyper(dat = obs.g10[[k]], emsd = list(ema = 10, emb = 10, alpha.sd = .5, beta.sd = .5)
                                                           , epsd = list(epa = 10, epb = 10, alpha.sd = .5, beta.sd = .5)
                                                           , thetasd = list(thetaa = 3, thetab = 9, alpha.sd = .5, beta.sd = .5)
                                                           , diag = FALSE, mode = "graph", model.checking = TRUE
                                                           , reps = 3, draws = n_draws, tinning = n_tinning, burntime = n_burntime
                                                           , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = FALSE))
    
    
    #save.image(file = paste(V, "nodes.thetaprior.modelcheck.33", n.group, true_np, paste(myemp, collapse = ""), note, "RData", sep = "."))
    #save.image(file = "thetaprior.similar.TECH03.1915.RData")
    
    if(i==1){
      theta1.em.density[[k]] <- density(theta1$em)
      theta1.ep.density[[k]] <- density(theta1$ep)
      theta1.em.CI[[k]] <- apply(theta1$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta1.ep.CI[[k]] <- apply(theta1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta1.b01[[k]] <- apply(theta1$net, c(3,4,5), mean)
      theta1.pred.y[[k]] <- pred.y.mean(theta1$pred.y, iter = n_draws, reps = 3, m = 6)
      theta1.np[[k]] <- theta1$np
      rm(theta1)
    }else if(i==2){
      theta2.em.density[[k]] <- density(theta2$em)
      theta2.ep.density[[k]] <- density(theta2$ep)
      theta2.em.CI[[k]] <- apply(theta2$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta2.ep.CI[[k]] <- apply(theta2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta2.b01[[k]] <- apply(theta2$net, c(3,4,5), mean)
      theta2.pred.y[[k]] <- pred.y.mean(theta2$pred.y, iter = n_draws, reps = 3, m = 6)
      theta2.np[[k]] <- theta2$np
    }else if(i==3){
      theta3.em.density[[k]] <- density(theta3$em)
      theta3.ep.density[[k]] <- density(theta3$ep)
      theta3.em.CI[[k]] <- apply(theta3$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta3.ep.CI[[k]] <- apply(theta3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta3.b01[[k]] <- apply(theta3$net, c(3,4,5), mean)
      theta3.pred.y[[k]] <- pred.y.mean(theta3$pred.y, iter = n_draws, reps = 3, m = 6)
      theta3.np[[k]] <- theta3$np
      rm(theta3)
    }else if(i==4){
      theta4.em.density[[k]] <- density(theta4$em)
      theta4.ep.density[[k]] <- density(theta4$ep)
      theta4.em.CI[[k]] <- apply(theta4$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta4.ep.CI[[k]] <- apply(theta4$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta4.b01[[k]] <- apply(theta4$net, c(3,4,5), mean)
      theta4.pred.y[[k]] <- pred.y.mean(theta4$pred.y, iter = n_draws, reps = 3, m = 6)
      theta4.np[[k]] <- theta4$np
      rm(theta4)
    }else if(i==5){
      theta5.em.density[[k]] <- density(theta5$em)
      theta5.ep.density[[k]] <- density(theta5$ep)
      theta5.em.CI[[k]] <- apply(theta5$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta5.ep.CI[[k]] <- apply(theta5$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta5.b01[[k]] <- apply(theta5$net, c(3,4,5), mean)
      theta5.pred.y[[k]] <- pred.y.mean(theta5$pred.y, iter = n_draws, reps = 3, m = 6)
      theta5.np[[k]] <- theta5$np
      rm(theta5)
    }
    print(k)
    
    save.image(file = paste("thetaprior.similar.again.rep50.", school, ".RData", sep=""))
    
    
  }
}
#b2.em.list <- c(100,100,100,100,100,100)
#b1.ep.list <- c(9,9,9,3,3,3)



if(2==3){
  sum((obs.g10[[1]][1,,]==0 & t(obs.g10[[1]][1,,])==0), na.rm = TRUE)/(V*(V-1))
  2*sum((obs.g10[[1]][1,,]==1 & t(obs.g10[[1]][1,,])==0), na.rm = TRUE)/(V*(V-1))
  sum((obs.g10[[1]][1,,]==1 & t(obs.g10[[1]][1,,])==1), na.rm = TRUE)/(V*(V-1))
  
  load("/Users/dongahkim/Downloads/thetaprior.similar.again.rep50.1.ARTS02.RData")
  mycol <- c("gray", "pink", "yellow", "green")
  mylegend1 = c("m1", "m2", "m3", "m4")
  
  xs <- seq(0,1,by = 0.001)
  emy <- 0.9*dbeta(xs, 2,5)+0.1*dbeta(xs, 10, 1)
  epy <- 0.9*dbeta(xs, 1,50)+0.1*dbeta(xs, 1, 10)
  
  plot(xs, emy, type = "l", lwd = 3, col = 2, xlab = "e-", ylab = "Density", main = "Posterior for all observations in 50 reps")
  
  for(i in 1:myrep){
    lines(theta1.em.density[[i]], col = mycol[1], lwd = 1)
    lines(theta2.em.density[[i]], col = mycol[2], lwd = 1)
    lines(theta3.em.density[[i]], col = mycol[3], lwd = 1)
    #lines(theta4.em.density[[i]], col = mycol[4], lwd = 2)
  }
  plot(xs, epy, type = "l", lwd = 3, col = 2, xlab = "e+", ylab = "Density", main = "Posterior for all observations in 50 reps")
  for(i in 1:myrep){
    lines(theta1.ep.density[[i]], col = mycol[1], lwd = 2)
    lines(theta2.ep.density[[i]], col = mycol[2], lwd = 2)
    lines(theta3.ep.density[[i]], col = mycol[3], lwd = 2)
    #lines(theta4.ep.density[[i]], col = mycol[4], lwd = 2)
  }
  
  
  
  mean.np <- matrix(NA, ncol = 6, nrow = myrep)
  mean.np1 <- mean.np2 <- mean.np3 <- mean.np4 <- matrix(NA, ncol = 6, nrow = myrep)
  q975.np1 <- q975.np2 <- q975.np3 <- q975.np4 <- matrix(NA, ncol = 6, nrow = myrep)
  q025.np1 <- q025.np2 <- q025.np3 <- q025.np4 <- matrix(NA, ncol = 6, nrow = myrep)
  
  for(i in 1:myrep){
    mean.np1[i,] <- apply(theta1.np[[i]], 2, mean)
    mean.np2[i,] <- apply(theta2.np[[i]], 2, mean)
    mean.np3[i,] <- apply(theta3.np[[i]], 2, mean)
    #mean.np4[i,] <- apply(theta4.np[[i]], 2, mean)
    q975.np1[i,] <- apply(theta1.np[[i]], 2, quantile, prob = 0.975)
    q975.np2[i,] <- apply(theta2.np[[i]], 2, quantile, prob = 0.975)
    q975.np3[i,] <- apply(theta3.np[[i]], 2, quantile, prob = 0.975)
    #q975.np4[i,] <- apply(theta4.np[[i]], 2, quantile, prob = 0.975)
    q025.np1[i,] <- apply(theta1.np[[i]], 2, quantile, prob = 0.025)
    q025.np2[i,] <- apply(theta2.np[[i]], 2, quantile, prob = 0.025)
    q025.np3[i,] <- apply(theta3.np[[i]], 2, quantile, prob = 0.025)
    #q025.np4[i,] <- apply(theta4.np[[i]], 2, quantile, prob = 0.025)
  }
  #boxplot(mean.np1)
  #boxplot(rbind(mean.np1, q975.np1, q025.np1))
  #boxplot(theta2.np[[1]])
  
  npbox <- cbind(rbind(mean.np1, q975.np1, q025.np1), rbind(mean.np2, q975.np2, q025.np2), rbind(mean.np3, q975.np3, q025.np3))
  
  boxplot(npbox[,c(1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12,18)], col = mycol, main = "ARTS02, Posterior for Theta", cex.main = 2, xaxt = "n")
  for (i in 1:V)  lines(c(3*(i-1)+0.5, 3*i+0.5), rep(true_np[i],2), col = 2, lwd = 2)
  abline(v = c(3.5, 6.5, 9.5, 12.5, 15.5), lty = 2)
  axis(side=1, at=c(2, 5, 8, 11, 14, 17), mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  axis(side=1, at = 1:18, mgp=c(0,1.5,0), labels = rep(paste0("m",1:3),6), tick = 0, las = 2)
  
  union.net <- union.intersection.net(obs.g10[[1]])$union.net
  intersection.net <- union.intersection.net(obs.g10[[1]])$intersection.net
  
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
    return(res)
  }
  
  uni.inter.xy <- function(obs, g){
    union.net <- union.intersection.net(obs)$union.net
    intersection.net <- union.intersection.net(obs)$intersection.net
    uni.xy <- roc.table(union.net, g)
    inter.xy <- roc.table(intersection.net, g)
    res <- c(uni.xy, inter.xy)
    names(res) <- c("uni.x", "uni.y", "inter.x", "inter.y")
    return(res)
  }
  
  
  
  #uni.inter.xy(obs.g10[[1]], g10[[1]])
  
  uni.inter.xy.table <- matrix(NA, ncol = 4, nrow = myrep)
  for(i in 1:myrep){
    uni.inter.xy.table[i,] <- uni.inter.xy(obs.g10[[i]], g10[[i]])
  }
  
  colnames(uni.inter.xy.table) <- c("uni.x", "uni.y", "inter.x", "inter.y")
  uni.inter.xy.table
  
  roc.theta1.b01 <- roc.theta2.b01 <- roc.theta3.b01 <- roc.theta4.b01 <- roc.g10 <- list()
  roc.theta1.b <- roc.theta2.b <- roc.theta3.b <- roc.theta4.b <- roc.g10.b <- list()
  
  for(i in 1:myrep){
    roc.theta1.b[[i]] <- c(theta1.b01[[i]])
    roc.theta2.b[[i]] <- c(theta2.b01[[i]])
    roc.theta3.b[[i]] <- c(theta3.b01[[i]])
    #roc.theta4.b[[i]] <- c(theta4.b01[[i]])
    roc.g10.b[[i]] <- c(g10[[i]])
    
    diagna <- which(is.na(roc.theta2.b[[i]]))
    roc.theta1.b01[[i]] <- roc.theta1.b[[i]][-diagna] 
    roc.theta2.b01[[i]] <- roc.theta2.b[[i]][-diagna] 
    roc.theta3.b01[[i]] <- roc.theta3.b[[i]][-diagna] 
    #roc.theta4.b01[[i]] <- roc.theta4.b[[i]][-diagna] 
    roc.g10[[i]] <- roc.g10.b[[i]][-diagna]
    rm(diagna)
  }
  
  theta1.pred <- theta2.pred <- theta3.pred <- theta4.pred <- list()
  theta1.perf <- theta2.perf <- theta3.perf <- theta4.perf <- list()
  for(i in 1:myrep){
    theta1.pred[[i]] <- prediction(roc.theta1.b01[[i]], roc.g10[[i]])
    theta2.pred[[i]] <- prediction(roc.theta2.b01[[i]], roc.g10[[i]])
    theta3.pred[[i]] <- prediction(roc.theta3.b01[[i]], roc.g10[[i]])
    #theta4.pred[[i]] <- prediction(roc.theta4.b01[[i]], roc.g10[[i]])
    
    theta1.perf[[i]] <- performance(theta1.pred[[i]],"tpr","fpr")
    theta2.perf[[i]] <- performance(theta2.pred[[i]],"tpr","fpr")
    theta3.perf[[i]] <- performance(theta3.pred[[i]],"tpr","fpr")
    #theta4.perf[[i]] <- performance(theta4.pred[[i]],"tpr","fpr")
    
    print(i)
  }
  
  myinteger <- seq(0,1,by = 0.01)
  mylocation <- c()
  ROC.rep.x <- ROC.rep.y <- matrix(NA, ncol = length(myinteger), nrow = myrep)
  for(i in 1:myrep){
    for(j in 1:length(myinteger)){
      mylocation[j] <- max(which((theta1.perf[[1]]@alpha.values)[[1]]>=myinteger[j]))
    }
    ROC.rep.x[i,] <- as.data.frame(theta1.perf[[i]]@x.values)[,1][mylocation]
    ROC.rep.y[i,] <- as.data.frame(theta1.perf[[i]]@y.values)[,1][mylocation]
    
  }
  
  res <- res1 <- c()
  for(i in 1:50){
    res[i] <- ifelse(ROC.rep.y[i,min(which(ROC.rep.x[i,] < uni.inter.xy.table[i,1]))]<uni.inter.xy.table[i,2],1,0)
    res1[i] <- ifelse(ROC.rep.y[i,max(which(ROC.rep.x[i,] > uni.inter.xy.table[i,1]))]<uni.inter.xy.table[i,2],1,0)
    
  }
  mean(res)
  mean(res1)
  
  
  plot(theta1.perf[[1]], col = "gray", main = paste("ARTS02,", "ROC curve", sep = " "), cex.main = 2, cex.lab = 2, cex.axis = 2)
  plot(theta2.perf[[1]], col = "pink", add = TRUE)
  plot(theta3.perf[[1]], col = "yellow", add = TRUE)
  
  for(i in 2:myrep){
    plot(theta1.perf[[i]], col = "gray", add = TRUE)
    plot(theta2.perf[[i]], col = "pink", add = TRUE)
    plot(theta3.perf[[i]], col = "yellow", add = TRUE)
  }
  
  points(uni.inter.xy.table[,1], uni.inter.xy.table[,2])
  points(uni.inter.xy.table[,3], uni.inter.xy.table[,4])
  
  points(mean(uni.inter.xy.table[,1]), mean(uni.inter.xy.table[,2]), col = 2, pch = 16)
  points(mean(uni.inter.xy.table[,3]), mean(uni.inter.xy.table[,4]), col = 2, pch = 15)
  
  
  lines(apply(ROC.rep.x,2,mean, na.rm = TRUE), apply(ROC.rep.y,2,mean, na.rm = TRUE), type = "l", col = 2, lwd = 2)
  legend("bottomright", legend = c(paste("Union : ", mean(res), sep = ""), "Intersection", "mean(ROC)"), pch = c(16,15,NA), lty = c(NA,NA,1), lwd = 2, col = c(2,2,2), cex = 2)
  
  
  dev.off()
  ###############################################
  if(2==3){
    rowMeans(theta1.b01[[1]], na.rm = TRUE)
    
    #theta1.perf[[1]]@x.values
    plot(theta1.perf[[i]], col = "gray", lwd = 2, main = paste("ARTS02,", "ROC curve", sep = " "), cex.main = 2, cex.lab = 2, cex.axis = 2)
    plot(as.data.frame(theta1.perf[[i]]@x.values)[,1], as.data.frame(theta1.perf[[i]]@y.values)[,1])
    
    for(i in 1:6){
      print(length(as.data.frame(theta1.perf[[i]]@y.values)[,1]))
    }
    
    fpr1 <- tpr1 <- matrix(NA, nrow = myrep, ncol = length(as.data.frame(theta1.perf[[i]]@x.values)[,1]))
    fpr2 <- tpr2 <- matrix(NA, nrow = myrep, ncol = length(as.data.frame(theta2.perf[[i]]@x.values)[,1]))
    fpr3 <- tpr3 <- matrix(NA, nrow = myrep, ncol = length(as.data.frame(theta3.perf[[i]]@x.values)[,1]))
    #fpr4 <- tpr4 <- matrix(NA, nrow = myrep, ncol = length(as.data.frame(theta4.perf[[i]]@x.values)[,1]))
    
    for(i in 1:myrep){
      fpr1[i,] <- as.data.frame(theta1.perf[[i]]@x.values)[,1]
      tpr1[i,] <- as.data.frame(theta1.perf[[i]]@y.values)[,1]
      fpr2[i,] <- as.data.frame(theta2.perf[[i]]@x.values)[,1]
      tpr2[i,] <- as.data.frame(theta2.perf[[i]]@y.values)[,1]
      fpr3[i,] <- as.data.frame(theta3.perf[[i]]@x.values)[,1]
      tpr3[i,] <- as.data.frame(theta3.perf[[i]]@y.values)[,1]
      #fpr4[i,] <- as.data.frame(theta4.perf[[i]]@x.values)[,1]
      #tpr4[i,] <- as.data.frame(theta4.perf[[i]]@y.values)[,1]
      print(i)
    }
    mean.fpr1 <- colMeans(fpr1, na.rm = TRUE)
    mean.tpr1 <- colMeans(tpr1, na.rm = TRUE)
    mean.fpr2 <- colMeans(fpr2, na.rm = TRUE)
    mean.tpr2 <- colMeans(tpr2, na.rm = TRUE)
    mean.fpr3 <- colMeans(fpr3, na.rm = TRUE)
    mean.tpr3 <- colMeans(tpr3, na.rm = TRUE)
    #mean.fpr4 <- colMeans(fpr4, na.rm = TRUE)
    #mean.tpr4 <- colMeans(tpr4, na.rm = TRUE)
    
    plot(mean.fpr2, mean.tpr2, type = "l", lwd = 3, col = "gray")
    points(mean(uni.inter.xy.table[,1], na.rm = TRUE), mean(uni.inter.xy.table[,2], na.rm = TRUE), col = 2, pch = 16, cex = 2)
    points(mean(uni.inter.xy.table[,3], na.rm = TRUE), mean(uni.inter.xy.table[,4], na.rm = TRUE), col = 2, pch = 15, cex = 2)
    
    plot(mean.fpr, mean.tpr, type = "l", lwd = 3, col = "gray", xlim = c(0,0.2), ylim = c(0.4, 1))
    points(mean(uni.inter.xy.table[,1], na.rm = TRUE), mean(uni.inter.xy.table[,2], na.rm = TRUE), col = 2, pch = 16, cex = 2)
    points(mean(uni.inter.xy.table[,3], na.rm = TRUE), mean(uni.inter.xy.table[,4], na.rm = TRUE), col = 2, pch = 15, cex = 2)
    
    
    plot(theta1.perf[[1]], col = "gray", lwd = 2, main = paste("ARTS02,", "ROC curve", sep = " "), cex.main = 2, cex.lab = 2, cex.axis = 2)
    for(i in 2:13){
      plot(theta1.perf[[i]], col = "gray", lwd = 1, add = TRUE)
      
    }
    
    
    
    
    
    
    
    
    
    
  }
}


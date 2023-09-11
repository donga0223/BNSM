rm(list=ls())

## choose school name that you want to summarize
school <- "ARTS03"
xs <- seq(0,1,by = 0.001)
setwd("/Users/dongahkim/OneDrive/BNSM/paper simulation/data_matched_sim")
if(school == "ARTS02"){
  emy <- 0.9*dbeta(xs, 2,5)+0.1*dbeta(xs, 10, 1)
  epy <- 0.9*dbeta(xs, 1,50)+0.1*dbeta(xs, 1, 10)
}else if(school == "ARTS03"){
  emy <- 0.9*dbeta(xs, 2,5)+0.1*dbeta(xs, 5, 2)
  epy <- 0.9*dbeta(xs, 1,50)+0.1*dbeta(xs, 1, 10)
}else if(school == "PREP02"){
  emy <- 0.5*dbeta(xs, 1,1)+0.5*dbeta(xs, 5, 2)
  epy <- 0.75*dbeta(xs, 1,100)+0.25*dbeta(xs, 1, 9)
}else if(school == "TECH03"){
  emy <- 0.8*dbeta(xs, 5,5)+0.2*dbeta(xs, 3,9)
  epy <- 0.9*dbeta(xs, 1,100)+0.1*dbeta(xs, 1, 1)
}


#load("/Users/dongahkim/Downloads/thetaprior.similar.again.rep50.1.ARTS03.RData")
#load("/Users/dongahkim/Downloads/thetaprior.similar.again.rep50.5.ARTS03.RData")
load("/Users/dongahkim/Downloads/thetaprior.similar.again.rep50.3.ARTS03.RData")

mylegend1 = c("m1", "m4", "m5")
simname <- paste("S-", school, sep = "")
col.list = c("gray", "pink", "yellow", "green", "skyblue")
mn <- length(mylegend1)
mycol =  col.list[1:mn]
t.n <- mn*6
myrep <- length(theta3.b01)

pdf("S-ARTS03.again.50reps.final.pdf")

par(mar = c(4.5, 5, 3.5, 1.5))

plot(xs, emy, type = "l", lwd = 3, col = 2, xlab = "e-", ylab = "Density", main = paste(simname, ", Posterior for all over 50 reps", sep = ""), cex.main = 2, cex.lab = 2, cex.axis = 2)

for(i in 1:myrep){
  lines(theta1.em.density[[i]], col = mycol[1], lwd = 1)
  lines(theta2.em.density[[i]], col = mycol[2], lwd = 1)
  lines(theta3.em.density[[i]], col = mycol[3], lwd = 1)
  #lines(theta4.em.density[[i]], col = mycol[4], lwd = 1)
  #lines(theta5.em.density[[i]], col = mycol[5], lwd = 1)
}
lines(xs, emy, type = "l", lwd = 3, col = 2)
legend("topright", legend = mylegend1, col = mycol, lwd = 3, cex = 2)

plot(xs, epy, type = "l", lwd = 3, col = 2, xlab = "e+", ylab = "Density", main = paste(simname, ", Posterior for all over 50 reps", sep = ""), cex.main = 2, cex.lab = 2, cex.axis = 2)
for(i in 1:myrep){
  lines(theta1.ep.density[[i]], col = mycol[1], lwd = 2)
  lines(theta2.ep.density[[i]], col = mycol[2], lwd = 2)
  lines(theta3.ep.density[[i]], col = mycol[3], lwd = 2)
  #lines(theta4.ep.density[[i]], col = mycol[4], lwd = 2)
  #lines(theta5.ep.density[[i]], col = mycol[4], lwd = 2)
}
lines(xs, epy, type = "l", lwd = 3, col = 2)
legend("topright", legend = mylegend1, col = mycol, lwd = 3, cex = 2)



mean.np <- matrix(NA, ncol = 6, nrow = myrep)
mean.np1 <- mean.np2 <- mean.np3 <- mean.np4 <- mean.np5 <- matrix(NA, ncol = 6, nrow = myrep)
q975.np1 <- q975.np2 <- q975.np3 <- q975.np4 <- q975.np5 <- matrix(NA, ncol = 6, nrow = myrep)
q025.np1 <- q025.np2 <- q025.np3 <- q025.np4 <- q025.np5 <- matrix(NA, ncol = 6, nrow = myrep)

for(i in 1:myrep){
  mean.np1[i,] <- apply(theta1.np[[i]], 2, mean)
  mean.np2[i,] <- apply(theta2.np[[i]], 2, mean)
  mean.np3[i,] <- apply(theta3.np[[i]], 2, mean)
  #mean.np4[i,] <- apply(theta4.np[[i]], 2, mean)
  #mean.np5[i,] <- apply(theta5.np[[i]], 2, mean)
  q975.np1[i,] <- apply(theta1.np[[i]], 2, quantile, prob = 0.975)
  q975.np2[i,] <- apply(theta2.np[[i]], 2, quantile, prob = 0.975)
  q975.np3[i,] <- apply(theta3.np[[i]], 2, quantile, prob = 0.975)
  #q975.np4[i,] <- apply(theta4.np[[i]], 2, quantile, prob = 0.975)
  #q975.np5[i,] <- apply(theta5.np[[i]], 2, quantile, prob = 0.975)
  q025.np1[i,] <- apply(theta1.np[[i]], 2, quantile, prob = 0.025)
  q025.np2[i,] <- apply(theta2.np[[i]], 2, quantile, prob = 0.025)
  q025.np3[i,] <- apply(theta3.np[[i]], 2, quantile, prob = 0.025)
  #q025.np4[i,] <- apply(theta4.np[[i]], 2, quantile, prob = 0.025)
  #q025.np5[i,] <- apply(theta5.np[[i]], 2, quantile, prob = 0.025)
}
#boxplot(mean.np1)
#boxplot(rbind(mean.np1, q975.np1, q025.np1))
#boxplot(theta2.np[[1]])

npbox <- cbind(rbind(mean.np1, q975.np1, q025.np1), rbind(mean.np2, q975.np2, q025.np2), rbind(mean.np3, q975.np3, q025.np3))
               #, rbind(mean.np4, q975.np4, q025.np4), rbind(mean.np5, q975.np5, q025.np5))

boxplot(npbox[,c(seq(1,t.n,6),seq(2,t.n,6),seq(3,t.n,6),seq(4,t.n,6),seq(5,t.n,6),seq(6,t.n,6))], col = mycol, main = paste(simname, ", Posterior for Theta over 50 reps", sep = ""), cex.main = 2, xaxt = "n")
abline(v = (1:m)*mn+0.5, lty = 2)
axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
axis(side=1, at = 1:(mn*m), mgp=c(0,1.5,0), labels = rep(paste0("m",c(1,4,5)),6), tick = 0, cex.axis = 0.7)
for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(true_np[i],2), col = 2, lwd = 2)
legend("topleft", legend = "True mean density", lwd = 3, cex = 1.5, col = 2)


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

roc.theta1.b01 <- roc.theta2.b01 <- roc.theta3.b01 <- roc.theta4.b01 <- roc.theta5.b01 <- roc.g10 <- list()
roc.theta1.b <- roc.theta2.b <- roc.theta3.b <- roc.theta4.b <- roc.theta5.b <- roc.g10.b <- list()

for(i in 1:myrep){
  roc.theta1.b[[i]] <- c(theta1.b01[[i]])
  roc.theta2.b[[i]] <- c(theta2.b01[[i]])
  roc.theta3.b[[i]] <- c(theta3.b01[[i]])
  #roc.theta4.b[[i]] <- c(theta4.b01[[i]])
  #roc.theta5.b[[i]] <- c(theta5.b01[[i]])
  roc.g10.b[[i]] <- c(g10[[i]])
  
  diagna <- which(is.na(roc.theta2.b[[i]]))
  roc.theta1.b01[[i]] <- roc.theta1.b[[i]][-diagna] 
  roc.theta2.b01[[i]] <- roc.theta2.b[[i]][-diagna] 
  roc.theta3.b01[[i]] <- roc.theta3.b[[i]][-diagna] 
  #roc.theta4.b01[[i]] <- roc.theta4.b[[i]][-diagna] 
  #roc.theta5.b01[[i]] <- roc.theta5.b[[i]][-diagna] 
  roc.g10[[i]] <- roc.g10.b[[i]][-diagna]
  rm(diagna)
}

theta1.pred <- theta2.pred <- theta3.pred <- theta4.pred <- theta5.pred <- list()
theta1.perf <- theta2.perf <- theta3.perf <- theta4.perf <- theta5.perf <- list()
for(i in 1:myrep){
  theta1.pred[[i]] <- prediction(roc.theta1.b01[[i]], roc.g10[[i]])
  theta2.pred[[i]] <- prediction(roc.theta2.b01[[i]], roc.g10[[i]])
  theta3.pred[[i]] <- prediction(roc.theta3.b01[[i]], roc.g10[[i]])
  #theta4.pred[[i]] <- prediction(roc.theta4.b01[[i]], roc.g10[[i]])
  #theta5.pred[[i]] <- prediction(roc.theta5.b01[[i]], roc.g10[[i]])
  
  theta1.perf[[i]] <- performance(theta1.pred[[i]],"tpr","fpr")
  theta2.perf[[i]] <- performance(theta2.pred[[i]],"tpr","fpr")
  theta3.perf[[i]] <- performance(theta3.pred[[i]],"tpr","fpr")
  #theta4.perf[[i]] <- performance(theta4.pred[[i]],"tpr","fpr")
  #theta5.perf[[i]] <- performance(theta5.pred[[i]],"tpr","fpr")
  
  print(i)
}

myinteger <- seq(0,1,by = 0.01)
mylocation <- c()
ROC.rep.x <- ROC.rep.y <- matrix(NA, ncol = length(myinteger), nrow = myrep)
for(i in 1:myrep){
  for(j in 1:length(myinteger)){
    mylocation[j] <- max(which((theta1.perf[[i]]@alpha.values)[[1]]>=myinteger[j]))
  }
  ROC.rep.x[i,] <- as.data.frame(theta1.perf[[i]]@x.values)[,1][mylocation]
  ROC.rep.y[i,] <- as.data.frame(theta1.perf[[i]]@y.values)[,1][mylocation]
  
}

res <- res1 <- c()
for(i in 1:myrep){
  res[i] <- ifelse(ROC.rep.y[i,min(which(ROC.rep.x[i,] < uni.inter.xy.table[i,1]))]<uni.inter.xy.table[i,2],1,0)
  res1[i] <- ifelse(ROC.rep.y[i,max(which(ROC.rep.x[i,] > uni.inter.xy.table[i,1]))]<uni.inter.xy.table[i,2],1,0)
  
}
mean(res, na.rm = TRUE)
mean(res1, na.rm = TRUE)

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

roc1.obsmean <- roc.threshold(theta1.b01, g10, threshold = obs.g10)
roc1.estmean <- roc.threshold(theta1.b01, g10, threshold = theta1.b01)

plot(theta1.perf[[1]], col = mycol[1], main = paste(simname, ", ROC curve", sep = ""), cex.main = 2, cex.lab = 2, cex.axis = 2)
plot(theta2.perf[[1]], col = mycol[2], add = TRUE)
plot(theta3.perf[[1]], col = mycol[3], add = TRUE)
#plot(theta4.perf[[1]], col = mycol[4], add = TRUE)
#plot(theta5.perf[[1]], col = mycol[5], add = TRUE)

for(i in 2:myrep){
  plot(theta1.perf[[i]], col = mycol[1], add = TRUE)
  plot(theta2.perf[[i]], col = mycol[2], add = TRUE)
  plot(theta3.perf[[i]], col = mycol[3], add = TRUE)
  #plot(theta4.perf[[i]], col = mycol[4], add = TRUE)
  #plot(theta5.perf[[i]], col = mycol[5], add = TRUE)
}


points(uni.inter.xy.table[,1], uni.inter.xy.table[,2], col = "gray")
points(uni.inter.xy.table[,3], uni.inter.xy.table[,4], pch = 0, col = "gray")


lines(apply(ROC.rep.x,2,mean, na.rm = TRUE), apply(ROC.rep.y,2,mean, na.rm = TRUE), type = "l", col = 2, lwd = 2)
points(roc1.obsmean[,1], roc1.obsmean[,2], pch = 2, col = "gray")
points(roc1.estmean[,1], roc1.estmean[,2], pch = 4, col = "gray")

points(mean(uni.inter.xy.table[,1]), mean(uni.inter.xy.table[,2]), col = 2, pch = 16, cex = 2)
points(mean(uni.inter.xy.table[,3]), mean(uni.inter.xy.table[,4]), col = 2, pch = 15, cex = 2)
points(mean(roc1.obsmean[,1]), mean(roc1.obsmean[,2]), pch = 17, col = 2, cex = 2)
points(mean(roc1.estmean[,1]), mean(roc1.estmean[,2]), pch = 8, col = 2, cex = 2)

legend("bottomright", pch = c(16,15,17,8,NA), lty = c(NA,NA,NA,NA,1), lwd = 2, col = 2, cex = 2
       , legend = c(paste("Union : ", round(mean(res, na.rm = TRUE),3), sep = ""), "Intersection", expression(t == bar(Y)), expression(t == bar(hat(theta))), "mean(ROC)"))





mean.pred.y <- matrix(NA, ncol = 6, nrow = myrep)
mean.pred.y1 <- mean.pred.y2 <- mean.pred.y3 <- mean.pred.y4 <- mean.pred.y5 <- matrix(NA, ncol = 6, nrow = myrep)
q975.pred.y1 <- q975.pred.y2 <- q975.pred.y3 <- q975.pred.y4 <- q975.pred.y5 <- matrix(NA, ncol = 6, nrow = myrep)
q025.pred.y1 <- q025.pred.y2 <- q025.pred.y3 <- q025.pred.y4 <- q025.pred.y5 <- matrix(NA, ncol = 6, nrow = myrep)

for(i in 1:myrep){
  mean.pred.y1[i,] <- apply(theta1.pred.y[[i]], 2, mean)
  mean.pred.y2[i,] <- apply(theta2.pred.y[[i]], 2, mean)
  mean.pred.y3[i,] <- apply(theta3.pred.y[[i]], 2, mean)
  #mean.pred.y4[i,] <- apply(theta4.pred.y[[i]]$res, 2, mean)
  #mean.pred.y5[i,] <- apply(theta5.pred.y[[i]]$res, 2, mean)
  
  q975.pred.y1[i,] <- apply(theta1.pred.y[[i]], 2, quantile, prob = 0.975)
  q975.pred.y2[i,] <- apply(theta2.pred.y[[i]], 2, quantile, prob = 0.975)
  q975.pred.y3[i,] <- apply(theta3.pred.y[[i]], 2, quantile, prob = 0.975)
  #q975.pred.y4[i,] <- apply(theta4.pred.y[[i]]$res, 2, quantile, prob = 0.975)
  #q975.pred.y5[i,] <- apply(theta5.pred.y[[i]]$res, 2, quantile, prob = 0.975)
  
  q025.pred.y1[i,] <- apply(theta1.pred.y[[i]], 2, quantile, prob = 0.025)
  q025.pred.y2[i,] <- apply(theta2.pred.y[[i]], 2, quantile, prob = 0.025)
  q025.pred.y3[i,] <- apply(theta3.pred.y[[i]], 2, quantile, prob = 0.025)
  #q025.pred.y4[i,] <- apply(theta4.pred.y[[i]]$res, 2, quantile, prob = 0.025)
  #q025.pred.y5[i,] <- apply(theta5.pred.y[[i]]$res, 2, quantile, prob = 0.025)
}

obs.netmean <- matrix(NA, ncol = 6, nrow = myrep)
for(i in 1:myrep){
  obs.netmean[i,] <- netmean(obs.g10[[i]])
}
#boxplot(mean.pred.y1)
#boxplot(rbind(mean.pred.y1, q975.pred.y1, q025.pred.y1))
#boxplot(theta2.pred.y[[1]])

npbox <- cbind(rbind(mean.pred.y1, q975.pred.y1, q025.pred.y1)
               , rbind(mean.pred.y2, q975.pred.y2, q025.pred.y2)
               , rbind(mean.pred.y3, q975.pred.y3, q025.pred.y3))
               #, rbind(mean.pred.y4, q975.pred.y4, q025.pred.y4)
               #, rbind(mean.pred.y5, q975.pred.y5, q025.pred.y5))
#npbox <- cbind(mean.pred.y1, mean.pred.y2, mean.pred.y3
boxplot(npbox[,c(seq(1,t.n,6),seq(2,t.n,6),seq(3,t.n,6),seq(4,t.n,6),seq(5,t.n,6),seq(6,t.n,6))], col = mycol
        , main = paste(simname, ", Posterior Predictive check", sep = ""), cex.main = 2, xaxt = "n")
abline(v = (1:m)*mn+0.5, lty = 2)
axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
axis(side=1, at = 1:(mn*m), mgp=c(0,1.5,0), labels = rep(paste0("m",c(1,4,5)),6), tick = 0, cex.axis = 0.7)
#for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(c(0.05, 0.09, 0.01, 0.025, 0.04, 0.13)[i],2), col = 2, lwd = 2)

for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(apply(obs.netmean,2,mean)[i],2), col = 2, lwd = 2)
legend("topleft", legend = "Observed mean density", lwd = 3, cex = 1.5, col = 2)

npbox <- cbind(mean.pred.y1, mean.pred.y2, mean.pred.y3, mean.pred.y4, mean.pred.y5)
#, rbind(mean.pred.y4, q975.pred.y4, q025.pred.y4)
#, rbind(mean.pred.y5, q975.pred.y5, q025.pred.y5))
#npbox <- cbind(mean.pred.y1, mean.pred.y2, mean.pred.y3
boxplot(npbox[,c(seq(1,t.n,6),seq(2,t.n,6),seq(3,t.n,6),seq(4,t.n,6),seq(5,t.n,6),seq(6,t.n,6))], col = mycol
        , main = paste(simname, ", Posterior Predictive check", sep = ""), cex.main = 2, xaxt = "n")
abline(v = (1:m)*mn+0.5, lty = 2)
axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
axis(side=1, at = 1:(mn*m), mgp=c(0,1.5,0), labels = rep(paste0("m",c(1,4,5)),6), tick = 0, cex.axis = 0.7)
#for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(c(0.05, 0.09, 0.01, 0.025, 0.04, 0.13)[i],2), col = 2, lwd = 2)

for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(apply(obs.netmean,2,mean)[i],2), col = 2, lwd = 2)
legend("topleft", legend = "Observed mean density", lwd = 3, cex = 1.5, col = 2)





ppc.pvalue <- function(pred.y, m, myrep){
  pvalue <- matrix(NA, ncol = m, nrow = myrep)
  for(i in 1:myrep){
    for(k in 1:m){
      aa <- mean(pred.y[[i]][,k]>=netmean(obs.g10[[i]])[k])
      if(aa > 0.5){
        pvalue[i,k] <- 1-aa
      }else{
        pvalue[i,k] <- aa
      }
      rm(aa)
    }
  }

  return(pvalue)
}

theta1.ppc.p <- ppc.pvalue(theta1.pred.y, m, myrep )
theta2.ppc.p <- ppc.pvalue(theta2.pred.y, m, myrep )
theta3.ppc.p <- ppc.pvalue(theta3.pred.y, m, myrep )
#theta4.ppc.p <- ppc.pvalue(theta4.pred.y, m, myrep = 4 , res = "res")
#theta5.ppc.p <- ppc.pvalue(theta5.pred.y, m, myrep = 4 , res = "res")

for(i in 1:6){
  print(mean(theta1.ppc.p[,i] > 0.05))
}

dev.off()


ppc.pvalue <- function(obs.data, pred.y.data, m, myrep){
  
  p.res <- obs.mean <- matrix(NA, nrow = myrep, ncol = m)
  
  for(i in 1:myrep){
    obs.mean[i,] <- netmean(obs.data[[i]])
    for(k in 1:m){
      aa <- mean((pred.y.data[[i]][,k]>= obs.mean[i,k])*1)
      p.res[i,k] <- ifelse(aa<0.5, aa, 1-aa)
    }
  }
  return(list(p.res = p.res, obs.mean = obs.mean))
}

ppc.pvalue(obs.data = obs.g10, pred.y.data = theta1.pred.y, m, myrep)

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

sym.rate(obs.g10[[2]],m=6)
source("/Users/dongahkim/OneDrive/BNSM/bnsm_summary_function.R")

plot.sociomatrix((obs.g10[[1]][4,,]+t(obs.g10[[1]][4,,]))/2, drawlab = FALSE)




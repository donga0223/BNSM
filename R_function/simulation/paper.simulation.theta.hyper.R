#source("bnsm_function.1.R")
source("thetaprior.function.R")


library(sna)
myrep <- 50
m <-  6
V <- 50

n.group <- 3
myemp <- c(1,9,1,9,9,1)
myepp <- c(1,9,9,1,1,9)
mymodel <- "theta.hyper.rep50"

true_np <- 0.05
mymodel_np <- "obs" 

g10 <- obs.g10 <- list()
trueep10 <- trueem10 <- matrix(NA, nrow = myrep, ncol = V)


for(k in 1:myrep){
  set.seed(1200+10*k)
  tmp_ind <- sample(1:V, V*0.1)
  tmp1 <- sample(1:length(tmp_ind), round(length(tmp_ind)/2))
  d_ind <- tmp_ind[tmp1]
  d_ind2 <- tmp_ind[-tmp1]
  
  tmp <- gen.emep.dep.net(V, m, deprho = -0.5, similarity = 0.7, netprob = true_np, d_ind = d_ind, d_ind2 = d_ind2, em = myemp, ep = myepp, myseed = (1200+k))
  
  g10[[k]] <- tmp$g
  obs.g10[[k]] <- tmp$obs.g
  trueem10[k,] <- tmp$trueem
  trueep10[k,] <- tmp$trueep
  
  
  #model_np <- netmean(tmp$obs.g)
  if(mymodel_np == "true"){
    model_np <- rep(true_np,m)
  }else if(mymodel_np == "obs"){
    model_np <- netmean(tmp$obs.g)
  }
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


for(k in 1:myrep){
  for(i in 1:length(a1.em.list)){
    
    a1.em = a1.em.list[i]; b1.em = b1.em.list[i]; a2.em = 1; b2.em = 5
    a1.ep = 1; b1.ep = 1; a2.ep = 1; b2.ep = 30 ### uniform prior
    a1.theta = a1.theta.list[i]; b1.theta = b1.theta.list[i]; a2.theta = 1; b2.theta = 5 
    
    alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
    alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
    alpha.acc.theta = 0; beta.acc.theta = 0; total.alpha.theta = 0; total.beta.theta = 0
    
    assign(paste("theta",i,sep=""), bnsm.hyper.theta.hyper(dat = obs.g10[[k]], emsd = list(ema = 10, emb = 10, alpha.sd = .5, beta.sd = .5)
                                                           , epsd = list(epa = 10, epb = 10, alpha.sd = .5, beta.sd = .5)
                                                           , thetasd = list(thetaa = 3, thetab = 9, alpha.sd = .5, beta.sd = .5)
                                                           , diag = FALSE, mode = "graph", model.checking = TRUE
                                                           , reps = 3, draws = 1500, tinning = 100, burntime = 15000
                                                           , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = FALSE))
    
    
    if(i==1){
      theta1.em.density[[k]] <- density(theta1$em)
      theta1.ep.density[[k]] <- density(theta1$ep)
      theta1.em.CI[[k]] <- apply(theta1$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta1.ep.CI[[k]] <- apply(theta1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta1.b01[[k]] <- apply(theta1$net, c(3,4,5), mean)
      theta1.pred.y[[k]] <- pred.y.mean(theta1$pred.y, iter = 1500, reps = 3, m = 6)
      theta1.np[[k]] <- theta1$np
      rm(theta1)
    }else if(i==2){
      theta2.em.density[[k]] <- density(theta2$em)
      theta2.ep.density[[k]] <- density(theta2$ep)
      theta2.em.CI[[k]] <- apply(theta2$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta2.ep.CI[[k]] <- apply(theta2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta2.b01[[k]] <- apply(theta2$net, c(3,4,5), mean)
      theta2.pred.y[[k]] <- pred.y.mean(theta2$pred.y, iter = 1500, reps = 3, m = 6)
      theta2.np[[k]] <- theta2$np
      rm(theta2)
    }else if(i==3){
      theta3.em.density[[k]] <- density(theta3$em)
      theta3.ep.density[[k]] <- density(theta3$ep)
      theta3.em.CI[[k]] <- apply(theta3$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta3.ep.CI[[k]] <- apply(theta3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta3.b01[[k]] <- apply(theta3$net, c(3,4,5), mean)
      theta3.pred.y[[k]] <- pred.y.mean(theta3$pred.y, iter = 1500, reps = 3, m = 6)
      theta3.np[[k]] <- theta3$np
      rm(theta3)
    }else if(i==4){
      theta4.em.density[[k]] <- density(theta4$em)
      theta4.ep.density[[k]] <- density(theta4$ep)
      theta4.em.CI[[k]] <- apply(theta4$em, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta4.ep.CI[[k]] <- apply(theta4$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
      theta4.b01[[k]] <- apply(theta4$net, c(3,4,5), mean)
      theta4.pred.y[[k]] <- pred.y.mean(theta4$pred.y, iter = 1500, reps = 3, m = 6)
      theta4.np[[k]] <- theta4$np
      rm(theta4)
    }
    print(k)
    save.image(paste("paper.simulation.final", true_np, mymodel, "RData", sep = "."))
  }
  
}



if(2==3){
myrep <- 1
mycol <- c("red", "gray", "pink", "yellow", "green")
mylegend <- c("True" , "M.Fle-Network1", "M.Fle-Network2", "M.Fle-Network3", "M.Fle-Network4")

em.reorder <- order(trueem10[myrep,])
ep.reorder <- order(trueep10[myrep,])
trueem1 <- trueem10[myrep,em.reorder]
trueep1 <- trueep10[myrep,ep.reorder]


theta1.em <- apply(theta1$em, 2, quantile, prob = c(0.025,0.975, 0.5))[,em.reorder]
theta1.ep <- apply(theta1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))[,ep.reorder]
theta2.em <- apply(theta2$em, 2, quantile, prob = c(0.025,0.975, 0.5))[,em.reorder]
theta2.ep <- apply(theta2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))[,ep.reorder]
theta3.em <- apply(theta3$em, 2, quantile, prob = c(0.025,0.975, 0.5))[,em.reorder]
theta3.ep <- apply(theta3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))[,ep.reorder]
theta4.em <- apply(theta4$em, 2, quantile, prob = c(0.025,0.975, 0.5))[,em.reorder]
theta4.ep <- apply(theta4$ep, 2, quantile, prob = c(0.025,0.975, 0.5))[,ep.reorder]



theta1.b01 <- apply(theta1$net, c(3,4,5), mean)
theta2.b01 <- apply(theta2$net, c(3,4,5), mean)
theta3.b01 <- apply(theta3$net, c(3,4,5), mean)
theta4.b01 <- apply(theta4$net, c(3,4,5), mean)


pdf("M.Fle-Network.pdf")
par(mar = c(4.5, 5, 3.5, 1.5))
plot(c(0, 1 + V), c(0,1), type = "n", ylab = "e- Credible Interval",xlab = "node number"
     , cex.lab = 2, cex.main = 2, main = "posterior for each node")
for (i in 1:V) lines(rep(i-0.15, 2), theta1.em[1:2, i], lwd = 3, col = "gray")
for (i in 1:V) lines(rep(i-0.05, 2), theta2.em[1:2, i], lwd = 3, col = "pink")
for (i in 1:V) lines(rep(i+0.05, 2), theta3.em[1:2, i], lwd = 3, col = "yellow")
for (i in 1:V) lines(rep(i+0.15, 2), theta4.em[1:2, i], lwd = 3, col = "green")

for (i in 1:V)  lines(c(i-0.2, i+0.2), rep(trueem1[i],2), col = "red", lwd = 2)
legend("topleft", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)


plot(c(0, 1 + V), c(0,1), type = "n", ylab = "e+ Credible Interval",xlab = "node number"
     , cex.lab = 2, cex.main = 2, main = "posterior for each node")
for (i in 1:V) lines(rep(i-0.15, 2), theta1.ep[1:2, i], lwd = 3, col = "gray")
for (i in 1:V) lines(rep(i-0.05, 2), theta2.ep[1:2, i], lwd = 3, col = "pink")
for (i in 1:V) lines(rep(i+0.05, 2), theta3.ep[1:2, i], lwd = 3, col = "yellow")
for (i in 1:V) lines(rep(i+0.15, 2), theta4.ep[1:2, i], lwd = 3, col = "green")
for (i in 1:V)  lines(c(i-0.2, i+0.2), rep(trueep1[i],2), col = "red", lwd = 2)
legend("topleft", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)


plot(density(trueem10[1,]), col = "red", lwd = 2, xlab = "e-",cex.main = 2, cex.lab = 2, "posterior for all observation", ylim = c(0,120))
lines(density(theta4$em, na.rm = TRUE), col = "green", lwd = 3)
lines(density(theta2$em, na.rm = TRUE), col = "pink", lwd = 3)
lines(density(theta1$em, na.rm = TRUE), col = "gray", lwd = 3)
lines(density(theta3$em, na.rm = TRUE), col = "yellow", lwd = 3)
legend("topright", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)

plot(density(trueep10[1,]), col = "red", lwd = 2, xlab = "e+",cex.main = 2, cex.lab = 2, "posterior for all observation", ylim = c(0,8))
lines(density(theta1$ep, na.rm = TRUE), col = "gray", lwd = 3)
lines(density(theta2$ep, na.rm = TRUE), col = "pink", lwd = 3)
lines(density(theta3$ep, na.rm = TRUE), col = "yellow", lwd = 3)
lines(density(theta4$ep, na.rm = TRUE), col = "green", lwd = 3)
legend("topright", legend = mylegend, col = mycol, lwd = 3, cex = 1.5)


m <- dim(theta1$np)[2]

netnp <- array(dim=c(4500, m*4))
for(i in 1:6){
  for(k in 1:4){
    #print(k+(i-1)*3)
    netnp[,(k+(i-1)*4)] <- get(paste("theta", k, sep=""))$np[,i]
  }
}

true_np <- netmean(g10[[1]])
boxplot(netnp, col = mycol[2:5], main = "Posterior mean density", cex.main = 2, xaxt = "n")
for (i in 1:V)  lines(c(4*(i-1)+0.5, 4*i+0.5), rep(true_np[i],2), col = 2, lwd = 2)
abline(v = c(4.5, 8.5, 12.5, 16.5, 20.5), lty = 2)
axis(side=1, at=c(2.5, 6.5, 10.5, 14.5, 18.5, 22.5), mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
axis(side=1, at = 1:24, mgp=c(0,1.5,0), labels = rep(mylegend[2:5],6), tick = 0, las = 2, cex.axis = 0.5)
legend("topright", col = 2, legend = "True mean density", lty = 1, lwd = 2, cex = 1.5)

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



roc.theta1.b01 <- c(theta1.b01)
roc.theta2.b01 <- c(theta2.b01)
roc.theta3.b01 <- c(theta3.b01)
roc.theta4.b01 <- c(theta4.b01)
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

threshold <-  round(netmean(obs.g10[[myrep]]),3)
theta1.bmean <- (theta1.b01>=threshold)*1
theta2.bmean <- (theta2.b01>=threshold)*1
theta3.bmean <- (theta3.b01>=threshold)*1
theta4.bmean <- (theta4.b01>=threshold)*1


theta1.bmean.xy <- roc.table(theta1.bmean, g10[[myrep]])
theta2.bmean.xy <- roc.table(theta2.bmean, g10[[myrep]])
theta3.bmean.xy <- roc.table(theta3.bmean, g10[[myrep]])
theta4.bmean.xy <- roc.table(theta4.bmean, g10[[myrep]])


threshold <-  0.5
theta1.05 <- (theta1.b01>=threshold)*1
theta2.05 <- (theta2.b01>=threshold)*1
theta3.05 <- (theta3.b01>=threshold)*1
theta4.05 <- (theta4.b01>=threshold)*1

theta1.05.xy <- roc.table(theta1.05, g10[[myrep]])
theta2.05.xy <- roc.table(theta2.05, g10[[myrep]])
theta3.05.xy <- roc.table(theta3.05, g10[[myrep]])
theta4.05.xy <- roc.table(theta4.05, g10[[myrep]])


plot(theta1.perf, col = "gray", lwd = 2, main = paste("M.Fle-Network,", "ROC curve", sep = " "), cex.main = 2, cex.lab = 2, cex.axis = 2)
plot(theta2.perf, col = "pink", lwd = 2, add = TRUE)
plot(theta3.perf, col = "yellow", lwd = 2, add = TRUE)
plot(theta4.perf, col = "green", lwd = 2, add = TRUE)
points(uni.xy[1], uni.xy[2], col = 2, pch = 16, cex = 2)
points(inter.xy[1], inter.xy[2], col = 2, pch = 15, cex = 2)
points(theta1.bmean.xy[1], theta1.bmean.xy[2], col = mycol[1], pch = 17, cex = 2)
points(theta2.bmean.xy[1], theta2.bmean.xy[2], col = mycol[2], pch = 17, cex = 2)
points(theta3.bmean.xy[1], theta3.bmean.xy[2], col = mycol[3], pch = 17, cex = 2)
points(theta4.bmean.xy[1], theta4.bmean.xy[2], col = mycol[4], pch = 17, cex = 2)
points(theta1.05.xy[1], theta1.05.xy[2], col = mycol[1], pch = 18, cex = 2)
points(theta2.05.xy[1], theta2.05.xy[2], col = mycol[2], pch = 18, cex = 2)
points(theta3.05.xy[1], theta3.05.xy[2], col = mycol[3], pch = 18, cex = 2)
points(theta4.05.xy[1], theta4.05.xy[2], col = mycol[4], pch = 18, cex = 2)
legend("bottomright", legend = c("Union", "Intersection", expression(paste("t=", bar(Y)[k], sep = "")), "t=0.5", mylegend[2:5])
       , pch = c(16,15,17,18,20,20,20,20), col = c(2,2,1,1,mycol[2:5]), cex = 1.5)

plot(theta1.perf, col = "gray", lwd = 2, main = paste("M.Fle-Network,", "ROC curve", sep = " "), cex.main = 2, cex.lab = 2, cex.axis = 2, xlim = c(0,0.3), ylim = c(0.7,1))
plot(theta2.perf, col = "pink", lwd = 2, add = TRUE)
plot(theta3.perf, col = "yellow", lwd = 2, add = TRUE)
plot(theta4.perf, col = "green", lwd = 2, add = TRUE)
points(uni.xy[1], uni.xy[2], col = 2, pch = 16, cex = 2)
points(inter.xy[1], inter.xy[2], col = 2, pch = 15, cex = 2)
points(theta1.bmean.xy[1], theta1.bmean.xy[2], col = mycol[1], pch = 17, cex = 2)
points(theta2.bmean.xy[1], theta2.bmean.xy[2], col = mycol[2], pch = 17, cex = 2)
points(theta3.bmean.xy[1], theta3.bmean.xy[2], col = mycol[3], pch = 17, cex = 2)
points(theta4.bmean.xy[1], theta4.bmean.xy[2], col = mycol[4], pch = 17, cex = 2)
points(theta1.05.xy[1], theta1.05.xy[2], col = mycol[1], pch = 18, cex = 2)
points(theta2.05.xy[1], theta2.05.xy[2], col = mycol[2], pch = 18, cex = 2)
points(theta3.05.xy[1], theta3.05.xy[2], col = mycol[3], pch = 18, cex = 2)
points(theta4.05.xy[1], theta4.05.xy[2], col = mycol[4], pch = 18, cex = 2)
legend("bottomright", legend = c("Union", "Intersection", expression(paste("t=", bar(Y)[k], sep = "")), "t=0.5", mylegend[2:5])
       , pch = c(16,15,17,18,20,20,20,20), col = c(2,2,1,1,mycol[2:5]), cex = 1.5)


plot.sociomatrix(g10[[1]][1,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE
                 , cex.main = 2, main = paste("True,", "Network 1", sep = " ") )

plot.sociomatrix(obs.g10[[1]][1,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE
                 , cex.main = 2, main = paste("Observed,", "Network 1", sep = " ") )




dev.off()

}




### rep summary plot
if(1==2){
  
  mycol <- c("red", "gray", "pink", "yellow", "green")
  mylegend <- c("True" , "M.Fle-Network1", "M.Fle-Network2", "M.Fle-Network3", "M.Fle-Network4")
  
  
  
  em.reorder <- order(trueem10[myrep,])
  ep.reorder <- order(trueep10[myrep,])
  trueem1 <- trueem10[myrep,em.reorder]
  trueep1 <- trueep10[myrep,ep.reorder]
  
  
  
  par(mar = c(4.5, 5, 3.5, 1.5))
  p = seq(0,1, length=1000)
  plot(p, 0.95*dbeta(p, 1, 9)+0.05*dbeta(p, 9, 1), ylab="Density", type ="l", col=2, lwd = 3, xlab = "e-", main = "Posterior for all observation")
  for(i in 1:30){
    lines(theta4.em.density[[i]], col = "green")
    lines(theta2.em.density[[i]], col = "pink")
    lines(theta1.em.density[[i]], col = "gray")
    lines(theta3.em.density[[i]], col = "yellow")
  }
  
  
  plot(p, 0.95*dbeta(p, 1, 9)+0.05*dbeta(p, 9, 1), ylab="Density", type ="l", col=2, lwd = 3, xlab = "e+", main = "Posterior for all observation")
  for(i in 1:30){
    lines(theta4.ep.density[[i]], col = "green")
    lines(theta2.ep.density[[i]], col = "pink")
    lines(theta1.ep.density[[i]], col = "gray")
    lines(theta3.ep.density[[i]], col = "yellow")
  }
  
  mean.np <- matrix(NA, nrow = 30, ncol = 6)
  for(i in 1:30){
    mean.np[i,] <- apply(theta1.np[[i]],2,mean)
    
  }
  boxplot(mean.np, ylim = c(0,0.1))
  abline(h = 0.05, col = 2)
  
  
  
  
}


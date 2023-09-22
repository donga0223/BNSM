pdf("Butts_PPC_plots.pdf")
schools <- c("ARTS02", "ARTS03", "PREP02", "TECH03")
mylegend1 = c("Shared", "Not Shared" )
mn <- length(mylegend1)
mycol = c("gray", "pink")
library(ROCR)
library(sna)


for(school in schools){
  print(school)
  load(paste("Butts_PPC_", school, ".RData", sep = ""))
  load(paste("Butts_PPC_one_", school, ".RData", sep=""))
  source("thetaprior_function_forpaper.R")
  
  m=6
  
  
  Butts_all <- pred.y.mean(Butts_all$pred.y, iter = 1500, reps = 3, m = 6) ## (4500,6)
  Butts_Net1 <- pred.y.mean(m1$pred.y, iter = 1500, reps = 3, m = 1) ## (4500,1)
  Butts_Net2 <- pred.y.mean(m2$pred.y, iter = 1500, reps = 3, m = 1) ## (4500,1)
  Butts_Net3 <- pred.y.mean(m3$pred.y, iter = 1500, reps = 3, m = 1) ## (4500,1)
  Butts_Net4 <- pred.y.mean(m4$pred.y, iter = 1500, reps = 3, m = 1) ## (4500,1)
  Butts_Net5 <- pred.y.mean(m5$pred.y, iter = 1500, reps = 3, m = 1) ## (4500,1)
  Butts_Net6 <- pred.y.mean(m6$pred.y, iter = 1500, reps = 3, m = 1) ## (4500,1)
  
  
  netnp <- array(dim=c(4500, m*2))
  netnp[, c(1,3,5,7,9,11)] <- Butts_all
  netnp[, 2] <- Butts_Net1
  netnp[, 4] <- Butts_Net2
  netnp[, 6] <- Butts_Net3
  netnp[, 8] <- Butts_Net4
  netnp[, 10] <- Butts_Net5
  netnp[, 12] <- Butts_Net6
  
  V <- dim(get(school))[2]
  
  
  boxplot(netnp, col = mycol, main = paste(school, " Posterior Predictive check \n From Butts methods"), cex.main = 2, xaxt = "n")
  #for (i in 1:V)  lines(c(4*(i-1)+0.5, 4*i+0.5), rep(true_np[i],2), col = 2, lwd = 2)
  #abline(v = c(4.5, 8.5, 12.5, 16.5, 20.5), lty = 2)
  abline(v = (1:m)*mn+0.5, lty = 2)
  
  #axis(side=1, at=c(2.5, 6.5, 10.5, 14.5, 18.5, 22.5), mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  #axis(side=1, at = 1:24, mgp=c(0,1.5,0), labels = rep(paste0("m",1:4),6), tick = 0)
  axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  #axis(side=1, at = 1:(mn*m), mgp=c(0,1.5,0), labels = rep(paste0("m",1:mn),6), tick = 0)
  
  for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(netmean(get(school))[i],2), col = 2, lwd = 2)
  legend("bottomright", legend = "observed density", lwd = 3, cex = 1.5, col = 2)
  legend("topleft", legend = mylegend1, col = mycol, pch = 19, cex = 1.5, horiz = TRUE)
}

dev.off()
#save.image(file = "Butts_PPC_plots.RData")



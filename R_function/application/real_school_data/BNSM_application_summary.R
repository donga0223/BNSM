rm(list=ls())
source("R_function/BNSM_summary_function.R")
library(sna)
library(RColorBrewer)

schools <- c("ARTS02", "ARTS03", "PREP02", "TECH03")
for(school in schools){
  load(paste("RDatas/application/BNSM_application_", school, ".RData", sep=""))
  
  pdfname <- paste("figures/BNSM_application_", school, ".pdf", sep="")
  
  #mylegend1 = c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9")
  mylegend1 = c("m1", "m2", "m3", "m4", "m9")
  
  mn <- length(mylegend1)
  mycol =  brewer.pal(mn, "Set1")
  
  V <- dim(get(school))[2]
  m <- dim(get(school))[1]
  
  
  pdf(pdfname)
  par(mar = c(4.5, 5, 3.5, 1.5))
  plot(c(0, 1 + V), c(0,1), type = "n", ylab = "e- Credible Interval",xlab = "node number"
       , cex.lab = 2, cex.main = 2, cex.axis = 2, main = paste(school, "Posterior for each node"))
  for (i in 1:V) lines(rep(i-0.2, 2), m1.res[[3]][1:2, i], lwd = 3, col = mycol[1])
  for (i in 1:V) lines(rep(i-0.1, 2), m2.res[[3]][1:2, i], lwd = 3, col = mycol[2])
  for (i in 1:V) lines(rep(i, 2), m3.res[[3]][1:2, i], lwd = 3, col = mycol[3])
  for (i in 1:V) lines(rep(i+.1, 2), m4.res[[3]][1:2, i], lwd = 3, col = mycol[4])
  #for (i in 1:V) lines(rep(i+.1, 2), m5.res[[3]][1:2, i], lwd = 3, col = mycol[5])
  #for (i in 1:V) lines(rep(i+.1, 2), m6.res[[3]][1:2, i], lwd = 3, col = mycol[6])
  #for (i in 1:V) lines(rep(i+.1, 2), m7.res[[3]][1:2, i], lwd = 3, col = mycol[7])
  #for (i in 1:V) lines(rep(i+.1, 2), m8.res[[3]][1:2, i], lwd = 3, col = mycol[8])
  for (i in 1:V) lines(rep(i+.1, 2), m9.res[[3]][1:2, i], lwd = 3, col = mycol[5])
  legend("topleft", legend = mylegend1, col = mycol, lwd = 3, cex = 2)
  
  
  plot(c(0, 1 + V), c(0,1), type = "n", ylab = "e+ Credible Interval",xlab = "node number"
       , cex.lab = 2, cex.main = 2, cex.axis = 2, main = paste(school, "Posterior for each node"))
  for (i in 1:V) lines(rep(i-0.2, 2), m1.res[[4]][1:2, i], lwd = 3, col = mycol[1])
  for (i in 1:V) lines(rep(i-0.1, 2), m2.res[[4]][1:2, i], lwd = 3, col = mycol[2])
  for (i in 1:V) lines(rep(i, 2), m3.res[[4]][1:2, i], lwd = 3, col = mycol[3])
  for (i in 1:V) lines(rep(i+.1, 2), m4.res[[4]][1:2, i], lwd = 3, col = mycol[4])
  #for (i in 1:V) lines(rep(i+.2, 2), m5.res[[4]][1:2, i], lwd = 3, col =mycol[5])
  #for (i in 1:V) lines(rep(i+.2, 2), m6.res[[4]][1:2, i], lwd = 3, col =mycol[6])
  #for (i in 1:V) lines(rep(i+.2, 2), m7.res[[4]][1:2, i], lwd = 3, col =mycol[7])
  #for (i in 1:V) lines(rep(i+.2, 2), m8.res[[4]][1:2, i], lwd = 3, col =mycol[8])
  for (i in 1:V) lines(rep(i+.2, 2), m9.res[[4]][1:2, i], lwd = 3, col =mycol[5])
  legend("topleft", legend = mylegend1, col = mycol, lwd = 3, cex = 2)
  
  ################################################################################
  ymax <- max(c(m1.res[[1]]$y, m2.res[[1]]$y, m3.res[[1]]$y, m4.res[[1]]$y, m9.res[[1]]$y) )
  plot(m1.res[[1]], col = mycol[1], lwd = 3, xlab = "e-"
       , cex.axis = 2, cex.main = 2, cex.lab = 2, main = paste(school, ", Posterior for all"), ylim = c(0,ymax))
  lines(m2.res[[1]], col = mycol[2], lwd = 3, lty = 2)
  lines(m3.res[[1]], col = mycol[3], lwd = 3)
  lines(m4.res[[1]], col = mycol[4], lwd = 3, lty = 2)
  #lines(m5.res[[1]], col = mycol[5], lwd = 3)
  #lines(m6.res[[1]], col = mycol[6], lwd = 1)
  #lines(m7.res[[1]], col = mycol[7], lwd = 3)
  #lines(m8.res[[1]], col = mycol[8], lwd = 1)
  lines(m9.res[[1]], col = mycol[5], lwd = 3)
  legend("topright", legend = mylegend1, col = mycol, lwd = 3, lty = c(1,2,1,2,1), cex = 2)
  
  
  ymax <- max(c(m1.res[[2]]$y, m2.res[[2]]$y, m3.res[[2]]$y, m4.res[[2]]$y, m9.res[[2]]$y) )
  
  plot(m9.res[[2]], col = mycol[5], lwd = 3, xlab = "e+"
       , cex.axis = 2, cex.main = 2, cex.lab = 2, main = paste(school, ", Posterior for all", sep = ""), ylim = c(0, ymax))
  lines(m1.res[[2]], col = mycol[1], lwd = 3)
  lines(m2.res[[2]], col = mycol[2], lwd = 3, lty = 2)
  lines(m3.res[[2]], col = mycol[3], lwd = 3)
  lines(m4.res[[2]], col = mycol[4], lwd = 3)
  #lines(m6.res[[2]], col = mycol[6], lwd = 1)
  #lines(m7.res[[2]], col = mycol[7], lwd = 3)
  #lines(m8.res[[2]], col = mycol[8], lwd = 1)
  #lines(m4.res[[2]], col = mycol[4], lwd = 3, lty = 2)
  legend("topright", legend = mylegend1, col = mycol, lwd = 3, lty = c(1,2,1,2,1), cex = 2)
  
  
  m5.res <- m9.res
  netnp <- array(dim=c(4500, mn*m))
  for(i in 1:m){
    for(k in c(1:mn)){
      netnp[,(k+(i-1)*mn)] <- get(paste("m", k,".res", sep=""))[[6]][,i]
    }
  }
  
  
  
  boxplot(netnp, col = mycol, main = paste(school, ", Posterior for Theta", sep=""), cex.main = 2, xaxt = "n")
  abline(v = (1:m)*mn+0.5, lty = 2)
  axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  
  for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(netmean(get(school))[i],2), col = 2, lwd = 2)
  legend("bottomright", legend = "observed density", lwd = 3, cex = 1.5, col = 2)
  legend("topleft", legend = mylegend1, col = mycol, pch = 19, cex = 1.5, horiz = TRUE)
  
  plot.sociomatrix(get(school)[1,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE
                   , main = paste(school, " Observed Network (", round(netmean(get(school))[1],4), ")", sep = ""), cex.main = 2)
  
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
  
  union.net <- union.intersection.net(get(school))$union.net
  intersection.net <- union.intersection.net(get(school))$intersection.net
  
  plot.sociomatrix(union.net[1,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE
                   , main = paste(school, " Union (", round(netmean(union.net)[1],4), ")", sep = ""), cex.main = 2)
  plot.sociomatrix(intersection.net[1,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE
                   , main = paste(school, " Intersection (", round(netmean(intersection.net)[1],4), ")", sep = ""), cex.main = 2)
  
  m1.res[[5]][which(is.na(m1.res[[5]]))] <- 0
  plot.sociomatrix(m1.res[[5]][1,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE, main = paste(school, "f2f"), cex.main = 2)
  plot.sociomatrix(m1.res[[5]][2,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE, main = paste(school, "online"), cex.main = 2)
  plot.sociomatrix(m1.res[[5]][3,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE, main = paste(school, "eatout"), cex.main = 2)
  plot.sociomatrix(m1.res[[5]][4,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE, main = paste(school, "pa"), cex.main = 2)
  plot.sociomatrix(m1.res[[5]][5,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE, main = paste(school, "screen"), cex.main = 2)
  plot.sociomatrix(m1.res[[5]][6,,], drawlab = FALSE, diaglab = FALSE, drawlines = FALSE, main = paste(school, "lunch"), cex.main = 2)
  
  netnp <- array(dim=c(4500, mn*m))
  for(i in 1:m){
    for(k in c(1:mn)){
      netnp[,(k+(i-1)*mn)] <- get(paste("m", k,".res", sep=""))[[7]][,i]
    }
  }
  
  
  boxplot(netnp, col = mycol, main = paste(school, " Posterior Predictive check",sep = ""), cex.main = 2, xaxt = "n")
  abline(v = (1:m)*mn+0.5, lty = 2)
  axis(side=1, at=(1:m)*mn-mn/2, mgp=c(0,0.5,0), labels = paste0("Net",1:6), cex.axis = 1.5, tick = 0)
  
  for (i in 1:V)  lines(c(mn*(i-1)+0.5, mn*i+0.5), rep(netmean(get(school))[i],2), col = 2, lwd = 2)
  legend("bottomright", legend = "observed density", lwd = 3, cex = 1.5, col = 2)
  legend("topleft", legend = mylegend1, col = mycol, pch = 19, cex = 1.5, horiz = TRUE)
  
  
  dev.off()
  
  
  
  
}


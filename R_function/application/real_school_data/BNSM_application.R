args <- commandArgs(trailingOnly = TRUE)
school <- args[1]


load("application/net.symmetric.application.RData")
source("thetaprior_function_forpaper.R")

m1.res <- m2.res <- m3.res <- m4.res <- m5.res <- list()
m6.res <- m7.res <- m8.res <- m9.res <- list()

a1.em.list <- c(3,3,1,1,3,3,1,1,1)
b1.em.list <- c(5,5,9,9,5,5,9,9,1)
b2.em.list <- c(5,5,5,5,5,5,5,5,30)

a1.ep.list <- c(rep(1,9))
b1.ep.list <- c(rep(1,4), rep(9,4), 1)
b2.ep.list <- c(rep(30,4), rep(5,4), 30)

a1.theta.list <- c(3,1,3,1,3,1,3,1,1)
b1.theta.list <- c(5,9,5,9,5,9,5,9,1)
b2.theta.list <- c(5,5,5,5,5,5,5,5,30)

load(paste("application/BNSM_application_", school, ".RData", sep=""))
#for(k in 1:length(a1.em.list)){
for(k in c(8,9)){
  
  a1.em = a1.em.list[k]; b1.em = b1.em.list[k]; a2.em = 1; b2.em = b2.em.list[k]
  #a1.ep = 1; b1.ep = 1; a2.ep = 1; b2.ep = 30 ### uniform prior
  a1.ep = a1.ep.list[k]; b1.ep = b1.ep.list[k]; a2.ep = 1; b2.ep = b2.ep.list[k] 
  a1.theta = a1.theta.list[k]; b1.theta = b1.theta.list[k]; a2.theta = 1; b2.theta = b2.theta.list[k]
  
  alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
  alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
  alpha.acc.theta = 0; beta.acc.theta = 0; total.alpha.theta = 0; total.beta.theta = 0
  
  assign(paste("m",k,sep=""), bnsm.hyper.theta.hyper(dat = get(school), emsd = list(ema = 10, emb = 10, alpha.sd = .5, beta.sd = .5)
                                                          , epsd = list(epa = 10, epb = 10, alpha.sd = .5, beta.sd = .5)
                                                            , thetasd = list(thetaa = 3, thetab = 9, alpha.sd = .5, beta.sd = .5)
                                                            , diag = FALSE, mode = "graph", model.checking = TRUE
                                                            , reps = 3, draws = 1500, tinning = 100, burntime = 15000
                                                            , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = FALSE))
  
  print(k)
  if(k==1){
    m1.res[[1]] <- density(m1$em)
    m1.res[[2]] <- density(m1$ep)
    m1.res[[3]] <- apply(m1$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m1.res[[4]] <- apply(m1$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m1.res[[5]] <- apply(m1$net, c(3,4,5), mean)
    m1.res[[6]] <- m1$np
    m1.res[[7]] <- pred.y.mean(m1$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m1)
  }else if(k==2){
    m2.res[[1]] <- density(m2$em)
    m2.res[[2]] <- density(m2$ep)
    m2.res[[3]] <- apply(m2$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m2.res[[4]] <- apply(m2$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m2.res[[5]] <- apply(m2$net, c(3,4,5), mean)
    m2.res[[6]] <- m2$np
    m2.res[[7]] <- pred.y.mean(m2$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m2)
  }else if(k==3){
    m3.res[[1]] <- density(m3$em)
    m3.res[[2]] <- density(m3$ep)
    m3.res[[3]] <- apply(m3$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m3.res[[4]] <- apply(m3$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m3.res[[5]] <- apply(m3$net, c(3,4,5), mean)
    m3.res[[6]] <- m3$np
    m3.res[[7]] <- pred.y.mean(m3$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m3)
  }else if(k==4){
    m4.res[[1]] <- density(m4$em)
    m4.res[[2]] <- density(m4$ep)
    m4.res[[3]] <- apply(m4$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m4.res[[4]] <- apply(m4$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m4.res[[5]] <- apply(m4$net, c(3,4,5), mean)
    m4.res[[6]] <- m4$np
    m4.res[[7]] <- pred.y.mean(m4$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m4)
  }else if(k==5){
    m5.res[[1]] <- density(m5$em)
    m5.res[[2]] <- density(m5$ep)
    m5.res[[3]] <- apply(m5$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m5.res[[4]] <- apply(m5$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m5.res[[5]] <- apply(m5$net, c(3,4,5), mean)
    m5.res[[6]] <- m5$np
    m5.res[[7]] <- pred.y.mean(m5$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m5)
  }
    if(k==6){
    m6.res[[1]] <- density(m6$em)
    m6.res[[2]] <- density(m6$ep)
    m6.res[[3]] <- apply(m6$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m6.res[[4]] <- apply(m6$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m6.res[[5]] <- apply(m6$net, c(3,4,5), mean)
    m6.res[[6]] <- m6$np
    m6.res[[7]] <- pred.y.mean(m6$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m6)
  }else if(k==7){
    m7.res[[1]] <- density(m7$em)
    m7.res[[2]] <- density(m7$ep)
    m7.res[[3]] <- apply(m7$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m7.res[[4]] <- apply(m7$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m7.res[[5]] <- apply(m7$net, c(3,4,5), mean)
    m7.res[[6]] <- m7$np
    m7.res[[7]] <- pred.y.mean(m7$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m7)
  }else if(k==8){
    m8.res[[1]] <- density(m8$em)
    m8.res[[2]] <- density(m8$ep)
    m8.res[[3]] <- apply(m8$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m8.res[[4]] <- apply(m8$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m8.res[[5]] <- apply(m8$net, c(3,4,5), mean)
    m8.res[[6]] <- m8$np
    m8.res[[7]] <- pred.y.mean(m8$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m8)
  }else if(k==9){
    m9.res[[1]] <- density(m9$em)
    m9.res[[2]] <- density(m9$ep)
    m9.res[[3]] <- apply(m9$em, 2, quantile, prob = c(0.025,0.975, 0.5))
    m9.res[[4]] <- apply(m9$ep, 2, quantile, prob = c(0.025,0.975, 0.5))
    m9.res[[5]] <- apply(m9$net, c(3,4,5), mean)
    m9.res[[6]] <- m9$np
    m9.res[[7]] <- pred.y.mean(m9$pred.y, iter = 1500, reps = 3, m = 6)
    rm(m9)
  }

  save.image(file = paste("application/BNSM_application_", school, ".RData", sep = ""))
}



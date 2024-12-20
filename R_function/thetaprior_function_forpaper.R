if(1==2){
  dat = tmp$obs.g
  np = model_np
  emsd = list(ema = 10, emb = 10, alpha.sd = 1, beta.sd = 1)
  epsd = list(epa = 10, epb = 10, alpha.sd = 1, beta.sd = 1)
  diag = FALSE
  mode = "graph"
  reps = 3
  draws = 100
  tinning = 1
  burntime = 150
  quiet = TRUE
  anames = NULL
  onames = NULL
  compute.sqrtrhat = FALSE
  
}

library(sna)
library(copula)
#library(MCMCpack)

gen.emep.dep.net <- function(V, m, deprho, similarity = NULL, netprob, d_ind = NULL, d_ind2 = NULL, em, ep, myseed){
  set.seed(myseed)
  if(is.null(similarity)){
    g <- rgraph(V, m = m, mode = "graph", tprob = netprob)
  }else{
    if(length(netprob) == 1){
      g1 <- rgraph(V, m = 1, mode = "graph", tprob = netprob) 
      g.p <- similarity*g1+(1-similarity)*(1-g1)
      g.p <- g.p*(netprob/mean(g.p))
      g <- rgraph(V, m = m, mode = "graph", tprob = g.p) 
    }else{
      g <- array(NA, dim = c(m,V,V))
      g1 <- rgraph(V, m = 1, mode = "graph", tprob = mean(netprob)) 
      g.p <- similarity*g1+(1-similarity)*(1-g1)
      for(h in 1:length(netprob)){
        g.p11 <- g.p*(netprob[h]/mean(g.p))
        g[h,,] <- rgraph(V, m = 1, mode = "graph", tprob = g.p11) 
      }
    }
  }
  
  myCop <- normalCopula(param=c(deprho), dim = 2, dispstr = "un")
  myMvd1 <- mvdc(copula=myCop, margins=c("beta", "beta"),
                 paramMargins=list(list(shape1=em[1], shape2=em[2]), list(shape1=ep[1], shape2=ep[2])))
  emep1 <- rMvdc(V, myMvd1)
  trueem <- emep1[,1]
  trueep <- emep1[,2]
  
  if(is.null(d_ind) == FALSE){
    myMvd2 <- mvdc(copula=myCop, margins=c("beta", "beta"),
                   paramMargins=list(list(shape1=em[3], shape2=em[4]), list(shape1=ep[3], shape2=ep[4])))
    emep2 <- rMvdc(length(d_ind), myMvd2)
    trueem[d_ind] <- emep2[,1] 
    trueep[d_ind] <- emep2[,2] 
  }
  
  if(is.null(d_ind2) == FALSE){
    myMvd3 <- mvdc(copula=myCop, margins=c("beta", "beta"),
                   paramMargins=list(list(shape1=em[5], shape2=em[6]), list(shape1=ep[5], shape2=ep[6])))
    emep3 <- rMvdc(length(d_ind2), myMvd3)
    trueem[d_ind2] <- emep3[,1] ## all
    trueep[d_ind2] <- emep3[,2] ## all
    
  }
  
  g.pp <- aperm(g, c(3,1,2))*(1-trueem)+(1-aperm(g, c(3,1,2)))*trueep
  g.p <- aperm(g.pp, c(2,1,3))
  
  obs.g <- array(NA, dim = c(m, V, V))
  for(i in 1:m){
    obs.g[i,,] <- rgraph(V, mode = "digraph", tprob = g.p[i,,], diag = FALSE)  
  }
  return(list(g = g, obs.g = obs.g, trueem = trueem, trueep = trueep))
} 


draw.alpha <- function(alpha,beta,theta,prop.sd, type) {
  V <- length(theta)
  alpha.star <- rnorm(1, alpha, prop.sd)
  num <- V*(lgamma(alpha.star+beta) - lgamma(alpha.star)) +alpha.star*sum(log(theta)) + log.prior(alpha.star,beta, type)
  den <- V*(lgamma(alpha+beta)- lgamma(alpha)) +alpha*sum(log(theta)) + log.prior(alpha,beta, type)
  acc <- ifelse((log(runif(1))<=num - den)&&(alpha.star>0)&&(!is.na(num-den))&&(alpha.star <= beta),1,0)
  if(type == "em"){
    alpha.acc.em <<- alpha.acc.em+acc
    total.alpha.em <<- total.alpha.em + 1
  }else if(type == "ep"){
    alpha.acc.ep <<- alpha.acc.ep+acc
    total.alpha.ep <<- total.alpha.ep + 1
  }else if(type == "theta"){
    alpha.acc.theta <<- alpha.acc.theta+acc
    total.alpha.theta <<- total.alpha.theta + 1
  }
  return(ifelse(acc,alpha.star,alpha))
}

draw.beta <- function(alpha,beta,theta,prop.sd, type) {
  V <- length(theta)
  beta.star <- rnorm(1, beta, prop.sd)
  num <- V*(lgamma(alpha+beta.star) - lgamma(beta.star)) +beta.star*sum(log(1-theta)) + log.prior(alpha,beta.star, type)
  den <- V*(lgamma(alpha+beta)      - lgamma(beta)) +beta     *sum(log(1-theta)) + log.prior(alpha,beta, type)
  acc <- ifelse((log(runif(1))<=num - den)&&(beta.star>0)&&(!is.na(num-den))&&(alpha <= beta.star),1,0)
  
  if(type == "em"){
    beta.acc.em <<- beta.acc.em+acc
    total.beta.em <<- total.beta.em + 1
  }else if(type == "ep"){
    beta.acc.ep <<- beta.acc.ep+acc
    total.beta.ep <<- total.beta.ep + 1
  }else if(type == "theta"){
    beta.acc.theta <<- beta.acc.theta+acc
    total.beta.theta <<- total.beta.theta + 1
  }
  return(ifelse(acc,beta.star,beta))
}

log.prior <- function(alpha, beta, type){
  if(type == "em"){
    a1 = a1.em; b1 = b1.em; a2 = a2.em; b2 = b2.em
    res <- (a1-1)*(log(alpha)-log(alpha+beta))+(b1-1)*(log(beta)-log(alpha+beta))-(a2-1)*log(alpha+beta+1)+(b2-1)*(log(alpha+beta)-log(alpha+beta+1))-log(alpha+beta)-2*log(alpha+beta+1)
  }else if(type == "ep"){
    a1 = a1.ep; b1 = b1.ep; a2 = a2.ep; b2 = b2.ep
    res <- (a1-1)*(log(alpha)-log(alpha+beta))+(b1-1)*(log(beta)-log(alpha+beta))-(a2-1)*log(alpha+beta+1)+(b2-1)*(log(alpha+beta)-log(alpha+beta+1))-log(alpha+beta)-2*log(alpha+beta+1)
  }else if(type == "theta"){
    a1 = a1.theta; b1 = b1.theta; a2 = a2.theta; b2 = b2.theta
    res <- -log(alpha+beta)-2*log(alpha+beta+1)
  }
  return(res)
}




sampletie_func <- function(ep, em, m, n, d, mode, diag, np, model.checking){
  ep.a <- em.a <- array(NA, dim = c(m,n,n))
  for(x in 1:n){
    ep.a[,x,] <- ep[x]
    em.a[,x,] <- em[x]
  }
  pygt <- d * (1 - em.a) + (1 - d) * em.a
  pygnt <- d * ep.a + (1 - d) * (1 - ep.a)
  pygt[which(is.na(pygt))] <- 1
  pygnt[which(is.na(pygnt))] <- 1
  a <- tieprob <- pred_y <- yprob <- array(NA, dim = c(m,n,n))
  for(l in 1:m){
    num <- np[l] * (pygt[l,,]*t(pygt[l,,]))
    den1 <- (1 - np[l]) * (pygnt[l,,]*t(pygnt[l,,]))
    tieprob[l,,] <- num/(num + den1)
  }
  a <- rgraph(n, m, tprob = tieprob, mode = mode, diag = diag)
  if (!diag) 
    a <- diag.remove(a)
  if (length(dim(a)) == 2){
    a <- array(a, dim = c(1, NROW(a), NCOL(a)))
  }
  
  if(model.checking == TRUE){
    yprob <- a*(1-em.a)+(1-a)*ep.a
    pred_y <- rgraph(n, m, tprob = yprob, mode = "digraph", diag = diag)
    pred_y[which(is.na(pygt))] <- NA
    if (!diag) 
      pred_y <- diag.remove(pred_y)
    return(list(a = a, pred_y = pred_y))
  }else if(model.checking == FALSE){
    return(a)
  }
}
 
netmean <- function(g){
  m <- dim(g)[1]
  V <- dim(g)[2]
  res <- apply(g,1,sum, na.rm = TRUE)/(V^2-V)
  return(res)
}




pred.y.mean <- function(pred.y, iter, reps, m){
  res <- matrix(NA, nrow = iter*reps, ncol = m)
  if(m>1){
    for(j in 1:reps){
      for(i in 1:iter){
        res[((j-1)*iter+i),] <- netmean(pred.y[j,i,,,])
      }
    }
  }else if(m==1){
    for(j in 1:reps){
      for(i in 1:iter){
        res[((j-1)*iter+i),] <- mean(pred.y[j,i,,,], na.rm = TRUE)
      }
    }
  }
  return(res)
}

bnsm.hyper.theta.hyper <- function(dat, emsd = list(ema = 10, emb = 10, alpha.sd = 0.25, beta.sd = 0.25)
                                   , epsd = list(epa = 3, epb = 9, alpha.sd = 0.25, beta.sd = 0.25)
                                   , thetasd = list(thetaa = 3, thetab = 9, alpha.sd = 0.25, beta.sd = 0.25)
                                   , diag = FALSE, mode = "graph", model.checking = TRUE
                                   , reps = 3, draws = 1500, tinning = 10, burntime = 500
                                   , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE) {
  
  start.time <- Sys.time()
  dat <- as.sociomatrix.sna(dat, simplify = TRUE)
  
  if (is.list(dat)) 
    stop("All bbnam input graphs must be of the same order.")
  if (length(dim(dat)) == 2) 
    dat <- array(dat, dim = c(1, NROW(dat), NCOL(dat)))
  if (reps == 1) 
    compute.sqrtrhat <- FALSE
  m <- dim(dat)[1]
  n <- dim(dat)[2]
  d <- dat
  slen <- burntime + draws*tinning
  n.draws <- burntime+tinning*(1:draws)
  
  out <- list()
  tieprob <- array(NA, dim=c(m,n,n))
  
  if (is.null(anames)) 
    anames <- paste("a", 1:n, sep = "")
  if (is.null(onames)) 
    onames <- paste("o", 1:m, sep = "")
  
  if (!diag) 
    d <- diag.remove(d)
  if (!quiet) 
    cat("Creating temporary variables and drawing initial conditions....\n")
  
  
  
  res.a <- array(dim = c(reps, draws, m, n, n))
  res.pred.y <- array(dim = c(reps, draws, m, n, n))
  res.em <- array(dim = c(reps, draws, n))
  res.ep <- array(dim = c(reps, draws, n))
  res.np <- array(dim = c(reps, draws, m))
  #res.np <- array(dim = c(reps, draws, n))
  res.emalpha <- array(dim = c(reps,draws))
  res.embeta <- array(dim = c(reps,draws))
  res.epalpha <- array(dim = c(reps,draws))
  res.epbeta <- array(dim = c(reps,draws))
  res.thetaalpha <- array(dim = c(reps,draws))
  res.thetabeta <- array(dim = c(reps,draws))
  
  res.alpha.acc.em <- rep(NA, reps)
  res.alpha.acc.ep <- rep(NA, reps)
  res.alpha.acc.theta <- rep(NA, reps)
  res.beta.acc.em <- rep(NA, reps)
  res.beta.acc.ep <- rep(NA, reps)
  res.beta.acc.theta <- rep(NA, reps)
  res.total.alpha.em <- rep(NA, reps)
  alpha.acc.rate.em <- rep(NA, reps)
  alpha.acc.rate.ep <- rep(NA, reps)
  alpha.acc.rate.theta <- rep(NA, reps)
  beta.acc.rate.em <- rep(NA, reps)
  beta.acc.rate.ep <- rep(NA, reps)
  beta.acc.rate.theta <- rep(NA, reps)
  
  for(j in 1:reps){
    emalpha <- emsd$ema
    embeta <- emsd$emb
    epalpha <- epsd$epa
    epbeta <- epsd$epb
    
    em <- rbeta(n, emalpha, embeta)
    ep <- rbeta(n, epalpha, epbeta)
    
    thetaalpha <- thetasd$thetaa
    thetabeta <- thetasd$thetab
    
    np <- rbeta(m, thetaalpha, thetabeta)
    #np <- rbeta(n, thetaalpha, thetabeta)
    
    ### updated values setting ######
    for(i in 2:slen){
      if(model.checking == TRUE){
        sample.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
        new.a <- sample.a$a
        pred.y <- sample.a$pred_y
      }else{
        new.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
      }
      
      cem <- matrix(nrow = n, ncol = 2)
      cep <- matrix(nrow = n, ncol = 2)
      
      cem[, 1] <- rowSums(colSums((1 - d) * new.a,  na.rm = TRUE))
      cem[, 2] <- rowSums(colSums(d * new.a,  na.rm = TRUE))
      cep[, 1] <- rowSums(colSums(d * (1-new.a),  na.rm = TRUE))
      cep[, 2] <- rowSums(colSums((1 - d) * (1-new.a),  na.rm = TRUE))
      
      new.emalpha <- draw.alpha(emalpha, embeta, theta = em, prop.sd = emsd$alpha.sd, type = "em")
      new.embeta <- draw.beta(new.emalpha, embeta, theta = em, prop.sd = emsd$beta.sd, type = "em")
      
      new.epalpha <- draw.alpha(epalpha, epbeta, theta = ep, prop.sd = epsd$alpha.sd, type = "ep")
      new.epbeta <- draw.beta(new.epalpha, epbeta, theta = ep, prop.sd = epsd$beta.sd, type = "ep")
      
      new.em <- rbeta(n, new.emalpha + cem[, 1], new.embeta + cem[, 2])
      new.ep <- rbeta(n, new.epalpha + cep[, 1], new.epbeta + cep[, 2])
      
      new.thetaalpha <- draw.alpha(thetaalpha, thetabeta, theta = np, prop.sd = thetasd$alpha.sd, type = "theta")
      new.thetabeta <- draw.beta(new.thetaalpha, thetabeta, theta = np, prop.sd = thetasd$beta.sd, type = "theta")
      
      new.np <- rbeta(m, rowSums(new.a, na.rm = TRUE)+new.thetaalpha, rowSums(1-new.a, na.rm = TRUE)+new.thetabeta)
      #new.np <- rbeta(n, rowSums(colSums(new.a,  na.rm = TRUE))+new.thetaalpha, rowSums(colSums(1-new.a,  na.rm = TRUE))+new.thetabeta)
      
      #### storing values after burnning and every tinning
      if(i %in% n.draws){
        k <- (i-burntime)/tinning
        res.a[j,k,,,] <- new.a
        res.pred.y[j,k,,,] <- pred.y
        res.emalpha[j,k] <- new.emalpha
        res.embeta[j,k] <- new.embeta
        res.epalpha[j,k] <- new.epalpha
        res.epbeta[j,k] <- new.epbeta
        res.thetaalpha[j,k] <- new.thetaalpha
        res.thetabeta[j,k] <- new.thetabeta
        
        res.em[j,k,] <- new.em
        res.ep[j,k,] <- new.ep
        res.np[j,k,] <- new.np
      }
      
      rm(emalpha, embeta, epalpha, epbeta, thetaalpha, thetabeta, em, ep, new.a, np, pred.y)
      emalpha <- new.emalpha
      embeta <- new.embeta
      epalpha <- new.epalpha
      epbeta <- new.epbeta
      thetaalpha <- new.thetaalpha
      thetabeta <- new.thetabeta
      
      em <- new.em
      ep <- new.ep
      np <- new.np
      rm(new.emalpha, new.embeta, new.epalpha, new.epbeta, new.thetaalpha, new.thetabeta, new.em, new.ep, new.np)
    }
    print(alpha.acc.em)
    print(total.alpha.em)
    print(alpha.acc.em/total.alpha.em)
    res.alpha.acc.em[j] <- alpha.acc.em
    res.alpha.acc.ep[j] <- alpha.acc.ep
    res.alpha.acc.theta[j] <- alpha.acc.theta
    res.beta.acc.em[j] <- beta.acc.em
    res.beta.acc.ep[j] <- beta.acc.ep
    res.beta.acc.theta[j] <- beta.acc.theta
    res.total.alpha.em[j] <- total.alpha.em
  }
  
  alpha.acc.rate.em[1] <- res.alpha.acc.em[1]/total.alpha.em[1]
  alpha.acc.rate.ep[1] <- res.alpha.acc.ep[1]/total.alpha.em[1]
  alpha.acc.rate.theta[1] <- res.alpha.acc.theta[1]/total.alpha.em[1]
  beta.acc.rate.em[1] <- res.beta.acc.em[1]/total.alpha.em[1]
  beta.acc.rate.ep[1] <- res.beta.acc.ep[1]/total.alpha.em[1]
  beta.acc.rate.theta[1] <- res.beta.acc.theta[1]/total.alpha.em[1]
  for(j in 2:reps){
    alpha.acc.rate.em[j] <- (res.alpha.acc.em[j]-res.alpha.acc.em[(j-1)])/(total.alpha.em[j]-total.alpha.em[(j-1)])
    alpha.acc.rate.ep[j] <- (res.alpha.acc.ep[j]-res.alpha.acc.ep[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    alpha.acc.rate.theta[j] <- (res.alpha.acc.theta[j]-res.alpha.acc.theta[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    
    beta.acc.rate.em[j] <- (res.beta.acc.em[j]-res.beta.acc.em[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    beta.acc.rate.ep[j] <- (res.beta.acc.ep[j]-res.beta.acc.ep[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    beta.acc.rate.theta[j] <- (res.beta.acc.theta[j]-res.beta.acc.theta[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
  }
  out$net <- res.a 
  out$pred.y <- res.pred.y
  out$em <- res.em[1, , ]
  out$ep <- res.ep[1, , ]
  out$np  <- res.np[1, ,]
  out$emalpha <- res.emalpha[1, ]
  out$embeta <- res.embeta[1, ]
  out$epalpha <- res.epalpha[1, ]
  out$epbeta  <- res.epbeta[1, ]
  out$thetaalpha <- res.thetaalpha[1, ]
  out$thetabeta  <- res.thetabeta[1, ]
  
  
  if (reps >= 2) 
    for (i in 2:reps) {
      out$em <- rbind(out$em, res.em[i, , ])
      out$ep <- rbind(out$ep, res.ep[i, , ])
      out$np  <- rbind(out$np, res.np[i, , ])
      out$emalpha <- c(out$emalpha, res.emalpha[i, ])
      out$epalpha <- c(out$epalpha, res.epalpha[i, ])
      out$thetaalpha <- c(out$thetaalpha, res.thetaalpha[i, ])
      out$thetabeta <- c(out$thetabeta, res.thetabeta[i, ])
      out$embeta <- c(out$embeta, res.embeta[i, ])
      out$epbeta <- c(out$epbeta, res.epbeta[i, ])
    }
  
  endtime <- Sys.time() - start.time
  if (!quiet) 
    cat("\tAggregated error parameters\n")
  out$anames <- anames
  out$onames <- onames
  out$nactors <- n
  out$nrelations <- m
  out$reps <- reps
  out$draws <- draws
  out$burntime <- burntime
  out$tinning <- tinning
  out$alpha.acc.em <- res.alpha.acc.em
  out$alpha.acc.ep <- res.alpha.acc.ep
  out$beta.acc.em <- res.beta.acc.em
  out$beta.acc.ep <- res.beta.acc.ep
  out$total.alpha.em <- res.total.alpha.em
  out$alpha.acc.rate.em <- alpha.acc.rate.em
  out$alpha.acc.rate.ep <- alpha.acc.rate.ep
  out$beta.acc.rate.em <- beta.acc.rate.em
  out$beta.acc.rate.ep <- beta.acc.rate.ep
  out$alpha.acc.rate.theta <- alpha.acc.rate.theta
  out$beta.acc.rate.theta <- beta.acc.rate.theta
  out$runtime <- endtime
  out$model <- "node"
  class(out) <- c("bnsm.node", "bnsm")
  out
}



bnsm.one.theta <- function(dat, emp = c(3,9), epp = c(3,9)
                           , thetaalpha = 1, thetabeta = 9
                           , diag = FALSE, mode = "graph", model.checking = TRUE
                           , reps = 3, draws = 1500, tinning = 10, burntime = 500
                           , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE){
  start.time <- Sys.time()
  dat <- as.sociomatrix.sna(dat, simplify = TRUE)
  if (is.list(dat)) 
    stop("All bbnam input graphs must be of the same order.")
  if (length(dim(dat)) == 2) 
    dat <- array(dat, dim = c(1, NROW(dat), NCOL(dat)))
  if (reps == 1) 
    compute.sqrtrhat <- FALSE
  m <- dim(dat)[1]
  n <- dim(dat)[2]
  d <- dat
  slen <- burntime + draws*tinning
  n.draws <- burntime+tinning*(1:draws)
  
  out <- list()
  tieprob <- array(NA, dim=c(m,n,n))
  
  if ((!is.matrix(emp)) || (NROW(emp) != n) || (NCOL(emp) != 2)) {
    if (length(emp) == 2) 
      emprior <- sapply(emp, rep, n)
    else emprior <- matrix(emp, n, 2)
  }
  if ((!is.matrix(epp)) || (NROW(epp) != n) || (NCOL(epp) != 2)) {
    if (length(epp) == 2) 
      epprior <- sapply(epp, rep, n)
    else epprior <- matrix(epp, n, 2)
  }
  if (is.null(anames)) 
    anames <- paste("a", 1:n, sep = "")
  if (is.null(onames)) 
    onames <- paste("o", 1:m, sep = "")
  
  if (!diag) 
    d <- diag.remove(d)
  if (!quiet) 
    cat("Creating temporary variables and drawing initial conditions....\n")
  
  res.a <- array(dim = c(reps, draws, m, n, n))
  res.pred.y <- array(dim = c(reps, draws, m, n, n))
  res.em <- array(dim = c(reps, draws, n))
  res.ep <- array(dim = c(reps, draws, n))
  res.np <- array(dim = c(reps, draws, m))
  emalpha <- emp[1]
  embeta <- emp[2]
  epalpha <- epp[1]
  epbeta <- epp[2]
  
  for(j in 1:reps){
    em <- rbeta(n, emalpha, embeta)
    ep <- rbeta(n, epalpha, epbeta)
    np <- rbeta(m, thetaalpha, thetabeta)
    
    ### updated values setting ###### 
    for (i in 2:slen){
      if(model.checking == TRUE){
        sample.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
        new.a <- sample.a$a
        pred.y <- sample.a$pred_y
      }else{
        new.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
      }
      
      cem <- matrix(nrow = n, ncol = 2)
      cep <- matrix(nrow = n, ncol = 2)
      
      cem[, 1] <- rowSums(colSums((1 - d) * new.a,  na.rm = TRUE))
      cem[, 2] <- rowSums(colSums(d * new.a,  na.rm = TRUE))
      cep[, 1] <- rowSums(colSums(d * (1-new.a),  na.rm = TRUE))
      cep[, 2] <- rowSums(colSums((1 - d) * (1-new.a),  na.rm = TRUE))
      
      new.em <- rbeta(n, emalpha + cem[, 1], embeta + cem[, 2])
      new.ep <- rbeta(n, epalpha + cep[, 1], epbeta + cep[, 2])
      
      new.np <- rbeta(m, rowSums(new.a, na.rm = TRUE)+thetaalpha, rowSums(1-new.a, na.rm = TRUE)+thetabeta)
      #### storing values after burnning and every tinning
      if(i %in% n.draws){
        k <- (i-burntime)/tinning
        res.a[j,k,,,] <- new.a
        res.pred.y[j,k,,,] <- pred.y
        res.em[j,k,] <- new.em
        res.ep[j,k,] <- new.ep
        res.np[j,k,] <- new.np
      }
      rm(em, ep, new.a)
      em <- new.em
      ep <- new.ep
      np <- new.np
      
      rm(new.em, new.ep, new.np)
    }
  }
  
  out$net <- res.a 
  out$pred.y <- res.pred.y
  out$em <- res.em[1, , ]
  out$ep <- res.ep[1, , ]
  out$np  <- res.np[1, ,]
  
  if (reps >= 2){ 
    for (i in 2:reps) {
      out$em <- rbind(out$em, res.em[i, , ])
      out$ep <- rbind(out$ep, res.ep[i, , ])
      out$np <- rbind(out$np, res.np[i, , ])
    }
  }
  
  endtime <- Sys.time() - start.time
  #out$nprior <- nprior
  out$anames <- anames
  out$onames <- onames
  out$nactors <- n
  out$nrelations <- m
  out$reps <- reps
  out$draws <- draws
  out$burntime <- burntime
  out$tinning <- tinning
  out$runtime <- endtime
  out$model <- "node"
  class(out) <- c("bnsm.node", "bnsm")
  out
}




bnsm.hyper <- function(dat, np = netmean(dat), emsd = list(ema = 10, emb = 10, alpha.sd = 0.25, beta.sd = 0.25)
                       , epsd = list(epa = 3, epb = 9, alpha.sd = 0.25, beta.sd = 0.25)
                       , diag = FALSE, mode = "graph", model.checking = TRUE
                       , reps = 3, draws = 1500, tinning = 10, burntime = 500
                       , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE) {
  
  start.time <- Sys.time()
  dat <- as.sociomatrix.sna(dat, simplify = TRUE)
  
  if (is.list(dat)) 
    stop("All bbnam input graphs must be of the same order.")
  if (length(dim(dat)) == 2) 
    dat <- array(dat, dim = c(1, NROW(dat), NCOL(dat)))
  if (reps == 1) 
    compute.sqrtrhat <- FALSE
  m <- dim(dat)[1]
  n <- dim(dat)[2]
  d <- dat
  slen <- burntime + draws*tinning
  n.draws <- burntime+tinning*(1:draws)
  
  out <- list()
  tieprob <- array(NA, dim=c(m,n,n))
  if((NROW(np) == m)||(length(np) == 1)) {
    nprior <- array(np, dim = c(m, n, n))
  }else if(NROW(np) == n) {
    nprior <- aperm(array(np, dim = c(n, m, n)), c(2,1,3))
  }else{
    stop("nprior is not correct.")
  }
  
  if (is.null(anames)) 
    anames <- paste("a", 1:n, sep = "")
  if (is.null(onames)) 
    onames <- paste("o", 1:m, sep = "")
  
  if (!diag) 
    d <- diag.remove(d)
  if (!quiet) 
    cat("Creating temporary variables and drawing initial conditions....\n")
  
  res.a <- array(dim = c(reps, draws, m, n, n))
  res.pred.y <- array(dim = c(reps, draws, m, n, n))
  res.em <- array(dim = c(reps, draws, n))
  res.ep <- array(dim = c(reps, draws, n))
  res.emalpha <- array(dim = c(reps,draws))
  res.embeta <- array(dim = c(reps,draws))
  res.epalpha <- array(dim = c(reps,draws))
  res.epbeta <- array(dim = c(reps,draws))
  res.alpha.acc.em <- rep(NA, reps)
  res.alpha.acc.ep <- rep(NA, reps)
  res.beta.acc.em <- rep(NA, reps)
  res.beta.acc.ep <- rep(NA, reps)
  res.total.alpha.em <- rep(NA, reps)
  alpha.acc.rate.em <- rep(NA, reps)
  alpha.acc.rate.ep <- rep(NA, reps)
  beta.acc.rate.em <- rep(NA, reps)
  beta.acc.rate.ep <- rep(NA, reps)
  
  for(j in 1:reps){
    #### initial values setting #####
    #alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
    #alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
    
    emalpha <- emsd$ema
    embeta <- emsd$emb
    epalpha <- epsd$epa
    epbeta <- epsd$epb
    
    em <- rbeta(n, emalpha, embeta)
    ep <- rbeta(n, epalpha, epbeta)
    
    ### updated values setting ######
    for(i in 2:slen){
      if(model.checking == TRUE){
        sample.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
        new.a <- sample.a$a
        pred.y <- sample.a$pred_y
      }else{
        new.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
      }
      
      cem <- matrix(nrow = n, ncol = 2)
      cep <- matrix(nrow = n, ncol = 2)
      
      cem[, 1] <- rowSums(colSums((1 - d) * new.a,  na.rm = TRUE))
      cem[, 2] <- rowSums(colSums(d * new.a,  na.rm = TRUE))
      cep[, 1] <- rowSums(colSums(d * (1-new.a),  na.rm = TRUE))
      cep[, 2] <- rowSums(colSums((1 - d) * (1-new.a),  na.rm = TRUE))
      
      new.emalpha <- draw.alpha(emalpha, embeta, theta = em, prop.sd = emsd$alpha.sd, type = "em")
      new.embeta <- draw.beta(new.emalpha, embeta, theta = em, prop.sd = emsd$beta.sd, type = "em")
      
      new.epalpha <- draw.alpha(epalpha, epbeta, theta = ep, prop.sd = epsd$alpha.sd, type = "ep")
      new.epbeta <- draw.beta(new.epalpha, epbeta, theta = ep, prop.sd = epsd$beta.sd, type = "ep")
      
      new.em <- rbeta(n, new.emalpha + cem[, 1], new.embeta + cem[, 2])
      new.ep <- rbeta(n, new.epalpha + cep[, 1], new.epbeta + cep[, 2])
      
      #### storing values after burnning and every tinning
      if(i %in% n.draws){
        k <- (i-burntime)/tinning
        res.a[j,k,,,] <- new.a
        res.pred.y[j,k,,,] <- pred.y
        res.emalpha[j,k] <- new.emalpha
        res.embeta[j,k] <- new.embeta
        res.epalpha[j,k] <- new.epalpha
        res.epbeta[j,k] <- new.epbeta
        
        res.em[j,k,] <- new.em
        res.ep[j,k,] <- new.ep
      }
      
      rm(emalpha, embeta, epalpha, epbeta, em, ep, new.a)
      emalpha <- new.emalpha
      embeta <- new.embeta
      epalpha <- new.epalpha
      epbeta <- new.epbeta
      
      em <- new.em
      ep <- new.ep
      
      rm(new.emalpha, new.embeta, new.epalpha, new.epbeta, new.em, new.ep)
    }
    print(alpha.acc.em)
    print(total.alpha.em)
    print(alpha.acc.em/total.alpha.em)
    res.alpha.acc.em[j] <- alpha.acc.em
    res.alpha.acc.ep[j] <- alpha.acc.ep
    res.beta.acc.em[j] <- beta.acc.em
    res.beta.acc.ep[j] <- beta.acc.ep
    res.total.alpha.em[j] <- total.alpha.em
  }
  
  alpha.acc.rate.em[1] <- res.alpha.acc.em[1]/total.alpha.em[1]
  alpha.acc.rate.ep[1] <- res.alpha.acc.ep[1]/total.alpha.em[1]
  beta.acc.rate.em[1] <- res.beta.acc.em[1]/total.alpha.em[1]
  beta.acc.rate.ep[1] <- res.beta.acc.ep[1]/total.alpha.em[1]
  for(j in 2:reps){
    alpha.acc.rate.em[j] <- (res.alpha.acc.em[j]-res.alpha.acc.em[(j-1)])/(total.alpha.em[j]-total.alpha.em[(j-1)])
    alpha.acc.rate.ep[j] <- (res.alpha.acc.ep[j]-res.alpha.acc.ep[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    beta.acc.rate.em[j] <- (res.beta.acc.em[j]-res.beta.acc.em[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    beta.acc.rate.ep[j] <- (res.beta.acc.ep[j]-res.beta.acc.ep[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
  }
  out$net <- res.a 
  out$pred.y <- res.pred.y
  out$em <- res.em[1, , ]
  out$ep <- res.ep[1, , ]
  out$emalpha <- res.emalpha[1, ]
  out$embeta <- res.embeta[1, ]
  out$epalpha <- res.epalpha[1, ]
  out$epbeta  <- res.epbeta[1, ]
  
  if (reps >= 2) 
    for (i in 2:reps) {
      out$em <- rbind(out$em, res.em[i, , ])
      out$ep <- rbind(out$ep, res.ep[i, , ])
      out$emalpha <- c(out$emalpha, res.emalpha[i, ])
      out$epalpha <- c(out$epalpha, res.epalpha[i, ])
      out$embeta <- c(out$embeta, res.embeta[i, ])
      out$epbeta <- c(out$epbeta, res.epbeta[i, ])
    }
  
  endtime <- Sys.time() - start.time
  if (!quiet) 
    cat("\tAggregated error parameters\n")
  out$nprior <- nprior
  out$anames <- anames
  out$onames <- onames
  out$nactors <- n
  out$nrelations <- m
  out$reps <- reps
  out$draws <- draws
  out$burntime <- burntime
  out$tinning <- tinning
  out$alpha.acc.em <- res.alpha.acc.em
  out$alpha.acc.ep <- res.alpha.acc.ep
  out$beta.acc.em <- res.beta.acc.em
  out$beta.acc.ep <- res.beta.acc.ep
  out$total.alpha.em <- res.total.alpha.em
  out$alpha.acc.rate.em <- alpha.acc.rate.em
  out$alpha.acc.rate.ep <- alpha.acc.rate.ep
  out$beta.acc.rate.em <- beta.acc.rate.em
  out$beta.acc.rate.ep <- beta.acc.rate.ep
  out$runtime <- endtime
  out$model <- "node"
  class(out) <- c("bnsm.node", "bnsm")
  out
}



bnsm.hyper.theta <- function(dat, emsd = list(ema = 10, emb = 10, alpha.sd = 0.25, beta.sd = 0.25)
                             , epsd = list(epa = 3, epb = 9, alpha.sd = 0.25, beta.sd = 0.25)
                             , thetaalpha = 1, thetabeta = 9, model.checking = TRUE
                             , diag = FALSE, mode = "graph"
                             , reps = 3, draws = 1500, tinning = 10, burntime = 500
                             , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE) {
  
  start.time <- Sys.time()
  dat <- as.sociomatrix.sna(dat, simplify = TRUE)
  
  if (is.list(dat)) 
    stop("All bbnam input graphs must be of the same order.")
  if (length(dim(dat)) == 2) 
    dat <- array(dat, dim = c(1, NROW(dat), NCOL(dat)))
  if (reps == 1) 
    compute.sqrtrhat <- FALSE
  m <- dim(dat)[1]
  n <- dim(dat)[2]
  d <- dat
  slen <- burntime + draws*tinning
  n.draws <- burntime+tinning*(1:draws)
  
  out <- list()
  tieprob <- array(NA, dim=c(m,n,n))
  
  if (is.null(anames)) 
    anames <- paste("a", 1:n, sep = "")
  if (is.null(onames)) 
    onames <- paste("o", 1:m, sep = "")
  
  if (!diag) 
    d <- diag.remove(d)
  if (!quiet) 
    cat("Creating temporary variables and drawing initial conditions....\n")
  
  
  
  res.a <- array(dim = c(reps, draws, m, n, n))
  res.pred.y <- array(dim = c(reps, draws, m, n, n))
  res.em <- array(dim = c(reps, draws, n))
  res.ep <- array(dim = c(reps, draws, n))
  res.np <- array(dim = c(reps, draws, m))
  #res.np <- array(dim = c(reps, draws, n))
  res.emalpha <- array(dim = c(reps,draws))
  res.embeta <- array(dim = c(reps,draws))
  res.epalpha <- array(dim = c(reps,draws))
  res.epbeta <- array(dim = c(reps,draws))
  res.alpha.acc.em <- rep(NA, reps)
  res.alpha.acc.ep <- rep(NA, reps)
  res.beta.acc.em <- rep(NA, reps)
  res.beta.acc.ep <- rep(NA, reps)
  res.total.alpha.em <- rep(NA, reps)
  alpha.acc.rate.em <- rep(NA, reps)
  alpha.acc.rate.ep <- rep(NA, reps)
  beta.acc.rate.em <- rep(NA, reps)
  beta.acc.rate.ep <- rep(NA, reps)
  
  for(j in 1:reps){
    #### initial values setting #####
    #alpha.acc.em = 0; beta.acc.em = 0; total.alpha.em = 0; total.beta.em = 0
    #alpha.acc.ep = 0; beta.acc.ep = 0; total.alpha.ep = 0; total.beta.ep = 0
    
    #new.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np)
    
    emalpha <- emsd$ema
    embeta <- emsd$emb
    epalpha <- epsd$epa
    epbeta <- epsd$epb
    
    em <- rbeta(n, emalpha, embeta)
    ep <- rbeta(n, epalpha, epbeta)
    
    np <- rbeta(m, thetaalpha, thetabeta)
    #np <- rbeta(n, thetaalpha, thetabeta)
    
    ### updated values setting ######
    for(i in 2:slen){
      if(model.checking == TRUE){
        sample.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
        new.a <- sample.a$a
        pred.y <- sample.a$pred_y
      }else{
        new.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np)
      }
      
      cem <- matrix(nrow = n, ncol = 2)
      cep <- matrix(nrow = n, ncol = 2)
      
      cem[, 1] <- rowSums(colSums((1 - d) * new.a,  na.rm = TRUE))
      cem[, 2] <- rowSums(colSums(d * new.a,  na.rm = TRUE))
      cep[, 1] <- rowSums(colSums(d * (1-new.a),  na.rm = TRUE))
      cep[, 2] <- rowSums(colSums((1 - d) * (1-new.a),  na.rm = TRUE))
      
      new.emalpha <- draw.alpha(emalpha, embeta, theta = em, prop.sd = emsd$alpha.sd, type = "em")
      new.embeta <- draw.beta(new.emalpha, embeta, theta = em, prop.sd = emsd$beta.sd, type = "em")
      
      new.epalpha <- draw.alpha(epalpha, epbeta, theta = ep, prop.sd = epsd$alpha.sd, type = "ep")
      new.epbeta <- draw.beta(new.epalpha, epbeta, theta = ep, prop.sd = epsd$beta.sd, type = "ep")
      
      new.em <- rbeta(n, new.emalpha + cem[, 1], new.embeta + cem[, 2])
      new.ep <- rbeta(n, new.epalpha + cep[, 1], new.epbeta + cep[, 2])
      
      new.np <- rbeta(m, rowSums(new.a, na.rm = TRUE)+thetaalpha, rowSums(1-new.a, na.rm = TRUE)+thetabeta)
      #new.np <- rbeta(n, rowSums(colSums(new.a,  na.rm = TRUE))+thetaalpha, rowSums(colSums(1-new.a,  na.rm = TRUE))+thetabeta)
      
      #### storing values after burnning and every tinning
      if(i %in% n.draws){
        k <- (i-burntime)/tinning
        res.a[j,k,,,] <- new.a
        res.pred.y[j,k,,,] <- pred.y
        res.emalpha[j,k] <- new.emalpha
        res.embeta[j,k] <- new.embeta
        res.epalpha[j,k] <- new.epalpha
        res.epbeta[j,k] <- new.epbeta
        
        res.em[j,k,] <- new.em
        res.ep[j,k,] <- new.ep
        
        res.np[j,k,] <- new.np
      }
      
      rm(emalpha, embeta, epalpha, epbeta, em, ep, new.a, pred.y)
      emalpha <- new.emalpha
      embeta <- new.embeta
      epalpha <- new.epalpha
      epbeta <- new.epbeta
      
      em <- new.em
      ep <- new.ep
      np <- new.np
      rm(new.emalpha, new.embeta, new.epalpha, new.epbeta, new.em, new.ep, new.np)
    }
    print(alpha.acc.em)
    print(total.alpha.em)
    print(alpha.acc.em/total.alpha.em)
    res.alpha.acc.em[j] <- alpha.acc.em
    res.alpha.acc.ep[j] <- alpha.acc.ep
    res.beta.acc.em[j] <- beta.acc.em
    res.beta.acc.ep[j] <- beta.acc.ep
    res.total.alpha.em[j] <- total.alpha.em
  }
  
  alpha.acc.rate.em[1] <- res.alpha.acc.em[1]/total.alpha.em[1]
  alpha.acc.rate.ep[1] <- res.alpha.acc.ep[1]/total.alpha.em[1]
  beta.acc.rate.em[1] <- res.beta.acc.em[1]/total.alpha.em[1]
  beta.acc.rate.ep[1] <- res.beta.acc.ep[1]/total.alpha.em[1]
  for(j in 2:reps){
    alpha.acc.rate.em[j] <- (res.alpha.acc.em[j]-res.alpha.acc.em[(j-1)])/(total.alpha.em[j]-total.alpha.em[(j-1)])
    alpha.acc.rate.ep[j] <- (res.alpha.acc.ep[j]-res.alpha.acc.ep[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    beta.acc.rate.em[j] <- (res.beta.acc.em[j]-res.beta.acc.em[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
    beta.acc.rate.ep[j] <- (res.beta.acc.ep[j]-res.beta.acc.ep[j-1])/(total.alpha.em[j]-total.alpha.em[j-1])
  }
  out$net <- res.a 
  out$pred.y <- res.pred.y
  out$em <- res.em[1, , ]
  out$ep <- res.ep[1, , ]
  out$np  <- res.np[1, ,]
  out$emalpha <- res.emalpha[1, ]
  out$embeta <- res.embeta[1, ]
  out$epalpha <- res.epalpha[1, ]
  out$epbeta  <- res.epbeta[1, ]
  
  
  if (reps >= 2) 
    for (i in 2:reps) {
      out$em <- rbind(out$em, res.em[i, , ])
      out$ep <- rbind(out$ep, res.ep[i, , ])
      out$np <- rbind(out$np, res.np[i, , ])
      out$emalpha <- c(out$emalpha, res.emalpha[i, ])
      out$epalpha <- c(out$epalpha, res.epalpha[i, ])
      out$embeta <- c(out$embeta, res.embeta[i, ])
      out$epbeta <- c(out$epbeta, res.epbeta[i, ])
    }
  
  endtime <- Sys.time() - start.time
  if (!quiet) 
    cat("\tAggregated error parameters\n")
  out$anames <- anames
  out$onames <- onames
  out$nactors <- n
  out$nrelations <- m
  out$reps <- reps
  out$draws <- draws
  out$burntime <- burntime
  out$tinning <- tinning
  out$alpha.acc.em <- res.alpha.acc.em
  out$alpha.acc.ep <- res.alpha.acc.ep
  out$beta.acc.em <- res.beta.acc.em
  out$beta.acc.ep <- res.beta.acc.ep
  out$total.alpha.em <- res.total.alpha.em
  out$alpha.acc.rate.em <- alpha.acc.rate.em
  out$alpha.acc.rate.ep <- alpha.acc.rate.ep
  out$beta.acc.rate.em <- beta.acc.rate.em
  out$beta.acc.rate.ep <- beta.acc.rate.ep
  out$runtime <- endtime
  out$model <- "node"
  class(out) <- c("bnsm.node", "bnsm")
  out
}



bnsm.one <- function(dat, np = netmean(dat), emp = c(3,9), epp = c(3,9)
                     , diag = FALSE, mode = "graph", model.checking = FALSE
                     , reps = 3, draws = 1500, tinning = 10, burntime = 500
                     , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE){
  start.time <- Sys.time()
  dat <- as.sociomatrix.sna(dat, simplify = TRUE)
  if (is.list(dat)) 
    stop("All bbnam input graphs must be of the same order.")
  if (length(dim(dat)) == 2) 
    dat <- array(dat, dim = c(1, NROW(dat), NCOL(dat)))
  if (reps == 1) 
    compute.sqrtrhat <- FALSE
  m <- dim(dat)[1]
  n <- dim(dat)[2]
  d <- dat
  slen <- burntime + draws*tinning
  n.draws <- burntime+tinning*(1:draws)
  
  
  out <- list()
  tieprob <- array(NA, dim=c(m,n,n))
  if((NROW(np) == m)||(length(np) == 1)) {
    nprior <- array(np, dim = c(m, n, n))
  }else if(NROW(np) == n) {
    nprior <- aperm(array(np, dim = c(n, m, n)), c(2,1,3))
  }else{
    stop("nprior is not correct.")
  }
  
  if ((!is.matrix(emp)) || (NROW(emp) != n) || (NCOL(emp) != 2)) {
    if (length(emp) == 2) 
      emprior <- sapply(emp, rep, n)
    else emprior <- matrix(emp, n, 2)
  }
  if ((!is.matrix(epp)) || (NROW(epp) != n) || (NCOL(epp) != 2)) {
    if (length(epp) == 2) 
      epprior <- sapply(epp, rep, n)
    else epprior <- matrix(epp, n, 2)
  }
  if (is.null(anames)) 
    anames <- paste("a", 1:n, sep = "")
  if (is.null(onames)) 
    onames <- paste("o", 1:m, sep = "")
  
  if (!diag) 
    d <- diag.remove(d)
  if (!quiet) 
    cat("Creating temporary variables and drawing initial conditions....\n")
  
  res.a <- array(dim = c(reps, draws, m, n, n))
  res.pred.y <- array(dim = c(reps, draws, m, n, n))
  res.em <- array(dim = c(reps, draws, n))
  res.ep <- array(dim = c(reps, draws, n))
  emalpha <- emp[1]
  embeta <- emp[2]
  epalpha <- epp[1]
  epbeta <- epp[2]
  
  for(j in 1:reps){
    em <- rbeta(n, emalpha, embeta)
    ep <- rbeta(n, epalpha, epbeta)

  ### updated values setting ######
    for (i in 2:slen){
      if(model.checking == TRUE){
        sample.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
        new.a <- sample.a$a
        pred.y <- sample.a$pred_y
      }else{
        new.a <- sampletie_func(ep = ep, em = em, m = m, n = n, d = d, mode = mode, diag = diag, np = np, model.checking = model.checking)
      }
      
      cem <- matrix(nrow = n, ncol = 2)
      cep <- matrix(nrow = n, ncol = 2)
      
      cem[, 1] <- rowSums(colSums((1 - d) * new.a,  na.rm = TRUE))
      cem[, 2] <- rowSums(colSums(d * new.a,  na.rm = TRUE))
      cep[, 1] <- rowSums(colSums(d * (1-new.a),  na.rm = TRUE))
      cep[, 2] <- rowSums(colSums((1 - d) * (1-new.a),  na.rm = TRUE))
      
      new.em <- rbeta(n, emalpha + cem[, 1], embeta + cem[, 2])
      new.ep <- rbeta(n, epalpha + cep[, 1], epbeta + cep[, 2])
      
      #### storing values after burnning and every tinning
      if(i %in% n.draws){
        k <- (i-burntime)/tinning
        res.a[j,k,,,] <- new.a
        res.pred.y[j,k,,,] <- pred.y
        res.em[j,k,] <- new.em
        res.ep[j,k,] <- new.ep
      }
      rm(em, ep, new.a)
      em <- new.em
      ep <- new.ep
      
      rm(new.em, new.ep)
    }
  }
  
  out$net <- res.a 
  out$pred.y <- res.pred.y
  out$em <- res.em[1, , ]
  out$ep <- res.ep[1, , ]
  
  if (reps >= 2) 
    for (i in 2:reps) {
      out$em <- rbind(out$em, res.em[i, , ])
      out$ep <- rbind(out$ep, res.ep[i, , ])
    }
  
  endtime <- Sys.time() - start.time
  out$nprior <- nprior
  out$anames <- anames
  out$onames <- onames
  out$nactors <- n
  out$nrelations <- m
  out$reps <- reps
  out$draws <- draws
  out$burntime <- burntime
  out$tinning <- tinning
  out$runtime <- endtime
  out$model <- "node"
  class(out) <- c("bnsm.node", "bnsm")
  out
}



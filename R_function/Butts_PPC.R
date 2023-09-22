
## load real data  (we cannot share the data)
load("net.symmetric.application.RData") 

## load function 
source("thetaprior_function_forpaper.R")
library(sna)

n_draws = 1500
n_tinning = 100
n_burntime = 15000

schools <- c("ARTS02", "ARTS03", "PREP02", "TECH03")

## Fitting real data to the Butts method with six different relations shared across the network 
## We have four datasets named ARTS02, ARTS03, PREP02, TECH03.

for(school in schools){
  Butts_all <- list()
  Butts_all <- bnsm.one.theta(dat = get(school)
                              , emp = c(1,11), epp = c(1,11)
                              , thetaalpha = 0.5, thetabeta = 0.5
                              , diag = FALSE, mode = "graph", model.checking = TRUE
                              , reps = 3, draws = n_draws, tinning = n_draws, burntime = n_burntime
                              , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE)
  
  print("Butts all done")
  save.image(file = paste("Butts_PPC_", school, ".RData", sep = ""))
  
  
  ## Fitting real data to the Butts method with not shared across the network 
  m1.res <- m2.res <- m3.res <- m4.res <- m5.res <- m6.res <- list()
  
  for(k in 1:6){
    print(k)
    assign(paste("m",k,sep=""), bnsm.one.theta(dat = get(school)[k,,]
                                               , emp = c(1,11), epp = c(1,11)
                                               , thetaalpha = 0.5, thetabeta = 0.5
                                               , diag = FALSE, mode = "graph", model.checking = TRUE
                                               , reps = 3, draws = n_draws, tinning = n_draws, burntime = n_burntime
                                               , quiet = TRUE, anames = NULL, onames = NULL, compute.sqrtrhat = TRUE))
    save.image(file = paste("Butts_PPC_one_", school, ".RData", sep=""))
  }
}

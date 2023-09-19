######################################################################
### calculate accuracy
## threshold : when \hat{theta}>threshold -> 1 , otherwise -> 0
## bnet : posterior distribution from result (b$net)
## g : original network data
## forall = TRUE : calculate accuracy for all network together.
## forall = FALSE : calculate accuracy for each network separately.
######################################################################


calaccuracy <- function(threshold, bnet, g, forall){
  g1 <- diag.remove(g)
  g1 <- lower.tri.remove(g1)
  b01 <- (apply(bnet, c(2,3,4), mean)>=threshold)*1
  b01re <- diag.remove(b01)
  b01re <- lower.tri.remove(b01re)
  
  if(forall == TRUE){
    res <- matrix(NA, nrow = 9, ncol = 1)
    res[1,] <- mean(g1 == b01re, na.rm = TRUE)
    res[2,] <- sum(g1 == b01re, na.rm = TRUE)
    res[3,] <- sum(g1 != b01re, na.rm = TRUE)
    res[4,] <- sum(g1==0 & b01re==0, na.rm = TRUE)
    res[5,] <- sum(g1==0 & b01re==1, na.rm = TRUE)
    res[6,] <- sum(g1==1 & b01re==0, na.rm = TRUE)
    res[7,] <- sum(g1==1 & b01re==1, na.rm = TRUE)
  }else{
    res <- matrix(NA, nrow = 9, ncol = dim(g1)[1])
    res[1,] <- apply(g1 == b01re, 1, mean, na.rm = TRUE)
    res[2,] <- apply(g1 == b01re, 1, sum, na.rm = TRUE)
    res[3,] <- apply(g1 != b01re, 1, sum, na.rm = TRUE)
    res[4,] <- apply(g1==0 & b01re==0, 1, sum, na.rm = TRUE)
    res[5,] <- apply(g1==0 & b01re==1, 1, sum, na.rm = TRUE)
    res[6,] <- apply(g1==1 & b01re==0, 1, sum, na.rm = TRUE)
    res[7,] <- apply(g1==1 & b01re==1, 1, sum, na.rm = TRUE)
  }
  res[8,] <- res[4,]/(res[4,]+res[5,])
  res[9,] <- res[7,]/(res[6,]+res[7,])
  rownames(res) <- c("all", "equ", "neq", "00", "01", "10", "11", "00/true0", "11/true1")
  return(res) 
}



######################################################################
### compare original network and observed network
## g : original network data
## obs.g : observed network data
######################################################################

compare_ori_obs <- function(g, obs.g, forall){
  g1 <- diag.remove(g)
  obs.g1 <- diag.remove(obs.g)
  if(forall == TRUE){
    res <- matrix(NA, nrow = 9, ncol = 1)
    res[1,] <- mean(g1 == obs.g1, na.rm = TRUE)
    res[2,] <- sum(g1 == obs.g1, na.rm = TRUE)
    res[3,] <- sum(g1 != obs.g1, na.rm = TRUE)
    res[4,] <- sum(g1==0 & obs.g1==0, na.rm = TRUE)
    res[5,] <- sum(g1==0 & obs.g1==1, na.rm = TRUE)
    res[6,] <- sum(g1==1 & obs.g1==0, na.rm = TRUE)
    res[7,] <- sum(g1==1 & obs.g1==1, na.rm = TRUE)
  }else{
    res <- matrix(NA, nrow = 9, ncol = dim(g)[1])
    res[1,] <- apply(g1 == obs.g1, 1, mean, na.rm = TRUE)
    res[2,] <- apply(g1 == obs.g1, 1, sum, na.rm = TRUE)
    res[3,] <- apply(g1 != obs.g1, 1, sum, na.rm = TRUE)
    res[4,] <- apply(g1==0 & obs.g1==0, 1, sum, na.rm = TRUE)
    res[5,] <- apply(g1==0 & obs.g1==1, 1, sum, na.rm = TRUE)
    res[6,] <- apply(g1==1 & obs.g1==0, 1, sum, na.rm = TRUE)
    res[7,] <- apply(g1==1 & obs.g1==1, 1, sum, na.rm = TRUE)
  }
  res[8,] <- res[4,]/(res[4,]+res[5,])
  res[9,] <- res[7,]/(res[6,]+res[7,])
  rownames(res) <- c("all", "equ", "neq", "00", "01", "10", "11", "00/true0", "11/true1")
  return(res) 
}


union.intersection.net <- function(observed){
  union.net <- intersection.net <- observed
  for(i in 1:dim(observed)[1]){
    union.net[i,,] <- ifelse(observed[i,,]+t(observed[i,,]) > 0, 1, 0)
    intersection.net[i,,] <- ifelse(observed[i,,]+t(observed[i,,]) > 1, 1, 0)
  }
  return(list(union.net = union.net, intersection.net = intersection.net))
}


######################################################################
### Sociomatrix plot with cell probabilities
### Edited code that already exists.
######################################################################


plot.sociomatrix <- function (x, labels = NULL, drawlab = TRUE, diaglab = TRUE, drawlines = TRUE, 
                              xlab = NULL, ylab = NULL, cex.lab = 1, font.lab = 1, col.lab = 1, cellprob = FALSE,
                              cellprob.cex = 1, cellprob.col = 2,
                              scale.values = TRUE, cell.col = gray, na.cell.col = rainbow, ...) 
{
  if ((!inherits(x, c("matrix", "array", "data.frame"))) || 
      (length(dim(x)) > 2)) 
    x <- as.sociomatrix.sna(x)
  if (is.list(x)) 
    x <- x[[1]]
  n <- dim(x)[1]
  o <- dim(x)[2]
  if (is.null(labels)) 
    labels <- list(NULL, NULL)
  if (is.null(labels[[1]])) {
    if (is.null(rownames(x))) 
      labels[[1]] <- 1:dim(x)[1]
    else labels[[1]] <- rownames(x)
  }
  if (is.null(labels[[2]])) {
    if (is.null(colnames(x))) 
      labels[[2]] <- 1:dim(x)[2]
    else labels[[2]] <- colnames(x)
  }
  if (scale.values) 
    d <- 1 - x
  else d <- x
  if (is.null(xlab)) 
    xlab <- ""
  if (is.null(ylab)) 
    ylab <- ""
  plot(1, 1, xlim = c(0, o + 1), ylim = c(n + 1, 0), type = "n", 
       axes = FALSE, xlab = xlab, ylab = ylab, ...)
  
  for (i in 1:n) for (j in 1:o){
    if(is.na(d[i,j])){
      d[i,j] <- 1
      rect(j - 0.5, i + 0.5, j + 
             0.5, i - 0.5, col = na.cell.col(d[i, j]), xpd = TRUE, border = drawlines)
    }else{
      rect(j - 0.5, i + 0.5, j + 0.5, i - 0.5, col = cell.col(d[i, j]), xpd = TRUE, border = drawlines)
    }
  }
  
  rect(0.5, 0.5, o + 0.5, n + 0.5, col = NA, xpd = TRUE)
  if (drawlab) {
    text(rep(0, n), 1:n, labels[[1]], cex = cex.lab, font = font.lab, 
         col = col.lab)
    text(1:o, rep(0, o), labels[[2]], cex = cex.lab, font = font.lab, 
         col = col.lab)
  }
  if ((n == o) & (drawlab) & (diaglab)) 
    if (all(labels[[1]] == labels[[2]])) 
      text(1:o, 1:n, labels[[1]], cex = cex.lab, font = font.lab, 
           col = col.lab)
  if(cellprob){
    for (i in 2:n) for (j in 1:(i-1)) text(i,j, round(x[i,j],2), cex = cellprob.cex, font = font.lab, 
                                           col = cellprob.col)
  }
}


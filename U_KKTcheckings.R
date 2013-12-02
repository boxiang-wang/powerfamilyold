hsvm <- function(v, varlist) {
  delta = varlist$delta
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- (1 - r - delta/2) else dl <- (r - 1)^2/delta/2
  dl
}

dhsvm <- function(v, varlist) {
  delta = varlist$delta
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- -1 else dl <- (r - 1) / delta
  dl
}

power <- function(v, varlist) {
  qv = varlist$qv
  decib = qv / (qv + 1)
  r <- v[1]
  dl = ifelse(r > decib, r ^ (-qv) * (qv ^ qv) / ((qv + 1) ^ (qv + 1)), 1 - r)
  dl
}

dpower <- function(v, varlist) {
  qv = varlist$qv
  decib = qv / (qv + 1)
  r <- v[1]
  dl = ifelse(r > decib, (-1) * r ^ (-qv - 1) * decib ^ (qv + 1), -1)
}


margin <- function(b0, beta, y, x, loss = c("hsvm", "power"), delta=2, qv=2) {
  loss <- match.arg(loss)
  nobs <- nrow(x)
  b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
  link <- x %*% beta + b0MAT
  if (loss %in% c("hsvm", "power")) {
    r <- y * link
  } else r <- y - link
  fun <- paste("d", loss, sep = "")
  varlist = list(delta=delta, qv=qv)
  dMat <- apply(r, c(1, 2), eval(fun), varlist) #dMat1
  if (loss %in% c("hsvm", "power")) {
    yxdMat <- t(x) %*% (dMat * y)/nobs
  } else yxdMat <- t(x) %*% dMat/nobs
  yxdMat
}
# l is the number of lambda
# p is nvars, n is nobs
# dim(dMat) = n by l
# dim(yxdMat) = p by l

KKT = function(b0, beta, y, x, lambda, lambda2, thr, 
                loss = c("hsvm", "power"), delta=2, qv=2) {
  loss = match.arg(loss)
  dl = margin(b0, beta, y, x, loss=loss, delta=delta, qv=qv)
  ctr = 0
  #ccounts = 0
  for (l in 1:length(lambda)) {
    p = nrow(beta)
    for(j in 1:p)
    {
      if(beta[j,l]==0)
      {
        BB = abs(dl[j,l]) - lambda[l]
        #ccounts = ccounts + 1
        if (BB > thr) 
        {
          cat("violate at b = 0", BB, "\n")
          ctr <- ctr + 1
        }
      } else{
        AA = dl[j,l] + lambda[l] * sign(beta[j,l]) + lambda2 * beta[j,l]
        #ccounts = ccounts + 1
        if (abs(sum(AA)) >= thr)
        {
          cat("violate at b != 0", abs(sum(AA)), "\n")
          ctr <- ctr + 1
        }
        
      }
    }
  }
  ccounts = length(lambda) * 
    p
  cat("# of violations", ctr/length(lambda), "\n")
  cat("% of violations", ctr/ccounts*100, "%", "\n")
  return(ctr/length(lambda))
} 

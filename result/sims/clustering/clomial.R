library(Clomial)
set.seed(1)

args = commandArgs(TRUE)
dirname = args[1]
clusters = as.integer(args[2])

Dc <- read.csv(paste0(dirname, "/total.tsv"), sep="\t", row.names = 1)
Dt <- read.csv(paste0(dirname, "/alt.tsv"), sep="\t", row.names = 1)

f = function (Dc, Dt, Mu, P) 
{
  result <- list()
  log.likelihood <- 0
  llS <- c()
  for (i in 1:nrow(Mu)) {
    p <- (Mu[i, ]/2) %*% P
    p[p > 1] <- 1
    #print(i)
    #print(Dt[i,])
    #print(Dc[i,])
    #print(p)
    lli <- dbinom(as.numeric(Dt[i, ]), as.numeric(Dc[i, ]), as.numeric(p), log = TRUE)
    log.likelihood <- log.likelihood + sum(lli)
    llS <- rbind(llS, lli)
  }
  rownames(llS) <- rownames(Mu)
  result[["ll"]] <- log.likelihood
  result[["llS"]] <- llS
  return(result)
}

bic = function (Dc, Dt, Mu, P) 
{
  result <- list()
  ll <- f(Dc = Dc, Dt = Dt, Mu = Mu, P = P)$ll
  result[["ll"]] <- ll
  N <- nrow(Mu)
  C <- ncol(Mu)
  M <- ncol(P)
  obsNum <- sum(Dc)
  obsNumAverage <- obsNum/(N * M)
  aic <- -2 * ll + (N * C + M * (C - 1)) * 2
  bic <- -2 * ll + (N * C + M * (C - 1)) * log(obsNum)
  result[["bic"]] <- bic
  result[["aic"]] <- aic
  result[["obsNum"]] <- obsNum
  return(result)
}

print(dirname)
bics = c()
for (i in 3:clusters) {
  clomial = Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=i,doParal=FALSE,binomTryNum=25)
  chosen <- choose.best(models=clomial$models, doTalk=TRUE)$bestModel
  bics[i] <- bic(Dc=Dc,Dt=Dt, Mu=chosen$Mu, P=chosen$P)$bic
  print(paste("bic", i, "clusters:", bics[i]))
  write.table(bics, paste0(dirname, "/bics.tsv"), sep="\t")
  write.table(round(chosen$Mu), paste0(dirname, "/", i, ".tsv"), sep="\t")
}

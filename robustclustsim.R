library(mvtnorm)
library(clusterGeneration)
library(mice)
library(MixtureMissing)
library(mclust)
library(LearnClust)

dims <- c(4, 7, 10)
nclust <- c(3, 4, 5)
size <- c(100, 500, 1000)
noiseprop <- c(0.05, 0.10, 0.15, 0.20)
missingprop <- c(0.1, 0.2, 0.3, 0.4, 0.5)
noisethresh <- c(15, 20, 25)
means3 <- c(-6, 0, 6)
means4 <- c(-9, -3, 3, 9)
means5 <- c(-12, -6, 0, 6, 12)
df3 <- c(3, 25, 50)
df4 <- c(3, 25, 50, 75)
df5 <- c(3, 25, 50, 75, 100)

# Empty data frame to store simulation results
resdf <- data.frame(seed=numeric(),
                    nClust=numeric(),
                    nrows=numeric(),
                    ncols=numeric(),
                    noisep=numeric(),
                    missingp=numeric(),
                    pidist=numeric(),
                    mudist=numeric(),
                    vardist=numeric(),
                    method=character(),
                    ARIFull=numeric(),
                    ARINoNoise=numeric(),
                    stringsAsFactors=FALSE)
restdf <- data.frame(seed=numeric(),
                     nClust=numeric(),
                     nrows=numeric(),
                     ncols=numeric(),
                     noisep=numeric(),
                     missingp=numeric(),
                     pidist=numeric(),
                     mudist=numeric(),
                     vardist=numeric(),
                     dfdist = numeric(),
                     method=character(),
                     ARIFull=numeric(),
                     ARINoNoise=numeric(),
                     stringsAsFactors=FALSE)

num_iters <- 5400
run <- 1
count <- 1

for (clusters in nclust){
  j <- 1
  for (dim in dims){
    j <- 1
    for (nrow in size){
      j <- 1
      for (noise in noiseprop){
        j <- 1
        for (missing in missingprop){
          j <- 1
          count <- 1
          while (count <= 10){
            set.seed(j)
            # create mixing proportions vector
            ifelse(nrow == 1000,
                   pivec <- c(round(runif(clusters-1, min=1/clusters-0.05, max=1/clusters+0.05), 3)),
                   pivec <- c(round(runif(clusters-1, min=1/clusters-0.05, max=1/clusters+0.05), 2)))
            pivec <- c(pivec, 1-sum(pivec))
            cat('Pi has been constructed \n')
            if (length(pivec)==length(means3)){
              means <- means3
              noisethr <- noisethresh[1]
              dfs <- df3
            }
            if (length(pivec)==length(means4)){
              means <- means4
              noisethr <- noisethresh[2]
              dfs <- df4
            }
            if (length(pivec)==length(means5)){
              means <- means5
              noisethr <- noisethresh[3]
              dfs <- df5
            }
            # generate mixture components
            varsvec <- c()
            for (i in 1:length(pivec)){
              sigma1 <- genPositiveDefMat(dim=dim, covMethod = "unifcorrmat", rangeVar = c(1, 5))$Sigma
              varsvec <- c(varsvec, sigma1)
              mix1 <- rmvnorm(n=round(nrow*pivec[i]), mean=rep(means[i], dim), sigma=sigma1)
              mix1t <- rmvt(n=round(nrow*pivec[i]), delta=rep(means[i], dim), sigma=sigma1, df=dfs[i], type='shifted')
              mix1 <- as.data.frame(mix1)
              mix1t <- as.data.frame(mix1t)
              mix1$label <- i
              mix1t$label <- i
              ifelse(i==1, mixdf <- mix1, mixdf <- rbind(mixdf, mix1))
              ifelse(i==1, mixtdf <- mix1t, mixtdf <- rbind(mixtdf, mix1t))
            }
            cat('Mixture components constructed \n')
            # convert some data points to noise for Gaussian mixtures
            aux <- FALSE
            while (!aux){
              noisevals <- sample.int(nrow, size=round(noise*nrow))
              if (length(unique(mixdf[noisevals, dim+1]))==clusters){
                aux <- TRUE
              }
            }
            # convert some data points to noise for t mixtures
            aux <- FALSE
            while (!aux){
              noisevalst <- sample.int(nrow, size=round(noise*nrow))
              if (length(unique(mixtdf[noisevalst, dim+1]))==clusters){
                aux <- TRUE
              }
            }
            cat('Uniqueness checks done \n')
            # remove some observations from the rest of the data points
            mixdfmissing <- ampute(mixdf[-c(noisevals), 1:dim],
                                   prop=missing, mech='MAR', type='MID',
                                   bycases=TRUE)$amp
            mixdfmissing$label <- mixdf[-c(noisevals), dim+1]
            mixtdfmissing <- ampute(mixtdf[-c(noisevalst), 1:dim],
                                    prop=missing, mech='MAR', type='MID',
                                    bycases=TRUE)$amp
            mixtdfmissing$label <- mixtdf[-c(noisevalst), dim+1]
            cat('Observations removed \n')
            # generate noise
            noisevec <- c()
            for (i in 1:dim){
              noisevec <- c(noisevec, runif(n=noise*nrow, min=-noisethr, max=noisethr))
            }
            noisetvec <- c()
            for (i in 1:dim){
              noisetvec <- c(noisetvec, runif(n=noise*nrow, min=-noisethr, max=noisethr))
            }
            noisedf <- cbind(matrix(noisevec, ncol=dim), mixdf[noisevals, dim+1])
            colnames(noisedf) <- colnames(mixdfmissing)
            finaldf <- rbind(mixdfmissing, noisedf)
            noisetdf <- cbind(matrix(noisetvec, ncol=dim), mixtdf[noisevalst, dim+1])
            colnames(noisetdf) <- colnames(mixtdfmissing)
            finaltdf <- rbind(mixtdfmissing, noisetdf)
            cat('Final data frames constructed \n')
            ### clustering
            aux <- TRUE # this should handle errors
            # MCNM CLUSTERING
            if (aux){
              tryCatch({mcnmclus <- MCNM(finaldf[,c(1:dim)], G=clusters, max_iter=500)
              mcnmari <- adjustedRandIndex(finaldf$label, mcnmclus$clusters)
              mcnmarin <- adjustedRandIndex(finaldf[c(1:nrow(mixdfmissing)), dim+1],
                                            mcnmclus$clusters[c(1:nrow(mixdfmissing))])},
              error = function(cond){
                print("NOT OK")
                message(cond)
                aux <<- FALSE
                return(aux)})
            }
            
            # MNM CLUSTERING
            if (aux){
              tryCatch({mnmclus <- MNM(finaldf[,c(1:dim)], G=clusters, max_iter=500)
              mnmari <- adjustedRandIndex(finaldf$label, mnmclus$clusters)
              mnmarin <- adjustedRandIndex(finaldf[c(1:nrow(mixdfmissing)), dim+1],
                                           mnmclus$clusters[c(1:nrow(mixdfmissing))])},
              error = function(cond){
                print("NOT OK")
                message(cond)
                aux <<- FALSE
                return(aux)})
            }
            # STUDENT-t CLUSTERING
            if (aux){
              tryCatch({mtmclus <- MtM(finaldf[,c(1:dim)], G=clusters, max_iter=500)
              mtmari <- adjustedRandIndex(finaldf$label, mtmclus$clusters)
              mtmarin <- adjustedRandIndex(finaldf[c(1:nrow(mixdfmissing)), dim+1],
                                           mtmclus$clusters[c(1:nrow(mixdfmissing))])},
              error = function(cond){
                print("NOT OK")
                message(cond)
                aux <<- FALSE
                return(aux)})
            }
            # MCNM CLUSTERING
            if (aux){
              tryCatch({mcnmtclus <- MCNM(finaltdf[,c(1:dim)], G=clusters, max_iter=500)
              mcnmtari <- adjustedRandIndex(finaltdf$label, mcnmtclus$clusters)
              mcnmtarin <- adjustedRandIndex(finaltdf[c(1:nrow(mixtdfmissing)), dim+1],
                                             mcnmtclus$clusters[c(1:nrow(mixtdfmissing))])},
              error = function(cond){
                print("NOT OK")
                message(cond)
                aux <<- FALSE
                return(aux)}
              )
            }
            # STUDENT-t CLUSTERING
            if (aux){
              tryCatch({mtmtclus <- MtM(finaltdf[,c(1:dim)], G=clusters, max_iter=500)
              mtmtari <- adjustedRandIndex(finaltdf$label, mtmtclus$clusters)
              mtmtarin <- adjustedRandIndex(finaltdf[c(1:nrow(mixtdfmissing)), dim+1],
                                            mtmtclus$clusters[c(1:nrow(mixtdfmissing))])},
              error = function(cond){
                print("NOT OK")
                message(cond)
                aux <<- FALSE
                return(aux)}
              )
            }
            if (aux){
              sumpimcnm <- 0
              sumpimnm <- 0
              sumpimtm <- 0
              summumcnm <- 0
              summumnm <- 0
              summumtm <- 0
              sumvarmcnm <- 0
              sumvarmnm <- 0
              sumvarmtm <- 0
              sumpimcnmt <- 0
              sumpimtmt <- 0
              summumcnmt <- 0
              summumtmt <- 0
              sumvarmcnmt <- 0
              sumvarmtmt <- 0
              sumdf <- 0
              for (i in 1:clusters){
                # values for Gaussians
                sumpimcnm <- sumpimcnm + abs(mcnmclus$pi[i] - pivec[i])
                sumpimnm <- sumpimnm + abs(mnmclus$pi[i] - pivec[i])
                sumpimtm <- sumpimtm + abs(mtmclus$pi[i] - pivec[i])
                summumcnm <- summumcnm + chebyshevDistance(mcnmclus$mu[i,], rep(means[i], dim))
                summumnm <- summumnm + chebyshevDistance(mnmclus$mu[i,], rep(means[i], dim))
                summumtm <- summumtm + chebyshevDistance(mtmclus$mu[i,], rep(means[i], dim))
                varmat <- matrix(varsvec[((i-1)*dim^2+1):(i*dim^2)], nrow=dim)
                sumvarmcnm <- sumvarmcnm + norm(mcnmclus$sigma[,,i] - varmat, "M")
                sumvarmnm <- sumvarmnm + norm(mnmclus$sigma[,,i] - varmat, "M")
                sumvarmtm <- sumvarmtm + norm(mtmclus$sigma[,,i] - varmat, "M")
                # values for t-distribution
                sumpimcnmt <- sumpimcnmt + abs(mcnmtclus$pi[i] - pivec[i])
                sumpimtmt <- sumpimtmt + abs(mtmtclus$pi[i] - pivec[i])
                summumcnmt <- summumcnmt + chebyshevDistance(mcnmtclus$mu[i,], rep(means[i], dim))
                summumtmt <- summumtmt + chebyshevDistance(mtmtclus$mu[i,], rep(means[i], dim))
                sumvarmcnmt <- sumvarmcnmt + norm(mcnmtclus$sigma[,,i] - varmat, "M")
                sumvarmtmt <- sumvarmtmt + norm(mtmtclus$sigma[,,i] - varmat, "M")
                sumdf <- sumdf + abs(mtmtclus$df[i] - dfs[i])
              }
              cat('Values calculated \n')
              mcnmrow <- data.frame("seed"=j, "nClust"=clusters, "nrows"=nrow,
                                    "ncols"=dim, "noisep"=noise, "missingp"=missing,
                                    "pidist"=sumpimcnm, "mudist"=summumcnm,
                                    "vardist"=sumvarmcnm, "method"="MCNM",
                                    "ARIFull"=mcnmari, "ARINoNoise"=mcnmarin)
              mnmrow <- data.frame("seed"=j, "nClust"=clusters, "nrows"=nrow,
                                   "ncols"=dim, "noisep"=noise, "missingp"=missing,
                                   "pidist"=sumpimnm, "mudist"=summumnm,
                                   "vardist"=sumvarmnm, "method"="MNM",
                                   "ARIFull"=mnmari, "ARINoNoise"=mnmarin)
              mtmrow <- data.frame("seed"=j, "nClust"=clusters, "nrows"=nrow,
                                   "ncols"=dim, "noisep"=noise, "missingp"=missing,
                                   "pidist"=sumpimtm, "mudist"=summumtm,
                                   "vardist"=sumvarmtm, "method"="MtM",
                                   "ARIFull"=mtmari, "ARINoNoise"=mtmarin)
              mcnmrowt <- data.frame("seed"=j, "nClust"=clusters, "nrows"=nrow,
                                     "ncols"=dim, "noisep"=noise, "missingp"=missing,
                                     "pidist"=sumpimcnmt, "mudist"=summumcnmt,
                                     "vardist"=sumvarmcnmt, "dfdist"=NA,
                                     "method"="MCNM", "ARIFull"=mcnmtari,
                                     "ARINoNoise"=mcnmtarin)
              mtmrowt <- data.frame("seed"=j, "nClust"=clusters, "nrows"=nrow,
                                    "ncols"=dim, "noisep"=noise, "missingp"=missing,
                                    "pidist"=sumpimtmt, "mudist"=summumtmt,
                                    "vardist"=sumvarmtmt, "dfdist"=sumdf,
                                    "method"="MtM", "ARIFull"=mtmtari,
                                    "ARINoNoise"=mtmtarin)
              resdf <- rbind(resdf, mcnmrow)
              resdf <- rbind(resdf, mnmrow)
              resdf <- rbind(resdf, mtmrow)
              save(resdf, file='resdf.RData')
              restdf <- rbind(restdf, mcnmrowt)
              restdf <- rbind(restdf, mtmrowt)
              save(restdf, file='restdf.RData')
              cat('Iteration', run, '/', num_iters,'complete \n')
              run <- run+1
              count <- count+1
            }
            j <- j+1
          }
          j <- j+1
        }
        j <- j+1
      }
      j <- j+1
    }
    j <- j+1
  }
  j <- j+1
}
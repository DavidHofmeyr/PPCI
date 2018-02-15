#### A modification of R's default optimisation to allow, in particular, verbosity levels used
#### in Hofmeyr and Pavlidis (2017)

### function ppclust.optim() finds the optimal projection vector for clustering using one of the
### objectives in Hofmeyr and Pavlidis (2017)
## arguments:
# v0 = vector. Initial projection vector
# f = function. objective function to be optimised
# df = function. gradient of objective
# X = matrix. dataset being clustered
# params = list. parameters used in f and df
# getb = function. determines the optimal hyperplane orthogonal to a vector v (the value of b in the description of a hyperplane in the paper)
# verbosity = integer. determines the verbosity level of the optimisation. see paper for details
# labels = vector. only used if verbosity > 0. plots of points with different labels are given different colours. the labels are
#           not used to aid the optimisation, only for illustrative purposes
# method = character. describes the objective being used. one of c('MCDC', 'MDH', 'NCutH')
# maxit = integer. maximum number of iterations in optimisation.
# ftol = double. realtive tolerance for termination of optimisation

## output is the optimal projection vector

ppclust.optim <- function(v0, f, df, X, params, getb, verbosity, labels, method, maxit, ftol){
  if(verbosity==0){
    if(method=='MCDC') optim(v0, f, df, X, params, method='BFGS', control = list(fnscale = -1e-4, maxit = maxit, reltol = ftol))
    else optim(v0, f, df, X, params, method='BFGS', control = list(maxit = maxit, reltol = ftol))
  }
  else{
    iter = 0
    if(verbosity==3){
      dirname <- paste(method, ' plots.', format(Sys.time(), "%m_%d_%y_%Ih%Mm%Ss"), sep = '')
      if(method=='MDH') dirname <- paste(dirname, '. alpha=', params$alpha, sep = '')
      dir.create(dirname)
      count <- 0
    }
    if(ncol(X)>2) v2 <- rARPACK::eigs_sym(cov(X), 2)$vectors[,2]
    else v2 <- eigen(cov(X))$vectors[,2]
    v0 <- list(par = v0)
    for(i in 1:30){
      v2 <- v2-v2%*%v0$par*v0$par
      xp <- X%*%cbind(v0$par, v2)
      if(method=='MDH') den <- density(xp[,1], bw = params$h)
      else den <- density(xp[,1])
      if(verbosity==3){
        count <- count + 1
        pdf(paste('~/', dirname, '/', method, ' plot-', count, '.pdf', sep = ''))
      }
      if(is.null(labels)) plot(xp, col = rgb(.5, .5, .5), xlab = '', ylab = '', main = '')
      else plot(xp, col = labels, xlab = '', ylab = '', main = '')
      lines(den$x, min(xp[,2])+.5*den$y/max(den$y)*(max(xp[,2])-min(xp[,2])), lwd = 2)
      if(method=='MDH'){
        s <- sd(xp[,1])
        pen <- den$y + params$C*((den$x<(-params$alpha*s))*(-params$alpha*s-den$x)^2 + (den$x>(params$alpha*s))*(den$x-params$alpha*s)^2)
        lines(den$x, min(xp[,2])+.5*pen/max(den$y)*(max(xp[,2])-min(xp[,2])), lty = 2)
      }
      bcrit <- getb(v0$par, X, params)
      abline(v=bcrit, lwd = 2, col = 2)
      if(verbosity>1){
        SR <- ifelse(is.null(labels), NA, round(success_ratio(xp[,1]<bcrit, labels), 3))
        if(method=='MDH'){
          reldep <- round(md_reldepth(v0$par, X, params), 3)
          denval <- round(f_md(v0$par, X, params), 3)
          title(main = paste('alpha = ', params$alpha, '. density = ', denval, '. relative depth = ', reldep, '. success ratio = ', SR, sep = ''))
        }
        else if(method=='NCutH'){
          NC <- round(f_ncut(v0$par, X, params), 3)
          title(main = paste('normalised cut = ', NC, '. success ratio = ', SR, sep = ''))
        }
        else{
          VR <- round(f_mc(v0$par, X, params), 3)
          title(main = paste('variance ratio = ', VR, '. success ratio = ', SR, sep = ''))
        }
      }
      if(verbosity==3) dev.off()
      if(method=='MCDC') v0 <- optim(v0$par, f, df, X, params, method = 'BFGS', control = list(fnscale = -1e-4, maxit=3, reltol = ftol))
      else v0 <- optim(v0$par, f, df, X, params, method = 'BFGS', control = list(maxit=3, reltol = ftol))
      Sys.sleep(.05)
      v0$par <- v0$par/norm_vec(v0$par)
      iter <- iter + v0$counts[2]
      if(v0$counts[2]==1 || i==30 || iter > maxit) break
    }
    v0
  }
}



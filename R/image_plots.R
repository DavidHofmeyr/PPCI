optidigits_mean_images <- function(clusters){
  op <- par(no.readonly = TRUE)
  T = table(PPCI::optidigits$c, clusters)
  f = function(p) -sum(diag(T[,p[1:min(10, max(clusters))]]))
  df = function(p){
    sm = sample(1:length(p), 2)
    rep1 = p[sm[1]]
    rep2 = p[sm[2]]
    p[sm[2]] = rep1
    p[sm[1]] = rep2
    p
  }
  p = optim(1:max(clusters), f, df, method = 'SANN')$par
  par(mfrow = c(1, max(clusters)))
  for(i in 1:max(clusters)){
    par(mar = c(0, 0, 0, 0))
    image(matrix(colMeans(PPCI::optidigits$x[which(clusters==p[i]),]), 8, 8)[,8:1], col = rgb((1:20)/20, (1:20)/20, (1:20)/20), xaxt = 'n', yaxt = 'n')
  }
  par(op)
}

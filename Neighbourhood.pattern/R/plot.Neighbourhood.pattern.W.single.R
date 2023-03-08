plot.Neighbourhood.pattern.W.single=function (a, b, n) 
{
  library(ggplot2)
  Neighbourhood.pattern.W.single = function(a, b, n) {
    Neighbourhood.single = function(a, b, n) {
      c = b
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, ] - a[1, ])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      d = cbind(b, d)
      d = subset(d, d > 0)
      d = d[order(d[, 3])[1:n], 1:3]
      colnames(d) = c("x", "y", "Distance")
      d
    }
    Neighbourhood = Neighbourhood.single(a, b, n)[, 1:2]
    threshold = matrix(cos(2 * pi/(n + 1)), 1, 2 * n)
    slope = matrix(NA, n, 1)
    e = matrix(NA, n, n)
    rb = matrix(NA, n, n)
    residual = matrix(NA, n, 1)
    stand1 = (Neighbourhood.single(a, b, n)[, 1] - a[1, 1])/Neighbourhood.single(a, 
                                                                                 b, n)[, 3]
    stand1[stand1 == 0] = 1e-08
    stand2 = (Neighbourhood.single(a, b, n)[, 2] - a[1, 2])/Neighbourhood.single(a, 
                                                                                 b, n)[, 3]
    standard = cbind(stand1, stand2)
    for (i in 1:n) {
      slope[i, ] = standard[i, 2]/standard[i, 1]
      for (j in 1:n) {
        residual[j, ] = slope[i, ] * standard[j, 1] - 
          standard[j, 2]
      }
      e[, i] = residual
      e[abs(e) < 1e-09] = 0
    }
    for (i in 3:(n + 2)) {
      standard.e = cbind(standard, e)
      pe = subset(standard.e, standard.e[, i] >= 0, c(1, 
                                                      2, i))
      pea = pe[, 1:2] %*% standard[(i - 2), ]
      pea = pea[order(pea, decreasing = T)]
      ne = subset(standard.e, standard.e[, i] < 0, c(1, 
                                                     2, i))
      nea = ne[, 1:2] %*% standard[(i - 2), ]
      nea = nea[order(nea)]
      rbe = rbind(as.matrix(pea), as.matrix(nea))
      rb[, (i - 2)] = rbe
    }
    rbn2 = c(rb[2, ], rb[n, ])
    angle = rbn2 - threshold
    ang = subset(matrix(angle), matrix(angle) > 0)
    num = length(ang)
    W = num/(2 * n)
    outcome = list(a = a, Neighbourhood = Neighbourhood, 
                   W = W)
    outcome
  }
  Ne = Neighbourhood.pattern.W.single(a, b, n)
  center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[, 
                                                          2], each = n)), c(1:n))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "group")
  Neigh = cbind(Ne$Neighbourhood, c(1:n))
  colnames(Neigh) = c("x", "y", "group")
  unitg = rbind(center, Neigh)
  max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[, 
                                                                         2] - Ne$a[, 2]))
  p = ggplot() + geom_line(data = unitg, aes(x = x, y = y, 
                                             group = group), linetype = 2, size = 1) + geom_point(data = Ne$a, 
                                                                                                  aes(x, y), size = 10, color = "red4", alpha = 0.7) + 
    geom_point(data = Ne$Neighbourhood, aes(x, y), size = 10, 
               color = "green4", alpha = 0.7)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Ne$Neighbourhood, 
                                                                    aes(x, y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - max, Ne$a[, 1] + max), y = c(Ne$a[, 
                                                                   2] - max, Ne$a[, 2] + max))
  p = p + annotate("text", x = Ne$a[, 1] - max + 0.125 * 
                     (max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("W=", 
                                                                                round(Ne$W, 3)))
  p + theme_bw()
}
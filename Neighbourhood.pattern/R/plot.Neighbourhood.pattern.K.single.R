plot.Neighbourhood.pattern.K.single=function (a, b, n) 
{
  library(ggplot2)
  Neighbourhood.pattern.K.single = function(a, b, n) {
    Neighbourhood.single1 = function(a, b, n) {
      c = b[, 1:2]
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      d = cbind(b, d)
      d = subset(d, d > 0)
      d = d[order(d[, 4])[1:n], 1:4]
      colnames(d) = c("x", "Y", "Size", 
                      "Distance")
      d
    }
    Nei.tree = Neighbourhood.single1(a, b, n)
    Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - a[, 3])
    colnames(Nei.dif.high) = c("x", "Y", "Size", 
                               "Distance", "Size_dif ")
    Neighbourhood = Nei.dif.high
    Neighbourhood1 = Neighbourhood
    Neighbourhood[which(Neighbourhood[, 5] <= 0), 5] = NA
    Neighbourhood1[which(Neighbourhood1[, 5] <= 0), 5] = 1e-08
    k1 = Neighbourhood1[, 4]/Neighbourhood1[, 5]
    k1 = k1/(k1 + 1)
    k2 = as.data.frame(k1)
    K = sum(k2, na.rm = T)
    K = K/n
    Neighbourhood_pattern_K = K
    outcome = list(a = a, Neighbourhood = Neighbourhood, 
                   Neighbourhood_pattern_K = Neighbourhood_pattern_K)
    outcome
  }
  Ne = Neighbourhood.pattern.K.single(a, b, n)
  center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[, 
                                                          2], each = n), rep(Ne$a[, 3], each = n)), c(1:n))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "size", 
                       "group")
  Neigh = cbind(Ne$Neighbourhood, c(1:n))
  Neigh = Neigh[, c(1, 2, 3, 6)]
  colnames(Neigh) = c("x", "y", "size", "group")
  unitg = rbind(center, Neigh)
  max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[, 
                                                                         2] - Ne$a[, 2]))
  p = ggplot()
  p = p + geom_line(data = unitg, aes(x = x, y = y, group = group), 
                    linetype = 2, size = 1)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 13/max(unitg[, 
                                                                 3]) * a[, 3] + 2, color = "red4", alpha = 0.7) + 
    geom_point(data = Neigh, aes(x, y), size = 13/max(unitg[, 
                                                            3]) * Neigh[, 3] + 2, color = "green4", alpha = 0.7)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Neigh, 
                                                                    aes(x, y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - 1.1 * max, Ne$a[, 1] + 1.1 * 
                       max), y = c(Ne$a[, 2] - 1.1 * max, Ne$a[, 2] + 1.1 * 
                                     max))
  p = p + annotate("text", x = Ne$a[, 1] - max + 0.3 * 
                     (max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("Neighbourhood_pattern_K=", 
                                                                                round(Ne$Neighbourhood_pattern_K, 3)))
  p + theme_bw()
}
plot.Neighbourhood.pattern.U.single=function (a, b, n) 
{
  library(ggplot2)
  Neighbourhood.pattern.U.single = function(a, b, n) {
    Neighbourhood.single = function(a, b, n) {
      a1 = a
      b1 = b
      a = a[, 1:2]
      b = b[, 1:2]
      c = b
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, ] - a[1, ])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      d = cbind(b, d, b1[, 3])
      d = subset(d, d > 0)
      d = d[order(d[, 3])[1:n], 1:4]
      d = d[, c(1:2, 4)]
      colnames(d) = c("x", "y", "Size")
      d
    }
    Neighbourhood = Neighbourhood.single(a, b, n)
    U = subset(Neighbourhood, Size > a[, 3])
    U = nrow(U)/n
    outcome = list(a = a, Neighbourhood = Neighbourhood, 
                   U = U)
    outcome
  }
  Ne = Neighbourhood.pattern.U.single(a, b, n)
  center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[, 
                                                          2], each = n), rep(Ne$a[, 3], each = n)), c(1:n))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "size", 
                       "group")
  Neigh = cbind(Ne$Neighbourhood, c(1:n))
  colnames(Neigh) = c("x", "y", "size", "group")
  unitg = rbind(center, Neigh)
  big = subset(Ne$Neighbourhood, Size > a[, 3])
  small = subset(Ne$Neighbourhood, Size <= a[, 3])
  max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[, 
                                                                         2] - Ne$a[, 2]))
  p = ggplot() + geom_line(data = unitg, aes(x = x, y = y, 
                                             group = group), linetype = 2, size = 1) + geom_point(data = Ne$a, 
                                                                                                  aes(x, y), size = 13/max(unitg[, 3]) * a[, 3] + 2, color = "red4", 
                                                                                                  alpha = 0.7) + geom_point(data = big, aes(x, y), size = 13/max(unitg[, 
                                                                                                                                                                       3]) * big[, 3] + 2, color = "green4", alpha = 0.7) + 
    geom_point(data = small, aes(x, y), size = 13/max(unitg[, 
                                                            3]) * small[, 3] + 2, color = "grey", alpha = 0.7)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = big, 
                                                                    aes(x, y), size = 2) + geom_point(data = small, aes(x, 
                                                                                                                        y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - max, Ne$a[, 1] + max), y = c(Ne$a[, 
                                                                   2] - max, Ne$a[, 2] + max))
  p = p + annotate("text", x = Ne$a[, 1] - max + 0.125 * 
                     (max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("U=", 
                                                                                round(Ne$U, 3)))
  p + theme_bw()
}
plot.Neighbourhood.pattern.M.single=function (a, b, n) 
{
  library(ggplot2)
  Neighbourhood.pattern.M.single = function(a, b, n) {
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
      colnames(d) = c("x", "y", "Species")
      rownames(d) = 1:n
      d
    }
    Neighbourhood = Neighbourhood.single(a, b, n)
    Neighbourhood1 = subset(Neighbourhood, Species != a[, 
                                                        3])
    M = nrow(Neighbourhood1)/n
    output = list(a = a, Neighbourhood = Neighbourhood, M = M)
    output
  }
  Ne = Neighbourhood.pattern.M.single(a, b, n)
  center = cbind(cbind(rep(a[, 1], each = n), rep(a[, 2], each = n)), 
                 c(1:n))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "group")
  center1 = center
  center1[, 4] = "red4"
    colnames(center1) = c("x", "y", "group", 
                          "col")
    Neigh = cbind(Ne$Neighbourhood[, 1:3], c(1:n))
    Neigh1 = Neigh
    for (i in 1:n) {
      if (Neigh[i, 3] == a[, 3]) {
        Neigh1[i, 5] = "red4"
      }
      else {
        Neigh1[i, 5] = "green4"
      }
    }
    colnames(Neigh) = c("x", "y", "Species", 
                        "group")
    colnames(Neigh1) = c("x", "y", "Species", 
                         "group", "col")
    Neigh = Neigh[, c(1:2, 4)]
    unitg = rbind(center, Neigh)
    max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[, 
                                                                           2] - Ne$a[, 2]))
    p = ggplot() + geom_line(data = unitg, aes(x = x, y = y, 
                                               group = group), linetype = 2, size = 1) + geom_point(data = a, 
                                                                                                    aes(x, y), size = 10, color = "red4", alpha = 0.7) + 
      geom_point(data = Neigh1, aes(x, y), size = 10, color = Neigh1$col, 
                 alpha = 0.7)
    p = p + geom_point(data = a, aes(x, y), size = 2, color = "black") + 
      geom_point(data = Neigh1, aes(x, y), size = 2, color = "black")
    p = p + lims(x = c(Ne$a[, 1] - max, Ne$a[, 1] + max), y = c(Ne$a[, 
                                                                     2] - max, Ne$a[, 2] + max))
    p = p + annotate("text", x = Ne$a[, 1] - max + 0.125 * 
                       (max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("M=", 
                                                                                  round(Ne$M, 3)))
    p + theme_bw()
}
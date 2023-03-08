plot.Neighbourhood.pattern.C.single=function (a, b, n) 
{
  library(ggplot2)
  library(ggforce)
  Neighbourhood.pattern.C.single = function(a, b, n) {
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
      colnames(d) = c("x", "y", "distance", 
                      "size")
      d
    }
    Neighbourhood = Neighbourhood.single(a, b, n)
    C = subset(Neighbourhood, distance < (Neighbourhood$size + 
                                            a$size)/2)
    C = nrow(C)/n
    outcome = list(a = a, Neighbourhood = Neighbourhood, 
                   C = C)
    outcome
  }
  Ne = Neighbourhood.pattern.C.single(a, b, n)
  center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[, 
                                                          2], each = n), rep(Ne$a[, 3], each = n)), c(1:n))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "size", 
                       "group")
  Neigh = cbind(Ne$Neighbourhood, c(1:n))
  Neigh = Neigh[, c(1, 2, 4, 5)]
  colnames(Neigh) = c("x", "y", "size", "group")
  unitg = rbind(center, Neigh)
  small = subset(Ne$Neighbourhood, distance < (Ne$Neighbourhood$size + 
                                                 a$size)/2)
  big = subset(Ne$Neighbourhood, distance >= (Ne$Neighbourhood$size + 
                                                a$size)/2)
  max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[, 
                                                                         2] - Ne$a[, 2]))
  p = ggplot()
  p = p + geom_line(data = unitg, aes(x = x, y = y, group = group), 
                    linetype = 2, size = 1) + geom_circle(data = Ne$a, aes(x0 = x, 
                                                                           y0 = y, r = size/2), fill = "red4", color = "black", 
                                                          alpha = 0.7) + geom_circle(data = big, aes(x0 = x, y0 = y, 
                                                                                                     r = size/2), color = "black", fill = "green4", 
                                                                                     alpha = 0.7) + geom_circle(data = small, aes(x0 = x, 
                                                                                                                                  y0 = y, r = size/2), fill = "grey", color = "black", 
                                                                                                                alpha = 0.7)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = big, 
                                                                    aes(x, y), size = 2) + geom_point(data = small, aes(x, 
                                                                                                                        y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - 1.5 * max, Ne$a[, 1] + 1.3 * 
                       max), y = c(Ne$a[, 2] - 1.5 * max, Ne$a[, 2] + 1.3 * 
                                     max))
  p = p + annotate("text", x = Ne$a[, 1] - 1.3 * max + 
                     0.125 * (max), y = Ne$a[, 2] + 1.3 * max - 0.125 * (max), 
                   label = paste0("C=", round(Ne$C, 3)))
  p + theme_bw()
}
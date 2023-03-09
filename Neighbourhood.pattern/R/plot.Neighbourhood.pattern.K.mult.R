plot.Neighbourhood.pattern.K.mult=function (a, b, n) 
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.K.mult = function(a, b, n) {
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
      Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - a[, 
                                                       3])
      colnames(Nei.dif.high) = c("x", "Y", 
                                 "Size", "Distance", "Size_dif ")
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
    d = matrix(NA, nrow(a), 4)
    pb = tkProgressBar("进度", "已完成 %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.pattern.K.single(a[j, 
      ], b, n)$Neighbourhood_pattern_K))
      info = sprintf("已完成 %d%%", round(j * 
                                         100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Size", "Neighbourhood_pattern_K")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  plot.data = Neighbourhood.pattern.K.mult(a, b, n)
  p = ggplot(plot.data, aes(x, y)) + geom_point(size = 12/max(plot.data[, 
                                                                        4]) * plot.data[, 4] + 3, shape = 1, color = "red")
  p1 = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw()
  p2 = p1 + xlab("") + ylab("") + labs(title = "The neighbourhood pattern K of forest")
  p2
}

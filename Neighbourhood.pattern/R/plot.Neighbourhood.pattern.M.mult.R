plot.Neighbourhood.pattern.M.mult=function (a, b, n) 
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.M.mult = function(a, b, n) {
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
        colnames(d) = c("x", "Y", "Species")
        rownames(d) = 1:n
        d
      }
      Neighbourhood = Neighbourhood.single(a, b, n)
      Neighbourhood1 = subset(Neighbourhood, Species != 
                                a[, 3])
      M = nrow(Neighbourhood1)/n
      output = list(a = a, Neighbourhood = Neighbourhood, 
                    M = M)
      output
    }
    d = matrix(NA, nrow(a), 4)
    pb = tkProgressBar("进度", "已完成 %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.pattern.M.single(a[j, 
      ], b, n)$M))
      info = sprintf("已完成 %d%%", round(j * 
                                         100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("X", "Y", "Species", 
                    "M")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  plot.data = Neighbourhood.pattern.M.mult(a, b, n)
  for (i in c(1:2, 4)) {
    plot.data[, i] = as.numeric(as.character(plot.data[, 
                                                       i]))
  }
  p = ggplot(plot.data, aes(X, Y)) + geom_point(size = 12/max(plot.data[, 
                                                                        4]) * plot.data[, 4] + 3, shape = 1, color = "red")
  p1 = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw()
  p2 = p1 + xlab("") + ylab("") + labs(title = "The neighbourhood pattern M of forest")
  p2
}

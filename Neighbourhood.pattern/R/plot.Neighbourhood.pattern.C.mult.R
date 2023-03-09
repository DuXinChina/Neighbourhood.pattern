plot.Neighbourhood.pattern.C.mult=function (a, b, n) 
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.C.mult = function(a, b, n) {
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
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("进度", "已完成 %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Neighbourhood.pattern.C.single(a[j, 
      ], b, n)$C))
      info = sprintf("已完成 %d%%", round(j * 
                                         100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "C")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  plot.data = Neighbourhood.pattern.C.mult(a, b, n)
  p = ggplot(plot.data, aes(x, y)) + geom_point(size = 12/max(plot.data[, 
                                                                        3]) * plot.data[, 3] + 3, shape = 1, color = "red")
  p1 = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw()
  p2 = p1 + xlab("") + ylab("") + labs(title = "The neighbourhood pattern C of forest")
  p2
}
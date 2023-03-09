plot.Neighbourhood.pattern.W.dif_sp=function (a, b, n) 
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.mult = function(a, b, n) {
    a = a[, 1:2]
    b = b[, 1:2]
    Neighbourhood.pattern.single = function(a, b, n) {
      Neighbourhood.single = function(a, b, n) {
        c = b
        for (i in 1:nrow(b)) {
          c[i, ] = (b[i, ] - a[1, ])^2
          d = (c[, 1] + c[, 2])^(1/2)
        }
        d = cbind(b, d)
        d = subset(d, d > 0)
        d = d[order(d[, 3])[1:n], 1:3]
        colnames(d) = c("x", "Y", "Distance")
        rownames(d) = 1:n
        d
      }
      threshold = matrix(cos(2 * pi/(n + 1)), 1, 2 * n)
      slope = matrix(NA, n, 1)
      e = matrix(NA, n, n)
      rb = matrix(NA, n, n)
      residual = matrix(NA, n, 1)
      stand1 = (Neighbourhood.single(a, b, n)[, 1] - a[1, 
                                                       1])/Neighbourhood.single(a, b, n)[, 3]
      stand1[stand1 == 0] = 1e-08
      stand2 = (Neighbourhood.single(a, b, n)[, 2] - a[1, 
                                                       2])/Neighbourhood.single(a, b, n)[, 3]
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
        pe = subset(standard.e, standard.e[, i] >= 0, 
                    c(1, 2, i))
        pea = pe[, 1:2] %*% standard[(i - 2), ]
        pea = pea[order(pea, decreasing = T)]
        ne = subset(standard.e, standard.e[, i] < 0, 
                    c(1, 2, i))
        nea = ne[, 1:2] %*% standard[(i - 2), ]
        nea = nea[order(nea)]
        rbe = rbind(as.matrix(pea), as.matrix(nea))
        rb[, (i - 2)] = rbe
      }
      rbn2 = c(rb[2, ], rb[n, ])
      angle = rbn2 - threshold
      ang = subset(matrix(angle), matrix(angle) > 0)
      num = length(ang)
      vaule = num/(2 * n)
      vaule
    }
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("进度", "已完成 %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.pattern.single(a[j, 
      ], b, n)))
      info = sprintf("已完成 %d%%", round(j * 
                                         100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("X", "Y", "Neighbourhood Pattern W")
    rownames(d) = 1:j
    d
  }
  plot.data = Neighbourhood.pattern.mult(a, b, n)
  plot.data = as.data.frame(plot.data)
  plot.data = cbind(plot.data, a$Species)
  colnames(plot.data) = c("X", "Y", "Neighbourhood Pattern W", 
                          "Species")
  p = ggplot() + geom_point(data = plot.data, aes(X, Y, color = Species), 
                            size = 12/max(plot.data[, 3]) * plot.data[, 3] + 1, shape = 1)
  p1 = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw()
  p2 = p1 + xlab("") + ylab("") + labs(title = "The neighbourhood pattern W of diffence species")
  p2
}

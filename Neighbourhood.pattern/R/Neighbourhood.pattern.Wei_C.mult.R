Neighbourhood.pattern.Wei_C.mult=function (a, b) 
{
  library(tcltk)
  Neighbourhood.pattern.Wei_C.single = function(a, b) {
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
    C = Neighbourhood.pattern.C.single(a, b, 4)$C
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
      W = num/(2 * n)
      outcome = list(a = a, Neighbourhood = Neighbourhood, 
                     W = W)
      outcome
    }
    wei = Neighbourhood.pattern.W.single(a, b, 4)$W
    W = wei
    if (wei == 0) {
      Wei = 1
    }
    if (wei == 0.25) {
      Wei = 0.75
    }
    if (wei == 0.5) {
      Wei = 0.5
    }
    if (wei == 0.75) {
      Wei = 0.375
    }
    if (wei == 1) {
      Wei = 0.25
    }
    Wei_C = Wei * C
    Neighbourhood = Neighbourhood.pattern.C.single(a, b, 
                                                   4)$Neighbourhood
    output = list(a = a, Neighbourhood = Neighbourhood, C = C, 
                  W = W, Wei_C = Wei_C)
    output
  }
  d = matrix(NA, nrow(a), 3)
  pb = tkProgressBar("进度", "已完成 %", 0, 
                     100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Neighbourhood.pattern.Wei_C.single(a[j, 
    ], b)$Wei_C))
    info = sprintf("已完成 %d%%", round(j * 100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "Wei_C")
  rownames(d) = 1:nrow(a)
  d = as.data.frame(d)
  d
}

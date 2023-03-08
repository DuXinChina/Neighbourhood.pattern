Neighbourhood.pattern.Wei_K.single=function (a, b, n) 
{
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
  K = Neighbourhood.pattern.K.single(a, b, n)
  W = Neighbourhood.pattern.W.single(a, b, n)
  a = K$a
  Neighbourhood = K$Neighbourhood
  Neighbourhood_pattern_K = K$Neighbourhood_pattern_K
  Neighbourhood_pattern_Wei_K = K$Neighbourhood_pattern_K * 
    (W$W^(1 - K$Neighbourhood_pattern_K))
  W = W$W
  outcome = list(a = a, Neighbourhood = Neighbourhood, Neighbourhood_pattern_K = Neighbourhood_pattern_K, 
                 W = W, Neighbourhood_pattern_Wei_K = Neighbourhood_pattern_Wei_K)
  outcome
}
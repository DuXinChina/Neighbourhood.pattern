K.ISAA=function (minx, maxx, miny, maxy, b, indis, lag) 
{
  library(ape)
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  K.mult = function(a, b) {
    library(tcltk)
    Neighbourhood.K.single = function(a, b) {
      b1 = subset(b, x >= a[, 1] & y > a[, 2] & size > 
                    a[, 3])
      b2 = subset(b, x > a[, 1] & y <= a[, 2] & size > 
                    a[, 3])
      b3 = subset(b, x <= a[, 1] & y < a[, 2] & size > 
                    a[, 3])
      b4 = subset(b, x < a[, 1] & y >= a[, 2] & size > 
                    a[, 3])
      Neighbourhood.K.single1 = function(a, b) {
        b1 = subset(b, x >= a[, 1] & y > a[, 2] & size > 
                      a[, 3])
        b2 = subset(b, x > a[, 1] & y <= a[, 2] & size > 
                      a[, 3])
        b3 = subset(b, x <= a[, 1] & y < a[, 2] & size > 
                      a[, 3])
        b4 = subset(b, x < a[, 1] & y >= a[, 2] & size > 
                      a[, 3])
        Neighbourhood.single = function(a, b) {
          c = b[, 1:2]
          for (i in 1:nrow(b)) {
            c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
            d = (c[, 1] + c[, 2])^(1/2)
          }
          d = cbind(b, d)
          d = subset(d, d > 0)
          d = d[order(d[, 4])[1], 1:4]
          colnames(d) = c("x", "Y", "Size", 
                          "Distance")
          d
        }
        b1.Nei.tree = Neighbourhood.single(a, b1)
        b2.Nei.tree = Neighbourhood.single(a, b2)
        b3.Nei.tree = Neighbourhood.single(a, b3)
        b4.Nei.tree = Neighbourhood.single(a, b4)
        Nei.tree = rbind(b1.Nei.tree, b2.Nei.tree, b3.Nei.tree, 
                         b4.Nei.tree)
        Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - 
                               a[, 3])
        colnames(Nei.dif.high) = c("x", "Y", 
                                   "Size", "Distance", "Size_dif ")
        Neighbourhood = Nei.dif.high
        k1 = Neighbourhood[, 4]/Neighbourhood[, 5]
        k2 = as.data.frame(k1)
        k2[, 1][is.infinite(k2[, 1])] = NA
        K = sum(k2, na.rm = T)
        key = 2
        outcome = list(a = a, Neighbourhood = Neighbourhood, 
                       K = K, key = key)
        outcome
      }
      if (nrow(b1) != 0 & nrow(b2) != 0 & nrow(b3) != 0 & 
          nrow(b4) != 0) {
        Neighbourhood.K.single1(a, b)
      }
      else {
        key1 = 1
        key = 1
        out = list(key = key, key1 = key1)
        out
      }
    }
    d = matrix(NA, nrow(a), 4)
    pb <- tkProgressBar("Progress", "Percent complete %", 
                        0, 100)
    star_time <- Sys.time()
    for (j in 1:nrow(a)) {
      KEY = Neighbourhood.K.single(a[j, ], b)
      sw = KEY$key
      if (sw == 1) {
        warning(paste("Please check", c(j), "row"))
      }
      else {
        d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.K.single(a[j, 
        ], b)$K))
      }
      info <- sprintf("Percent complete %d%%", round(j * 
                                                       100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)", 
                                                    info), info)
    }
    end_time <- Sys.time()
    close(pb)
    run_time <- end_time - star_time
    colnames(d) = c("X", "Y", "Size", "K")
    rownames(d) = 1:nrow(a)
    if (is.na(sum(d)) == TRUE) {
      print("Please ensure that there are trees above the reference tree in all four quadrants")
      d
    }
    else {
      d
    }
  }
  nx = (maxx - minx - indis)%/%lag
  ny = (maxy - miny - indis)%/%lag
  n = min(nx, ny)
  I = matrix(NA, n, 1)
  Z = matrix(NA, n, 1)
  distance = matrix(NA, n, 1)
  x = seq(minx, maxx, by = indis)
  y = seq(miny, maxy, by = indis)
  xy = expand.grid(x, y)
  xy$size = 0
  colnames(xy) = c("x", "y", "size")
  vaule = K.mult(xy, b)
  vaule = vaule[complete.cases(vaule), ]
  vaule = vaule[, c(1, 2, 4)]
  colnames(vaule) = c("x", "y", "K")
  vaule = as.data.frame(vaule)
  for (i in 1:n) {
    vaule.dists = as.matrix(dist(cbind(vaule$x, vaule$y)))
    vaule.dists = (vaule.dists - indis)%/%(i * lag) * i * 
      lag + 0.5 * (i * lag + indis)
    vaule.dists.inv = 1/vaule.dists
    diag(vaule.dists.inv) = 0
    Moran = Moran.I(vaule$K, vaule.dists.inv)
    I[i] = Moran$observed
    Z[i] = (Moran$observed - Moran$expected)/Moran$sd
    distance[i] = lag * (i) + indis
    ISAA = cbind(distance, I, Z)
    colnames(ISAA) = c("lag", "Moran_I", "Z_Score")
  }
  yBL = max(ISAA[, 3])/max(ISAA[, 2])
  yBL = 1.1 * yBL
  ISAA = as.data.frame(ISAA)
  ISAAI = ISAA[, c(1, 2)]
  ISAAI$Index = c(" Moran_I")
  colnames(ISAAI) = c("Lag", "Moran_I", "Index")
  ISAAZ = ISAA[, c(1, 3)]
  ISAAZ$Index = c("Z_Score")
  colnames(ISAAZ) = c("Lag", "Moran_I", "Index")
  ISAAZ$Moran_I = ISAAZ$Moran_I/yBL
  ISAAIZ = rbind(ISAAI, ISAAZ)
  Z = ISAA[, 3]
  Zbef = matrix(NA, length(ISAA[, 3]) - 2, 1)
  Zbac = matrix(NA, length(ISAA[, 3]) - 2, 1)
  for (i in 1:(length(Z) - 2)) {
    Zbef[i] = Z[i + 1] - Z[i]
    Zbac[i] = Z[i + 2] - Z[i + 1]
  }
  ZBB = cbind(Zbef, Zbac)
  ZBB = cbind(ISAA$lag[-c(1, length(Z))], ZBB)
  ZBB = ZBB[which(ZBB[, 2] > 0 & ZBB[, 3] < 0), 1]
  p = ggplot(ISAAIZ, aes(x = Lag)) + geom_line(aes(y = Moran_I, 
                                                   group = Index)) + geom_point(aes(y = Moran_I, shape = Index), 
                                                                                size = 2)
  p = p + scale_y_continuous(sec.axis = sec_axis(~. * yBL, 
                                                 name = "Z_Score")) + theme_bw()
  p = p + geom_hline(aes(yintercept = 1.96/yBL), linetype = 5, 
                     col = "black") + geom_hline(aes(yintercept = -1.96/yBL), 
                                                 linetype = 5, col = "black") + geom_hline(aes(yintercept = 0/yBL), 
                                                                                           linetype = 1) + geom_vline(xintercept = ZBB, linetype = 5, 
                                                                                                                      col = "red") + theme(legend.title = element_blank(), 
                                                                                                                                           legend.position = c(0.8, 0.8))
  p
}
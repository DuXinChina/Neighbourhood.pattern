Scale.dependence.Wei_S.ISAA=function (minx, maxx, miny, maxy, b, seq, indis, lag, scale, MI) 
{
  library(ape)
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Scale.dependence.Wei_S.mult = function(a, b, scale, MI) {
    library(tcltk)
    Scale.dependence.Wei_S.single = function(a, b, scale, 
                                             MI) {
      Scale.dependence.S.single = function(a, b, scale, 
                                           MI) {
        Neighbourhood.single1 = function(a, b, scale) {
          c = b[, 1:2]
          for (i in 1:nrow(b)) {
            c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
            d = (c[, 1] + c[, 2])^(1/2)
          }
          d = cbind(b, d)
          d = subset(d, d < scale)
          colnames(d) = c("x", "Y", "Size", 
                          "Distance")
          d
        }
        Nei.tree = Neighbourhood.single1(a, b, scale)
        n = nrow(Nei.tree)
        Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - 
                               a[, 3])
        colnames(Nei.dif.high) = c("x", "Y", 
                                   "Size", "Distance", "Size_dif ")
        n.dif = subset(Nei.dif.high, Nei.dif.high$Size_dif > 
                         0)
        n.dif = nrow(n.dif)
        if (n.dif == 0) {
          Scale_dependence_S = 0
          Neighbourhood = Nei.dif.high
        }
        if (n.dif != 0) {
          Neighbourhood = Nei.dif.high
          Neighbourhood1 = Neighbourhood
          Neighbourhood[which(Neighbourhood[, 5] <= 0), 
                        5] = NA
          Neighbourhood1[which(Neighbourhood1[, 5] <= 
                                 0), 5] = 0
          Neighbourhood1[which(Neighbourhood1[, 4] == 
                                 0), 4] = 1e-08
          S1 = Neighbourhood1[, 5]/Neighbourhood1[, 4]
          S2 = as.data.frame(S1)
          S = sum(S2, na.rm = T)
          Scale_dependence_S = S
        }
        Scale_dependence_S = (atan(Scale_dependence_S/(MI * 
                                                         pi))/pi * 2)
        outcome = list(a = a, Neighbourhood = Neighbourhood, 
                       Scale_dependence_S = Scale_dependence_S)
        outcome
      }
      Scale.dependence.Simpson = function(a, b, scale) {
        Neighbourhood.single1 = function(a, b, scale) {
          c = b[, 1:2]
          for (i in 1:nrow(b)) {
            c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
            d = (c[, 1] + c[, 2])^(1/2)
          }
          d = cbind(b, d)
          d = subset(d, d > 0)
          d = subset(d, d < scale)
          colnames(d) = c("x", "y", "Size", 
                          "Distance")
          d
        }
        Nei.tree = Neighbourhood.single1(a, b, scale)
        Nei.tree
        Nei.tree[, 3] = Nei.tree[, 3] - a[, 3]
        Nei.tree = subset(Nei.tree, Size > 0)
        Nei.tree = subset(Nei.tree, Distance > 0)
        Nei.tree
        Simpn = sum(Nei.tree[, 3])
        Nei.tree1 = subset(Nei.tree, x >= a[, 1] & y > 
                             a[, 2])
        Nei.tree2 = subset(Nei.tree, x > a[, 1] & y <= 
                             a[, 2])
        Nei.tree3 = subset(Nei.tree, x <= a[, 1] & y < 
                             a[, 2])
        Nei.tree4 = subset(Nei.tree, x < a[, 1] & y >= 
                             a[, 2])
        Simpn1 = sum(Nei.tree1[, 3])
        Simpn2 = sum(Nei.tree2[, 3])
        Simpn3 = sum(Nei.tree3[, 3])
        Simpn4 = sum(Nei.tree4[, 3])
        simp_wei = sum((Simpn1/Simpn)^2, (Simpn2/Simpn)^2, 
                       (Simpn3/Simpn)^2, (Simpn4/Simpn)^2)
        if (is.nan(simp_wei) == T) 
          (simp_wei = 0.25)
        simp_wei
      }
      S = Scale.dependence.S.single(a, b, scale, MI)
      Levins_Simpson = Scale.dependence.Simpson(a, b, scale)
      a = S$a
      Neighbourhood = S$Neighbourhood
      Scale_dependence_S = S$Scale_dependence_S
      Scale_dependence_Wei_S = Scale_dependence_S * (0.25/Levins_Simpson)
      outcome = list(a = a, Neighbourhood = Neighbourhood, 
                     Scale_dependence_S = Scale_dependence_S, Levins_Simpson = Levins_Simpson, 
                     Scale_dependence_Wei_S = Scale_dependence_Wei_S)
      outcome
    }
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("Progress", "Percent complete %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Scale.dependence.Wei_S.single(a[j, 
      ], b, scale, MI)$Scale_dependence_Wei_S))
      info = sprintf("Percent complete %d%%", round(j * 
                                                      100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Scale_dependence_Wei_S")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }

  nx = (maxx - minx - indis)%/%lag
  ny = (maxy - miny - indis)%/%lag
  n = min(nx, ny)
  I = matrix(NA, n, 1)
  Z = matrix(NA, n, 1)
  distance = matrix(NA, n, 1)
  x = seq(minx, maxx, length= seq+1)
  y = seq(miny, maxy, length= seq+1)
  xy = expand.grid(x, y)
  xy$size = 0
  colnames(xy) = c("x", "y", "size")
  vaule = Scale.dependence.Wei_S.mult(xy, b, scale, MI)
  vaule = vaule[complete.cases(vaule), ]
  vaule = vaule[, c(1, 2, 3)]
  colnames(vaule) = c("x", "y", "Scale_dependence_Wei_S")
  vaule = as.data.frame(vaule)
  for (i in 1:n) {
    vaule.dists = as.matrix(dist(cbind(vaule$x, vaule$y)))
    vaule.dists = (vaule.dists - indis)%/%(i * lag) * i * lag + 0.5 * (i * lag + indis)
    vaule.dists[vaule.dists<=0]=0.5 * (i * lag + indis)
    vaule.dists.inv = 1/vaule.dists
    Moran = Moran.I(vaule$Scale_dependence_Wei_S, vaule.dists.inv)
    I[i] = Moran$observed
    Z[i] = (Moran$observed - Moran$expected)/Moran$sd
    distance[i] = lag * (i) + indis
    ISAA = cbind(distance, I, Z)
    colnames(ISAA) = c("lag", "Moran_I", "Z_Score")

  }
  vaule.dists = as.matrix(dist(cbind(vaule$x, vaule$y)))
  vaule.dists [vaule.dists <indis]=0.5*indis
  vaule.dists.inv = 1/vaule.dists
  Moran = Moran.I(vaule$Scale_dependence_Wei_S, vaule.dists.inv)
  I = Moran$observed
  Z = (Moran$observed - Moran$expected)/Moran$sd
  distance =indis
  ISAA_1 = cbind(distance, I, Z)
  colnames(ISAA_1) = c("lag", "Moran_I", "Z_Score")
  ISAA=rbind(ISAA_1,ISAA)
  
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
  Z_positive=subset(Z,Z>0)
  Zbef = matrix(NA, length(ISAA[, 3]) - 2, 1)
  Zbac = matrix(NA, length(ISAA[, 3]) - 2, 1)
  
  for (i in 1:(length(Z) - 2)) {
    Zbef[i] = Z[i + 1] - Z[i]
    Zbac[i] = Z[i + 2] - Z[i + 1]
  }
  n=length(vaule$Scale_dependence_Wei_S) 
  ZBB = cbind(Zbef, Zbac)
  ZBB = cbind(ISAA$lag[-c(1, length(Z))], ZBB)
  ZBB = ZBB[which(ZBB[, 2] > 0 & ZBB[, 3] < 0), 1]
  Z_positive=cbind(ISAA$lag[-c(1, length(Z))], Z[-c(1, length(Z))])
  ZBB_positive=intersect(ZBB, Z_positive[which(Z_positive[,2]>0),1])
  p = ggplot(ISAAIZ, aes(x = Lag)) + 
    geom_hline(aes(yintercept = 1.96/yBL), linetype = 5, col = "black",size=1) + 
    geom_hline(aes(yintercept = -1.96/yBL), linetype = 5, col = "black",size=1) + 
    geom_hline(aes(yintercept = (-1/(n-1))), linetype = 2,col = "blue",size=1) + 
    geom_vline(xintercept = ZBB_positive, linetype = 5, col = "red",size=1) + 
    geom_hline(aes(yintercept = 0), linetype = 1,col = "gray",size=1)+
    geom_line(aes(y = Moran_I, group = Index)) +
    geom_point(aes(y = Moran_I, shape = Index), size = 2)
  p = p + scale_y_continuous(sec.axis = sec_axis(~. * yBL, name = "Z_Score")) +
    theme_bw()
  p = p +
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))
  p
}
 

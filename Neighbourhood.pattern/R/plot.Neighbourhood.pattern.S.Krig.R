plot.Neighbourhood.pattern.S.Krig=function (minx, maxx, miny, maxy, b, seq, n) 
{
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Neighbourhood.pattern.S.mult = function(a, b, n) {
    Neighbourhood.pattern.S.single = function(a, b, n) {
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
      Neighbourhood1[which(Neighbourhood1[, 5] <= 0), 5] = 0
      Neighbourhood1[which(Neighbourhood1[, 4] == 0), 4] = 1e-08
      S1 = Neighbourhood1[, 5]/Neighbourhood1[, 4]
      S1 = S1/(S1 + 1)
      S2 = as.data.frame(S1)
      S = sum(S2, na.rm = T)
      S = S/n
      Neighbourhood_pattern_S = S
      outcome = list(a = a, Neighbourhood = Neighbourhood, 
                     Neighbourhood_pattern_S = Neighbourhood_pattern_S)
      outcome
    }
    d = matrix(NA, nrow(a), 4)
    pb = tkProgressBar("进度", "已完成 %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.pattern.S.single(a[j, 
      ], b, n)$Neighbourhood_pattern_S))
      info = sprintf("已完成 %d%%", round(j * 
                                         100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Size", "Neighbourhood_pattern_S")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  xgrid = seq(minx, maxx, length.out = seq + 1)
  ygrid = seq(miny, maxy, length.out = seq + 1)
  basexy = expand.grid(xgrid, ygrid)
  xgrid.cen = seq(minx + 0.5 * (maxx - minx) - 0.5/seq * (maxx - 
                                                            minx), minx + 0.5 * (maxx - minx) + 0.5/seq * (maxx - 
                                                                                                             minx), length.out = 3)
  ygrid.cen = seq(miny + 0.5 * (maxy - miny) - 0.5/seq * (maxy - 
                                                            miny), miny + 0.5 * (maxy - miny) + 0.5/seq * (maxy - 
                                                                                                             miny), length.out = 3)
  basexy.cen = expand.grid(xgrid.cen, ygrid.cen)
  xgrid.left.top = seq(minx + 0.25 * (maxx - minx) - 0.5/seq * 
                         (maxx - minx), minx + 0.25 * (maxx - minx) + 0.5/seq * 
                         (maxx - minx), length.out = 3)
  ygrid.left.top = seq(miny + 0.75 * (maxy - miny) - 0.5/seq * 
                         (maxy - miny), miny + 0.75 * (maxy - miny) + 0.5/seq * 
                         (maxy - miny), length.out = 3)
  basexy.left.top = expand.grid(xgrid.left.top, ygrid.left.top)
  xgrid.left.bottom = seq(minx + 0.25 * (maxx - minx) - 0.5/seq * 
                            (maxx - minx), minx + 0.25 * (maxx - minx) + 0.5/seq * 
                            (maxx - minx), length.out = 3)
  ygrid.left.bottom = seq(miny + 0.25 * (maxy - miny) - 0.5/seq * 
                            (maxy - miny), miny + 0.25 * (maxy - miny) + 0.5/seq * 
                            (maxy - miny), length.out = 3)
  basexy.left.bottom = expand.grid(xgrid.left.bottom, ygrid.left.bottom)
  xgrid.right.top = seq(minx + 0.75 * (maxx - minx) - 0.5/seq * 
                          (maxx - minx), minx + 0.75 * (maxx - minx) + 0.5/seq * 
                          (maxx - minx), length.out = 3)
  ygrid.right.top = seq(miny + 0.75 * (maxy - miny) - 0.5/seq * 
                          (maxy - miny), miny + 0.75 * (maxy - miny) + 0.5/seq * 
                          (maxy - miny), length.out = 3)
  basexy.right.top = expand.grid(xgrid.right.top, ygrid.right.top)
  xgrid.right.bottom = seq(minx + 0.75 * (maxx - minx) - 0.5/seq * 
                             (maxx - minx), minx + 0.75 * (maxx - minx) + 0.5/seq * 
                             (maxx - minx), length.out = 3)
  ygrid.right.bottom = seq(miny + 0.25 * (maxy - miny) - 0.5/seq * 
                             (maxy - miny), miny + 0.25 * (maxy - miny) + 0.5/seq * 
                             (maxy - miny), length.out = 3)
  basexy.right.bottom = expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  xgrid = seq(minx, maxx, length.out = 200)
  ygrid = seq(miny, maxy, length.out = 200)
  basexy1 = expand.grid(xgrid, ygrid)
  Basexy1 = basexy1
  basexy = rbind(basexy, basexy.cen, basexy.left.top, basexy.left.bottom, 
                 basexy.right.top, basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x + y
  gridded(basexy1) <- TRUE
  basexy[, 3] = 0
  colnames(basexy) <- c("x", "y", "Size")
  basexy = dplyr::distinct(basexy)
  data = Neighbourhood.pattern.S.mult(basexy, b, n)
  data1 = data
  coordinates(data) <- c("x", "y")
  spplot(data, "Neighbourhood_pattern_S")
  vgm1 <- variogram(Neighbourhood_pattern_S ~ 1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1, vgm("Sph"))
  print(m)
  p1 = plot(vgm1, model = m)
  print(p1)
  krige_res <- krige(Neighbourhood_pattern_S ~ 1, data, basexy1, 
                     model = m)
  z = krige_res$var1.pred
  Basexyz = cbind(Basexy1, z)
  colnames(Basexyz) = c("x", "y", "Vaule")
  p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y, 
                                                  fill = Vaule)) + theme_bw() + scale_fill_gradientn(colours = terrain.colors(10))
  p2 = p2 + labs(title = "The neighbourhood pattern S") + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  print(p2)
}
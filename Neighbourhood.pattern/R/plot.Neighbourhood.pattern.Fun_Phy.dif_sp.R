plot.Neighbourhood.pattern.Fun_Phy.dif_sp=function (a, b, eco_distance, n) 
{
  library(ggplot2)
  Neighbourhood.pattern.Fun_Phy.mult = function(a, b, eco_distance, 
                                                n) {
    library(tcltk)
    Neighbourhood.pattern.Fun_Phy.single = function(a, b, 
                                                    eco_distance, n) {
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
      dis = matrix(NA, nrow(Neighbourhood), 1)
      for (i in 1:nrow(Neighbourhood)) {
        dis[i, ] = unique(eco_distance[eco_distance[, 
                                                    1] == as.character(a[, 3]) & eco_distance[, 
                                                                                              2] == as.character(Neighbourhood[i, 3]) | eco_distance[, 
                                                                                                                                                     1] == as.character(Neighbourhood[i, 3]) & eco_distance[, 
                                                                                                                                                                                                            2] == as.character(a[, 3]), 3])
      }
      Neighbourhood = cbind(Neighbourhood, dis)
      colnames(Neighbourhood) = c("x", "Y", 
                                  "Species", "Fun|Phy_Distance")
      Neighbourhood_pattern_Fun.phy_Richness = sum(Neighbourhood[, 
                                                                 4]/n)
      output = list(a = a, Neighbourhood = Neighbourhood, 
                    Neighbourhood_pattern_Fun.phy_Richness = Neighbourhood_pattern_Fun.phy_Richness)
      output
    }
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("进度", "已完成 %", 
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Neighbourhood.pattern.Fun_Phy.single(a[j, 
      ], b, eco_distance, n)$Neighbourhood_pattern_Fun.phy_Richness))
      info = sprintf("已完成 %d%%", round(j * 
                                         100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Neighbourhood_pattern_Fun.phy_Richness")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  plot.data = Neighbourhood.pattern.Fun_Phy.mult(a, b, eco_distance, 
                                                 n)
  plot.data = cbind(plot.data, a$Species)
  plot.data = as.data.frame(plot.data)
  colnames(plot.data) = c("x", "y", "Neighbourhood_pattern_Fun.phy_Richness", 
                          "Species")
  p = ggplot() + geom_point(data = plot.data, aes(x, y, color = Species), 
                            size = 12/max(plot.data[, 3]) * plot.data[, 3] + 3, shape = 1)
  p1 = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw()
  p2 = p1 + xlab("") + ylab("") + labs(title = "The neighbourhood pattern Fun.phy_Richness of diffence species")
  p2
}

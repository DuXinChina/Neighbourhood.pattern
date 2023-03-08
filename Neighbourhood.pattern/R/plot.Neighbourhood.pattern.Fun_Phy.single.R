plot.Neighbourhood.pattern.Fun_Phy.single=function (a, b, eco_distance, n) 
{
  library(ggplot2)
  Neighbourhood.pattern.Fun_Phy.single = function(a, b, eco_distance, 
                                                  n) {
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
      dis[i, ] = unique(eco_distance[eco_distance[, 1] == 
                                       as.character(a[, 3]) & eco_distance[, 2] == as.character(Neighbourhood[i, 
                                                                                                              3]) | eco_distance[, 1] == as.character(Neighbourhood[i, 
                                                                                                                                                                    3]) & eco_distance[, 2] == as.character(a[, 3]), 
                                     3])
    }
    Neighbourhood = cbind(Neighbourhood, dis)
    colnames(Neighbourhood) = c("x", "Y", "Species", 
                                "Fun|Phy_Distance")
    Neighbourhood_pattern_Fun.phy_Richness = sum(Neighbourhood[, 
                                                               4]/n)
    output = list(a = a, Neighbourhood = Neighbourhood, Neighbourhood_pattern_Fun.phy_Richness = Neighbourhood_pattern_Fun.phy_Richness)
    output
  }
  Ne = Neighbourhood.pattern.Fun_Phy.single(a, b, eco_distance, 
                                            n)
  center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[, 
                                                          2], each = n), rep(0, each = n), c(1:n)))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "Fun.phy_distance", 
                       "group")
  Neigh = cbind(Ne$Neighbourhood[, c(1:2, 4)], c(1:n))
  center1 = cbind(cbind(c(-100, -100), c(-100, -100), c(1:0), 
                        c(1:2)))
  center1 = as.data.frame(center1)
  colnames(center1) = c("x", "y", "Fun.phy_distance", 
                        "group")
  colnames(Neigh) = c("x", "y", "Fun.phy_distance", 
                      "group")
  Neigh1 = rbind(Neigh, center1)
  unitg = rbind(center, Neigh)
  options(warn = -1)
  max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[, 
                                                                         2] - Ne$a[, 2]))
  p = ggplot() + geom_line(data = unitg, aes(x = x, y = y, 
                                             group = group), linetype = 2, size = 1) + geom_point(data = Ne$a, 
                                                                                                  aes(x, y), size = 10, color = "red4", alpha = 0.7) + 
    geom_point(data = Neigh1, aes(x, y, color = Fun.phy_distance), 
               size = 10, alpha = 0.7)
  p = p + scale_colour_gradient(low = "red4", high = "green4")
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Neigh, 
                                                                    aes(x, y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - max, Ne$a[, 1] + max), y = c(Ne$a[, 
                                                                   2] - max, Ne$a[, 2] + max))
  p = p + annotate("text", x = Ne$a[, 1] - max + 0.65 * 
                     (max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("Neighbourhood_pattern_Fun.phy_Richness=", 
                                                                                round(Ne$Neighbourhood_pattern_Fun.phy_Richness, 3)))
  p = p + theme_bw()
  print(p)
  options(warn = 1)
}
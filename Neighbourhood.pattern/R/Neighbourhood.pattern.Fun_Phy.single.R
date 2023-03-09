Neighbourhood.pattern.Fun_Phy.single=function (a, b, eco_distance, n) 
{
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
    dis[i, ] = unique(eco_distance[eco_distance[, 1] == as.character(a[, 
                                                                       3]) & eco_distance[, 2] == as.character(Neighbourhood[i, 
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
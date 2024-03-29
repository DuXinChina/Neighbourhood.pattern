Neighbourhood.pattern.S.single=function (a, b, n) 
{
  Neighbourhood.single1 = function(a, b, n) {
    c = b[, 1:2]
    for (i in 1:nrow(b)) {
      c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
      d = (c[, 1] + c[, 2])^(1/2)
    }
    d = cbind(b, d)
    d = subset(d, d > 0)
    d = d[order(d[, 4])[1:n], 1:4]
    colnames(d) = c("x", "Y", "Size", "Distance")
    d
  }
  Nei.tree = Neighbourhood.single1(a, b, n)
  Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - a[, 3])
  colnames(Nei.dif.high) = c("x", "Y", "Size", 
                             "Distance", "Size_dif ")
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
  outcome = list(a = a, Neighbourhood = Neighbourhood, Neighbourhood_pattern_S = Neighbourhood_pattern_S)
  outcome
}
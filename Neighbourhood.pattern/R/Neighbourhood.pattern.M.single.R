Neighbourhood.pattern.M.single=function (a, b, n) 
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
    colnames(d) = c("x", "y", "Species")
    d
  }
  Neighbourhood = Neighbourhood.single(a, b, n)
  Neighbourhood1 = subset(Neighbourhood, Species != a[, 3])
  M = nrow(Neighbourhood1)/n
  output = list(a = a, Neighbourhood = Neighbourhood, M = M)
  output
}
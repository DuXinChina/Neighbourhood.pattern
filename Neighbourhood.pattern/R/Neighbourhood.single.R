Neighbourhood.single=function (a, b, n) 
{
  c = b
  for (i in 1:nrow(b)) {
    c[i, ] = (b[i, ] - a[1, ])^2
    d = (c[, 1] + c[, 2])^(1/2)
  }
  d = cbind(b, d)
  d = subset(d, d > 0)
  d = d[order(d[, 3])[1:n], 1:3]
  colnames(d) = c("x", "Y", "Distance")
  Neighbourhood = d
  a = a
  out = list(a = a, Neighbourhood = Neighbourhood)
  out
}
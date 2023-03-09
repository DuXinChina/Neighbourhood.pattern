Neighbourhood.pattern.C.single=function (a, b, n) 
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
    colnames(d) = c("x", "y", "distance", 
                    "size")
    d
  }
  Neighbourhood = Neighbourhood.single(a, b, n)
  C = subset(Neighbourhood, distance < (Neighbourhood$size + 
                                          a$size)/2)
  C = nrow(C)/n
  outcome = list(a = a, Neighbourhood = Neighbourhood, C = C)
  outcome
}
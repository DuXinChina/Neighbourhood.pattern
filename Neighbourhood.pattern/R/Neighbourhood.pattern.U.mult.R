Neighbourhood.pattern.U.mult=function (a, b, n) 
{
  library(tcltk)
  Neighbourhood.pattern.U.single = function(a, b, n) {
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
      colnames(d) = c("x", "y", "Size")
      d
    }
    Neighbourhood = Neighbourhood.single(a, b, n)
    U = subset(Neighbourhood, Size > a[, 3])
    U = nrow(U)/n
    outcome = list(a = a, Neighbourhood = Neighbourhood, 
                   U = U)
    outcome
  }
  d = matrix(NA, nrow(a), 4)
  pb = tkProgressBar("进度", "已完成 %", 0, 
                     100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.pattern.U.single(a[j, 
    ], b, n)$U))
    info = sprintf("已完成 %d%%", round(j * 100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "size", "U")
  rownames(d) = 1:nrow(a)
  d = as.data.frame(d)
  d
}

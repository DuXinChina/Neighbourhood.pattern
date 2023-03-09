Neighbourhood.mult=function (a, b, n) 
{
  library(tcltk)
  Neighbourhood.single = function(a, b, n) {
    c = b
    for (i in 1:nrow(b)) {
      c[i, ] = (b[i, ] - a[1, ])^2
      d = (c[, 1] + c[, 2])^(1/2)
    }
    d = cbind(b, d)
    d = subset(d, d > 0)
    d = d[order(d[, 3])[1:n], 1:3]
    colnames(d) = c("x", "Y", "Distance")
    rownames(d) = 1:n
    Neighbourhood = d
    a = a
    out = list(a = a, Neighbourhood = Neighbourhood)
    out
  }
  e = matrix(NA, nrow(a), 3 * n + 2)
  pb = tkProgressBar("进度", "已完成 %", 0, 
                     100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d = cbind(a[j, ], Neighbourhood.single(a[j, ], b, n)$Neighbourhood[1, 
    ])
    for (i in 2:n) {
      d = cbind(d, Neighbourhood.single(a[j, ], b, n)$Neighbourhood[i, 
      ])
      d = as.matrix(d)
      d
    }
    e[j, ] = d
    info = sprintf("已完成 %d%%", round(j * 100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)", 
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  rownames(e) = 1:j
  colnames(e) = c("X", "Y", rep(c("X", "Y", 
                                  "Distance"), n))
  e
}
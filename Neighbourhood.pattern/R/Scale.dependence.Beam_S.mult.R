
Scale.dependence.Beam_S.mult=function(a,b,scale,mi,MI)
{
library(tcltk)
Scale.dependence.Beam_S.single=function(a,b,scale,mi,MI)
{
  Neighbourhood.single1 = function(a, b, scale) {
    c = b[, 1:2]
    for (i in 1:nrow(b)) {
      c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
      d = (c[, 1] + c[, 2])^(1/2)
    }
    d = cbind(b, d)
    d = subset(d, d < scale)
    colnames(d) = c("x", "Y", "Size", 
                    "Distance")
    d
  }
  Nei.tree = Neighbourhood.single1(a, b, scale)
  n = nrow(Nei.tree)
  Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - a[, 3])
  colnames(Nei.dif.high) = c("x", "Y",  "Size", "Distance", "Size_dif ")
  n.dif = subset(Nei.dif.high, Nei.dif.high$Size_dif >  0)
  n.dif = nrow(n.dif)
  if (n.dif == 0) {
    Scale_dependence_S_Beam = 0
    Neighbourhood = Nei.dif.high
  }
  if (n.dif != 0) {
    Neighbourhood = Nei.dif.high
    a=a[rep(1,nrow(Neighbourhood)),]
    Weight=Neighbourhood[,1:2]-a[,1:2]
    Weight[,3]=atan(Weight[,2]/abs(Weight[,1]))
    Weight[,3]=Weight[,3]
    Weight[,3]=(pi/2)+Weight[,3]
    ###上面weight的第三列是邻近树木到a的方向与正南方向的夹�?
    Weight[,4]=((pi)-mi*Weight[,3])/(pi)
    Weight[which(Weight[,4]<0),4]=0
    Neighbourhood1 = Neighbourhood
    Neighbourhood[which(Neighbourhood[, 5] <= 0),  5] = NA
    Neighbourhood1[which(Neighbourhood1[, 5] <= 0),  5] = 0
    Neighbourhood1[which(Neighbourhood1[, 4] == 0),  4] = 1e-08
    S1 = Neighbourhood1[, 5]/Neighbourhood1[, 4]*Weight[,4]
    S2 = as.data.frame(S1)
    S = sum(S2, na.rm = T)
    Scale_dependence_S_Beam = S
  }
  Scale_dependence_S_Beam = (atan(Scale_dependence_S_Beam/(MI * pi))/pi * 2)
  outcome = Scale_dependence_S_Beam
  outcome
}
d = matrix(NA, nrow(a), 3)  
pb = tkProgressBar("", "Percent complete %", 0, 100)
star_time = Sys.time()
for (j in 1:nrow(a)) {
  d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Scale.dependence.Beam_S.single(a[j,],b,scale,mi,MI)))
  info = sprintf("Percent complete %d%%", round(j * 100/nrow(a)))
  setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)", 
                                                info), info)
}
end_time = Sys.time()
close(pb)
run_time = end_time - star_time
colnames(d) = c("x", "y", "Scale_dependence_Beam_S")
rownames(d) = 1:nrow(a)
d = as.data.frame(d)
d
}
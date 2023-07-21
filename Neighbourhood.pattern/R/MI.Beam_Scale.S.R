
MI.Beam_Scale.S=function (minx, maxx, miny, maxy, b, seq, scale, mi) 
{
library(sp)
library(gstat)
library(tcltk)
xgrid = seq(minx, maxx, length.out = seq + 1)
ygrid = seq(miny, maxy, length.out = seq + 1)
basexy = expand.grid(xgrid, ygrid)
basexy[, 3] = 0
colnames(basexy) <- c("x", "y", "Size")
a = dplyr::distinct(basexy)
Scale.dependence.Beam_S.mult=function(a,b,scale,mi)
{
  library(tcltk)
  Scale.dependence.Beam_S.single=function(a,b,scale,mi)
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
    Scale_dependence_S_Beam = Scale_dependence_S_Beam
    outcome = Scale_dependence_S_Beam
    outcome
  }
  d = matrix(NA, nrow(a), 3)  
  pb = tkProgressBar("", "Percent complete %", 0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Scale.dependence.Beam_S.single(a[j,],b,scale,mi)))
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
data = Scale.dependence.Beam_S.mult(a,b,scale,mi)
Scale.dependence.Beam_S=data$Scale_dependence_Beam_S
Rank = rank(Scale.dependence.Beam_S)
Percentile_Rank_of_Scale.dependence.Beam_S = Rank/length(Rank)
nls = nls(Percentile_Rank_of_Scale.dependence.Beam_S ~ (atan(Scale.dependence.Beam_S/(MI * 
                                                                              pi))/pi * 2), start = list(MI = 1), algorithm = "port")
res = summary(nls)
MI = res$parameters[1, 1]
Revise_Scale.dependence.Beam_S = (atan(Scale.dependence.Beam_S/(MI * 
                                                        pi))/pi * 2)
plot(Percentile_Rank_of_Scale.dependence.Beam_S,Revise_Scale.dependence.Beam_S,  
     main = "Q-Q plot", xlim = c(0, 1), ylim = c(0, 
                                                 1)) + abline(0, 1, col = "red4")
R_square = 1 - sum((Percentile_Rank_of_Scale.dependence.Beam_S - 
                      Revise_Scale.dependence.Beam_S)^2)/sum((Percentile_Rank_of_Scale.dependence.Beam_S - mean(Percentile_Rank_of_Scale.dependence.Beam_S))^2)
output = list(MI = MI, R_square = R_square)
output
}

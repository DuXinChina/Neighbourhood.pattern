

MI.Beam_Tan.S=function (minx, maxx, miny, maxy, b, seq, tan, mi) 
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

Tan.dependence.Beam_S.mult=function (a, b, tan, mi) 
{
  library(tcltk)
  Tan.dependence.Beam_S.single=function (a, b, tan, mi) 
  {
    Neighbourhood.single1 = function(a, b, tan) {
      c = b[, 1:2]
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      d = cbind(b, d)
      d = subset(d, d > 0)
      d$tangent = (d$size - a[, 3])/d$d
      d = subset(d, tangent > tan)
      colnames(d) = c("x", "Y", "Size", 
                      "Distance", "tangent")
      d
    }
    Nei.tree = Neighbourhood.single1(a, b, tan)
    n.dif = nrow(Nei.tree)
    if (n.dif == 0) {
      Tan_dependence_Beam_S = 0
      Neighbourhood = Nei.tree
    }
    if (n.dif != 0) {
      Neighbourhood = Nei.tree
      Neighbourhood1 = Neighbourhood
      Neighbourhood[which(Neighbourhood[, 5] <= 0), 5] = NA
      Neighbourhood1[which(Neighbourhood1[, 5] <= 0), 5] = 0
      a=a[rep(1,nrow(Neighbourhood)),]
      Weight=Neighbourhood[,1:2]-a[,1:2]
      Weight[,3]=atan(Weight[,2]/abs(Weight[,1]))
      Weight[,3]=Weight[,3]
      Weight[,3]=(pi/2)+Weight[,3]
      Weight[,4]=((pi)-mi*Weight[,3])/(pi)
      Weight[which(Weight[,4]<0),4]=0
      S1 = Neighbourhood1[, 5]*Weight[,4]
      S2 = as.data.frame(S1)
      S = sum(S2, na.rm = T)
      Tan_dependence_Beam_S = S
    }
    Tan_dependence_Beam_S = Tan_dependence_Beam_S
    if (is.nan(Tan_dependence_Beam_S) == T) 
      (Tan_dependence_Beam_S = 1)
    outcome =  Tan_dependence_Beam_S
    outcome
  }
  
  d = matrix(NA, nrow(a), 3)  
  pb = tkProgressBar("", "Percent complete %", 0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Beam_S.single(a[j,], b, tan, mi)))
    info = sprintf("Percent complete %d%%", round(j * 100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)", 
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "Tan_dependence_Beam_S")
  rownames(d) = 1:nrow(a)
  d = as.data.frame(d)
  d
}
data = Tan.dependence.Beam_S.mult(a, b, tan, mi) 
Tan_dependence_Beam_S=data$Tan_dependence_Beam_S
Rank = rank(Tan_dependence_Beam_S)
Percentile_Rank_of_Tan_dependence_Beam_S = Rank/length(Rank)
nls = nls(Percentile_Rank_of_Tan_dependence_Beam_S ~ (atan(Tan_dependence_Beam_S/(MI * 
                                                                                        pi))/pi * 2), start = list(MI = 1), algorithm = "port")
res = summary(nls)
MI = res$parameters[1, 1]
Revise_Tan_dependence_Beam_S = (atan(Tan_dependence_Beam_S/(MI * 
                                                                  pi))/pi * 2)
plot(Percentile_Rank_of_Tan_dependence_Beam_S,Revise_Tan_dependence_Beam_S,  
     main = "Q-Q plot", xlim = c(0, 1), ylim = c(0, 
                                                 1)) + abline(0, 1, col = "red4")
R_square = 1 - sum((Percentile_Rank_of_Tan_dependence_Beam_S - 
                      Revise_Tan_dependence_Beam_S)^2)/sum((Percentile_Rank_of_Tan_dependence_Beam_S - mean(Percentile_Rank_of_Tan_dependence_Beam_S))^2)
output = list(MI = MI, R_square = R_square)
output
}

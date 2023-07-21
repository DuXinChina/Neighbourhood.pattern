a=data.frame(x=50,y=50,size=0)
Tan.dependence.Beam_S.single=function (a, b, tan, mi, MI) 
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
    Tan_dependence_Beam_S = (atan(Tan_dependence_Beam_S/(MI * pi))/pi * 
                          2)
    if (is.nan(Tan_dependence_Beam_S) == T) 
      (Tan_dependence_Beam_S = 1)
    outcome =  Tan_dependence_Beam_S
    outcome
  }

Tan.dependence.Beam_S.single(a[1,], b, tan=3, mi=0.8, MI=5) 

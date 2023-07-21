plot.Scale.dependence.Beam_S.Krig=function(minx, maxx, miny, maxy, b, seq, scale, mi, MI)
{
library(sp)
library(gstat)
library(tcltk)
library(ggplot2)
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
xgrid = seq(minx, maxx, length.out = seq + 1)
ygrid = seq(miny, maxy, length.out = seq + 1)
basexy = expand.grid(xgrid, ygrid)
xgrid.cen = seq(minx + 0.5 * (maxx - minx) - 0.5/seq * (maxx - 
                                                          minx), minx + 0.5 * (maxx - minx) + 0.5/seq * (maxx - 
                                                                                                           minx), length.out = 3)
ygrid.cen = seq(miny + 0.5 * (maxy - miny) - 0.5/seq * (maxy - 
                                                          miny), miny + 0.5 * (maxy - miny) + 0.5/seq * (maxy - 
                                                                                                           miny), length.out = 3)
basexy.cen = expand.grid(xgrid.cen, ygrid.cen)
xgrid.left.top = seq(minx + 0.25 * (maxx - minx) - 0.5/seq * 
                       (maxx - minx), minx + 0.25 * (maxx - minx) + 0.5/seq * 
                       (maxx - minx), length.out = 3)
ygrid.left.top = seq(miny + 0.75 * (maxy - miny) - 0.5/seq * 
                       (maxy - miny), miny + 0.75 * (maxy - miny) + 0.5/seq * 
                       (maxy - miny), length.out = 3)
basexy.left.top = expand.grid(xgrid.left.top, ygrid.left.top)
xgrid.left.bottom = seq(minx + 0.25 * (maxx - minx) - 0.5/seq * 
                          (maxx - minx), minx + 0.25 * (maxx - minx) + 0.5/seq * 
                          (maxx - minx), length.out = 3)
ygrid.left.bottom = seq(miny + 0.25 * (maxy - miny) - 0.5/seq * 
                          (maxy - miny), miny + 0.25 * (maxy - miny) + 0.5/seq * 
                          (maxy - miny), length.out = 3)
basexy.left.bottom = expand.grid(xgrid.left.bottom, ygrid.left.bottom)
xgrid.right.top = seq(minx + 0.75 * (maxx - minx) - 0.5/seq * 
                        (maxx - minx), minx + 0.75 * (maxx - minx) + 0.5/seq * 
                        (maxx - minx), length.out = 3)
ygrid.right.top = seq(miny + 0.75 * (maxy - miny) - 0.5/seq * 
                        (maxy - miny), miny + 0.75 * (maxy - miny) + 0.5/seq * 
                        (maxy - miny), length.out = 3)
basexy.right.top = expand.grid(xgrid.right.top, ygrid.right.top)
xgrid.right.bottom = seq(minx + 0.75 * (maxx - minx) - 0.5/seq * 
                           (maxx - minx), minx + 0.75 * (maxx - minx) + 0.5/seq * 
                           (maxx - minx), length.out = 3)
ygrid.right.bottom = seq(miny + 0.25 * (maxy - miny) - 0.5/seq * 
                           (maxy - miny), miny + 0.25 * (maxy - miny) + 0.5/seq * 
                           (maxy - miny), length.out = 3)
basexy.right.bottom = expand.grid(xgrid.right.bottom, ygrid.right.bottom)
xgrid = seq(minx, maxx, length.out = 200)
ygrid = seq(miny, maxy, length.out = 200)
basexy1 = expand.grid(xgrid, ygrid)
Basexy1 = basexy1
basexy = rbind(basexy, basexy.cen, basexy.left.top, basexy.left.bottom, 
               basexy.right.top, basexy.right.bottom)
colnames(basexy1) <- c("x", "y")
coordinates(basexy1) <- ~x + y
gridded(basexy1) <- TRUE
basexy[, 3] = 0
colnames(basexy) <- c("x", "y", "Size")
basexy = dplyr::distinct(basexy)
data = Scale.dependence.Beam_S.mult(basexy,b,scale,mi,MI)
data1 = data
coordinates(data) <- c("x", "y")
spplot(data, "Scale_dependence_Beam_S")
vgm1 <- variogram(Scale_dependence_Beam_S ~ 1, data)
plot(vgm1, plot.numbers = TRUE)
m <- fit.variogram(vgm1, vgm("Sph"))
print(m)
sd = subset(vgm1$dist, vgm1$dist < m$range[2])
spre = m$psill[1] + (m$psill[2]) * ((3 * sd)/(2 * m$range[2]) - 
                                      sd^3/(2 * (m$range[2])^3))
bd = subset(vgm1$dist, vgm1$dist > m$range[2])
bpre = rep((m$psill[1] + m$psill[2]), length(bd))
pre = rbind(as.matrix(spre), as.matrix(bpre))
Coefficient_of_Determination = 1 - sum((pre - vgm1$gamma)^2)/sum((vgm1$gamma - 
                                                                    mean(vgm1$gamma))^2)
p1 = plot(vgm1, model = m)
print(p1)
krige_res <- krige(Scale_dependence_Beam_S ~ 1, data, basexy1, 
                   model = m)
z = krige_res$var1.pred
Basexyz = cbind(Basexy1, z)
colnames(Basexyz) = c("x", "y", "Vaule")
Basexyz[which(Basexyz$Vaule < 0), 3] = 0
p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y, fill = Vaule)) + theme_bw() + scale_fill_gradientn(limits = c(min(Basexyz$Vaule),max(Basexyz$Vaule)), colours = terrain.colors(10))
p2 = p2 + labs(title = "Scale_dependence_Beam_S") + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
print(p2)
}

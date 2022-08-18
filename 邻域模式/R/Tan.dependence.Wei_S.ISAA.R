
Tan.dependence.Wei_S.ISAA=function(minx,maxx,miny,maxy,b,indis,lag,tan,MI)
{
library(ape)
library(sp)
library(gstat)
library(tcltk)  
library(ggplot2)
Tan.dependence.Wei_S.mult=function(a,b,tan)
{
  library(tcltk)  
  Tan.dependence.Wei_S.single=function(a,b,tan)
  {
    
    Tan.dependence.S.single=function(a,b,tan)
    {
      Neighbourhood.single1=function(a,b,tan)###找出一定尺度内的邻体         
      {#3
        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4 
        d=cbind(b,d)
        #d=subset(d,d>0)
        d$tangent=(d$size-a[,3])/d$d
        d=subset(d,tangent>tan)
        colnames(d) = c("x","Y","Size","Distance","tangent")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,tan)
      n.dif=nrow(Nei.tree)
      if(n.dif==0)
      {Tan_dependence_S=0
      Neighbourhood=Nei.tree}
      if(n.dif!=0)
      {
        Neighbourhood=Nei.tree
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
        Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
        S1=Neighbourhood1[,5]/Neighbourhood1[,4]
        S2=as.data.frame(S1)
        S=sum(S2,na.rm=T)
        Tan_dependence_S=S
      }
      Tan_dependence_S=(atan(Tan_dependence_S/(MI*pi))/pi*2)
      if(is.nan(Tan_dependence_S)==T)(Tan_dependence_S=1)
      outcome=list(a=a,Neighbourhood=Neighbourhood,Tan_dependence_S=Tan_dependence_S)
      outcome
    }
    
    Tan.dependence.Simpson=function(a,b,tan)
    {
      Neighbourhood.single1=function(a,b,tan)###找出一定尺度内的邻体         
      {#3
        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4 
        d=cbind(b,d)
        #d=subset(d,d>0)
        d$tangent=(d$size-a[,3])/d$d
        d=subset(d,tangent>tan)
        colnames(d) = c("x","y","Size","Distance","tangent")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,tan)
      Nei.tree[,3]=Nei.tree[,3]-a[,3]
      Nei.tree=subset(Nei.tree,Size>0)
      Nei.tree
      Simpn=sum(Nei.tree[,3])
      Nei.tree=subset(Nei.tree,Distance>0)
      Nei.tree1=subset(Nei.tree,x>=a[,1]&y>a[,2])####将Nei.tree划分到四个象限
      Nei.tree2=subset(Nei.tree,x>a[,1]&y<=a[,2])
      Nei.tree3=subset(Nei.tree,x<=a[,1]&y<a[,2])
      Nei.tree4=subset(Nei.tree,x<a[,1]&y>=a[,2])
      Simpn1=sum(Nei.tree1[,3])
      Simpn2=sum(Nei.tree2[,3])
      Simpn3=sum(Nei.tree3[,3])
      Simpn4=sum(Nei.tree4[,3])
      simp_wei=sum((Simpn1/Simpn)^2,(Simpn2/Simpn)^2,(Simpn3/Simpn)^2,(Simpn4/Simpn)^2)
      if(is.nan(simp_wei)==T)(simp_wei=0.25)
      simp_wei
    }
    #####计算尺度依赖加权荫蔽度
    S=Tan.dependence.S.single(a,b,tan)
    Levins_Simpson=Tan.dependence.Simpson(a,b,tan)
    a=S$a
    Neighbourhood=S$Neighbourhood
    Tan_dependence_S=S$Tan_dependence_S
    Tan_dependence_Wei_S=Tan_dependence_S*(0.25/Levins_Simpson)
    #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
    outcome=list(a=a,Neighbourhood=Neighbourhood,Tan_dependence_S=Tan_dependence_S,Levins_Simpson=Levins_Simpson,Tan_dependence_Wei_S=Tan_dependence_Wei_S)
    outcome
  }
  d=matrix(NA,nrow(a),3)
  pb=tkProgressBar("Progress","Percent complete %", 0, 100) 
  star_time=Sys.time()
  ## 记录程序开始时间 
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Tan.dependence.Wei_S.single(a[j,],b,tan)$Tan_dependence_Wei_S))
    info=sprintf("Percent complete %d%%", round(j*100/nrow(a))) 
    ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("Progress (%s)", info),info) 
    ## 设置进度条         
  }
  end_time =Sys.time()  
  ## 记录程序结束时间
  close(pb)  
  run_time = end_time - star_time 
  ## 计算程序运行时间
  colnames(d)=c("x","y","Tan_dependence_Wei_S")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
#计算均匀布点的合理分割分数
nx=(maxx-minx-indis)%/%lag
ny=(maxy-miny-indis)%/%lag
n=min(nx,ny)


I=matrix(NA,n,1)
Z=matrix(NA,n,1)
distance=matrix(NA,n,1)
####生成均匀布点

x=seq(minx,maxx,by=indis)
y=seq(miny,maxy,by=indis)
xy=expand.grid(x,y)
xy$size=0
colnames(xy)=c("x","y","size")
vaule=Tan.dependence.Wei_S.mult(xy,b,tan)

for(i in 1:n)
{
  vaule.dists=as.matrix(dist(cbind(vaule$x,vaule$y)))
  vaule.dists=(vaule.dists-indis)%/%(i*lag)*i*lag+0.5*(i*lag+indis)
  vaule.dists.inv = 1/vaule.dists
  diag(vaule.dists.inv)=0
  Moran=Moran.I(vaule$Tan_dependence_Wei_S,vaule.dists.inv)
  I[i]=Moran$observed
  Z[i]=(Moran$observed-Moran$expected)/Moran$sd
  distance[i]=lag*(i)+indis
  ISAA=cbind(distance,I,Z)
  colnames(ISAA)=c("lag","Moran_I","Z_Score")
}

yBL=max(ISAA[,3])/max(ISAA[,2])
yBL=1.1*yBL
ISAA=as.data.frame(ISAA)
ISAAI=ISAA[,c(1,2)]
ISAAI$Index=c(" Moran_I")
colnames(ISAAI)=c("Lag","Moran_I","Index")
ISAAZ=ISAA[,c(1,3)]
ISAAZ$Index=c("Z_Score")
colnames(ISAAZ)=c("Lag","Moran_I","Index")
ISAAZ$Moran_I=ISAAZ$Moran_I/yBL
ISAAIZ=rbind(ISAAI,ISAAZ)

Z=ISAA[,3]
Zbef=matrix(NA,length(ISAA[,3])-2,1)
Zbac=matrix(NA,length(ISAA[,3])-2,1)
for(i in 1:(length(Z)-2))
{
  Zbef[i]=Z[i+1]-Z[i]
  Zbac[i]=Z[i+2]-Z[i+1]
}
ZBB=cbind(Zbef,Zbac)
ZBB=cbind(ISAA$lag[-c(1,length(Z))],ZBB)
ZBB=ZBB[which(ZBB[,2]>0 & ZBB[,3]<0),1]

p=ggplot(ISAAIZ,aes(x=Lag)) + geom_line(aes(y=Moran_I,group=Index)) + geom_point(aes(y=Moran_I,shape=Index),size=2)
p=p+scale_y_continuous(sec.axis = sec_axis(~.*yBL, name = "Z_Score"))+theme_bw()
p=p+geom_hline(aes(yintercept=1.96/yBL),linetype=5,col="black")+
  geom_hline(aes(yintercept=-1.96/yBL),linetype=5,col="black")+
  geom_hline(aes(yintercept=0/yBL),linetype=1)+
  geom_vline(xintercept=ZBB,linetype=5,col="red")+
  theme(legend.title=element_blank(),legend.position = c(0.8,0.8))
p
}
Tan.dependence.Wei_S.ISAA(minx=5,maxx=45,miny=5,maxy=45,b=b,indis=1,lag=1,tan=4,MI=2.5)

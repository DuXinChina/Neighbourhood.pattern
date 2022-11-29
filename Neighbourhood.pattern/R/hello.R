####边界效应function~
boundary.eff=function(a,minx,maxx,miny,maxy,frame)
{
  a=subset(a,x>minx+frame&x<maxx-frame)
  a=subset(a,y>miny+frame&y<maxy-frame)
  a
}
####全林分开敞度function~
K.mult=function(a,b)####计算林分的开敞度
{
  library(tcltk)
  Neighbourhood.K.single=function(a,b)
  {
    b1=subset(b,x>=a[,1]&y>a[,2]&size>a[,3])####将b划分成四个象限
    b2=subset(b,x>a[,1]&y<=a[,2]&size>a[,3])
    b3=subset(b,x<=a[,1]&y<a[,2]&size>a[,3])
    b4=subset(b,x<a[,1]&y>=a[,2]&size>a[,3])
    Neighbourhood.K.single1=function(a,b)####计算单株林木a的开敞度
    {
      b1=subset(b,x>=a[,1]&y>a[,2]&size>a[,3])####将b划分成四个象限
      b2=subset(b,x>a[,1]&y<=a[,2]&size>a[,3])
      b3=subset(b,x<=a[,1]&y<a[,2]&size>a[,3])
      b4=subset(b,x<a[,1]&y>=a[,2]&size>a[,3])
      Neighbourhood.single=function(a,b)###找出a1林木的临近个体
      {#3
        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1],1:4]

        colnames(d) = c("x","Y","Size","Distance")
        d
      }#3
      b1.Nei.tree=Neighbourhood.single(a,b1)###分别找出四个象限中最近的树
      b2.Nei.tree=Neighbourhood.single(a,b2)
      b3.Nei.tree=Neighbourhood.single(a,b3)
      b4.Nei.tree=Neighbourhood.single(a,b4)
      Nei.tree=rbind(b1.Nei.tree,b2.Nei.tree,b3.Nei.tree,b4.Nei.tree)####葫芦娃合体
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      k1=Neighbourhood[,4]/Neighbourhood[,5]
      k2=as.data.frame(k1)
      k2[,1][is.infinite(k2[,1])] = NA
      K=sum(k2,na.rm=T)
      key=2
      outcome=list(a=a,Neighbourhood=Neighbourhood,K=K,key=key)
      outcome
    }

    if(nrow(b1)!=0&nrow(b2)!=0&nrow(b3)!=0&nrow(b4)!=0) {Neighbourhood.K.single1(a,b) }else{key1=1;key=1;out=list(key=key,key1=key1);out}
  }



  d=matrix(NA,nrow(a),4)
  pb <- tkProgressBar("进度","已完成 %", 0, 100)
  star_time <- Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    KEY=Neighbourhood.K.single(a[j,],b)
    sw=KEY$key
    if (sw==1){warning(paste("检查第",c(j),"行"))} else {d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.K.single(a[j,],b)$K))}
    info <- sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time <- Sys.time()  ## 记录程序结束时间

  ## 第三个位置关闭进度条
  close(pb)

  run_time <- end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("X","Y","Size","K")
  rownames(d)=1:nrow(a)
  if (is.na(sum(d))==TRUE) {print("请确保以参考林木为中心的四个象限内均有高于其的立木");d} else {d}

}
####单一林木的开敞度function~

K.single=function(a,b)
{
  b1=subset(b,x>=a[,1]&y>a[,2]&size>a[,3])####将b划分成四个象限
  b2=subset(b,x>a[,1]&y<=a[,2]&size>a[,3])
  b3=subset(b,x<=a[,1]&y<a[,2]&size>a[,3])
  b4=subset(b,x<a[,1]&y>=a[,2]&size>a[,3])
  Neighbourhood.K.single1=function(a,b)####计算单株林木a的开敞度
  {
    b1=subset(b,x>=a[,1]&y>a[,2]&size>a[,3])####将b划分成四个象限
    b2=subset(b,x>a[,1]&y<=a[,2]&size>a[,3])
    b3=subset(b,x<=a[,1]&y<a[,2]&size>a[,3])
    b4=subset(b,x<a[,1]&y>=a[,2]&size>a[,3])
    Neighbourhood.single=function(a,b)###找出a1林木的临近个体
    {#3
      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1],1:4]

      colnames(d) = c("x","Y","Size","Distance")
      d
    }#3
    b1.Nei.tree=Neighbourhood.single(a,b1)###分别找出四个象限中最近的树
    b2.Nei.tree=Neighbourhood.single(a,b2)
    b3.Nei.tree=Neighbourhood.single(a,b3)
    b4.Nei.tree=Neighbourhood.single(a,b4)
    Nei.tree=rbind(b1.Nei.tree,b2.Nei.tree,b3.Nei.tree,b4.Nei.tree)####葫芦娃合体
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    k1=Neighbourhood[,4]/Neighbourhood[,5]
    k2=as.data.frame(k1)
    k2[,1][is.infinite(k2[,1])] = NA
    K=sum(k2,na.rm=T)
    outcome=list(a=a,Neighbourhood=Neighbourhood,K=K)
    outcome
  }

  if(nrow(b1)!=0&nrow(b2)!=0&nrow(b3)!=0&nrow(b4)!=0) {Neighbourhood.K.single1(a,b) }else{warning("请确保以参考林木为中心的四个象限内均有高于其的立木")}
}
####寻找全林分的相邻木function~
Neighbourhood.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.single=function(a,b,n)
  {
    c=b
    for (i in 1:nrow(b))
    {
      c[i,]=(b[i,]-a[1,])^2
      d=(c[,1]+c[,2])^(1/2)
    }

    d=cbind(b,d)
    d=subset(d,d>0)
    d=d[order(d[,3])[1:n],1:3]

    colnames(d) = c("x","Y","Distance")
    rownames(d)=1:n
    Neighbourhood=d
    a=a
    out=list(a=a,Neighbourhood=Neighbourhood)
    out
  }

  e=matrix(NA,nrow(a),3*n+2)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d=cbind(a[j,],Neighbourhood.single(a[j,],b,n)$Neighbourhood[1,])

    for(i in 2:n)
    {
      d=cbind(d,Neighbourhood.single(a[j,],b,n)$Neighbourhood[i,])
      d=as.matrix(d)
      d
    }
    e[j,]=d
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  rownames(e)=1:j
  colnames(e) = c("X","Y",rep(c("X","Y","Distance"),n))
  e
}
####寻找单一参照树的相邻木function~
Neighbourhood.single=function(a,b,n)
{
  c=b
  for (i in 1:nrow(b))
  {
    c[i,]=(b[i,]-a[1,])^2
    d=(c[,1]+c[,2])^(1/2)
  }

  d=cbind(b,d)
  d=subset(d,d>0)
  d=d[order(d[,3])[1:n],1:3]

  colnames(d) = c("x","Y","Distance")
  #rownames(d)=1:n
  Neighbourhood=d
  a=a
  out=list(a=a,Neighbourhood=Neighbourhood)
  out
}
####计算全林分密集度function~

Neighbourhood.pattern.C.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.C.single=function(a,b,n)
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      colnames(d) = c("x","y","distance","size")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
    C=nrow(C)/n
    outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
    outcome
  }#2
  d=matrix(NA,nrow(a),3)
  pb = tkProgressBar("进度","已完成 %", 0, 100)
  star_time = Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.C.single(a[j,],b,n)$C))
    info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time = Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time= end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","C")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照木密集度function~
Neighbourhood.pattern.C.single=function(a,b,n)
{#2
  Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
  {#3
    a1=a
    b1=b
    a=a[,1:2]
    b=b[,1:2]
    c=b
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,]-a[1,])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d,b1[,3])
    d=subset(d,d>0)
    d=d[order(d[,3])[1:n],1:4]
    colnames(d) = c("x","y","distance","size")
    #rownames(d)=1:n
    d
  }#3
  Neighbourhood=Neighbourhood.single(a,b,n)
  C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
  C=nrow(C)/n
  outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
  outcome
}#2
####计算全林分的邻体模式系统发育多样性function~

Neighbourhood.pattern.Fun_Phy.mult=function(a,b,eco_distance,n)
{
  library(tcltk)
  Neighbourhood.pattern.Fun_Phy.single=function(a,b,eco_distance,n)
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      d=d[,c(1:2,4)]
      colnames(d) = c("x","Y","Species")
      rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    dis=matrix(NA,nrow(Neighbourhood),1)
    for (i in 1:nrow(Neighbourhood))
    {
      dis[i,]=unique(eco_distance[eco_distance[,1]==as.character(a[,3])&eco_distance[,2]==as.character(Neighbourhood[i,3])|eco_distance[,1]==as.character(Neighbourhood[i,3]) & eco_distance[,2]==as.character(a[,3]),3])

    }
    Neighbourhood=cbind(Neighbourhood,dis)
    colnames(Neighbourhood) = c("x","Y","Species","Fun|Phy_Distance")
    Neighbourhood_pattern_Fun.phy_Richness=sum(Neighbourhood[,4]/n)
    output=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_Fun.phy_Richness=Neighbourhood_pattern_Fun.phy_Richness)
    output
  }#2
  d=matrix(NA,nrow(a),3)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Fun_Phy.single(a[j,],b,eco_distance,n)$Neighbourhood_pattern_Fun.phy_Richness))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Neighbourhood_pattern_Fun.phy_Richness")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树的邻体模式系统发育多样性
Neighbourhood.pattern.Fun_Phy.single=function(a,b,eco_distance,n)
{#2
  Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
  {#3
    a1=a
    b1=b
    a=a[,1:2]
    b=b[,1:2]
    c=b
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,]-a[1,])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d,b1[,3])
    d=subset(d,d>0)
    d=d[order(d[,3])[1:n],1:4]
    d=d[,c(1:2,4)]
    colnames(d) = c("x","Y","Species")
    rownames(d)=1:n
    d
  }#3
  Neighbourhood=Neighbourhood.single(a,b,n)
  dis=matrix(NA,nrow(Neighbourhood),1)
  for (i in 1:nrow(Neighbourhood))
  {
    dis[i,]=unique(eco_distance[eco_distance[,1]==as.character(a[,3])&eco_distance[,2]==as.character(Neighbourhood[i,3])|eco_distance[,1]==as.character(Neighbourhood[i,3]) & eco_distance[,2]==as.character(a[,3]),3])

  }
  Neighbourhood=cbind(Neighbourhood,dis)
  colnames(Neighbourhood) = c("x","Y","Species","Fun|Phy_Distance")
  Neighbourhood_pattern_Fun.phy_Richness=sum(Neighbourhood[,4]/n)
  output=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_Fun.phy_Richness=Neighbourhood_pattern_Fun.phy_Richness)
  output
}#2
####计算全林分邻体模式开敞度function~
Neighbourhood.pattern.K.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.K.single=function(a,b,n)
  {
    Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1:n],1:4]
      colnames(d) = c("x","Y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,n)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    Neighbourhood1=Neighbourhood
    Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
    Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
    k1=Neighbourhood1[,4]/Neighbourhood1[,5]
    k1=k1/(k1+1)
    k2=as.data.frame(k1)
    K=sum(k2,na.rm=T)
    K=K/n
    Neighbourhood_pattern_K=K
    outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
    outcome
  }

  d=matrix(NA,nrow(a),4)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.K.single(a[j,],b,n)$Neighbourhood_pattern_K))
    info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time = Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time= end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Size","Neighbourhood_pattern_K")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树邻体模式开敞度function~
Neighbourhood.pattern.K.single=function(a,b,n)
{
  Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
  {#3

    c=b[,1:2]
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,1:2]-a[1,1:2])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d)
    d=subset(d,d>0)
    d=d[order(d[,4])[1:n],1:4]
    colnames(d) = c("x","Y","Size","Distance")
    d
  }
  Nei.tree=Neighbourhood.single1(a,b,n)
  Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
  colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
  Neighbourhood=Nei.dif.high
  Neighbourhood1=Neighbourhood
  Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
  Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
  k1=Neighbourhood1[,4]/Neighbourhood1[,5]
  k1=k1/(k1+1)
  k2=as.data.frame(k1)
  K=sum(k2,na.rm=T)
  K=K/n
  Neighbourhood_pattern_K=K
  outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
  outcome
}
####计算全林分混交度function~
Neighbourhood.pattern.M.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.M.single=function(a,b,n)###计算a1林木的混交度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      d=d[,c(1:2,4)]
      colnames(d) = c("x","Y","Species")
      rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    Neighbourhood1=subset(Neighbourhood,Species!=a[,3])
    M=nrow(Neighbourhood1)/n
    output=list(a=a,Neighbourhood=Neighbourhood,M=M)
    output
  }#2

  d=matrix(NA,nrow(a),4)
  pb = tkProgressBar("进度","已完成 %", 0, 100)
  star_time = Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.M.single(a[j,],b,n)$M))
    info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time = Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time= end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Species","M")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参数树混交度
Neighbourhood.pattern.M.single=function(a,b,n)###计算a1林木的混交度
{#2
  Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
  {#3
    a1=a
    b1=b
    a=a[,1:2]
    b=b[,1:2]
    c=b
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,]-a[1,])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d,b1[,3])
    d=subset(d,d>0)
    d=d[order(d[,3])[1:n],1:4]
    d=d[,c(1:2,4)]
    colnames(d) = c("x","y","Species")
    #rownames(d)=1:n
    d
  }#3
  Neighbourhood=Neighbourhood.single(a,b,n)
  Neighbourhood1=subset(Neighbourhood,Species!=a[,3])
  M=nrow(Neighbourhood1)/n
  output=list(a=a,Neighbourhood=Neighbourhood,M=M)
  output
}#2
####计算全林分邻体模式荫蔽度function~
Neighbourhood.pattern.S.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.S.single=function(a,b,n)
  {
    Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1:n],1:4]
      colnames(d) = c("x","Y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,n)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    Neighbourhood1=Neighbourhood
    Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
    Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
    Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
    S1=Neighbourhood1[,5]/Neighbourhood1[,4]
    S1=S1/(S1+1)
    S2=as.data.frame(S1)
    S=sum(S2,na.rm=T)
    S=S/n
    Neighbourhood_pattern_S=S
    outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
    outcome
  }

  d=matrix(NA,nrow(a),4)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.S.single(a[j,],b,n)$Neighbourhood_pattern_S))
    info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time = Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time= end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Size","Neighbourhood_pattern_S")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树的邻体模式荫蔽度function~
Neighbourhood.pattern.S.single=function(a,b,n)
{
  Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
  {#3

    c=b[,1:2]
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,1:2]-a[1,1:2])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d)
    d=subset(d,d>0)
    d=d[order(d[,4])[1:n],1:4]
    colnames(d) = c("x","Y","Size","Distance")
    d
  }
  Nei.tree=Neighbourhood.single1(a,b,n)
  Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
  colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
  Neighbourhood=Nei.dif.high
  Neighbourhood1=Neighbourhood
  Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
  Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
  Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
  S1=Neighbourhood1[,5]/Neighbourhood1[,4]
  S1=S1/(S1+1)
  S2=as.data.frame(S1)
  S=sum(S2,na.rm=T)
  S=S/n
  Neighbourhood_pattern_S=S
  outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
  outcome
}
####计算全林分的大小比function~
Neighbourhood.pattern.U.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.U.single=function(a,b,n)###计算a1林木的大小比
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      d=d[,c(1:2,4)]
      colnames(d) = c("x","y","Size")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    U=subset(Neighbourhood,Size>a[,3])
    U=nrow(U)/n
    outcome=list(a=a,Neighbourhood=Neighbourhood,U=U)
    outcome
  }#2
  d=matrix(NA,nrow(a),4)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.U.single(a[j,],b,n)$U))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","size","U")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树的大小比function~
Neighbourhood.pattern.U.single=function(a,b,n)###计算a1林木的大小比
{#2
  Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
  {#3
    a1=a
    b1=b
    a=a[,1:2]
    b=b[,1:2]
    c=b
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,]-a[1,])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d,b1[,3])
    d=subset(d,d>0)
    d=d[order(d[,3])[1:n],1:4]
    d=d[,c(1:2,4)]
    colnames(d) = c("x","y","Size")
    #rownames(d)=1:n
    d
  }#3
  Neighbourhood=Neighbourhood.single(a,b,n)
  U=subset(Neighbourhood,Size>a[,3])
  U=nrow(U)/n
  outcome=list(a=a,Neighbourhood=Neighbourhood,U=U)
  outcome
}#2
####计算全林分角尺度function~
Neighbourhood.pattern.W.mult=function(a,b,n)
{#1
  library(tcltk)
  Neighbourhood.pattern.single=function(a,b,n)###计算a1林木的角尺度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:3]

      colnames(d) = c("x","Y","Distance")
      rownames(d)=1:n
      d
    }#3
    threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
    slope=matrix(NA,n,1)
    e=matrix(NA,n,n)
    rb=matrix(NA,n,n)
    residual=matrix(NA,n,1)
    stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]####校正为以a为原点，单位长度的向量
    stand1[stand1==0]=0.00000001####防止slope的分母为0
    stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]####校正为以a为原点，单位长度的向量
    standard=cbind(stand1,stand2)####校正为以a为原点，单位长度的向量

    for(i in 1:n)
    {#5
      slope[i,]=standard[i,2]/standard[i,1]####求tan角

      for(j in 1:n)
      {
        residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
      }
      e[,i]=residual
      e[abs(e)<0.000000001]=0
    }#5
    for(i in 3:(n+2))
    {
      standard.e=cbind(standard,e)
      pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
      pea=pe[,1:2]%*%standard[(i-2),]
      pea=pea[order(pea,decreasing = T)]
      ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
      nea=ne[,1:2]%*%standard[(i-2),]
      nea=nea[order(nea)]
      rbe=rbind(as.matrix(pea),as.matrix(nea))
      rb[,(i-2)]=rbe
    }
    rbn2=c(rb[2,],rb[n,])
    angle=rbn2-threshold
    ang=subset(matrix(angle),matrix(angle)>=0)
    num=length(ang)
    vaule=num/(2*n)

    vaule
  }#2
  ########
  d=matrix(NA,nrow(a),3)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.single(a[j,],b,n)))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","W")
  rownames(d)=1:j
  d
}#1
####计算单一参照树角尺度function~
Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
{#2
  Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
  {#3
    c=b
    for (i in 1:nrow(b))
    {#4
      c[i,]=(b[i,]-a[1,])^2
      d=(c[,1]+c[,2])^(1/2)
    }#4

    d=cbind(b,d)
    d=subset(d,d>0)
    d=d[order(d[,3])[1:n],1:3]

    colnames(d) = c("x","y","Distance")
    #rownames(d)=1:n
    d
  }#3
  Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
  threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
  slope=matrix(NA,n,1)
  e=matrix(NA,n,n)
  rb=matrix(NA,n,n)
  residual=matrix(NA,n,1)
  stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
  stand1[stand1==0]=0.00000001
  stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
  standard=cbind(stand1,stand2)

  for(i in 1:n)
  {#5
    slope[i,]=standard[i,2]/standard[i,1]

    for(j in 1:n)
    {
      residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
    }
    e[,i]=residual
    e[abs(e)<0.000000001]=0
  }#5
  for(i in 3:(n+2))
  {
    standard.e=cbind(standard,e)
    pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
    pea=pe[,1:2]%*%standard[(i-2),]
    pea=pea[order(pea,decreasing = T)]
    ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
    nea=ne[,1:2]%*%standard[(i-2),]
    nea=nea[order(nea)]
    rbe=rbind(as.matrix(pea),as.matrix(nea))
    rb[,(i-2)]=rbe
  }
  rbn2=c(rb[2,],rb[n,])
  angle=rbn2-threshold
  ang=subset(matrix(angle),matrix(angle)>0)
  num=length(ang)
  W=num/(2*n)
  outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
  outcome
}#2
####计算全林分加权密集度function~
Neighbourhood.pattern.Wei_C.mult=function(a,b)
{#
  library(tcltk)
  Neighbourhood.pattern.Wei_C.single=function(a,b)
  {#
    #####密集度
    Neighbourhood.pattern.C.single=function(a,b,n)###计算a1林木的密集度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        colnames(d) = c("x","y","distance","size")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
      C=nrow(C)/n
      outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
      outcome
    }#2

    C=Neighbourhood.pattern.C.single(a,b,4)$C

    ####角尺度
    Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","y","Distance")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      W=num/(2*n)
      outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
      outcome
    }#2
    wei=Neighbourhood.pattern.W.single(a,b,4)$W
    W=wei
    if(wei==0) {Wei=1}
    if(wei==0.25) {Wei=0.75}
    if(wei==0.5) {Wei=0.5}
    if(wei==0.75) {Wei=0.375}
    if(wei==1) {Wei=0.25}

    Wei_C=Wei*C
    Neighbourhood=Neighbourhood.pattern.C.single(a,b,4)$Neighbourhood
    output=list(a=a,Neighbourhood=Neighbourhood,C=C,W=W,Wei_C=Wei_C)
    output
  }
  d=matrix(NA,nrow(a),3)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_C.single(a[j,],b)$Wei_C))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Wei_C")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树的加权密集度function~
Neighbourhood.pattern.Wei_C.single=function(a,b)
{#
  #####密集度
  Neighbourhood.pattern.C.single=function(a,b,n)###计算a1林木的密集度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      colnames(d) = c("x","y","distance","size")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
    C=nrow(C)/n
    outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
    outcome
  }#2

  C=Neighbourhood.pattern.C.single(a,b,4)$C

  ####角尺度
  Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:3]

      colnames(d) = c("x","y","Distance")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
    threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
    slope=matrix(NA,n,1)
    e=matrix(NA,n,n)
    rb=matrix(NA,n,n)
    residual=matrix(NA,n,1)
    stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
    stand1[stand1==0]=0.00000001
    stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
    standard=cbind(stand1,stand2)

    for(i in 1:n)
    {#5
      slope[i,]=standard[i,2]/standard[i,1]

      for(j in 1:n)
      {
        residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
      }
      e[,i]=residual
      e[abs(e)<0.000000001]=0
    }#5
    for(i in 3:(n+2))
    {
      standard.e=cbind(standard,e)
      pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
      pea=pe[,1:2]%*%standard[(i-2),]
      pea=pea[order(pea,decreasing = T)]
      ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
      nea=ne[,1:2]%*%standard[(i-2),]
      nea=nea[order(nea)]
      rbe=rbind(as.matrix(pea),as.matrix(nea))
      rb[,(i-2)]=rbe
    }
    rbn2=c(rb[2,],rb[n,])
    angle=rbn2-threshold
    ang=subset(matrix(angle),matrix(angle)>0)
    num=length(ang)
    W=num/(2*n)
    outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
    outcome
  }#2
  wei=Neighbourhood.pattern.W.single(a,b,4)$W
  W=wei
  if(wei==0) {Wei=1}
  if(wei==0.25) {Wei=0.75}
  if(wei==0.5) {Wei=0.5}
  if(wei==0.75) {Wei=0.375}
  if(wei==1) {Wei=0.25}

  Wei_C=Wei*C
  Neighbourhood=Neighbourhood.pattern.C.single(a,b,4)$Neighbourhood
  output=list(a=a,Neighbourhood=Neighbourhood,C=C,W=W,Wei_C=Wei_C)
  output
}
####计算全林分加权邻体模式开敞度function~
Neighbourhood.pattern.Wei_K.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.Wei_K.single=function(a,b,n)
  {
    #####加载结构单元开敞度
    Neighbourhood.pattern.K.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
      k1=Neighbourhood1[,4]/Neighbourhood1[,5]
      k1=k1/(k1+1)
      k2=as.data.frame(k1)
      K=sum(k2,na.rm=T)
      K=K/n
      Neighbourhood_pattern_K=K
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
      outcome
    }
    #####加载角尺度
    Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","y","Distance")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      W=num/(2*n)
      outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
      outcome
    }#2

    ######计算角尺度与结构单元开敞度
    K=Neighbourhood.pattern.K.single(a,b,n)
    W=Neighbourhood.pattern.W.single(a,b,n)
    a=K$a
    Neighbourhood=K$Neighbourhood
    Neighbourhood_pattern_K=K$Neighbourhood_pattern_K
    Neighbourhood_pattern_Wei_K=K$Neighbourhood_pattern_K*(W$W^(1-K$Neighbourhood_pattern_K))
    W=W$W
    outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_K=Neighbourhood_pattern_K,W=W,  Neighbourhood_pattern_Wei_K=  Neighbourhood_pattern_Wei_K)
    outcome
  }
  d=matrix(NA,nrow(a),3)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_K.single(a[j,],b,n)$Neighbourhood_pattern_Wei_K))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Neighbourhood_pattern_Wei_K")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树的加权邻体模式开敞度function~
Neighbourhood.pattern.Wei_K.single=function(a,b,n)
{
  #####加载结构单元开敞度
  Neighbourhood.pattern.K.single=function(a,b,n)
  {
    Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1:n],1:4]
      colnames(d) = c("x","Y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,n)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    Neighbourhood1=Neighbourhood
    Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
    Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
    k1=Neighbourhood1[,4]/Neighbourhood1[,5]
    k1=k1/(k1+1)
    k2=as.data.frame(k1)
    K=sum(k2,na.rm=T)
    K=K/n
    Neighbourhood_pattern_K=K
    outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
    outcome
  }
  #####加载角尺度
  Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:3]

      colnames(d) = c("x","y","Distance")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
    threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
    slope=matrix(NA,n,1)
    e=matrix(NA,n,n)
    rb=matrix(NA,n,n)
    residual=matrix(NA,n,1)
    stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
    stand1[stand1==0]=0.00000001
    stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
    standard=cbind(stand1,stand2)

    for(i in 1:n)
    {#5
      slope[i,]=standard[i,2]/standard[i,1]

      for(j in 1:n)
      {
        residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
      }
      e[,i]=residual
      e[abs(e)<0.000000001]=0
    }#5
    for(i in 3:(n+2))
    {
      standard.e=cbind(standard,e)
      pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
      pea=pe[,1:2]%*%standard[(i-2),]
      pea=pea[order(pea,decreasing = T)]
      ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
      nea=ne[,1:2]%*%standard[(i-2),]
      nea=nea[order(nea)]
      rbe=rbind(as.matrix(pea),as.matrix(nea))
      rb[,(i-2)]=rbe
    }
    rbn2=c(rb[2,],rb[n,])
    angle=rbn2-threshold
    ang=subset(matrix(angle),matrix(angle)>0)
    num=length(ang)
    W=num/(2*n)
    outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
    outcome
  }#2

  ######计算角尺度与结构单元开敞度
  K=Neighbourhood.pattern.K.single(a,b,n)
  W=Neighbourhood.pattern.W.single(a,b,n)
  a=K$a
  Neighbourhood=K$Neighbourhood
  Neighbourhood_pattern_K=K$Neighbourhood_pattern_K
  Neighbourhood_pattern_Wei_K=K$Neighbourhood_pattern_K*(W$W^(1-K$Neighbourhood_pattern_K))
  W=W$W
  outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_K=Neighbourhood_pattern_K,W=W,  Neighbourhood_pattern_Wei_K=  Neighbourhood_pattern_Wei_K)
  outcome
}
####计算全林分的加权邻体模式荫蔽度function~
Neighbourhood.pattern.Wei_S.mult=function(a,b,n)
{
  library(tcltk)
  Neighbourhood.pattern.Wei_S.single=function(a,b,n)
  {
    ######加载邻体模式荫蔽度函数
    Neighbourhood.pattern.S.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
      Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
      S1=Neighbourhood1[,5]/Neighbourhood1[,4]
      S1=S1/(S1+1)
      S2=as.data.frame(S1)
      S=sum(S2,na.rm=T)
      S=S/n
      Neighbourhood_pattern_S=S
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
      outcome
    }

    #####加载角尺度计算函数
    Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","y","Distance")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      W=num/(2*n)
      outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
      outcome
    }#2

    #####计算邻体模式加权荫蔽度
    S=Neighbourhood.pattern.S.single(a,b,n)
    W=Neighbourhood.pattern.W.single(a,b,n)
    a=S$a
    Neighbourhood=S$Neighbourhood
    Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
    Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
    Neighbourhood_pattern_Wei_S=Neighbourhood_pattern_S*(1/n)^(W$W)
    W=W$W
    outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_S=Neighbourhood_pattern_S,W=W,  Neighbourhood_pattern_Wei_S=  Neighbourhood_pattern_Wei_S)
    outcome
  }
  d=matrix(NA,nrow(a),3)
  pb=tkProgressBar("进度","已完成 %", 0, 100)
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_S.single(a[j,],b,n)$Neighbourhood_pattern_Wei_S))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time =Sys.time()  ## 记录程序结束时间
  close(pb)
  run_time = end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","Neighbourhood_pattern_Wei_S")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  d
}
####计算单一参照树的加权邻体模式荫蔽度function~
Neighbourhood.pattern.Wei_S.single=function(a,b,n)
{
  ######加载邻体模式荫蔽度函数
  Neighbourhood.pattern.S.single=function(a,b,n)
  {
    Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1:n],1:4]
      colnames(d) = c("x","y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,n)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    Neighbourhood1=Neighbourhood
    Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
    Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
    Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
    S1=Neighbourhood1[,5]/Neighbourhood1[,4]
    S1=S1/(S1+1)
    S2=as.data.frame(S1)
    S=sum(S2,na.rm=T)
    S=S/n
    Neighbourhood_pattern_S=S
    outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
    outcome
  }

  #####加载角尺度计算函数
  Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:3]

      colnames(d) = c("x","y","Distance")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
    threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
    slope=matrix(NA,n,1)
    e=matrix(NA,n,n)
    rb=matrix(NA,n,n)
    residual=matrix(NA,n,1)
    stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
    stand1[stand1==0]=0.00000001
    stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
    standard=cbind(stand1,stand2)

    for(i in 1:n)
    {#5
      slope[i,]=standard[i,2]/standard[i,1]

      for(j in 1:n)
      {
        residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
      }
      e[,i]=residual
      e[abs(e)<0.000000001]=0
    }#5
    for(i in 3:(n+2))
    {
      standard.e=cbind(standard,e)
      pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
      pea=pe[,1:2]%*%standard[(i-2),]
      pea=pea[order(pea,decreasing = T)]
      ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
      nea=ne[,1:2]%*%standard[(i-2),]
      nea=nea[order(nea)]
      rbe=rbind(as.matrix(pea),as.matrix(nea))
      rb[,(i-2)]=rbe
    }
    rbn2=c(rb[2,],rb[n,])
    angle=rbn2-threshold
    ang=subset(matrix(angle),matrix(angle)>0)
    num=length(ang)
    W=num/(2*n)
    outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
    outcome
  }#2

  #####计算邻体模式加权荫蔽度
  S=Neighbourhood.pattern.S.single(a,b,n)
  W=Neighbourhood.pattern.W.single(a,b,n)
  a=S$a
  Neighbourhood=S$Neighbourhood
  Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
  Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
  Neighbourhood_pattern_Wei_S=Neighbourhood_pattern_S*(1/n)^(W$W)
  W=W$W
  outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_S=Neighbourhood_pattern_S,W=W,  Neighbourhood_pattern_Wei_S=  Neighbourhood_pattern_Wei_S)
  outcome
}
####绘制单一参照树开敞度图function~
plot.K.single=function(a,b)
{
  library(ggplot2)
  b1=subset(b,x>=a[,1]&y>a[,2]&size>a[,3])####将b划分成四个象限
  b2=subset(b,x>a[,1]&y<=a[,2]&size>a[,3])
  b3=subset(b,x<=a[,1]&y<a[,2]&size>a[,3])
  b4=subset(b,x<a[,1]&y>=a[,2]&size>a[,3])
  options(warn=-1)
  K.single=function(a,b)
  {
    Neighbourhood.K.single1=function(a,b)####计算单株林木a的开敞度
    {
      b1=subset(b,x>=a[,1]&y>a[,2]&size>a[,3])####将b划分成四个象限
      b2=subset(b,x>a[,1]&y<=a[,2]&size>a[,3])
      b3=subset(b,x<=a[,1]&y<a[,2]&size>a[,3])
      b4=subset(b,x<a[,1]&y>=a[,2]&size>a[,3])
      Neighbourhood.single=function(a,b)###找出a1林木的临近个体
      {#3
        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1],1:4]

        colnames(d) = c("x","Y","Size","Distance")
        d
      }#3
      b1.Nei.tree=Neighbourhood.single(a,b1)###分别找出四个象限中最近的树
      b2.Nei.tree=Neighbourhood.single(a,b2)
      b3.Nei.tree=Neighbourhood.single(a,b3)
      b4.Nei.tree=Neighbourhood.single(a,b4)
      Nei.tree=rbind(b1.Nei.tree,b2.Nei.tree,b3.Nei.tree,b4.Nei.tree)####葫芦娃合体
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      k1=Neighbourhood[,4]/Neighbourhood[,5]
      k2=as.data.frame(k1)
      k2[,1][is.infinite(k2[,1])] = NA
      K=sum(k2,na.rm=T)
      outcome=list(a=a,Neighbourhood=Neighbourhood,K=K)
      outcome
    }

    if(nrow(b1)!=0&nrow(b2)!=0&nrow(b3)!=0&nrow(b4)!=0) {Neighbourhood.K.single1(a,b) }else{}
  }
  options(warn=1)
  Ne=K.single(a,b)

  if(nrow(b1)!=0&nrow(b2)!=0&nrow(b3)!=0&nrow(b4)!=0) {
    center=cbind(cbind(rep(a[,1],each=n),rep(a[,2],each=n),rep(a[,3],each=n)),c(1:n))
    center=as.data.frame(center)
    colnames(center)=c("x","y","size","group")
    Neigh=cbind(Ne$Neighbourhood[,1:3],c(1:n))
    colnames(Neigh)=c("x","y","size","group")
    unitg=rbind(center,Neigh)
    max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
    p=ggplot()
    p=p+geom_hline(data=Ne$a,aes(yintercept=y),linetype=5,col="red")+geom_vline(data=Ne$a,aes(xintercept=x),linetype=5,col="red")
    p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_point(data=Ne$a,aes(x, y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Ne$Neighbourhood,aes(x, y),size=13/max(unitg[,3])*Ne$Neighbourhood[,3]+2,color="green4",alpha=0.7)
    p=p+geom_point(data=Ne$a,aes(x, y),size=2)+geom_point(data=Ne$Neighbourhood,aes(x, y),size=2)
    p=p+lims(x=c(Ne$a[,1]-max,Ne$a[,1]+max),y=c(Ne$a[,2]-max,Ne$a[,2]+max))
    p=p+ annotate("text",x=Ne$a[,1]-max+0.125*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("K=",round(Ne$K,3)))
    p+theme_bw()
  }
  else{warning("请确保以参考林木为中心的四个象限内均有高于其的立木")}
}
####区分物种绘制全林分密集度图function~
plot.Neighbourhood.pattern.C.dif_sp=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)
  colnames(a)=c("x","y","size")
  colnames(b)=c("x","y","size")
  Neighbourhood.pattern.C.mult=function(a,b,n)
  {
    Neighbourhood.pattern.C.single=function(a,b,n)
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        colnames(d) = c("x","y","distance","size")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
      C=nrow(C)/n
      outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
      outcome
    }#2
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.C.single(a[j,],b,n)$C))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","C")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.C.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","Neighbourhood Pattern C","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern C of diffence species")
  p2
}
####绘制全林分密集度图function~
plot.Neighbourhood.pattern.C.mult=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.C.mult=function(a,b,n)
  {
    Neighbourhood.pattern.C.single=function(a,b,n)
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        colnames(d) = c("x","y","distance","size")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
      C=nrow(C)/n
      outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
      outcome
    }#2
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.C.single(a[j,],b,n)$C))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","C")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }

  plot.data=Neighbourhood.pattern.C.mult(a,b,n)###计算画图所需数据
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern C of forest")
  p2
}
####绘制单一参照树的密集度function~
plot.Neighbourhood.pattern.C.single=function(a,b,n)
{
  library(ggplot2)
  library(ggforce)
  Neighbourhood.pattern.C.single=function(a,b,n)###计算a1林木的大小比
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      colnames(d) = c("x","y","distance","size")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
    C=nrow(C)/n
    outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
    outcome
  }#2
  Ne=Neighbourhood.pattern.C.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  Neigh=Neigh[,c(1,2,4,5)]
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  small=subset(Ne$Neighbourhood,distance<(Ne$Neighbourhood$size+a$size)/2)
  big=subset(Ne$Neighbourhood,distance>=(Ne$Neighbourhood$size+a$size)/2)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_circle(data=Ne$a,aes(x0=x, y0=y,r=size/2),fill="red4",color="black",alpha=0.7)+geom_circle(data=big,aes(x0=x, y0=y,r=size/2),color="black",fill="green4",alpha=0.7)+geom_circle(data=small,aes(x0=x, y0=y,r=size/2),fill="grey",color="black",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=big,aes(x,y),size=2)+geom_point(data=small,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.5*max,Ne$a[,1]+1.3*max),y=c(Ne$a[,2]-1.5*max,Ne$a[,2]+1.3*max))
  p=p+ annotate("text",x=Ne$a[,1]-1.3*max+0.125*(max),y=Ne$a[,2]+1.3*max-0.125*(max),label = paste0("C=",round(Ne$C,3)))
  p+theme_bw()
}
####区分物种绘制全林分邻体模式系统发育丰富度图function~
plot.Neighbourhood.pattern.Fun_Phy.dif_sp=function(a,b,eco_distance,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Fun_Phy.mult=function(a,b,eco_distance,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Fun_Phy.single=function(a,b,eco_distance,n)
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        d=d[,c(1:2,4)]
        colnames(d) = c("x","Y","Species")
        rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      dis=matrix(NA,nrow(Neighbourhood),1)
      for (i in 1:nrow(Neighbourhood))
      {
        dis[i,]=unique(eco_distance[eco_distance[,1]==as.character(a[,3])&eco_distance[,2]==as.character(Neighbourhood[i,3])|eco_distance[,1]==as.character(Neighbourhood[i,3]) & eco_distance[,2]==as.character(a[,3]),3])

      }
      Neighbourhood=cbind(Neighbourhood,dis)
      colnames(Neighbourhood) = c("x","Y","Species","Fun|Phy_Distance")
      Neighbourhood_pattern_Fun.phy_Richness=sum(Neighbourhood[,4]/n)
      output=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_Fun.phy_Richness=Neighbourhood_pattern_Fun.phy_Richness)
      output
    }#2
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Fun_Phy.single(a[j,],b,eco_distance,n)$Neighbourhood_pattern_Fun.phy_Richness))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Fun.phy_Richness")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Fun_Phy.mult(a,b,eco_distance,n)
  plot.data=cbind(plot.data,a$Species)
  plot.data=as.data.frame(plot.data)
  colnames(plot.data)=c("x","y","Neighbourhood_pattern_Fun.phy_Richness","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Fun.phy_Richness of diffence species")
  p2
}
#####绘制全林分的邻体模式系统发育丰富度function~
plot.Neighbourhood.pattern.Fun_Phy.mult=function(a,b,eco_distance,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Fun_Phy.mult=function(a,b,eco_distance,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Fun_Phy.single=function(a,b,eco_distance,n)
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        d=d[,c(1:2,4)]
        colnames(d) = c("x","Y","Species")
        rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      dis=matrix(NA,nrow(Neighbourhood),1)
      for (i in 1:nrow(Neighbourhood))
      {
        dis[i,]=unique(eco_distance[eco_distance[,1]==as.character(a[,3])&eco_distance[,2]==as.character(Neighbourhood[i,3])|eco_distance[,1]==as.character(Neighbourhood[i,3]) & eco_distance[,2]==as.character(a[,3]),3])

      }
      Neighbourhood=cbind(Neighbourhood,dis)
      colnames(Neighbourhood) = c("x","Y","Species","Fun|Phy_Distance")
      Neighbourhood_pattern_Fun.phy_Richness=sum(Neighbourhood[,4]/n)
      output=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_Fun.phy_Richness=Neighbourhood_pattern_Fun.phy_Richness)
      output
    }#2
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Fun_Phy.single(a[j,],b,eco_distance,n)$Neighbourhood_pattern_Fun.phy_Richness))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Fun.phy_Richness")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Fun_Phy.mult(a,b,eco_distance,n)###计算画图所需数据
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Fun.phy_Richness of forest")
  p2
}
####绘制单一参照树的邻体模式系统发育丰富度图function~
plot.Neighbourhood.pattern.Fun_Phy.single=function(a,b,eco_distance,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Fun_Phy.single=function(a,b,eco_distance,n)
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      d=d[,c(1:2,4)]
      colnames(d) = c("x","Y","Species")
      rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    dis=matrix(NA,nrow(Neighbourhood),1)
    for (i in 1:nrow(Neighbourhood))
    {
      dis[i,]=unique(eco_distance[eco_distance[,1]==as.character(a[,3])&eco_distance[,2]==as.character(Neighbourhood[i,3])|eco_distance[,1]==as.character(Neighbourhood[i,3]) & eco_distance[,2]==as.character(a[,3]),3])

    }
    Neighbourhood=cbind(Neighbourhood,dis)
    colnames(Neighbourhood) = c("x","Y","Species","Fun|Phy_Distance")
    Neighbourhood_pattern_Fun.phy_Richness=sum(Neighbourhood[,4]/n)
    output=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_Fun.phy_Richness=Neighbourhood_pattern_Fun.phy_Richness)
    output
  }#2
  Ne=Neighbourhood.pattern.Fun_Phy.single(a,b,eco_distance,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(0,each=n),c(1:n)))
  center=as.data.frame(center)
  colnames(center)=c("x","y","Fun.phy_distance","group")
  Neigh=cbind(Ne$Neighbourhood[,c(1:2,4)],c(1:n))
  center1=cbind(cbind(c(-100,-100),c(-100,-100),c(1:0),c(1:2)))
  center1=as.data.frame(center1)
  colnames(center1)=c("x","y","Fun.phy_distance","group")
  colnames(Neigh)=c("x","y","Fun.phy_distance","group")
  Neigh1=rbind(Neigh,center1)
  unitg=rbind(center,Neigh)

  options(warn=-1)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_point(data=Ne$a,aes(x, y),size=10,color="red4",alpha=0.7)+geom_point(data=Neigh1,aes(x, y,color=Fun.phy_distance),size=10,alpha=0.7)
  p=p+scale_colour_gradient(low="red4",high="green4")
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-max,Ne$a[,1]+max),y=c(Ne$a[,2]-max,Ne$a[,2]+max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.65*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("Neighbourhood_pattern_Fun.phy_Richness=",round(Ne$Neighbourhood_pattern_Fun.phy_Richness,3)))
  p=p+theme_bw()
  print(p)
  options(warn=1)
}
####区分物种绘制全林分邻体模式开敞度图function~
plot.Neighbourhood.pattern.K.dif_sp=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)
  Neighbourhood.pattern.K.mult=function(a,b,n)
  {
    Neighbourhood.pattern.K.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
      k1=Neighbourhood1[,4]/Neighbourhood1[,5]
      k1=k1/(k1+1)
      k2=as.data.frame(k1)
      K=sum(k2,na.rm=T)
      K=K/n
      Neighbourhood_pattern_K=K
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
      outcome
    }

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.K.single(a[j,],b,n)$Neighbourhood_pattern_K))
      info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time = Sys.time()  ## 记录程序结束时间
    close(pb)

    run_time <- end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Size","Neighbourhood_pattern_K")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.K.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","Size","Neighbourhood Pattern K","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern K of diffence species")
  p2
}
####绘制邻体模式开敞度克里金图function~
plot.Neighbourhood.pattern.K.Krig=function(minx,maxx,miny,maxy,b,seq,n)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  ####加载Neighbourhood.pattern.K.mult功能
  Neighbourhood.pattern.K.mult=function(a,b,n)
  {
    Neighbourhood.pattern.K.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        #d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
      k1=Neighbourhood1[,4]/Neighbourhood1[,5]
      k1=k1/(k1+1)
      k2=as.data.frame(k1)
      K=sum(k2,na.rm=T)
      K=K/n
      Neighbourhood_pattern_K=K
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
      outcome
    }

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.K.single(a[j,],b,n)$Neighbourhood_pattern_K))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Size","Neighbourhood_pattern_K")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Neighbourhood.pattern.K
  data=Neighbourhood.pattern.K.mult(basexy,b,n)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Neighbourhood_pattern_K")
  vgm1 <- variogram(Neighbourhood_pattern_K~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Neighbourhood_pattern_K~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title = "The neighbourhood pattern K")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
####绘制全林分邻体模式开敞度图function~
plot.Neighbourhood.pattern.K.mult=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.K.mult=function(a,b,n)
  {
    Neighbourhood.pattern.K.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
      k1=Neighbourhood1[,4]/Neighbourhood1[,5]
      k1=k1/(k1+1)
      k2=as.data.frame(k1)
      K=sum(k2,na.rm=T)
      K=K/n
      Neighbourhood_pattern_K=K
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
      outcome
    }

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.K.single(a[j,],b,n)$Neighbourhood_pattern_K))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Size","Neighbourhood_pattern_K")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.K.mult(a,b,n)###计算画图所需数据
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern K of forest")
  p2
}
####绘制单一参照树邻体模式开敞度图function~
plot.Neighbourhood.pattern.K.single=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.K.single=function(a,b,n)
  {
    Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1:n],1:4]
      colnames(d) = c("x","Y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,n)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    Neighbourhood1=Neighbourhood
    Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
    Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
    k1=Neighbourhood1[,4]/Neighbourhood1[,5]
    k1=k1/(k1+1)
    k2=as.data.frame(k1)
    K=sum(k2,na.rm=T)
    K=K/n
    Neighbourhood_pattern_K=K
    outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
    outcome
  }
  Ne=Neighbourhood.pattern.K.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  Neigh=Neigh[,c(1,2,3,6)]
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)
  p=p+geom_point(data=Ne$a,aes(x,y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Neigh,aes(x,y),size=13/max(unitg[,3])*Neigh[,3]+2,color="green4",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.1*max,Ne$a[,1]+1.1*max),y=c(Ne$a[,2]-1.1*max,Ne$a[,2]+1.1*max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.3*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("Neighbourhood_pattern_K=",round(Ne$Neighbourhood_pattern_K,3)))
  p+theme_bw()
}
####区分物种绘制全林分混交度图function~
plot.Neighbourhood.pattern.M.dif_spe=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.M.mult=function(a,b,n)
  {
    Neighbourhood.pattern.M.single=function(a,b,n)###计算a1林木的混交度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        d=d[,c(1:2,4)]
        colnames(d) = c("x","Y","Species")
        rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      Neighbourhood1=subset(Neighbourhood,Species!=a[,3])
      M=nrow(Neighbourhood1)/n
      output=list(a=a,Neighbourhood=Neighbourhood,M=M)
      output
    }#2

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.M.single(a[j,],b,n)$M))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("X","Y","Species","M")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }

  plot.data=Neighbourhood.pattern.M.mult(a,b,n)###计算画图所需数据

  for (i in c(1:2,4))
  {
    plot.data[,i]=as.numeric(as.character(plot.data[,i]))
  }
  p=ggplot(plot.data, aes(X, Y,color=Species))+geom_point(size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern M of different species")
  p2
}
####绘制全林分混交度图function~
plot.Neighbourhood.pattern.M.mult=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.M.mult=function(a,b,n)
  {
    Neighbourhood.pattern.M.single=function(a,b,n)###计算a1林木的混交度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        d=d[,c(1:2,4)]
        colnames(d) = c("x","Y","Species")
        rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      Neighbourhood1=subset(Neighbourhood,Species!=a[,3])
      M=nrow(Neighbourhood1)/n
      output=list(a=a,Neighbourhood=Neighbourhood,M=M)
      output
    }#2

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.M.single(a[j,],b,n)$M))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("X","Y","Species","M")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }

  plot.data=Neighbourhood.pattern.M.mult(a,b,n)###计算画图所需数据
  for (i in c(1:2,4))
  {
    plot.data[,i]=as.numeric(as.character(plot.data[,i]))
  }
  p=ggplot(plot.data, aes(X, Y))+geom_point(size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern M of forest")
  p2
}
####绘制单一参数树的混交度function~
plot.Neighbourhood.pattern.M.single=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.M.single=function(a,b,n)###计算a1林木的混交度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      d=d[,c(1:2,4)]
      colnames(d) = c("x","y","Species")
      rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    Neighbourhood1=subset(Neighbourhood,Species!=a[,3])
    M=nrow(Neighbourhood1)/n
    output=list(a=a,Neighbourhood=Neighbourhood,M=M)
    output
  }#2
  Ne=Neighbourhood.pattern.M.single(a,b,n)
  center=cbind(cbind(rep(a[,1],each=n),rep(a[,2],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","group")
  center1=center
  center1[,4]="red4"
  colnames(center1)=c("x","y","group","col")
  Neigh=cbind(Ne$Neighbourhood[,1:3],c(1:n))
  Neigh1=Neigh
  for(i in 1:n)
  {
    if(Neigh[i,3]==a[,3]){Neigh1[i,5]="red4"}
    else
    {Neigh1[i,5]="green4"}
  }
  colnames(Neigh)=c("x","y","Species","group")
  colnames(Neigh1)=c("x","y","Species","group","col")
  Neigh=Neigh[,c(1:2,4)]

  unitg=rbind(center,Neigh)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_point(data=a,aes(x, y),size=10,color="red4",alpha=0.7)+geom_point(data=Neigh1,aes(x, y),size=10,color=Neigh1$col,alpha=0.7)
  p=p+geom_point(data=a,aes(x, y),size=2,color="black")+geom_point(data=Neigh1,aes(x, y),size=2,color="black")
  p=p+lims(x=c(Ne$a[,1]-max,Ne$a[,1]+max),y=c(Ne$a[,2]-max,Ne$a[,2]+max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.125*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("M=",round(Ne$M,3)))
  p+theme_bw()
}
####区分物种绘制全林分邻体模式荫蔽度
plot.Neighbourhood.pattern.S.dif_sp=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)
  Neighbourhood.pattern.S.mult=function(a,b,n)
  {
    Neighbourhood.pattern.S.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
      Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
      S1=Neighbourhood1[,5]/Neighbourhood1[,4]
      S1=S1/(S1+1)
      S2=as.data.frame(S1)
      S=sum(S2,na.rm=T)
      S=S/n
      Neighbourhood_pattern_S=S
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
      outcome
    }
    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.S.single(a[j,],b,n)$Neighbourhood_pattern_S))
      info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time = Sys.time()  ## 记录程序结束时间
    close(pb)

    run_time <- end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Size","Neighbourhood_pattern_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.S.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","Size","Neighbourhood Pattern S","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern S of diffence species")
  p2
}
####绘制邻体模式荫蔽度克里金图function~
plot.Neighbourhood.pattern.S.Krig=function(minx,maxx,miny,maxy,b,seq,n)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  ####加载Neighbourhood.pattern.S.mult功能
  Neighbourhood.pattern.S.mult=function(a,b,n)
  {

    Neighbourhood.pattern.S.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
      Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
      S1=Neighbourhood1[,5]/Neighbourhood1[,4]
      S1=S1/(S1+1)
      S2=as.data.frame(S1)
      S=sum(S2,na.rm=T)
      S=S/n
      Neighbourhood_pattern_S=S
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
      outcome
    }

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.S.single(a[j,],b,n)$Neighbourhood_pattern_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Size","Neighbourhood_pattern_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Neighbourhood.pattern.S
  data=Neighbourhood.pattern.S.mult(basexy,b,n)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Neighbourhood_pattern_S")
  vgm1 <- variogram(Neighbourhood_pattern_S~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Neighbourhood_pattern_S~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title = "The neighbourhood pattern S")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
####绘制全林分邻体模式荫蔽度图function~
plot.Neighbourhood.pattern.S.mult=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.S.mult=function(a,b,n)
  {

    Neighbourhood.pattern.S.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
      Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
      S1=Neighbourhood1[,5]/Neighbourhood1[,4]
      S1=S1/(S1+1)
      S2=as.data.frame(S1)
      S=sum(S2,na.rm=T)
      S=S/n
      Neighbourhood_pattern_S=S
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
      outcome
    }

    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.S.single(a[j,],b,n)$Neighbourhood_pattern_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Size","Neighbourhood_pattern_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.S.mult(a,b,n)###计算画图所需数据
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern S of forest")
  p2
}
####绘制单一参照树邻体模式荫蔽度图function~
plot.Neighbourhood.pattern.S.single=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.S.single=function(a,b,n)
  {
    Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,4])[1:n],1:4]
      colnames(d) = c("x","Y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,n)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    Neighbourhood=Nei.dif.high
    Neighbourhood1=Neighbourhood
    Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
    Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
    Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
    S1=Neighbourhood1[,5]/Neighbourhood1[,4]
    S1=S1/(S1+1)
    S2=as.data.frame(S1)
    S=sum(S2,na.rm=T)
    S=S/n
    Neighbourhood_pattern_S=S
    outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
    outcome
  }

  Ne=Neighbourhood.pattern.S.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  Neigh=Neigh[,c(1,2,3,6)]
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)
  p=p+geom_point(data=Ne$a,aes(x,y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Neigh,aes(x,y),size=13/max(unitg[,3])*Neigh[,3]+2,color="green4",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.1*max,Ne$a[,1]+1.1*max),y=c(Ne$a[,2]-1.1*max,Ne$a[,2]+1.1*max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.3*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("Neighbourhood_pattern_S=",round(Ne$Neighbourhood_pattern_S,3)))
  p+theme_bw()
}
####区分物种绘制全林分大小比图function~
plot.Neighbourhood.pattern.U.dif_sp=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)
  colnames(a)=c("x","y","size")
  colnames(b)=c("x","y","size")
  Neighbourhood.pattern.U.mult=function(a,b,n)
  {
    Neighbourhood.pattern.U.single=function(a,b,n)###计算a1林木的大小比
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        d=d[,c(1:2,4)]
        colnames(d) = c("x","y","Size")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      U=subset(Neighbourhood,Size>a[,3])
      U=nrow(U)/n
      outcome=list(a=a,Neighbourhood=Neighbourhood,U=U)
      outcome
    }#2
    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.U.single(a[j,],b,n)$U))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","size","U")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.U.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","size","Neighbourhood Pattern U","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern U of diffence species")
  p2
}
####绘制全林分大小比图function~
plot.Neighbourhood.pattern.U.mult=function(a,b,n)
{
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.U.mult=function(a,b,n)
  {
    Neighbourhood.pattern.U.single=function(a,b,n)###计算a1林木的大小比
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        d=d[,c(1:2,4)]
        colnames(d) = c("x","y","Size")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      U=subset(Neighbourhood,Size>a[,3])
      U=nrow(U)/n
      outcome=list(a=a,Neighbourhood=Neighbourhood,U=U)
      outcome
    }#2
    d=matrix(NA,nrow(a),4)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.U.single(a[j,],b,n)$U))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","size","U")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }

  plot.data=Neighbourhood.pattern.U.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,4])*plot.data[,4]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern U of forest")
  p2
}
####绘制单一参照树大小比图function~
plot.Neighbourhood.pattern.U.single=function(a,b,n)###计算a1林木的大小比
{#2
  library(ggplot2)
  Neighbourhood.pattern.U.single=function(a,b,n)###计算a1林木的大小比
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      a1=a
      b1=b
      a=a[,1:2]
      b=b[,1:2]
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d,b1[,3])
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:4]
      d=d[,c(1:2,4)]
      colnames(d) = c("x","y","Size")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)
    U=subset(Neighbourhood,Size>a[,3])
    U=nrow(U)/n
    outcome=list(a=a,Neighbourhood=Neighbourhood,U=U)
    outcome
  }#2
  Ne=Neighbourhood.pattern.U.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  big=subset(Ne$Neighbourhood,Size>a[,3])
  small=subset(Ne$Neighbourhood,Size<=a[,3])
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_point(data=Ne$a,aes(x, y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=big,aes(x, y),size=13/max(unitg[,3])*big[,3]+2,color="green4",alpha=0.7)+geom_point(data=small,aes(x, y),size=13/max(unitg[,3])*small[,3]+2,color="grey",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=big,aes(x,y),size=2)+geom_point(data=small,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-max,Ne$a[,1]+max),y=c(Ne$a[,2]-max,Ne$a[,2]+max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.125*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("U=",round(Ne$U,3)))
  p+theme_bw()
}
####区分物种绘制全林分角尺度图function~
plot.Neighbourhood.pattern.W.dif_sp=function(a,b,n)
{####画角尺度图
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.mult=function(a,b,n)
  {#1
    ########
    ########
    ########
    a=a[,1:2]
    b=b[,1:2]
    Neighbourhood.pattern.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","Y","Distance")
        rownames(d)=1:n
        d
      }#3
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      vaule=num/(2*n)

      vaule
    }#2
    ########


    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.single(a[j,],b,n)))
      info= sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time = Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time= end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("X","Y","Neighbourhood Pattern W")
    rownames(d)=1:j
    d
  }#1
  plot.data=Neighbourhood.pattern.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a$Species)
  colnames(plot.data)=c("X","Y","Neighbourhood Pattern W","Species")
  p=ggplot()+geom_point(data=plot.data, aes(X, Y,color=Species),size=12/max(plot.data[,3])*plot.data[,3]+1,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern W of diffence species")
  p2
}
####绘制全林分角尺度图function~
plot.Neighbourhood.pattern.W.mult=function(a,b,n)
{####画角尺度图
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.mult=function(a,b,n)
  {#1
    ########
    ########
    ########
    Neighbourhood.pattern.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","Y","Distance")
        rownames(d)=1:n
        d
      }#3
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      vaule=num/(2*n)

      vaule
    }#2
    ########


    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(Neighbourhood.pattern.single(a[j,],b,n)))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("X","Y","Neighbourhood Pattern W")
    rownames(d)=1:j
    d
  }#1
  plot.data=Neighbourhood.pattern.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  p=ggplot(plot.data, aes(X, Y))+geom_point(size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern W of forest")
  p2
}
####绘制单一参照树的角尺度图function~
plot.Neighbourhood.pattern.W.single=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
  {#2
    Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
    {#3
      c=b
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,]-a[1,])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4

      d=cbind(b,d)
      d=subset(d,d>0)
      d=d[order(d[,3])[1:n],1:3]

      colnames(d) = c("x","y","Distance")
      #rownames(d)=1:n
      d
    }#3
    Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
    threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
    slope=matrix(NA,n,1)
    e=matrix(NA,n,n)
    rb=matrix(NA,n,n)
    residual=matrix(NA,n,1)
    stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
    stand1[stand1==0]=0.00000001
    stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
    standard=cbind(stand1,stand2)

    for(i in 1:n)
    {#5
      slope[i,]=standard[i,2]/standard[i,1]

      for(j in 1:n)
      {
        residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
      }
      e[,i]=residual
      e[abs(e)<0.000000001]=0
    }#5
    for(i in 3:(n+2))
    {
      standard.e=cbind(standard,e)
      pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
      pea=pe[,1:2]%*%standard[(i-2),]
      pea=pea[order(pea,decreasing = T)]
      ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
      nea=ne[,1:2]%*%standard[(i-2),]
      nea=nea[order(nea)]
      rbe=rbind(as.matrix(pea),as.matrix(nea))
      rb[,(i-2)]=rbe
    }
    rbn2=c(rb[2,],rb[n,])
    angle=rbn2-threshold
    ang=subset(matrix(angle),matrix(angle)>0)
    num=length(ang)
    W=num/(2*n)
    outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
    outcome
  }#2
  Ne=Neighbourhood.pattern.W.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  colnames(Neigh)=c("x","y","group")
  unitg=rbind(center,Neigh)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_point(data=Ne$a,aes(x, y),size=10,color="red4",alpha=0.7)+geom_point(data=Ne$Neighbourhood,aes(x, y),size=10,color="green4",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x, y),size=2)+geom_point(data=Ne$Neighbourhood,aes(x, y),size=2)
  p=p+lims(x=c(Ne$a[,1]-max,Ne$a[,1]+max),y=c(Ne$a[,2]-max,Ne$a[,2]+max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.125*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("W=",round(Ne$W,3)))
  p+theme_bw()
}
####区分物种绘制全林分加权密集度图function~
plot.Neighbourhood.pattern.Wei_C.dif_sp=function(a,b)
{#
  library(ggplot2)
  library(tcltk)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)
  colnames(a)=c("x","y","size")
  colnames(b)=c("x","y","size")
  Neighbourhood.pattern.Wei_C.mult=function(a,b)
  {#
    Neighbourhood.pattern.Wei_C.single=function(a,b)
    {#
      #####密集度
      Neighbourhood.pattern.C.single=function(a,b,n)###计算a1林木的密集度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          a1=a
          b1=b
          a=a[,1:2]
          b=b[,1:2]
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d,b1[,3])
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:4]
          colnames(d) = c("x","y","distance","size")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)
        C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
        C=nrow(C)/n
        outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
        outcome
      }#2

      C=Neighbourhood.pattern.C.single(a,b,4)$C

      ####角尺度
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2
      wei=Neighbourhood.pattern.W.single(a,b,4)$W
      W=wei
      if(wei==0) {Wei=1}
      if(wei==0.25) {Wei=0.75}
      if(wei==0.5) {Wei=0.5}
      if(wei==0.75) {Wei=0.375}
      if(wei==1) {Wei=0.25}

      Wei_C=Wei*C
      Neighbourhood=Neighbourhood.pattern.C.single(a,b,4)$Neighbourhood
      output=list(a=a,Neighbourhood=Neighbourhood,C=C,W=W,Wei_C=Wei_C)
      output
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_C.single(a[j,],b)$Wei_C))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Wei_C")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Wei_C.mult(a,b)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","Neighbourhood Pattern Wei_C","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Wei_C of diffence species")
  p2
}
####绘制全林分加权密集度图function~
plot.Neighbourhood.pattern.Wei_C.mult=function(a,b)
{#
  library(ggplot2)
  library(tcltk)
  Neighbourhood.pattern.Wei_C.mult=function(a,b)
  {#
    Neighbourhood.pattern.Wei_C.single=function(a,b)
    {#
      #####密集度
      Neighbourhood.pattern.C.single=function(a,b,n)###计算a1林木的密集度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          a1=a
          b1=b
          a=a[,1:2]
          b=b[,1:2]
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d,b1[,3])
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:4]
          colnames(d) = c("x","y","distance","size")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)
        C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
        C=nrow(C)/n
        outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
        outcome
      }#2

      C=Neighbourhood.pattern.C.single(a,b,4)$C

      ####角尺度
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2
      wei=Neighbourhood.pattern.W.single(a,b,4)$W
      W=wei
      if(wei==0) {Wei=1}
      if(wei==0.25) {Wei=0.75}
      if(wei==0.5) {Wei=0.5}
      if(wei==0.75) {Wei=0.375}
      if(wei==1) {Wei=0.25}

      Wei_C=Wei*C
      Neighbourhood=Neighbourhood.pattern.C.single(a,b,4)$Neighbourhood
      output=list(a=a,Neighbourhood=Neighbourhood,C=C,W=W,Wei_C=Wei_C)
      output
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_C.single(a[j,],b)$Wei_C))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Wei_C")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Wei_C.mult(a,b)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Wei_C of forest")
  p2
}
####绘制单一参照树的加权密集度图function~
plot.Neighbourhood.pattern.Wei_C.single=function(a,b)
{
  library(ggplot2)
  library(ggforce)
  Neighbourhood.pattern.Wei_C.single=function(a,b)
  {#
    #####密集度
    Neighbourhood.pattern.C.single=function(a,b,n)###计算a1林木的密集度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        a1=a
        b1=b
        a=a[,1:2]
        b=b[,1:2]
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d,b1[,3])
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:4]
        colnames(d) = c("x","y","distance","size")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)
      C=subset(Neighbourhood,distance<(Neighbourhood$size+a$size)/2)
      C=nrow(C)/n
      outcome=list(a=a,Neighbourhood=Neighbourhood,C=C)
      outcome
    }#2

    C=Neighbourhood.pattern.C.single(a,b,4)$C

    ####角尺度
    Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","y","Distance")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      W=num/(2*n)
      outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
      outcome
    }#2
    wei=Neighbourhood.pattern.W.single(a,b,4)$W
    W=wei
    if(wei==0) {Wei=1}
    if(wei==0.25) {Wei=0.75}
    if(wei==0.5) {Wei=0.5}
    if(wei==0.75) {Wei=0.375}
    if(wei==1) {Wei=0.25}

    Wei_C=Wei*C
    Neighbourhood=Neighbourhood.pattern.C.single(a,b,4)$Neighbourhood
    output=list(a=a,Neighbourhood=Neighbourhood,C=C,W=W,Wei_C=Wei_C)
    output
  }
  Ne=Neighbourhood.pattern.Wei_C.single(a,b)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  Neigh=Neigh[,c(1,2,4,5)]
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  small=subset(Ne$Neighbourhood,distance<(Ne$Neighbourhood$size+a$size)/2)
  big=subset(Ne$Neighbourhood,distance>=(Ne$Neighbourhood$size+a$size)/2)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_circle(data=Ne$a,aes(x0=x, y0=y,r=size/2),fill="red4",color="black",alpha=0.7)+geom_circle(data=big,aes(x0=x, y0=y,r=size/2),color="black",fill="green4",alpha=0.7)+geom_circle(data=small,aes(x0=x, y0=y,r=size/2),fill="grey",color="black",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=big,aes(x,y),size=2)+geom_point(data=small,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.5*max,Ne$a[,1]+1.3*max),y=c(Ne$a[,2]-1.5*max,Ne$a[,2]+1.3*max))
  p=p+ annotate("text",x=Ne$a[,1]-1.3*max+0.125*(max),y=Ne$a[,2]+1.3*max-0.125*(max),label = paste0("Weight_C=",round(Ne$Wei_C,3)))
  p+theme_bw()
}
####区分物种绘制全林分加权邻体模式开敞度function~
plot.Neighbourhood.pattern.Wei_K.dif_sp=function(a,b,n)
{
  library(ggplot2)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)
  Neighbourhood.pattern.Wei_K.mult=function(a,b,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Wei_K.single=function(a,b,n)
    {
      #####加载结构单元开敞度
      Neighbourhood.pattern.K.single=function(a,b,n)
      {
        Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,4])[1:n],1:4]
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,n)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
        k1=Neighbourhood1[,4]/Neighbourhood1[,5]
        k1=k1/(k1+1)
        k2=as.data.frame(k1)
        K=sum(k2,na.rm=T)
        K=K/n
        Neighbourhood_pattern_K=K
        outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
        outcome
      }
      #####加载角尺度
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2

      ######计算角尺度与结构单元开敞度
      K=Neighbourhood.pattern.K.single(a,b,n)
      W=Neighbourhood.pattern.W.single(a,b,n)
      a=K$a
      Neighbourhood=K$Neighbourhood
      Neighbourhood_pattern_K=K$Neighbourhood_pattern_K
      Neighbourhood_pattern_Wei_K=K$Neighbourhood_pattern_K*(W$W^(1-K$Neighbourhood_pattern_K))
      W=W$W
      outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_K=Neighbourhood_pattern_K,W=W,  Neighbourhood_pattern_Wei_K=  Neighbourhood_pattern_Wei_K)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_K.single(a[j,],b,n)$Neighbourhood_pattern_Wei_K))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Wei_K")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Wei_K.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","Neighbourhood Pattern Wei_K","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Wei_K of diffence species")
  p2
}
####绘制加权邻体模式开敞度克里金图function~
plot.Neighbourhood.pattern.Wei_K.Krig=function(minx,maxx,miny,maxy,b,seq,n)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Neighbourhood.pattern.Wei_K.mult=function(a,b,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Wei_K.single=function(a,b,n)
    {
      #####加载结构单元开敞度
      Neighbourhood.pattern.K.single=function(a,b,n)
      {
        Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,4])[1:n],1:4]
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,n)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
        k1=Neighbourhood1[,4]/Neighbourhood1[,5]
        k1=k1/(k1+1)
        k2=as.data.frame(k1)
        K=sum(k2,na.rm=T)
        K=K/n
        Neighbourhood_pattern_K=K
        outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
        outcome
      }
      #####加载角尺度
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2

      ######计算角尺度与结构单元开敞度
      K=Neighbourhood.pattern.K.single(a,b,n)
      W=Neighbourhood.pattern.W.single(a,b,n)
      a=K$a
      Neighbourhood=K$Neighbourhood
      Neighbourhood_pattern_K=K$Neighbourhood_pattern_K
      Neighbourhood_pattern_Wei_K=K$Neighbourhood_pattern_K*(W$W^(1-K$Neighbourhood_pattern_K))
      W=W$W
      outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_K=Neighbourhood_pattern_K,W=W,  Neighbourhood_pattern_Wei_K=  Neighbourhood_pattern_Wei_K)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_K.single(a[j,],b,n)$Neighbourhood_pattern_Wei_K))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Wei_K")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Neighbourhood.pattern.Wei_K
  data=Neighbourhood.pattern.Wei_K.mult(basexy,b,n)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Neighbourhood_pattern_Wei_K")
  vgm1 <- variogram(Neighbourhood_pattern_Wei_K~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Neighbourhood_pattern_Wei_K~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title = "The neighbourhood pattern Wei_K")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
####绘制全林分加权邻体模式开敞度图function~
plot.Neighbourhood.pattern.Wei_K.mult=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Wei_K.mult=function(a,b,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Wei_K.single=function(a,b,n)
    {
      #####加载结构单元开敞度
      Neighbourhood.pattern.K.single=function(a,b,n)
      {
        Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,4])[1:n],1:4]
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,n)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
        k1=Neighbourhood1[,4]/Neighbourhood1[,5]
        k1=k1/(k1+1)
        k2=as.data.frame(k1)
        K=sum(k2,na.rm=T)
        K=K/n
        Neighbourhood_pattern_K=K
        outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
        outcome
      }
      #####加载角尺度
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2

      ######计算角尺度与结构单元开敞度
      K=Neighbourhood.pattern.K.single(a,b,n)
      W=Neighbourhood.pattern.W.single(a,b,n)
      a=K$a
      Neighbourhood=K$Neighbourhood
      Neighbourhood_pattern_K=K$Neighbourhood_pattern_K
      Neighbourhood_pattern_Wei_K=K$Neighbourhood_pattern_K*(W$W^(1-K$Neighbourhood_pattern_K))
      W=W$W
      outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_K=Neighbourhood_pattern_K,W=W,  Neighbourhood_pattern_Wei_K=  Neighbourhood_pattern_Wei_K)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_K.single(a[j,],b,n)$Neighbourhood_pattern_Wei_K))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Wei_K")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Wei_K.mult(a,b,n)###计算画图所需数据
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Wei_K of forest")
  p2
}
####绘制单一参照树的加权邻体模式开敞度图function~
plot.Neighbourhood.pattern.Wei_K.single=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Wei_K.single=function(a,b,n)
  {
    #####加载结构单元开敞度
    Neighbourhood.pattern.K.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0.00000001
      k1=Neighbourhood1[,4]/Neighbourhood1[,5]
      k1=k1/(k1+1)
      k2=as.data.frame(k1)
      K=sum(k2,na.rm=T)
      K=K/n
      Neighbourhood_pattern_K=K
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_K=Neighbourhood_pattern_K)
      outcome
    }
    #####加载角尺度
    Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","y","Distance")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      W=num/(2*n)
      outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
      outcome
    }#2

    ######计算角尺度与结构单元开敞度
    K=Neighbourhood.pattern.K.single(a,b,n)
    W=Neighbourhood.pattern.W.single(a,b,n)
    a=K$a
    Neighbourhood=K$Neighbourhood
    Neighbourhood_pattern_K=K$Neighbourhood_pattern_K
    Neighbourhood_pattern_Wei_K=K$Neighbourhood_pattern_K*(W$W^(1-K$Neighbourhood_pattern_K))
    W=W$W
    outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_K=Neighbourhood_pattern_K,W=W,  Neighbourhood_pattern_Wei_K=  Neighbourhood_pattern_Wei_K)
    outcome
  }
  Ne=Neighbourhood.pattern.Wei_K.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  Neigh=Neigh[,c(1,2,3,6)]
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)
  p=p+geom_point(data=Ne$a,aes(x,y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Neigh,aes(x,y),size=13/max(unitg[,3])*Neigh[,3]+2,color="green4",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.1*max,Ne$a[,1]+1.1*max),y=c(Ne$a[,2]-1.1*max,Ne$a[,2]+1.1*max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.6*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("Neighbourhood_pattern_Wei_K=",round(Ne$Neighbourhood_pattern_Wei_K,3)))
  p+theme_bw()
}
####区分物种绘制全林分加权邻体模式荫蔽度图function~
plot.Neighbourhood.pattern.Wei_S.dif_sp=function(a,b,n)
{
  library(ggplot2)
  a1=a
  a=cbind(a$x,a$y,a$size)
  b=cbind(b$x,b$y,b$size)
  a=as.data.frame(a)
  b=as.data.frame(b)

  Neighbourhood.pattern.Wei_S.mult=function(a,b,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Wei_S.single=function(a,b,n)
    {
      ######加载邻体模式荫蔽度函数
      Neighbourhood.pattern.S.single=function(a,b,n)
      {
        Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,4])[1:n],1:4]
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,n)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
        Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
        S1=Neighbourhood1[,5]/Neighbourhood1[,4]
        S1=S1/(S1+1)
        S2=as.data.frame(S1)
        S=sum(S2,na.rm=T)
        S=S/n
        Neighbourhood_pattern_S=S
        outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
        outcome
      }

      #####加载角尺度计算函数
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2

      #####计算邻体模式加权荫蔽度
      S=Neighbourhood.pattern.S.single(a,b,n)
      W=Neighbourhood.pattern.W.single(a,b,n)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
      Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
      Neighbourhood_pattern_Wei_S=Neighbourhood_pattern_S*(1/n)^(W$W)
      W=W$W
      outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_S=Neighbourhood_pattern_S,W=W,  Neighbourhood_pattern_Wei_S=  Neighbourhood_pattern_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_S.single(a[j,],b,n)$Neighbourhood_pattern_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Wei_S.mult(a,b,n)###计算画图所需数据
  plot.data=as.data.frame(plot.data)
  plot.data=cbind(plot.data,a1$Species)
  colnames(plot.data)=c("x","y","Neighbourhood Pattern Wei_S","Species")
  p=ggplot()+geom_point(data=plot.data, aes(x, y,color=Species),size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1)
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Wei_S of diffence species")
  p2
}
####绘制加权邻体模式荫蔽度克里金图function~
plot.Neighbourhood.pattern.Wei_S.Krig=function(minx,maxx,miny,maxy,b,seq,n)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Neighbourhood.pattern.Wei_S.mult=function(a,b,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Wei_S.single=function(a,b,n)
    {
      ######加载邻体模式荫蔽度函数
      Neighbourhood.pattern.S.single=function(a,b,n)
      {
        Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          #d=subset(d,d>0)
          d=d[order(d[,4])[1:n],1:4]
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,n)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
        Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
        S1=Neighbourhood1[,5]/Neighbourhood1[,4]
        S1=S1/(S1+1)
        S2=as.data.frame(S1)
        S=sum(S2,na.rm=T)
        S=S/n
        Neighbourhood_pattern_S=S
        outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
        outcome
      }

      #####加载角尺度计算函数
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          # d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2

      #####计算邻体模式加权荫蔽度
      S=Neighbourhood.pattern.S.single(a,b,n)
      W=Neighbourhood.pattern.W.single(a,b,n)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
      Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
      Neighbourhood_pattern_Wei_S=Neighbourhood_pattern_S*(1/n)^(W$W)
      W=W$W
      outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_S=Neighbourhood_pattern_S,W=W,  Neighbourhood_pattern_Wei_S=  Neighbourhood_pattern_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_S.single(a[j,],b,n)$Neighbourhood_pattern_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Neighbourhood.pattern.Wei_S
  data=Neighbourhood.pattern.Wei_S.mult(basexy,b,n)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Neighbourhood_pattern_Wei_S")
  vgm1 <- variogram(Neighbourhood_pattern_Wei_S~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Neighbourhood_pattern_Wei_S~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title = "The neighbourhood pattern Wei_S")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
####绘制全林分加权邻体模式荫蔽度图function~
plot.Neighbourhood.pattern.Wei_S.mult=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Wei_S.mult=function(a,b,n)
  {
    library(tcltk)
    Neighbourhood.pattern.Wei_S.single=function(a,b,n)
    {
      ######加载邻体模式荫蔽度函数
      Neighbourhood.pattern.S.single=function(a,b,n)
      {
        Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,4])[1:n],1:4]
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,n)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
        Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
        S1=Neighbourhood1[,5]/Neighbourhood1[,4]
        S1=S1/(S1+1)
        S2=as.data.frame(S1)
        S=sum(S2,na.rm=T)
        S=S/n
        Neighbourhood_pattern_S=S
        outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
        outcome
      }

      #####加载角尺度计算函数
      Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
      {#2
        Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
        {#3
          c=b
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,]-a[1,])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4

          d=cbind(b,d)
          d=subset(d,d>0)
          d=d[order(d[,3])[1:n],1:3]

          colnames(d) = c("x","y","Distance")
          #rownames(d)=1:n
          d
        }#3
        Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
        threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
        slope=matrix(NA,n,1)
        e=matrix(NA,n,n)
        rb=matrix(NA,n,n)
        residual=matrix(NA,n,1)
        stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
        stand1[stand1==0]=0.00000001
        stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
        standard=cbind(stand1,stand2)

        for(i in 1:n)
        {#5
          slope[i,]=standard[i,2]/standard[i,1]

          for(j in 1:n)
          {
            residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
          }
          e[,i]=residual
          e[abs(e)<0.000000001]=0
        }#5
        for(i in 3:(n+2))
        {
          standard.e=cbind(standard,e)
          pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
          pea=pe[,1:2]%*%standard[(i-2),]
          pea=pea[order(pea,decreasing = T)]
          ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
          nea=ne[,1:2]%*%standard[(i-2),]
          nea=nea[order(nea)]
          rbe=rbind(as.matrix(pea),as.matrix(nea))
          rb[,(i-2)]=rbe
        }
        rbn2=c(rb[2,],rb[n,])
        angle=rbn2-threshold
        ang=subset(matrix(angle),matrix(angle)>0)
        num=length(ang)
        W=num/(2*n)
        outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
        outcome
      }#2

      #####计算邻体模式加权荫蔽度
      S=Neighbourhood.pattern.S.single(a,b,n)
      W=Neighbourhood.pattern.W.single(a,b,n)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
      Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
      Neighbourhood_pattern_Wei_S=Neighbourhood_pattern_S*(1/n)^(W$W)
      W=W$W
      outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_S=Neighbourhood_pattern_S,W=W,  Neighbourhood_pattern_Wei_S=  Neighbourhood_pattern_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Neighbourhood.pattern.Wei_S.single(a[j,],b,n)$Neighbourhood_pattern_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Neighbourhood_pattern_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  plot.data=Neighbourhood.pattern.Wei_S.mult(a,b,n)###计算画图所需数据
  p=ggplot(plot.data, aes(x, y))+geom_point(size=12/max(plot.data[,3])*plot.data[,3]+3,shape=1,color="red")
  p1=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()
  p2=p1+xlab("")+ylab("")+labs(title = "The neighbourhood pattern Wei_S of forest")
  p2
}
####绘制单一参照树加权邻体模式荫蔽度图function~
plot.Neighbourhood.pattern.Wei_S.single=function(a,b,n)
{
  library(ggplot2)
  Neighbourhood.pattern.Wei_S.single=function(a,b,n)
  {
    ######加载邻体模式荫蔽度函数
    Neighbourhood.pattern.S.single=function(a,b,n)
    {
      Neighbourhood.single1=function(a,b,n)###找出a1林木的临近个体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,4])[1:n],1:4]
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,n)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
      Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
      S1=Neighbourhood1[,5]/Neighbourhood1[,4]
      S1=S1/(S1+1)
      S2=as.data.frame(S1)
      S=sum(S2,na.rm=T)
      S=S/n
      Neighbourhood_pattern_S=S
      outcome=list(a=a,Neighbourhood=Neighbourhood,Neighbourhood_pattern_S=Neighbourhood_pattern_S)
      outcome
    }

    #####加载角尺度计算函数
    Neighbourhood.pattern.W.single=function(a,b,n)###计算a1林木的角尺度
    {#2
      Neighbourhood.single=function(a,b,n)###找出a1林木的临近个体
      {#3
        c=b
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,]-a[1,])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4

        d=cbind(b,d)
        d=subset(d,d>0)
        d=d[order(d[,3])[1:n],1:3]

        colnames(d) = c("x","y","Distance")
        #rownames(d)=1:n
        d
      }#3
      Neighbourhood=Neighbourhood.single(a,b,n)[,1:2]
      threshold=matrix(cos(2*pi/(n+1)),1,2*n)###计算角度临界值
      slope=matrix(NA,n,1)
      e=matrix(NA,n,n)
      rb=matrix(NA,n,n)
      residual=matrix(NA,n,1)
      stand1=(Neighbourhood.single(a,b,n)[,1]-a[1,1])/Neighbourhood.single(a,b,n)[,3]
      stand1[stand1==0]=0.00000001
      stand2=(Neighbourhood.single(a,b,n)[,2]-a[1,2])/Neighbourhood.single(a,b,n)[,3]
      standard=cbind(stand1,stand2)

      for(i in 1:n)
      {#5
        slope[i,]=standard[i,2]/standard[i,1]

        for(j in 1:n)
        {
          residual[j,]=slope[i,]*standard[j,1]-standard[j,2]
        }
        e[,i]=residual
        e[abs(e)<0.000000001]=0
      }#5
      for(i in 3:(n+2))
      {
        standard.e=cbind(standard,e)
        pe=subset(standard.e,standard.e[,i]>=0,c(1,2,i))
        pea=pe[,1:2]%*%standard[(i-2),]
        pea=pea[order(pea,decreasing = T)]
        ne=subset(standard.e,standard.e[,i]<0,c(1,2,i))
        nea=ne[,1:2]%*%standard[(i-2),]
        nea=nea[order(nea)]
        rbe=rbind(as.matrix(pea),as.matrix(nea))
        rb[,(i-2)]=rbe
      }
      rbn2=c(rb[2,],rb[n,])
      angle=rbn2-threshold
      ang=subset(matrix(angle),matrix(angle)>0)
      num=length(ang)
      W=num/(2*n)
      outcome=list(a=a,Neighbourhood=Neighbourhood,W=W)
      outcome
    }#2

    #####计算邻体模式加权荫蔽度
    S=Neighbourhood.pattern.S.single(a,b,n)
    W=Neighbourhood.pattern.W.single(a,b,n)
    a=S$a
    Neighbourhood=S$Neighbourhood
    Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
    Neighbourhood_pattern_S=S$Neighbourhood_pattern_S
    Neighbourhood_pattern_Wei_S=Neighbourhood_pattern_S*(1/n)^(W$W)
    W=W$W
    outcome=list(a=a,Neighbourhood=Neighbourhood, Neighbourhood_pattern_S=Neighbourhood_pattern_S,W=W,  Neighbourhood_pattern_Wei_S=  Neighbourhood_pattern_Wei_S)
    outcome
  }
  Ne=Neighbourhood.pattern.Wei_S.single(a,b,n)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(Ne$Neighbourhood,c(1:n))
  Neigh=Neigh[,c(1,2,3,6)]
  colnames(Neigh)=c("x","y","size","group")
  unitg=rbind(center,Neigh)
  max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  p=ggplot()
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)
  p=p+geom_point(data=Ne$a,aes(x,y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Neigh,aes(x,y),size=13/max(unitg[,3])*Neigh[,3]+2,color="green4",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.1*max,Ne$a[,1]+1.1*max),y=c(Ne$a[,2]-1.1*max,Ne$a[,2]+1.1*max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.6*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("Neighbourhood_pattern_Wei_S=",round(Ne$Neighbourhood_pattern_Wei_S,3)))
  p+theme_bw()
}
####计算尺度化加权荫蔽度修正系数function~
MI.Scale.S=function(minx,maxx,miny,maxy,b,seq,scale,MI)
{
  library(sp)
  library(gstat)
  library(tcltk)
  #################
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  a= dplyr::distinct(basexy)####删除重复样点
  ################
  Scale.dependence.S.mult=function(a,b,scale)
  {
    library(tcltk)
    Scale.dependence.S.single=function(a,b,scale,MI)
    {
      ######加载尺度依赖荫蔽度函数
      Scale.dependence.S.single=function(a,b,scale,MI)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3
          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          #d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        n=nrow(Nei.tree)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        n.dif=subset(Nei.dif.high,Nei.dif.high$Size_dif>0)
        n.dif=nrow(n.dif)
        if(n.dif==0)
        {Scale_dependence_S=0
        Neighbourhood=Nei.dif.high}
        if(n.dif!=0)
        {
          Neighbourhood=Nei.dif.high
          Neighbourhood1=Neighbourhood
          Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
          Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
          Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
          S1=Neighbourhood1[,5]/Neighbourhood1[,4]
          S2=as.data.frame(S1)
          S=sum(S2,na.rm=T)
          Scale_dependence_S=S
        }
        Scale_dependence_S=Scale_dependence_S
        outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S)
        outcome
      }

      #####计算尺度依赖加权荫蔽度
      S=Scale.dependence.S.single(a,b,scale)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Scale_dependence_S=S$Scale_dependence_S
      Scale_dependence_Wei_S=Scale_dependence_S
      #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
      outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S,Levins_Simpson=Levins_Simpson,Scale_dependence_Wei_S=Scale_dependence_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Scale.dependence.S.single(a[j,],b,scale)$Scale_dependence_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Scale_dependence_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  data=Scale.dependence.S.mult(a,b,scale)
  data=data[,3]
  n_=length(data)
  x=matrix()
  data=data[order(data)]
  for(i in 1:n_){x[i]=(atan(data[i]/(MI*pi))/pi*2)}
  data=1:n_
  lm=lm(data~x)
  print(plot(x,data)+abline(lm))
  print(summary(lm))
}
####计算单一参照树的尺度化加权荫蔽度function~
Scale.dependence.Wei_S.single=function(a,b,scale,MI)
{
  ######加载尺度依赖荫蔽度函数
  Scale.dependence.S.single=function(a,b,scale,MI)
  {
    Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
    {#3
      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4
      d=cbind(b,d)
      d=subset(d,d>0)
      d=subset(d,d<scale)
      colnames(d) = c("x","Y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,scale)
    n=nrow(Nei.tree)
    Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
    colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
    n.dif=subset(Nei.dif.high,Nei.dif.high$Size_dif>0)
    n.dif=nrow(n.dif)
    if(n.dif==0)
    {Scale_dependence_S=0
    Neighbourhood=Nei.dif.high}
    if(n.dif!=0)
    {
      Neighbourhood=Nei.dif.high
      Neighbourhood1=Neighbourhood
      Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
      Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
      Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
      S1=Neighbourhood1[,5]/Neighbourhood1[,4]
      S2=as.data.frame(S1)
      S=sum(S2,na.rm=T)
      Scale_dependence_S=S
    }
    Scale_dependence_S=(atan( Scale_dependence_S/(MI*pi))/pi*2)
    outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S)
    outcome
  }
  Scale.dependence.Simpson=function(a,b,scale)
  {
    Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
    {#3

      c=b[,1:2]
      for (i in 1:nrow(b))
      {#4
        c[i,]=(b[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }#4
      d=cbind(b,d)
      #d=subset(d,d>0)
      d=subset(d,d<scale)
      colnames(d) = c("x","y","Size","Distance")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,scale)
    Nei.tree
    Nei.tree[,3]=Nei.tree[,3]-a[,3]
    Nei.tree=subset(Nei.tree,Size>0)
    Nei.tree=subset(Nei.tree,Distance>0)
    Nei.tree
    Simpn=sum(Nei.tree[,3])
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
  S=Scale.dependence.S.single(a,b,scale,MI)
  Levins_Simpson=Scale.dependence.Simpson(a,b,scale)
  a=S$a
  Neighbourhood=S$Neighbourhood
  Scale_dependence_S=S$Scale_dependence_S
  Scale_dependence_Wei_S=Scale_dependence_S*(0.25/Levins_Simpson)
  #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
  outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S,Levins_Simpson=Levins_Simpson,Scale_dependence_Wei_S=Scale_dependence_Wei_S)
  outcome
}
#####绘制单一参照树的尺度化加权荫蔽度图function~
plot.Scale.dependence.Wei_S.single=function(a,b,scale,MI)
{
  library(ggplot2)
  library(ggforce)

  Scale.dependence.Wei_S.single=function(a,b,scale,MI)
  {
    ######加载尺度依赖荫蔽度函数
    Scale.dependence.S.single=function(a,b,scale,MI)
    {
      Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
      {#3
        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4
        d=cbind(b,d)
        d=subset(d,d>0)
        d=subset(d,d<scale)
        colnames(d) = c("x","Y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,scale)
      n=nrow(Nei.tree)
      Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
      colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
      n.dif=subset(Nei.dif.high,Nei.dif.high$Size_dif>0)
      n.dif=nrow(n.dif)
      if(n.dif==0)
      {Scale_dependence_S=0
      Neighbourhood=Nei.dif.high}
      if(n.dif!=0)
      {
        Neighbourhood=Nei.dif.high
        Neighbourhood1=Neighbourhood
        Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
        Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
        Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
        S1=Neighbourhood1[,5]/Neighbourhood1[,4]
        S2=as.data.frame(S1)
        S=sum(S2,na.rm=T)
        Scale_dependence_S=S
      }
      Scale_dependence_S=(atan( Scale_dependence_S/(MI*pi))/pi*2)
      outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S)
      outcome
    }
    Scale.dependence.Simpson=function(a,b,scale)
    {
      Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
      {#3

        c=b[,1:2]
        for (i in 1:nrow(b))
        {#4
          c[i,]=(b[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }#4
        d=cbind(b,d)
        #d=subset(d,d>0)
        d=subset(d,d<scale)
        colnames(d) = c("x","y","Size","Distance")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,scale)
      Nei.tree
      Nei.tree[,3]=Nei.tree[,3]-a[,3]
      Nei.tree=subset(Nei.tree,Size>0)
      Nei.tree=subset(Nei.tree,Distance>0)
      Nei.tree
      Simpn=sum(Nei.tree[,3])
      Nei.tree1=subset(Nei.tree,x>=a[,1]&y>a[,2])####将Nei.tree划分到四个象限
      Nei.tree2=subset(Nei.tree,x>a[,1]&y<=a[,2])
      Nei.tree3=subset(Nei.tree,x<=a[,1]&y<a[,2])
      Nei.tree4=subset(Nei.tree,x<a[,1]&y>=a[,2])
      Simpn1=sum(Nei.tree1[,3])
      Simpn2=sum(Nei.tree2[,3])
      Simpn3=sum(Nei.tree3[,3])
      Simpn4=sum(Nei.tree4[,3])
      simp_wei=sum((Simpn1/Simpn)^2,(Simpn2/Simpn)^2,(Simpn3/Simpn)^2,(Simpn4/Simpn)^2)
      if(is.nan(simp_wei)==T)(simp_wei=1)
      simp_wei
    }
    #####计算尺度依赖加权荫蔽度
    S=Scale.dependence.S.single(a,b,scale,MI)
    Levins_Simpson=Scale.dependence.Simpson(a,b,scale)
    a=S$a
    Neighbourhood=S$Neighbourhood
    Scale_dependence_S=S$Scale_dependence_S
    Scale_dependence_Wei_S=Scale_dependence_S*(0.25/Levins_Simpson)
    #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
    outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S,Levins_Simpson=Levins_Simpson,Scale_dependence_Wei_S=Scale_dependence_Wei_S)
    outcome
  }
  Ne=Scale.dependence.Wei_S.single(a,b,scale,MI)
  big=subset(Ne$Neighbourhood,Size>a[,3])
  small=subset(Ne$Neighbourhood,Size<=a[,3])
  n=nrow(big)
  center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
  center=as.data.frame(center)
  colnames(center)=c("x","y","size","group")
  Neigh=cbind(big,c(1:n))
  Neigh=Neigh[,c(1,2,3,6)]
  colnames(Neigh)=c("x","y","size","group")
  colnames(small)=c("x","y","size")
  unitg=rbind(center,Neigh)
  p=ggplot()
  p=p+geom_hline(data=Ne$a,aes(yintercept=y),linetype=5,col="red")+geom_vline(data=Ne$a,aes(xintercept=x),linetype=5,col="red")
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)+geom_circle(data=Ne$a,aes(x0=x, y0=y,r=scale),linetype=2,color="black",alpha=1)
  p=p+geom_point(data=Ne$a,aes(x,y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Neigh,aes(x,y),size=13/max(unitg[,3])*Neigh[,3]+2,color="green4",alpha=0.7)
  p=p+geom_point(data=small,aes(x,y),size=13/max(unitg[,3])*small[,3]+2,color="grey",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)+geom_point(data=small,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.1*scale,Ne$a[,1]+1.1*scale),y=c(Ne$a[,2]-1.1*scale,Ne$a[,2]+1.1*scale))
  p=p+ annotate("text",x=Ne$a[,1]-scale+0.5*(scale),y=Ne$a[,2]+scale-0.125*(scale),label = paste0("Scale_dependence_Wei_S=",round(Ne$Scale_dependence_Wei_S,3)))
  p+theme_bw()
}
####绘制尺度化加权荫蔽度克里金图function~
plot.Scale.dependence.Wei_S.Krig=function(minx,maxx,miny,maxy,b,seq,scale,MI)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Scale.dependence.Wei_S.mult=function(a,b,scale,MI)
  {
    library(tcltk)
    Scale.dependence.Wei_S.single=function(a,b,scale,MI)
    {
      ######加载尺度依赖荫蔽度函数
      Scale.dependence.S.single=function(a,b,scale,MI)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3
          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          #d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        n=nrow(Nei.tree)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        n.dif=subset(Nei.dif.high,Nei.dif.high$Size_dif>0)
        n.dif=nrow(n.dif)
        if(n.dif==0)
        {Scale_dependence_S=0
        Neighbourhood=Nei.dif.high}
        if(n.dif!=0)
        {
          Neighbourhood=Nei.dif.high
          Neighbourhood1=Neighbourhood
          Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
          Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
          Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
          S1=Neighbourhood1[,5]/Neighbourhood1[,4]
          S2=as.data.frame(S1)
          S=sum(S2,na.rm=T)
          Scale_dependence_S=S
        }
        Scale_dependence_S=(atan( Scale_dependence_S/(MI*pi))/pi*2)
        outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S)
        outcome
      }
      Scale.dependence.Simpson=function(a,b,scale)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        Nei.tree
        Nei.tree[,3]=Nei.tree[,3]-a[,3]
        Nei.tree=subset(Nei.tree,Size>0)
        Nei.tree=subset(Nei.tree,Distance>0)
        Nei.tree
        Simpn=sum(Nei.tree[,3])
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
      S=Scale.dependence.S.single(a,b,scale,MI)
      Levins_Simpson=Scale.dependence.Simpson(a,b,scale)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Scale_dependence_S=S$Scale_dependence_S
      Scale_dependence_Wei_S=Scale_dependence_S*(0.25/Levins_Simpson)
      #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
      outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S,Levins_Simpson=Levins_Simpson,Scale_dependence_Wei_S=Scale_dependence_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Scale.dependence.Wei_S.single(a[j,],b,scale,MI)$Scale_dependence_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Scale_dependence_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Scale_dependence_Wei_S
  data=Scale.dependence.Wei_S.mult(basexy,b,scale,MI)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Scale_dependence_Wei_S")
  vgm1 <- variogram(Scale_dependence_Wei_S~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Scale_dependence_Wei_S~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title="The Scale dependence Wei_S")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
####计算加权进界邻体荫蔽度的修正系数function~
#MI.Tan.S=function(minx,maxx,miny,maxy,b,seq,tan,MI)
{
  library(sp)
  library(gstat)
  library(tcltk)
  #################
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  a= dplyr::distinct(basexy)####删除重复样点
  ################
  Tan.dependence.S.mult=function(a,b,tan)
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

        if(is.nan(Tan_dependence_S)==T)(Tan_dependence_S=100000000)
        outcome=list(a=a,Neighbourhood=Neighbourhood,Tan_dependence_S=Tan_dependence_S)
        outcome
      }


      #####计算尺度依赖加权荫蔽度
      S=Tan.dependence.S.single(a,b,tan)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Tan_dependence_S=S$Tan_dependence_S

      Tan_dependence_Wei_S=Tan_dependence_S
      outcome=list(a=a,Neighbourhood=Neighbourhood,Tan_dependence_S=Tan_dependence_S,Levins_Simpson=Levins_Simpson,Tan_dependence_Wei_S=Tan_dependence_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Tan.dependence.Wei_S.single(a[j,],b,tan)$Tan_dependence_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Tan_dependence_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  data=Tan.dependence.S.mult(a,b,tan)
  # data=subset(data,data[,3]>0)
  data=data[,3]
  n_=length(data)
  x=matrix()
  data=data[order(data)]
  for(i in 1:n_){x[i]=(atan(data[i]/(MI*pi))/pi*2)}
  data=1:n_
  lm=lm(data~x)
  print(plot(x,data)+abline(lm))
  print(summary(lm))
}
####计算单一参照树的加权进界邻体荫蔽度function~
Tan.dependence.Wei_S.single=function(a,b,tan,MI)
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
      d=subset(d,d>0)
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
      d=subset(d,d>0)
      d$tangent=(d$size-a[,3])/d$d
      d=subset(d,tangent>tan)
      colnames(d) = c("x","y","Size","Distance","tangent")
      d
    }
    Nei.tree=Neighbourhood.single1(a,b,tan)
    Nei.tree[,3]=Nei.tree[,3]-a[,3]
    Nei.tree=subset(Nei.tree,Size>0)
    Nei.tree=subset(Nei.tree,Distance>0)
    Nei.tree
    Simpn=sum(Nei.tree[,3])
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
####绘制单一参照树的加权进界邻体荫蔽度图function~
plot.Tan.dependence.Wei_S.single=function(a,b,tan,MI)
{
  library(ggplot2)
  Tan.dependence.Wei_S.single=function(a,b,tan,MI)
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
        d=subset(d,d>0)
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
        d=subset(d,d>0)
        d$tangent=(d$size-a[,3])/d$d
        d=subset(d,tangent>tan)
        colnames(d) = c("x","y","Size","Distance","tangent")
        d
      }
      Nei.tree=Neighbourhood.single1(a,b,tan)
      Nei.tree[,3]=Nei.tree[,3]-a[,3]
      Nei.tree=subset(Nei.tree,Size>0)
      Nei.tree=subset(Nei.tree,Distance>0)
      Nei.tree
      Simpn=sum(Nei.tree[,3])
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
  Ne= Tan.dependence.Wei_S.single(a,b,tan,MI)
  big=subset(Ne$Neighbourhood,Size>a[,3])
  small=subset(Ne$Neighbourhood,Size<=a[,3])
  n=nrow(big)
  if(n==0)
  {
    center=cbind(a[,1:3],1)
    colnames(center)=c("x","y","size","group")
    unitg=center
    Neigh=subset(unitg,size<0)
    colnames(Neigh)=c("x","y","size","group")
    colnames(small)=c("x","y","size")
    max=10
  }
  if(n!=0)
  {
    center=cbind(cbind(rep(Ne$a[,1],each=n),rep(Ne$a[,2],each=n),rep(Ne$a[,3],each=n)),c(1:n))
    center=as.data.frame(center)
    colnames(center)=c("x","y","size","group")
    Neigh=cbind(big,c(1:n))
    Neigh=Neigh[,c(1,2,3,6)]
    colnames(Neigh)=c("x","y","size","group")
    colnames(small)=c("x","y","size")
    unitg=rbind(center,Neigh)
    max=max(abs(Ne$Neighbourhood[,1]-Ne$a[,1]),abs(Ne$Neighbourhood[,2]-Ne$a[,2]))
  }
  p=ggplot()
  p=p+geom_hline(data=Ne$a,aes(yintercept=y),linetype=5,col="red")+geom_vline(data=Ne$a,aes(xintercept=x),linetype=5,col="red")
  p=p+geom_line(data=unitg,aes(x=x,y=y,group=group),linetype=2,size=1)
  p=p+geom_point(data=Ne$a,aes(x,y),size=13/max(unitg[,3])*a[,3]+2,color="red4",alpha=0.7)+geom_point(data=Neigh,aes(x,y),size=13/max(unitg[,3])*Neigh[,3]+2,color="green4",alpha=0.7)
  p=p+geom_point(data=small,aes(x,y),size=13/max(unitg[,3])*small[,3]+2,color="grey",alpha=0.7)
  p=p+geom_point(data=Ne$a,aes(x,y),size=2)+geom_point(data=Neigh,aes(x,y),size=2)+geom_point(data=small,aes(x,y),size=2)
  p=p+lims(x=c(Ne$a[,1]-1.1*max,Ne$a[,1]+1.1*max),y=c(Ne$a[,2]-1.1*max,Ne$a[,2]+1.1*max))
  p=p+ annotate("text",x=Ne$a[,1]-max+0.3*(max),y=Ne$a[,2]+max-0.125*(max),label = paste0("Tan_dependence_Wei_S=",round(Ne$Tan_dependence_Wei_S,3)))
  p=p+theme_bw()
  suppressMessages(print(p))
}
####绘制加权进界邻体荫蔽度克里金图function~
plot.Scale.dependence.Wei_S.Krig=function(minx,maxx,miny,maxy,b,seq,scale,MI)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Scale.dependence.Wei_S.mult=function(a,b,scale,MI)
  {
    library(tcltk)
    Scale.dependence.Wei_S.single=function(a,b,scale,MI)
    {
      ######加载尺度依赖荫蔽度函数
      Scale.dependence.S.single=function(a,b,scale,MI)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3
          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          #d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        n=nrow(Nei.tree)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        n.dif=subset(Nei.dif.high,Nei.dif.high$Size_dif>0)
        n.dif=nrow(n.dif)
        if(n.dif==0)
        {Scale_dependence_S=0
        Neighbourhood=Nei.dif.high}
        if(n.dif!=0)
        {
          Neighbourhood=Nei.dif.high
          Neighbourhood1=Neighbourhood
          Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
          Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
          Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
          S1=Neighbourhood1[,5]/Neighbourhood1[,4]
          S2=as.data.frame(S1)
          S=sum(S2,na.rm=T)
          Scale_dependence_S=S
        }
        Scale_dependence_S=(atan( Scale_dependence_S/(MI*pi))/pi*2)
        outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S)
        outcome
      }
      Scale.dependence.Simpson=function(a,b,scale)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        Nei.tree
        Nei.tree[,3]=Nei.tree[,3]-a[,3]
        Nei.tree=subset(Nei.tree,Size>0)
        Nei.tree=subset(Nei.tree,Distance>0)
        Nei.tree
        Simpn=sum(Nei.tree[,3])
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
      S=Scale.dependence.S.single(a,b,scale,MI)
      Levins_Simpson=Scale.dependence.Simpson(a,b,scale)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Scale_dependence_S=S$Scale_dependence_S
      Scale_dependence_Wei_S=Scale_dependence_S*(0.25/Levins_Simpson)
      #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
      outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S,Levins_Simpson=Levins_Simpson,Scale_dependence_Wei_S=Scale_dependence_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Scale.dependence.Wei_S.single(a[j,],b,scale,MI)$Scale_dependence_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Scale_dependence_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Scale_dependence_Wei_S
  data=Scale.dependence.Wei_S.mult(basexy,b,scale,MI)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Scale_dependence_Wei_S")
  vgm1 <- variogram(Scale_dependence_Wei_S~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Scale_dependence_Wei_S~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title="The Scale dependence Wei_S")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
plot.Tan.dependence.Wei_S.Krig=function(minx,maxx,miny,maxy,b,seq,tan,MI)
{
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
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Tan.dependence.Wei_S.single(a[j,],b,tan)$Tan_dependence_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Tan_dependence_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Scale_dependence_Wei_S
  data=Tan.dependence.Wei_S.mult(basexy,b,tan)
  #data=na.omit(data)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Tan_dependence_Wei_S")
  vgm1 <- variogram(Tan_dependence_Wei_S~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  sd=subset(vgm1$dist,vgm1$dist<m$range[2])
  spre=m$psill[1]+(m$psill[2])*((3*sd)/(2*m$range[2])-sd^3/(2*(m$range[2])^3))
  bd=subset(vgm1$dist,vgm1$dist>m$range[2])
  bpre=rep((m$psill[1]+m$psill[2]),length(bd))
  pre=rbind(as.matrix(spre),as.matrix(bpre))
  Coefficient_of_Determination=1- sum((pre-vgm1$gamma)^2)/sum((vgm1$gamma-mean(vgm1$gamma))^2)
  print(paste("Coefficient_of_Determination=",Coefficient_of_Determination))
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Tan_dependence_Wei_S~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours = terrain.colors(10))
  p2=p2+labs(title="The Tan dependence Wei_S")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
plot.Scale.dependence.Wei_S.Krig=function(minx,maxx,miny,maxy,b,seq,scale,MI)
{
  #### Krig准备工作开始, 导入程序包
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Scale.dependence.Wei_S.mult=function(a,b,scale,MI)
  {
    library(tcltk)
    Scale.dependence.Wei_S.single=function(a,b,scale,MI)
    {
      ######加载尺度依赖荫蔽度函数
      Scale.dependence.S.single=function(a,b,scale,MI)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3
          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          #d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","Y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        n=nrow(Nei.tree)
        Nei.dif.high=cbind(Nei.tree,Nei.tree[,3]-a[,3])
        colnames(Nei.dif.high) = c("x","Y","Size","Distance","Size_dif ")
        n.dif=subset(Nei.dif.high,Nei.dif.high$Size_dif>0)
        n.dif=nrow(n.dif)
        if(n.dif==0)
        {Scale_dependence_S=0
        Neighbourhood=Nei.dif.high}
        if(n.dif!=0)
        {
          Neighbourhood=Nei.dif.high
          Neighbourhood1=Neighbourhood
          Neighbourhood[which(Neighbourhood[,5]<=0),5]=NA
          Neighbourhood1[which(Neighbourhood1[,5]<=0),5]=0
          Neighbourhood1[which(Neighbourhood1[,4]==0),4]=0.00000001
          S1=Neighbourhood1[,5]/Neighbourhood1[,4]
          S2=as.data.frame(S1)
          S=sum(S2,na.rm=T)
          Scale_dependence_S=S
        }
        Scale_dependence_S=(atan( Scale_dependence_S/(MI*pi))/pi*2)
        outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S)
        outcome
      }
      Scale.dependence.Simpson=function(a,b,scale)
      {
        Neighbourhood.single1=function(a,b,scale)###找出一定尺度内的邻体
        {#3

          c=b[,1:2]
          for (i in 1:nrow(b))
          {#4
            c[i,]=(b[i,1:2]-a[1,1:2])^2
            d=(c[,1]+c[,2])^(1/2)
          }#4
          d=cbind(b,d)
          d=subset(d,d>0)
          d=subset(d,d<scale)
          colnames(d) = c("x","y","Size","Distance")
          d
        }
        Nei.tree=Neighbourhood.single1(a,b,scale)
        Nei.tree
        Nei.tree[,3]=Nei.tree[,3]-a[,3]
        Nei.tree=subset(Nei.tree,Size>0)
        Nei.tree=subset(Nei.tree,Distance>0)
        Nei.tree
        Simpn=sum(Nei.tree[,3])
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
      S=Scale.dependence.S.single(a,b,scale,MI)
      Levins_Simpson=Scale.dependence.Simpson(a,b,scale)
      a=S$a
      Neighbourhood=S$Neighbourhood
      Scale_dependence_S=S$Scale_dependence_S
      Scale_dependence_Wei_S=Scale_dependence_S*(0.25/Levins_Simpson)
      #Scale_dependence_Wei_S=Scale_dependence_Wei_S/(Scale_dependence_Wei_S+1)
      outcome=list(a=a,Neighbourhood=Neighbourhood,Scale_dependence_S=Scale_dependence_S,Levins_Simpson=Levins_Simpson,Scale_dependence_Wei_S=Scale_dependence_Wei_S)
      outcome
    }
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("进度","已完成 %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,1:2]),as.matrix(Scale.dependence.Wei_S.single(a[j,],b,scale,MI)$Scale_dependence_Wei_S))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
    }
    end_time =Sys.time()  ## 记录程序结束时间
    close(pb)
    run_time = end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","Scale_dependence_Wei_S")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)

  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)

  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)

  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)

  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  basexy[,3]=0
  colnames(basexy) <- c("x", "y","Size")
  basexy= dplyr::distinct(basexy)####删除重复样点
  ######计算插值点的Scale_dependence_Wei_S
  data=Scale.dependence.Wei_S.mult(basexy,b,scale,MI)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标
  spplot(data,"Scale_dependence_Wei_S")
  vgm1 <- variogram(Scale_dependence_Wei_S~1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  sd=subset(vgm1$dist,vgm1$dist<m$range[2])
  spre=m$psill[1]+(m$psill[2])*((3*sd)/(2*m$range[2])-sd^3/(2*(m$range[2])^3))
  bd=subset(vgm1$dist,vgm1$dist>m$range[2])
  bpre=rep((m$psill[1]+m$psill[2]),length(bd))
  pre=rbind(as.matrix(spre),as.matrix(bpre))
  Coefficient_of_Determination=1- sum((pre-vgm1$gamma)^2)/sum((vgm1$gamma-mean(vgm1$gamma))^2)
  p1=plot(vgm1, model=m)
  print(p1)
  krige_res <- krige(Scale_dependence_Wei_S~1, data, basexy1, model = m)

  ### 查看克里格插值的结果

  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  Basexyz[which(Basexyz$Vaule<0),3]=0
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(limits=c(0,1),colours = terrain.colors(10))
  p2=p2+labs(title="The Scale dependence Wei_S")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}

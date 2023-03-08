boundary.eff=function (a, minx, maxx, miny, maxy, frame) 
{
  a = subset(a, x > minx + frame & x < maxx - frame)
  a = subset(a, y > miny + frame & y < maxy - frame)
  a
}
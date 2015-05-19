system.time(mod<-glmer(data=sdat[[2]][,],formula = as.factor(P_A) ~ 1 + poly(Phylo.Relatedness,2,raw=TRUE) + (1 + poly(Phylo.Relatedness,2,raw=TRUE)|Species), family = "binomial"))
wald<-confint.merMod(mod,method="Wald")
system.time(prof<-confint.merMod(mod,method="profile",parm=3,.progress="txt"))
system.time(boo<-confint(mod,method="boot"))
system.time(boo2<-confint(mod,method="boot",nsim=1000))

m<-melt(boo2)
m2<-split(m,m$Var2)
d<-melt(lapply(m2,function(x){
  y<-trajF(x[1,3],x[2,3],x[3,3],sdat[[2]][,]$Phylo.Relatedness)
  x<-sdat[[2]][,]$Phylo.Relatedness
  d<-data.frame(x,y)
}),id.vars=c("x","y"))



ggplot(d,aes(x=x,y=y,col=L1)) + geom_line()
30177/60/60


system.time(modg<-glm(data=sdat[[2]][,],formula = as.factor(P_A) ~ poly(Phylo.Relatedness,2,raw=TRUE), family = "binomial"))
m<-coefficients(modg)
y<-trajF(m[[1]],m[[2]],m[[3]],sdat[[2]][,]$Phylo.Relatedness)
x<-sdat[[2]][,]$Phylo.Relatedness
d<-data.frame(x,y)

ggplot(d,aes(x=x,y=y)) + geom_line()

confint(modg)

specSdf$param<-as.factor(str_match(specSdf$parameter,"\\w+"))
gdf<-group_by(specSdf,Species,param) %>% 
  summarize(Lower=quantile(estimate,0.025),Upper=quantile(estimate,0.975),mean=mean(estimate))

min(gdf$mean)  

vp<-function(intercept,linear,polynomial,x=seq(0,1,.01)){

#define trajectory function
trajF<-function(intercept,linear,polynomial,x){
  p<-inv.logit(intercept + linear * x  + polynomial * x^2)
  return(p)
}

y<-trajF(intercept,linear,polynomial,x)
d<-data.frame(x=x,y=y)
}
d<-vp(-.66,6.7,-12.2)
ggplot(d,aes(x=x,y=y)) + geom_line()

intercept=-1
linear=10
polynomial=-7

pol<-seq(5,7,.1)
d<-lapply(pol,function(x){vp(-1,x,-18)})
names(d)<-pol
d<-melt(d,id.vars=c("x","y"))
ggplot(d,aes(x=x,y=y,col=as.numeric(L1))) + geom_line(aes(group=L1)) + scale_color_continuous(low="blue",high="red")


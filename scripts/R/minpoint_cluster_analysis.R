library(tidyverse)

dir_path="data/cluster_analysis/"

files=dir(dir_path,full.names = TRUE)

df=data.frame(NULL)
for (f in files){
  df<-rbind(df,read.csv(f))
  
}

df<-df %>% 
  group_by(annotation) %>% 
  mutate(clusters=clusters/max(clusters))

lm.fit=lm(log10(clusters)~log10(minpoints),data=df)
lm.predict=predict(lm.fit,data.frame(minpoints=seq(min(df$minpoints),max(df$minpoints))))

df.predict=data.frame(minpoints=seq(min(df$minpoints),max(df$minpoints)),clusters=10^lm.predict)
df.d2dx2<-df.predict %>% 
  distinct(minpoints,clusters) %>% 
  arrange(minpoints)
n=nrow(df.d2dx2)
d2dx2=(df.d2dx2$clusters[2:n]-df.d2dx2$clusters[1:(n-1)])/df.d2dx2$clusters[1:(n-1)]*100

indx<-which(d2dx2>-5)[1]
minpoint_opt<-df.predict$minpoints[indx-1]

ggplot()+
  geom_point(data=df,aes(x=log10(minpoints),y=log10(clusters)),alpha=0.1)+
  geom_smooth(data=df,aes(x=log10(minpoints),y=log10(clusters)),method="lm") +
  geom_line(data=data.frame(minpoints=c(log10(minpoint_opt),log10(minpoint_opt)),clusters=c(log10(min(df$clusters)),log10(max(df$clusters)))),aes(x=minpoints,y=clusters),color="red") 


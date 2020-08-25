library(tidyverse)

df=read.csv("data/mds_cluster_dataframes/ATP_synthase_subunit_delta.csv")

df_good<-df %>% filter(Cluster > -1)
df_noise<-df %>% filter(Cluster == -1)
ggplot()+
  geom_point(data=df_good,aes(x=Component_1,y=Component_2,color=as.factor(Cluster))) +
  geom_point(data=df_noise,aes(x=Component_1,y=Component_2),color='black')

# library("bios2mds")
# yll_align <- import.fasta("/home/acwach/Homology/yll_align.fas")
# 
# yll_align<-c(yll_align[-c(35,53,63)],yll_align[c(35,53,63)])
# 
# yll_dist<-mat.dis(yll_align, yll_align)
# 
# rnames<-paste(colnames(yll_dist),runif(length(colnames(yll_dist))),sep="_")
# colnames(yll_dist)<-rnames
# rownames(yll_dist)<-rnames
# 
# #####2D MDS
# yll_mmds<-mmds(yll_dist, pc = 2, group.file = NULL)
# 
# km1<-kmeans.run(yll_mmds$coord, nb.clus = 3)
# clus1<-substr(names(unlist(km1[[2]][1])),3,1000)
# clus2<-substr(names(unlist(km1[[2]][2])),3,1000)
# clus3<-substr(names(unlist(km1[[2]][3])),3,1000)
# 
# yll058w_colors<-rep("black",length(rnames))
# yll058w_colors[which(rnames %in% clus1)]<-"Cluster 1"
# yll058w_colors[which(rnames %in% clus2)]<-"Cluster 2"
# yll058w_colors[which(rnames %in% clus3)]<-"Cluster 3"
# 
# yll_taxa<-array()
# yll_taxa<-rep("Ascomycota",length(rnames))
# yll_taxa[c(3,10,13,103)]<-"Bacteria"
# yll_taxa[c(109,113,114,287,288,289,297,298,299)]<-"Fungal outgroups"
# 
# df_yll<-data.frame(PC1=yll_mmds$coord$PC1,PC2=yll_mmds$coord$PC2,km_clust=yll058w_colors,spp_name=rnames,taxa=yll_taxa)
# df_yll$km_clust<-factor(df_yll$km_clust,levels=c("Cluster 1","Cluster 2","Cluster 3"))#,"Cluster 4"))
# 
# 
# require("ggplot2")
# png("yll_mds_ggplot2d.png",height=2000,width=3000,res=300)
# 
# fig5b <- ggplot(df_yll)+
# 	theme_linedraw()+
# 	geom_point(aes(x=PC1,y=PC2,col=km_clust,shape=yll_taxa), size = 2)+
# 	scale_color_manual(values=c("#E64A19","#303F9F","#7B1FA2"))+
# 	scale_shape_manual(values=c(20, 0, 2))+
# 	labs(x="PC1",y="PC2")+
# 	annotate("segment", x = df_yll[331,]$PC1, xend = df_yll[331,]$PC1, y = df_yll[331,]$PC2, yend = df_yll[331,]$PC2+.7)+
# 	annotate(geom="text", x=df_yll[331,]$PC1-.5, y=df_yll[331,]$PC2+.9, label="YJR130C", color="black",size = 2)+
# 	annotate("segment", x = df_yll[330,]$PC1, xend = df_yll[330,]$PC1-.5, y = df_yll[330,]$PC2, yend = df_yll[330,]$PC2+.5)+
# 	annotate(geom="text", x=df_yll[330,]$PC1-.7, y=df_yll[330,]$PC2+.5, label="YML082W", color="black",size = 2)+
# 	annotate("segment", x = df_yll[329,]$PC1, xend = df_yll[329,]$PC1+.5, y = df_yll[329,]$PC2, yend = df_yll[329,]$PC2)+
# 	annotate(geom="text", x=df_yll[329,]$PC1+1.1, y=df_yll[329,]$PC2, label="YLL058W", color="black",size = 2)+
#   theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_blank(),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm")) +
#   guides(shape = guide_legend(nrow = 3, byrow=F, order = 1,
#                               override.aes=list(size = 2)),
#          color = guide_legend(nrow = 3, byrow=F, order = 2,
#                               override.aes=list(size = 2)))
# 
# cluster.legend <- g_legend(fig5b + guides(shape = F, color = guide_legend(nrow = 1, byrow=F, override.aes=list(size = 2))))
# 
# dev.off()
# 
# 
# ##making tree
# 
# library("ggtree")
# 
# tree0<-read.tree(text='(Lipomycetaceae,(Trigonopsidaceae,(Dipodascaceae,
# ((CUG-Ser1 clade,Pichiaceae),((Saccharomycopsis,Ascoidea),((Cyberlindnera,Barnettozyma),(Hanseniaspora,((Lachancea,(Eremothecium,Kluyveromyces)),((Torulaspora,Zygotorulaspora),(Tetrapisispora,((Naumovozyma,Kazachstania),(Nakaseomyces,Saccharomyces))))))))))));
# ')
# tree1<-read.tree(text='(Lipomycetaceae,(Trigonopsidaceae,(Dipodascaceae,
# (((Cephaloascus,((Meyerozyma,Yamadazyma),((Metschnikowia,Hyphopichia),(Debaryomyces,(Priceomyces,((Spathaspora,Scheffersomyces),(Teunomyces,Suhomyces) ) )))))
# ,Pichiaceae),((Saccharomycopsis,Ascoidea),((Cyberlindnera,Barnettozyma),(Hanseniaspora,((Lachancea,(Eremothecium,Kluyveromyces)),((Torulaspora,Zygotorulaspora),(Tetrapisispora,((Naumovozyma,Kazachstania),(Nakaseomyces,Saccharomyces))))))))))));
# ')
# 
# cluster2_present<-c("Lachancea","Zygotorulaspora","Torulaspora","Saccharomyces")
# cluster1_present<-c("Lachancea","Tetrapisispora","Nakaseomyces","Kazachstania","Hanseniaspora","Torulaspora","Zygotorulaspora","Kluyveromyces","Eremothecium","Naumovozyma","Saccharomyces")
# cluster3_present<-c("Trigonopsidaceae","Lipomycetaceae","Dipodascaceae","Pichiaceae","Saccharomycopsis","Debaryomyces","Yamadazyma",
# "Sporopachydermia","Cyberlindnera","Alloascoidea","Starmera","Nakazawaea","Cephaloascus","Priceomyces","Barnettozyma","Metschnikowia","Peterozyma",
# "Cephaloascus","Meyerozyma","Spathaspora","Scheffersomyces","Hyphopichia")
#  
# clus1colors<-rep("white",30)
# clus1colors[c(1,2,3,4,5,6,7,8,9,10,11)]<-"#E64A19"
# 
# clus2colors<-rep("white",30)
# clus2colors[c(1,6,7,10)]<-"#303F9F"
# 
# clus3colors<-rep("white",30)
# clus3colors[c(12,13,15,18,19,20,21,22,23,24,25,26,27,28,29,30)]<-"#7B1FA2"
# 
# png("tree.png",height=2000,width=2000,res=300)
# fig5c <- ggtree(tree1)+
#   geom_tiplab(size = 2.2)+
#   xlim(0, 20)+
#   annotate("point",x=17,y=(1:30)-.1,col=rev(clus1colors)) +
#   annotate("point",x=18,y=(1:30)-.1,col=rev(clus2colors)) +
#   annotate("point",x=19,y=(1:30)-.1,col=rev(clus3colors)) 
# 
# #geom_point(x=30,y=(1:20)-.1)
# dev.off()

##### NEW
library(bios2mds)
library(ggplot2)
library(ggtree)

yll_align<-import.fasta("/home/acwach/Homology/yll_align.fas")

yll_align<-c(yll_align[-c(35,53,63)],yll_align[c(35,53,63)])

yll_dist<-mat.dis(yll_align, yll_align)

rnames<-paste(colnames(yll_dist),runif(length(colnames(yll_dist))),sep="_")
colnames(yll_dist)<-rnames
rownames(yll_dist)<-rnames

#####2D MDS
yll_mmds<-mmds(yll_dist, pc = 2, group.file = NULL)

km1<-kmeans.run(yll_mmds$coord, nb.clus = 3)
clus1<-substr(names(unlist(km1[[2]][1])),3,1000)
clus2<-substr(names(unlist(km1[[2]][2])),3,1000)
clus3<-substr(names(unlist(km1[[2]][3])),3,1000)

yll058w_colors<-rep("black",length(rnames))
yll058w_colors[which(rnames %in% clus1)]<-"Cluster 1"
yll058w_colors[which(rnames %in% clus2)]<-"Cluster 2"
yll058w_colors[which(rnames %in% clus3)]<-"Cluster 3"

yll_taxa<-array()
yll_taxa<-rep("Ascomycota",length(rnames))
yll_taxa[c(3,10,13,103)]<-"Bacteria"
yll_taxa[c(109,113,114,287,288,289,297,298,299)]<-"Fungal outgroups"

df_yll<-data.frame(PC1=yll_mmds$coord$PC1,PC2=yll_mmds$coord$PC2,km_clust=yll058w_colors,spp_name=rnames,taxa=yll_taxa)
df_yll$km_clust<-factor(df_yll$km_clust,levels=c("Cluster 1","Cluster 2","Cluster 3"))#,"Cluster 4"))

unique(df_yll$km_clust)

fig.mds <- ggplot(df_yll)+
  geom_point(aes(x=PC1,y=PC2,col=km_clust,shape=yll_taxa),
             size = 2)+
  theme(legend.title=element_blank())+
  scale_color_manual(values=c("Cluster 1"="#006666",
                              "Cluster 2"="#CC33FF",
                              "Cluster 3"="#99CCFF"),
                     labels = c('Cluster 1' = 'Ancestral class',
                                'Cluster 2' = '*YLL058W* class',
                                'Cluster 3' = '*STR2/YML082W* class'))+ #
  scale_shape_manual(values=c(20, 0, 2))+
  labs(x="PC1",y="PC2")+
  annotate("segment", x = df_yll[331,]$PC1, xend = df_yll[331,]$PC1, y = df_yll[331,]$PC2, yend = df_yll[331,]$PC2+.7)+
  annotate(geom="text", x=df_yll[331,]$PC1, y=df_yll[331,]$PC2+.9, label="STR2", color="black", size = 2,fontface = "italic")+
  annotate("segment", x = df_yll[330,]$PC1, xend = df_yll[330,]$PC1-.5, y = df_yll[330,]$PC2, yend = df_yll[330,]$PC2+.5)+
  annotate(geom="text", x=df_yll[330,]$PC1-.5, y=df_yll[330,]$PC2+.7, label="YML082W", color="black", size = 2,fontface = "italic")+
  annotate("segment", x = df_yll[329,]$PC1, xend = df_yll[329,]$PC1+.5, y = df_yll[329,]$PC2, yend = df_yll[329,]$PC2)+
  annotate(geom="text", x=df_yll[329,]$PC1+0.8, y=df_yll[329,]$PC2, label="YLL058W", color="black", size = 2,fontface = "italic") +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

ggsave(sprintf("%s/FigureYLLMDS.jpg",fig_path), fig.mds,
       height = one.5c, width = one.5c, units = 'mm',
       bg = 'white',
       dpi = 600)

##making tree
tree0<-read.tree(text='(Lipomycetaceae,(Trigonopsidaceae,(Dipodascaceae,
((CUG-Ser1 clade,Pichiaceae),((Saccharomycopsis,Ascoidea),((Cyberlindnera,Barnettozyma),(Hanseniaspora,((Lachancea,(Eremothecium,Kluyveromyces)),((Torulaspora,Zygotorulaspora),(Tetrapisispora,((Naumovozyma,Kazachstania),(Nakaseomyces,Saccharomyces))))))))))));
')
tree1<-read.tree(text='(Lipomycetaceae,(Trigonopsidaceae,(Dipodascaceae,
(((Cephaloascus,((Meyerozyma,Yamadazyma),((Metschnikowia,Hyphopichia),(Debaryomyces,(Priceomyces,((Spathaspora,Scheffersomyces),(Teunomyces,Suhomyces) ) )))))
,Pichiaceae),((Saccharomycopsis,Ascoidea),((Cyberlindnera,Barnettozyma),(Hanseniaspora,((Lachancea,(Eremothecium,Kluyveromyces)),((Torulaspora,Zygotorulaspora),(Tetrapisispora,((Naumovozyma,Kazachstania),(Nakaseomyces,Saccharomyces))))))))))));
')

cluster2_present<-c("Lachancea","Zygotorulaspora","Torulaspora","Saccharomyces")
cluster1_present<-c("Lachancea","Tetrapisispora","Nakaseomyces","Kazachstania","Hanseniaspora","Torulaspora","Zygotorulaspora","Kluyveromyces","Eremothecium","Naumovozyma","Saccharomyces")
cluster3_present<-c("Trigonopsidaceae","Lipomycetaceae","Dipodascaceae","Pichiaceae","Saccharomycopsis","Debaryomyces","Yamadazyma",
                    "Sporopachydermia","Cyberlindnera","Alloascoidea","Starmera","Nakazawaea","Cephaloascus","Priceomyces","Barnettozyma","Metschnikowia","Peterozyma",
                    "Cephaloascus","Meyerozyma","Spathaspora","Scheffersomyces","Hyphopichia")

clus1colors<-rep("white",30)
clus1colors[c(1,2,3,4,5,6,7,8,9,10,11)]<-"#006666"
clus1colors2<-rep("white",30)
clus1colors2[c(1,2,3,4,5,6,7,8,9,10,11)]<-"#9E9E9E"

clus2colors<-rep("white",30)
clus2colors[c(1,6,7,10)]<-"#CC33FF"
clus2colors2<-rep("white",30)
clus2colors2[c(1,6,7,10)]<-"#212121"

clus3colors<-rep("white",30)
clus3colors[c(12,13,15,18,19,20,21,22,23,24,25,26,27,28,29,30)]<-"#99CCFF"
clus3colors2<-rep("white",30)
clus3colors2[c(12,13,15,18,19,20,21,22,23,24,25,26,27,28,29,30)]<-"#9E9E9E"

fig.yll.tree <- ggtree(tree1)+
  geom_tiplab(size = 2.8)+
  xlim(-9, 36-2)+
  annotate("rect",xmin=29.4-2.2,xmax=30.6-1.8,ymin=(1:30)-.6,ymax=(1:30)+.4,fill=rev(clus1colors2),col='white') +
  annotate("rect",xmin=31.4-2.2,xmax=32.6-1.8,ymin=(1:30)-.6,ymax=(1:30)+.4,fill=rev(clus2colors2),col='white') +
  annotate("rect",xmin=33.4-2.2,xmax=34.6-1.8,ymin=(1:30)-.6,ymax=(1:30)+.4,fill=rev(clus3colors2),col='white') +
  annotate("point",x=30-2,y=(1:30)-.1,col=rev(clus1colors),size=1) +
  annotate("point",x=32-2,y=(1:30)-.1,col=rev(clus2colors),size=1) +
  annotate("point",x=34-2,y=(1:30)-.1,col=rev(clus3colors),size=1) #+
  # annotate("text",x = -5, y = 25,
  #          label = "atop(atop('Whole genome duplication', 'creates paralogs '*italic('STR2')), 'and '*italic('YML082W'))",
  #          size = 3, parse = T)
#geom_point(x=30,y=(1:20)-.1)
# atop(textstyle('Whole genome duplication'), textsyle('creates paralogs')*italic('STR2'))
# 
# ggsave(sprintf("%s/FigureYLLTree.jpg",fig_path), fig.yll.tree,
#        height = one.5c, width = 200, units = 'mm',
#        bg = 'white',
#        dpi = 600)


fig.yll.tree.leg <- data.frame(levels = c('a','b','c','d','e')) %>%
  ggplot(aes(x = 1, y = levels, col = levels, shape = levels)) +
  geom_point() +
  scale_color_manual(name = 'Classes of *YLL058W* homologs',
                     values = c('a' = '#99CCFF',
                                'b' = '#CC33FF',
                                'c' = '#006666',
                                'd' = '#212121',
                                'e' = '#9E9E9E'),
                     labels = c('a' = 'Ancestral class',
                                'b' = '*YLL058W* class',
                                'c' = '*STR2/YML082W* class',
                                'd' = 'In sulfur-related gene cluster',
                                'e' = 'Not in sulfur gene cluster')) +
  scale_shape_manual(name = 'Classes of *YLL058W* homologs',
                     values = c('a' = 16,
                                'b' = 16,
                                'c' = 16,
                                'd' = 15,
                                'e' = 15),
                     labels = c('a' = 'Ancestral class',
                                'b' = '*YLL058W* class',
                                'c' = '*STR2/YML082W* class',
                                'd' = 'In sulfur-related gene cluster',
                                'e' = 'Not in sulfur gene cluster')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title =  ggtext::element_markdown(size = titles),
        legend.text =  ggtext::element_markdown(size = txt),
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(4, "mm")) +
  guides(col = guide_legend(nrow=3, override.aes=list(size = 2)))

fig.yll.tree.leg <- g_legend(fig.yll.tree.leg)

fig5b <- plot_grid(fig.yll.tree, fig.yll.tree.leg,
          ncol = 1, rel_heights = c(5,1))


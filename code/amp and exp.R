#fly to mouse
gene_list <- read.csv("E:/GO_layers/blast/restart/ensembl/mouse_fly.csv",header=TRUE)
mouse_gene <- read.csv("E:/GO_layers/blast/circidian_gene/mouse_gene.csv")
mouse_gene <- mouse_gene[,1]
fly_gene <- read.csv("E:/GO_layers/blast/circidian_gene/fly_gene.csv")
fly_gene <- fly_gene[,1]
fly2mouse <- unique(gene_list[gene_list[,2] %in% fly_gene, 1])
flynot2mouse <- unique(gene_list[gene_list[,2] %in% fly_gene, 2])
#mouse to fly
mouse2fly <- gene_list[gene_list[,1] %in% mouse_gene, 2]

#mouse to fly and fly to mouse overlap
mouse2fly2mouse <- gene_list[gene_list[,2] %in% mouse2fly, 1]
fly_orthology_gene <- unique(gene_list[gene_list[,2] %in% mouse2fly, 2])


length(intersect(mouse2fly2mouse, fly2mouse))
#orthology fly2mouse overlap  not tissue
length(intersect(fly2mouse, mouse_gene))
temp_b <- draw.pairwise.venn(area1=6393,area2=6267,cross.area=1563,category = c("mouse", "fly"),filename=NULL,  scaled=FALSE,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')


#orthology fly2mouse overlap  tissue average
# mouse_gene <- read.csv("E:/GO_layers/photo/result/mouse/mouse_gene.csv")
# fly_gene <- read.csv("E:/GO_layers/circadian_data/fly/fly_gene.csv")

# 
# gene_id <- list()
# for(i in 1:ncol(mouse_gene))
# {
#   gene_id[[i]] <- bitr(mouse_gene[,i], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL
# }
# temp <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
# colnames(temp) <- mouse_tissue
# write.csv(temp, "E:/GO_layers/photo/result/mouse/mouse_genesymbol.csv")
mouse_gene <- read.csv("E:/GO_layers/photo/result/mouse/mouse_genesymbol.csv")
fly_gene <- read.csv("E:/GO_layers/circadian_data/fly/fly_gene.csv")
mouse_tissue <- colnames(mouse_gene)
fly_tissue <- colnames(fly_gene)
# gene_id <- list()
# for(i in 1:ncol(fly_gene))
# {
#   gene_id[[i]] <- as.character(gene_list[gene_list[,2] %in% fly_gene[,i], 1])
# }
# temp <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
# write.csv(temp, "E:/GO_layers/circadian_data/fly/fly2mouse_gene.csv")
fly2mouse_gene <- read.csv("E:/GO_layers/circadian_data/fly/fly2mouse_gene.csv")
common_term_num <- 0
total_num <- 0
data1_term_num <- 0
data2_term_num <- 0
not_fly_gene <- 0
plot_list = list()
for(i in 1:length(mouse_tissue))
{
  current_mouse_gene <- unique(mouse_gene[,i])
  for(j in 1:length(fly_tissue))
  {
    current_fly_gene <- na.omit(fly_gene[,j])
    total_num <- total_num + 1
    orth_fly_gene <- unique(gene_list[gene_list[,2] %in% current_fly_gene, 1])
    not_fly_gene <- not_fly_gene + length(setdiff(current_fly_gene, gene_list[,2]))
    temp_b <- draw.pairwise.venn(area1=length(current_mouse_gene),area2=(length(orth_fly_gene)+length(setdiff(current_fly_gene, gene_list[,2]))),cross.area=length(intersect(current_mouse_gene, orth_fly_gene)),category = c("mouse", "fly"),rotation.degree = 180,filename=NULL,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')
    plot_list <- c(plot_list, list(temp_b))
    common_term_num <- common_term_num + length(intersect(current_mouse_gene, orth_fly_gene))
    data1_term_num <- data1_term_num + length(setdiff(current_mouse_gene, orth_fly_gene))
    data2_term_num <- data2_term_num + length(setdiff(orth_fly_gene, current_mouse_gene))
  }
}
total_num
common_num = common_term_num/total_num
data1_all = data1_term_num/total_num + common_num
data2_all = data2_term_num/total_num + common_num
not_fly_gene = not_fly_gene/total_num
common_num <- round(common_num, digits = 2)
data1_all <- round(data1_all, digits = 2)
data2_all <- round(data2_all, digits = 2)
temp_b <- draw.pairwise.venn(area1=data1_all,area2=(data2_all + not_fly_gene),cross.area=common_num,category = c("mouse", "fly"),rotation.degree = 180,filename=NULL,  scaled=FALSE,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')
grid.draw(temp_b)
percent_num <- common_num / (data1_all + data2_all - common_num + not_fly_gene)
plot_grid(plotlist = plot_list,nrow=12, ncol=2)





#term   tissue average
#fly not ortholgy to mouse
old_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv",header=T,row.names=1)
young_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_youngARSER_output.csv",header=T,row.names=1)
expr_gene <- intersect(rownames(old_data),rownames(young_data))
expr_gene <- bitr(expr_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID
r1 <- read.csv("E:/GO_layers/circadian_data/fly/fly_gene_entrezid.csv")

#write.csv(gene, "E:/GO_layers/photo/result/fly_gene.csv")
#gene_id <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID
#r1 <- gene_id
ck1<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Dm.eg.db",keyType = "ENTREZID",ont="BP",readable=TRUE,universe=expr_gene)
ck2<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Dm.eg.db",keyType = "ENTREZID",ont="MF",readable=TRUE,universe=expr_gene)
ck3<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Dm.eg.db",keyType = "ENTREZID",ont="CC",readable=TRUE,universe=expr_gene)
write.csv(ck1@compareClusterResult, "E:/GO_layers/circadian_data/fly/fly_notorthology_nolevel_bp.csv")
write.csv(ck2@compareClusterResult, "E:/GO_layers/circadian_data/fly/fly_notorthology_nolevel_mf.csv")
write.csv(ck3@compareClusterResult, "E:/GO_layers/circadian_data/fly/fly_notorthology_nolevel_cc.csv")


#mouse
b=list.files("E:/GO_layers/circadian_data/mouse/allgenes12tissues/")
expr_gene <- c()
for(i in 1:length(b))
{
  file_name <- paste("E:/GO_layers/circadian_data/mouse/allgenes12tissues/", b[i],sep="")
  temp_data <- read.csv(file_name, header = TRUE, row.names = 1)
  if(i == 1)
  {
    expr_gene <- temp_data[,1]
  }else
  {
    expr_gene <- intersect(expr_gene, temp_data[,1])
  }
}

expr_gene <- bitr(expr_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
r1<-read.csv("E:/GO_layers/photo/mouse.csv")
mouse_ck1<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",ont="BP",readable=TRUE,universe=expr_gene)
mouse_ck2<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",ont="MF",readable=TRUE,universe=expr_gene)
mouse_ck3<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",ont="CC",readable=TRUE,universe=expr_gene)
write.csv(mouse_ck1@compareClusterResult, "E:/GO_layers/circadian_data/mouse/mouse_notorthology_nolevel_bp.csv")
write.csv(mouse_ck2@compareClusterResult, "E:/GO_layers/circadian_data/mouse/mouse_notorthology_nolevel_mf.csv")
write.csv(mouse_ck3@compareClusterResult, "E:/GO_layers/circadian_data/mouse/mouse_notorthology_nolevel_cc.csv")


#fly mouse not orthology photo not average
length(unique(ck1@compareClusterResult$ID))
length(unique(mouse_ck1@compareClusterResult$ID))
length(intersect(ck1@compareClusterResult$ID, mouse_ck1@compareClusterResult$ID))
temp_b <- draw.pairwise.venn(area1=1425,area2=514,cross.area=200,category = c("mouse", "fly"),filename=NULL,  scaled=FALSE,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')

#fly mouse not orthology photo average tissue
common_term_num <- 0
total_num <- 0
data1_term_num <- 0
data2_term_num <- 0
not_fly_gene <- 0
plot_list = list()
for(i in 1:length(mouse_tissue))
{
  current_mouse_gene <- unique(mouse_ck3@compareClusterResult[mouse_ck3@compareClusterResult$Cluster==mouse_tissue[i],]$ID)
  for(j in 1:length(fly_tissue))
  {
    current_fly_gene <- unique(ck3@compareClusterResult[ck3@compareClusterResult$Cluster==fly_tissue[j],]$ID)
    total_num <- total_num + 1
    orth_fly_gene <- current_fly_gene
    if(length(current_mouse_gene) > length(orth_fly_gene))
    {
      number_flag <- 0
    }else
    {
      number_flag <- 180
    }
    temp_b <- draw.pairwise.venn(area1=length(current_mouse_gene),area2=length(orth_fly_gene),cross.area=length(intersect(current_mouse_gene, orth_fly_gene)),rotation.degree = number_flag ,category = c("mouse", "fly"),filename=NULL,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')
    plot_list <- c(plot_list, list(temp_b))
    common_term_num <- common_term_num + length(intersect(current_mouse_gene, orth_fly_gene))
    data1_term_num <- data1_term_num + length(setdiff(current_mouse_gene, orth_fly_gene))
    data2_term_num <- data2_term_num + length(setdiff(orth_fly_gene, current_mouse_gene))
  }
}
total_num
common_num = common_term_num/total_num
data1_all = data1_term_num/total_num + common_num
data2_all = data2_term_num/total_num + common_num
not_fly_gene = not_fly_gene/total_num
common_num <- round(common_num, digits = 2)
data1_all <- round(data1_all, digits = 2)
data2_all <- round(data2_all, digits = 2)
temp_b <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,category = c("mouse", "fly"),rotation.degree = 180,filename=NULL,  scaled=FALSE,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')
grid.draw(temp_b)
percent_num <- common_num / (data1_all + data2_all - common_num + not_fly_gene)
plot_grid(plotlist = plot_list,nrow=12, ncol=2)
















#fly ortholgy to mouse
b=list.files("E:/GO_layers/circadian_data/mouse/allgenes12tissues/")
expr_gene <- c()
for(i in 1:length(b))
{
  file_name <- paste("E:/GO_layers/circadian_data/mouse/allgenes12tissues/", b[i],sep="")
  temp_data <- read.csv(file_name, header = TRUE, row.names = 1)
  if(i == 1)
  {
    expr_gene <- temp_data[,1]
  }else
  {
    expr_gene <- intersect(expr_gene, temp_data[,1])
  }
}

expr_gene <- bitr(expr_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
fly_gene <- read.csv("E:/GO_layers/circadian_data/fly/fly2mouse_gene.csv")
gene_id <- list()
for(i in 1:ncol(fly_gene))
{
  gene_id[[i]] <- bitr(fly_gene[,i], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
}
temp <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
write.csv(temp, "E:/GO_layers/circadian_data/fly/fly2mouse_gene_geneid.csv")
r1 <- read.csv("E:/GO_layers/circadian_data/fly/fly2mouse_gene_geneid.csv")

#write.csv(gene, "E:/GO_layers/photo/result/fly_gene.csv")
#gene_id <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID
#r1 <- gene_id
ck1<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",keyType = "ENTREZID",ont="BP",readable=TRUE,universe=expr_gene)
ck2<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",keyType = "ENTREZID",ont="MF",readable=TRUE,universe=expr_gene)
ck3<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",keyType = "ENTREZID",ont="CC",readable=TRUE,universe=expr_gene)
write.csv(ck1@compareClusterResult, "E:/GO_layers/circadian_data/fly/fly_orthology_nolevel_bp.csv")
write.csv(ck2@compareClusterResult, "E:/GO_layers/circadian_data/fly/fly_orthology_nolevel_mf.csv")
write.csv(ck3@compareClusterResult, "E:/GO_layers/circadian_data/fly/fly_orthology_nolevel_cc.csv")


fly_ck1 <- read.csv("E:/GO_layers/circadian_data/fly/fly_notorthology_nolevel_bp.csv")
fly_ck2 <- read.csv("E:/GO_layers/circadian_data/fly/fly_notorthology_nolevel_mf.csv")
fly_ck3 <- read.csv("E:/GO_layers/circadian_data/fly/fly_notorthology_nolevel_cc.csv")

length(unique(ck3@compareClusterResult$ID))
length(unique(fly_ck3$ID))
length(intersect(ck3@compareClusterResult$ID, fly_ck3$ID))


#mouse
b=list.files("E:/GO_layers/circadian_data/mouse/allgenes12tissues/")
expr_gene <- c()
for(i in 1:length(b))
{
  file_name <- paste("E:/GO_layers/circadian_data/mouse/allgenes12tissues/", b[i],sep="")
  temp_data <- read.csv(file_name, header = TRUE, row.names = 1)
  if(i == 1)
  {
    expr_gene <- temp_data[,1]
  }else
  {
    expr_gene <- intersect(expr_gene, temp_data[,1])
  }
}

expr_gene <- bitr(expr_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
r1<-read.csv("E:/GO_layers/photo/mouse.csv")
mouse_ck1<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",ont="BP",readable=TRUE,universe=expr_gene)
mouse_ck2<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",ont="MF",readable=TRUE,universe=expr_gene)
mouse_ck3<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",ont="CC",readable=TRUE,universe=expr_gene)
write.csv(mouse_ck1@compareClusterResult, "E:/GO_layers/circadian_data/mouse/mouse_notorthology_nolevel_bp.csv")
write.csv(mouse_ck2@compareClusterResult, "E:/GO_layers/circadian_data/mouse/mouse_notorthology_nolevel_mf.csv")
write.csv(mouse_ck3@compareClusterResult, "E:/GO_layers/circadian_data/mouse/mouse_notorthology_nolevel_cc.csv")
dotplot(mouse_ck1, showCategory=5, includeAll=FALSE)
dotplot(ck1, showCategory=5, includeAll=FALSE)


#fly mouse orthology photo not average
length(unique(ck1@compareClusterResult$ID))
length(unique(mouse_ck1@compareClusterResult$ID))
length(intersect(ck1@compareClusterResult$ID, mouse_ck1@compareClusterResult$ID))
temp_b <- draw.pairwise.venn(area1=1425,area2=514,cross.area=200,category = c("mouse", "fly"),filename=NULL,  scaled=FALSE,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')

#fly mouse not orthology photo average tissue
common_term_num <- 0
total_num <- 0
data1_term_num <- 0
data2_term_num <- 0
not_fly_gene <- 0
plot_list = list()
for(i in 1:length(mouse_tissue))
{
  current_mouse_gene <- unique(mouse_ck3@compareClusterResult[mouse_ck3@compareClusterResult$Cluster==mouse_tissue[i],]$ID)
  for(j in 1:length(fly_tissue))
  {
    current_fly_gene <- unique(ck3@compareClusterResult[ck3@compareClusterResult$Cluster==fly_tissue[j],]$ID)
    total_num <- total_num + 1
    orth_fly_gene <- current_fly_gene
    if(length(current_mouse_gene) > length(orth_fly_gene))
    {
      number_flag <- 0
    }else
    {
      number_flag <- 180
    }
    temp_b <- draw.pairwise.venn(area1=length(current_mouse_gene),area2=length(orth_fly_gene),cross.area=length(intersect(current_mouse_gene, orth_fly_gene)),rotation.degree = number_flag ,category = c("mouse", "fly"),filename=NULL,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')
    plot_list <- c(plot_list, list(temp_b))
    common_term_num <- common_term_num + length(intersect(current_mouse_gene, orth_fly_gene))
    data1_term_num <- data1_term_num + length(setdiff(current_mouse_gene, orth_fly_gene))
    data2_term_num <- data2_term_num + length(setdiff(orth_fly_gene, current_mouse_gene))
  }
}
total_num
common_num = common_term_num/total_num
data1_all = data1_term_num/total_num + common_num
data2_all = data2_term_num/total_num + common_num
not_fly_gene = not_fly_gene/total_num
common_num <- round(common_num, digits = 2)
data1_all <- round(data1_all, digits = 2)
data2_all <- round(data2_all, digits = 2)
temp_b <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,category = c("mouse", "fly"),rotation.degree = 180,filename=NULL,  scaled=FALSE,ext.text=FALSE,offset=0,resolution =300, imagetype="png", col=c('white','white'), cat.pos=c(0,0),fill=c(colors()[38], colors()[616]),  cex='white',cat.cex='white')
grid.draw(temp_b)
percent_num <- common_num / (data1_all + data2_all - common_num + not_fly_gene)
plot_grid(plotlist = plot_list,nrow=12, ncol=2)





#orthology phase

gene_list <- read.csv("E:/GO_layers/blast/restart/ensembl/mouse_fly.csv",header=TRUE)
b=list.files("E:/GO_layers/circadian_data/mouse/allgenes12tissues/")
expr_gene <- c()
for(i in 1:length(b))
{
  file_name <- paste("E:/GO_layers/circadian_data/mouse/allgenes12tissues/", b[i],sep="")
  temp_data <- read.csv(file_name, header = TRUE, row.names = 1)
  if(i == 1)
  {
    expr_gene <- temp_data[,1]
  }else
  {
    expr_gene <- intersect(expr_gene, temp_data[,1])
  }
}

mouse_expr_gene <- expr_gene



old_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv",header=TRUE,row.names=1)
young_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_youngARSER_output.csv",header=TRUE,row.names=1)
expr_gene <- intersect(rownames(old_data),rownames(young_data))
fly_expr_gene <- expr_gene
orth_mouse_gene <- unique(gene_list[,1])
orth_fly_gene <- unique(gene_list[,2])
fly_overlap <- intersect(fly_expr_gene, orth_fly_gene)
mouse_overlap <- intersect(mouse_expr_gene,gene_list[,1])
new_gene_list <- gene_list[gene_list[,1] %in% mouse_overlap,]
expr_gene_list <- new_gene_list[new_gene_list[,2] %in% fly_overlap,]




phase_plot_data <- data.frame(0, 0)
mouse_gene <- read.csv("E:/GO_layers/photo/result/mouse/mouse_genesymbol.csv")

mouse_gene <- read.csv("E:/GO_layers/blast/circidian_gene/mouse_gene.csv")
mouse_gene <- mouse_gene[,1]
fly_gene <- read.csv("E:/GO_layers/circadian_data/fly/fly_gene.csv")
fly_gene <- fly_gene[,1]





mouse_jtk_data <- read.csv("E:/GO_layers/circadian_data/mouse/allgenes12tissues/Adr_Jtk.csv")

mous_symbol <- mouse_anno[mouse_anno$probeset_id %in% probe_id,]$gene_symbol_mouse
{
  current_fly_gene <- fly_gene[fly_gene %in% expr_gene_list[,2]]
  fly_jtk_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv")
  for(i in 1:length(current_fly_gene))
  {
    loop_fly_phase <- fly_jtk_data[fly_jtk_data$Gene.Name == as.character(current_fly_gene[i]),"phase"]
    fly2mouse_gene <- expr_gene_list[expr_gene_list[,2] %in% current_fly_gene[i], 1]
    for(j in 1:length(fly2mouse_gene))
    {
      if(fly2mouse_gene[j] %in% mouse_jtk_data$X.1 == TRUE)
      {
        #print("true")
        loop_mouse_phase <- mouse_jtk_data[which(mouse_jtk_data$X.1 %in% fly2mouse_gene[j]), "LAG"]
        phase_plot_data <- rbind(phase_plot_data, c(loop_fly_phase, loop_mouse_phase))
      }
    }
  }
}



plot(x=phase_plot_data[,1], y=phase_plot_data[,2],xlab="fly phase", ylab = "mouse phase")
for(i in 1:length(phase_plot_data[,1]))
{
  if (is.na(phase_plot_data[i,1]))
  {
    phase_plot_data[i,1] <- 0
  }
}
for(i in 1:length(phase_plot_data[,2]))
{
  if (is.na(phase_plot_data[i,2]))
  {
    phase_plot_data[i,2] <- 0
  }
}
cor.test(as.numeric(phase_plot_data[,1]), as.numeric(phase_plot_data[,2]))







gene_list <- read.csv("E:/GO_layers/blast/restart/ensembl/mouse_fly.csv",header=TRUE)
mouse_gene <- read.csv("E:/GO_layers/blast/circidian_gene/mouse_gene.csv")
mouse_gene <- mouse_gene[,1]
fly_gene <- read.csv("E:/GO_layers/blast/circidian_gene/fly_gene.csv")
fly_gene <- fly_gene[,1]
fly2mouse <- unique(gene_list[gene_list[,2] %in% fly_gene, 1])
flynot2mouse <- unique(gene_list[gene_list[,2] %in% fly_gene, 2])
#mouse to fly
mouse2fly <- gene_list[gene_list[,1] %in% mouse_gene, 2]

#mouse to fly and fly to mouse overlap
mouse2fly2mouse <- gene_list[gene_list[,2] %in% mouse2fly, 1]

same_gene <- intersect(mouse2fly2mouse, fly2mouse)
fly_orthology_gene <- unique(gene_list[gene_list[,1] %in% same_gene, 2])
mouse_orthology_gene <- same_gene


phase_plot_data <- data.frame(0, 0)
mouse_jtk_data <- read.csv("E:/GO_layers/circadian_data/mouse/allgenes12tissues/Kidney_Jtk.csv")
fly_jtk_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv")
for(i in 1:length(same_gene))
{
  loop_fly_phase <- fly_jtk_data[fly_jtk_data$Gene.Name %in% fly_orthology_gene[i],"phase"]
  loop_mouse_phase <- mouse_jtk_data[which(mouse_jtk_data$X.1 %in% mouse_orthology_gene[i]), "LAG"]
  phase_plot_data <- rbind(phase_plot_data, c(loop_fly_phase, loop_mouse_phase))
}


plot(x=phase_plot_data[,1], y=phase_plot_data[,2],xlab="fly phase", ylab = "mouse phase")
for(i in 1:length(phase_plot_data[,1]))
{
  if (is.na(phase_plot_data[i,1]))
  {
    phase_plot_data[i,1] <- 0
  }
}
for(i in 1:length(phase_plot_data[,2]))
{
  if (is.na(phase_plot_data[i,2]))
  {
    phase_plot_data[i,2] <- 0
  }
}
cor.test(as.numeric(phase_plot_data[,1]), as.numeric(phase_plot_data[,2]))










number <- c(313,379,55,47,475,43,9,417,20,284, 500,592,72,23,7, 15, 39,9,6,31,2,42,62,58,13,4,5,10,55,2,2,65,1,37,118,67,6,12)
names(number) <- c("OLD","YOUNG","ADR","AOR","BFAT","BST","CERE","HEA","HYPO","KID","LIV","LUN","MUS","WAT","ADR&OLD","AOR&OLD","BFAT&OLD","BST&OLD","CERE&OLD","HEA&OLD","HYPO&OLD","KID&OLD","LIV&OLD","LUN&OLD","MUS&OLD","WAT&OLD","ADR&OLD","AOR&YOUNG","BFAT&YOUNG","BST&YOUNG","CERE&YOUNG","HEA&YOUNG","HYPO&YOUNG","KID&YOUNG","LIV&YOUNG","LUN&YOUNG","MUS&YOUNG","WAT&YOUNG")
data <- fromExpression(number)
upset(data,nsets = 14)



number <- c(4607,4510,445,316,1055,156,191,787,149,2059, 2631,1949,275,286,93, 62, 205,37,49,152,34,426,514,375,60,61,105,66,202,36,42,171,37,415,530,352,55,66)
names(number) <- c("OLD","YOUNG","ADR","AOR","BFAT","BST","CERE","HEA","HYPO","KID","LIV","LUN","MUS","WAT","ADR&OLD","AOR&OLD","BFAT&OLD","BST&OLD","CERE&OLD","HEA&OLD","HYPO&OLD","KID&OLD","LIV&OLD","LUN&OLD","MUS&OLD","WAT&OLD","ADR&OLD","AOR&YOUNG","BFAT&YOUNG","BST&YOUNG","CERE&YOUNG","HEA&YOUNG","HYPO&YOUNG","KID&YOUNG","LIV&YOUNG","LUN&YOUNG","MUS&YOUNG","WAT&YOUNG")
data <- fromExpression(number)
upset(data,nsets = 14)






arabid_gene <- read.csv("E:/GO_layers/circadian_data/拟南芥/circadian_gene_arabidopsis.csv")
arabid_anno <- read.csv("E:/GO_layers/circadian_data/拟南芥/gse5612_annotation.csv")
arabid_symbol <- arabid_anno[arabid_anno$ID %in% arabid_gene[,1], ]$Gene.symbol
write.csv(arabid_symbol,"E:/GO_layers/circadian_data/拟南芥/symbol.csv")



human_cell_gene <- read.csv("E:/GO_layers/circadian_data/human/4tissue/cell/circadian_gene.csv")
human_anno <- read.csv("E:/GO_layers/circadian_data/human/4tissue/cell/2.csv")
gene_name <- human_cell_gene$cycID













#################################################amp and exp new
library(plotrix)
# Create the input vectors.
colors = c("green","orange","brown")
months <- c("Mar","Apr","May","Jun","Jul")
regions <- c("East","West","North")

# Create the matrix of the values.
Values <- matrix(c(2,9,3,11,9,4,8,7,3,12,5,2,8,10,11), nrow = 3, ncol = 5, byrow = TRUE)

# Give the chart file a name
png(file = "barchart_stacked.png")

# Create the bar chart
barplot(Values, main = "total revenue", names.arg = months, xlab = "month", ylab = "revenue", col = colors)

# Add the legend to the chart
legend("topleft", regions, cex = 1.3, fill = colors)





colors <- c(colors()[38], colors()[616])



#mosue
tissuename <- c("Adr","Aorta","BFAT","BS","Cere","Heart","Hypo","Kidney","Liver","Lung","Mus","WFAT")
aa= matrix(0,nrow = 12,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 12,ncol = 2)
for(i in 1:12){
  expr = read.csv(paste0("E:/GO_layers/circadian_data/mouse/allgenes12tissues/",tissuename[i],'_Jtk.csv'))
  cir = expr[,'ADJ.P']<0.05
  expr_cir = expr[cir,]
  expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
  expr_cirmean= rowMeans(expr_cir[,8:31])
  aa[i,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
  aa[i,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
  aa[i,3]= aa[i,1]^2
  aaa[,1]= aa[,'R2']
  aaa[,2]= 1- aa[,'R2']
  pie3D(aaa[i,],main = tissuename[i],pty="s",col=colors)

}
pie3D(aaa[i,],main = tissuename[i],pty="s",col=colors)
colors <- c(colors()[38], colors()[616])
barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)


#human
tissuename<- c("blood","cell","IPS","gse63546","gse76368","gse2703")
aa= matrix(0,nrow = 6,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 6,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/human/blood.csv")
cir = expr[,'meta3d_Pvalue']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[-which(expr_cir$X==0),]
expr_cir <- expr_cir[expr_cir[,'meta3d_AMP']!=0,]
aa[1,1]=cor.test(log(as.numeric(expr_cir[,'meta3d_Base'])),log(as.numeric(expr_cir[,'meta3d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cir[,'meta3d_Base'])),log(as.numeric(expr_cir[,'meta3d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2


aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)

expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/human/cell.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[-which(expr_cir$Gene.symbol==""),]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:48])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[2,1]= aa[1,'R2']
aaa[2,2]= 1- aa[1,'R2']
pie3D(aaa[2,],main = tissuename[2],pty="s",col=colors)



expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/human/ips.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[-which(expr_cir$X ==0),]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= expr_cir[,"X.1"]
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2
aaa[3,1]= aa[1,'R2']
aaa[3,2]= 1- aa[1,'R2']
pie3D(aaa[3,],main = tissuename[2],pty="s",col=colors)



expr = read.csv("E:/GO_layers/circadian_data/human/new_data/JTKresult_data_gse63546 (1).csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,7:20])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2
aaa[4,1]= aa[1,'R2']
aaa[4,2]= 1- aa[1,'R2']
pie3D(aaa[3,],main = tissuename[2],pty="s",col=colors)


expr = read.csv("E:/GO_layers/circadian_data/human/new_data/JTKGSE76368_data1 (1).csv")
cir = expr[,'JTK_pvalue']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'JTK_amplitude']!=0,]
expr_cirmean= rowMeans(expr_cir[,10:17])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'JTK_amplitude'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'JTK_amplitude'])))$p.value
aa[1,3]= aa[1,1]^2
aaa[5,1]= aa[1,'R2']
aaa[5,2]= 1- aa[1,'R2']
pie3D(aaa[3,],main = tissuename[2],pty="s",col=colors)


expr = read.csv("E:/GO_layers/circadian_data/human/new_data/JTKresult_data1_gse2703 (1).csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:13])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2
aaa[6,1]= aa[1,'R2']
aaa[6,2]= 1- aa[1,'R2']
pie3D(aaa[3,],main = tissuename[2],pty="s",col=colors)














#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#fly
tissuename<- c("young","old")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyyoung.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:19])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)

expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyold.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:19])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[2,1]= aa[1,'R2']
aaa[2,2]= 1- aa[1,'R2']
pie3D(aaa[2,],main = tissuename[2],pty="s",col=colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#yeast
tissuename<- c("Sample2","Sample6")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample2.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,7:26])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)

expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample6.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,7:30])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[2,1]= aa[1,'R2']
aaa[2,2]= 1- aa[1,'R2']
pie3D(aaa[2,],main = tissuename[2],pty="s",col=colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#arabidopsis
tissuename<- c("flower","seed")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/flower.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,7:24])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)

expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/JTK.Agse5612.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:20])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[2,1]= aa[1,'R2']
aaa[2,2]= 1- aa[1,'R2']
pie3D(aaa[2,],main = tissuename[2],pty="s",col=colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#neurospora
tissuename<- c("neurospora")
aa= matrix(0,nrow = 1,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 1,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/neurospora/JTKresult_RNA.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,7:18])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#cyanobacteria
tissuename<- c("cyanobacteria")
aa= matrix(0,nrow = 1,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 1,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/cyanobacteria/JTK.LLREPLICATE2.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,9:28])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)


#archaea
tissuename<- c("archaea")
aa= matrix(0,nrow = 1,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 1,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/archaea/meta2d_AR2_achaea.csv")
cir = expr[,'JTK_pvalue']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'JTK_amplitude']!=0,]
expr_cirmean= rowMeans(expr_cir[,20:55])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'JTK_amplitude'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'JTK_amplitude'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
pie3D(aaa[1,],main = tissuename[1],pty="s",col=colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)



################################################percent plot
#mouse
tissuename= c("Adr","Aorta","BFAT","BS","Cere","Heart","Hypo","Kidney","Liver","Lung","Mus","WFAT")
aa= matrix(-1,nrow = 12,ncol = 5)
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 12,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)
for(i in 1:12){
  expr = read.csv(paste0("E:/GO_layers/circadian_data/mouse/allgenes12tissues/",tissuename[i],'_Jtk.csv'))
  bhp= expr[,'BH.Q']
  bhp= bhp[order(rowMeans(expr[,7:30]))]
  smallp= bhp<0.05
  len = length(smallp)
  n=5
  quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
  quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
  prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
  aa[i,]=prop
  numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
  bb[i,]=numb
  cat(tissuename[i],'   ',bb[i,],'\n')
  bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])
  
  diffextrnum=c(0,0,0,0)
  for(j in 1:10000){
    aab=rep(0,length(bhp))
    aab[sample(1:length(bhp),sum(smallp))]=1
    aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
    aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
    diffextrnum= diffextrnum+(aad>=bbdiff)
  }
  diffpvalue[i,]=diffextrnum/10000
}

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/rising_trend_percentage.csv')
tissuenamefull= c("Adrenal Gland","Aorta","Brown Fat","Brainstem","Cerebellum","Heart","Hypothalamus",
                  "Kidney","Liver","Lung","Skeletal Muscle","White Fat")
{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/rising_trend_percentage.csv',row.names=1)
  aa= aa/rowSums(aa)
  
  filcolor <-  c(colors[589],colors()[590],colors()[591],colors()[592],colors()[593])
  # png("fig1.png",width=10,height=10,units="in",res=350)
  for(i in 1:12){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # #plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  parmai=par()$mai
  par(mai=rep(0,4))
  plot.new()
  plot.new()
  lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  legend("topleft",lbls,
         fill= filcolor,
         cex=1.5,bty = 'n')
  par(mai=parmai)
  # dev.off()
}





#human
aa= matrix(-1,nrow = 6,ncol = 5)
rownames(aa)=c("blood","cell","IPS","gse63546","gse76368","gse2703")
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 6,ncol = 4)
rownames(diffpvalue)=c("blood","cell","IPS","gse63546","gse76368","gse2703")
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)


#blood
blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/human/blood.csv")
mean_exp <- order(as.numeric(blood_data[,'meta3d_Base']))
bhp = blood_data[,'meta3d_Pvalue']
bhp=bhp[order(mean_exp)]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000




#bcell
bcell_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/human/cell.csv")
expr_cirmean <- c()
for(i in 1:nrow(bcell_data))
{
  expr_cirmean[length(expr_cirmean) + 1] <- mean(as.numeric(bcell_data[i,8:48]))
}
i=2
bhp = bcell_data[,'ADJ.P']
bhp=bhp[order(rowMeans(bcell_data[,8:48]))]
smallp= bhp<0.05
len = length(smallp)
n=5
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[2,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[2,]=numb
cat(tissuename[2],'   ',bb[2,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000






#ips
bcell_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/human/ips.csv")
#expr_cir <- bcell_data
expr_cir <- bcell_data[-which(bcell_data$X ==0),]
expr_cirmean= expr_cir[,"X.1"]

bcell_data <- expr_cir

i=3
bhp = bcell_data[,'ADJ.P']
bhp=bhp[order(expr_cir[,"X.1"])]
smallp= bhp<0.05
len = length(smallp)
n=5
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[3,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[3,]=numb
cat(tissuename[3],'   ',bb[3,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000



#gse63546
bcell_data <- read.csv("E:/GO_layers/circadian_data/human/new_data/JTKresult_data_gse63546 (1).csv")
#expr_cir <- bcell_data
expr_cir <- bcell_data
expr_cirmean= rowMeans(expr_cir[,7:20])

bcell_data <- expr_cir

i=4
bhp = bcell_data[,'ADJ.P']
bhp=bhp[order(expr_cirmean)]
smallp= bhp<0.05
len = length(smallp)
n=5
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[i,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[i,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000



#gse76368
bcell_data <- read.csv("E:/GO_layers/circadian_data/human/new_data/JTKGSE76368_data1 (1).csv")
expr_cir <- bcell_data
expr_cirmean= rowMeans(expr_cir[,10:17])

bcell_data <- expr_cir

i=5
bhp = bcell_data[,'JTK_pvalue']
bhp=bhp[order(expr_cirmean)]
smallp= bhp<0.05
len = length(smallp)
n=5
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[i,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[i,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000




#gse2703
bcell_data <- read.csv("E:/GO_layers/circadian_data/human/new_data/JTKresult_data1_gse2703 (1).csv")
expr_cir <- bcell_data
expr_cirmean= rowMeans(expr_cir[,8:13])

bcell_data <- expr_cir

i=6
bhp = bcell_data[,'ADJ.P']
bhp=bhp[order(expr_cirmean)]
smallp= bhp<0.05
len = length(smallp)
n=5
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[i,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[i,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000












tissuenamefull <- c("blood","cell","IPS","gse63546","gse76368","gse2703")

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/humanrising_trend_percentage_2.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/humanrising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  for(i in 1:6){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  parmai=par()$mai
  par(mai=rep(0,4))
  plot.new()
  plot.new()
  lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  legend("topleft",lbls,
         fill= filcolor,
         cex=1.5,bty = 'n')
  par(mai=parmai)
  # dev.off()
}








#fly
aa= matrix(-1,nrow = 2,ncol = 5)
rownames(aa)=c('Old','Young')
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 2,ncol = 4)
rownames(diffpvalue)=c('Old','Young')
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)





tissuename<- c("young","old")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyyoung.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:19])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']

expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyold.csv")
cir = expr[,'ADJ.P']<0.05
expr_cir = expr[cir,]
expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,8:19])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
aa[1,3]= aa[1,1]^2


















#old
blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyold.csv")
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,8:19]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000


#young
blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyyoung.csv")
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,8:19]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=2
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[2,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[2,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000



tissuenamefull <- c('Old','Young')

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/fly_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/fly_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  #par(mfrow=c(4,4))
  for(i in 1:2){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  # parmai=par()$mai
  # par(mai=rep(0,4))
  # plot.new()
  # plot.new()
  # lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  # legend("topleft",lbls,
  #        fill= filcolor,
  #        cex=1.5,bty = 'n')
  # par(mai=parmai)
  # dev.off()
}








#yeast
aa= matrix(-1,nrow = 2,ncol = 5)
tissuename <- c("sample2","sample6")
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 2,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)




#sample2
blood_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/JTK.Sample2.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,6:25]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000


#sample6
blood_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/JTK.Sample6.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,6:25]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=2
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[2,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[2,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000



tissuenamefull <- tissuename

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/yeast_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/yeast_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  #par(mfrow=c(4,4))
  for(i in 1:2){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  # parmai=par()$mai
  # par(mai=rep(0,4))
  # plot.new()
  # plot.new()
  # lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  # legend("topleft",lbls,
  #        fill= filcolor,
  #        cex=1.5,bty = 'n')
  # par(mai=parmai)
  # dev.off()
}




#arabidopsis
aa= matrix(-1,nrow = 2,ncol = 5)
tissuename <- c("flower","seed")
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 2,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)




blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/flower.csv")
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,7:24]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000

blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/JTK.Agse5612.csv")
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,8:20]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=2
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[2,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[2,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000



tissuenamefull <- tissuename

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/arabidopsis_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/arabidopsis_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  #par(mfrow=c(4,4))
  for(i in 1:2){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  # parmai=par()$mai
  # par(mai=rep(0,4))
  # plot.new()
  # plot.new()
  # lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  # legend("topleft",lbls,
  #        fill= filcolor,
  #        cex=1.5,bty = 'n')
  # par(mai=parmai)
  # dev.off()
}




#neurospora
aa= matrix(-1,nrow = 1,ncol = 5)
tissuename <- c("Neurospora")
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 1,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)



blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/neurospora/JTKresult_RNA.csv")
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,7:18]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000


tissuenamefull <- tissuename

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/neurospora_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/neurospora_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  #par(mfrow=c(4,4))
  for(i in 1:1){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  # parmai=par()$mai
  # par(mai=rep(0,4))
  # plot.new()
  # plot.new()
  # lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  # legend("topleft",lbls,
  #        fill= filcolor,
  #        cex=1.5,bty = 'n')
  # par(mai=parmai)
  # dev.off()
}




#cyanobacteria
aa= matrix(-1,nrow = 1,ncol = 5)
tissuename <- c("Cyanobacteria")
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 1,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)




blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/cyanobacteria/JTK.LLREPLICATE2.csv")
bhp = blood_data[,'ADJ.P']
bhp=bhp[order(rowMeans(blood_data[,9:28]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000


tissuenamefull <- tissuename

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/cyanobacteria_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/cyanobacteria_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  #par(mfrow=c(4,4))
  for(i in 1:1){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  parmai=par()$mai
  par(mai=rep(0,4))
  plot.new()
  plot.new()
  lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  legend("topleft",lbls,
         fill= filcolor,
         cex=1.5,bty = 'n')
  par(mai=parmai)
  # dev.off()
}


#archaea
aa= matrix(-1,nrow = 1,ncol = 5)
tissuename <- c("archaea")
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 1,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)




blood_data <- read.csv("E:/GO_layers/final_photo/exp_amp_data/archaea/meta2d_AR2_achaea.csv")
bhp = blood_data[,'JTK_pvalue']
bhp=bhp[order(rowMeans(blood_data[,20:55]))]
smallp= bhp<0.05
len = length(smallp)
n=5
i=1
quantilepoint= apply(as.matrix(1:n),1,function(x)floor(len/n*x))
quantilepoint= cbind(c(0,quantilepoint[1:(n-1)])+1,quantilepoint)
prop= apply(quantilepoint,1,function(x)mean(smallp[x[1]:x[2]]))
aa[1,]=prop
numb= apply(quantilepoint,1,function(x)sum(smallp[x[1]:x[2]]))
bb[1,]=numb
cat(tissuename[i],'   ',bb[i,],'\n')
bbdiff=c(bb[i,2]-bb[i,1],bb[i,3]-bb[i,2],bb[i,4]-bb[i,3],bb[i,5]-bb[i,4])

diffextrnum=c(0,0,0,0)
for(j in 1:10000){
  aab=rep(0,length(bhp))
  aab[sample(1:length(bhp),sum(smallp))]=1
  aac= apply(quantilepoint,1,function(x)sum(aab[x[1]:x[2]]))
  aad= c(aac[2]-aac[1],aac[3]-aac[2],aac[4]-aac[3],aac[5]-aac[4])
  diffextrnum= diffextrnum+(aad>=bbdiff)
}
diffpvalue[1,]=diffextrnum/10000


tissuenamefull <- tissuename

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/archaea_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/archaea_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  #filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  #par(mfrow=c(4,4))
  for(i in 1:1){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    # parmai=par()$mai
    # par(new= T)
    # par(mai=rep(0,4))
    # plot.new()
    # text(0,1,letters[i],cex=2)
    # par(mai=parmai)
  }
  # legend
  parmai=par()$mai
  par(mai=rep(0,4))
  plot.new()
  plot.new()
  lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  legend("topleft",lbls,
         fill= filcolor,
         cex=1.5,bty = 'n')
  par(mai=parmai)
  # dev.off()
}







######################################################orthoFinder othology gene photo
species_sort <- c('human','mouse','fly','arabidopsis','yeast','neurospora','cyanobacteria','archaea')
species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','cyanobacteria','archaea')
gene_num_list <- data.frame(0,0)
colnames(gene_num_list) <- c("species","gene_num")

venn_data <- data.frame(0,0,0,0,0)
colnames(venn_data) <- c("species1","species2","num1","num2","overlap")
for(i in 1:(length(species) - 1))
{
  for(j in (i+1):length(species))
  {
    current_file_name <- paste(paste(paste(paste("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/",species[i],sep=""),'__v__',sep=""), species[j],sep=""),".txt",sep="")
    if (file.exists(current_file_name))
    {
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[i],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[j],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      
      
      
      gene1_num <- 0
      gene2_num <- 0
      overlap <- 0
      num <- 0
      gene1_total_num <- c()
      gene2_total_num <- 0
      for(k in 1:ncol(gene1))
      {
        for(l in 1:ncol(gene2))
        {
          num <- num + 1
          if(which(species_sort == species[i]) < which(species_sort == species[j]))
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V1
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
          }else
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V2
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
          }
          overlap <- overlap + length(intersect(current_gene1, current_gene2))
          
        }
      }
      gene1_num <- gene1_num / num
      gene2_num <- gene2_num / num
      overlap <- overlap / num
      
      # gene1 <- unique(gene1[,1])
      # gene2 <- unique(gene2[,1])
      # gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene1))))
      # gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene2))))
      # orth_gene1 <- gene_list[gene_list[,1]%in%gene1 ,2]
      # orth_gene2 <- gene_list[gene_list[,2]%in%gene2 ,2]
      # orth_gene1 <- na.omit(orth_gene1)
      # orth_gene2 <- na.omit(orth_gene2)
      # temp_row <- c(species[i], species[j], length(unique(orth_gene1)), length(unique(orth_gene2)), length(unique(intersect(orth_gene1, orth_gene2))))
      temp_row <- c(species[i], species[j], gene1_num, gene2_num, overlap)
      venn_data <- rbind(venn_data, temp_row)
      gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene1_total_num))))
      gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene2_total_num))))
    }else
    {
      current_file_name <- paste(paste(paste(paste("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/",species[j],sep=""),'__v__',sep=""), species[i],sep=""),".txt",sep="")
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[j],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[i],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      gene1_num <- 0
      gene2_num <- 0
      overlap <- 0
      num <- 0
      for(k in 1:ncol(gene1))
      {
        for(l in 1:ncol(gene2))
        {
          num <- num + 1
          if(which(species_sort == species[i]) > which(species_sort == species[j]))
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V1
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
          }else
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V2
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
          }
          overlap <- overlap + length(intersect(current_gene1, current_gene2))
          
        }
      }
      gene1_num <- gene1_num / num
      gene2_num <- gene2_num / num
      overlap <- overlap / num
      temp_row <- c(species[j], species[i], gene1_num, gene2_num, overlap)
      venn_data <- rbind(venn_data, temp_row)
      
      
      
      
      
      
      
      
      # gene1 <- gene1[,1]
      # gene2 <- gene2[,1]
      # gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene1))))
      # gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene2))))
      # orth_gene1 <- gene_list[ gene_list[,1]%in%gene1 ,2]
      # orth_gene2 <- gene_list[ gene_list[,2]%in%gene2 ,2]
      # orth_gene1 <- na.omit(orth_gene1)
      # orth_gene2 <- na.omit(orth_gene2)
      # temp_row <- c(species[i], species[j], length(unique(orth_gene2)), length(unique(orth_gene1)), length(unique(intersect(orth_gene1, orth_gene2))))
      # venn_data <- rbind(venn_data, temp_row)
    }
  }
}





gene_num_list <- data.frame(0,0)
colnames(gene_num_list) <- c("species","gene_num")
for(i in 1:length(species))
{
  current_file_gene_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[i],sep=""),"_gene.csv",sep="")
  current_data_gene <- read.csv(current_file_gene_name)
  for(j in 1:ncol(current_data_gene))
  {
    if(j == 1)
    {
      temp_gene <- as.character(current_data_gene[,j])
    }else
    {
      temp_gene <- c(temp_gene, as.character(current_data_gene[,j]))
    }
  }
  gene_num_list <- rbind(gene_num_list, c(species[i], as.character(length(unique(temp_gene)))))
}




gene_num_list <- gene_num_list[-1,]
gene_num_list$gene_num <- as.numeric(gene_num_list$gene_num)
venn_data_new <- venn_data[2:nrow(venn_data),]
plot_list = list()
for(i in 1:nrow(venn_data_new))
{
  data1_all = as.numeric(venn_data_new[i,3])
  data2_all = as.numeric(venn_data_new[i,4])
  common_num = as.numeric(venn_data_new[i,5])
  if(data1_all > data2_all)
  {
    number_flag = 0
  }else
  {
    number_flag = 180
  }
  temp_plot <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, ext.text=FALSE, height = 450, width = 450, rotation.degree = number_flag, resolution =300, imagetype="png", col="black", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
  plot_list <- c(plot_list, list(temp_plot))
}


overlap_data <- data.frame(nrow=21, ncol=2)
for(i in 1:nrow(venn_data_new))
{
  if(i == 1)
  { 
    name = paste(venn_data_new[i,1], venn_data_new[i,2],sep="_")
    overlap_data[i,] <- c(name,as.numeric(venn_data_new$overlap[i]))
  }else
  {
    name = paste(venn_data_new[i,1], venn_data_new[i,2],sep="_")
    overlap_data[i,] <- c(name,as.numeric(venn_data_new$overlap[i]))
  }
}

overlap_data$ncol <- as.numeric(overlap_data$ncol)
#plot_data1 <- list()
#plot_data1[c(2,3,4,5,6,7,10,11,12,13,14,18,19,20,21,26,27,28,34,35,42)] = plot_list 
#plot_grid(plotlist = plot_data1)
plot_data1 <- list()
temp_matrix <- matrix(1:64,byrow=TRUE,nrow=8)
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = plot_list 
grid.newpage()
species1 <- c('archaea','cyanobacteria','arabidopsis','neurospora','yeast','fly','mouse','human')
p.hist.len <- ggplot(gene_num_list, mapping=aes(x=species,y=gene_num,fill=species,group=factor(1))) + geom_bar(stat="identity")
#p.hist.wid <- ggplot(iris) + geom_histogram(aes(x=Sepal.Width)) + coord_flip()
p.hist.wid <- ggplot(overlap_data, mapping=aes(x=nrow,y=ncol,fill=overlap_data[,1], group=factor(1))) + geom_bar(stat="identity") + coord_flip() + theme(legend.position = "none")
p.scatter<-ggplot(data=point_data, aes(x=point_data[,1], y=point_data[,2]))+geom_point(color='1',size=0) + scale_x_discrete(species) + theme_bw()  + theme(panel.border = element_blank() ,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

vp.len <- viewport(x=0, y=0.66, width=0.66, height=0.34, just=c("left", "bottom"))
vp.wid <- viewport(x=0.66, y=0, width=0.34, height=0.66, just=c("left", "bottom"))
vp.scatter <- viewport(x=0, y=0, width=0.66, height=0.66, just=c("left", "bottom"))
vp.scatter1 <- viewport(x=0.11, y=0.09, width=0.50, height=0.50, just=c("left", "bottom"))

print(p.hist.len, vp=vp.len)
print(p.hist.wid, vp=vp.wid)
species_test <- c('arabidopsis','archaea','cyanobacteria','fly','mouse','neurospora','yeast')
species_test1 <- c('yeast','neurospora','mouse','fly','cyanobacteria','archaea','arabidopsis')
print(as.ggplot(expression(plot(1:8,1:8,type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:8,labels = species) + axis(side = 4,at = 1:8,labels = species1))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=8, ncol=8, axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)






######################################################orthoFinder othology term photo
library(org.Mm.eg.db)#Mouse
library(org.Hs.eg.db)#Human   babbon also
library(org.Dm.eg.db)#Fly
library(org.At.tair.db)#Arabidopsis
library(org.Sc.sgd.db)#yeast
library(AnnotationDbi)
library(GO.db)
species_sort <- c('human','mouse','fly','arabidopsis','yeast','neurospora','cyanobacteria','archaea')
species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','cyanobacteria','archaea')
mouse_expr_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/mouse_gene.csv")
mouse_expr_gene <- mouse_expr_gene$x
human_expr_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/human_gene.csv")
human_expr_gene <- human_expr_gene$x
fly_expr_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/fly_gene.csv")
fly_expr_gene <- fly_expr_gene$x
yeast_expr_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/yeast_gene.csv")
yeast_expr_gene <- yeast_expr_gene$x
arabid_expr_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/arabidopsis_gene.csv")
arabid_expr_gene <- arabid_expr_gene$x
expr_list <- list(mouse_expr_gene,human_expr_gene,fly_expr_gene,yeast_expr_gene,arabid_expr_gene)
names(expr_list) <- c("mouse","human","fly","yeast","arabidopsis")
go_list <- list(org.Mm.eg.db, org.Hs.eg.db, org.Dm.eg.db, org.Sc.sgd.db, org.At.tair.db)
names(go_list) <- c("mouse","human","fly","yeast","arabidopsis")


for(i in 1:(length(species) - 1))
{
  for(j in (i+1): length(species))
  {
    current_file_name <- paste(paste(paste(paste("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/",species[i],sep=""),'__v__',sep=""), species[j],sep=""),".txt",sep="")
    if (file.exists(current_file_name))
    {
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[i],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[j],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      
      
      
      species1_new_gene_list <- list()
      species2_new_gene_list <- list()
      for(k in 1:ncol(gene1))
      {
        if(which(species_sort == species[i]) < which(species_sort == species[j]))
        {
          current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
          species1_new_gene_list[[k]] <- unique(as.character(current_gene1))
          
        }else
        {
          current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
          species1_new_gene_list[[k]] <- unique(as.character(current_gene1))
        }
      }
      for(l in 1:ncol(gene2))
      {
        if(which(species_sort == species[i]) < which(species_sort == species[j]))
        {
          current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
          species2_new_gene_list[[l]] <- unique(as.character(current_gene2))
          
        }else
        {
          current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
          species2_new_gene_list[[l]] <- unique(as.character(current_gene2))
        }
      }
      temp1 <- do.call(cbind, lapply(lapply(species1_new_gene_list, unlist), `length<-`, max(lengths(species1_new_gene_list))))
      temp2 <- do.call(cbind, lapply(lapply(species2_new_gene_list, unlist), `length<-`, max(lengths(species2_new_gene_list))))
      file1_name <- paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],'_vs_',species[j],'_',species[i],".csv",sep="")
      file2_name <- paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],'_vs_',species[j],'_',species[j],".csv",sep="")
      colnames(temp1) <- colnames(gene1)
      colnames(temp2) <- colnames(gene2)
      write.csv(temp1, file1_name, col.names = TRUE,row.names=FALSE)
      write.csv(temp2, file2_name, col.names = TRUE,row.names=FALSE)
    }else
    {
      current_file_name <- paste(paste(paste(paste("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/",species[j],sep=""),'__v__',sep=""), species[i],sep=""),".txt",sep="")
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[j],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("E:/GO_layers/final_photo/circadian_gene_tissue/",species[i],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      species1_new_gene_list <- list()
      species2_new_gene_list <- list()
      for(k in 1:ncol(gene1))
      {
        if(which(species_sort == species[i]) > which(species_sort == species[j]))
        {
          current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
          species1_new_gene_list[[k]] <- unique(as.character(current_gene1))
          
        }else
        {
          current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
          species1_new_gene_list[[k]] <- unique(as.character(current_gene1))
        }
      }
      for(l in 1:ncol(gene2))
      {
        if(which(species_sort == species[i]) > which(species_sort == species[j]))
        {
          current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
          species2_new_gene_list[[l]] <- unique(as.character(current_gene2))
          
        }else
        {
          current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
          species2_new_gene_list[[l]] <- unique(as.character(current_gene2))
        }
      }
      temp1 <- do.call(cbind, lapply(lapply(species1_new_gene_list, unlist), `length<-`, max(lengths(species1_new_gene_list))))
      temp2 <- do.call(cbind, lapply(lapply(species2_new_gene_list, unlist), `length<-`, max(lengths(species2_new_gene_list))))
      file1_name <- paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],'_vs_',species[j],'_',species[j],".csv",sep="")
      file2_name <- paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],'_vs_',species[j],'_',species[i],".csv",sep="")
      colnames(temp1) <- colnames(gene1)
      colnames(temp2) <- colnames(gene2)
      write.csv(temp1, file1_name, col.names = TRUE,row.names=FALSE)
      write.csv(temp2, file2_name, col.names = TRUE,row.names=FALSE)
    }
  }
}




for(i in 1:(length(species) - 1))
{
  for(j in (i+1): length(species))
  {
    current_file_name <- paste(paste(paste(paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],sep=""),'_vs_',sep=""), species[j],'_',species[i],sep=""),".csv",sep="")
    if (file.exists(current_file_name))
    {
      enrich_gene1 <- read.csv(current_file_name)
      current_file_name2 <- paste(paste(paste(paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],sep=""),'_vs_',sep=""), species[j],'_',species[j],sep=""),".csv",sep="")
      enrich_gene2 <- read.csv(current_file_name2)
      gene_id <- list()
      for(n in 1:ncol(enrich_gene1))
      {
        gene_id[[n]] <- bitr(as.character(enrich_gene1[,n]), fromType="SYMBOL", toType="ENTREZID", OrgDb=go_list[[species[i]]])$ENTREZID
      }
      temp1 <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
      colnames(temp1) <- colnames(enrich_gene1)
      
      gene_id <- list()
      for(n in 1:ncol(enrich_gene2))
      {
        gene_id[[n]] <- bitr(as.character(enrich_gene2[,n]), fromType="SYMBOL", toType="ENTREZID", OrgDb=go_list[[species[j]]])$ENTREZID
      }
      temp2 <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
      colnames(temp2) <- colnames(enrich_gene2)
      
      ck11 <- compareCluster(geneClusters = as.data.frame(temp1),fun = "enrichGO",OrgDb=go_list[[species[i]]],ont="BP",readable=TRUE,universe=as.character(expr_list[[species[i]]]))
      ck12 <- compareCluster(geneClusters = as.data.frame(temp2),fun = "enrichGO",OrgDb=go_list[[species[j]]],ont="BP",readable=TRUE,universe=as.character(expr_list[[species[j]]]))
      ck12 <- enrichGO(as.character(temp1[,1]),OrgDb=go_list[[species[i]]],ont="BP",readable=TRUE,universe=as.character(expr_list[[species[i]]]))
      write.csv(c)
    }
  }
}









######################################################adjust phase photo ordered by pvalue





file_name <- list.files(path = "E:/GO_layers/final_photo/exp_amp_data/mouse/", pattern = NULL, all.files = FALSE,full.names = FALSE, recursive = T,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
for(k in 1:length(file_name))
{
  current_file <- paste("E:/GO_layers/final_photo/exp_amp_data/mouse/", file_name[k], sep="")
  data <- read.csv(current_file)
  data <- na.omit(data)
  for(i in 1:nrow(data))
  {
    current_phase <- data[i,]$LAG
    current_phase <- current_phase %% 24
    diff_phase <- floor(current_phase / 2)
    day1 = data[i, 8:19]
    new_day1 <- day1
    for(j in 1:length(day1))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 12
      }
      new_day1[num] <- day1[j]
    }
    data[i, 8:19] <- new_day1
    
    day2 = data[i, 20:31]
    new_day2 <- day2
    for(j in 1:length(day2))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 12
      }
      new_day2[num] <- day2[j]
    }
    data[i, 20:31] <- new_day2
    
  }
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/mouse/",file_name[k],sep="")
  write.csv(data, file_name_new)
  # data <- data[,6:29]
  # rt2<-apply(data,1,scale)
  # r<-t(rt2)
  # pheatmap(r,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  
  
}















for(i in 1:nrow(data))
{
  current_phase <- data[i,]$JTK_adjphase
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 20:31]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase * 2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day1[num] <- day1[j]
  }
  data[i, 20:31] <- new_day1
  
  day2 = data[i, 32:43]
  new_day2 <- day2
  for(j in 1:length(day2))
  {
    num = j - diff_phase * 2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day2[num] <- day2[j]
  }
  data[i, 32:43] <- new_day2
  
  
  day3 = data[i, 44:55]
  new_day3 <- day3
  for(j in 1:length(day3))
  {
    num = j - diff_phase * 2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day3[num] <- day3[j]
  }
  data[i, 44:55] <- new_day3
  
}
write.csv(data, "E:/new.csv")
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/mouse/Adr_Jtk.csv")
data <- rt[,8:31]
rt2<-apply(data,1,scale) ##按行scale
#write.csv(rt2,"D:/cyclingdata/archaea Halobacterium salinarum NRC-1/A/古菌/cycling6.csv")
r<-t(rt2)
p001 <- which(rt$ADJ.P < 0.01 & rt$ADJ.P >0.009)
p001 <- p001[length(p001)]
p005 <- which(rt$ADJ.P < 0.05 & rt$ADJ.P >0.04)
p005 <- p005[length(p005)]
p01 <- which(rt$ADJ.P < 0.1 & rt$ADJ.P >0.09)
p01 <- p01[length(p01)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#archaea
data <- read.csv("C:/Users/wanglab/Downloads/meta2d_AR2_achaea.csv")
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$JTK_adjphase
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 20:31]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase * 2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day1[num] <- day1[j]
  }
  data[i, 20:31] <- new_day1
  
  day2 = data[i, 32:43]
  new_day2 <- day2
  for(j in 1:length(day2))
  {
    num = j - diff_phase * 2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day2[num] <- day2[j]
  }
  data[i, 32:43] <- new_day2
  
  
  day3 = data[i, 44:55]
  new_day3 <- day3
  for(j in 1:length(day3))
  {
    num = j - diff_phase * 2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day3[num] <- day3[j]
  }
  data[i, 44:55] <- new_day3
  
}
write.csv(data, "E:/new.csv")
rt<-read.csv("E:/new.csv")
data <- rt[,20:55]
pvalue <- rt[,2]
rt2<-apply(data,1,scale) ##按行scale
#write.csv(rt2,"D:/cyclingdata/archaea Halobacterium salinarum NRC-1/A/古菌/cycling6.csv")
r<-t(rt2)
a <- pheatmap(r,cluster_rows = FALSE,cluster_cols = pvalue,show_colnames = TRUE,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))



























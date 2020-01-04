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
upset(data,nsets = 14,keep.order=TRUE)



number <- c(4607,4510,445,316,1055,156,191,787,149,2059, 2631,1949,275,286,93, 62, 205,37,49,152,34,426,514,375,60,61,105,66,202,36,42,171,37,415,530,352,55,66)
names(number) <- c("OLD","YOUNG","ADR","AOR","BFAT","BST","CERE","HEA","HYPO","KID","LIV","LUN","MUS","WAT","ADR&OLD","AOR&OLD","BFAT&OLD","BST&OLD","CERE&OLD","HEA&OLD","HYPO&OLD","KID&OLD","LIV&OLD","LUN&OLD","MUS&OLD","WAT&OLD","ADR&OLD","AOR&YOUNG","BFAT&YOUNG","BST&YOUNG","CERE&YOUNG","HEA&YOUNG","HYPO&YOUNG","KID&YOUNG","LIV&YOUNG","LUN&YOUNG","MUS&YOUNG","WAT&YOUNG")
data <- fromExpression(number)
upset(data,nsets = 14,keep.order=TRUE)






arabid_gene <- read.csv("E:/GO_layers/circadian_data/拟南芥/circadian_gene_arabidopsis.csv")
arabid_anno <- read.csv("E:/GO_layers/circadian_data/拟南芥/gse5612_annotation.csv")
arabid_symbol <- arabid_anno[arabid_anno$ID %in% arabid_gene[,1], ]$Gene.symbol
write.csv(arabid_symbol,"E:/GO_layers/circadian_data/拟南芥/symbol.csv")



human_cell_gene <- read.csv("E:/GO_layers/circadian_data/human/4tissue/cell/circadian_gene.csv")
human_anno <- read.csv("E:/GO_layers/circadian_data/human/4tissue/cell/2.csv")
gene_name <- human_cell_gene$cycID













#################################################amp and exp new
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

# Save the file
dev.off()








#data
tissuename <- c("Adr","Aorta","BFAT","BS","Cere","Heart","Hypo","Kidney","Liver","Lung","Mus","WFAT")
par(mfcol=c(4,3))
aa= matrix(0,nrow = 12,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
bb=1:12

for(i in 1:12){
  expr = read.csv(paste0("E:/GO_layers/circadian_data/mouse/allgenes12tissues/",tissuename[i],'_Jtk.csv'))
  cir = expr[,'BH.Q']<0.05
  expr_cir = expr[cir,]
  expr_cir <- expr_cir[expr_cir[,'AMP']!=0,]
  expr_cirmean= rowMeans(expr_cir[,8:31])
  
  
  
  
  
  plot(log(expr_cirmean),log(expr_cir[,'AMP']),xlab = '',ylab = '',main=tissuename[i])
  aa[i,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$estimate
  aa[i,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'AMP'])))$p.value
  aa[i,3]= aa[i,1]^2
  bb[i] = sum(cir)
}
#write.csv(aa,"correlation_between_amp_and_expr.csv")

{
  #aa= read.csv("correlation_between_amp_and_expr.csv",row.names = 1)
  aaa= matrix(0,nrow = 12,ncol = 2)
  rownames(aaa)= tissuename
  colnames(aaa)=c('expression','others')
  aaa[,1]= aa[,'R2']
  aaa[,2]= 1- aa[,'R2']
  
  par(mfcol=c(4,3))
  parmar= par()$mar
  par(mar = c(0,0, 2, 0),las = 1)
  for(i in 1:12){
    pie(aaa[i,],labels = round(aaa[i,]*100)/100,main = tissuename[i],col= c('skyblue','white'),init.angle = 90)
  }
  par(mar=parmar,las= 0)
  
  par(mfcol=c(1,1))
  plot.new()
  legend('top',c('expression','others'),fill = c('skyblue','white'),cex=3)
}
tissuenamefull= c("Adrenal Gland","Aorta","Brown Fat","Brainstem","Cerebellum","Heart","Hypothalamus",
                  "Kidney","Liver","Lung","Skeletal Muscle","White Fat")


{
  #aa= read.csv("correlation_between_amp_and_expr.csv",row.names = 1)
  aaa= matrix(0,nrow = 12,ncol = 2)
  rownames(aaa)= tissuename
  colnames(aaa)=c('expression','others')
  aaa[,1]= aa[,'R2']
  aaa[,2]= 1- aa[,'R2']
  # png("fig2.png",width=10,height=12,units="in",res=350)
  # par(mfrow=c(4,4))
  par(plt=c(0.15,0.85,0.25,0.65))
  par(mai=c(0.6732,0.5412,0.5412,0.2772))
  
  {
    mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
    widths <-c(9,9,9,9)  
    heights <-c(6,6,6,2)  
    layout(mat,widths = widths, heights = heights)  
    nf <-layout(mat,  widths = widths, heights = heights)  
    layout.show(nf)
  }
  
  for(i in 1:12){
    expr = read.csv(paste0("E:/GO_layers/circadian_data/mouse/allgenes12tissues/",tissuename[i],'_Jtk.csv'))
    cir = expr[,'BH.Q']<0.05
    expr_cir = expr[cir,]
    expr_cirmean= rowMeans(expr_cir[,7:30])
    plot(log(expr_cirmean),log(expr_cir[,'AMP']),
         xlab = 'expression',ylab = 'amplitude',main=tissuenamefull[i],
         xlim = c(2,10),ylim = c(-2,9),
         pch=20,cex=0.5)
    text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)
    
    par(plt=c(0.60,0.95,0.20,0.55))
    par(new=T)
    pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
    text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
    par(plt=c(0.15,0.85,0.25,0.65))
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
  }
  {
    # parmai=par()$mai
    # par(mai=rep(0,4))
    # plot.new()
    # plot.new()
    # legend('top',c('variation explained by expression','variation explained by other factors'),
    #        fill = c('skyblue','white'),cex=1.5,bty='n')
    # par(mai=parmai)
  }
  
  # dev.off()
  
  {
    parraw_mar= par()$mar
    par(mar = c(0,0, 0, 0),las = 1)  
    plot(c(0, 4),c(0, 4),  
         xlim = c(0, 4), ylim = c(0, 4),  
         xaxs = 'i', yaxs = 'i',  
         xaxt = 'n', yaxt = 'n', type = 'n', ann =F, axes = F) 
    legend('top',c('variation explained by expression','variation explained by other factors'),
           fill = c('skyblue','white'),cex=1.5,bty='n')
    par(mar=parraw_mar,las= 0)
  }
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
write.csv(venn_data_new,"E:/GO_layers/final_photo/venn_data.csv")
venn_data_new <- read.csv("E:/GO_layers/final_photo/venn_data.csv")
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
  temp_plot <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, ext.text=FALSE, height = 450, width = 450, rotation.degree = number_flag, resolution =300, imagetype="png", col="black", fill=c(colors()[38], 'steelblue2'), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
  plot_list <- c(plot_list, list(temp_plot))
}
grid.newpage()


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
species11 <- c('H.sap','M.mu','D.me','S.ce','N.cr','A.th','Syn','H.sal')
species12 <- c('H.sal','Syn','A.th','N.cr','S.ce','D.me','M.mu','H.sap')
print(as.ggplot(expression(plot(1:8,1:8,type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:8,labels = species11) + axis(side = 4,at = 1:8,labels = species12))), vp=vp.scatter)
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
neurospora_expr_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/neurospora_gene.csv")
neurospora_expr_gene <- neurospora_expr_gene$CycID
neuro_ortho <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/neurospora__v__yeast.txt")
neurospora_expr_gene <- neuro_ortho[neuro_ortho[,1] %in% neurospora_expr_gene,2]
cyanobacteria_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/cyanobacteria_gene.csv")
cyanobacteria_expr_gene <- cyanobacteria_gene$X7942_ID
cyano_ortho <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/cyanobacteria__v__yeast.txt")
cyanobacteria_expr_gene <- cyano_ortho[cyano_ortho[,1] %in% cyanobacteria_expr_gene,2]
archaea_gene <- read.csv("E:/GO_layers/final_photo/expression_gene/archaea_gene.csv")
archaea_expr_gene <- archaea_gene$CycID
arch_ortho <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/archaea__v__yeast.txt")
archaea_expr_gene <- arch_ortho[arch_ortho[,1] %in% archaea_expr_gene, 2]

expr_list <- list(mouse_expr_gene,human_expr_gene,fly_expr_gene,yeast_expr_gene,arabid_expr_gene,neurospora_expr_gene,cyanobacteria_expr_gene,archaea_expr_gene)
names(expr_list) <- c("mouse","human","fly","yeast","arabidopsis","neurospora","cyanobacteria","archaea")
go_list <- list(org.Mm.eg.db, org.Hs.eg.db, org.Dm.eg.db, org.Sc.sgd.db, org.At.tair.db,org.Sc.sgd.db,org.Sc.sgd.db,org.Sc.sgd.db)
names(go_list) <- c("mouse","human","fly","yeast","arabidopsis","neurospora","cyanobacteria","archaea")


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
      print(species[i])
      print(species[j])
      enrich_gene1 <- read.csv(current_file_name)
      current_file_name2 <- paste(paste(paste(paste("E:/GO_layers/final_photo/enrich_orthology/",species[i],sep=""),'_vs_',sep=""), species[j],'_',species[j],sep=""),".csv",sep="")
      enrich_gene2 <- read.csv(current_file_name2)
      gene_id <- list()
      
      
      # if(species[i] == "neurospora")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/neurospora__v__yeast.txt")
      #   enrich_gene1 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[i] == "cyanobacteria")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/cyanobacteria__v__yeast.txt")
      #   enrich_gene1 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[i] == "archaea")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/archaea__v__yeast.txt")
      #   enrich_gene1 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # 
      # 
      # if(species[j] == "neurospora")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/neurospora__v__yeast.txt")
      #   enrich_gene2 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[j] == "cyanobacteria")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/cyanobacteria__v__yeast.txt")
      #   enrich_gene2 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[j] == "archaea")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/archaea__v__yeast.txt")
      #   enrich_gene2 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      
      if(species[i] == "yeast" | species[i] == "neurospora" | species[i] == "cyanobacteria"|species[i] == "archaea")
      {
        current_fromtype1 <- "ORF"
      }else
      {
        current_fromtype1 <- "SYMBOL"
      }
      if(species[j] == "yeast" | species[j] == "neurospora" | species[j] == "cyanobacteria"|species[j] == "archaea")
      {
        current_fromtype2 <- "ORF"
      }else
      {
        current_fromtype2 <- "SYMBOL"
      }
      
      if(species[i] != "arabidopsis" & species[i] !="neurospora" & species[i] !="cyanobacteria" &species[i] !="archaea")
      {
        gene_id <- list()
        for(n in 1:ncol(enrich_gene1))
        {
          gene_id[[n]] <- bitr(as.character(enrich_gene1[,n]), fromType=current_fromtype1, toType="ENTREZID", OrgDb=go_list[[species[i]]])$ENTREZID
        }
        temp1 <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
        colnames(temp1) <- colnames(enrich_gene1)
        enrich_gene1 <- temp1
      }
      if(species[j] != "arabidopsis" & species[j] !="neurospora" & species[j] !="cyanobacteria" &species[j] !="archaea")
      {
        gene_id <- list()
        for(n in 1:ncol(enrich_gene2))
        {
          gene_id[[n]] <- bitr(as.character(enrich_gene2[,n]), fromType=current_fromtype2, toType="ENTREZID", OrgDb=go_list[[species[j]]])$ENTREZID
        }
        temp2 <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
        colnames(temp2) <- colnames(enrich_gene2)
        enrich_gene2 <- temp2
      }
      if(species[i] == "arabidopsis")
      {
        current_type1 = "TAIR"
      }else
      {
        current_type1 = "ENTREZID"
      }
      
      if(species[j] == "arabidopsis")
      {
        current_type2 = "TAIR"
      }else
      {
        current_type2 = "ENTREZID"
      }
      

      

      
      if(species[i] !="neurospora" & species[i] !="cyanobacteria" &species[i] !="archaea")
      {
        ck11 <- compareCluster(geneClusters = as.data.frame(enrich_gene1),fun = "enrichGO",keyType = current_type1,OrgDb=go_list[[species[i]]],ont="BP",universe=as.character(expr_list[[species[i]]]))
        species1_filename <- paste("E:/GO_layers/final_photo/enrich_result/",species[i],"_vs_",species[j],"_",species[i],"_bp",".csv",sep="")
        write.csv(ck11@compareClusterResult, species1_filename)
      }
      if(species[j] !="neurospora" & species[j] !="cyanobacteria" &species[j] !="archaea")
      {
        ck12 <- compareCluster(geneClusters = as.data.frame(enrich_gene2),fun = "enrichGO",keyType = current_type2,OrgDb=go_list[[species[j]]],ont="BP",universe=as.character(expr_list[[species[j]]]))
        species2_filename <- paste("E:/GO_layers/final_photo/enrich_result/",species[i],"_vs_",species[j],"_",species[j],"_bp",".csv",sep="")
        write.csv(ck12@compareClusterResult, species2_filename)
      }
      
      
    
      
      
    }else
    {
      current_file_name <- paste(paste(paste(paste("E:/GO_layers/final_photo/enrich_orthology/",species[j],sep=""),'_vs_',sep=""), species[i],'_',species[j],sep=""),".csv",sep="")
      print(species[i])
      print(species[j])
      enrich_gene1 <- read.csv(current_file_name)
      current_file_name2 <- paste(paste(paste(paste("E:/GO_layers/final_photo/enrich_orthology/",species[j],sep=""),'_vs_',sep=""), species[i],'_',species[i],sep=""),".csv",sep="")
      enrich_gene2 <- read.csv(current_file_name2)
      gene_id <- list()
      
      
      # if(species[j] == "neurospora")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/neurospora__v__yeast.txt")
      #   enrich_gene1 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[j] == "cyanobacteria")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/cyanobacteria__v__yeast.txt")
      #   enrich_gene1 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[j] == "archaea")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/archaea__v__yeast.txt")
      #   enrich_gene1 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      
      
      # if(species[i] == "neurospora")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/neurospora__v__yeast.txt")
      #   enrich_gene2 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[i] == "cyanobacteria")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/cyanobacteria__v__yeast.txt")
      #   enrich_gene2 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      # if(species[i] == "archaea")
      # {
      #   current_orthology_list <- read.table("E:/GO_layers/final_photo/orthoFinder_result/protein2gene/archaea__v__yeast.txt")
      #   enrich_gene2 <- current_orthology_list[current_orthology_list[,1] %in% enrich_gene1, 2]
      # }
      
      if(species[j] == "yeast" | species[j] == "neurospora" | species[j] == "cyanobacteria"|species[j] == "archaea")
      {
        current_fromtype1 <- "ORF"
      }else
      {
        current_fromtype1 <- "SYMBOL"
      }
      if(species[i] == "yeast" | species[i] == "neurospora" | species[i] == "cyanobacteria"|species[i] == "archaea")
      {
        current_fromtype2 <- "ORF"
      }else
      {
        current_fromtype2 <- "SYMBOL"
      }
      
      if(species[j] != "arabidopsis" & species[j] !="neurospora" & species[j] !="cyanobacteria" &species[j] !="archaea")
      {
        gene_id <- list()
        for(n in 1:ncol(enrich_gene1))
        {
          gene_id[[n]] <- bitr(as.character(enrich_gene1[,n]), fromType=current_fromtype1, toType="ENTREZID", OrgDb=go_list[[species[j]]])$ENTREZID
        }
        temp1 <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
        colnames(temp1) <- colnames(enrich_gene1)
        enrich_gene1 <- temp1
      }
      if(species[i] != "arabidopsis" & species[i] !="neurospora" & species[i] !="cyanobacteria" &species[i] !="archaea")
      {
        gene_id <- list()
        for(n in 1:ncol(enrich_gene2))
        {
          gene_id[[n]] <- bitr(as.character(enrich_gene2[,n]), fromType=current_fromtype2, toType="ENTREZID", OrgDb=go_list[[species[i]]])$ENTREZID
        }
        temp2 <- do.call(cbind, lapply(lapply(gene_id, unlist), `length<-`, max(lengths(gene_id))))
        colnames(temp2) <- colnames(enrich_gene2)
        enrich_gene2 <- temp2
      }
      if(species[j] == "arabidopsis")
      {
        current_type1 = "TAIR"
      }else
      {
        current_type1 = "ENTREZID"
      }
      
      if(species[i] == "arabidopsis")
      {
        current_type2 = "TAIR"
      }else
      {
        current_type2 = "ENTREZID"
      }
      
      
      if(species[j] !="neurospora" & species[j] !="cyanobacteria" &species[j] !="archaea")
      {
        ck11 <- compareCluster(geneClusters = as.data.frame(enrich_gene1),fun = "enrichGO",keyType = current_type1,OrgDb=go_list[[species[j]]],ont="BP",universe=as.character(expr_list[[species[j]]]))
        species1_filename <- paste("E:/GO_layers/final_photo/enrich_result/",species[j],"_vs_",species[i],"_",species[j],"_bp",".csv",sep="")
        write.csv(ck11@compareClusterResult, species1_filename)
      }

      if(species[i] !="neurospora" & species[i] !="cyanobacteria" &species[i] !="archaea")
      {
        ck12 <- compareCluster(geneClusters = as.data.frame(enrich_gene2),fun = "enrichGO",keyType = current_type2,OrgDb=go_list[[species[i]]],ont="BP",universe=as.character(expr_list[[species[i]]]))
        species2_filename <- paste("E:/GO_layers/final_photo/enrich_result/",species[j],"_vs_",species[i],"_",species[i],"_bp",".csv",sep="")
        write.csv(ck12@compareClusterResult, species2_filename)
      }
      
    }
  
  }
}


#upset

number <- c(3.7096031,3.6762335,4.074952,2.075227,3.8537829,4.0515654,3.0769231,9.4161361,14.4481356,3.8220919, 12.6066947,12.7404447,2.9028436,13.4369038,4.6260601,12.2164868, 11.5183246, 4.0160643,5.0022036,19.7991392,31.2217195,5.3475936,5.6559076,3.9568345,4.4117647,18.4515531,4.0332147,2.7027027)
names(number) <- c("human&mouse","human&fly","human&yeast","human&neurospora","human&arabidopsis","human&cyanobacteria","human&archaea","mouse&fly","mouse&yeast","mouse&neurospora","mouse&arabidopsis","mouse&cyanobacteria","mouse&archaea","fly&yeast","fly&neurospora","fly&arabidopsis","fly&cyanobacteria","fly&archaea","yeast&neurospora","yeast&arabidopsis","yeast&cyanobacteria","yeast&archaea","neurospora&arabidopsis","neurospora&cyanobacteria","neurospora&archaea","arabidopsis&cyanobacteria","arabidopsis&archaea","archaea&cyanobacteria")
data <- fromExpression(number)
upset(data,nsets = 28,keep.order=TRUE,order.by="degree")







##############################read data and draw term photo

file_dir <- "E:/GO_layers/final_photo/enrich_result/"
species <- c("human","mouse","fly","yeast","neurospora","arabidopsis","cyanobacteria","archaea")
term_venn_data <- data.frame(0,0,0, 0, 0)
for(i in 1:(length(species)-1))
{
  for(j in (i + 1):length(species))
  {
    current_file1_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[i],"_bp.csv",sep="")
    if (file.exists(current_file1_name))
    {
      current_file2_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[j],"_bp.csv",sep="")
      current_file1 <- read.csv(current_file1_name)
      current_file2 <- read.csv(current_file2_name)
      species1_term_num <- 0
      species2_term_num <- 0
      overlap_num <- 0
      num <- 0
      cluster1 <- unique(current_file1$Cluster)
      cluster2 <- unique(current_file2$Cluster)
      for(k in 1:length(cluster1))
      {
        for(l in 1:length(cluster2))
        {
          current_data1 <- current_file1[current_file1$Cluster %in% cluster1[k],]
          current_data2 <- current_file2[current_file2$Cluster %in% cluster2[l],]
          num <- num + 1
          species1_term_num <- species1_term_num + length(unique(current_data1$ID))
          species2_term_num <- species2_term_num + length(unique(current_data2$ID))
          overlap_num<- overlap_num + length(intersect(current_data1$ID,current_data2$ID))
          
          
        }
      }
      temp_row <- c(species[i], species[j], species1_term_num/num, species2_term_num/num, overlap_num/num)
      term_venn_data <- rbind(term_venn_data, temp_row)
    }else
    {
      current_file1_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[j],"_bp.csv",sep="")
      current_file2_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[i],"_bp.csv",sep="")
      current_file1 <- read.csv(current_file1_name)
      current_file2 <- read.csv(current_file2_name)
      species1_term_num <- 0
      species2_term_num <- 0
      overlap_num <- 0
      num <- 0
      cluster1 <- unique(current_file1$Cluster)
      cluster2 <- unique(current_file2$Cluster)
      for(k in 1:length(cluster1))
      {
        for(l in 1:length(cluster2))
        {
          current_data1 <- current_file1[current_file1$Cluster %in% cluster1[k],]
          current_data2 <- current_file2[current_file2$Cluster %in% cluster2[l],]
          num <- num + 1
          species1_term_num <- species1_term_num + length(unique(current_data1$ID))
          species2_term_num <- species2_term_num + length(unique(current_data2$ID))
          overlap_num<- overlap_num + length(intersect(current_data1$ID,current_data2$ID))
          
          
        }
      }
      temp_row <- c(species[j], species[i], species1_term_num/num, species2_term_num/num, overlap_num/num)
      term_venn_data <- rbind(term_venn_data, temp_row)
    }
  }
}

write.csv(term_venn_data,"E:/GO_layers/final_photo/venn_term_data.csv")



species <- c("human","mouse","fly","yeast","neurospora","arabidopsis","cyanobacteria","archaea")
venn_data_new <- read.csv("E:/GO_layers/final_photo/venn_term_data.csv")
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
  temp_plot <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, ext.text=FALSE, height = 450, width = 450, rotation.degree = number_flag, resolution =300, imagetype="png", col="black", fill=c(colors()[38], 'steelblue2'), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
  plot_list <- c(plot_list, list(temp_plot))
}
grid.newpage()


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
species11 <- c('H.sap','M.mu','D.me','S.ce','N.cr','A.th','Syn','H.sal')
species12 <- c('H.sal','Syn','A.th','N.cr','S.ce','D.me','M.mu','H.sap')
print(as.ggplot(expression(plot(1:8,1:8,type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:8,labels = species11) + axis(side = 4,at = 1:8,labels = species12))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=8, ncol=8, axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)



#upset term
number <- c(3.020008,5.9047619,2.4054983,1.9417476,5.4669704,5.8641975,0,5.8521999,11.4348143,11.952191, 10.22499,11.770073,0,24.3325706,3.2634033,14.6464646, 25.0591017, 1.910828,2.4604569,14.6226415,25.2144082,1,1,16.1073826,15.0943396,11.4774115,1.07,4.1666667)
names(number) <- c("human&mouse","human&fly","human&yeast","human&neurospora","human&arabidopsis","human&cyanobacteria","human&archaea","mouse&fly","mouse&yeast","mouse&neurospora","mouse&arabidopsis","mouse&cyanobacteria","mouse&archaea","fly&yeast","fly&neurospora","fly&arabidopsis","fly&cyanobacteria","fly&archaea","yeast&neurospora","yeast&arabidopsis","yeast&cyanobacteria","yeast&archaea","neurospora&arabidopsis","neurospora&cyanobacteria","neurospora&archaea","arabidopsis&cyanobacteria","arabidopsis&archaea","archaea&cyanobacteria")
data <- fromExpression(number)
upset(data,nsets = 28)






number_gene <- c(3.7096031,3.6762335,4.074952,2.075227,3.8537829,4.0515654,3.0769231,9.4161361,14.4481356,3.8220919, 12.6066947,12.7404447,2.9028436,13.4369038,4.6260601,12.2164868, 11.5183246, 4.0160643,5.0022036,19.7991392,31.2217195,5.3475936,5.6559076,3.9568345,4.4117647,18.4515531,4.0332147,2.7027027)
number_term <- c(3.020008,5.9047619,2.4054983,1.9417476,5.4669704,5.8641975,0,5.8521999,11.4348143,11.952191, 10.22499,11.770073,0,24.3325706,3.2634033,14.6464646, 25.0591017, 1.910828,2.4604569,14.6226415,25.2144082,1,1,16.1073826,15.0943396,11.4774115,1.07,4.1666667)
dataset <- data.frame(value = c(number_gene, number_term), group = factor(rep(c("Overlap between circadian gene","Overlap between GO term"),each=28)))
boxplot(number_gene,number_term,ylab="Percentage",border=c(colors()[38], 'steelblue2'))



#barplot
data <- read.csv("E:/GO_layers/final_photo/photo/new1.csv")
for(i in 1:nrow(data))
{
  temp <- data[i,]
  temp <- reshape2::melt(temp)
  temp_name <- paste("E:/GO_layers/final_photo/photo/percent/plot",i,".pdf",sep="")
  pdf(temp_name, width=5.22, height=5.22)
  barplot(temp$value,horiz=TRUE,col="mediumblue")
  dev.off()
  #a <- ggplot(temp,aes(x=variable, y=value)) + geom_col(width=0.5,fill="mediumblue")  + scale_y_continuous(expand = c(0, 0))+ coord_flip() + theme_bw()
  #ggsave(a, file=temp_name, width=5.22, height=5.22)

}


##########################same GO
# file_dir <- "E:/GO_layers/final_photo/enrich_result/"
# species <- c("human","mouse","fly","yeast","neurospora","arabidopsis","cyanobacteria","archaea")
# term_venn_data <- data.frame(0,0,0, 0, 0)
# for(i in 1:(length(species)-1))
# {
#   for(j in (i + 1):length(species))
#   {
#     current_file1_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[i],"_bp.csv",sep="")
#     if (file.exists(current_file1_name))
#     {
#       current_file2_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[j],"_bp.csv",sep="")
#       current_file1 <- read.csv(current_file1_name)
#       current_file2 <- read.csv(current_file2_name)
#       if(i==1 & j==1)
#       {
#         same_go <- c(as.character(unique(current_file1$ID)))
#         same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
#       }else
#       {
#         same_go <- c(same_go, c(as.character(unique(current_file1$ID))))
#         same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
#       }
#     }else
#     {
#       current_file1_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[j],"_bp.csv",sep="")
#       current_file2_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[i],"_bp.csv",sep="")
#       current_file1 <- read.csv(current_file1_name)
#       current_file2 <- read.csv(current_file2_name)
#       
#       if(i==1 & j==1)
#       {
#         same_go <- c(as.character(unique(current_file1$ID)))
#         same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
#       }else
#       {
#         same_go <- c(same_go, c(as.character(unique(current_file1$ID))))
#         same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
#       }
#     }
#   }
# }
# write.csv(table(same_go),"E:/GO_layers/final_photo/samego.csv")



file_dir <- "E:/GO_layers/final_photo/enrich_result/"
species <- c("human","mouse","fly","yeast","neurospora","arabidopsis","cyanobacteria","archaea")
term_venn_data <- data.frame(0,0,0, 0, 0)
num <- 0
for(i in 1:(length(species)-1))
{
  for(j in (i + 1):length(species))
  {
    
    current_file1_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[i],"_bp.csv",sep="")
    if (file.exists(current_file1_name))
    {
      current_file2_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[j],"_bp.csv",sep="")
      current_file1 <- read.csv(current_file1_name)
      current_file2 <- read.csv(current_file2_name)
      tissue1 <- unique(current_file1$Cluster)
      tissue2 <- unique(current_file2$Cluster)
      for(k in 1:length(tissue1))
      {
        for(l in 1:length(tissue2))
        {
          num <- num + 1
          current_data1 <- current_file1[current_file1$Cluster %in% tissue1[k],]
          current_data2 <- current_file2[current_file2$Cluster %in% tissue2[l],]
          if(num == 1 & k == 1 & l == 1)
          {
            same_go <- c(as.character(unique(current_file1$ID)))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
          }else
          {
            same_go <- c(same_go, c(as.character(unique(current_file1$ID))))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
          }
        }
      }

    }else
    {
      current_file1_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[j],"_bp.csv",sep="")
      current_file2_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[i],"_bp.csv",sep="")
      current_file1 <- read.csv(current_file1_name)
      current_file2 <- read.csv(current_file2_name)
      tissue1 <- unique(current_file1$Cluster)
      tissue2 <- unique(current_file2$Cluster)
      for(k in 1:length(tissue1))
      {
        for(l in 1:length(tissue2))
        {
          num <- num + 1
          current_data1 <- current_file1[current_file1$Cluster %in% tissue1[k],]
          current_data2 <- current_file2[current_file2$Cluster %in% tissue2[l],]
          if(num == 1 & k == 1 & l == 1)
          {
            same_go <- c(as.character(unique(current_file1$ID)))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
          }else
          {
            same_go <- c(same_go, c(as.character(unique(current_file1$ID))))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
          }
        }
      }
    }
  }
}
write.csv(table(same_go),"E:/GO_layers/final_photo/samego_new.csv")












# go_data <- read.csv("E:/GO_layers/final_photo/topgo.csv")
# 
# ggplot(go_data, aes(x=a,y=ID,size=Number)) + 
#   geom_dotplot(binaxis='y', stackdir='center')
# 
# dotplot(go_data)





data <- read.csv("E:/GO_layers/final_photo/samego_new.csv")
data <- data[,-1]
data <- data[1:30,]
# barplot(data$Freq,names.arg = data$name,ylim=c(0,50))
# text(cex = 1,data$name,srt = 45)
ggplot(data, aes(y=Freq,x=reorder(name,X=Freq))) + geom_col(fill = "royalblue2") + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),panel.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), panel.grid.major = element_blank(),plot.background = element_rect(fill = "transparent",colour = NA)) + coord_flip()











##########################same GO new




file_dir <- "E:/GO_layers/final_photo/enrich_result/"
species <- c("human","mouse","fly","yeast","neurospora","arabidopsis","cyanobacteria","archaea")
term_venn_data <- data.frame(0,0,0, 0, 0)
num <- 0

for(i in 1:(length(species)-1))
{
  for(j in (i + 1):length(species))
  {
    
    current_file1_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[i],"_bp.csv",sep="")
    if (file.exists(current_file1_name))
    {
      current_file2_name <- paste(file_dir, species[i],"_vs_",species[j],"_",species[j],"_bp.csv",sep="")
      current_file1 <- read.csv(current_file1_name)
      current_file2 <- read.csv(current_file2_name)
      tissue1 <- unique(current_file1$Cluster)
      tissue2 <- unique(current_file2$Cluster)
      for(k in 1:length(tissue1))
      {
        for(l in 1:length(tissue2))
        {
          num <- num + 1
          current_data1 <- current_file1[current_file1$Cluster %in% tissue1[k],]
          current_data2 <- current_file2[current_file2$Cluster %in% tissue2[l],]
          if(num == 1 & k == 1 & l == 1)
          {
            same_go <- c(as.character(unique(current_file1$ID)))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
            expect_data <- as.numeric(1000* (1 / length(tissue1)) * (1 / length(tissue2)))
            expect_list <- rep(expect_data, length(unique(current_file1$ID)))
            expect_list <- c(expect_list, rep(expect_data, length(unique(current_file2$ID))))
          }else
          {
            same_go <- c(same_go, c(as.character(unique(current_file1$ID))))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
            expect_data <- 1000* (1 / length(tissue1)) * (1 / length(tissue2))
            expect_list <- c(expect_list, rep(expect_data, length(unique(current_file1$ID))))
            expect_list <- c(expect_list, rep(expect_data, length(unique(current_file2$ID))))
          }
        }
      }
      
    }else
    {
      current_file1_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[j],"_bp.csv",sep="")
      current_file2_name <- paste(file_dir, species[j],"_vs_",species[i],"_",species[i],"_bp.csv",sep="")
      current_file1 <- read.csv(current_file1_name)
      current_file2 <- read.csv(current_file2_name)
      tissue1 <- unique(current_file1$Cluster)
      tissue2 <- unique(current_file2$Cluster)
      for(k in 1:length(tissue1))
      {
        for(l in 1:length(tissue2))
        {
          num <- num + 1
          current_data1 <- current_file1[current_file1$Cluster %in% tissue1[k],]
          current_data2 <- current_file2[current_file2$Cluster %in% tissue2[l],]
          if(num == 1 & k == 1 & l == 1)
          {
            same_go <- c(as.character(unique(current_file1$ID)))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
            expect_list <- rep(expect_data, length(unique(current_file1$ID)))
            expect_list <- c(expect_list, rep(expect_data, length(unique(current_file2$ID))))
          }else
          {
            same_go <- c(same_go, c(as.character(unique(current_file1$ID))))
            same_go <- c(same_go, c(as.character(unique(current_file2$ID))))
            expect_list <- c(expect_list, rep(expect_data, length(unique(current_file1$ID))))
            expect_list <- c(expect_list, rep(expect_data, length(unique(current_file2$ID))))
          }
        }
      }
    }
  }
}

new_data <- as.data.frame(cbind(same_go, expect_list))
term_id <- unique(new_data$same_go)
new_data_output <- data.frame(term_id, 0)
for(i in 1:length(term_id))
{
  idx <- which(new_data[,1] == term_id[i])
  new_data_output[i,2] <- sum(as.numeric(as.vector((new_data[idx,2]))))
}


write.csv(new_data_output,"E:/GO_layers/final_photo/samego_new_expected.csv")

write.csv(table(same_go),"E:/GO_layers/final_photo/samego_new.csv")












# go_data <- read.csv("E:/GO_layers/final_photo/topgo.csv")
# 
# ggplot(go_data, aes(x=a,y=ID,size=Number)) + 
#   geom_dotplot(binaxis='y', stackdir='center')
# 
# dotplot(go_data)





data <- read.csv("E:/GO_layers/final_photo/samego_new.csv")
data <- read.csv("E:/GO_layers/final_photo/samego_new_expected.csv")
data <- data[,-1]
data <- data[1:30,]
# barplot(data$Freq,names.arg = data$name,ylim=c(0,50))
# text(cex = 1,data$name,srt = 45)
ggplot(data, aes(y=Freq,x=reorder(name,X=Freq))) + geom_col(fill = "royalblue2") + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),panel.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(), panel.grid.major = element_blank(),plot.background = element_rect(fill = "transparent",colour = NA)) + coord_flip()





















######################################################adjust phase photo ordered by pvalue
#mouse
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

#############################mouse 4h
tissuename1 <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(k in 1:length(tissuename))
{
  file_name2 <- paste("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename1[k],".csv",sep="")
  data <- read.csv(file_name2)
  data <- na.omit(data)
  for(i in 1:nrow(data))
  {
    current_phase <- data[i,]$meta2d_phase
    current_phase <- current_phase %% 24
    diff_phase <- floor(current_phase / 4)
    day1 = data[i, 24:29]
    new_day1 <- day1
    for(j in 1:length(day1))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 6
      }
      new_day1[num] <- day1[j]
    }
    data[i, 24:29] <- new_day1
    
    day2 = data[i, 30:35]
    new_day2 <- day2
    for(j in 1:length(day2))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 6
      }
      new_day2[num] <- day2[j]
    }
    data[i, 30:35] <- new_day2
    
  }
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/mouse_4h/",tissuename1[k],".csv",sep="")
  write.csv(data, file_name_new)
  # data <- data[,6:29]
  # rt2<-apply(data,1,scale)
  # r<-t(rt2)
  # pheatmap(r,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  
  
}





rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/mouse/WFAT_Jtk.csv")
rt <- rt[order(rt$BH.Q),]
data <- rt[,9:32]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$BH.Q < 0.01 & rt$BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$BH.Q < 0.05 & rt$BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$BH.Q < 0.1 & rt$BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$BH.Q < 0.2 & rt$BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))









########################################

tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
for(j in 1:length(tissuename))
{
  {
  file_name1 <- paste("E:/GO_layers/final_photo/adjust_phase/mouse_4h/",tissuename[j],".csv",sep="")
  rt<-read.csv(file_name1)
  rt <- rt[order(rt$meta2d_BH.Q),]
  data <- rt[,25:36]
  rt2<-apply(data,1,scale) ##按行scale
  r<-t(rt2)
  p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
  p001 <- p001[length(p001)]
  p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
  p005 <- p005[length(p005)]
  p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
  p01 <- p01[length(p01)]
  p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
  p02 <- p02[length(p02)]
  pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  }
}














########################################yeast
library(pheatmap)

  current_file <- paste("E:/circidian_algorithm/result/yeast/meta2d/", "meta2d_JTK.Sample2.csv", sep="")
  data <- read.csv(current_file,check.names = FALSE)
  data <- na.omit(data)
  for(i in 1:nrow(data))
  {
    current_phase <- data$meta2d_phase[i]
    current_phase <- current_phase %% 130
    diff_phase <- floor(current_phase/13)
    day1 = data[i, 24:33]
    new_day1 <- day1
    for(j in 1:length(day1))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 10
      }
      new_day1[num] <- day1[j]
    }
    data[i, 24:33] <- new_day1
    
    day2 = data[i, 34:43]
    new_day2 <- day2
    for(j in 1:length(day2))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 10
      }
      new_day2[num] <- day2[j]
    }
    data[i, 34:43] <- new_day2
    
  }
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/yeast/","meta2d_JTK.Sample2.csv",sep="")
  write.csv(data, file_name_new)
  # data <- data[,6:29]
  # rt2<-apply(data,1,scale)
  # r<-t(rt2)
  # pheatmap(r,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  
  

#meta2d
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/yeast/","meta2d_JTK.Sample2.csv",sep="")
yeast_data <- read.csv(file_name_new)
rt <- yeast_data
data <- yeast_data[,24:43]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))

#jtk
yeast_data <- read.csv(file_name_new)
rt <- yeast_data
data <- yeast_data[,24:43]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/yeast/meta2d_JTK.Sample2.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/yeast/rain/sample2.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,24:43]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))








  current_file <- paste("E:/circidian_algorithm/result/yeast/meta2d/", "meta2d_JTK.Sample6.csv", sep="")
  data <- read.csv(current_file,check.names = FALSE)
  data <- na.omit(data)
  for(i in 1:nrow(data))
  {
    current_phase <- data$meta2d_phase[i]
    current_phase <- current_phase %% 432
    diff_phase <- floor(current_phase/36)
    day1 = data[i, 24:35]
    new_day1 <- day1
    for(j in 1:length(day1))
    {
      num = j - diff_phase + 3
      if(num <= 0)
      {
        num <- num + 12
      }
      new_day1[num] <- day1[j]
    }
    data[i, 24:35] <- new_day1
    
    day2 = data[i, 36:47]
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
    data[i, 36:47] <- new_day2
    
  }
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/yeast/","meta2d_JTK.Sample6.csv",sep="")
  write.csv(data, file_name_new)
  # data <- data[,6:29]
  # rt2<-apply(data,1,scale)
  # r<-t(rt2)
  # pheatmap(r,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  
  
#meta2d
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/yeast/","meta2d_JTK.Sample6.csv",sep="")
{
yeast_data <- read.csv(file_name_new)
rt <- yeast_data
data <- yeast_data[,25:44]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("darkblue", "darkblue","darkblue","purple","red", "red","red"))(20))

}

  
#jtk
yeast_data <- read.csv(file_name_new)
rt <- yeast_data
data <- yeast_data[,25:44]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("darkblue", "darkblue","darkblue","purple","red", "red","red"))(10))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/yeast/meta2d_JTK.Sample6.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/yeast/rain/sample6.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,24:43]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


  
  
  
  
  
  
  
  
  
  
  
  
  
#fly

  current_file <- paste("E:/GO_layers/final_photo/exp_amp_data/fly/", file_name[k], sep="")
  data <- read.csv(current_file)
  data <- na.omit(data)
  for(i in 1:nrow(data))
  {
    current_phase <- data[i,]$LAG
    diff_phase <- floor(current_phase / 4)
    day1 = data[i, 8:13]
    new_day1 <- day1
    for(j in 1:length(day1))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 6
      }
      new_day1[num] <- day1[j]
    }
    data[i, 8:13] <- new_day1
    
    day2 = data[i, 14:19]
    new_day2 <- day2
    for(j in 1:length(day2))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 6
      }
      new_day2[num] <- day2[j]
    }
    data[i, 14:19] <- new_day2
    
  }
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/fly/",file_name[k],sep="")
  write.csv(data, file_name_new)
  # data <- data[,6:29]
  # rt2<-apply(data,1,scale)
  # r<-t(rt2)
  # pheatmap(r,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  
  



#meta2d

rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/fly/JTK.flyyoung.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyyoung.csv")
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
rt <- cbind(rt,meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)])
data <- rt[,8:19]
rt <- rt[-which(rowMeans(data)==0),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.01 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.05 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.1 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.2 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#jtk
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/fly/JTK.flyyoung.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyyoung.csv")
rt <- rt[order(meta2d_data$JTK_pvalue),]
rt <- cbind(rt,meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)])
data <- rt[,8:19]
rt <- rt[-which(rowMeans(data)==0),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.01 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.05 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.1 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.2 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))

#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/fly/JTK.flyyoung.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/rain/young.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,8:19]
rt <- rt[-which(rowMeans(data)==0),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))







rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/fly/JTK.flyold.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyold.csv")
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
rt <- cbind(rt,meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)])
data <- rt[,8:19]
rt <- rt[-which(rowMeans(data)==0),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.01 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.05 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.1 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` < 0.2 & rt$`meta2d_data$meta2d_BH.Q[order(meta2d_data$meta2d_BH.Q)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))

#jtk
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/fly/JTK.flyold.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyold.csv")
rt <- rt[order(meta2d_data$JTK_pvalue),]
rt <- cbind(rt,meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)])
data <- rt[,8:19]
rt <- rt[-which(rowMeans(data)==0),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.01 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.05 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.1 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` < 0.2 & rt$`meta2d_data$JTK_pvalue[order(meta2d_data$JTK_pvalue)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/fly/JTK.flyold.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/fly/rain/old.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,8:19]
rt <- rt[-which(rowMeans(data)==0),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#arabidopsis
file_name <- list.files(path = "E:/GO_layers/final_photo/exp_amp_data/arabidopsis/", pattern = NULL, all.files = FALSE,full.names = FALSE, recursive = T,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

current_file <- paste("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/", file_name[1], sep="")
data <- read.csv(current_file)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$LAG
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 7:24]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase * 3
    if(num <= 0)
    {
      num <- num + 18
    }
    new_day1[num] <- day1[j]
  }
  data[i, 7:24] <- new_day1
  

  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/arabidopsis/",file_name[1],sep="")
write.csv(data, file_name_new)
  
  


current_file <- paste("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/", file_name[2], sep="")
data <- read.csv(current_file)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$LAG
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 7:24]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase * 3
    if(num <= 0)
    {
      num <- num + 18
    }
    new_day1[num] <- day1[j]
  }
  data[i, 7:24] <- new_day1
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/arabidopsis/",file_name[2],sep="")
write.csv(data, file_name_new)




###############arabidopsis meta2d
file_name2 <- paste("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_","flower.csv",sep="")
data <- read.csv(file_name2)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$meta2d_phase
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 19:36]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase*3
    if(num <= 0)
    {
      num <- num + 18
    }
    new_day1[num] <- day1[j]
  }
  data[i, 19:36] <- new_day1
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/arabidopsis/","flower.csv",sep="")
write.csv(data, file_name_new)



rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/JTK.flyyoung.csv")
rt <- rt[order(rt$ADJ.P),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$ADJ.P < 0.01 & rt$ADJ.P >0)
p001 <- p001[length(p001)]
p005 <- which(rt$ADJ.P < 0.05 & rt$ADJ.P >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$ADJ.P < 0.1 & rt$ADJ.P >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$ADJ.P < 0.2 & rt$ADJ.P >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#meta2d


yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/flower.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
data <- yeast_data[,20:37]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#jtk
yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/flower.csv")
rt <- yeast_data

meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$JTK_pvalue),]
data <- yeast_data[,20:37]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/flower.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/rain/flower.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,20:37]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))












file_name2 <- paste("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_","JTK.Agse5612.csv",sep="")
data <- read.csv(file_name2)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$meta2d_phase
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 24:29]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase
    if(num <= 0)
    {
      num <- num + 6
    }
    new_day1[num] <- day1[j]
  }
  data[i, 24:29] <- new_day1
  
  day2 = data[i, 30:35]
  new_day2 <- day2
  for(j in 1:length(day2))
  {
    num = j - diff_phase
    if(num <= 0)
    {
      num <- num + 6
    }
    new_day2[num] <- day2[j]
  }
  data[i, 30:35] <- new_day2
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/arabidopsis/","seed.csv",sep="")
write.csv(data, file_name_new)




#meta2d


yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/seed.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
data <- yeast_data[,25:36]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#jtk
yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/seed.csv")
rt <- yeast_data

meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$JTK_pvalue),]
data <- yeast_data[,25:36]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/seed.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/arabidopsis/rain/seed.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,25:36]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("blue", "blue","darkblue", "red","red"))(60))











rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/arabidopsis/JTK.flyyoung.csv")
rt <- rt[order(rt$ADJ.P),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$ADJ.P < 0.01 & rt$ADJ.P >0)
p001 <- p001[length(p001)]
p005 <- which(rt$ADJ.P < 0.05 & rt$ADJ.P >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$ADJ.P < 0.1 & rt$ADJ.P >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$ADJ.P < 0.2 & rt$ADJ.P >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#human
file_name <- list.files(path = "E:/GO_layers/final_photo/exp_amp_data/human/", pattern = NULL, all.files = FALSE,full.names = FALSE, recursive = T,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

current_file <- paste("E:/GO_layers/final_photo/exp_amp_data/human/", "cell.csv", sep="")
data <- read.csv(current_file)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$LAG
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 1)
  day1 = data[i, 7:30]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase
    if(num <= 0)
    {
      num <- num + 24
    }
    new_day1[num] <- day1[j]
  }
  data[i, 7:30] <- new_day1
  
  day2 = data[i, 31:54]
  new_day2 <- day2
  for(j in 1:length(day2))
  {
    num = j - diff_phase
    if(num <= 0)
    {
      num <- num + 24
    }
    new_day2[num] <- day2[j]
  }
  data[i, 31:54] <- new_day2
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/human/","cell.csv",sep="")
write.csv(data, file_name_new)




current_file <- paste("E:/GO_layers/final_photo/exp_amp_data/human/", "human_blood_jtk", sep="")
data <- read.csv(current_file)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$LAG
  diff_phase <- floor(current_phase / 3)
  day1 = data[i, 7:24]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase * 3
    if(num <= 0)
    {
      num <- num + 18
    }
    new_day1[num] <- day1[j]
  }
  data[i, 7:24] <- new_day1
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/human/",file_name[2],sep="")
write.csv(data, file_name_new)





rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/human/JTK.flyyoung.csv")
rt <- rt[order(rt$ADJ.P),]
data <- rt[,8:19]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$ADJ.P < 0.01 & rt$ADJ.P >0)
p001 <- p001[length(p001)]
p005 <- which(rt$ADJ.P < 0.05 & rt$ADJ.P >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$ADJ.P < 0.1 & rt$ADJ.P >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$ADJ.P < 0.2 & rt$ADJ.P >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))







#######################

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






############neurospora
file_name2 <- paste("E:/circidian_algorithm/result/neurospora/meta2d/meta2d_","JTKresult_RNA.csv",sep="")
data <- read.csv(file_name2)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$meta2d_phase
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 2)
  day1 = data[i, 24:35]
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
  data[i, 24:35] <- new_day1
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/neurospora/","neurospora.csv",sep="")
write.csv(data, file_name_new)

#meta2d


yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/neurospora/neurospora.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
data <- yeast_data[,25:36]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#jtk
yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/neurospora/neurospora.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$JTK_pvalue),]
data <- yeast_data[,25:36]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/neurospora/neurospora.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/neurospora/rain/neurospora.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,25:36]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))





#################cyanobacteria
file_name2 <- paste("E:/circidian_algorithm/result/cyanobacteria/meta2d/meta2d_","JTK.LLREPLICATE2.csv",sep="")
data <- read.csv(file_name2)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$meta2d_phase
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 4)
  day1 = data[i, 19:30]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase*2
    if(num <= 0)
    {
      num <- num + 12
    }
    new_day1[num] <- day1[j]
  }
  data[i, 19:30] <- new_day1
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/cyanobacteria/","cyanobacteria.csv",sep="")
write.csv(data, file_name_new)


#meta2d


yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/cyanobacteria/cyanobacteria.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
data <- yeast_data[,20:31]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#jtk
yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/cyanobacteria/cyanobacteria.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$JTK_pvalue),]
data <- yeast_data[,20:31]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/cyanobacteria/cyanobacteria.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/cyanobacteria/rain/cyanobacteria.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,20:31]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




###############human
file_name2 <- paste("E:/circidian_algorithm/result/human/meta2d/meta2d_","data.csv",sep="")
data <- read.csv(file_name2)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$meta2d_phase
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 1)
  day1 = data[i, 25:48]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase
    if(num <= 0)
    {
      num <- num + 24
    }
    new_day1[num] <- day1[j]
  }
  data[i, 25:48] <- new_day1
  
  
  day2 = data[i, 49:72]
  new_day2 <- day2
  for(j in 1:length(day2))
  {
    num = j - diff_phase 
    if(num <= 0)
    {
      num <- num + 24
    }
    new_day2[num] <- day2[j]
  }
  data[i, 49:72] <- new_day2
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/human/","cell.csv",sep="")
write.csv(data, file_name_new)


#meta2d


yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/human/cell.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
data <- yeast_data[,26:73]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#jtk
yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/human/cell.csv")
rt <- yeast_data

meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$JTK_pvalue),]
data <- yeast_data[,26:73]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/human/cell.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/human/rain/cell.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,26:73]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))







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
pheatmap(r,cluster_rows = FALSE,cluster_cols = pvalue,show_colnames = TRUE,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))



############################archaea 4 data
tissuename <- c("hal_1","hal_2","hal_4")
for(k in 1:length(tissuename))
{
  file_name2 <- paste("E:/circidian_algorithm/result/archaea/meta2d/meta2d_",tissuename[k],".csv",sep="")
  data <- read.csv(file_name2)
  data <- na.omit(data)
  for(i in 1:nrow(data))
  {
    current_phase <- data[i,]$meta2d_phase
    current_phase <- current_phase %% 24
    diff_phase <- floor(current_phase / 3)
    day1 = data[i, 24:31]
    new_day1 <- day1
    for(j in 1:length(day1))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 8
      }
      new_day1[num] <- day1[j]
    }
    data[i, 24:31] <- new_day1
    
    day2 = data[i, 32:39]
    new_day2 <- day2
    for(j in 1:length(day2))
    {
      num = j - diff_phase
      if(num <= 0)
      {
        num <- num + 8
      }
      new_day2[num] <- day2[j]
    }
    data[i, 32:39] <- new_day2
    
  }
  file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/archaea/",tissuename[k],".csv",sep="")
  write.csv(data, file_name_new)
  # data <- data[,6:29]
  # rt2<-apply(data,1,scale)
  # r<-t(rt2)
  # pheatmap(r,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))
  
  
}

file_name2 <- paste("E:/circidian_algorithm/result/archaea/meta2d/meta2d_","hal_3.csv",sep="")
data <- read.csv(file_name2)
data <- na.omit(data)
for(i in 1:nrow(data))
{
  current_phase <- data[i,]$meta2d_phase
  current_phase <- current_phase %% 24
  diff_phase <- floor(current_phase / 3)
  day1 = data[i, 24:31]
  new_day1 <- day1
  for(j in 1:length(day1))
  {
    num = j - diff_phase
    if(num <= 0)
    {
      num <- num + 8
    }
    new_day1[num] <- day1[j]
  }
  data[i, 24:31] <- new_day1
  
  
  
}
file_name_new <- paste("E:/GO_layers/final_photo/adjust_phase/archaea/","hal_3.csv",sep="")
write.csv(data, file_name_new)



#meta2d


yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/archaea/hal_4.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$meta2d_BH.Q),]
data <- yeast_data[,25:40]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$meta2d_BH.Q < 0.01 & rt$meta2d_BH.Q >0)
p001 <- p001[length(p001)]
p005 <- which(rt$meta2d_BH.Q < 0.05 & rt$meta2d_BH.Q >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$meta2d_BH.Q < 0.1 & rt$meta2d_BH.Q >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$meta2d_BH.Q < 0.2 & rt$meta2d_BH.Q >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))




#jtk
yeast_data <- read.csv("E:/GO_layers/final_photo/adjust_phase/archaea/hal_4.csv")
rt <- yeast_data
meta2d_data <- yeast_data
rt <- rt[order(meta2d_data$JTK_pvalue),]
data <- yeast_data[,25:40]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$JTK_pvalue < 0.01 & rt$JTK_pvalue >0)
p001 <- p001[length(p001)]
p005 <- which(rt$JTK_pvalue < 0.05 & rt$JTK_pvalue >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$JTK_pvalue < 0.1 & rt$JTK_pvalue >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$JTK_pvalue < 0.2 & rt$JTK_pvalue >0.1)
p02 <- p02[length(p02)]
pheatmap(r,labels_row=rt$ADJ.P, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = T,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))


#rain
rt<-read.csv("E:/GO_layers/final_photo/adjust_phase/archaea/hal_4.csv")
meta2d_data <- read.csv("E:/circidian_algorithm/result/archaea/rain/hal4.csv")
rt <- rt[order(meta2d_data$bhq),]
rt <- cbind(rt,meta2d_data$bhq[order(meta2d_data$bhq)])
data <- rt[,25:40]
rt2<-apply(data,1,scale) ##按行scale
r<-t(rt2)
p001 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.01 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0)
p001 <- p001[length(p001)]
p005 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.05 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.01)
p005 <- p005[length(p005)]
p01 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.1 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.05)
p01 <- p01[length(p01)]
p02 <- which(rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` < 0.2 & rt$`meta2d_data$bhq[order(meta2d_data$bhq)]` >0.1)
p02 <- p02[length(p02)]
pheatmap(r, gaps_row = c(p001,p005,p01,p02),cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,color = colorRampPalette(c("green", "darkgreen","black", "red","red"))(60))





#####################top exp variation range















#mosue
tissuename <- c("Adr","Aorta","BFAT","BS","Cere","Heart","Hypo","Kidney","Liver","Lung","Mus","WFAT")
cir_gene <- read.csv("E:/GO_layers/final_photo/circadian_gene_tissue/mouse_gene.csv")
new_data <- data.frame(0,0,0,0,0)
for(i in 1:12){
  expr = read.csv(paste0("E:/GO_layers/circadian_data/mouse/allgenes12tissues/",tissuename[i],'_Jtk.csv'))
  cir = expr[,'BH.Q']<0.05
  current_cir_gene <- unique(expr[cir,]$X.1)
  naidx <- which(expr$X.1=="#N/A")
  expr <- expr[-naidx,]
  rownames(expr) <- c(1:length(expr[,1]))
  expr_amp <- c()
  expr_cirmean <- c()
  expr_sd <- c()
  expr_cirrange <- c()
  for(j in 1:nrow(expr))
  {
    expr_cirmean[j] <- mean(as.numeric(expr[j,8:31]))
    expr_sd[j] <- sd(as.numeric(expr[j,8:31]))
    expr_cirrange[j] <- max(as.numeric(expr[j,8:31])) - min(as.numeric(expr[j,8:31]))
    expr_amp[j] <- expr[j,"AMP"]
  }
  expr_cirramp <- expr_amp / expr_cirmean
  expr_cirvaration <-  expr_sd /expr_cirmean
  expr_cirrange <- expr_cirrange / expr_cirmean
  expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$X.1
  expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$X.1
  expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$X.1
  expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$X.1
  mean_expr <- expr_mean[1:length(current_cir_gene)]
  vara_expr <- expr_varation[1:length(current_cir_gene)]
  range_expr <- expr_range[1:length(current_cir_gene)]
  ramp_expr <- expr_ramp[1:length(current_cir_gene)]
  new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
  print(length(intersect(current_cir_gene, mean_expr)))
  
}

write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/mouse.csv")

#human
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/human/blood.csv")
expr <- expr[-which(expr$X==0),]
cir = expr[,'meta3d_Pvalue']<0.05
current_cir_gene <- unique(expr[cir,]$X)
rownames(expr) <- c(1:length(expr[,1]))
expr_amp <- c()
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_cirramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,10:19]))
  expr_sd[j] <- sd(as.numeric(expr[j,10:19]))
  expr_cirrange[j] <- max(as.numeric(expr[j,10:19])) - min(as.numeric(expr[j,10:19]))
  expr_cirramp[j] <- expr[j,"meta3d_rAMP"]
}
expr_cirvaration <-  expr_sd /expr_cirmean
expr_cirrange <- expr_cirrange / expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$X
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$X
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$X
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$X
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))







expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/human/cell.csv")
expr <- expr[-which(expr$Gene.symbol==""),]
cir = expr[,'ADJ.P']<0.05
current_cir_gene= expr[cir,]$Gene.symbol
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,7:54]))
  expr_sd[j] <- sd(as.numeric(expr[j,7:54]))
  expr_cirrange[j] <- max(as.numeric(expr[j,7:54])) - min(as.numeric(expr[j,7:54]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_cirrange <- expr_cirrange / expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$Gene.symbol
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$Gene.symbol
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$Gene.symbol
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$Gene.symbol
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))
write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/human.csv")



#fly
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyyoung.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$Gene.Name)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,8:19]))
  expr_sd[j] <- sd(as.numeric(expr[j,8:19]))
  expr_cirrange[j] <- max(as.numeric(expr[j,8:19])) - min(as.numeric(expr[j,8:19]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_cirrange <- expr_cirrange/expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$Gene.Name
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$Gene.Name
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$Gene.Name
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$Gene.Name
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))


expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/fly/JTK.flyold.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$Gene.Name)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,8:19]))
  expr_sd[j] <- sd(as.numeric(expr[j,8:19]))
  expr_cirrange[j] <- max(as.numeric(expr[j,8:19])) - min(as.numeric(expr[j,8:19]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_cirrange <- expr_cirrange/ expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$Gene.Name
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$Gene.Name
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$Gene.Name
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$Gene.Name
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))

write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/fly.csv")


#yeast
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample2.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$Transcript_ID)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,7:26]))
  expr_sd[j] <- sd(as.numeric(expr[j,7:26]))
  expr_cirrange[j] <- max(as.numeric(expr[j,7:26])) - min(as.numeric(expr[j,7:26]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirrange <- expr_cirrange / expr_cirmean
expr_cirramp <- expr_ciramp /expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$Transcript_ID
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$Transcript_ID
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$Transcript_ID
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$Transcript_ID
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))


expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/yeast/JTK.Sample6.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$Transcript_ID)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,7:30]))
  expr_sd[j] <- sd(as.numeric(expr[j,7:30]))
  expr_cirrange[j] <- max(as.numeric(expr[j,7:30])) - min(as.numeric(expr[j,7:30]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_cirrange <- expr_cirrange / expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$Transcript_ID
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$Transcript_ID
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$Transcript_ID
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$Transcript_ID
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))
write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/yeast.csv")



#arabidopsis
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/flower.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$CycID)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,7:24]))
  expr_sd[j] <- sd(as.numeric(expr[j,7:24]))
  expr_cirrange[j] <- max(as.numeric(expr[j,7:24])) - min(as.numeric(expr[j,7:24]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_cirrange <- expr_cirrange / expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$CycID
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$CycID
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$CycID
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$CycID
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))



expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/arabidopsis/JTK.Agse5612.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$Gene.symbol)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,8:20]))
  expr_sd[j] <- sd(as.numeric(expr[j,8:20]))
  expr_cirrange[j] <- max(as.numeric(expr[j,8:20])) - min(as.numeric(expr[j,8:20]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirrange <- expr_cirrange / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$Gene.symbol
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$Gene.symbol
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$Gene.symbol
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$Gene.symbol
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))
write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/arabidopsis.csv")



#neurospora
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/neurospora/JTKresult_RNA.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$CycID)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,7:18]))
  expr_sd[j] <- sd(as.numeric(expr[j,7:18]))
  expr_cirrange[j] <- max(as.numeric(expr[j,7:18])) - min(as.numeric(expr[j,7:18]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirrange <- expr_cirrange / expr_cirmean
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$CycID
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$CycID
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$CycID
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$CycID
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))

write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/neurospora.csv")




#cyanobacteria
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/cyanobacteria/JTK.LLREPLICATE2.csv")
cir = expr[,'ADJ.P']<0.05
current_cir_gene <- unique(expr[cir,]$X7942_ID)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,9:28]))
  expr_sd[j] <- sd(as.numeric(expr[j,9:28]))
  expr_cirrange[j] <- max(as.numeric(expr[j,9:28])) - min(as.numeric(expr[j,9:28]))
  expr_ciramp[j] <- as.numeric(expr[j,"AMP"])
}
expr_cirrange <- expr_cirrange / expr_cirmean
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$X7942_ID
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$X7942_ID
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$X7942_ID
expr_ramp<- expr[order(expr_cirramp,decreasing = TRUE),]$X7942_ID
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))

write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/cyanobacteria.csv")




#archaea
new_data <- data.frame(0,0,0,0,0)
expr = read.csv("E:/GO_layers/final_photo/exp_amp_data/archaea/meta2d_AR2_achaea.csv")
cir = expr[,'JTK_pvalue']<0.05
current_cir_gene <- unique(expr[cir,]$CycID)
rownames(expr) <- c(1:length(expr[,1]))
expr_cirmean <- c()
expr_sd <- c()
expr_cirrange <- c()
expr_ciramp <- c()
for(j in 1:nrow(expr))
{
  expr_cirmean[j] <- mean(as.numeric(expr[j,20:55]))
  expr_sd[j] <- sd(as.numeric(expr[j,20:55]))
  expr_cirrange[j] <- max(as.numeric(expr[j,20:55])) - min(as.numeric(expr[j,20:55]))
  expr_ciramp[j] <- as.numeric(expr[j,"JTK_amplitude"])
}
expr_cirrange <- expr_cirrange / expr_cirmean
expr_cirramp <- expr_ciramp / expr_cirmean
expr_cirvaration <-  expr_sd /expr_cirmean
expr_mean <- expr[order(expr_cirmean,decreasing = TRUE),]$CycID
expr_varation <- expr[order(expr_cirvaration,decreasing = TRUE),]$CycID
expr_range <- expr[order(expr_cirrange,decreasing = TRUE),]$CycID
expr_ramp <- expr[order(expr_cirramp,decreasing = TRUE),]$CycID
mean_expr <- expr_mean[1:length(current_cir_gene)]
vara_expr <- expr_varation[1:length(current_cir_gene)]
range_expr <- expr_range[1:length(current_cir_gene)]
ramp_expr <- expr_ramp[1:length(current_cir_gene)]
new_data <- rbind(new_data, c(length(current_cir_gene),length(intersect(current_cir_gene, mean_expr)),length(intersect(current_cir_gene, vara_expr)),length(intersect(current_cir_gene, range_expr)),length(intersect(current_cir_gene, ramp_expr))))
print(length(intersect(current_cir_gene, mean_expr)))

write.csv(new_data, "E:/GO_layers/final_photo/mean_sd_range/archaea.csv")



##############################draw mean sd range photo
# 
# temp <- head(mpg)
# ggplot(temp,aes(x=class)) + geom_bar(aes(fill=factor(cyl))) # 按照cyl分组# 其他参数
# ggplot(mpg,aes(x=class)) + geom_bar(aes(col=factor(cyl)))
# ggplot(mpg,aes(x=class)) + geom_bar(aes(alpha=factor(cyl)))ggplot(mpg,aes(x=class)) + geom_bar(fill="blue") # 显示蓝色
# ggplot(mpg,aes(x=class)) + geom_bar(aes(fill="blue")) # 显示粉红色
# ggplot(mpg,aes(x=class)) + geom_bar(aes(fill="a")) # 显示粉红色
# 
# 
# 
# 
# 
# barplot(as.matrix(new_data), col = c('blue', 'orange', 'green', 'yellow', 'red', 'hotpink', 'cyan','purple', 'burlywood1', 'skyblue', 'gray'),
#         legend = rownames(data), 
#         cex.axis = 2, cex.names = 2,  las = 1, width = 0.5, space = 0.5, beside = FALSE,
#         args.legend = list(x = 'right', bty = 'n', inset = -0.18, cex = 2, y.intersp = 1.2, x.intersp = 0.7, text.width = 1))
# mtext('Relative Abundance(%)', cex = 2, side = 2, line = 4)
# dev.off()



# s1=c(1:3)
# s2=c(2:4)
# s3=c(3:5)
# #建立数据集
# d=data.frame(s1,s2,s3,row.names=c('miffy','kitty','tommy'))
# #转换数据集结构
# library(reshape2)
# df=melt(as.matrix(d))
# #使用ggplot2绘图
# library(ggplot2)
# ggplot(df,aes(Var2,value,fill=Var1))+geom_bar(stat="identity",position="fill")



data <- read.csv("E:/GO_layers/final_photo/mean_sd_range/data.csv",row.names=1)
filename <- colnames(data)
for(i in 1:ncol(data))
{
  temp_data <- data[,i]
  new_data <- data.frame(ncol=4,nrow=2)
  new_data[1,1] <- temp_data[2]
  new_data[2,1] <- temp_data[1] - temp_data[2]
  new_data[1,2] <- temp_data[3]
  new_data[2,2] <- temp_data[1] - temp_data[3]
  new_data[1,3] <- temp_data[4]
  new_data[2,3] <- temp_data[1] - temp_data[4]
  new_data[1,4] <- temp_data[5]
  new_data[2,4] <- temp_data[1] - temp_data[5]
  rownames(new_data) <- c("current","total-current")
  colnames(new_data) <- c("mean","variation","range","rAMP")
  
  new_data <- melt(as.matrix(new_data))
  a<- ggplot(new_data,aes(x=Var2,y=as.numeric(value),fill=Var1)) + geom_bar(stat="identity",position="fill",fill=c("red1","royalblue2","limegreen","royalblue2","mediumblue","royalblue2","gold","royalblue2"))+ theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+labs(x = filename[i],y = "percentage")

  temp_name <- paste("E:/GO_layers/final_photo/photo/mean_sd_range/plot",i,".pdf",sep="")
  ggsave(a, file=temp_name, width=5.22, height=5.22)
}


data <- read.csv("E:/GO_layers/final_photo/mean_sd_range/data.csv",row.names=1)
data <- t(data)
data <- as.matrix(data)
percent_mean <- data[,2] / data[,1]
percent_variation <- data[,3] / data[,1]
percent_range <- data[,4] / data[,1]
percent_ramp <- data[,5] / data[,1]
#dataset <- data.frame(value = t(c(percent_mean, percent_variation,percent_range,percent_ramp)), group = factor(c("Mean","Variation","Range","rAMP")))
#boxplot(c(percent_mean, percent_variation,percent_range),ylab="Percentage",border=c(colors()[38], 'steelblue2',colors()[37]))


dataset <- data.frame(value = c(percent_mean, percent_variation,percent_range,percent_ramp), group = factor(rep(c("Mean","Variation","Range","rAMP"),each=23)))
aov(dataset$value~dataset$group)
boxplot(percent_mean,percent_variation,percent_range,percent_ramp,ylab="Percentage",col=c("red1", 'limegreen',"mediumblue","gold"))



data <- read.csv("E:/GO_layers/proteinhuman.csv",header=TRUE)
melt(data)
for(i in 1:nrow(data))
{
  if(t.test(data[i,2:7],data[i,8:14])$p.value < 0.05)
  {
    print(i)
  }
}





####################################mean range ramp venn plot
#ramp
data <- read.csv("E:/GO_layers/final_photo/mean_sd_range/data.csv")
ramp_data <- data[5,]
venn_data_new <- read.csv("E:/GO_layers/final_photo/mean_sd_range/plotdata.csv")
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
  temp_plot <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, ext.text=FALSE, rotation.degree = number_flag, resolution =300, imagetype="png", col="black", fill=c("violet", 'steelblue2'), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
  plot_list <- c(plot_list, list(temp_plot))
}
grid.newpage()
plot_grid(plotlist = plot_list,nrow=6, ncol=4, axis="bl",rel_widths=c(1,1,1))









#range

data <- read.csv("E:/GO_layers/final_photo/mean_sd_range/data.csv")
ramp_data <- data[6,]
venn_data_new <- read.csv("E:/GO_layers/final_photo/mean_sd_range/plotdata.csv")
plot_list = list()
for(i in 1:nrow(venn_data_new))
{
  data1_all = as.numeric(venn_data_new[i,3])
  data2_all = as.numeric(venn_data_new[i,4])
  common_num = as.numeric(venn_data_new[i,6])
  if(data1_all > data2_all)
  {
    number_flag = 0
  }else
  {
    number_flag = 180
  }
  temp_plot <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, ext.text=FALSE, height = 450, width = 450, rotation.degree = number_flag, resolution =300, imagetype="png", col="black", fill=c("violet", 'steelblue2'), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
  plot_list <- c(plot_list, list(temp_plot))
}
grid.newpage()
plot_grid(plotlist = plot_list,nrow=6, ncol=4, axis="bl",rel_widths=c(1,1,1))




















########################kegg
library(clusterProfiler)

library(org.Mm.eg.db)
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

mus_kegg <- download_KEGG("mus")


for(i in 1:12)
{
  a <- enrichKEGG(gene=r1[,i],organism="mmu")
  fileanme <- paste("E:/GO_layers/kegg/noback/mouse",i,".csv",sep="")
  write.csv(a@result,fileanme)
}


for(i in 1:12)
{
  a <- enrichKEGG(gene=r1[,i],organism="mmu",universe=expr_gene)
  fileanme <- paste("E:/GO_layers/kegg/back/mouse",i,".csv",sep="")
  write.csv(a@result,fileanme)
}





ck1<-compareCluster(geneClusters = r1,fun = "enrichKEGG",organism="mmu",pvalueCutoff=1,qvalueCutoff=1)
ck1<-compareCluster(geneClusters = r1,fun = "enrichKEGG",organism="mmu",pvalueCutoff=1,qvalueCutoff=1,universe=expr_gene)


ck2<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",pvalueCutoff=1,qvalueCutoff=1,ont="MF",readable=TRUE,universe=expr_gene)
ck3<-compareCluster(geneClusters = r1,fun = "enrichGO",OrgDb="org.Mm.eg.db",pvalueCutoff=1,qvalueCutoff=1,ont="CC",readable=TRUE,universe=expr_gene)





















#################################################amp and exp new meta2d and rain
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

# Save the file
dev.off()








#data
tissuename <- c("adr","AOR","BFAT","bst","CER","HEA","HYP","KID","LIV","LUN","MUS","WAT")
par(mfcol=c(4,3))
aa= matrix(0,nrow = 12,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
bb=1:12




for(i in 1:12){
  expr = read.csv(paste0("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[i],'.csv'))
  file_name1 <- paste("E:/circidian_algorithm/circidian_gene/mouse/",tissuename[i],'.csv',sep="")
  gene_data <- read.csv(file_name1)
  cir = gene_data$x
  expr_cir = expr[expr$CycID %in% cir,]
  expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
  expr_cirmean= rowMeans(expr_cir[,24:35])
  
  
  
  
  
  plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),xlab = '',ylab = '',main=tissuename[i])
  aa[i,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
  aa[i,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
  aa[i,3]= aa[i,1]^2
  bb[i] = sum(cir)
}
#write.csv(aa,"correlation_between_amp_and_expr.csv")

{
  #aa= read.csv("correlation_between_amp_and_expr.csv",row.names = 1)
  aaa= matrix(0,nrow = 12,ncol = 2)
  rownames(aaa)= tissuename
  colnames(aaa)=c('expression','others')
  aaa[,1]= aa[,'R2']
  aaa[,2]= 1- aa[,'R2']
  
  par(mfcol=c(4,3))
  parmar= par()$mar
  par(mar = c(0,0, 2, 0),las = 1)
  for(i in 1:12){
    pie(aaa[i,],labels = round(aaa[i,]*100)/100,main = tissuename[i],col= c('skyblue','white'),init.angle = 90)
  }
  par(mar=parmar,las= 0)
  
  par(mfcol=c(1,1))
  plot.new()
  legend('top',c('expression','others'),fill = c('skyblue','white'),cex=3)
}
tissuenamefull= c("Adrenal Gland","Aorta","Brown Fat","Brainstem","Cerebellum","Heart","Hypothalamus",
                  "Kidney","Liver","Lung","Skeletal Muscle","White Fat")


{
  #aa= read.csv("correlation_between_amp_and_expr.csv",row.names = 1)
  aaa= matrix(0,nrow = 12,ncol = 2)
  rownames(aaa)= tissuename
  colnames(aaa)=c('expression','others')
  aaa[,1]= aa[,'R2']
  aaa[,2]= 1- aa[,'R2']
  # png("fig2.png",width=10,height=12,units="in",res=350)
  # par(mfrow=c(4,4))
  par(plt=c(0.15,0.85,0.25,0.65))
  par(mai=c(0.6732,0.5412,0.5412,0.2772))
  
  {
    mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
    widths <-c(9,9,9,9)  
    heights <-c(6,6,6,2)  
    layout(mat,widths = widths, heights = heights)  
    nf <-layout(mat,  widths = widths, heights = heights)  
    layout.show(nf)
  }
  
  for(i in 1:12){
    expr = read.csv(paste0("E:/circidian_algorithm/result/mouse/4h/JTK/meta2d_",tissuename[i],'.csv'))
    file_name1 <- paste("E:/circidian_algorithm/circidian_gene/mouse/",tissuename[i],'.csv',sep="")
    gene_data <- read.csv(file_name1)
    cir = gene_data$x
    expr_cir = expr[expr$CycID %in% cir,]
    expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
    expr_cirmean= rowMeans(expr_cir[,24:35])
    
    
    
    
    plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
         xlab = 'expression',ylab = 'amplitude',main=tissuenamefull[i],
         xlim = c(2,10),ylim = c(-2,9),
         pch=20,cex=0.5)
    text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)
    
    par(plt=c(0.60,0.95,0.20,0.55))
    par(new=T)
    pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
    text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
    par(plt=c(0.15,0.85,0.25,0.65))
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
  }
  {
    # parmai=par()$mai
    # par(mai=rep(0,4))
    # plot.new()
    # plot.new()
    # legend('top',c('variation explained by expression','variation explained by other factors'),
    #        fill = c('skyblue','white'),cex=1.5,bty='n')
    # par(mai=parmai)
  }
  
  # dev.off()
  
  {
    parraw_mar= par()$mar
    par(mar = c(0,0, 0, 0),las = 1)  
    plot(c(0, 4),c(0, 4),  
         xlim = c(0, 4), ylim = c(0, 4),  
         xaxs = 'i', yaxs = 'i',  
         xaxt = 'n', yaxt = 'n', type = 'n', ann =F, axes = F) 
    legend('top',c('variation explained by expression','variation explained by other factors'),
           fill = c('skyblue','white'),cex=1.5,bty='n')
    par(mar=parraw_mar,las= 0)
  }
}




#fly
tissuename<- c("young","old")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyyoung.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/fly/young.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:35])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))



expr = read.csv("E:/circidian_algorithm/result/fly/meta2d/meta2d_JTK.flyold.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/fly/old.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:35])
aa[2,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[2,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[2,3]= aa[2,1]^2

aaa[2,1]= aa[2,'R2']
aaa[2,2]= 1- aa[2,'R2']
i=2
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#yeast
tissuename<- c("Sample2","Sample6")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/yeast/meta2d/meta2d_JTK.Sample2.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/yeast/sample2.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:43])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

expr = read.csv("E:/circidian_algorithm/result/yeast/meta2d/meta2d_JTK.Sample6.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/yeast/sample6.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:47])
aa[2,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[2,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[2,3]= aa[2,1]^2

aaa[2,1]= aa[2,'R2']
aaa[2,2]= 1- aa[2,'R2']
i=2
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)

#arabidopsis
tissuename<- c("flower","seed")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_flower.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/arabidopsis/flower.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,19:36])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

expr = read.csv("E:/circidian_algorithm/result/arabidopsis/meta2d/meta2d_JTK.Agse5612.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/arabidopsis/seed.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:36])
aa[2,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[2,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[2,3]= aa[2,1]^2

aaa[2,1]= aa[2,'R2']
aaa[2,2]= 1- aa[2,'R2']
i=2
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))


#neurospora
tissuename<- c("neurospora")
aa= matrix(0,nrow = 1,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 1,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/neurospora/meta2d/meta2d_JTKresult_RNA.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/neurospora/neurospora.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:35])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))


#cyanobacteria
tissuename<- c("cyanobacteria")
aa= matrix(0,nrow = 1,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 1,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/cyanobacteria/meta2d/meta2d_JTK.LLREPLICATE2.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/cyanobacteria/cyanobacteria.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,19:38])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))
#barplot(t(aaa), main = "R2", names.arg = tissuename, xlab = "tissue", ylab = "percent", col = colors)


#archaea
tissuename<- c("hal1","hal2","hal3","hal4")
aa= matrix(0,nrow = 4,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 4,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_1.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/archaea/hal1.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:42])
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],

     pch=20,cex=0.5)
text(-0.5,0,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

expr = read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_2.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/archaea/hal2.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:42])
aa[2,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[2,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[2,3]= aa[2,1]^2

aaa[2,1]= aa[2,'R2']
aaa[2,2]= 1- aa[2,'R2']
i=2
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     
     pch=20,cex=0.5)
text(-0.2,0,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))


expr = read.csv("E:/circidian_algorithm/result/archaea/meta2d/meta2d_hal_4.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/archaea/hal4.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= rowMeans(expr_cir[,24:42])
aa[3,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[3,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[3,3]= aa[3,1]^2

aaa[3,1]= aa[3,'R2']
aaa[3,2]= 1- aa[3,'R2']
i=3
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[4],
     
     pch=20,cex=0.5)
text(-0.6,-0.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))



#human
tissuename<- c("cell","epistem")
aa= matrix(0,nrow = 2,ncol = 3)
rownames(aa)= tissuename
colnames(aa)= c('cor','pvalue','R2')
aaa= matrix(0,nrow = 2,ncol = 2)
expr = read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_data.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/human/cell.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= expr_cir$meta2d_Base
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[1,3]= aa[1,1]^2

aaa[1,1]= aa[1,'R2']
aaa[1,2]= 1- aa[1,'R2']
i=1
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

expr = read.csv("E:/circidian_algorithm/result/human/meta2d/meta2d_RAW DATA.csv")
gene_data <- read.csv("E:/circidian_algorithm/circidian_gene/human/epistem.csv")
cir = gene_data$x
expr_cir = expr[expr$CycID %in% cir,]
expr_cir <- expr_cir[expr_cir[,'meta2d_AMP']!=0,]
expr_cirmean= expr_cir$meta2d_Base
aa[2,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$estimate
aa[2,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'meta2d_AMP'])))$p.value
aa[2,3]= aa[1,1]^2

aaa[2,1]= aa[2,'R2']
aaa[2,2]= 1- aa[2,'R2']
i=2
plot(log(expr_cirmean),log(expr_cir[,'meta2d_AMP']),
     xlab = 'expression',ylab = 'amplitude',main=tissuename[i],
     xlim = c(2,10),ylim = c(-2,9),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[i,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)

par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[i,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[i,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

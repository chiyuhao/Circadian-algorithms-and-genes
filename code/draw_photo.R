library(UpSetR)
library(grid)
library(futile.logger)
library(VennDiagram)
library(gridGraphics)
library(grid)
library(ggplotify)
library(cowplot)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
movies <- read.csv(system.file("extdata","movies.csv",package = "UpSetR"), header = TRUE, sep=";")
upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
between <- function(row, min, max){
  newData <- (row["ReleaseDate"] < max) & (row["ReleaseDate"] > min)
}

upset(movies, sets=c("Drama","Comedy","Action","Thriller","Western","Documentary"),
      queries = list(list(query = intersects, params = list("Drama", "Thriller")),
                     list(query = between, params=list(1970,1980), color="red", active=TRUE)))
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
peakAnno



# A simple two-set diagram （根据数据量多少而确定圆的大小）
venn.plot <- draw.pairwise.venn(100, 70, 30, c("First", "Second"));
grid.draw(venn.plot);
grid.newpage();

# Same diagram as above, but without scaling（不会根据数据的多少而自动适应圆的大小）
venn.plot <- draw.pairwise.venn(100, 70, 30, c("First", "Second"), scaled = FALSE);
grid.draw(venn.plot);
grid.newpage();

venn.plot <- draw.pairwise.venn(
  area1 = 100,#区域1的数
  area2 = 70,#区域2的数
  cross.area = 68,#交叉数
  category = c("First", "Second"),#分类名称
  fill = c("blue", "red"),#区域填充颜色
  lty = "blank",#区域边框线类型
  cex = 2, #区域内部数字的字体大小
  cat.cex = 2, #分类名称字体大小
  cat.pos = c(285, 105), #分类名称在圆的位置，默认正上方，通过角度进行调整
  cat.dist = 0.09,#分类名称距离边的距离（可以为负数）
  cat.just = list(c(-1, -1), c(1, 1)),#分类名称的位置
  ext.pos = 30,#线的角度 默认是正上方12点位置
  ext.dist = -0.05,#外部线的距离
  ext.length = 0.85, #外部线长度
  ext.line.lwd = 2, #外部线的宽度
  ext.line.lty = "dashed"#外部线为虚线
);
grid.draw(venn.plot);


venn.plot <- draw.triple.venn(
  area1 = 65,
  area2 = 75,
  area3 = 85,
  n12 = 35,
  n23 = 15,
  n13 = 25,
  n123 = 5,
  category = c("First", "Second", "Third"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
);
grid.draw(venn.plot);#画图展示
# Writing to file
tiff(filename = "Triple_Venn_diagram.tiff", compression = "lzw");?? #保存图片
dev.off();


A <- sample(1:1000, 400, replace = FALSE);
B <- sample(1:1000, 600, replace = FALSE);
C <- sample(1:1000, 350, replace = FALSE);
D <- sample(1:1000, 550, replace = FALSE);
E <- sample(1:1000, 375, replace = FALSE);
venn.plot <- venn.diagram(
  x = list(
    A = A,
    B = B,
    C = C,
    D = D,
    E = E
  ),
  filename = NULL,
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
);
grid.draw(venn.plot);#画图展示












library(VennDiagram)

files=list.files(path = ".", pattern = "type")

for (i in files){
  
  a=read.table(i)
  
  individual=strsplit(i,"\\.")[[1]][1]
  
  image_name=paste(individual,".tiff",sep="")
  
  IGM=which(a[,3]>0)
  
  IGA=which(a[,4]>0)
  
  IGG=which(a[,5]>0)
  
  venn.diagram(list(IGM=IGM,IGA=IGA,IGG=IGG), fill=c("red","green","blue"), alpha=c(0.5,0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename=image_name)
  
}


par(mfrow = c(2, 2))
x <- seq(-pi,pi,by=0.1)
plot(x,sin(x),typ="l")
plot(x,cos(x))
plot(x,2*sin(x)*cos(x))
plot(x,tan(x))










###############################################################################start  term photo
human_data_all <- read.csv("E:/GO_layers/expressed_gene_enrich/human_result_padj.csv")
#human_term <- unique(human_data$ID)
mouse_data_all <- read.csv("E:/GO_layers/expressed_gene_enrich/mouse_result_padj.csv")
#mouse_term <- unique(mouse_data$ID)
fly_data <- read.csv("E:/GO_layers/expressed_gene_enrich/fly_result_padj.csv")
#fly_term <- unique(fly_data$ID)
yeast_data <- read.csv("E:/GO_layers/expressed_gene_enrich/yeast_result_padj.csv")
#yeast_term <- unique(yeast_data$ID)
arabidopsis_data <- read.csv("E:/GO_layers/expressed_gene_enrich/arabidopsis_result_padj.csv")
#arabidopsis_term <- unique(arabidopsis_data$ID)
cyanobacteria_data <- read.csv("E:/GO_layers/expressed_gene_enrich/cyanobacteria_result_padj.csv")
archaea_result_data <- read.csv("E:/GO_layers/expressed_gene_enrich/archaea_result_padj.csv")
neurospra_result_data <- read.csv("E:/GO_layers/expressed_gene_enrich/neurospra_result_padj.csv")
species_data <- list(human_data_all, mouse_data_all,fly_data,yeast_data,arabidopsis_data,cyanobacteria_data,archaea_result_data,neurospra_result_data)

plot_list_new <- list()
for(k in 1:7)
{
  for(l in (k+1):8)
  {
    human_data <- species_data[[k]][species_data[[k]]$ont=="BP",]
    mouse_data <- species_data[[l]][species_data[[l]]$ont=="BP",]
    temp_tissue1 <- unique(human_data$Cluster)
    temp_tissue2 <- unique(mouse_data$Cluster)
    common_term_num <- 0
    total_num <- 0
    data1_term_num <- 0
    data2_term_num <- 0
    fly_tissue <- temp_tissue1
    human_tissue <- fly_tissue
    mouse_tissue <- temp_tissue2
    # human_data <- species_data[[k]][species_data[[k]]$ont=="CC",]
    # mouse_data <- species_data[[k]][species_data[[l]]$ont=="CC",]
    for(i in 1:length(fly_tissue))
    {
      for(j in 1:length(mouse_tissue))
      {
        total_num <- total_num + 1
        common_term_num <- common_term_num + length(intersect(human_data[human_data$Cluster==human_tissue[i],]$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
        data1_term_num <- data1_term_num + length(setdiff(human_data[human_data$Cluster==human_tissue[i],]$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
        data2_term_num <- data2_term_num + length(setdiff(mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID, human_data[human_data$Cluster==human_tissue[i],]$ID))
      }
    }
    total_num
    common_num = common_term_num/total_num
    data1_all = data1_term_num/total_num + common_num
    data2_all = data2_term_num/total_num + common_num
    temp_b <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
    plot_list_new <- c(plot_list_new, list(temp_b))
  }
  
}
species_new <- c('human', 'mouse','fly','yeast','arabidopsis','cyanobacteria','archaea','neurospra')
species_new1 <- c('neurospra','archaea','cyanobacteria','arabidopsis','yeast','fly','mouse','human')

plot_data1 <- list()
plot_data1[c(2,3,4,5,6,7,8,11,12,13,14,15,16,20,21,22,23,24,29,30,31,32,38,39,40,47,48,56)] = plot_list_new 
vp.scatter <- viewport(x=0, y=0, width=0.66, height=0.66, just=c("left", "bottom"))
vp.scatter1 <- viewport(x=0.11, y=0.09, width=0.50, height=0.50, just=c("left", "bottom"))

print(as.ggplot(expression(plot(1:8,1:8,type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:8,labels = species_new) + axis(side = 4,at = 1:8,labels = species_new1))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=8, ncol=8, axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)








for(i in 1:length(species_data))
{
  print(species_new[i])
  print(length(unique(species_data[[i]][species_data[[i]]$ont=="BP",'ID'])))
}





#####################not average term photo
human_data_all <- read.csv("E:/GO_layers/expressed_gene_enrich/human_result_padj.csv")
#human_term <- unique(human_data$ID)
mouse_data_all <- read.csv("E:/GO_layers/expressed_gene_enrich/mouse_result_padj.csv")
#mouse_term <- unique(mouse_data$ID)
fly_data <- read.csv("E:/GO_layers/expressed_gene_enrich/fly_result_padj.csv")
#fly_term <- unique(fly_data$ID)
yeast_data <- read.csv("E:/GO_layers/expressed_gene_enrich/yeast_result_padj.csv")
#yeast_term <- unique(yeast_data$ID)
arabidopsis_data <- read.csv("E:/GO_layers/expressed_gene_enrich/arabidopsis_result_padj.csv")
#arabidopsis_term <- unique(arabidopsis_data$ID)
cyanobacteria_data <- read.csv("E:/GO_layers/expressed_gene_enrich/cyanobacteria_result_padj.csv")
archaea_result_data <- read.csv("E:/GO_layers/expressed_gene_enrich/archaea_result_padj.csv")
neurospra_result_data <- read.csv("E:/GO_layers/expressed_gene_enrich/neurospra_result_padj.csv")
species_data <- list(human_data_all, mouse_data_all,fly_data,yeast_data,arabidopsis_data,cyanobacteria_data,archaea_result_data,neurospra_result_data)
species_new <- c('human', 'mouse','fly','yeast','arabidopsis','cyanobacteria','archaea','neurospra')



plot_list_new <- list()
for(k in 1:7)
{
  for(l in (k+1):8)
  {
    human_data <- species_data[[k]][species_data[[k]]$ont=="BP",]
    mouse_data <- species_data[[l]][species_data[[l]]$ont=="BP",]
    temp_tissue1 <- unique(human_data$Cluster)
    temp_tissue2 <- unique(mouse_data$Cluster)
    common_term_num <- 0
    total_num <- 0
    data1_term_num <- 0
    data2_term_num <- 0
    fly_tissue <- temp_tissue1
    human_tissue <- fly_tissue
    mouse_tissue <- temp_tissue2
    # human_data <- species_data[[k]][species_data[[k]]$ont=="CC",]
    # mouse_data <- species_data[[k]][species_data[[l]]$ont=="CC",]

    total_num <- total_num + 1
    common_term_num <- common_term_num + length(intersect(human_data$ID, mouse_data$ID))
    data1_term_num <- data1_term_num + length(setdiff(human_data$ID, mouse_data$ID))
    data2_term_num <- data2_term_num + length(setdiff(mouse_data$ID, human_data$ID))

    total_num
    common_num = common_term_num/total_num
    data1_all = data1_term_num/total_num + common_num
    data2_all = data2_term_num/total_num + common_num
    temp_b <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
    plot_list_new <- c(plot_list_new, list(temp_b))
  }
  
}
species_new <- c('human', 'mouse','fly','yeast','arabidopsis','cyanobacteria','archaea','neurospra')
species_new1 <- c('neurospra','archaea','cyanobacteria','arabidopsis','yeast','fly','mouse','human')

plot_data1 <- list()
plot_data1[c(2,3,4,5,6,7,8,11,12,13,14,15,16,20,21,22,23,24,29,30,31,32,38,39,40,47,48,56)] = plot_list_new 
vp.scatter <- viewport(x=0, y=0, width=0.66, height=0.66, just=c("left", "bottom"))
vp.scatter1 <- viewport(x=0.11, y=0.09, width=0.50, height=0.50, just=c("left", "bottom"))

print(as.ggplot(expression(plot(1:8,1:8,type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:8,labels = species_new) + axis(side = 4,at = 1:8,labels = species_new1))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=8, ncol=8, axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)

######################################################################################

fly_tissue <- unique(fly_data$Cluster)
mouse_tissue <- unique(mouse_data$Cluster)
common_term_num <- 0
total_num <- 0
data1_term_num <- 0
data2_term_num <- 0
human_data <- fly_data
human_tissue <- fly_tissue
for(i in 1:length(fly_tissue))
{
  for(j in 1:length(mouse_tissue))
  {
    total_num <- total_num + 1
    common_term_num <- common_term_num + length(intersect(human_data[human_data$Cluster==human_tissue[i],]$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
    data1_term_num <- data1_term_num + length(setdiff(human_data[human_data$Cluster==human_tissue[i],]$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
    data2_term_num <- data2_term_num + length(setdiff(mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID, human_data[human_data$Cluster==human_tissue[i],]$ID))
  }
}
total_num
common_num = common_term_num/total_num
data1_all = data1_term_num/total_num + common_num
data2_all = data2_term_num/total_num + common_num
common_num
data1_all
data2_all



b2 <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
grid.draw(b2)




yeast_tissue <- unique(yeast_data$Cluster)
mouse_tissue <- unique(mouse_data$Cluster)
common_term_num <- 0
total_num <- 0
data1_term_num <- 0
data2_term_num <- 0
human_data <- yeaast_data
for(i in 1:1)
{
  for(j in 1:length(mouse_tissue))
  {
    total_num <- total_num + 1
    common_term_num <- common_term_num + length(intersect(human_data$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
    data1_term_num <- data1_term_num + length(setdiff(human_data$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
    data2_term_num <- data2_term_num + length(setdiff(mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID, human_data$ID))
  }
}
total_num
common_num = common_term_num/total_num
data1_all = data1_term_num/total_num + common_num
data2_all = data2_term_num/total_num + common_num

b3 <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
grid.draw(b3)






human_tissue <- unique(human_data$Cluster)
mouse_tissue <- unique(mouse_data$Cluster)
common_term_num <- 0
total_num <- 0
data1_term_num <- 0
data2_term_num <- 0
for(i in 1:length(human_tissue))
{
  for(j in 1:length(mouse_tissue))
  {
    total_num <- total_num + 1
    common_term_num <- common_term_num + length(intersect(human_data[human_data$Cluster==human_tissue[i],]$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
    data1_term_num <- data1_term_num + length(setdiff(human_data[human_data$Cluster==human_tissue[i],]$ID, mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID))
    data2_term_num <- data2_term_num + length(setdiff(mouse_data[mouse_data$Cluster==mouse_tissue[j],]$ID, human_data[human_data$Cluster==human_tissue[i],]$ID))
  }
}
total_num
common_num = common_term_num/total_num
data1_all = data1_term_num/total_num + common_num
data2_all = data2_term_num/total_num + common_num

b1 <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
grid.draw(b1)



plot_grid(plotlist = list(b1,b2,b3),labels = c("human","fly","yeast"))




















species <- c('mouse','human','fly','yeast','arabidopsis')
term_num <- c(length(mouse_term), length(human_term), length(fly_term), length(yeast_term), length(arabidopsis_term))
term_num_data <- data.frame(species, term_num)

photo_input <- list(mouse_term, human_term, fly_term, yeast_term, arabidopsis_term)
intersect_name <- c()
intersect_num <- c()
for(i in 1:5)
{
  if (i != 5)
  {
    for(j in (i + 1) : 5)
    {
      intersect_name[length(intersect_name) + 1] <- paste(species[i], species[j], sep="_")
      intersect_num[length(intersect_num) + 1] <- length(intersect(photo_input[[i]], photo_input[[j]]))
    }
  }
} 

right_photo_input <- data.frame(intersect_name, intersect_num)



b1 <- venn.diagram(x=list(x=mouse_term, y=human_term),filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
#grid.draw(b1)
b2 <- venn.diagram(x=list(x=mouse_term, y=fly_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),  cat.cex=0)
#grid.draw(b2)
b3 <- venn.diagram(x=list(x=mouse_term, y=yeast_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)
#grid.draw(b3)
b4 <- venn.diagram(x=list(x=mouse_term,y=arabidopsis_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)

b5 <- venn.diagram(x=list(x=human_term, y=fly_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)

b6 <- venn.diagram(x=list(x=human_term, y=yeast_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)
b7 <- venn.diagram(x=list(x=human_term, y=arabidopsis_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[616], colors()[38]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)
b8 <- venn.diagram(x=list(x=fly_term, y=yeast_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)
b9 <- venn.diagram(x=list(x=fly_term, y=arabidopsis_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[616], colors()[38]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)
b10 <- venn.diagram(x=list(x=yeast_term,y=arabidopsis_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[616], colors()[38]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300), cat.cex=0)
plot_data <- list(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
#plot_grid(plotlist = plot_data,labels = LETTERS[1:10])
plot_data1 <- list()
plot_data1[c(1,2,3,4,6,7,8,11,12,16)] = plot_data 

# plot(1:3)
# grid(NA, 5, lwd = 2)

#p.scatter <- ggplot(data.frame(0,0,0,0,0)) +geom_histogram(aes())

#p.scatter <- ggplot() + geom_point(colour="0", size = 4) + scale_x_discrete('A')

point_data <- data.frame(species, species)




#p.hist.len <- ggplot(term_num_data) + geom_histogram(aes(x=species))
p.hist.len <- ggplot(term_num_data, mapping=aes(x=species,y=term_num,fill=species,group=factor(1))) + geom_bar(stat="identity")
#p.hist.wid <- ggplot(iris) + geom_histogram(aes(x=Sepal.Width)) + coord_flip()
p.hist.wid <- ggplot(right_photo_input, mapping=aes(x=intersect_name,y=intersect_num,fill=intersect_name, group=factor(1))) + geom_bar(stat="identity") + coord_flip() + theme(legend.position = "none")
p.scatter<-ggplot(data=point_data, aes(x=point_data[,1], y=point_data[,2]))+geom_point(color='1',size=0) + scale_x_discrete(species) + theme_bw()  + theme(panel.border = element_blank() ,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))



grid.newpage()
vp.len <- viewport(x=0, y=0.66, width=0.66, height=0.34, just=c("left", "bottom"))
vp.wid <- viewport(x=0.66, y=0, width=0.34, height=0.66, just=c("left", "bottom"))
vp.scatter <- viewport(x=0, y=0, width=0.66, height=0.66, just=c("left", "bottom"))
vp.scatter1 <- viewport(x=0.18, y=0.13, width=0.40, height=0.40, just=c("left", "bottom"))

# direct the charts into the specified viewport
print(p.hist.len, vp=vp.len)
print(p.hist.wid, vp=vp.wid)
#print(p.scatter, vp=vp.scatter)

species_test <- c('mouse','human','fly','yeast','arabid','')
species_test1 <- c('arabid','yeast','fly','human','mouse','')
print(as.ggplot(expression(plot(1:6,1:6,type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:6,labels = species_test) + axis(side = 4,at = 1:6,labels = species_test1))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=4, ncol=4, axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)

#print(plot_grid(plotlist = plot_data,labels = LETTERS[1:6]), vp=vp.scatter)

#grid.newpage()












upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
between <- function(row, min, max){
  newData <- (row["ReleaseDate"] < max) & (row["ReleaseDate"] > min)
}






intersect_name <- c()
intersect_num <- c()
for(i in 1:5)
{
  intersect_name[length(intersect_name) + 1] <- species[i]
  intersect_num[length(intersect_num) + 1] <- length(unique(photo_input[[i]]))
  if (i != 5)
  {
    for(j in (i + 1) : 5)
    {
      intersect_name[length(intersect_name) + 1] <- paste(species[i], species[j], sep="&")
      intersect_num[length(intersect_num) + 1] <- length(intersect(photo_input[[i]], photo_input[[j]]))
    }
  }
}

input <- intersect_num
names(input) <- intersect_name
data <- fromExpression(input)
upset(data)
###############################################################################################################end


label.col = c("darkred", "white", "darkblue", "white","white", "white", "darkgreen")


b2 <- venn.diagram(x=list(human_term=human_term, mouse_term=mouse_term), filename=NULL,  col="white", fill=c(colors()[616], colors()[38]),cat.fontfamily = "serif",fontface = "bold",lwd=c(1, 1),scale=FALSE)
grid.draw(b2)



#############################################
human_data <- read.csv("E:/GO_layers/photo/result/human/human_result_padj.csv")
human_term <- human_data$ID
mouse_data <- read.csv("E:/GO_layers/photo/result/mouse/mouse_result_padj.csv")
mouse_term <- mouse_data$ID
fly_data <- read.csv("E:/GO_layers/photo/result/fly/fly_result_padj.csv")
fly_term <- fly_data$ID
yeast_data <- read.csv("E:/GO_layers/photo/result/yeast/yeast_result_padj.csv")
yeast_term <- yeast_data$ID
arabidopsis_data <- read.csv("E:/GO_layers/photo/result/Arabidopsis/arabidopsis_result_padj.csv")
arabidopsis_term <- arabidopsis_data$ID
photo_input <- list(human_term, fly_term, mouse_term,yeast_term,arabidopsis_term)
par(mfrow = c(3, 4))

draw.pairwise.venn(area1=7,area2=8,cross.area=5
                   ,category=c('A','B'),lwd=rep(1,1),lty=rep(2,2)
                   ,col=c('red','green'),fill=c('red','green')
                   ,cat.col=c('red','green')
                   ,rotation.degree=9)

a <- venn.diagram(list(A=1:10,B=3:8,C=6:9), fill=c("red","green","blue"), alpha=c(0.5,0.5,0.5),   filename=NULL)
grid.draw(a,record=TRUE);
for(i in 1:5)
{
  for(j in 1:5)
  {
    if(i == j)
    {
      next
    }
    current_data <- list(photo_input[i], photo_)
    venn.diagram(x=list( Fly=fly_term,Human=human_term), filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[616], colors()[38]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07), cat.pos=c(60, 300),cat.cex=1,scale=FALSE)
  }
  
}






################################################



species <- c('human','fly','mouse','yeast','arabidopsis')
a<-species
b<-species



library(plotrix)
plot(factor(a,levels=species),factor(b,levels=species))
draw.circle(30,50,30,border="white",col=rgb(0,1,0,0.6))
draw.circle(70,50,30,border="white",col=rgb(0,0,1,0.6))
abline(v=30,col=rgb(64/255,244/255,208/255,1))
abline(v=70,col=rgb(64/255,244/255,208/255,1))
abline(h=50,col=rgb(64/255,244/255,208/255,1))
dev.off()

m = 0
sigma = 1
t = seq(-10,10,by=0.1)
n1 = 1 / sqrt(2 * pi * sigma) * exp(-(t - m)^2/(2*sigma) )
plot(t,n1)
par(new=TRUE)
sigma = 2
n2 = 1 / sqrt(2 * pi * sigma) * exp(-(t - m)^2/(2*sigma) )
plot(t,n2)







library("grid")
library("ggplotify")

p1 <- as.grob(~barplot(1:10))
p2 <- as.grob(expression(plot(rnorm(10))))
p3 <- as.grob(function() plot(sin))

library("vcd")
data(Titanic)
p4 <- as.grob(~mosaic(Titanic))

library("lattice")
data(mtcars)
p5 <- as.grob(densityplot(~mpg|cyl, data=mtcars))
grid.newpage()
grid.draw(p1)
vp = viewport(x=.35, y=.75, width=.35, height=.3)
pushViewport(vp)
grid.draw(p2)
upViewport()
library(ggplot2)
p1 <- as.ggplot(~barplot(1:10)) +
  annotate("text", x = .6, y = .5,
           label = "Hello Base Plot", size = 5,
           color = 'firebrick', angle=45)

p2 <- as.ggplot(expression(plot(rnorm(10))))
p3 <- as.ggplot(function() plot(sin))

p4 <- as.ggplot(~mosaic(Titanic))

p5 <- as.ggplot(densityplot(~mpg|cyl, data=mtcars))
library(cowplot)

library(colorspace)
col <- rainbow_hcl(3)
names(col) <- unique(iris$Species)

color <- col[iris$Species]
p6 <- as.ggplot(~plot(iris$Sepal.Length, iris$Sepal.Width, col=color, pch=15))

p7 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, color=Species)) +
  geom_point(shape=15) + scale_color_manual(values=col, name="")

legend <- get_legend(p7)

## also able to annotate base or other plots using ggplot2
library(ggimage)
p8 <- p6 + geom_subview(x=.7, y=.78, subview=legend)


p9 <- as.ggplot(~image(volcano))

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3, labels=LETTERS[1:9])








###########################################gene
mouse_data <- read.csv("E:/GO_layers/photo/result/mouse/mouse_gene.csv")
mouse_gene <- c()
for(i in 1:ncol(mouse_data))
{
  mouse_gene <- c(mouse_gene, mouse_data[,i])
}
mouse_gene <- unique(mouse_gene) #6398

human_data <- read.csv("E:/GO_layers/photo/result/human/human_gene.csv")
human_gene <- c()
for(i in 1:ncol(human_data))
{
  human_gene <- c(human_gene, human_data[,i])
}
human_gene <- unique(human_gene) #7941


library(biomaRt)
mart = useMart('ensembl')
listDatasets(mart)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
data <- read.csv("E:/GO_layers/tongyuan/CYCLING GENES.csv")
genes <- as.vector(data$CYCLING.GENE)

#human to mouse
genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
               
               values = genes ,mart = mouse,
               
               attributesL = c("hgnc_symbol","chromosome_name", "start_position"),
               
               martL = human,
               
               uniqueRows=T)

genes
write.csv(genes, "E:/GO_layers/tongyuan/result.csv")
human_gene <- read.csv("E:/GO_layers/tongyuan/result.csv")
human_gene <- human_gene$MGI.symbol
human_gene <- unique(human_gene)#6346
mouse_gene <- bitr(mouse_gene, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL
mouse_gene <- unique(mouse_gene) #6393


fly_data <- read.csv("E:/GO_layers/photo/result/fly/fly_gene.csv")
fly_gene <- fly_data[,2]
fly <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")
genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
               
               values = toupper(fly_gene) ,mart = mouse,
               
               attributesL = c("chromosome_name", "start_position"),
               
               martL = fly,
               
               uniqueRows=T)
fly_gene <- genes$MGI.symbol#3127->317

yeast_data <- read.csv("E:/GO_layers/photo/result/yeast/yeast_gene.csv.csv")
yeast_gene <- yeast_data$x
yeast_gene1 <- bitr(yeast_gene, fromType="ORF", toType="ENTREZID", OrgDb="org.Sc.sgd.db")$ENTREZID
yeast = useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
genes = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id",
               
               values = yeast_gene1 ,mart = mouse,
               
               attributesL = c("chromosome_name", "start_position"),
               
               martL = yeast,
               
               uniqueRows=T)
yeast_gene <- genes$NCBI.gene.ID #4026 ->0


arabid_data <- read.csv("E:/GO_layers/photo/result/Arabidopsis/arabidopsis_gene.csv.csv")
arabid_gene <- arabid_data$x
arabid_gene1 <- bitr(arabid_gene, fromType="TAIR", toType="ENTREZID", OrgDb="org.At.tair.db")$ENTREZID

mart = useMart(biomart="plants_mart",host="plants.ensembl.org")
arabid = useDataset("athaliana_eg_gene",mart= mart)

genes = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id",
               
               values = toupper(yeast_gene1) ,mart = mouse,
               
               attributesL = c("chromosome_name", "start_position"),
               
               martL = arabid,
               
               uniqueRows=T)




#######################################################Supply
mouse_data <- read.csv("E:/GO_layers/photo/result/mouse/mouse_gene.csv")
mouse_list <- list()
col_name <- colnames(mouse_data)
for(i in 1:(ncol(mouse_data) - 1))
{
  for(j in (i + 1):ncol(mouse_data))
  {
    input_list <- list(na.omit(mouse_data[,i]), na.omit(mouse_data[,j]))
    names(input_list) <- c(col_name[i], col_name[j])
    
    b1 <- venn.diagram(x=input_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.cex=1)
    #grid.draw(b1)
    mouse_list <- c(mouse_list, list(b1))
  }
}
plot_grid(plotlist = mouse_list)




human_data <- read.csv("E:/GO_layers/photo/result/human/human_gene.csv")
human_list <- list()
col_name <- colnames(human_data)
for(i in 1:(ncol(human_data) - 1))
{
  for(j in (i + 1):ncol(human_data))
  {
    input_list <- list(na.omit(human_data[,i]), na.omit(human_data[,j]))
    names(input_list) <- c(col_name[i], col_name[j])
    
    b1 <- venn.diagram(x=input_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.cex=1)
    #grid.draw(b1)
    human_list <- c(human_list, list(b1))
  }
}

plot_grid(plotlist = human_list)

























# ###################################exp amp
# 
# #mouse
# file_name <- list.files(path = "E:/GO_layers/circadian_data/mouse/allgenes12tissues/", pattern = NULL, all.files = FALSE,full.names = FALSE, recursive = T,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
# plot_list <- list()
# for(i in 1:length(file_name))
# {
#   current_file <- paste("E:/GO_layers/circadian_data/mouse/allgenes12tissues/", file_name[i], sep="")
#   temp_data <- read.csv(current_file)
#   temp_data <- temp_data[temp_data$BH.Q<0.05,]
#   input_data <- data.frame(temp_data[,7],rowMeans(temp_data[,8:30]))
#   colnames(input_data) <- c("amp","exp")
#   temp_plot <- ggplot(input_data,aes(x=exp,y=amp)) + geom_point()
#   plot_list <- c(plot_list, list(temp_plot))
#   
# }
# 
# plot_grid(plotlist = plot_list,labels = file_name)
# 
# 
# 
# #human
# bcell_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/4tissue_exp/Bcell.csv")
# skin_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/4tissue_exp/skin.csv")
# bcell_input_data <- bcell_data[,c(3,10)]
# skin_input_data <- skin_data[,c(6,3)]
# colnames(bcell_input_data) <- c("amp","exp")
# colnames(skin_input_data) <- c("amp","exp")
# bcell_plot <- ggplot(bcell_input_data,aes(x=exp,y=amp)) + geom_point()
# skin_plot <- ggplot(skin_input_data,aes(x=exp,y=amp)) + geom_point()
# plot_grid(plotlist = list(bcell_plot, skin_plot),labels = c("Bcell","skin"))
# 
# 
# #fly
# old_data <- read.csv("E:/GO_layers/circadian_data/fly/old_pvalue.csv")
# young_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_youngARSER_output.csv")
# old_input_data <- data.frame(old_data[,7],rowMeans(old_data[,14:25]))
# colnames(old_input_data) <- c("amp","exp")
# old_plot <- ggplot(old_input_data,aes(x=exp,y=amp)) + geom_point()
# young_input_data <- data.frame(young_data[,7],rowMeans(young_data[,14:25]))
# colnames(young_input_data) <- c("amp","exp")
# young_plot <- ggplot(young_input_data,aes(x=exp,y=amp)) + geom_point()
# plot_grid(plotlist = list(old_plot, young_plot),labels = c("old","young"))
# 
# 
# 
# 
# #yeast
# yeast_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/amp_exp.csv")
# yeast_input_data <- data.frame(yeast_data[,6],rowMeans(yeast_data[,7:26]))
# colnames(yeast_input_data) <- c("amp","exp")
# ggplot(yeast_input_data,aes(x=exp,y=amp)) + geom_point()
# 
# 
# #arabidopsis
# arabid_data <- read.csv("E:/GO_layers/circadian_data/Arabidopsis/GSEPLANTS/arab_amp_exp.csv")
# arabid_input_data <- data.frame(arabid_data[,6],rowMeans(arabid_data[,7:24]))
# colnames(arabid_input_data) <- c("amp","exp")
# ggplot(arabid_input_data,aes(x=exp,y=amp)) + geom_point()
# 
# 
# 
# 
# 
# mouse_data <- read.csv("E:/GO_layers/blast/circidian_gene/mouse_gene.csv")
# write.csv(mouse_data, "E:/GO_layers/blast/circidian_gene/mouse_gene_id.csv")
# mouse_new_data <- bitr(mouse_Data$ADR, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL
# write.csv(mouse_new_data, "E:/GO_layers/blast/circidian_gene/mouse_gene.csv")
#######################################################################orthlogy gene
species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','cyanobacteria','archaea')
gene_num_list <- data.frame(0,0)
colnames(gene_num_list) <- c("species","gene_num")

venn_data <- data.frame(0,0,0,0,0)
colnames(venn_data) <- c("species1","species2","num1","num2","overlap")
for(i in 1:(length(species) - 1))
{
  for(j in (i+1):length(species))
  {
    current_file_name <- paste(paste(paste(paste("E:/GO_layers/blast/ensembl_gene/",species[i],sep=""),'_',sep=""), species[j],sep=""),".sw.blastp",sep="")
    if (file.exists(current_file_name))
    {
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("E:/GO_layers/blast/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("E:/GO_layers/blast/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      
      gene1 <- unique(gene1[,1])
      gene2 <- unique(gene2[,1])
      gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene1))))
      gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene2))))
      orth_gene1 <- gene_list[gene_list[,1]%in%gene1 ,2]
      orth_gene2 <- gene_list[gene_list[,2]%in%gene2 ,2]
      orth_gene1 <- na.omit(orth_gene1)
      orth_gene2 <- na.omit(orth_gene2)
      temp_row <- c(species[i], species[j], length(unique(orth_gene1)), length(unique(orth_gene2)), length(unique(intersect(orth_gene1, orth_gene2))))
      venn_data <- rbind(venn_data, temp_row)
    }else
    {
      current_file_name <- paste(paste(paste(paste("E:/GO_layers/blast/ensembl_gene/",species[j],sep=""),'_',sep=""), species[i],sep=""),".sw.blastp",sep="")
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("E:/GO_layers/blast/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("E:/GO_layers/blast/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      gene1 <- gene1[,1]
      gene2 <- gene2[,1]
      gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene1))))
      gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene2))))
      orth_gene1 <- gene_list[ gene_list[,1]%in%gene1 ,2]
      orth_gene2 <- gene_list[ gene_list[,2]%in%gene2 ,2]
      orth_gene1 <- na.omit(orth_gene1)
      orth_gene2 <- na.omit(orth_gene2)
      temp_row <- c(species[i], species[j], length(unique(orth_gene2)), length(unique(orth_gene1)), length(unique(intersect(orth_gene1, orth_gene2))))
      venn_data <- rbind(venn_data, temp_row)
    }
  }
}
  


gene_num_list <- gene_num_list[-1,]
gene_num_list <- unique(gene_num_list)
gene_num_list$gene_num <- as.numeric(gene_num_list$gene_num)
venn_data_new <- venn_data[2:nrow(venn_data),]
plot_list = list()
for(i in 1:nrow(venn_data_new))
{
  data1_all = as.numeric(venn_data_new[i,3])
  data2_all = as.numeric(venn_data_new[i,4])
  common_num = as.numeric(venn_data_new[i,5])
  temp_plot <- draw.pairwise.venn(area1=data1_all,area2=data2_all,cross.area=common_num,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
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










##################################################circadian gene in different tissues


b1 <- venn.diagram(x=list(A=data1,B=data2),filename=NULL, height = 450, width = 450, ext.text=FALSE,offset=0,resolution =300, euler.d = TRUE,imagetype="png", col=c('white','white'), fill=c("cornflowerblue", "green"), alpha=c(0.6, 0.6), cex='white',cat.cex='white')
T<-draw.pairwise.venn(area1=10,area2=100,cross.area = 5,filename=NULL,inverted = 1,
                       ext.text=FALSE, imagetype="png", col=c('white','white'), fill=c("cornflowerblue", "green"))

T<-draw.pairwise.venn(area1=100,area2=10,cross.area = 5,filename=NULL,inverted = 0,
                      height = 450, width = 450, ext.text=FALSE,offset=0,resolution =300, euler.d = TRUE,imagetype="png", col=c('white','white'), fill=c("cornflowerblue", "green"), alpha=c(0.6, 0.6), cex=1,cat.cex=1)



















#mouse
mouse_data <- read.csv("E:/GO_layers/photo/result/mouse/mouse_gene_phototdata.csv")
mouse_data <- read.csv("E:/GO_layers/photo/result/mouse/mouse_gene.csv")
tissue <- colnames(mouse_data)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- na.omit(mouse_data[,i])
    data2 <- na.omit(mouse_data[,j])
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    small_flag = -1
    if(length(data1) > length(data2))
    {
      small_flag = 0
    }else
    {
      small_flag = 1
    }
    b1 <- venn.diagram(x=list(A=data1,B=data2),filename=NULL, height = 450, width = 450, inverted = small_flag, ext.text=FALSE,offset=0,resolution =300, imagetype="png", cex=0, cat.cex=0, col=c('white','white'), fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6))
    #b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c("cornflowerblue", "green"), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
    #b1 <- draw.pairwise.venn(area1=length(temp_list[[1]]),area2=length(temp_list[[2]]),cross.area=length(intersect(temp_list[[1]],temp_list[[2]])),filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c("cornflowerblue", "green"), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=0)
    
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:144,byrow=TRUE,nrow=12)
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=12, ncol=12)





                                                                                                                                                                                 
p.len <- viewport(x=0, y=0.66, width=0.66, height=0.34, just=c("left", "bottom"))
vp.wid <- viewport(x=0.66, y=0, width=0.34, height=0.66, just=c("left", "bottom"))
vp.scatter <- viewport(x=0, y=0, width=0.66, height=0.66, just=c("left", "bottom"))
vp.scatter1 <- viewport(x=0.11, y=0.09, width=0.50, height=0.50, just=c("left", "bottom"))


species = tissue
species1 = rev(tissue)
print(as.ggplot(expression(plot(1:length(tissue),1:length(tissue),type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:length(tissue),labels = species) + axis(side = 4,at = 1:length(tissue),labels = species1))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=length(tissue), ncol=length(tissue), axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)








#human
mouse_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/4tissue_gene.csv")
tissue <- colnames(mouse_data)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- na.omit(mouse_data[,i])
    data2 <- na.omit(mouse_data[,j])
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=1)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(ncol(mouse_data)*ncol(mouse_data)),byrow=TRUE,nrow=ncol(mouse_data))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=ncol(mouse_data), ncol=ncol(mouse_data))



#fly
mouse_data <- read.csv("E:/GO_layers/circadian_data/fly/fly_gene.csv")
tissue <- colnames(mouse_data)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- na.omit(mouse_data[,i])
    data2 <- na.omit(mouse_data[,j])
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=1)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(ncol(mouse_data)*ncol(mouse_data)),byrow=TRUE,nrow=ncol(mouse_data))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=ncol(mouse_data), ncol=ncol(mouse_data))



#yeast
mouse_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/gene_two.csv")
tissue <- colnames(mouse_data)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- na.omit(mouse_data[,i])
    data2 <- na.omit(mouse_data[,j])
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=1)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(ncol(mouse_data)*ncol(mouse_data)),byrow=TRUE,nrow=ncol(mouse_data))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=ncol(mouse_data), ncol=ncol(mouse_data))


##################################################################term in different tissues


#mouse
#mouse_data <- read.csv("E:/GO_layers/expressed_gene_enrich/mouse_result_padj.csv")
mouse_data <- read.csv("E:/GO_layers/circadian_data/mouse/mosue_notlevel_bp.csv")
temp_plot_list <- list()
num <- 0
mouse_data <- mouse_data[mouse_data$ont=='MF',]
tissue <- unique(mouse_data$Cluster)
tissue
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- mouse_data[mouse_data$Cluster==tissue[i],'ID']
    data2 <- mouse_data[mouse_data$Cluster==tissue[j],'ID']
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    small_flag = -1
    if(length(data1) > length(data2))
    {
      small_flag = 0
    }else
    {
      small_flag = 1
    }
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", ext.text=FALSE,offset=0,inverted = small_flag, col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=0, cat.dist=c(-0.07, -0.07),cat.pos=c(300, 300),cat.cex=0)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(length(tissue)*length(tissue)),byrow=TRUE,nrow=length(tissue))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 

p.len <- viewport(x=0, y=0.66, width=0.66, height=0.34, just=c("left", "bottom"))
vp.wid <- viewport(x=0.66, y=0, width=0.34, height=0.66, just=c("left", "bottom"))
vp.scatter <- viewport(x=0, y=0, width=0.66, height=0.66, just=c("left", "bottom"))
vp.scatter1 <- viewport(x=0.11, y=0.09, width=0.50, height=0.50, just=c("left", "bottom"))


species = tissue
species1 = rev(tissue)
print(as.ggplot(expression(plot(1:length(tissue),1:length(tissue),type = 'n',xaxt= 'n',yaxt='n',bty = '7',xlab = NA,ylab = NA) + axis(side = 3,at = 1:length(tissue),labels = species) + axis(side = 4,at = 1:length(tissue),labels = species1))), vp=vp.scatter)
print(plot_grid(plotlist = plot_data1,nrow=length(tissue), ncol=length(tissue), axis="bl",rel_widths=c(1,1,1)), vp=vp.scatter1)






#human
mouse_data <- read.csv("E:/GO_layers/expressed_gene_enrich/human_result_padj.csv")
tissue <- unique(mouse_data$Cluster)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- mouse_data[mouse_data$Cluster==tissue[i],'ID']
    data2 <- mouse_data[mouse_data$Cluster==tissue[j],'ID']
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=1)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(length(tissue)*length(tissue)),byrow=TRUE,nrow=length(tissue))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=length(tissue), ncol=length(tissue))



#fly
mouse_data <- read.csv("E:/GO_layers/expressed_gene_enrich/fly_result_padj.csv")
mouse_data <- mouse_data[mouse_data$ont=="BP",]
tissue <- unique(mouse_data$Cluster)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- mouse_data[mouse_data$Cluster==tissue[i],'ID']
    data2 <- mouse_data[mouse_data$Cluster==tissue[j],'ID']
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=1)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(length(tissue)*length(tissue)),byrow=TRUE,nrow=length(tissue))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=length(tissue), ncol=length(tissue))



#yeast
mouse_data <- read.csv("E:/GO_layers/expressed_gene_enrich/yeast_result_padj.csv")
mouse_data <- mouse_data[mouse_data$ont=="BP",]
tissue <- unique(mouse_data$Cluster)
temp_plot_list <- list()
num <- 0
for(i in 1:(length(tissue) - 1))
{
  for(j in (i + 1):length(tissue))
  {
    num <- num + 1
    data1 <- mouse_data[mouse_data$Cluster==tissue[i],'ID']
    data2 <- mouse_data[mouse_data$Cluster==tissue[j],'ID']
    temp_list  = list(data1, data2)
    names(temp_list) <- c(tissue[i],tissue[j])
    b1 <- venn.diagram(x=temp_list,filename=NULL, height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[38], colors()[616]), alpha=c(0.6, 0.6), lwd=c(1, 1), cex=1, cat.dist=c(-0.07, -0.07),cat.pos=c(60, 300),cat.cex=1)
    temp_plot_list <- c(temp_plot_list , list(b1))
  }
}
plot_data1 <- list()
temp_matrix <- matrix(1:(length(tissue)*length(tissue)),byrow=TRUE,nrow=length(tissue))
temp_num <- temp_matrix[upper.tri(temp_matrix)]
temp_num <- sort(temp_num)
plot_data1[c(temp_num)] = temp_plot_list 
plot_grid(plotlist=plot_data1,nrow=length(tissue), ncol=length(tissue))

































































#####################################################amp and exp
#######################################mouse
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
           fill = c(colors()[38], colors()[616]),cex=1.5,bty='n')
    par(mar=parraw_mar,las= 0)
  }
}



legend('top',c("Mean","Variation","Range","rAMP"),
      fill = c("red1", 'limegreen',"mediumblue","gold"),cex=1.5,bty='n')





#######################################################human
aa= matrix(0,nrow = 2,ncol = 3)





#blood
blood_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/blood/DATA.csv",header=TRUE,row.names = 1)
mean_exp <- as.numeric(rowMeans(blood_data))
blood_data1 <- read.csv("E:/GO_layers/circadian_data/human/4tissue/blood/1/meta3dSubjectID_AF0089SleepExtension_DESIGN.csv",header=TRUE,row.names = 1)
blood_data <- cbind(blood_data, mean_exp)
cir = blood_data1[,'JTK_pvalue']<0.05
expr_cir = blood_data1[cir,]
expr_cir <- expr_cir[expr_cir[,'JTK_amplitude']!=0,]
#expr_cirmean= expr_cir[,'mean_exp']
expr_cirmean <- rowMeans(expr_cir[,9:18])
plot(log(expr_cirmean),log(expr_cir[,'JTK_amplitude']),xlab = '',ylab = '',main='Blood')
aa[1,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'JTK_amplitude'])))$estimate
aa[1,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(expr_cir[,'JTK_amplitude'])))$p.value
aa[1,3]= aa[1,1]^2
bb[1] = sum(cir)
#bcell
bcell_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/b-cell/bcell_data.csv",header=TRUE,row.names = 1)
expr_cirmean <- c()
for(i in 1:nrow(bcell_data))
{
  expr_cirmean[length(expr_cirmean) + 1] <- mean(as.numeric(bcell_data[i,16:21]))
}
plot(log(expr_cirmean),log(bcell_data[,'Max_Amp']),xlab = '',ylab = '',main='Bcell')
aa[2,1]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(bcell_data[,'Max_Amp'])))$estimate
aa[2,2]=cor.test(log(as.numeric(expr_cirmean)),log(as.numeric(bcell_data[,'Max_Amp'])))$p.value
aa[2,3]= aa[2,1]^2
#bb[2] = sum(cir)

colnames(aa)= c('cor','pvalue','R2')









aaa= matrix(0,nrow = 2,ncol = 2)
rownames(aaa)= c("Blood", "Bcell")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']





par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/blood/DATA.csv",header=TRUE,row.names = 1)
mean_exp <- as.numeric(rowMeans(blood_data))
blood_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/blood/1/meta3dGroupID_SleepExtension_DESIGN.csv",header=TRUE,row.names = 1)
blood_data <- cbind(blood_data, mean_exp)
cir = blood_data[,'meta3d_Pvalue']<0.05
expr_cir = blood_data[cir,]
expr_cir <- expr_cir[expr_cir[,'meta3d_AMP']!=0,]
expr_cirmean= expr_cir[,'mean_exp']
plot(log(expr_cirmean),log(expr_cir[,'meta3d_AMP']),xlab = 'expression',ylab = 'amplitude',main='Blood',
     xlim = c(0,4),ylim = c(-2,2),
     pch=20,cex=0.5)

text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)




bcell_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/b-cell/bcell_data.csv",header=TRUE,row.names = 1)
expr_cirmean <- c()
for(i in 1:nrow(bcell_data))
{
  expr_cirmean[length(expr_cirmean) + 1] <- mean(as.numeric(bcell_data[i,16:21]))
}
plot(log(expr_cirmean),log(bcell_data[,'Max_Amp']),xlab = 'expression',ylab = 'amplitude',main='Bcell',
     xlim = c(2,10),ylim = c(-2,2),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[2,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)
par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[2,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[2,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)





  
  
# {
#   parraw_mar= par()$mar
#   par(mar = c(0,0, 0, 0),las = 1)  
#   plot(c(0, 4),c(0, 4),  
#        xlim = c(0, 4), ylim = c(0, 4),  
#        xaxs = 'i', yaxs = 'i',  
#        xaxt = 'n', yaxt = 'n', type = 'n', ann =F, axes = F) 
#   legend('top',c('variation explained by expression','variation explained by other factors'),
#          fill = c('skyblue','white'),cex=1.5,bty='n')
#   par(mar=parraw_mar,las= 0)
# }

####################################################fly
aa= matrix(0,nrow = 2,ncol = 3)





#old
blood_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$fdr_BH<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,15:23]))
amp <- blood_data$amplitude
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='old')
aa[1,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[1,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[1,3]= aa[1,1]^2

#young
bcell_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_youngARSER_output.csv",header=TRUE,row.names = 1)
bcell_data <- bcell_data[bcell_data$fdr_BH<0.05,]
mean_exp <- as.numeric(rowMeans(bcell_data[,15:23]))
amp <- bcell_data$amplitude
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='young')
aa[2,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[2,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[2,3]= aa[1,1]^2

colnames(aa)= c('cor','pvalue','R2')









aaa= matrix(0,nrow = 2,ncol = 2)
rownames(aaa)= c("Old", "Young")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']





par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$fdr_BH<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,15:23]))
amp <- blood_data$amplitude
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Old',
     xlim = c(0,4),ylim = c(-2,2),
     pch=20,cex=0.5)


text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)




bcell_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_youngARSER_output.csv",header=TRUE,row.names = 1)
bcell_data <- bcell_data[bcell_data$fdr_BH<0.05,]
mean_exp <- as.numeric(rowMeans(bcell_data[,15:23]))
amp <- bcell_data$amplitude
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='young',
     xlim = c(2,10),ylim = c(-2,2),
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[2,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)
par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[2,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[2,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)







####################################################yeast
aa= matrix(0,nrow = 2,ncol = 3)
#sample2
#blood_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/JTK.Sample2.csv",header=TRUE,row.names = 1)
blood_data <- read.csv("E:/GO_layers/circadian_data/Yeast/old/JTK.yeast_gb_rpkm_nonCycNorm_Sample2.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$ADJ.P<0.05,]
blood_data <- blood_data[blood_data[,'AMP']!=0,]
#mean_exp <- as.numeric(rowMeans(blood_data[,12:31]))
mean_exp <- as.numeric(rowMeans(blood_data[,6:25]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='Sample2')
aa[1,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[1,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[1,3]= aa[1,1]^2

#sample6
#bcell_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/JTK.Sample6.csv",header=TRUE,row.names = 1)
bcell_data <- read.csv("E:/GO_layers/circadian_data/Yeast/old/data.csv",header=TRUE,row.names = 1)
bcell_data <- bcell_data[bcell_data$ADJ.P<0.05,]
bcell_data <- bcell_data[bcell_data[,'AMP']!=0,]
#mean_exp <- as.numeric(rowMeans(bcell_data[,6:29]))
mean_exp <- as.numeric(rowMeans(bcell_data[,9:33]))
amp <- bcell_data$AMP
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='Sample6')
aa[2,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[2,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[2,3]= aa[2,1]^2

colnames(aa)= c('cor','pvalue','R2')


aaa= matrix(0,nrow = 2,ncol = 2)
rownames(aaa)= c("Old", "Young")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']


par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/JTK.Sample2.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$ADJ.P<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:25]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Sample2',xlim=c(-2,10),ylim=c(-2,4),
     pch=20,cex=0.5)


text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)




bcell_data <- read.csv("E:/GO_layers/circadian_data/Yeast/JTK/JTK.Sample6.csv",header=TRUE,row.names = 1)
bcell_data <- bcell_data[bcell_data$ADJ.P<0.05,]
mean_exp <- as.numeric(rowMeans(bcell_data[,6:29]))
amp <- bcell_data$AMP
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Sample6',
     pch=20,cex=0.5)
text(2,7.5,paste0('r = ',round(aa[2,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)
par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[2,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[2,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)




####################################################arabidopsis
aa= matrix(0,nrow = 1,ncol = 3)
blood_data <- read.csv("E:/GO_layers/circadian_data/Arabidopsis/GSEPLANTS/JTKresult_expr.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$BH.Q<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:23]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='Arabidopsis')
aa[1,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[1,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[1,3]= aa[1,1]^2



colnames(aa)= c('cor','pvalue','R2')


aaa= matrix(0,nrow = 1,ncol = 2)
rownames(aaa)= c("Arabidopsis")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']


par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/Arabidopsis/GSEPLANTS/JTKresult_expr.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$BH.Q<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:23]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Arabidopsis',
     pch=20,cex=0.5)


text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)


####################################################neurospora
aa= matrix(0,nrow = 1,ncol = 3)
blood_data <- read.csv("E:/GO_layers/circadian_data/Neurospora/12915_2015_126_MOESM3_ESM.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$pvalue<0.05,]
mean_exp <- blood_data$mean
amp <- blood_data$amplitude
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='Neurospora')
aa[1,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[1,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[1,3]= aa[1,1]^2



colnames(aa)= c('cor','pvalue','R2')


aaa= matrix(0,nrow = 1,ncol = 2)
rownames(aaa)= c("Neurospora")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']


par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/Neurospora/12915_2015_126_MOESM3_ESM.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$pvalue<0.05,]
mean_exp <- blood_data$mean
amp <- blood_data$amplitude
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Neurospora',
     pch=20,cex=0.5)


text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)


####################################################cyanobacteria
aa= matrix(0,nrow = 1,ncol = 3)
blood_data <- read.csv("E:/GO_layers/circadian_data/cyanobacteria/JTK.LLREPLICATE2.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$BH.Q<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:25]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='Cyanobacteria')
aa[1,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[1,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[1,3]= aa[1,1]^2



colnames(aa)= c('cor','pvalue','R2')


aaa= matrix(0,nrow = 1,ncol = 2)
rownames(aaa)= c("Cyanobacteria")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']


par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/cyanobacteria/JTK.LLREPLICATE2.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$BH.Q<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:25]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Cyanobacteria',
     pch=20,cex=0.5)


text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)





####################################################archaea
aa= matrix(0,nrow = 1,ncol = 3)
blood_data <- read.csv("E:/GO_layers/circadian_data/Archaea/LSresult_data.csv",header=TRUE)
blood_data <- blood_data[blood_data$p<0.2,]
blood_data <- read.csv("E:/GO_layers/circadian_data/Archaea/meta2d_data.csv",header=TRUE)
gene_name <- read.csv("E:/GO_layers/circadian_data/Archaea/gene.csv")
gene_name <- gene_name[,1]
blood_data1 <- blood_data[blood_data[,1] %in% gene_name,]


blood_data1 <- blood_data1[blood_data1$LS_pvalue<0.2,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:25]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = '',ylab = '',main='Cyanobacteria')
aa[1,1]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$estimate
aa[1,2]=cor.test(log(as.numeric(mean_exp)),log(as.numeric(amp)))$p.value
aa[1,3]= aa[1,1]^2



colnames(aa)= c('cor','pvalue','R2')


aaa= matrix(0,nrow = 1,ncol = 2)
rownames(aaa)= c("Cyanobacteria")
colnames(aaa)=c('expression','others')
aaa[,1]= aa[,'R2']
aaa[,2]= 1- aa[,'R2']


par(mfcol=c(4,3))


bb=1:12
{
  mat <-rbind(matrix(1:12, 3, 4,byrow = T),13)  
  widths <-c(9,9,9,9)  
  heights <-c(6,6,6,2)  
  layout(mat,widths = widths, heights = heights)  
  nf <-layout(mat,  widths = widths, heights = heights)  
  layout.show(nf)
}



blood_data <- read.csv("E:/GO_layers/circadian_data/cyanobacteria/JTK.LLREPLICATE2.csv",header=TRUE,row.names = 1)
blood_data <- blood_data[blood_data$BH.Q<0.05,]
mean_exp <- as.numeric(rowMeans(blood_data[,6:25]))
amp <- blood_data$AMP
plot(log(mean_exp),log(amp),xlab = 'expression',ylab = 'amplitude',main='Cyanobacteria',
     pch=20,cex=0.5)


text(2,7.5,paste0('r = ',round(aa[1,1]*1000)/1000,'\n','p < 0.001'),pos = 4,cex=1.5)


par(plt=c(0.60,0.95,0.20,0.55))
par(new=T)
pie(aaa[1,],labels = NA,col= c('skyblue','white'),init.angle = 90)
text(-0.7,-0.3,paste0(round(aaa[1,1]*10000)/100,'%'),cex = 1,pos=4)
par(plt=c(0.15,0.85,0.25,0.65))

parmai=par()$mai
par(new= T)
par(mai=rep(0,4))
plot.new()
text(0,1,letters[i],cex=2)
par(mai=parmai)
















# proportion and expression levels -------------------
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

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:12){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
  }
  # legend
  parmai=par()$mai
  par(mai=rep(0,4))
  plot.new()
  plot.new()
  lbls= paste((0:4)*2,"0%-",(0:4)*2+2,"0%",sep = "")
  legend("topleft",lbls,
         fill= c(colors[589],colors()[590],colors()[591],colors()[592],colors()[593]),
         cex=1.5,bty = 'n')
  par(mai=parmai)
  # dev.off()
}





#human
aa= matrix(-1,nrow = 2,ncol = 5)
rownames(aa)=c('blood','bcell')
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 2,ncol = 4)
rownames(diffpvalue)=c('blood','bcell')
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)

#

#blood
blood_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/blood/DATA.csv",header=TRUE,row.names = 1)
mean_exp <- order(as.numeric(rowMeans(blood_data)))
blood_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/blood/1/meta3dGroupID_SleepExtension_DESIGN.csv",header=TRUE,row.names = 1)
blood_data <- cbind(blood_data, mean_exp)
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
bcell_data <- read.csv("E:/GO_layers/circadian_data/human/4tissue/b-cell/bcell_data.csv",header=TRUE,row.names = 1)
expr_cirmean <- c()
for(i in 1:nrow(bcell_data))
{
  expr_cirmean[length(expr_cirmean) + 1] <- mean(as.numeric(bcell_data[i,16:21]))
}
i=2
bhp = bcell_data[,'P']
bhp=bhp[order(rowMeans(bcell_data[,16:21]))]
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



tissuenamefull <- c("blood","bcell")

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/humanrising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/humanrising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:2){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
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




#old
blood_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_old_ARSER_output.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'fdr_BH']
bhp=bhp[order(rowMeans(blood_data[,15:23]))]
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
blood_data <- read.csv("E:/GO_layers/circadian_data/fly/GSE81100_youngARSER_output.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'fdr_BH']
bhp=bhp[order(rowMeans(blood_data[,15:23]))]
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
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:2){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
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
bhp = blood_data[,'BH.Q']
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
bhp = blood_data[,'BH.Q']
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
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:2){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
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




#arabidopsis
aa= matrix(-1,nrow = 1,ncol = 5)
tissuename <- c("Arabidopsis")
rownames(aa)=tissuename
colnames(aa)= c('level1','level2','level3','level4','level5')
bb= aa
diffpvalue=matrix(-1,nrow = 1,ncol = 4)
rownames(diffpvalue)=tissuename
colnames(diffpvalue)=c('diff12','diff23','diff34','diff45')
set.seed(1)




blood_data <- read.csv("E:/GO_layers/circadian_data/Arabidopsis/GSEPLANTS/JTKresult_expr.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'BH.Q']
bhp=bhp[order(rowMeans(blood_data[,6:23]))]
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

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/arabidopsis_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/arabidopsis_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:1){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
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



blood_data <- read.csv("E:/GO_layers/circadian_data/Neurospora/12915_2015_126_MOESM3_ESM.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'pvalue']
bhp=bhp[order(blood_data[,'mean'])]
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
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:1){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
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




blood_data <- read.csv("E:/GO_layers/circadian_data/cyanobacteria/JTK.LLREPLICATE2.csv",header=TRUE,row.names = 1)
bhp = blood_data[,'BH.Q']
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


tissuenamefull <- tissuename

write.csv(aa,'E:/GO_layers/circadian_data/trend_percentage/cyanobacteria_rising_trend_percentage.csv')

{
  aa= read.csv('E:/GO_layers/circadian_data/trend_percentage/cyanobacteria_rising_trend_percentage.csv',row.names = 1)
  aa= aa/rowSums(aa)
  filcolor= c("yellow","green","green3","skyblue","blue")
  # png("fig1.png",width=10,height=10,units="in",res=350)
  par(mfrow=c(4,4))
  for(i in 1:1){
    aan= as.numeric(aa[i,])
    barplot(aan,
            main = tissuenamefull[i],
            ylim = c(0,0.6),
            ylab = "distribution of cycling genes",
            col = filcolor,
            xlab = "cost of mRNA")
    
    parmai=par()$mai
    par(new= T)
    par(mai=rep(0,4))
    plot.new()
    text(0,1,letters[i],cex=2)
    par(mai=parmai)
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
